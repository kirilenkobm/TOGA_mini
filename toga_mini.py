#!/usr/bin/env python3
"""Master script for the TOGA pipeline.

Perform all operations from the beginning to the end.
If you need to call TOGA: most likely this is what you need.
"""
import argparse
import json
import os
import shutil
import sys
from datetime import datetime as dt
from typing import Optional

from toga_modules.constants import Constants
from toga_modules.bed_hdf5_index import bed_hdf5_index
from toga_modules.chain_bst_index import chain_bst_index
from toga_modules.classify_chains import classify_chains
from toga_modules.common import call_process
from toga_modules.common import get_fst_col
from toga_modules.common import make_symlink
from toga_modules.common import setup_logger
from toga_modules.common import to_log
from toga_modules.common import get_shared_lib_extension

from toga_modules.filter_bed import prepare_bed_file
from toga_modules.make_pr_pseudogenes_annotation import create_processed_pseudogenes_track
from toga_modules.merge_chains_output import merge_chains_output
from toga_modules.create_orthologous_loci_table import create_orthologous_loci_table
from toga_modules.stitch_fragments import stitch_scaffolds
from toga_modules.toga_sanity_checks import TogaSanityChecker

__author__ = "Bogdan M. Kirilenko"
__github__ = "https://github.com/kirilenkobm"
__version__ = "Mini"

LOCATION = os.path.dirname(__file__)


class Toga:
    """TOGA manager class."""
    def __init__(self, args):
        """Initiate toga class."""
        self.t0 = dt.now()
        # define the project name
        if args.project_dir:
            _dirname = os.path.dirname(args.project_dir)
            self.project_name = os.path.basename(_dirname)
        else:
            self.project_name = self.generate_project_name()
        # create a project dir
        self.wd: str = (  # had to add annotation to supress typing warnings in PyCharm 2023.3
            os.path.abspath(args.project_dir)
            if args.project_dir
            else os.path.join(os.getcwd(), self.project_name)
        )
        os.mkdir(self.wd) if not os.path.isdir(self.wd) else None

        # manage logfiles
        _log_filename = self.t0.strftime("%Y_%m_%d_at_%H_%M")
        self.quiet = args.quiet
        self.log_file = os.path.join(self.wd, f"toga_{_log_filename}.log")
        self.log_dir = os.path.join(self.wd, "temp_logs")  # temp file to collect logs from processes
        os.mkdir(self.log_dir) if not os.path.isdir(self.log_dir) else None
        setup_logger(self.log_file, write_to_console=not self.quiet)

        # check if all files TOGA needs are here
        self.temp_files = []  # remove at the end, list of temp files
        to_log("#### Initiating TOGA class ####")
        self.max_workers = args.max_workers

        self.toga_exe_path = os.path.dirname(__file__)
        self.version = "TOGA-mini-dev"
        TogaSanityChecker.check_args_correctness(self, args)
        self.__modules_addr()
        self.temp_wd = os.path.join(self.wd, Constants.TEMP)
        self.project_name = self.project_name.replace("/", "")
        os.mkdir(self.temp_wd) if not os.path.isdir(self.temp_wd) else None

        # to avoid crash on filesystem without locks:
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        chain_basename = os.path.basename(args.chain_input)
        # dir to collect log files with rejected reference genes:
        self.rejected_dir = os.path.join(self.wd, "skipped_transcripts_logs")
        os.mkdir(self.rejected_dir) if not os.path.isdir(self.rejected_dir) else None

        # filter chain in this folder
        g_ali_basename = "genome_alignment"
        self.chain_file = os.path.join(self.temp_wd, f"{g_ali_basename}.chain")
        # there is an assumption that chain file has .chain extension
        # chain indexing was a bit problematic: (i) bsddb3 fits perfectly but is very
        # painful to install, (ii) sqlite is also fine but might be dysfunctional on some
        # cluster file systems, so we create chain_ID: (start_byte, offset) dictionary for
        # instant extraction of a particular chain from the chain file
        # we save these dictionaries into two files: a text file (tsv) and binary file with BST
        # depending on the case we will use both (for maximal performance)
        self.chain_index_file = os.path.join(self.temp_wd, f"{g_ali_basename}.bst")
        self.chain_index_txt_file = os.path.join(
            self.temp_wd, f"{g_ali_basename}.chain_ID_position"
        )

        # make the command, prepare the chain file
        if not os.path.isfile(args.chain_input):
            self.die(f"Error! File {args.chain_input} doesn't exist!")

        if chain_basename.endswith(".gz"):  # version for gz
            chain_filter_cmd = (
                f"gzip -dc {args.chain_input} | "
                f"{self.CHAIN_SCORE_FILTER} stdin "
                f"{args.min_score} > {self.chain_file}"
                # Tried to replace C binary with AWK, something to think about
                # f"awk -f {self.CHAIN_SCORE_FILTER_AWK} {args.min_score} "
                # f"> {self.chain_file}"
            )
        elif args.no_chain_filter:  # it is .chain and score filter is not required
            chain_filter_cmd = f"rsync -a {args.chain_input} {self.chain_file}"
        else:  # it is .chain | a score filter required
            chain_filter_cmd = (
                f"{self.CHAIN_SCORE_FILTER} {args.chain_input} "
                f"{args.min_score} > {self.chain_file}"
            )

        # filter chains with score < threshold
        call_process(
            chain_filter_cmd, "Please check if you use a proper chain file."
        )

        # bed define bed files addresses
        self.ref_bed = os.path.join(self.temp_wd, "toga_filtered_reference_annotation.bed")
        self.index_bed_file = os.path.join(self.temp_wd, "toga_filtered_reference_annotation.hdf5")
        
        # chain processing directory
        self.results_dir = os.path.join(self.temp_wd, "chain_classification_results")

        # filter bed file
        bed_filter_rejected_file = "skipped_on_bed_filter_stage.txt"
        bed_filter_rejected = os.path.join(self.rejected_dir, bed_filter_rejected_file)
        # keeping UTRs!
        prepare_bed_file(
            args.bed_input,
            self.ref_bed,
            # TODO: check whether we like to include this parameter
            save_rejected=bed_filter_rejected,
            only_chrom=args.limit_to_ref_chrom,
        )

        # mics things
        self.isoforms_arg = args.isoforms if args.isoforms else None
        self.isoforms = None  # will be assigned after a completeness check
        self.rejected_log = os.path.join(self.wd, "genes_rejection_reason.tsv")

        # define to call CESAR or not to call
        self.t_2bit = self.__find_two_bit(args.tDB)
        self.q_2bit = self.__find_two_bit(args.qDB)

        self.orthologous_chain_limit = args.orthologous_chain_limit
        self.o2o_only = args.o2o_only
        self.use_parallel_loci = not args.disable_parallel_loci

        self.fragmented_genome = False if args.disable_fragments_joining else True
        self.orthology_score_threshold = args.orthology_score_threshold
        if self.orthology_score_threshold < 0.0 or args.orthology_score_threshold > 1.0:
            self.die(
                "orth_score_threshold parameter must be in range [0..1], got "
                f"{self.orthology_score_threshold}; Abort"
            )

        self.chain_results_df = os.path.join(self.wd, "chain_features.tsv")

        self.bed_fragmented_exons_data = os.path.join(
            self.temp_wd, "bed_fragments_to_exons.tsv"
        )

        # genes to be classified as missing
        self._transcripts_not_intersected = []
        self._transcripts_not_classified = []
        self.predefined_glp_cesar_split = os.path.join(
            self.temp_wd, "predefined_glp_cesar_split.tsv"
        )

        self.__check_param_files()

        # create symlinks to 2bits: let the user know what 2bits were used
        self.t_2bit_link = os.path.join(self.wd, "t2bit.link")
        self.q_2bit_link = os.path.join(self.wd, "q2bit.link")
        make_symlink(self.t_2bit, self.t_2bit_link)
        make_symlink(self.q_2bit, self.q_2bit_link)

        # dump input parameters, object state
        self.toga_params_file = os.path.join(self.temp_wd, "toga_init_state.json")
        self.toga_args_file = os.path.join(self.wd, "project_args.json")
        self.version_file = os.path.join(self.wd, "version.txt")

        with open(self.toga_params_file, "w") as f:
            # default=string is a workaround to serialize a datetime object
            json.dump(self.__dict__, f, default=str)
        with open(self.toga_args_file, "w") as f:
            json.dump(vars(args), f, default=str)
        with open(self.version_file, "w") as f:
            f.write(self.version)

        to_log(f"Saving output to {self.wd}")
        to_log(f"Arguments stored in {self.toga_args_file}")

    def __check_param_files(self):
        """Check that all parameter files exist."""
        files_to_check = [
            self.t_2bit,
            self.q_2bit,
            self.ref_bed,
            self.chain_file,
            self.isoforms_arg,
        ]
        for item in files_to_check:
            if not item:
                # this file just isn't given
                continue
            elif not os.path.isfile(item):
                self.die(f"Error! File {item} not found!")

        # sanity checks: check that bed file chromosomes match reference 2bit
        with open(self.ref_bed, "r") as f:
            lines = [line.rstrip().split("\t") for line in f]
            t_in_bed = set(x[3] for x in lines)
            chromosomes_in_bed = set(x[0] for x in lines)
            # 2bit check function accepts a dict chrom: size
            # from a bed12 file we cannot infer sequence length
            # None is just a placeholder that indicated that we don't need
            # to compare chrom lengths with 2bit
            chrom_sizes_in_bed = {x: None for x in chromosomes_in_bed}
        self.isoforms = TogaSanityChecker.check_isoforms_file(self.isoforms_arg, t_in_bed, self.temp_wd)
        TogaSanityChecker.check_2bit_file_completeness(self.t_2bit, chrom_sizes_in_bed, self.ref_bed)
        # need to check that chain chromosomes and their sizes match 2bit file data
        with open(self.chain_file, "r") as f:
            header_lines = [x.rstrip().split() for x in f if x.startswith("chain")]
            t_chrom_to_size = {x[2]: int(x[3]) for x in header_lines}
            q_chrom_to_size = {x[7]: int(x[8]) for x in header_lines}
        f.close()
        TogaSanityChecker.check_2bit_file_completeness(self.t_2bit, t_chrom_to_size, self.chain_file)
        TogaSanityChecker.check_2bit_file_completeness(self.q_2bit, q_chrom_to_size, self.chain_file)

    def die(self, msg, rc=1):
        """Show msg in stderr, exit with the rc given."""
        to_log(msg)
        to_log(f"Program finished with exit code {rc}\n")
        self.__mark_crashed()
        sys.exit(rc)

    def __modules_addr(self):
        """Define addresses of modules."""
        self.LOCATION = os.path.dirname(__file__)  # folder containing pipeline scripts
        lib_ext = get_shared_lib_extension()
        self.CHAIN_SCORE_FILTER = os.path.join(
            self.LOCATION, Constants.C_BIN_DIR, "chain_score_filter"
        )
        self.CHAIN_COORDS_CONVERT_LIB = os.path.join(
            self.LOCATION, Constants.C_LIB_DIR, f"chain_coords_converter_slib{lib_ext}"
        )

        self.CHAIN_INDEX_SLIB = os.path.join(
            self.LOCATION, Constants.C_LIB_DIR, f"chain_bst_lib{lib_ext}"
        )
        self.BED_BDB_INDEX = os.path.join(
            self.LOCATION, Constants.MODULES_DIR, "bed_hdf5_index.py"
        )

    def __find_two_bit(self, db):
        """Find a 2bit file."""
        if os.path.isfile(db):
            return os.path.abspath(db)
        self.die(f"Two bit file {db} not found! Abort")
        return None

    def run(self, up_to_and_incl: Optional[int] = None) -> None:
        """Run TOGA from start to finish"""
        # 0) preparation:
        # define the project name and mkdir for it,
        # move a chain file filtered, define initial values,
        # make indexed files for the chain
        self.__mark_start()
        to_log("\n\n#### STEP 0: making chain and bed file indexes\n")
        self.__make_indexed_chain()
        self.__make_indexed_bed()
        self.__time_mark("Made indexes")
        if up_to_and_incl is not None and up_to_and_incl == 0: return None

        # 1) make a joblist for chain features extraction
        to_log("\n\n#### STEP 1: Generate extract chain features jobs\n")
        self.__split_chain_jobs()
        self.__time_mark("Split chain jobs")
        if up_to_and_incl is not None and up_to_and_incl == 1: return None

        # 2) extract chain features: parallel process
        to_log("\n\n#### STEP 2: Extract chain features: parallel step\n")
        self.__extract_chain_features()
        self.__time_mark("Chain jobs done")
        to_log(f"Logs from individual chain runner jobs are show below")
        self.__collapse_logs("chain_runner_")
        if up_to_and_incl is not None and up_to_and_incl == 2: return None

        # 3) create chain features dataset
        to_log("\n\n#### STEP 3: Merge step 2 output\n")
        self.__merge_chains_output()
        self.__time_mark("Chains output merged")
        if up_to_and_incl is not None and up_to_and_incl == 3: return None

        # 4) classify chains as orthologous, paralogous, etc. using xgboost
        to_log("\n\n#### STEP 4: Classify chains using gradient boosting model\n")
        self.__classify_chains()
        self.__time_mark("Chains classified")
        if up_to_and_incl is not None and up_to_and_incl == 4: return None

        # 5) create cluster jobs for CESAR2.0
        to_log("\n\n#### STEP 5: Create orthologous loci table")
        self.__create_orthologous_loci_table()
        self.__time_mark("Orthologous loci table created")
        if up_to_and_incl is not None and up_to_and_incl == 5: return None

        # 6) Create a bed track for processed pseudogenes
        to_log("\n\n#### STEP 6: Create processed pseudogenes track\n")
        self.__create_processed_pseudogenes_track()

        # Cleanup
        shutil.rmtree(self.temp_wd) if os.path.isdir(self.temp_wd) else None
        shutil.rmtree(self.log_dir) if os.path.isdir(self.log_dir) else None
        return None

    def __collapse_logs(self, prefix):
        """Merge logfiles starting with the prefix into a single log."""
        log_filenames_with_prefix = [x for x in os.listdir(self.log_dir) if x.startswith(prefix)]
        log_f = open(self.log_file, "a")
        for log_filename in log_filenames_with_prefix:
            full_path = os.path.join(self.log_dir, log_filename)
            clipped_filename = log_filename.split(".")[0]  # remove .log
            in_f = open(full_path, "r")
            for line in in_f:
                log_f.write(f"{clipped_filename}: {line}")
            in_f.close()
        log_f.close()

    def __mark_start(self):
        """Indicate that a TOGA process has started."""
        p_ = os.path.join(self.wd, Constants.RUNNING)
        f = open(p_, "w")
        now_ = str(dt.now())
        f.write(f"TOGA process started at {now_}\n")
        f.close()

    def __mark_crashed(self):
        """Indicate that the TOGA process died."""
        running_f = os.path.join(self.wd, Constants.RUNNING)
        crashed_f = os.path.join(self.wd, Constants.CRASHED)
        os.remove(running_f) if os.path.isfile(running_f) else None
        f = open(crashed_f, "w")
        now_ = str(dt.now())
        f.write(f"TOGA CRASHED AT {now_}\n")
        f.close()

    def __make_indexed_chain(self):
        """Make a chain index file."""
        # make *.bb file
        to_log("Started chain indexing...")
        self.temp_files.append(self.chain_index_file)
        self.temp_files.append(self.chain_file)
        self.temp_files.append(self.chain_index_txt_file)

        if os.path.isfile(self.chain_file) and os.path.isfile(self.chain_index_txt_file):
            return

        chain_bst_index(
            self.chain_file, self.chain_index_file, txt_index=self.chain_index_txt_file
        )

    def __time_mark(self, msg):
        """Left time mark."""
        if self.time_log is None:
            return
        t = dt.now() - self.t0
        with open(self.time_log, "a") as f:
            f.write(f"{msg} at {t}\n")

    def __make_indexed_bed(self):
        """Create gene_ID: bed line bdb indexed file."""
        to_log("Started bed file indexing...")
        bed_hdf5_index(self.ref_bed, self.index_bed_file)
        self.temp_files.append(self.index_bed_file)

    def __split_chain_jobs(self):
        """Create the results directory for parallel processing."""
        os.makedirs(self.results_dir, exist_ok=True)

    def __extract_chain_features(self):
        """Execute extract chain features jobs using direct parallel processing."""
        from toga_modules.parallel_chain_runner import run_parallel_chain_extraction
        
        to_log("Extracting chain features using direct parallel processing")
        
        # Run the parallel chain extraction
        try:
            result_files = run_parallel_chain_extraction(
                chain_file=self.chain_file,
                bed_file=self.ref_bed,
                bed_index=self.index_bed_file,
                results_dir=self.results_dir,
                max_workers=self.max_workers,
                log_file=self.log_file,
                quiet=self.quiet
            )
        except KeyboardInterrupt:
            to_log("Chain feature extraction interrupted by user")
            raise
        
        to_log(f"Chain feature extraction completed. Results in: {result_files}")

        # Store the path to the results for the merge step
        if result_files:
            self.chain_results_file = result_files[0]
        else:
            self.die("No results generated from chain feature extraction!")

    def __merge_chains_output(self):
        """Call parse results."""
        # define where to save intermediate table
        merge_chains_output(
            self.ref_bed, self.isoforms, self.results_dir, self.chain_results_df
        )
        # .append(self.chain_results_df)  -> UCSC plugin needs that

    def __classify_chains(self):
        """Run the decision tree."""
        # define input and output."""
        to_log("Classifying chains")
        self.transcript_to_chain_classes = os.path.join(self.temp_wd, "trans_to_chain_classes.tsv")
        self.pred_scores = os.path.join(self.wd, "orthology_scores.tsv")
        self.se_model = os.path.join(self.LOCATION, "chain_class_models", "se_model.dat")
        self.me_model = os.path.join(self.LOCATION, "chain_class_models", "me_model.dat")
        cl_rej_log = os.path.join(self.rejected_dir, "classify_chains_rejected.txt")

        if not os.path.isfile(self.se_model) or not os.path.isfile(self.me_model):
            self.die(f"Cannot find {self.se_model} or {self.me_model} models!")

        classify_chains(
            self.chain_results_df,
            self.transcript_to_chain_classes,
            self.se_model,
            self.me_model,
            rejected=cl_rej_log,
            raw_out=self.pred_scores,
            annot_threshold=self.orthology_score_threshold,
        )
        # extract not classified transcripts
        # first column in the rejected log
        self._transcripts_not_classified = get_fst_col(cl_rej_log)
        TogaSanityChecker.check_chains_classified(self.chain_results_df)

    def __create_processed_pseudogenes_track(self):
        """Create annotation of processed genes in the query."""
        to_log("Creating processed pseudogenes track.")
        processed_pseudogenes_track = os.path.join(self.wd, "proc_pseudogenes.bed")
        create_processed_pseudogenes_track(
            self.pred_scores, self.chain_file, self.index_bed_file, processed_pseudogenes_track
        )

    def __create_orthologous_loci_table(self):
        """Call create_orthologous_loci_table.py."""
        if self.fragmented_genome:
            to_log("Detecting fragmented transcripts")
            # need to stitch fragments together
            gene_fragments = stitch_scaffolds(
                self.chain_file, self.pred_scores, self.ref_bed, True
            )
            fragmented_dict_file = os.path.join(self.temp_wd, "gene_fragments.txt")
            f = open(fragmented_dict_file, "w")
            for k, v in gene_fragments.items():
                v_str = ",".join(map(str, v))
                f.write(f"{k}\t{v_str}\n")
            f.close()
            to_log(f"Fragments data saved to {fragmented_dict_file}")
        else:
            # no fragment file: ok
            to_log("Skip fragmented genes detection")
            fragmented_dict_file = None

        # create orthologous loci table
        to_log("Creating orthologous loci table")
        
        skipped_path = os.path.join(self.rejected_dir, "skipped.txt")
        self.paralogs_log = os.path.join(self.temp_wd, "paralogs.txt")

        create_orthologous_loci_table(
            orthologs_file=self.transcript_to_chain_classes,
            bed_file=self.ref_bed,
            bdb_bed_file=self.index_bed_file,
            bdb_chain_file=self.chain_index_file,
            target_two_bit=self.t_2bit,
            query_two_bit=self.q_2bit,
            toga_out_dir=self.wd,
            chains_limit=self.orthologous_chain_limit,
            skipped_genes=skipped_path,
            o2o_only=self.o2o_only,
            fragments_data=fragmented_dict_file,
            log_file=self.log_file,
            quiet=self.quiet,
            max_workers=self.max_workers,
            use_parallel=self.use_parallel_loci
        )

    @staticmethod
    def generate_project_name():
        """Generate project name automatically."""
        today_and_now = dt.now().strftime("%Y.%m.%d_at_%H:%M:%S")
        project_name = f"TOGA_project_on_{today_and_now}"
        return project_name


def parse_args(arg_strs: list[str] = None):
    """Parse arguments from the command line, or from a list of strings."""
    app = argparse.ArgumentParser()
    app.add_argument(
        "chain_input",
        type=str,
        help="Chain file. Extensions like FILE.chain or FILE.chain.gz are applicable."
    )
    app.add_argument(
        "bed_input", type=str, help="Bed file with annotations for the target genome."
    )
    app.add_argument(
        "tDB", default=None, help="Reference genome sequence in 2bit format."
    )
    app.add_argument("qDB", default=None, help="Query genome sequence in 2bit format.")

    # global ops
    app.add_argument(
        "project_dir",
        default=None,
        help=(
            "Project directory. TOGA will save all intermediate and output files "
            "exactly in this directory."
        )
    )
    app.add_argument(
        "--min_score",
        "--msc",
        type=int,
        default=15000,
        help=(
            "Chain score threshold. Exclude chains that have a lower score "
            "from the analysis. Default value is 15000."
        )
    )
    app.add_argument(
        "--isoforms", "-i", type=str, default="", help="Path to isoforms data file"
    )
    app.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Don't print to console"
    )
    app.add_argument(
        "--limit_to_ref_chrom",
        default=None,
        help="Find orthologs " "for a single reference chromosome only",
    )
    app.add_argument(
        "--max_workers",
        "--mw",
        type=int,
        default=10,
        help=(
            "Maximum number of workers for parallel processing. "
            "Default is 10."
        )
    )
    app.add_argument(
        "--no_chain_filter",
        "--ncf",
        action="store_true",
        dest="no_chain_filter",
        help=(
            "A flag. Do not filter the chain file (make sure you specified a "
            ".chain but not .gz file in this case)"
        ),
    )
    app.add_argument(
        "--orthology_score_threshold",
        "--ost",
        default=0.5,
        type=float,
        help="Score threshold to distinguish orthologs from paralogs, default 0.5",
    )
    app.add_argument(
        "--orthologous_chain_limit",
        type=int,
        default=100,
        help=(
            "Skip genes that have more that ORTHOLOGOUS_CHAIN_LIMIT orthologous "
            "chains. Recommended values are a 50-100."
        )
    )
    app.add_argument(
        "--o2o_only",
        "--o2o",
        action="store_true",
        dest="o2o_only",
        help="Process only the genes that have a single orthologous chain",
    )
    app.add_argument(
        "--disable_fragments_joining",
        "--dfj",
        dest="disable_fragments_joining",
        action="store_true",
        help="Disable assembling query genes from pieces",
    )
    app.add_argument(
        "--disable_parallel_loci",
        "--dpl",
        dest="disable_parallel_loci",
        action="store_true",
        help="Disable parallel processing for orthologous loci table creation",
    )
    # print help if there are no args
    if not arg_strs and len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)

    args = app.parse_args(arg_strs)
    return args


def main():
    """Entry point."""
    args = parse_args()
    toga_manager = Toga(args)
    toga_manager.run()


if __name__ == "__main__":
    main()
