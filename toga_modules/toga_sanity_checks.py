"""This module contains functions to ensure TOGA arguments correctness."""
import os
import shutil
from itertools import islice
from twobitreader import TwoBitFile

from toga_modules.constants import Constants
from toga_modules.common import to_log
from toga_modules.common import read_isoforms_file

SANITY_CHECKER_PREFIX = "SANITY_CHECKER"
U12_FILE_COLS = 3


class TogaSanityChecker:
    """Utility class to ensure TOGA arguments correctness."""
    @staticmethod
    def check_dir_args_safety(toga_cls, location):
        # No specific directory safety checks needed for ProcessPoolExecutor strategy
        pass
        # TODO: consider other dangerous scenario
        to_log("Does it work?")
        return

    @staticmethod
    def check_args_correctness(toga_cls, args):
        """Check that arguments are correct.

        Error exit if any argument is wrong.
        """
        if not os.path.isfile(args.chain_input):
            toga_cls.die(f"Error! Chain file {args.chain_input} does not exist!")
        if not os.path.isfile(args.bed_input):
            toga_cls.die(f"Error! Bed file {args.bed_input} does not exist!")
        # TODO: consider other scenario

    @staticmethod
    def check_2bit_file_completeness(two_bit_file, chromosome_sizes, chrom_file):
        """Check that the 2bit file is readable."""
        try:  # try to catch EOFError: if 2bitreader cannot read file
            two_bit_reader = TwoBitFile(two_bit_file)
            # check what sequences are in the file:
            twobit_seq_to_size = two_bit_reader.sequence_sizes()
            twobit_sequences = set(twobit_seq_to_size.keys())
            to_log(f"Found {len(twobit_sequences)} sequences in {two_bit_file}")
        except EOFError as err:  # this is a file but twobit reader couldn't read it
            raise ValueError(str(err))
        # another check: that bed or chain chromosomes intersect 2bit file sequences
        check_chromosomes = set(chromosome_sizes.keys())  # chroms in the input file
        intersection = twobit_sequences.intersection(check_chromosomes)
        # chroms_not_in_2bit = check_chroms.difference(twobit_sequences)

        # if len(chroms_not_in_2bit) > 0:
        #     missing_top_100 = list(chroms_not_in_2bit)[:100]
        #     missing_str = "\n".join(missing_top_100)
        #     err = (
        #         f"Error! 2bit file: {two_bit_file}; chain/bed file: {chrom_file}; "
        #         f"Some chromosomes present in the chain/bed file are not found in the "
        #         f"Two bit file. First <=100: {missing_str}"
        #     )
            # raise ValueError(err)
        # check that sizes also match
        for chrom in intersection:
            twobit_seq_len = twobit_seq_to_size[chrom]
            comp_file_seq_len = chromosome_sizes[chrom]
            # if None: this is from bed file: cannot compare
            if comp_file_seq_len is None:
                continue
            if twobit_seq_len == comp_file_seq_len:
                continue
            # got different sequence length in chain and 2bit files
            # which means these chains come from something different
            err = (
                f"Error! 2bit file: {two_bit_file}; chain_file: {chrom_file} "
                f"Chromosome: {chrom}; Sizes don't match! "
                f"Size in twobit: {twobit_seq_len}; size in chain: {comp_file_seq_len}"
            )
            to_log(err)
            raise ValueError(err)

    @staticmethod
    def check_isoforms_file(isoforms_arg, t_in_bed, temp_wd):
        """Sanity checks for isoforms file."""
        if not isoforms_arg:
            to_log("Continue without isoforms file: not provided")
            return None  # not provided: nothing to check
        # isoform file provided: need to check correctness and completeness
        # then check isoform file itself
        _, isoform_to_gene, header = read_isoforms_file(isoforms_arg)
        header_maybe_gene = header[0]  # header is optional, if not the case: first field is a gene
        header_maybe_trans = header[1]  # and the second is the isoform
        # save filtered isoforms file here:  (without unused transcripts)
        isoforms_file = os.path.join(temp_wd, "isoforms.tsv")
        # this set contains isoforms found in the isoforms file
        t_in_i = set(isoform_to_gene.keys())
        # there are transcripts that appear in bed but not in the isoform file
        # if this set is non-empty: raise an error
        u_in_b = t_in_bed.difference(t_in_i)

        if len(u_in_b) != 0:  # isoform file is incomplete
            extra_t_list = "\n".join(
                list(u_in_b)[:100]
            )  # show first 100 (or maybe show all?)
            err_msg = (
                f"Error! There are {len(u_in_b)} transcripts in the bed "
                f"file absent in the isoforms file! "
                f"There are the transcripts (first 100):\n{extra_t_list}"
            )
            raise ValueError(err_msg)

        t_in_both = t_in_bed.intersection(t_in_i)  # isoforms data that we save
        # if header absent: second field found in the bed file
        # then we don't need to write the original header
        # if present -> let keep it
        # there is not absolutely correct: MAYBE there is no header at all, but
        # the first line of the isoforms file is not in the bed file, so we still will write it
        skip_header = header_maybe_trans in t_in_bed

        # write isoforms file
        f = open(isoforms_file, "w")
        if not skip_header:
            f.write(f"{header_maybe_gene}\t{header_maybe_trans}\n")
        to_log(f"Writing isoforms data for {len(t_in_both)} transcripts.")
        for trans in t_in_both:
            gene = isoform_to_gene[trans]
            f.write(f"{gene}\t{trans}\n")

        f.close()
        return isoforms_file

    @staticmethod
    def check_chains_classified(chain_results_df):
        """Check whether a chain classification result is non-empty."""
        def has_more_than_one_line(file_path):
            with open(file_path, 'r') as f:
                return sum(1 for _ in islice(f, 2)) > 1

        is_complete = has_more_than_one_line(chain_results_df)
        if not is_complete:
            msg = f"Chain results file {chain_results_df} is empty! Abort."
            to_log(msg)
            raise ValueError(msg)

    @staticmethod
    def check_dependencies(toga_cls):
        """Check all dependencies."""
        # ProcessPoolExecutor is built-in to Python, no need to check for external dependencies
        # Only check C compilation status
        c_not_compiled = any(
            os.path.isfile(f) is False
            for f in [
                toga_cls.CHAIN_SCORE_FILTER,
                toga_cls.CHAIN_COORDS_CONVERT_LIB,
                toga_cls.CHAIN_FILTER_BY_ID,
                toga_cls.EXTRACT_SUBCHAIN_LIB,
                toga_cls.CHAIN_INDEX_SLIB,
            ]
        )
        if c_not_compiled:
            to_log("Warning! C code is not compiled, trying to compile...")
