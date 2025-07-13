#!/usr/bin/env python3
"""Create orthologous loci table.

According to predicted orthologous chains create a table with
transcript ID, chain ID, and query region coordinates.
"""
import argparse
import os
import sys
from collections import defaultdict
from datetime import datetime as dt
import ctypes
from common import chain_extract_id
from common import make_cds_track
from common import die
from common import setup_logger
from common import to_log

__author__ = "Bogdan M. Kirilenko"

LOCATION = os.path.dirname(__file__)

REL_LENGTH_THR = 50
ABS_LENGTH_TRH = 1_000_000

ORTHOLOG = "ORTH"
PARALOG = "PARA"
SPAN = "SPAN"

MODULE_NAME_FOR_LOG = "create_orthologous_loci_table"

# connect shared lib; define input and output data types
chain_coords_conv_lib_path = os.path.join(
    LOCATION, "..", "util_c", "lib", "libchain_coords_converter_slib.dylib"
)

ch_lib = ctypes.CDLL(chain_coords_conv_lib_path)
ch_lib.chain_coords_converter.argtypes = [
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_char_p),
]
ch_lib.chain_coords_converter.restype = ctypes.POINTER(ctypes.c_char_p)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("orthologs_file", help="Output of the chain classifier.")
    app.add_argument("bed_file", type=str, help="BED FILE")
    app.add_argument("bdb_bed_file", type=str, help="BDB BED FILE")
    app.add_argument("bdb_chain_file", type=str, help="BDB CHAIN FILE")
    app.add_argument("tDB", type=str, help="target 2 bit")
    app.add_argument("qDB", type=str, help="query 2 bit")
    app.add_argument("toga_out_dir", type=str, help="Toga output directory")

    app.add_argument(
        "--chains_limit",
        type=int,
        default=15,
        help="Skip genes with amount of orthologs more than the limit.",
    )
    app.add_argument(
        "--skipped_genes",
        default=None,
        help="If a gene was skipped due to number of chains limit, save it into a file.",
    )
    app.add_argument(
        "--o2o_only",
        "--o2o",
        action="store_true",
        dest="o2o_only",
        help="Process only the genes that have a single orthologous chain",
    )
    app.add_argument(
        "--annotate_paralogs",
        "--ap",
        action="store_true",
        dest="annotate_paralogs",
        help="Annotate paralogs instead of orthologs.",
    )
    app.add_argument(
        "--fragments_data", help="Gene: fragments file for fragmented genomes."
    )
    app.add_argument(
        "--log_file",
        default=None,
        help="Log file"
    )
    app.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Don't print to console"
    )
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_orthologs(orthologs_file, only_o2o=False, annotate_paralogs=False):
    """Read orthologs file."""
    # convert fields param string to list
    # fields = [x.upper() for x in fields_raw.split(",") if x != ""]
    to_log(f"{MODULE_NAME_FOR_LOG}: reading orthology data...")
    to_log(f"{MODULE_NAME_FOR_LOG}: for each transcript, find chains to produce annotations")

    genes_chains = {}
    chain_gene_field = {}
    skipped = []  # genes skipped at this stage
    _no_chains_intersecting = []
    f = open(orthologs_file, "r")  # open the file
    f.__next__()  # skip header
    # first column: transcript identifier
    # then: chain class fields (like column 2 - orthologous chains, 3 - paralogous)
    for line in f:
        # parse line
        line_info = line[:-1].split("\t")
        # "0" is a placeholder meaning "no chains there"
        transcript = line_info[0]
        selected, chains = [], {}

        chains[ORTHOLOG] = [x for x in line_info[1].split(",") if x != "0"]
        chains[PARALOG] = [x for x in line_info[2].split(",") if x != "0"]
        chains[SPAN] = [x for x in line_info[3].split(",") if x != "0"]
        # Processed pseudogenes column ignored -> they are processed separately
        all_chains = chains[ORTHOLOG] + chains[PARALOG] + chains[SPAN]

        if len(all_chains) == 0:
            # no orthologous loci can be created for this gene
            # because there are no chains we could use
            to_log(f"* skipping transcript {transcript}: no chains intersect the transcript")
            skipped.append((transcript, "0", "No chains intersecting the gene"))
            _no_chains_intersecting.append(transcript)
            continue

        # user can ask to process only the genes that have a single orthologous chain
        # here we check that this is the case
        not_one2one = len(chains[ORTHOLOG]) == 0 or len(chains[ORTHOLOG]) > 1
        if only_o2o and not_one2one:  # we requested only a single orthologous chain
            to_log(
                f"* skipping transcript {transcript}: only one2one requested, "
                f"got {len(chains[ORTHOLOG])} orthologs"
            )
            skipped.append((transcript, "0", "Only one2one requested, this gene didn't pass"))
            continue

        # use orthologous chains by default,
        # if no orthologous chains -> use spanning chains (SPAN)
        # no spanning chains -> use paralogous
        if annotate_paralogs:
            selected_field = PARALOG
        elif len(chains[ORTHOLOG]) > 0:
            selected_field = ORTHOLOG
        elif len(chains[SPAN]) > 0:
            selected_field = SPAN
        else:
            selected_field = PARALOG

        to_log(f"* selected chain class to annotate transcript {transcript}: {selected_field}")

        selected = chains[selected_field].copy()
        # mark used field
        for chain in selected:
            key = (chain, transcript)
            chain_gene_field[key] = selected_field

        # write to the dict, gene to chains we will use
        genes_chains[transcript] = selected

    f.close()
    die(
        "Error! No gene:chains pairs selected! Probably --fields parameter is wrong!"
    ) if len(genes_chains) == 0 else None

    to_log(f"{MODULE_NAME_FOR_LOG}: number of transcripts to create orthologous loci: {len(genes_chains)}")
    to_log(f"{MODULE_NAME_FOR_LOG}: total number of {len(chain_gene_field)} transcript/chain pairs")
    to_log(f"{MODULE_NAME_FOR_LOG}: skipped total of {len(skipped)} transcripts")
    to_log(
        f"{MODULE_NAME_FOR_LOG}: out of them, transcripts "
        f"not intersected by chains: {len(_no_chains_intersecting)}"
    )
    return genes_chains, chain_gene_field, skipped, _no_chains_intersecting


def read_bed(bed):
    """Read bed 12 file.

    For each transcript extract genetic coordinates and exon sizes.
    """
    bed_data = {}
    to_log(f"{MODULE_NAME_FOR_LOG}: reading bed file {bed}")
    f = open(bed, "r")
    for line in f:
        cds_track = make_cds_track(line).split("\t")
        bed_info = line[:-1].split("\t")
        chrom = bed_info[0]
        chrom_start = int(bed_info[1])
        chrom_end = int(bed_info[2])
        name = bed_info[3]
        block_sizes = [int(x) for x in cds_track[10].split(",") if x != ""]
        bed_data[name] = (chrom, chrom_start, chrom_end, block_sizes)
    f.close()
    to_log(f"{MODULE_NAME_FOR_LOG}: got data for {len(bed_data)} transcripts")
    return bed_data


def precompute_regions(
    batch, bed_data, bdb_chain_file, chain_gene_field, limit, q_2bit
):
    """Precompute region for each chain: bed pair."""
    to_log(f"{MODULE_NAME_FOR_LOG}: precomputing query regions for each transcript/chain pair")
    to_log(f"{MODULE_NAME_FOR_LOG}: batch size: {len(batch)}")
    chain_to_genes, skipped = defaultdict(list), []

    # revert the dict, from gene-to-chain to chain-to-genes
    to_log(f"{MODULE_NAME_FOR_LOG}: first, invert gene-to-chains dict to chain-to-genes")
    for transcript, chains_not_sorted in batch.items():
        if len(chains_not_sorted) == 0:
            to_log(f"* skip transcript {transcript}: no chains to create annotation in the query")
            skipped.append((transcript, "no orthologous chains"))
            continue

        chains = sorted(chains_not_sorted, key=lambda x: int(x))
        chains = chains[:limit]

        if len(chains_not_sorted) > limit:
            to_log(
                f"* !!for transcript {transcript} there is "
                f"{len(chains_not_sorted)} chains to produce the annotation. "
                f"Only the first {limit} chains will be considered, the remaining "
                f"chains will be skipped."
            )
            # skip genes that have > limit orthologous chains
            chains_skipped = chains[limit:]
            skipped.append(
                (
                    transcript,
                    ",".join(chains_skipped),
                    f"number of chains ({limit} chains) limit exceeded",
                )
            )

        for chain_id in chains:
            chain_to_genes[chain_id].append(transcript)

    # read regions themselves
    gene_chain_grange = defaultdict(dict)
    chains_num, iter_num = len(chain_to_genes.keys()), 0
    task_size = len(chain_to_genes)
    to_log(f"{MODULE_NAME_FOR_LOG}: for each of {task_size} involved chains, precompute regions")

    for chain_id, genes in chain_to_genes.items():
        # extract chain itself
        chain_body = chain_extract_id(bdb_chain_file, chain_id).encode()
        all_gene_ranges = []
        for transcript in genes:
            # get genomic coordinates for each gene
            gene_data = bed_data.get(transcript)
            grange = f"{gene_data[0]}:{gene_data[1]}-{gene_data[2]}"
            all_gene_ranges.append(grange)

        # we need to get corresponding regions in the query
        # for now we have chain blocks coordinates and gene
        # regions in the reference genome
        # use chain_coords_converter shared library to
        # convert target -> query coordinates via chain
        # first need to convert to C-types
        c_chain = ctypes.c_char_p(chain_body)
        c_shift = ctypes.c_int(1)
        granges_bytes = [s.encode("utf-8") for s in all_gene_ranges]
        granges_num = len(all_gene_ranges)
        c_granges_num = ctypes.c_int(granges_num)
        granges_arr = (ctypes.c_char_p * (granges_num + 1))()
        granges_arr[:-1] = granges_bytes
        granges_arr[granges_num] = None

        # then call the function
        raw_ch_conv_out = ch_lib.chain_coords_converter(
            c_chain, c_shift, c_granges_num, granges_arr
        )
        chain_coords_conv_out = []  # keep lines here
        # convert C output to python-readable type
        for i in range(granges_num + 1):
            chain_coords_conv_out.append(raw_ch_conv_out[i].decode("utf-8"))

        for line in chain_coords_conv_out[1:]:
            # then parse the output
            # line contains information about transcript range in the query
            # and the corresponding locus in the reference
            line_info = line.rstrip().split()
            # line info is: region num, region in reference, region in query
            # one line per one gene, in the same order
            num = int(line_info[0])
            # regions format is chrom:start-end
            q_chrom = line_info[1].split(":")[0]
            q_grange = line_info[1].split(":")[1].split("-")
            q_start, q_end = int(q_grange[0]), int(q_grange[1])
            que_len = q_end - q_start
            t_grange = line_info[2].split(":")[1].split("-")
            t_start, t_end = int(t_grange[0]), int(t_grange[1])
            tar_len = t_end - t_start
            len_delta = abs(tar_len - que_len)
            delta_gene_times = len_delta / tar_len
            transcript = genes[num]  # shared lib returns data per gene in the same order
            field = chain_gene_field.get((chain_id, transcript))
            # check that corresponding region in the query is not too long
            # for instance query locus is 50 times longer than the gene
            # or it's longer than 1M base and also this is a SPAN chain
            high_rel_len = delta_gene_times > REL_LENGTH_THR
            high_abs_len = len_delta > ABS_LENGTH_TRH
            long_loci_field = field == SPAN
            if (high_rel_len or high_abs_len) and long_loci_field:
                to_log(
                    f" * !!skipping transcript {transcript} / chain "
                    f"{chain_id} orthologous locus: too long query region"
                )
                skipped.append((transcript, chain_id, "too long query locus"))
                continue
            # for each chain-gene pair save query region length and coordinates
            # need this for required memory estimation and for orthologous loci table
            query_region = f"{q_chrom}:{q_start}-{q_end}"
            gene_chain_grange[transcript][chain_id] = {
                "query_length": que_len,
                "search_locus": query_region
            }

        del raw_ch_conv_out  # not sure if necessary but...
        iter_num += 1  # verbosity
        if iter_num % 10_000 == 0:
            to_log(f"PROCESSED {iter_num} CHAINS OUT OF {task_size}")
        # eprint(f"Chain {iter_num} / {chains_num}", end="\r")
    to_log(f"{MODULE_NAME_FOR_LOG}: precomputed regions for {len(gene_chain_grange)} transcripts")
    to_log(f"{MODULE_NAME_FOR_LOG}: skipped {len(skipped)} orthologous loci")
    return gene_chain_grange, skipped


def read_fragments_data(in_file):
    """Read gene: fragments file."""
    to_log(f"{MODULE_NAME_FOR_LOG}: reading transcript fragments data from {in_file}")
    ret = {}
    f = open(in_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        gene = line_data[0]
        chain_str = line_data[1]
        # chains = [int(x) for x in line_data.split(",") if x != ""]
        # actually there are strings:
        chains = [x for x in chain_str.split(",") if x != ""]
        ret[gene] = chains
    f.close()
    to_log(f"{MODULE_NAME_FOR_LOG}: got data for {len(ret)} transcripts "
           f"potentially fragmented in the query genome")
    return ret


def main():
    """Entry point."""
    t0 = dt.now()
    args = parse_args()
    setup_logger(args.log_file, write_to_console=not args.quiet)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash
    args_dict = vars(args)
    to_log(f"{MODULE_NAME_FOR_LOG}: the arguments list is:")
    for k, v in args_dict.items():
        to_log(f"* {k}: {v}")

    # get lists of orthologous chains per each gene
    batch, chain_gene_field, skipped_1, m_ = read_orthologs(
        args.orthologs_file, only_o2o=args.o2o_only, annotate_paralogs=args.annotate_paralogs
    )
    
    # load reference bed file data; coordinates and exon sizes
    bed_data = read_bed(args.bed_file)

    # if this is a fragmented genome: we need to change procedure for split genes
    if args.fragments_data:
        gene_fragments_dict = read_fragments_data(args.fragments_data)
    else:
        gene_fragments_dict = dict()

    # pre-compute chain : gene : region data - this is the main functionality now
    regions, skipped_2 = precompute_regions(
        batch,
        bed_data,
        args.bdb_chain_file,
        chain_gene_field,
        args.chains_limit,
        args.qDB,
    )

    # Create orthologous loci table directly from precomputed regions
    orthologous_loci = []
    to_log(f"{MODULE_NAME_FOR_LOG}: creating orthologous loci table for {len(regions)} transcripts")

    for transcript, chains_data in regions.items():
        gene_fragments = gene_fragments_dict.get(transcript, False)
        
        # Handle fragmented genes - only use chains that are fragments
        if gene_fragments:
            chains_data = {
                k: v for k, v in chains_data.items() if k in gene_fragments
            }

        # For each chain, add entry to orthologous loci table
        for chain_id, chain_data in chains_data.items():
            query_region = chain_data["search_locus"]
            orthologous_loci.append((transcript, chain_id, query_region))
            to_log(f" * added locus: {transcript} -> chain {chain_id} -> {query_region}")

    # Save orthologous loci table
    output_file = os.path.join(args.toga_out_dir, "orthologous_loci.tsv")
    to_log(f"{MODULE_NAME_FOR_LOG}: saving orthologous loci table to {output_file}")
    
    with open(output_file, "w") as f:
        f.write("transcript_id\tchain_id\tquery_region\n")
        for transcript_id, chain_id, query_region in orthologous_loci:
            f.write(f"{transcript_id}\t{chain_id}\t{query_region}\n")
    
    to_log(f"{MODULE_NAME_FOR_LOG}: saved {len(orthologous_loci)} orthologous loci")

    # save skipped genes if required
    if args.skipped_genes:
        skipped = skipped_1 + skipped_2
        to_log(f"{MODULE_NAME_FOR_LOG}: saving {len(skipped)} skipped transcripts to {args.skipped_genes}")
        f = open(args.skipped_genes, "w")
        # usually we have gene + reason why skipped
        # we split them with tab
        f.write("\n".join(["\t".join(x) for x in skipped]) + "\n")
        f.close()

    runtime = dt.now() - t0
    to_log(f"{MODULE_NAME_FOR_LOG}: orthologous loci table creation done in {runtime}")


if __name__ == "__main__":
    main()
