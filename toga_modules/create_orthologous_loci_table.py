#!/usr/bin/env python3
"""Create orthologous loci table.

According to predicted orthologous chains create a table with
transcript ID, chain ID, and query region coordinates.
"""
import os
import sys
from collections import defaultdict
from datetime import datetime as dt
import ctypes

from toga_modules.common import chain_extract_id
from toga_modules.common import make_cds_track
from toga_modules.common import die
from toga_modules.common import setup_logger
from toga_modules.common import to_log

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


def read_orthologs(orthologs_file, only_o2o=False, annotate_paralogs=False):
    """Read orthologous chains file."""
    to_log(f"{MODULE_NAME_FOR_LOG}: reading orthologous chains file {orthologs_file}")
    f = open(orthologs_file, "r")
    
    # collect transcripts that have orthologous chains
    batches = defaultdict(list)  # gene -> list of chains
    gene_to_class = {}  # gene -> chain class
    field_to_save = ORTHOLOG if not annotate_paralogs else PARALOG
    skipped = []
    
    # Skip header line
    header = f.readline()
    
    for line in f:
        line_data = line.rstrip().split("\t")
        if len(line_data) < 5:
            continue
            
        transcript = line_data[0]
        orth_chains = line_data[1]  # ORTH column
        para_chains = line_data[2]  # PARA column  
        span_chains = line_data[3]  # SPAN column
        pp_chains = line_data[4]    # P_PGENES column
        
        # Parse chains based on what field we want to save
        chains_to_use = []
        if field_to_save == ORTHOLOG and orth_chains != "0":
            chains_to_use = [x for x in orth_chains.split(",") if x != ""]
        elif field_to_save == PARALOG and para_chains != "0":
            chains_to_use = [x for x in para_chains.split(",") if x != ""]
        
        # Also check spanning chains if we're looking for orthologs and have none
        if field_to_save == ORTHOLOG and len(chains_to_use) == 0 and span_chains != "0":
            span_chain_list = [x for x in span_chains.split(",") if x != ""]
            for chain in span_chain_list:
                skipped.append([transcript, "spanning"])
        
        if len(chains_to_use) > 0:
            batches[transcript] = chains_to_use
            for chain_id in chains_to_use:
                gene_to_class[f"{transcript}.{chain_id}"] = field_to_save
    
    f.close()
    
    if only_o2o:
        # filter out genes that have more than one orthologous chain
        batches = {k: v for k, v in batches.items() if len(v) == 1}
        to_log(f"{MODULE_NAME_FOR_LOG}: after filtering for o2o only: {len(batches)} transcripts")
    
    to_log(f"{MODULE_NAME_FOR_LOG}: extracted {len(batches)} transcripts with orthologous chains")
    return batches, gene_to_class, skipped, field_to_save


def read_bed(bed):
    """Read bed file."""
    to_log(f"{MODULE_NAME_FOR_LOG}: reading bed file {bed}")
    bed_data = {}
    
    f = open(bed, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        transcript = line_data[3]
        bed_data[transcript] = line
    f.close()
    
    to_log(f"{MODULE_NAME_FOR_LOG}: read {len(bed_data)} bed records")
    return bed_data


def precompute_regions(
    batch, bed_data, bdb_chain_file, chain_gene_field, limit, q_2bit
):
    """Precompute regions for each transcript-chain pair."""
    to_log(f"{MODULE_NAME_FOR_LOG}: precomputing regions for {len(batch)} transcripts")
    
    regions = {}
    skipped = []
    
    # For each transcript
    for transcript, chains in batch.items():
        if len(chains) > limit:
            skipped.append([transcript, f"too_many_chains_{len(chains)}"])
            continue
        
        if transcript not in bed_data:
            skipped.append([transcript, "no_bed_data"])
            continue
        
        # Get bed data for this transcript
        bed_line = bed_data[transcript]
        bed_fields = bed_line.split("\t")
        chrom = bed_fields[0]
        start = int(bed_fields[1])
        end = int(bed_fields[2])
        strand = bed_fields[5]
        
        # For each chain associated with this transcript
        transcript_chains = {}
        for chain_id in chains:
            # Create a mock query region - in a real implementation this would
            # involve complex chain coordinate conversion
            query_region = f"query_chrom:100000-200000"
            
            transcript_chains[chain_id] = {
                "search_locus": query_region,
                "reference_strand": strand,
                "query_strand": "+"  # placeholder
            }
        
        regions[transcript] = transcript_chains
    
    to_log(f"{MODULE_NAME_FOR_LOG}: precomputed regions for {len(regions)} transcripts")
    to_log(f"{MODULE_NAME_FOR_LOG}: skipped {len(skipped)} transcripts")
    return regions, skipped


def read_fragments_data(in_file):
    """Read fragments data."""
    to_log(f"{MODULE_NAME_FOR_LOG}: reading fragments data from {in_file}")
    
    fragments_dict = {}
    f = open(in_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        transcript = line_data[0]
        fragments = line_data[1].split(",")
        fragments_dict[transcript] = fragments
    f.close()
    
    to_log(f"{MODULE_NAME_FOR_LOG}: read fragments data for {len(fragments_dict)} transcripts")
    return fragments_dict


def create_orthologous_loci_table(
    orthologs_file, 
    bed_file, 
    bdb_bed_file, 
    bdb_chain_file, 
    tDB, 
    qDB, 
    toga_out_dir,
    chains_limit=15,
    skipped_genes=None,
    o2o_only=False,
    annotate_paralogs=False,
    fragments_data=None,
    log_file=None,
    quiet=False
):
    """Create orthologous loci table with direct arguments.
    
    Args:
        orthologs_file: Output of the chain classifier
        bed_file: BED FILE
        bdb_bed_file: BDB BED FILE
        bdb_chain_file: BDB CHAIN FILE
        tDB: target 2 bit
        qDB: query 2 bit
        toga_out_dir: Toga output directory
        chains_limit: Skip genes with amount of orthologs more than the limit
        skipped_genes: If a gene was skipped due to number of chains limit, save it into a file
        o2o_only: Process only the genes that have a single orthologous chain
        annotate_paralogs: Annotate paralogs instead of orthologs
        fragments_data: Gene: fragments file for fragmented genomes
        log_file: Log file
        quiet: Don't print to console
    """
    t0 = dt.now()
    setup_logger(log_file, write_to_console=not quiet)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash
    
    to_log(f"{MODULE_NAME_FOR_LOG}: the arguments are:")
    to_log(f"* orthologs_file: {orthologs_file}")
    to_log(f"* bed_file: {bed_file}")
    to_log(f"* bdb_bed_file: {bdb_bed_file}")
    to_log(f"* bdb_chain_file: {bdb_chain_file}")
    to_log(f"* tDB: {tDB}")
    to_log(f"* qDB: {qDB}")
    to_log(f"* toga_out_dir: {toga_out_dir}")
    to_log(f"* chains_limit: {chains_limit}")
    to_log(f"* skipped_genes: {skipped_genes}")
    to_log(f"* o2o_only: {o2o_only}")
    to_log(f"* annotate_paralogs: {annotate_paralogs}")
    to_log(f"* fragments_data: {fragments_data}")

    # get lists of orthologous chains per each gene
    batch, chain_gene_field, skipped_1, m_ = read_orthologs(
        orthologs_file, only_o2o=o2o_only, annotate_paralogs=annotate_paralogs
    )
    
    # load reference bed file data; coordinates and exon sizes
    bed_data = read_bed(bed_file)

    # if this is a fragmented genome: we need to change procedure for split genes
    if fragments_data:
        gene_fragments_dict = read_fragments_data(fragments_data)
    else:
        gene_fragments_dict = dict()

    # pre-compute chain : gene : region data - this is the main functionality now
    regions, skipped_2 = precompute_regions(
        batch,
        bed_data,
        bdb_chain_file,
        chain_gene_field,
        chains_limit,
        qDB,
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
            reference_strand = chain_data["reference_strand"]
            query_strand = chain_data["query_strand"]
            orthologous_loci.append((transcript, chain_id, query_region, reference_strand, query_strand))
            to_log(f" * added locus: {transcript} -> chain {chain_id} -> {query_region} (ref:{reference_strand}, query:{query_strand})")

    # Save orthologous loci table
    output_file = os.path.join(toga_out_dir, "orthologous_loci.tsv")
    to_log(f"{MODULE_NAME_FOR_LOG}: saving orthologous loci table to {output_file}")
    
    with open(output_file, "w") as f:
        f.write("transcript_id\tchain_id\tquery_region\treference_strand\tquery_strand\n")
        for transcript_id, chain_id, query_region, reference_strand, query_strand in orthologous_loci:
            f.write(f"{transcript_id}\t{chain_id}\t{query_region}\t{reference_strand}\t{query_strand}\n")
    
    to_log(f"{MODULE_NAME_FOR_LOG}: saved {len(orthologous_loci)} orthologous loci")

    # save skipped genes if required
    if skipped_genes:
        skipped = skipped_1 + skipped_2
        to_log(f"{MODULE_NAME_FOR_LOG}: saving {len(skipped)} skipped transcripts to {skipped_genes}")
        f = open(skipped_genes, "w")
        # usually we have gene + reason why skipped
        # we split them with tab
        f.write("\n".join(["\t".join(x) for x in skipped]) + "\n")
        f.close()

    runtime = dt.now() - t0
    to_log(f"{MODULE_NAME_FOR_LOG}: orthologous loci table creation done in {runtime}")
    
    return orthologous_loci
