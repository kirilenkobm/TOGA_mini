#!/usr/bin/env python3
"""Create an orthologous loci table.

According to predicted orthologous chains create a table with
transcript ID, chain ID, and query region coordinates.
"""
import os
from collections import defaultdict
from datetime import datetime as dt
import ctypes
import concurrent.futures
import multiprocessing

from toga_modules.common import chain_extract_id
from toga_modules.common import load_chain_dict
from toga_modules.common import make_cds_track
from toga_modules.common import die
from toga_modules.common import setup_logger
from toga_modules.common import to_log
from toga_modules.common import get_shared_lib_extension

__author__ = "Bogdan M. Kirilenko"

LOCATION = os.path.dirname(__file__)

REL_LENGTH_THR = 14  # for shorter RNAs
ABS_LENGTH_TRH = 1_000_000

ORTHOLOG = "ORTH"
PARALOG = "PARA"
SPAN = "SPAN"

MODULE_NAME_FOR_LOG = "create_orthologous_loci_table"

# connect shared lib; define input and output data types
chain_coords_conv_lib_path = os.path.join(
    LOCATION, "..", "util_c", "lib", f"libchain_coords_converter_slib{get_shared_lib_extension()}"
)

ch_lib = ctypes.CDLL(chain_coords_conv_lib_path)
ch_lib.chain_coords_converter.argtypes = [
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_char_p),
]
ch_lib.chain_coords_converter.restype = ctypes.POINTER(ctypes.c_char_p)

# Global variables for worker processes
_chain_file = None
_chain_dict = None
_bed_data = None
_chain_gene_field = None
_limit = None


def _worker_init(chain_file, bed_data, chain_gene_field, limit):
    """Initialize each worker process with shared data."""
    global _chain_file, _chain_dict, _bed_data, _chain_gene_field, _limit
    
    # Set environment for this process
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    
    # Store shared data
    _chain_file = chain_file
    _bed_data = bed_data
    _chain_gene_field = chain_gene_field
    _limit = limit
    
    # Load chain dictionary once per worker
    index_file = chain_file.replace(".chain", ".chain_ID_position")
    if os.path.isfile(index_file):
        from toga_modules.common import load_chain_dict
        _chain_dict = load_chain_dict(index_file)
    else:
        _chain_dict = None

    # Load shared library for this worker
    from toga_modules.common import get_shared_lib_extension
    location = os.path.dirname(os.path.abspath(__file__))
    lib_path = os.path.join(
        location, "..", "util_c", "lib", f"libchain_coords_converter_slib{get_shared_lib_extension()}"
    )
    
    global _ch_lib
    _ch_lib = ctypes.CDLL(lib_path)
    _ch_lib.chain_coords_converter.argtypes = [
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_char_p),
    ]
    _ch_lib.chain_coords_converter.restype = ctypes.POINTER(ctypes.c_char_p)


def extract_chain_fast(chain_file, chain_dict, chain_id):
    """Extract chain string using chain dictionary - faster than chain_extract_id.
    
    Similar to the approach used in chain_runner.py and other modules.
    """
    chain_id_int = int(chain_id)
    if chain_id_int not in chain_dict:
        # Fall back to the original method if not found
        return None
    
    start_byte, offset = chain_dict[chain_id_int]
    with open(chain_file, "rb") as f:
        f.seek(start_byte)
        chain_data = f.read(offset).decode("utf-8")
    return chain_data


def read_orthologs(orthologs_file, only_o2o=False, annotate_paralogs=False):
    """Read the orthologs file."""
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
            to_log(f"* skipping transcript {transcript}: no chains intersecting the gene")
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

        selected = chains[selected_field].copy()
        # mark a used field
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
    batch, bed_data, bdb_chain_file, chain_gene_field, limit
):
    """Precompute a region for each chain: bed pair."""
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

    # Try to load chain dictionary for faster extraction
    chain_file = bdb_chain_file.replace(".bst", ".chain")
    index_file = chain_file.replace(".chain", ".chain_ID_position")
    chain_dict = None
    
    if os.path.isfile(index_file):
        chain_dict = load_chain_dict(index_file)  # Load once

    # read regions themselves
    gene_chain_grange = defaultdict(dict)
    chains_num, iter_num = len(chain_to_genes.keys()), 0
    task_size = len(chain_to_genes)
    to_log(f"{MODULE_NAME_FOR_LOG}: for each of {task_size} involved chains, precompute regions")

    for chain_id, genes in chain_to_genes.items():
        # extract chain itself - try fast method first, fall back to original if needed
        chain_body = None
        if chain_dict is not None:
            try:
                chain_text = extract_chain_fast(chain_file, chain_dict, chain_id)
                if chain_text is not None:
                    chain_body = chain_text.encode()
            except Exception:
                pass  # Fall back to original

        if chain_body is None:
            # Use original method as fallback
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
        # convert target -> query coordinates via a chain
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

        # Extract strand information from the first line (chain header)
        chain_header = chain_coords_conv_out[0].rstrip().split()
        # Header format: "chain tName tStrand tSize tStart tEnd qName qStrand qSize qStart qEnd"
        reference_strand = chain_header[2]  # tStrand
        query_strand = chain_header[7]      # qStrand

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
                skipped.append((transcript, chain_id, "too long query locus"))
                continue
            # for each chain-gene pair save query region length, coordinates, and strand info
            # need this for required memory estimation and for orthologous loci table
            query_region = f"{q_chrom}:{q_start}-{q_end}"
            gene_chain_grange[transcript][chain_id] = {
                "query_length": que_len,
                "search_locus": query_region,
                "reference_strand": reference_strand,
                "query_strand": query_strand
            }

        del raw_ch_conv_out  # not sure if necessary but...
        iter_num += 1  # verbosity
        if iter_num % 2500 == 0:
            to_log(f"PROCESSED {iter_num} CHAINS OUT OF {task_size}")
        # eprint(f"Chain {iter_num} / {chains_num}", end="\r")
    to_log(f"{MODULE_NAME_FOR_LOG}: precomputed regions for {len(gene_chain_grange)} transcripts")
    to_log(f"{MODULE_NAME_FOR_LOG}: skipped {len(skipped)} orthologous loci")
    return gene_chain_grange, skipped


def _process_single_chain(args):
    """Process a single chain in worker process."""
    chain_id, genes = args
    
    try:
        # Use global variables set by worker initialization
        global _chain_file, _chain_dict, _bed_data, _chain_gene_field, _limit, _ch_lib
        
        # extract chain itself - try fast method first, fall back to original if needed
        chain_body = None
        if _chain_dict is not None:
            try:
                chain_text = extract_chain_fast(_chain_file, _chain_dict, chain_id)
                if chain_text is not None:
                    chain_body = chain_text.encode()
            except Exception:
                pass  # Fall back to original

        if chain_body is None:
            # Use original method as fallback
            bdb_chain_file = _chain_file.replace(".chain", ".bst")
            chain_body = chain_extract_id(bdb_chain_file, chain_id).encode()
        
        all_gene_ranges = []
        for transcript in genes:
            # get genomic coordinates for each gene
            gene_data = _bed_data.get(transcript)
            if gene_data is None:
                continue
            grange = f"{gene_data[0]}:{gene_data[1]}-{gene_data[2]}"
            all_gene_ranges.append(grange)

        if not all_gene_ranges:
            return (chain_id, {}, [])

        # we need to get corresponding regions in the query
        # for now we have chain blocks coordinates and gene
        # regions in the reference genome
        # use chain_coords_converter shared library to
        # convert target -> query coordinates via a chain
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
        raw_ch_conv_out = _ch_lib.chain_coords_converter(
            c_chain, c_shift, c_granges_num, granges_arr
        )
        chain_coords_conv_out = []  # keep lines here
        # convert C output to python-readable type
        for i in range(granges_num + 1):
            chain_coords_conv_out.append(raw_ch_conv_out[i].decode("utf-8"))

        # Extract strand information from the first line (chain header)
        chain_header = chain_coords_conv_out[0].rstrip().split()
        # Header format: "chain tName tStrand tSize tStart tEnd qName qStrand qSize qStart qEnd"
        reference_strand = chain_header[2]  # tStrand
        query_strand = chain_header[7]      # qStrand

        gene_chain_grange = {}
        skipped = []
        
        for line in chain_coords_conv_out[1:]:
            # then parse the output
            # line contains information about transcript range in the query
            # and the corresponding locus in the reference
            line_info = line.rstrip().split()
            if len(line_info) < 3:
                continue
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
            delta_gene_times = len_delta / tar_len if tar_len > 0 else 0
            transcript = genes[num]  # shared lib returns data per gene in the same order
            field = _chain_gene_field.get((chain_id, transcript))
            # check that corresponding region in the query is not too long
            # for instance query locus is 50 times longer than the gene
            # or it's longer than 1M base and also this is a SPAN chain
            high_rel_len = delta_gene_times > REL_LENGTH_THR
            high_abs_len = len_delta > ABS_LENGTH_TRH
            long_loci_field = field == SPAN
            if (high_rel_len or high_abs_len) and long_loci_field:
                skipped.append((transcript, chain_id, "too long query locus"))
                continue
            # for each chain-gene pair save query region length, coordinates, and strand info
            # need this for required memory estimation and for orthologous loci table
            query_region = f"{q_chrom}:{q_start}-{q_end}"
            gene_chain_grange[transcript] = {
                "query_length": que_len,
                "search_locus": query_region,
                "reference_strand": reference_strand,
                "query_strand": query_strand
            }

        del raw_ch_conv_out  # not sure if necessary but...
        
        return (chain_id, gene_chain_grange, skipped)
        
    except Exception as e:
        # Can't use to_log here as it might not be available in worker
        print(f"Error processing chain {chain_id}: {e}")
        return (chain_id, {}, [(str(genes), chain_id, f"processing error: {e}")])


def precompute_regions_parallel(
    batch, bed_data, bdb_chain_file, chain_gene_field, limit, max_workers=None
):
    """Precompute regions for each chain:bed pair using parallel processing."""
    to_log(f"{MODULE_NAME_FOR_LOG}: precomputing query regions for each transcript/chain pair using parallel processing")
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

    chain_file = bdb_chain_file.replace(".bst", ".chain")
    task_size = len(chain_to_genes)
    to_log(f"{MODULE_NAME_FOR_LOG}: for each of {task_size} involved chains, precompute regions")

    # Set up parallel processing
    if max_workers is None:
        max_workers = min(multiprocessing.cpu_count(), task_size)
    to_log(f"{MODULE_NAME_FOR_LOG}: using {max_workers} workers for parallel processing")

    # Prepare work items
    chain_items = list(chain_to_genes.items())
    
    # Process all chains in parallel
    gene_chain_grange = defaultdict(dict)
    all_skipped = []
    
    try:
        # Use ProcessPoolExecutor with worker initialization
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_worker_init,
            initargs=(chain_file, bed_data, chain_gene_field, limit)
        ) as executor:
            # Submit all jobs
            futures = [
                executor.submit(_process_single_chain, (chain_id, genes))
                for chain_id, genes in chain_items
            ]
            
            # Collect results as they complete
            completed = 0
            for future in concurrent.futures.as_completed(futures):
                try:
                    chain_id, chain_gene_data, chain_skipped = future.result()
                    
                    # Merge results for this chain
                    for transcript, data in chain_gene_data.items():
                        gene_chain_grange[transcript][chain_id] = data
                    
                    # Collect skipped entries
                    all_skipped.extend(chain_skipped)
                    
                    completed += 1
                    
                    if completed % 2500 == 0:  # Log progress every 2500 chains
                        to_log(f"{MODULE_NAME_FOR_LOG}: PROCESSED {completed} CHAINS OUT OF {task_size}")
                        
                except Exception as exc:
                    to_log(f"{MODULE_NAME_FOR_LOG}: Chain processing failed with exception: {exc}")
            
            to_log(f"{MODULE_NAME_FOR_LOG}: completed all {task_size} chains")
    
    except KeyboardInterrupt:
        to_log(f"{MODULE_NAME_FOR_LOG}: parallel processing interrupted by user")
        raise

    # Combine original skipped with processing skipped
    all_skipped = skipped + all_skipped
    
    to_log(f"{MODULE_NAME_FOR_LOG}: precomputed regions for {len(gene_chain_grange)} transcripts")
    to_log(f"{MODULE_NAME_FOR_LOG}: skipped {len(all_skipped)} orthologous loci")
    return gene_chain_grange, all_skipped


def read_fragments_data(in_file):
    """Read gene: fragments file."""
    to_log(f"{MODULE_NAME_FOR_LOG}: reading transcript fragments data from {in_file}")
    ret = {}
    f = open(in_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        gene = line_data[0]
        chain_str = line_data[1]
        chains = [x for x in chain_str.split(",") if x != ""]
        ret[gene] = chains
    f.close()
    to_log(f"{MODULE_NAME_FOR_LOG}: got data for {len(ret)} transcripts "
           f"potentially fragmented in the query genome")
    return ret


def create_orthologous_loci_table(
    orthologs_file, 
    bed_file, 
    bdb_bed_file, 
    bdb_chain_file, 
    target_two_bit,
    query_two_bit,
    toga_out_dir,
    chains_limit=15,
    skipped_genes=None,
    o2o_only=False,
    annotate_paralogs=False,
    fragments_data=None,
    log_file=None,
    quiet=False,
    max_workers=None,
    use_parallel=True
):
    """Create orthologous loci table with direct arguments.
    
    Args:
        orthologs_file: Output of the chain classifier
        bed_file: BED FILE
        bdb_bed_file: BDB BED FILE
        bdb_chain_file: BDB CHAIN FILE
        target_two_bit: target 2 bit
        query_two_bit: query 2 bit
        toga_out_dir: Toga output directory
        chains_limit: Skip genes with amount of orthologs more than the limit
        skipped_genes: If a gene was skipped due to number of chains limit, save it into a file
        o2o_only: Process only the genes that have a single orthologous chain
        annotate_paralogs: Annotate paralogs instead of orthologs
        fragments_data: Gene: fragments file for fragmented genomes
        log_file: Log file
        quiet: Don't print to console
        max_workers: Maximum number of workers for parallel processing
        use_parallel: Use parallel processing for precomputing regions
    """
    t0 = dt.now()
    setup_logger(log_file, write_to_console=not quiet)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash
    
    to_log(f"{MODULE_NAME_FOR_LOG}: the arguments are:")
    to_log(f"* orthologs_file: {orthologs_file}")
    to_log(f"* bed_file: {bed_file}")
    to_log(f"* bdb_bed_file: {bdb_bed_file}")
    to_log(f"* bdb_chain_file: {bdb_chain_file}")
    to_log(f"* tDB: {target_two_bit}")
    to_log(f"* qDB: {query_two_bit}")
    to_log(f"* toga_out_dir: {toga_out_dir}")
    to_log(f"* chains_limit: {chains_limit}")
    to_log(f"* skipped_genes: {skipped_genes}")
    to_log(f"* o2o_only: {o2o_only}")
    to_log(f"* annotate_paralogs: {annotate_paralogs}")
    to_log(f"* fragments_data: {fragments_data}")
    to_log(f"* use_parallel: {use_parallel}")
    to_log(f"* max_workers: {max_workers}")

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
    if use_parallel:
        regions, skipped_2 = precompute_regions_parallel(
            batch,
            bed_data,
            bdb_chain_file,
            chain_gene_field,
            chains_limit,
            max_workers,
        )
    else:
        regions, skipped_2 = precompute_regions(
            batch,
            bed_data,
            bdb_chain_file,
            chain_gene_field,
            chains_limit,
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
