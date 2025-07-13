"""TOGA common functions."""
import os
import subprocess
import sys
import ctypes
from collections import defaultdict
import logging
import h5py

from toga_modules.toga_util import TogaUtil

SLIB_NAME = f"libchain_bst_lib{TogaUtil.get_shared_lib_extension()}"
ISOFORMS_FILE_COLS = 2


def parts(lst, n=3):
    """Split an iterable into parts with size n."""
    return [lst[i : i + n] for i in iter(range(0, len(lst), n))]


def die(msg, rc=1):
    """Show msg in stderr, exit with the rc given."""
    sys.stderr.write(msg + "\n")
    sys.exit(rc)


def setup_logger(log_file, write_to_console=True):
    # Set up logging
    logger = logging.getLogger('toga')
    logger.setLevel(logging.INFO)

    # Check if logger already has handlers to prevent duplicates
    if logger.handlers:
        # Logger already configured, don't add duplicate handlers
        return

    if write_to_console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        logger.addHandler(console_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        logger.addHandler(file_handler)


def to_log(msg):
    # Get the 'toga' logger
    logger = logging.getLogger('toga')
    # Log to both file and console
    logger.info(msg)


def bed_extract_id(index_file, gene_ids):
    """Extract a bed track from a BDB file."""
    h = h5py.File(index_file, "r")
    # accept both a list of gene_ids (type=list) and a single gene_id (type=str)
    if type(gene_ids) != str:
        keys = [str(gene_id) for gene_id in gene_ids]
    else:
        keys = [
            str(gene_ids),
        ]
    bed_lines = []
    for key in keys:
        try:  # catch key error
            b_bed_line = h[key][()]
            bed_line = b_bed_line.decode("utf-8")
            bed_lines.append(bed_line)
        except KeyError:
            # a single item not found
            # crush only if ALL tracks not found
            continue
    h.close()
    if len(bed_lines) == 0:
        # if nothing found -> there must be an error
        sys.stderr.write("Error! Bed tracks not found! (bed_extract_id)")
        sys.exit(1)
    return "".join(bed_lines)


def make_cds_track(line):
    """Trim UTRs from a bed track."""
    line_data = line.rstrip().split("\t")
    if len(line_data) != 12:
        sys.exit(f"Error! Bed line:\n{line}\nis a not bed-12 formatted line!")
    # parse bed12 line according to the specification
    chrom = line_data[0]
    chrom_start = int(line_data[1])
    # chromEnd = int(line_data[2])
    name = line_data[3]  # gene_name usually
    name += "_CDS"  # mark that UTRs are trimmed
    bed_score = int(line_data[4])  # never used
    strand = line_data[5]
    thick_start = int(line_data[6])
    thick_end = int(line_data[7])
    item_rgb = line_data[8]  # never used
    block_count = int(line_data[9])
    # chrom start and end define the entire transcript location
    # this includes both UTRs and CDS
    # thick start and end limit the CDS only
    block_sizes = [int(x) for x in line_data[10].split(",") if x != ""]
    block_starts = [int(x) for x in line_data[11].split(",") if x != ""]
    block_ends = [block_starts[i] + block_sizes[i] for i in range(block_count)]
    # block starts are given in the relative coordinates -> need to convert them
    # into absolute coordinates using chrom start
    block_abs_starts = [block_starts[i] + chrom_start for i in range(block_count)]
    block_abs_ends = [block_ends[i] + chrom_start for i in range(block_count)]
    # arrays for blocks with trimmed UTRs
    block_new_starts, block_new_ends = [], []

    for block_num in range(block_count):
        # go block-by-block
        block_start = block_abs_starts[block_num]
        block_end = block_abs_ends[block_num]

        # skip the block if it is entirely UTR
        if block_end <= thick_start:
            continue
        elif block_start >= thick_end:
            continue

        # if we are here: this is not an entirely UTR exon
        # it might intersect the CDS border or to be in the CDS entirely
        # remove UTRs: block start must be >= CDS_start (thick_start)
        # the block end must be <= CDS_end (thick_end)
        block_new_start = block_start if block_start >= thick_start else thick_start
        block_new_end = block_end if block_end <= thick_end else thick_end
        # save blocks with updated coordinates
        # also convert them back to relative coordinates with - thick_start
        # after the update thick_start/End are equal to chrom_start/End
        block_new_starts.append(block_new_start - thick_start)
        block_new_ends.append(block_new_end - thick_start)

    # block_count could change due to entirely UTR exons
    block_new_count = len(block_new_starts)
    # this is also a subject to change
    new_block_sizes = [
        block_new_ends[i] - block_new_starts[i] for i in range(block_new_count)
    ]

    # save the updated bed line with trimmed UTRs
    new_track = [
        chrom,
        thick_start,
        thick_end,
        name,
        bed_score,
        strand,
        thick_start,
        thick_end,
        item_rgb,
        block_new_count,
        ",".join([str(x) for x in new_block_sizes]) + ",",
        ",".join([str(x) for x in block_new_starts]) + ",",
    ]
    new_line = "\t".join([str(x) for x in new_track])
    return new_line


def chain_extract_id(index_file, chain_id, chain_file=None):
    """Extract chain text using an index file."""
    # within TOGA should be fine:
    chain_file = chain_file if chain_file else index_file.replace(".bst", ".chain")
    if not os.path.isfile(chain_file):
        # need this check anyway
        sys.exit(f"chain_extract_id error: cannot find {chain_file} file")
    # connect a shared library
    # .so must be there: in the modules/ dir
    script_location = os.path.dirname(__file__)
    slib_location = os.path.join(script_location, "..", "util_c", "lib", SLIB_NAME)
    sh_lib = ctypes.CDLL(slib_location)
    sh_lib.get_s_byte.argtypes = [
        ctypes.c_char_p,
        ctypes.c_uint64,
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
    ]
    sh_lib.get_s_byte.restype = ctypes.c_uint64

    # call library: find chain start byte and offset
    c_index_path = ctypes.c_char_p(index_file.encode())
    c_chain_id = ctypes.c_uint64(int(chain_id))
    c_sb = ctypes.c_uint64(0)  # write results in c_sb and c_of
    c_of = ctypes.c_uint64(0)  # provide them byref -> like pointers

    _ = sh_lib.get_s_byte(
        c_index_path, c_chain_id, ctypes.byref(c_sb), ctypes.byref(c_of)
    )

    if c_sb.value == c_of.value == 0:
        # if they are 0: nothing found then, raise Error
        sys.stderr.write(f"Error, chain {chain_id} ")
        sys.stderr.write("not found\n")
        sys.exit(1)

    # we got start byte and offset, extract chain from the file
    f = open(chain_file, "rb")
    f.seek(c_sb.value)  # jump to start_byte_position
    chain = f.read(c_of.value).decode("utf-8")  # read OFFSET bytes
    f.close()
    return chain


def flatten(lst):
    """Flat list out of a list of lists."""
    return [item for sublist in lst for item in sublist]


def load_chain_dict(chain_index_file):
    """Load dict chain ID: position in the file."""
    ans = {}
    if not os.path.isfile(chain_index_file):
        sys.exit(f"Error! File {chain_index_file} not found.")
    f = open(chain_index_file, "r")
    # tab-separated file: chain_ID, start_byte, offset
    for line in f:
        line_data = line.rstrip().split("\t")
        chain_id = int(line_data[0])
        start_byte = int(line_data[1])
        offset = int(line_data[2])
        val = (start_byte, offset)
        ans[chain_id] = val
    f.close()
    return ans


def read_isoforms_file(isoforms_file, pre_def_trans_list=None):
    """Read isoform file.

    Return gene: [isoforms] dict and
    isoforms: gene dict.
    Also returns a header (str, str)."""
    isoform_to_gene = {}
    gene_to_isoforms = defaultdict(list)
    if not os.path.isfile(isoforms_file):
        die(f"Error! Isoforms file {isoforms_file} not found")
    f = open(isoforms_file, "r")
    # the first line could be a header or may be not
    # process it separately just in case
    header_fields = f.__next__().rstrip().split("\t")
    header_gene = header_fields[0]
    header_trans = header_fields[1]
    header = (header_gene, header_trans)
    if pre_def_trans_list and header_trans not in pre_def_trans_list:
        pass
    else:
        isoform_to_gene[header_trans] = header_gene
        gene_to_isoforms[header_gene].append(header_trans)
    # process the rest of the file
    for l_num, line in enumerate(f):
        line_data = line.rstrip().split("\t")
        if len(line_data) != ISOFORMS_FILE_COLS:
            err_msg = (
                f"Isoforms file {isoforms_file} line num {l_num} corrupted: "
                f"Expected {ISOFORMS_FILE_COLS} lines, got {len(line_data)}\n"
            )
            die(err_msg)
        gene = line_data[0]
        transcript = line_data[1]
        # ore defined list of accepted transcripts: if provided, skip
        # transcripts that do not appear in this list
        if pre_def_trans_list and transcript not in pre_def_trans_list:
            continue
        isoform_to_gene[transcript] = gene
        gene_to_isoforms[gene].append(transcript)
    f.close()
    return gene_to_isoforms, isoform_to_gene, header


def make_symlink(src, dest):
    """Create a symlink.

    os.symlink expects that dest doesn't exist.
    Need to make a couple of checks before calling it.
    """
    if os.path.islink(dest):
        return
    elif os.path.isfile(dest):
        return
    os.symlink(src, dest)


def get_fst_col(path):
    """Just extract the first file's column."""
    with open(path, "r") as f:
        ret = [line.rstrip().split("\t")[0] for line in f]
    return ret


def call_process(cmd, extra_msg=None):
    """Call a subprocess and catch errors."""
    to_log(f"Calling cmd:\n{cmd}\n")
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        to_log(extra_msg) if extra_msg else None
        sys.exit(f"Error! Process:\n{cmd}\ndied! Abort.")
    to_log("Command finished with exit code 0.")
