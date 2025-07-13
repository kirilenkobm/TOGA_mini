"""Class holding all project-wide constants."""
import os

__author__ = "Bogdan M. Kirilenko"
__github__ = "https://github.com/kirilenkobm"


class Constants:
    LOCATION = os.path.dirname(__file__)

    U12_FILE_COLS = 3
    U12_AD_FIELD = {"A", "D"}
    ISOFORMS_FILE_COLS = 2
    NF_DIR_NAME = "nextflow_logs"
    NEXTFLOW = "nextflow"
    CESAR_PUSH_INTERVAL = 30  # CESAR jobs push interval
    ITER_DURATION = 60  # CESAR jobs check interval
    MEMLIM_ARG = "--memlim"
    FRAGM_ARG = "--fragments"

    CESAR_RUNNER = os.path.abspath(os.path.join(LOCATION, "cesar_runner.py"))
    CESAR_RUNNER_TMP = "{0} {1} {2} --check_loss {3} --rejected_log {4}"
    CESAR_PRECOMPUTED_REGIONS_DIRNAME = "cesar_precomputed_regions"
    CESAR_PRECOMPUTED_MEMORY_DIRNAME = "cesar_precomputed_memory"
    CESAR_PRECOMPUTED_ORTHO_LOCI_DIRNAME = "cesar_precomputed_orthologous_loci"

    CESAR_PRECOMPUTED_MEMORY_DATA = "cesar_precomputed_memory.tsv"
    CESAR_PRECOMPUTED_REGIONS_DATA = "cesar_precomputed_regions.tsv"
    CESAR_PRECOMPUTED_ORTHO_LOCI_DATA = "cesar_precomputed_orthologous_loci.tsv"

    NUM_CESAR_MEM_PRECOMP_JOBS = 500
    PARA_STRATEGIES = ["nextflow", "para", "custom"]  # TODO: add snakemake

    TEMP_CHAIN_CLASS = "temp_chain_trans_class"
    MODULES_DIR = "toga_modules"
    C_LIB_DIR = "util_c/lib"
    C_BIN_DIR = "util_c/bin"
    RUNNING = "RUNNING"
    CRASHED = "CRASHED"
    TEMP = "temp"

    # from CESAR_wrapper.py #
    FRAGMENT_CHAIN_ID = -1
    ORTH_LOC_LINE_SUFFIX = "#ORTHLOC"
    UNDEF_REGION = "None:0-0"

    # Sequence related #
    ATG_CODON = "ATG"
    XXX_CODON = "XXX"
    GAP_CODON = "---"
    NNN_CODON = "NNN"
    STOP_CODONS = {"TAG", "TGA", "TAA"}

    ACCEPTOR_SITE = ("ag",)
    DONOR_SITE = (
        "gt",
        "gc",
    )

class ConstColors:
    BLUE = "0,0,200"
    LIGHT_BLUE = "0,200,255"
    LIGHT_RED = "255,50,50"
    SALMON = "255,160,120"
    GREY = "130,130,130"
    BROWN = "159,129,112"
    BLACK = "10,10,10"


class InactMutClassesConst:
    MISS_EXON = "Missing exon"
    DEL_EXON = "Deleted exon"
    DEL_MISS = {MISS_EXON, DEL_EXON}
    COMPENSATION = "COMPENSATION"
    SSM = "SSM"
    # (ag)acceptor-EXON-donor(gt)
    SSM_D = "SSMD"  # Donor, right, GT,GC
    SSM_A = "SSMA"  # Acceptor, left, AG

    START_MISSING = "START_MISSING"
    ATG = "ATG"
    FS_DEL = "FS_DEL"
    FS_INS = "FS_INS"
    BIG_DEL = "BIG_DEL"
    BIG_INS = "BIG_INS"
    STOP = "STOP"

    STOPS = {"TAG", "TAA", "TGA"}
    D_M = {"D", "M"}
    LEFT_SPLICE_CORR = ("ag",)  # acceptor
    RIGHT_SPLICE_CORR = (
        "gt",
        "gc",
    )  # donor
    LEFT_SSID = 0
    RIGHT_SSID = 1
    ACCEPTOR = 0
    DONOR = 1

    BIG_INDEL_SIZE = 50
    SAFE_EXON_DEL_SIZE = 40  # actually 39
    FIRST_LAST_DEL_SIZE = 20
    BIG_EXON_THR = BIG_INDEL_SIZE * 5


# Standalone constants #
COMPLEMENT_BASE = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "G",
    "n": "n",
}


GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "---": "-",
    "NNN": "X",
}
