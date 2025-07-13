#!/usr/bin/env python3
import argparse
import sys
from datetime import datetime as dt
from collections import defaultdict

from toga_modules.common import setup_logger
from toga_modules.common import to_log
from toga_modules.common import make_cds_track
from toga_modules.common import flatten

SOURCE = "SOURCE"
SINK = "SINK"
SCORE_THRESHOLD = 0.5
EXON_COV_THRESHOLD = 1.33
MAX_OVERLAP = 250


class Vertex:
    """Graph vertex."""

    def __init__(self, name, start, end, score):
        self.name = name
        self.start = start
        self.end = end
        self.score = score
        self.children = list()

    def add_child(self, vertex):
        if vertex not in self.children:
            self.children.append(vertex)


class Graph:
    """Build a directed graph using adjacency list."""

    def __init__(self):
        self.vertices = {}

    def add_vertex(self, vertex):
        """Add vertex if it's not in the vertices list."""
        if isinstance(vertex, Vertex) and vertex.name not in self.vertices:
            self.vertices[vertex.name] = vertex
            return True
        else:
            return False

    def add_edge(self, parent, child):
        """Find vertices with parent and child names."""
        if parent in self.vertices and child in self.vertices:
            self.vertices[parent].add_child(child)
            return True
        else:
            return False

    def topological_sort_util(self, vertex, visited, stack):
        """Perform Depth First Search

        Mark the current node as visited.
        """
        visited[vertex] = True
        # check all children of this vertex if they're visited
        for i in self.vertices[vertex].children:
            if visited[i] is False:
                self.topological_sort_util(i, visited, stack)
        # add current vertex to stack
        stack.insert(0, vertex)

    def topological_sort(self):
        """Perform topological sort.

        Use recursive function topological_sort_util().
        Mark all the vertices as not visited.
        """
        visited = {vert: False for vert in self.vertices}
        # initiate stack to store sorted vertices
        stack = []
        for vertex in self.vertices:
            if not visited[vertex]:
                self.topological_sort_util(vertex, visited, stack)
        # return a sorted list of vertices
        return stack

    def __repr__(self):
        lines = []
        for elem in self.vertices.keys():
            line = f"{elem}\t{self.vertices[elem].children}\n"
            lines.append(line)
        return "".join(lines)


def parse_args():
    """Parse CMD args."""
    app = argparse.ArgumentParser()
    app.add_argument("chain_file", help="Chain file")
    app.add_argument(
        "chain_scores_file", help="XGBoost output: chain orthology probabilities"
    )
    app.add_argument("bed_file", help="Bed file containing gene loci.")
    app.add_argument(
        "--only_fragmented",
        "--of",
        action="store_true",
        dest="only_fragmented",
        help="Output fragmented genes only.",
    )
    app.add_argument("--log_file", help="Log file", default=None)
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    return app.parse_args()


def read_gene_scores(score_file):
    """Read orthology_score.tsv file into a dict.
    Dict structure is:
    {transcript_id : [(chain, score), (chain2, score2), ..] }.
    """
    ret = defaultdict(list)
    f = open(score_file, "r")
    f.__next__()  # skip header
    for line in f:
        line_data = line.rstrip().split()
        transcript = line_data[0]
        chain_id = int(line_data[1])
        chain_score = float(line_data[2])
        if chain_score < SCORE_THRESHOLD:
            continue
        item = (chain_id, chain_score)
        ret[transcript].append(item)
    f.close()
    to_log(
        f"stitch fragments: processing {len(ret)} transcripts with scores >= {SCORE_THRESHOLD}"
    )
    return ret


def read_chain_file(chain_file):
    """Read chain file.

    Create dict chain_id: (start, end)."""
    ret = {}
    to_log(
        f"stitch fragments: parsing chain file {chain_file} to get a mapping "
        f"between chain ID and coordinates in the query genome"
    )
    f = open(chain_file, "r")
    for line in f:
        if not line.startswith("chain"):
            continue
        line_data = line.rstrip().split()
        start = int(line_data[5])
        end = int(line_data[6])
        chain_id = int(line_data[12])
        ret[chain_id] = (start, end)
    f.close()
    to_log(f"stitch fragments: parsed {len(ret)} chains")
    return ret


def read_gene_loci(bed_file):
    """For each bed entry get coding locus."""
    # TODO: not the most optimal solution, fix it
    ret = {}
    f = open(bed_file, "r")
    for line in f:
        cds_line = make_cds_track(line).split("\t")
        # extract absolute exon coordinates
        chrom_start = int(cds_line[1])
        name = cds_line[3]
        if name.endswith("_CDS"):
            name = name[:-4]
        # TODO: fix duplicated code fragment
        block_count = int(cds_line[9])
        block_sizes = [int(x) for x in cds_line[10].split(",") if x != ""]
        block_starts = [int(x) for x in cds_line[11].split(",") if x != ""]
        block_ends = [block_starts[i] + block_sizes[i] for i in range(block_count)]
        block_abs_starts = [block_starts[i] + chrom_start for i in range(block_count)]
        block_abs_ends = [block_ends[i] + chrom_start for i in range(block_count)]
        exon_nums = list(range(block_count))
        exon_coords = list(zip(exon_nums, block_abs_starts, block_abs_ends))
        ret[name] = exon_coords
    f.close()
    return ret


def build_chain_graph(chain_id_to_loc, intersecting_chains_wscores):
    """Build chain graph.

    Read chains and corresponding scores into
    a chain dictionary {chain: (start, end, score)}.
    """
    chain_graph = Graph()
    # add all vertices to the chain graph
    for chain_id, score in intersecting_chains_wscores:
        start, end = chain_id_to_loc.get(chain_id, (None, None))
        if start is None:
            raise ValueError(f"Cannot find chain {chain_id}")
        vertex = Vertex(chain_id, start, end, -1 * score)
        chain_graph.add_vertex(vertex)

    # add edges to the chain graph
    for i in chain_graph.vertices:
        for j in chain_graph.vertices:
            if i == j:
                # no need to connect the point to itself
                continue
            i_vertex = chain_graph.vertices[i]
            j_vertex = chain_graph.vertices[j]
            # allow some overlap between chains
            # defined in the MAX_OVERLAP constant
            i_vertex_and_flank = i_vertex.end - MAX_OVERLAP
            if i_vertex_and_flank <= j_vertex.start:
                chain_graph.add_edge(i_vertex.name, j_vertex.name)
    return chain_graph


def add_source_sink_graph(graph_name):
    """Add artificial Source and Sink vertices to the chain graph.

    Assign them zero length and zero score.
    """
    source_end = min(
        [graph_name.vertices[vertex].start for vertex in graph_name.vertices]
    )
    source_start = source_end
    sink_start = max(
        [graph_name.vertices[vertex].end for vertex in graph_name.vertices]
    )
    sink_end = sink_start
    graph_name.add_vertex(Vertex(SOURCE, source_start, source_end, 0))
    graph_name.add_vertex(Vertex(SINK, sink_start, sink_end, 0))

    # add edges from Source to each vertex
    for vertex in graph_name.vertices:
        if vertex != SOURCE:
            graph_name.add_edge(SOURCE, vertex)

    # add edges from each vertex to Sink
    for vertex in graph_name.vertices:
        if vertex != SINK:
            graph_name.add_edge(vertex, SINK)
    return  # all


def find_shortest_path(graph_name, source, sink, sorted_vertices):
    """Find the shortest path in directed acyclic graph.

    Initiate dictionary with the shortest paths to each node:
    {vertex: (value, path itself)}.
    """
    shortest_paths = {}
    for sorted_vertex in sorted_vertices:
        shortest_paths[sorted_vertex] = (0, [source])

    # check each child of the current vertex
    # and update shortest path to this vertex in the dictionary
    for sorted_vertex in sorted_vertices:
        children = graph_name.vertices[sorted_vertex].children
        for child in children:
            current_score = shortest_paths[child][0]
            sp_sv_0 = shortest_paths[sorted_vertex][0]
            gn_sv_s = graph_name.vertices[child].score
            score_if_updated = sp_sv_0 + gn_sv_s
            if score_if_updated < current_score:
                new_path = list(shortest_paths[sorted_vertex][1])
                if sorted_vertex not in new_path:
                    new_path.append(sorted_vertex)
                shortest_paths[child] = (score_if_updated, new_path)
    return shortest_paths[sink]


def check_exon_coverage(chains, chain_id_to_loc, exons_loci):
    """For each chain check whether it intersects all gene exons."""
    exon_num = len(exons_loci)
    chain_id_coverage = {}
    for chain_id in chains:
        # for each chain create a bool list indicating what exon
        # it intersects
        # remove exons where end < chain_start
        # and exon start > chain_end
        chain_loc = chain_id_to_loc[chain_id]
        intersect_exon_nums = [
            x[0] for x in exons_loci if x[2] > chain_loc[0] and x[1] < chain_loc[1]
        ]
        bool__exon_cov = [False for _ in range(exon_num)]
        for i in intersect_exon_nums:
            bool__exon_cov[i] = True
        # bool__exon_cov = [intersect(e, chain_loc) > 0 for e in exons_loci]
        chain_id_coverage[chain_id] = bool__exon_cov
    return chain_id_coverage


def get_average_exon_cov(chain_to_exon_cov, exon_num):
    """Compute average exon coverage."""
    exon_cov = [0 for _ in range(exon_num)]
    for coverage in chain_to_exon_cov.values():
        # there are bool values in coverage
        # covert them to int such as True = 1 and False = 0
        coverage_numeric = [1 if c else 0 for c in coverage]
        for i in range(exon_num):
            exon_cov[i] += coverage_numeric[i]
    average_cov = sum(exon_cov) / exon_num
    return average_cov


def stitch_scaffolds(chain_file, chain_scores_file, bed_file, fragments_only=False):
    """Stitch chains of fragmented orthologs."""
    to_log("stitch_fragments: started stitching fragmented orthologous loci (if any)")
    transcript_score_dict = read_gene_scores(chain_scores_file)
    # func read_chain_file returns data about all chains in the file
    # however, we need only orthologous ones
    # to avoid contamination with paralogous chains we further filter the
    # chain_id_to_loc dictionary
    # transcript score dict: gene_id: [(chain, score), (chain, score), ...]
    # Iterate over dict values (lists of tuples), get the 1st elem of each tuple (chain_id)
    orthologous_chains = set(
        flatten([val[0] for val in vals] for vals in transcript_score_dict.values())
    )
    to_log(f"stitch fragments: processing total of {len(orthologous_chains)} chains with scores")
    chain_id_to_loc__no_filt = read_chain_file(chain_file)
    chain_id_to_loc = {
        key: val for key, val in chain_id_to_loc__no_filt.items() if key in orthologous_chains
    }
    genes_to_exon_coords = read_gene_loci(bed_file)
    transcript_to_path = {}
    task_size = len(transcript_score_dict.keys())
    to_log(f"stitch fragments: processing {task_size} transcripts")
    count = 1

    for transcript, intersecting_chains_wscores in transcript_score_dict.items():
        count += 1
        if len(intersecting_chains_wscores) <= 1:
            continue
        # intersecting chains: list of tuples
        # [(chain, score), (chain, score), ...]
        # chains that intersect this transcript
        exon_coords = genes_to_exon_coords.get(transcript)
        if exon_coords is None:
            # must never happen
            err_msg = f"Cannot find a bed track for {transcript}"
            to_log(f"stitch fragments: FATAL ERROR {err_msg}")
            raise ValueError(err_msg)
        # extract some extra information about exon coverage
        intersecting_chains = [x[0] for x in intersecting_chains_wscores]
        chain_id_to_exon_cov = check_exon_coverage(
            intersecting_chains, chain_id_to_loc, exon_coords
        )
        chain_id_covers_all = {
            key: all(val for val in vals) for key, vals in chain_id_to_exon_cov.items()
        }
        if any(chain_id_covers_all.values()):
            # if there is a chain that covers the transcript entirely: skip this
            continue
        average_exon_coverage = get_average_exon_cov(
            chain_id_to_exon_cov, len(exon_coords)
        )
        if average_exon_coverage > EXON_COV_THRESHOLD:
            # skip if each exon is covered > EXON_COV_THRESHOLD times in average
            continue
        # Initiate chain graph
        chain_graph = build_chain_graph(chain_id_to_loc, intersecting_chains_wscores)
        add_source_sink_graph(chain_graph)

        # Topologically sort chain graph
        sorted_vertices = chain_graph.topological_sort()

        # Find 'longest' (=highest scoring) path in the graph =
        # find the shortest path in the graph with negative scoring vertices.
        longest_path_chain_graph = find_shortest_path(
            chain_graph, SOURCE, SINK, sorted_vertices
        )
        _, _path = longest_path_chain_graph
        path = _path[1:]  # starts with [SOURCE, ... ]
        if fragments_only and len(path) < 2:
            # this transcript is covered entirely by a single chain
            continue
        to_log(f"stitch fragments: transcript {transcript} is potentially fragmented")
        transcript_to_path[transcript] = path
        del chain_graph
    to_log(f"stitch fragments: identified {len(transcript_to_path)} fragmented transcripts")
    return transcript_to_path


if __name__ == "__main__":
    t0 = dt.now()
    args = parse_args()
    setup_logger(args.log_file)
    transcript_to_path = stitch_scaffolds(
        args.chain_file,
        args.chain_scores_file,
        args.bed_file,
        fragments_only=args.only_fragmented,
    )
    # save output
    for k, v in transcript_to_path.items():
        v_str = ",".join(map(str, v))
        print(f"{k}\t{v_str}")
    elapsed = dt.now() - t0
    to_log(f"# Elapsed: {elapsed}")
