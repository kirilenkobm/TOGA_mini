#!/usr/bin/env python3
"""Direct parallel execution of chain feature extraction without joblist files."""

import os
import concurrent.futures
import multiprocessing
import random
from datetime import datetime as dt

from toga_modules.chain_bed_intersect import chain_bed_intersect
from toga_modules.chain_bst_index import chain_bst_index
from toga_modules.common import load_chain_dict, setup_logger, to_log

# Global variables for worker processes
_chain_file = None
_bed_index = None  
_chain_dict = None


def _worker_init(chain_file, bed_index):
    """Initialize each worker process with shared data."""
    global _chain_file, _bed_index, _chain_dict
    
    # Set environment for this process
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    
    # Store file paths
    _chain_file = chain_file
    _bed_index = bed_index
    
    # Load chain dictionary once per worker
    index_file = chain_file.replace(".chain", ".chain_ID_position")
    if os.path.isfile(index_file):
        from toga_modules.common import load_chain_dict
        _chain_dict = load_chain_dict(index_file)
    else:
        _chain_dict = {}


def _process_single_chain(args):
    """Process a single chain-transcripts pair in worker process."""
    chain_id, transcripts = args
    
    try:
        # Import here since it's done once per worker startup
        from toga_modules.chain_runner import chain_feat_extractor
        
        # Use global variables set by worker initialization
        global _chain_file, _bed_index, _chain_dict
        
        # Call the chain feature extractor
        result = chain_feat_extractor(
            chain_id,
            transcripts,
            _chain_file,
            _bed_index,
            _chain_dict,
            extended=False
        )
        
        return (chain_id, result)
        
    except Exception as e:
        # Can't use to_log here as it might not be available in worker
        print(f"Error processing chain {chain_id}: {e}")
        return (chain_id, None)


def run_parallel_chain_extraction(
    chain_file,
    bed_file,
    bed_index,
    results_dir,
    max_workers=None,
    log_file=None,
    quiet=False
):
    """Run chain feature extraction in parallel without joblist files.
    
    Args:
        chain_file: Chain file for local alignments
        bed_file: Bed file, gene annotations  
        bed_index: Indexed bed file
        results_dir: Directory to save results
        max_workers: Maximum number of workers for parallel processing
        log_file: Path to logfile
        quiet: Don't print to console
    
    Returns:
        List of result files created
    """
    t0 = dt.now()
    setup_logger(log_file, write_to_console=not quiet)
    
    to_log("Starting direct parallel chain feature extraction")
    
    # Create results directory
    os.makedirs(results_dir, exist_ok=True)
    
    # Get chain index file
    index_file = chain_file.replace(".chain", ".chain_ID_position")
    if not os.path.isfile(index_file):
        to_log(f"Creating chain index file {index_file}")
        try:
            chain_bst_index(chain_file, index_file)
        except Exception as e:
            to_log(f"Warning: Failed to create index file: {e}")
    
    # Get chain-gene intersections
    to_log("Finding intersections between reference transcripts and chains")
    chain_genes_raw, skipped = chain_bed_intersect(chain_file, bed_file)
    
    # Convert to the format expected by chain_feat_extractor (comma-separated with trailing comma)
    chain_genes = {k: ",".join(v) + "," for k, v in chain_genes_raw.items()}
    del chain_genes_raw
    
    to_log(f"Found {len(chain_genes)} chain-transcript pairs to process")
    to_log(f"Skipped {len(skipped)} transcripts that don't intersect any chain")
    
    # Randomize order for better load balancing
    chain_items = list(chain_genes.items())
    random.shuffle(chain_items)
    
    # Set up parallel processing
    if max_workers is None:
        max_workers = min(multiprocessing.cpu_count(), len(chain_items))
    to_log(f"Using {max_workers} workers for parallel processing")
    
    # Process all chains in parallel
    all_results = []
    result_files = []
    
    try:
        # Use ProcessPoolExecutor with worker initialization
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_worker_init,
            initargs=(chain_file, bed_index)
        ) as executor:
            # Submit all jobs - now just passing (chain_id, transcripts) tuples
            futures = [
                executor.submit(_process_single_chain, (chain_id, transcripts))
                for chain_id, transcripts in chain_items
            ]
            
            # Collect results as they complete
            completed = 0
            for future in concurrent.futures.as_completed(futures):
                try:
                    chain_id, result = future.result()
                    if result is not None:
                        all_results.append(result)
                    completed += 1
                    
                    if completed % 5000 == 0:  # Log progress every 50 jobs
                        to_log(f"Completed {completed}/{len(chain_items)} chains")
                        
                except Exception as exc:
                    to_log(f"Chain processing failed with exception: {exc}")
            
            to_log(f"Completed all {len(chain_items)} chains")
    
    except KeyboardInterrupt:
        to_log("Parallel processing interrupted by user")
        raise
    
    # Write all results to a single output file
    output_file = os.path.join(results_dir, "chain_features.txt")
    with open(output_file, 'w') as f:
        for result in all_results:
            if result:
                chain_output, genes_output, time_output = result
                f.write(chain_output)
                f.write(genes_output) 
                f.write(time_output)
    
    result_files.append(output_file)
    
    total_time = dt.now() - t0
    to_log(f"Parallel chain feature extraction completed in {total_time}")
    to_log(f"Results written to {output_file}")
    
    return result_files 