# Created by alexandra at 25/08/2023
import argparse
import multiprocessing
import time

import FeatureEngineering.CONSTANTS as CONSTANTS
from FragmentsManipulations import FragmentSizeCountMultiprocess as fscm

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Calculate fragment size counts for cfDNA WGS bam file sample')
parser.add_argument('-b', '--bam-file', dest='bam_file', help='Bam file for which the fragemnt sizes are calculated',
                    required=True)
parser.add_argument('-o', '--output-path', dest='output_path', help='Path were the output files are saved',
                    required=True)
parser.add_argument('-mp', '--mapping-quality', dest='mapping_quality', help='Min alignment mapping quality threshold',
                    default=30, type=int)
parser.add_argument('-sfs', '--start-size', dest='size_fragment_start',
                    help='What is the fragment size to start the calculation for', default=30, type=int)
parser.add_argument('-sfe', '--end-size', dest='size_fragment_end',
                    help='What is the fragment size to end the calculation for', default=700, type=int)
parser.add_argument('-ws', '--window-size', dest='window_size',
                    help='The window for fragment size calculation to tile the genome', default=1000000, type=int)
parser.add_argument('-wf', '--window-flag', dest='by_window_flag',
                    help='Flag to calculate fragment sizes on a tiled genome by a given window ', default=False,
                    type=bool)
parser.add_argument('-snp', '--mutations-path', dest='mutations_path',
                    help='Path to the csv files with the SNP mutations', default=None)
parser.add_argument('-gv', '--genome-version', dest='genome_version',
                    help='Genome version pf the alignment', default="hg19", required=False)
parser.add_argument('-t', '--threads', dest='no_threads', help='Number of threads', default="8", type=int)
parser.add_argument('-cna', '--cnas-path', dest='cnas_path',
                    help='Path to the csv files with the CNAs', default=None)
parser.add_argument('-gcCorrection', '--GC-corection', dest='gc_correction_type',
                    help='GC Correction method: None, GCParagon, Griffin, DELFI', default=None)
parser.add_argument('-srf', '--specific-regions-flag', dest='specific_regions_flag',
                    help='Flag to calculate fragment sizes on specific genomic intervals ', default=False,
                    type=bool)
parser.add_argument('-sr', '--specific-regions', dest='specific_regions',
                    help='Specific regions to calculate fragment sizes', default=None)


args = parser.parse_args()

bam_file = args.bam_file
output_path = args.output_path
mapping_quality = args.mapping_quality
size_fragment_start = args.size_fragment_start
size_fragment_end = args.size_fragment_end
window_size = args.window_size
by_window_flag = args.by_window_flag
mutations_path = args.mutations_path
genome_version = args.genome_version
threads = args.no_threads
cnas_path = args.cnas_path
gcCorrection = args.gc_correction_type
specific_regions_flag = args.specific_regions_flag
specific_regions_file = args.specific_regions

fragmentation_patterns = fscm.FragmentSizeCount(bam_file=bam_file,
                                                output_path=output_path,
                                                mapping_quality=mapping_quality,
                                                size_fragment_start=size_fragment_start,
                                                size_fragment_end=size_fragment_end,
                                                window_size=window_size,
                                                by_window_flag=by_window_flag,
                                                mutations_path=mutations_path,
                                                genome_version=genome_version,
                                                cnas_path=cnas_path,
                                                gc_correction_type=gcCorrection,
                                                specific_regions_flag=specific_regions_flag,
                                                specific_regions_file=specific_regions_file)


def task(chromosome_name):
    fragmentation_patterns.calculate_fragment_sizes(chromosome=chromosome_name)


if __name__ == "__main__":
    start_time = time.perf_counter()
    processes = []
    # CONSTANTS.chromosomes
    # Creates 10 processes then starts them
    for chromosome in CONSTANTS.chromosomes:
        p = multiprocessing.Process(target=task, args=(chromosome,))
        p.start()
        processes.append(p)

    # Joins all the processes
    for p in processes:
        p.join()

    fragmentation_patterns.merge_chromosomes()
    finish_time = time.perf_counter()

    print(f"Program finished in {finish_time - start_time} seconds")