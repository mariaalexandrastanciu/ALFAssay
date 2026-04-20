# Created by alexandra at 25/08/2023
import multiprocessing
import time
import os
import glob
from FragmentsManipulations import FragmentSizeCountMultiprocess as fscm
import FeatureEngineering.CONSTANTS as CONSTANTS

bam_file_hg38 = "FragmentsManipulations/data/input/*.bam"
mapping_quality = 40
size_fragment_start = 40
size_fragment_end = 700
window_size = 5000000
by_window_flag = True
gc_correction_type = "GCParagon"
specific_regions_flag = False
specific_regions_file = ""


def get_output_prefix(bam_file, output_dir):
    sample_name = os.path.basename(bam_file)
    sample_name = sample_name.replace(".markDup.GCtagged.bam", "")
    sample_name = sample_name.replace(".bam", "")
    return os.path.join(output_dir, sample_name)


def task(bam_file):
    output_dir = "FragmentsManipulations/data/output"
    output_path = get_output_prefix(bam_file, output_dir)

    print(f"Processing {bam_file}")
    print(f"Output prefix: {output_path}")

    fragmentation_patterns = fscm.FragmentSizeCount(
        bam_file=bam_file,
        output_path=output_path,
        mapping_quality=mapping_quality,
        size_fragment_start=size_fragment_start,
        size_fragment_end=size_fragment_end,
        window_size=window_size,
        by_window_flag=by_window_flag,
        genome_version="hg38",
        gc_correction_type=gc_correction_type
    )

    for chromosome in CONSTANTS.chromosomes:
        fragmentation_patterns.calculate_fragment_sizes(chromosome=chromosome)

    fragmentation_patterns.merge_chromosomes()


if __name__ == "__main__":
    start_time = time.perf_counter()

    input_dir = os.path.dirname(bam_file_hg38)
    bam_files = glob.glob(os.path.join(input_dir, "*.bam"))

    processes = []

    for bam_file in bam_files:
        p = multiprocessing.Process(target=task, args=(bam_file,))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time - start_time} seconds")