# Created by alexandra at 25/08/2023
import multiprocessing
import time
from FragmentsManipulations import FragmentSizeCountMultiprocess as fscm
import helpers.CONSTANTS as CONSTANTS

# bam_file = "/Users/alexandra/PhD/PearlStudy/bam_files/NIPT-PRL-002-25-BL-Px_S4.bam"
# bam_file = "/Users/alexandra/PhD/NeoRheaStudy/bams/NIPT-NRH-002-PL2-Nx_S5.bam"
bam_file_hg38 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/bams/Genome-IJB-HP-31-xx_S19.markDup.GCtagged.bam" # "/Users/alexandra/PhD/PearlStudy/bam_files/NIPT-PRL-001-25-BL-Px_S1.markDup.GCtagged.bam" #"/Users/alexandra/PhD/NeoRheaStudy/bams/NIPT-NRH-061-PL3-Nx_S69.bam"
bam_file_hg19 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/bams/IJB-HP-10.bam"
output_path = "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/test/testdata/healthyDataMultiP"
mapping_quality = 40
size_fragment_start = 40
size_fragment_end = 700
window_size = 1000000
by_window_flag = True
gc_correction_type="DELFI"
specific_regions_flag=False
specific_regions_file= "/refdata/amp_segments.bed"


fragmentation_patterns = fscm.FragmentSizeCount(bam_file=bam_file_hg38, output_path=output_path, mapping_quality=mapping_quality,
                                                size_fragment_start=size_fragment_start, size_fragment_end=size_fragment_end,
                                                window_size=window_size, by_window_flag=by_window_flag,
                                                # mutations_path="/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/MutScanParssed",
                                                genome_version="hg38", gc_correction_type=gc_correction_type,
                                                specific_regions_flag=specific_regions_flag,
                                                specific_regions_file=specific_regions_file)


def task(chromosome_name):
    fragmentation_patterns.calculate_fragment_sizes(chromosome=chromosome_name)


if __name__ == "__main__":
    start_time = time.perf_counter()
    processes = []

    for chromosome in CONSTANTS.chromosomes:
        p = multiprocessing.Process(target=task, args=(chromosome,))
        p.start()
        # p.join()
        processes.append(p)
    # Joins all the processes
    for p in processes:
        p.join()

    fragmentation_patterns.merge_chromosomes()
    finish_time = time.perf_counter()

    print(f"Program finished in {finish_time - start_time} seconds")