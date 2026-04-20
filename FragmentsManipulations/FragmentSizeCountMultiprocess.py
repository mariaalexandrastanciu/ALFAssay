# Created by alexandra at 25/08/2023
import glob
import os
import time

import numpy as np
import pandas as pd
import pyranges as pr
import pysam as pysam
from FeatureEngineering import CONSTANTS
import subprocess
try:
    from Bio.SeqUtils import gc_fraction
    def GC(sequence):
        return 100 * gc_fraction(sequence, ambiguous="ignore")
except ImportError:
    # Older versions have this:
    from Bio.SeqUtils import GC

from DataHarmonization import DataHarmonization as dh


class FragmentSizeCount:
    def __init__(self, bam_file, output_path, mapping_quality=30, size_fragment_start=30, size_fragment_end=700,
                 window_size=100000, by_window_flag=True, mutations_path=None, genome_version="hg19",
                 cnas_path=None, gc_correction_type=None, specific_regions_flag=False, specific_regions_file=None):
        self.bam_file = bam_file
        self.output_path = output_path
        self.mapping_quality = mapping_quality
        self.size_fragment_start = size_fragment_start
        self.size_fragment_end = size_fragment_end
        self.bed_fragment_size = output_path  # + "_fragment_size_expanded.bed"
        self.csv_fragment_size = output_path  # + "_fragment_size_summary.csv"
        self.bed_window_fragment_size = output_path  # + "_fragment_size_window.bed"
        self.window_size = window_size
        self.by_window_flag = by_window_flag
        self.mutations_path = mutations_path
        self.genome_version = genome_version
        self.cnas_path = cnas_path
        self.gc_correction_type = gc_correction_type
        self.specific_regions_flag = specific_regions_flag
        self.specific_regions_file = specific_regions_file
        # local_dir = "/Users/alexandra/PhD/PyCharmProjects/ALFAssay_/"
        local_dir = ""
        if self.genome_version == "hg19":
            self.gaps = local_dir + "refdata/hg19_gaps"
            self.delfi_high_confidence_regions = local_dir + "refdata/HighConfidenceGenomicRegions.csv"
            self.blacklisted_regions = local_dir + "refdata/wgEncodeDukeMapabilityRegionsExcludable.bed"
            self.genome_without_gaps = None
        elif self.genome_version == "hg38":
            self.gaps = local_dir + "refdata/hg38_gaps"
            self.delfi_high_confidence_regions = local_dir + "refdata/HighConfidenceGenomicRegions_hg38_v2.csv"
            self.blacklisted_regions = local_dir + "refdata/wgEncodeDukeMapabilityRegionsExcludable_hg38.bed"
            self.genome_without_gaps = "refdata/hg38_no_gaps"
        else:
            print("Genome version undefined")
        if window_size==5000000:
            self.hg38_bins = "DataHarmonization/data/5mb_bin_hg38.csv"
        elif window_size==1000000:
            self.hg38_bins = "DataHarmonization/data/1mb_bin_hg38.csv"
        elif window_size==100000:
            self.hg38_bins = "DataHarmonization/data/100k_bin_hg38.csv"

        self.target_population = "DataHarmonization/data/target20_novaseq.csv"
        self.filename_nn = output_path

    def calculate_fragment_sizes(self, chromosome):
        """
        Calculates the fragment size for each read in the bam file given, for a specific interval
        default is fragment size between [30,700] and for a specific mapping quality(default 30);
        it creates a file per samples with all chr start position end position and the fragment size
        which is saved as a given filename
        :return:
        """

        if self.mutations_path is not None:
            sampleId = os.path.basename(self.bam_file).replace(".bam", "")
            mutations = pd.read_csv(self.mutations_path + "/" + sampleId + ".csv", sep="\t")
            withMutationsFlag = True
        else:
            mutations = pd.DataFrame(columns=["Chromosome", "Position"])
            withMutationsFlag = False

        if self.cnas_path is not None:
            sampleId = os.path.basename(self.bam_file).replace(".bam", "")
            cnas = pd.read_csv(self.cnas_path + "/" + sampleId + ".csv", sep="\t")
            cnas = cnas[cnas["CNA"] > 0]
        else:
            cnas = pd.DataFrame(columns=["Chromosome", "Start", "End", "CNA"])
        specific_regions_chr_df = []
        bam = pysam.AlignmentFile(self.bam_file, "rb")
        col_names = ['Chromosome', 'Start', 'End', "fragment_size", "NoSNPs", "isReverse", "GC", "gcParagonWeight",
                     "end4BP", "end6BP"]
        count_range = list(range(self.size_fragment_start, self.size_fragment_end + 1, 1))
        if self.specific_regions_flag:
            specific_regions_df = pd.read_csv(self.specific_regions_file, sep="\t")
            specific_regions_chr_df = specific_regions_df[specific_regions_df["Chromosome"]==chromosome]
        # for chromosome in CONSTANTS.chromosomes:
        task_per_chr(bam, chromosome, mutations, cnas, col_names, count_range,
                     self.csv_fragment_size, self.size_fragment_start, self.size_fragment_end,
                     self.gaps, self.blacklisted_regions, self.delfi_high_confidence_regions,
                     self.by_window_flag, self.mapping_quality, self.bed_fragment_size,
                     self.bed_window_fragment_size, self.window_size,
                     self.genome_without_gaps, withMutationsFlag, self.gc_correction_type,
                     self.specific_regions_flag, specific_regions_chr_df,
                     self.target_population, self.hg38_bins, self.filename_nn)

    def merge_chromosomes(self):
        chr_files_summary = glob.glob(self.csv_fragment_size + "*chr*_fragment_size_summary.csv")
        all_chr = pd.DataFrame()
        for i, chr_file in enumerate(chr_files_summary):
            chr_file = pd.read_csv(chr_file, sep="\t")
            if i == 0:
                all_chr = chr_file
            else:
                all_chr = pd.concat([all_chr, chr_file]).groupby(['fragment_size']).sum().reset_index()
        filename_all = self.csv_fragment_size + "_fragment_size_summary.csv"
        all_chr.to_csv(filename_all, sep="\t", index=False)

        if self.by_window_flag:
            input_files = self.bed_window_fragment_size + "*chr*_fragment_size_window.bed"
            output_file = self.csv_fragment_size + "_fragment_size_summary_window.bed"
            subprocess.check_call(["./bashScripts/concatFragmentSizesOnWindows.sh", input_files, output_file])
            print("File created: " + output_file)

        input_files_exp = self.bed_window_fragment_size + "*chr*_fragment_size_expanded.bed"
        output_file_exp = self.csv_fragment_size + "_fragment_size_expanded.bed"
        subprocess.check_call(["./bashScripts/concatFragmentSizesExpaned.sh", input_files_exp, output_file_exp])
        print("File created: " + output_file_exp)

        input_files_nn = self.filename_nn + "*chr*_nn_data.bed"
        output_file_nn = self.filename_nn + "_nn_data.bed"
        subprocess.check_call(["./bashScripts/concatFragmentSizesOnWindows.sh", input_files_nn, output_file_nn])
        print("File created: " + output_file_nn)


        ### for specific regions ##
        if self.specific_regions_flag:
            input_files_sr = self.bed_window_fragment_size + "*chr*_fragment_size_window_sr.bed"
            output_file_sr = self.csv_fragment_size + "_fragment_size_summary_window_sr.bed"
            subprocess.check_call(["./bashScripts/concatFragmentSizesOnWindows.sh", input_files_sr, output_file_sr])
            print("File created: " + output_file_sr)

            input_files_exp_sr = self.filename_nn + "*chr*_fragment_size_expanded_sr.bed"
            output_file_exp_sr = self.filename_nn + "_fragment_size_expanded_sr.bed"
            subprocess.check_call(["./bashScripts/concatFragmentSizesOnWindows.sh", input_files_exp_sr,
                                   output_file_exp_sr])
            print("File created: " + output_file_exp_sr)

            input_files_nn_sr = self.filename_nn + "*chr*_nn_data_sr.bed"
            output_file_nn_sr = self.filename_nn + "_nn_data_sr.bed"
            subprocess.check_call(["./bashScripts/concatFragmentSizesOnWindows.sh", input_files_nn_sr, output_file_nn_sr])
            print("File created: " + output_file_nn_sr)

def create_tiled_genome_with_arm(gaps, window_size=100000):
    gap_regions = pd.read_table(gaps).rename(columns={"chrom": "Chromosome",
                                                      "chromStart": "Start",
                                                      "chromEnd": "End"})
    gap_regions = gap_regions[gap_regions['type'].isin(["centromere", "telomere"])]

    gap_regions_genomic_regions = pr.PyRanges(gap_regions)
    chrom_sizes = pr.data.chromsizes()
    chromosomes_outside_gaps = chrom_sizes.subtract(gap_regions_genomic_regions).as_df()
    chromosomes_outside_gaps = chromosomes_outside_gaps[chromosomes_outside_gaps["Chromosome"].isin(
        CONSTANTS.chromosomes)]
    chromosomes_outside_gaps["arm"] = CONSTANTS.chr_arms
    genome_per_arm = pr.PyRanges(chromosomes_outside_gaps)

    tiled_genome = pr.genomicfeatures.tile_genome(chrom_sizes, window_size)
    tiled_genome_wg = tiled_genome.subtract(tiled_genome.overlap(gap_regions_genomic_regions)).overlap(
        genome_per_arm)

    ## add arm values
    genome_arm_counts = genome_per_arm.count_overlaps(tiled_genome_wg).as_df()
    genome_arm_counts["chr_arm_index"] = genome_arm_counts.index

    arms = []
    for genome_arm_count in genome_arm_counts.values:
        arm = [CONSTANTS.chr_arms[genome_arm_count[5]] for x in range(genome_arm_count[4])]
        arms.append(arm)
    arms_flat = [item for sublist in arms for item in sublist]
    tiled_genome_wg.arm = arms_flat
    return tiled_genome_wg


def create_tiled_genome_with_arm_hg38(gaps, window_size=1000000):
    genome = pr.PyRanges(pd.read_csv(gaps, sep="\t"))

    tiled_genome = pr.genomicfeatures.tile_genome(genome, window_size)
    tiled_genome_wa = tiled_genome.intersect(genome)

    return tiled_genome_wa


def window_fragment_size_calculation(fragments, chromosome, count_range, window_size, delfi_high_confidence_regions,
                                     gaps, genome_without_gaps, withMutationsFlag, gc_correction_type):
    if genome_without_gaps is None:
        tiled_genome_no_filters = create_tiled_genome_with_arm(gaps, window_size)
    else:
        tiled_genome_no_filters = create_tiled_genome_with_arm_hg38(genome_without_gaps, window_size)
    tiled_genome = pr.PyRanges(pd.read_table(delfi_high_confidence_regions, sep=","))
    tiled_genome_by_window = tiled_genome.overlap(tiled_genome_no_filters)
    tiled_genome_chr = tiled_genome_by_window[tiled_genome_by_window.Chromosome == chromosome]
    fragments_data = pd.read_csv(fragments, sep="\t")
    genome_fragment_lengths = [y for x, y in fragments_data.groupby('fragment_size')]
    fragment_size_counts_window_per_chr = []
    fragment_size_mutations_window_per_chr = []
    # if there are fragment lengths not available ,put ones with value 0
    SNPs_range = ["NoSNPs_" + str(s) for s in count_range]

    if len(genome_fragment_lengths) != len(count_range):
        diff = np.setdiff1d(np.asarray(count_range),
                            np.sort(np.asarray(fragments_data["fragment_size"].unique())))
        i = 0
        for counts in count_range:
            if counts in diff:
                genome_fragment_length_counts = pd.DataFrame(np.column_stack((np.asarray([0] * len(tiled_genome_chr)).T,
                                                                              np.asarray(
                                                                                  [0] * len(tiled_genome_chr)).T)),
                                                             columns=["fragment_size_exists", "NoSNPs"])
                fragment_sizes_counts = genome_fragment_length_counts["fragment_size_exists"].T.values.tolist()
                fragment_sizes_mutations = genome_fragment_length_counts["NoSNPs"].T.values.tolist()
                # i = max(0, i - 1)
            else:
                genome_fragment_length = genome_fragment_lengths[i][
                    genome_fragment_lengths[i]["fragment_size"] == counts]
                # genome_fragment_length_pr = pr.PyRanges(genome_fragment_length)
                # genome_fragment_length_counts = tiled_genome_chr.count_overlaps(
                #     genome_fragment_length_pr).as_df()

                genome_fragment_length_pr = pr.PyRanges(genome_fragment_length)
                genome_fragment_per_window = tiled_genome_chr.join(genome_fragment_length_pr, how="left").as_df()
                genome_fragment_per_window["fragment_size_exists"] = np.where(
                    genome_fragment_per_window["fragment_size"]
                    == -1, 0, 1)
                genome_fragment_per_window["NoSNPs"] = np.where(genome_fragment_per_window["NoSNPs"]
                                                                == -1, 0, genome_fragment_per_window["NoSNPs"])
                genome_fragment_length_counts = genome_fragment_per_window.groupby(["Start", "End", "arm"]). \
                    agg(
                    {"fragment_size_exists": [("FragCounts", 'sum')], "NoSNPs": [("SNPs_sum", 'sum')]}).reset_index()

                i += 1
                [fragment_sizes_counts] = genome_fragment_length_counts["fragment_size_exists"].T.values.tolist()
                [fragment_sizes_mutations] = genome_fragment_length_counts["NoSNPs"].T.values.tolist()
            fragment_size_counts_window_per_chr.append(fragment_sizes_counts)
            fragment_size_mutations_window_per_chr.append(fragment_sizes_mutations)
    else:
        for genome_fragment_length in genome_fragment_lengths:
            genome_fragment_length_pr = pr.PyRanges(genome_fragment_length)
            genome_fragment_per_window = tiled_genome_chr.join(genome_fragment_length_pr, how="left").as_df()
            if gc_correction_type is None:
                genome_fragment_per_window["fragment_size_exists"] = np.where(
                    genome_fragment_per_window["fragment_size"] == -1, 0,  1)
            elif gc_correction_type == "DELFI":
                genome_fragment_per_window["fragment_size_exists"] = np.where(genome_fragment_per_window["fragment_size"]
                                                                          == -1, 0,
                                                                          genome_fragment_per_window["weight"])
            elif gc_correction_type == "GCParagon":
                genome_fragment_per_window["fragment_size_exists"] = np.where(genome_fragment_per_window["fragment_size"]
                                                                          == -1, 0,
                                                                          genome_fragment_per_window["gcParagonWeight"])

            genome_fragment_per_window["NoSNPs"] = np.where(genome_fragment_per_window["NoSNPs"]
                                                            == -1, 0, genome_fragment_per_window["NoSNPs"])
            genome_fragment_length_counts = genome_fragment_per_window.groupby(["Start", "End", "arm"]). \
                agg({"fragment_size_exists": [("FragCounts", 'sum')], "NoSNPs": [("SNPs_sum", 'sum')]}).reset_index()

            # genome_fragment_length_counts = tiled_genome_chr.count_overlaps(genome_fragment_length_pr).as_df()
            [fragment_sizes_counts] = genome_fragment_length_counts["fragment_size_exists"].T.values.tolist()
            [fragment_sizes_mutations] = genome_fragment_length_counts["NoSNPs"].T.values.tolist()
            fragment_size_counts_window_per_chr.append(fragment_sizes_counts)
            fragment_size_mutations_window_per_chr.append(fragment_sizes_mutations)
    print("Done for chr " + chromosome)
    fragment_size_counts_per_chr_df = pd.DataFrame(np.asarray(fragment_size_counts_window_per_chr).T,
                                                   columns=count_range)
    fragment_size_counts_per_chr_df["NoReads"] = fragment_size_counts_per_chr_df[count_range]. \
        sum(axis=1, numeric_only=True)

    if withMutationsFlag:
        fragment_size_mutations_per_chr_df = pd.DataFrame(np.asarray(fragment_size_mutations_window_per_chr).T,
                                                          columns=SNPs_range)
        fragment_size_window_per_chr_df = pd.concat([tiled_genome_chr.as_df(), fragment_size_counts_per_chr_df,
                                                     fragment_size_mutations_per_chr_df], axis=1)
    else:
        fragment_size_window_per_chr_df = pd.concat([tiled_genome_chr.as_df(), fragment_size_counts_per_chr_df],
                                                    axis=1)

    return fragment_size_window_per_chr_df


def task_per_chr(bam, chromosome, mutations, cnas, col_names, count_range, csv_fragment_size,
                 size_fragment_start, size_fragment_end, gaps, blacklisted_regions,
                 delfi_high_confidence_regions,
                 by_window_flag, mapping_quality, bed_fragment_size, bed_window_fragment_size,
                 window_size, genome_without_gaps, withMutationsFlag, gc_correction_type,
                 specific_regions_flag, specific_regions_chr_df, target_population, hg38_bins,
                 filename_nn):
    mutations_on_chr = mutations[mutations["Chromosome"] == chromosome]
    no_mutations_on_chr = len(mutations_on_chr.index)
    nr_lines = 100000
    write_mode = "w"
    fragment_size_per_100000 = []
    first_100000 = True
    pos_no = 0
    filename_chr = bed_fragment_size + "_" + chromosome + "_fragment_size_expanded.bed"
    start = time.time()
    print(" start time  " + str(start) + " for chromosome: " + chromosome)
    nr_reads = bam.count(chromosome)
    for read in bam.fetch(chromosome):
        pos_no = pos_no + 1
        mutation_at_this_pos = 0
        # keep only mapped reads with above threshold quality and consider only forward reads and keep fragments
        # between 30 and 700 - we cannot guarantee the qualities outside this interval(??)
        if read.mapping_quality >= mapping_quality and read.is_proper_pair:
            if size_fragment_start <= abs(read.template_length) <= size_fragment_end:

                start_frag = read.pos
                end_frag = read.pos - read.alen + read.template_length

                if no_mutations_on_chr > 0:
                    mutation_at_this_pos = len(mutations_on_chr[(start_frag <= mutations_on_chr["Position"]) &
                                                                (mutations_on_chr[
                                                                     "Position"] <= end_frag)].index)

                gcContent = GC(read.query_sequence)
                try:
                    gcParagonWeight = read.get_tag("GC")
                except KeyError as e:
                    gcParagonWeight = 0
                    print("No GC weight for GCParagon correction")
                end_4bp = read.get_forward_sequence()[-4:]
                end_6bp = read.get_forward_sequence()[-6:]
                fragment_size_per_100000.append([chromosome, start_frag, end_frag,
                                             abs(read.template_length), mutation_at_this_pos, read.is_reverse,
                                             gcContent, gcParagonWeight, end_4bp, end_6bp])
        if pos_no % 100000 == 0 or nr_reads == pos_no or nr_reads < 100000:
            end = time.time()
            print("\n" + str(pos_no) + " lines down in " + str(end - start) + " and chr " + chromosome)

            fragment_size_per_100000_df = pd.DataFrame(fragment_size_per_100000, columns=col_names)
            # print("\n end time  " + str(time.time()) + " for chromosome: " + chromosome)
            ## filter black listed and gap genomic regions ##
            # fragment_size_per_100000_df = filter_genome(fragment_size_per_100000_df, gaps, blacklisted_regions,
            #                                             delfi_high_confidence_regions)
            if len(fragment_size_per_100000_df.index) == 0 and first_100000:
                nr_lines = nr_lines + 100000

            if len(fragment_size_per_100000_df.index) > 0:
                # fragment_size_per_100000_df.to_csv(filename_chr, sep="\t", mode=write_mode, index=False,
                #                                    header=not os.path.exists(filename_chr))
                fragment_size_per_100000_df = \
                    fragment_size_per_100000_df[fragment_size_per_100000_df["isReverse"]==False][
                        ['Chromosome', 'Start', 'End', "fragment_size", "NoSNPs", "GC", "gcParagonWeight",
                         "end4BP", "end6BP"]]
                fragment_size_per_100000_df.to_csv(filename_chr, sep="\t", mode=write_mode, index=False,
                                                   header=first_100000)
                if pos_no > nr_lines:
                    first_100000 = False
                    write_mode = "a"

                ## calculate fragment counts summary ###

                counts_of_fragment_sizes(fragment_size_per_100000_df, first_100000, chromosome, csv_fragment_size,
                                         size_fragment_start, size_fragment_end,gc_correction_type)
            fragment_size_per_100000 = []
            ## calculate fragment counts on window ###
    dataHarmonization = dh.DataHarmonization(filename_chr,
                                             hg38_bins,
                                             target_population)

    dataHarmonization.data_harmonization()

    if by_window_flag:
        fragment_size_window_per_chr_df = window_fragment_size_calculation(filename_chr,
                                                                           chromosome, count_range,
                                                                           window_size,
                                                                           delfi_high_confidence_regions,
                                                                           gaps,
                                                                           genome_without_gaps, withMutationsFlag,
                                                                           gc_correction_type)
        filename_chr_window = bed_window_fragment_size + "_" + chromosome + "_fragment_size_window.bed"
        # fragment_size_per_chr_df = self.filter_genome(fragment_size_window_per_chr_df)
        fragment_size_window_per_chr_df.to_csv(filename_chr_window, sep="\t", index=False)

    result = create_genome_for_NN(filename_chr, chromosome, hg38_bins)
    filename_nn_data = filename_nn + "_" + chromosome + "_nn_data.bed"
    result.to_csv(filename_nn_data, sep="\t", index=False)

    if specific_regions_flag:
        filename_chr_sr = bed_fragment_size + "_" + chromosome + "_fragment_size_expanded_sr.bed"
        select_specific_regions(filename_chr, chromosome, filename_chr_sr, specific_regions_chr_df)
        fragment_size_window_per_chr_df = window_fragment_size_calculation(filename_chr_sr,
                                                                           chromosome, count_range,
                                                                           window_size,
                                                                           delfi_high_confidence_regions,
                                                                           gaps,
                                                                           genome_without_gaps, withMutationsFlag,
                                                                           gc_correction_type)
        filename_chr_window = bed_window_fragment_size + "_" + chromosome + "_fragment_size_window_sr.bed"
        # fragment_size_per_chr_df = self.filter_genome(fragment_size_window_per_chr_df)
        fragment_size_window_per_chr_df.to_csv(filename_chr_window, sep="\t", index=False)
        result = create_genome_for_NN(filename_chr_sr, chromosome, hg38_bins)
        filename_nn_sr = filename_nn + "_" + chromosome + "_nn_data_sr.bed"
        result.to_csv(filename_nn_sr, sep="\t", index=False)

    end = time.time()
    print(" Finished in  " + str(end - start) + "seconds for chromosome: " + chromosome)


def counts_of_fragment_sizes(extended_frag_size_df, iter_no, chromosome, csv_fragment_size,
                             size_fragment_start, size_fragment_end, gc_correction_type):

    filename_summary = csv_fragment_size + "_" + chromosome + "_fragment_size_summary.csv"
    if iter_no:
        fragment_size = np.asarray(range(size_fragment_start, size_fragment_end + 1, 1)).T
        counts = pd.Series((np.ones(len(fragment_size)))).T
        previous_frag_counts = pd.DataFrame({'fragment_size': fragment_size, 'fragment_counts': counts.values,
                                             'NoSNPs': 0})
    else:
        previous_frag_counts = pd.read_csv(filename_summary, sep="\t")

    if gc_correction_type=="GCParagon":
        new_fragment_counts = extended_frag_size_df.groupby(['fragment_size'])['gcParagonWeight'].sum()
    else:
        new_fragment_counts = extended_frag_size_df.groupby(['fragment_size'])['fragment_size'].count()
    new_nr_mutations = extended_frag_size_df.groupby(['fragment_size'])['NoSNPs'].sum()
    new_fragment_counts_df = pd.DataFrame({'fragment_size': new_fragment_counts.index,
                                           'fragment_counts': new_fragment_counts.values,
                                           'NoSNPs': new_nr_mutations.values})

    fragment_size_summary = previous_frag_counts.set_index('fragment_size'). \
        add(new_fragment_counts_df.set_index('fragment_size'), fill_value=0).reset_index()

    fragment_size_summary.to_csv(filename_summary, sep="\t", index=False)


def filter_genome(df_to_filter, gaps, blacklisted_regions, delfi_high_confidence_regions):
    genomic_regions_prefiltered = pr.PyRanges(df_to_filter)
    ## gap regions ###
    gap_regions = pd.read_table(gaps)[["chrom", "chromStart", "chromEnd"]]. \
        rename(columns={"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End"})
    gap_genomic_regions = pr.PyRanges(gap_regions)
    ## black listed ducks ###
    black_listed_regions = pd.read_table(blacklisted_regions, header=None)[
        [0, 1, 2]].rename(columns={0: "Chromosome", 1: "Start", 2: "End"})
    black_listed_genomic_regions = pr.PyRanges(black_listed_regions)

    highConfidenceGenomicRegions = pd.read_table(delfi_high_confidence_regions, sep=",")
    highConfidenceGenomicRegions = pr.PyRanges(highConfidenceGenomicRegions)
    ## check overlaps
    genomic_regions_filtered = genomic_regions_prefiltered.overlap(gap_genomic_regions, invert=True). \
        overlap(black_listed_genomic_regions, invert=True)
    # genomic_regions_filtered = genomic_regions_filtered.overlap(black_listed_genomic_regions, invert=True)
    genomic_regions_filtered_df = genomic_regions_filtered.overlap(highConfidenceGenomicRegions).as_df()

    return genomic_regions_filtered_df

def create_genome_for_NN(filename_chr, chromosome, hg38_bins):
    sample_data = pr.PyRanges(pd.read_csv(filename_chr, sep="\t"))
    high_confidence_regions = pr.PyRanges(pd.read_table(hg38_bins, sep="\t"))
    high_confidence_regions = high_confidence_regions[high_confidence_regions.Chromosome == chromosome]
    binned_sample_data_hcr = high_confidence_regions.intersect(sample_data).as_df()
    binned_sample_data_frag = sample_data.intersect(high_confidence_regions).as_df()
    all_data = binned_sample_data_frag.merge(binned_sample_data_hcr, on=["Chromosome", "Start", "End"])
    all_data["GC"] = all_data["GC"]
    noReadsGCParagonFrag = all_data.groupby(["bin"])['gcParagonWeight'].sum().reset_index(name="noReadsGCParagonFrag")
    noReads = all_data.groupby(["bin"])["fragment_size"].count().reset_index(name="noReads")
    gcContent = all_data.groupby(["bin"])["GC"].mean().reset_index(name="gcContent")
    ultraShortFragGCPara = all_data[all_data["fragment_size"] < 90].groupby(["bin"])[
        "gcParagonWeight"].sum().reset_index(name="ultraShortFragGCPara")
    shortFragGCPara = all_data[(all_data["fragment_size"] >= 90) & (all_data["fragment_size"] < 151)].groupby(["bin"])[
        "gcParagonWeight"].sum().reset_index(name="shortFragGCPara")
    longFragGCPara = \
        all_data[(all_data["fragment_size"] >= 151) & (all_data["fragment_size"] < 221)].groupby(["bin"])[
            "gcParagonWeight"].sum().reset_index(name="longFragGCPara")
    ultraLongFragGCPara = \
        all_data[all_data["fragment_size"] > 220].groupby(["bin"])[
            "gcParagonWeight"].sum().reset_index(name="ultraLongFragGCPara")

    ultraShortFrag = all_data[all_data["fragment_size"] < 90].groupby(["bin"])[
        "fragment_size"].count().reset_index(name="ultraShortFrag")
    shortFrag = \
        all_data[(all_data["fragment_size"] >= 90) & (all_data["fragment_size"] < 151)].groupby(["bin"])[
            "fragment_size"].count().reset_index(name="shortFrag")
    longFrag = \
        all_data[(all_data["fragment_size"] >= 151) & (all_data["fragment_size"] < 221)].groupby(["bin"])[
            "fragment_size"].count().reset_index(name="longFrag")
    ultraLongFrag = \
        all_data[all_data["fragment_size"] > 220].groupby(["bin"])[
            "fragment_size"].count().reset_index(name="ultraLongFrag")

    shortOverLong = shortFrag["shortFrag"]/longFrag["longFrag"]
    shortOverLongGCPara = shortFragGCPara["shortFragGCPara"] / longFragGCPara["longFragGCPara"]
    bin_data = high_confidence_regions[["Chromosome", "Start", "End", "bin", "arm"]].as_df()
    dfs = [bin_data, gcContent, noReads, noReadsGCParagonFrag, ultraShortFrag, ultraShortFragGCPara, shortFrag,
           shortFragGCPara, longFrag, longFragGCPara, ultraLongFrag, ultraLongFragGCPara]

    result = pd.concat(dfs).groupby('bin', as_index=False).first()
    result["shortOverLong"] = shortOverLong
    result["shortOverLongGCPara"] = shortOverLongGCPara

    return result


def select_specific_regions(filename_chr, chromosome, filename_chr_sr, specific_regions_df):
    sample_data = pr.PyRanges(pd.read_csv(filename_chr, sep="\t"))
    specific_regions_pr = pr.PyRanges(specific_regions_df)
    result = sample_data.overlap(specific_regions_pr)
    result.to_csv(filename_chr_sr, sep="\t")

