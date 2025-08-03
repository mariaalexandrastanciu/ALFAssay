# Created by alexandra at 08/08/2024
import pandas as pd
import pyranges as pr

class DataHarmonization:

    def __init__(self, df, high_qualtity_regions, gc_target):
        self.df_filename = df
        self.df = pd.read_csv(df, sep="\t")
        self.high_quality_regions = pd.read_csv(high_qualtity_regions, sep="\t")
        self.gc_target = pd.read_csv(gc_target, sep=",")

    def filter_regions(self):

        sample_genomic_regions = pr.PyRanges(self.df)
        highConfidenceGenomicRegions = pr.PyRanges(self.high_quality_regions)

        genomic_regions_hc = sample_genomic_regions.overlap(highConfidenceGenomicRegions).as_df()

        return genomic_regions_hc

    def gc_Correction(self, hc_sample_genomic_regions):
        # sample_genomic_regions = pd.read_csv(self.df, sep="\t")
        # gc_target_df = pd.read_csv(self.gc_target, sep="\t")
        hc_sample_genomic_regions["GC"] = hc_sample_genomic_regions["GC"]/100
        sample_genomic_regions = hc_sample_genomic_regions.round({"GC": 2})
        gc_counts = sample_genomic_regions.groupby(["Chromosome", "GC"]).size().reset_index(name="n")
        gc_counts = gc_counts.loc[(gc_counts['GC'] >= 0.20) & (gc_counts['GC'] <= 0.80)].\
            sort_values(by=['GC', 'Chromosome'])

        gc_target_df = self.gc_target.rename(columns={"gcmed": "target", "seqnames": "Chromosome",
                                                    "gc": "GC"})

        gc_counts = pd.merge(gc_counts, gc_target_df, how="left", on=["GC", "Chromosome"])

        gc_counts["weight"] = gc_counts["target"]/gc_counts["n"]
        sample_genomic_regions = pd.merge(sample_genomic_regions, gc_counts,
                                                              how="left", on=["GC", "Chromosome"])

        return sample_genomic_regions

    def save_file(self, clean_df):
        clean_df.to_csv(self.df_filename, sep="\t", index=False)

    def data_harmonization(self):
        genomic_regions_hc = self.filter_regions()
        sample_genomic_regions = self.gc_Correction(genomic_regions_hc)

        self.save_file(sample_genomic_regions)

        return sample_genomic_regions