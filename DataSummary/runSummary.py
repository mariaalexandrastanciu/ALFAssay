# Created by alexandra at 18/08/2023
import summaryPerSample as sS
from scipy.stats import ks_2samp
import DataWrapper as DW
from helpers import CONSTANTS as C

#create data
genome_version = "hg38"
analysis = "NN_1MB_summary"  #FragSize_hg38_30_700_mp60
v_start_fragment = 30
v_end_fragment = 700
task = C.ctDNADetected
gather_samples = True

DataWrapperNMF = DW.DataWrapper(
        cancer_paths=["/Users/alexandra/PhD/NeoRheaStudy", "/Users/alexandra/PhD/PearlStudy", "/Users/alexandra/PhD/SynergyStudy"], # "/Users/alexandra/PhD/PearlStudy"
        healthy_path="/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/" + analysis,
        output_path="/Users/alexandra/PhD/FragmentationPatterns/Data/" + analysis,
        gather_samples=gather_samples, ichorCNAThreshold=0.1, start_fragment_count_size=30, start_fragment=v_start_fragment,
        end_fragment=v_end_fragment, genome_version=genome_version, task=task)
DataWrapperNMF.gather_data_for_analysis()

# create summary #########
input_files = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_1MB_summary/"
output_summary_dir = "/Users/alexandra/PhD/FragmentationPatterns/Summary/NN_1MB_summary/"
sS.create_summary(input_files, output_summary_dir)

## plot fragment density ###
output_density_plot_dir = "/Users/alexandra/PhD/FragmentationPatterns/Summary/NN_5MB_summary/Plots/FragmentDensityPlots/"
input_files = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_1MB_summary/*"
sS.plot_fragment_density(output_density_plot_dir, input_files)


# plot correlation ##
root = "/Users/alexandra/PhD/FragmentationPatterns/Summary/NN_1MB_summary"
sS.plot_correlations(root)


