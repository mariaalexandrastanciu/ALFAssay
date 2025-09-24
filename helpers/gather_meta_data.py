# Created by alexandra at 15/02/2023
import pyranges as pr
import pandas as pd
import glob
import os.path
import numpy as np
import helpers.utils as u


def meta_data_Neorhea():

    print("Gather data for Neorhea")

    sample_names_vs_plasma = pd.read_csv("/Users/alexandra/PhD/NeoRheaStudy/plasma_name_vs_sample_name.csv", sep=",").drop_duplicates()

    fp_samples = pd.read_csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/fp_samples.csv", header=None)
    fp_samples.columns=["SampleID_FP"]
    fp_samples["SampleID"]=fp_samples["SampleID_FP"].str[5:16].str.replace("-","_")

    fp_patients = pd.merge(sample_names_vs_plasma, fp_samples, on="SampleID")
    conditions = [
        fp_samples['SampleID'].str.contains('PL1'),
        fp_samples['SampleID'].str.contains('PL2'),
        fp_samples['SampleID'].str.contains('PL3')]
    choices = ['Pre-treatment', 'C1D28', 'Surgery']
    fp_patients["timepoint"] = np.select(conditions, choices, default='Pre-treatment')

    ### create inivata results for fragmentation patterns
    inivata_data = pd.read_csv("/Users/alexandra/PhD/NeoRheaStudy/INIVATA/NeoRHEA_results_15Dec2022.csv", sep=",") 
    inivata_data.rename(columns={"patient_id":"Patient"}, inplace=True)

    meta_data_inivata = pd.merge(fp_patients, inivata_data, how="left", on=["Patient","timepoint"])

    ### create survival data for fragmentation patterns
    survival_data = pd.read_csv("/Users/alexandra/PhD/NeoRheaStudy/INIVATA/RaDaRanalysis/clinicalData.csv", sep=",")
    survival_data["Patient"] = "NRH" + survival_data["Patient"].astype(str)
    meta_survival_data = pd.merge(fp_patients, survival_data, how="left", on=["Patient"]) 
    meta_survival_data["SampleID"]=meta_survival_data["SampleID_FP"].str.replace("_fragment_size_summary.csv","")
    meta_survival_data = meta_survival_data[["SampleID", "Patient","timepoint","status","BCFS"]]
    meta_survival_data.to_csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/Neorhea_survival.csv", sep="\t", index=False, header=True)


    ### save the required info  for inivata
    meta_data_inivata["SampleID"]=meta_data_inivata["SampleID_FP"].str.replace("_fragment_size_summary.csv","")
    meta_data_inivata = meta_data_inivata[["SampleID", "Patient","timepoint","mean_VAF","ctDNA_detected"]]
    conditions = [
        meta_data_inivata["ctDNA_detected"]=="TRUE",
        meta_data_inivata["ctDNA_detected"]=="FALSE",
        meta_data_inivata["ctDNA_detected"]=="Fail"]
    choices = [True, False, "NA"] 
    meta_data_inivata["ctDNA_detected"] = np.select(conditions, choices, default="NA")
    meta_data_inivata.to_csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/Neorhea_INIVATA_results.csv",
                             sep="\t", index=False, header=True)


    print("Done gathering data for Neorhea")


###TO do - FIX this

def meta_data_Pearl():
    print("Gather data for Pearl")

    fp_samples = pd.read_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fp_samples.csv", header=None)
    survival_data = pd.read_csv("/Users/alexandra/PhD/PearlStudy/survivalInfo.csv", sep=",")
    fp_samples.columns=["SampleID_FP"]
    fp_samples["Patient"]=fp_samples["SampleID_FP"].str[5:15]

    meta_survival_data = pd.merge(fp_samples,survival_data,how="left",on="Patient")
    meta_survival_data["SampleID"]=meta_survival_data["SampleID_FP"].str.replace("_fragment_size_summary.csv","")

    meta_survival_data = meta_survival_data[["SampleID", "status", "BPFS"]]
    meta_survival_data.to_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/Pearl_survival_bl_pd.csv", sep="\t",
                              index=False, header=True)

    pearl = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/PearlStudyMetaData.csv", sep="\t")
    pearl_cols = pearl.columns
    pearl = pearl[["PatientId", "Label", "ichorTF", "timepoint", "VAF", "ctDNADetected", "VAFg0p001"]]
    meta_survival_data["PatientId"] = meta_survival_data["SampleID"]

    pearl = pd.merge(pearl, meta_survival_data, on="PatientId")
    pearl["survivalStatus"] = pearl["status"]
    pearl["survivalTime"] = pearl["BPFS"]
    pearl = pearl[pearl_cols]
    pearl.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/PearlStudyMetaData_.csv", sep="\t",
                 index=False, header=True)

    print("Done gathering data for Pearl")




def meta_data_Pearl_new():
    print("Gather data for Pearl")
    #
    # fp_samples = pd.read_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fp_samples.csv", header=None)
    # survival_data = pd.read_csv("/Users/alexandra/PhD/PearlStudy/pearl_followup.csv", sep=",")
    # survival_data_dates = pd.read_csv("/Users/alexandra/PhD/PearlStudy/survivalInfo.csv", sep=",")
    #
    # survival_data_dates['FirstRecururenceDt'] = pd.to_datetime(survival_data_dates['FirstRecururenceDt'])
    # survival_data_dates['FollowupDate'] = pd.to_datetime(survival_data_dates['FollowupDate'])
    # survival_data_dates['DeathDate'] = pd.to_datetime(survival_data_dates['DeathDate'])
    #
    # fp_samples.columns=["SampleID"]
    # fp_samples["SampleID"] = fp_samples["SampleID"].str.replace("_fragment_size_summary.csv", "")
    # fp_samples["Patient"]=fp_samples["SampleID"].str[9:12].str.lstrip('0')
    # fp_samples["Patient"] = fp_samples["Patient"].apply(int)
    #
    #
    # survival_data["Patient"] = survival_data["id"].values
    # survival_data["Patient"] = survival_data["Patient"].apply(int)
    #
    # meta_survival_data = pd.merge(fp_samples,survival_data, how="left", on="Patient")
    # meta_survival_data["BPFS"] = meta_survival_data["PFStime"] * 12 * 30
    # meta_survival_data["status"] = meta_survival_data["PFSevent"]
    # meta_survival_data = meta_survival_data[["SampleID", "status", "BPFS"]]
    # meta_survival_data["BPFS"] = meta_survival_data["BPFS"]*12*30
    # meta_survival_data.to_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/Pearl_survival_bl_pd.csv", sep=",",
    #                           index=False, header=True)

    meta_survival_data = pd.read_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/Pearl_survival_david.csv", sep=",")

    pearl = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/PearlStudyMetaData.csv", sep="\t")
    pearl_cols = pearl.columns
    pearl = pearl[["PatientId",	"Label"	,"ichorTF"	,"timepoint",	"VAF"	,"ctDNADetected", "VAFg0p001"]]
    meta_survival_data["PatientId"] = meta_survival_data["SampleID"]

    pearl = pd.merge(pearl, meta_survival_data, on="PatientId")
    pearl["survivalStatus"] = pearl["status"]
    pearl["survivalTime"] = pearl["BPFS"]
    pearl = pearl[pearl_cols]
    pearl.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/PearlStudyMetaData_.csv", sep="\t",
                              index=False, header=True)



    print("Done gathering data for Pearl")

meta_data_Pearl_new()

def meta_data_healthy():
    print("Gather data for healthy")
    # meta_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/healthyMetaData.csv", sep="\t")
    fp_samples = pd.read_csv("/Users/alexandra/PhD/healthy_sWGS/healthy_samples.csv", sep="\t", header=None)

    # meta_data["PatientId"] = "Genome-" + meta_data["PatientId"]
    fp_samples["PatientId"] = fp_samples[0].replace({"_nn_data.bed":""}, regex=True)
    fp_samples["Label"] = fp_samples["PatientId"].apply(u.get_label_healthy)
    fp_samples["ichorTF"] = 0
    fp_samples["timepoint"] = "BL"
    fp_samples["VAF"] = 0
    fp_samples["ctDNADetected"] = 0
    fp_samples["survivalStatus"] = 0
    fp_samples["survivalTime"] = "NA"
    fp_samples["VAFg0p001"] = "NA"
    columns = ["PatientId", "Label", "ichorTF", "timepoint", "VAF", "ctDNADetected",
                                       "survivalStatus", "survivalTime", "VAFg0p001"]
    healthy = pd.DataFrame(fp_samples[columns], columns=columns)

    # fp_samples = pd.read_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fp_samples.csv", header=None)
    # survival_data = pd.read_csv("/Users/alexandra/PhD/PearlStudy/survivalInfo.csv", sep=",")
    # fp_samples.columns=["SampleID_FP"]
    # fp_samples["Patient"]=fp_samples["SampleID_FP"].str[5:15]
    #
    # meta_survival_data = pd.merge(fp_samples,survival_data,how="left",on="Patient")
    # meta_survival_data["SampleID"]=meta_survival_data["SampleID_FP"].str.replace("_fragment_size_summary.csv","")
    #
    # meta_survival_data = meta_survival_data[["SampleID", "status", "BPFS"]]
    healthy.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/healthyMetaData.csv", sep="\t",
                       index=False, header=True)

    print("Done")

# meta_data_healthy();
# meta_data_Neorhea();


# import MetaDataParsing as MTP
# # cancer_paths="/Users/alexandra/PhD/NeoRheaStudy"
# # healthy_path = "/Users/alexandra/PhD/healthy_sWGS"
# # output_path = "/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData"
# # add_meta_obj = MTP.MetaDataParsing(cancer_paths, healthy_path, output_path)
# # add_meta_obj.dummyMetaDataHealthy()
# #
# # meta_data_healthy()
#
# MetaDataParsing = MTP.MetaDataParsing(
#     cancer_paths=["/Users/alexandra/PhD/NeoRheaStudy", "/Users/alexandra/PhD/PearlStudy"],
#     healthy_path="/Users/alexandra/PhD/healthy_sWGS",
#     output_path="/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData")
# MetaDataParsing.dummyMetaDataHealthy()
# meta_data_healthy()
# MetaDataParsing.metaDataPerStudy()

def meta_data_delfi():
    metadata = pd.read_table("/Users/alexandra/PhD/DELFIStudy/EGAD00001007796-metadata/samples.tsv")
    metadata_females = metadata[metadata["biological_sex"] == "female"]
    metadata_cancer_healthy = metadata_females[metadata_females["phenotype"].isin(["lung cancer", "cancer",
                                                                                   "no lung cancer, other cancer",
                                                                                   "healthy"])]

    metadata_cancer_healthy['phenotype'] = np.where(metadata_cancer_healthy['phenotype'] == "healthy",
             "Healthy", "Cancer")
    delfi_metadata = pd.DataFrame(data=metadata_cancer_healthy[["subject_id", "phenotype"]].values,
                                  columns=["PatientId", "Label"])
    delfi_metadata["ichorTF"] = "NA"
    delfi_metadata["timepoint"] = "BL"
    delfi_metadata["VAF"] = "NA"
    delfi_metadata["ctDNADetected"] = "NA"
    delfi_metadata["survivalStatus"] = "NA"
    delfi_metadata["survivalTime"] = "NA"
    delfi_metadata.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/DELFIMetaData.csv", sep="\t",
                          index=False, header=True)
    print("ok")

# meta_data_delfi()


def first_healthy():
    print("Gather data for healthy")
    # meta_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/healthyMetaData.csv", sep="\t")
    fp_samples = pd.read_csv("/Users/alexandra/PhD/healthy_sWGS/first_batch_healthy_samples.csv", sep="\t", header=None)

    # meta_data["PatientId"] = "Genome-" + meta_data["PatientId"]
    fp_samples["PatientId"] = fp_samples[0].replace({"_fragment_size_summary_window.csv":""}, regex=True)
    # fp_samples["Label"] = "Healthy"
    # fp_samples["ichorTF"] = 0
    # fp_samples["timepoint"] = "BL"
    # fp_samples["VAF"] = 0
    # fp_samples["ctDNADetected"] = False
    # fp_samples["survivalStatus"] = 0
    # fp_samples["survivalTime"] = "NA"
    columns = ["PatientId", "Label", "ichorTF", "timepoint", "VAF", "ctDNADetected",
                                       "survivalStatus", "survivalTime"]

    healthy = pd.DataFrame(fp_samples[["PatientId"]], columns=["PatientId"])

    # fp_samples = pd.read_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fp_samples.csv", header=None)
    # survival_data = pd.read_csv("/Users/alexandra/PhD/PearlStudy/survivalInfo.csv", sep=",")
    # fp_samples.columns=["SampleID_FP"]
    # fp_samples["Patient"]=fp_samples["SampleID_FP"].str[5:15]
    #
    # meta_survival_data = pd.merge(fp_samples,survival_data,how="left",on="Patient")
    # meta_survival_data["SampleID"]=meta_survival_data["SampleID_FP"].str.replace("_fragment_size_summary.csv","")
    #
    # meta_survival_data = meta_survival_data[["SampleID", "status", "BPFS"]]
    healthy.to_csv("/Users/alexandra/PhD/healthy_sWGS/first_batch_healthy_sample_names.csv", sep="\t",
                       index=False, header=True)

    print("Done")
# first_healthy()

def add_vafg01():
    neorhea_meta_data_real = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/NeoRheaStudyMetaData.csv",
                                         sep="\t")
    # neorhea_meta_data_VAFg0p001 = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/NeoRheaStudyMetaData.csv",
    #                                           sep="\t")

    neorhea_meta_data_real["VAFg0p001"] = np.where(neorhea_meta_data_real["VAF"]>0.001, True, False)

    neorhea_meta_data_real.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/NeoRheaStudyMetaData.csv",
                                  sep="\t", index=False, header=True)

    files_to_update = glob.glob("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/*.csv")
    for file in files_to_update:
        if "NeoRheaStudy" in file:
            None
        elif "AllStudiesMetaData" in file:
            meta_data = pd.read_csv(file, sep="\t")
            VAFg0p001 = np.where(meta_data["VAF"]>0.001, True, False)
            meta_data["VAFg0p001"] = np.where(meta_data["study"]=="NeoRheaStudy", VAFg0p001,
                                              meta_data["ctDNADetected"])
            meta_data.to_csv(file, sep="\t", index=False, header=True)
        else:
            print(file)
            meta_data = pd.read_csv(file, sep="\t")
            meta_data["VAFg0p001"] = meta_data["ctDNADetected"]
            meta_data.to_csv(file, sep="\t", index=False, header=True)

# add_vafg01()


def meta_data_synergy():
    print("Gather data for synergy")
    # meta_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/healthyMetaData.csv", sep="\t")
    fp_samples = pd.read_csv("/Users/alexandra/PhD/SynergyStudy/synergy_samples.csv", sep="\t")
    # week1_ctDNA_pos = pd.read_csv("/Users/alexandra/PhD/SynergyStudy/metadata/week1_ctDNA_David.csv", sep=",")
    # meta_data["PatientId"] = "Genome-" + meta_data["PatientId"]
    ichorCNA=pd.read_csv("/Users/alexandra/PhD/SynergyStudy/FragmentationPatterns/ichor_CNA_results.csv", sep="\t")

    survival = pd.read_csv("/Users/alexandra/PhD/SynergyStudy/FragmentationPatterns/Synergy_survival.csv", sep="\t")

    fp_samples["PatientId"] = fp_samples["Sample"]
    fp_samples["Label"] = "Cancer"
    fp_samples["timepoint"] = fp_samples["Sample"].apply(u.get_timepoint_from_path)
    fp_samples["VAFg0p001"] = "NA"

    fp_samples = pd.merge(fp_samples,ichorCNA, left_on="Sample", right_on="SampleID" )
    fp_samples = fp_samples.rename(columns={"ichorCNATF": "ichorTF"})
    fp_samples["VAF"] = fp_samples["ichorTF"]
    fp_samples["ctDNADetected"] = np.where(fp_samples["VAF"]>=0.03, True, False)

    fp_samples = pd.merge(fp_samples, survival, left_on="Sample", right_on="SampleID")

    fp_samples = fp_samples.rename(columns={"status": "survivalStatus", "BPFS":"survivalTime"})

    columns = ["PatientId", "Label", "ichorTF", "timepoint", "VAF", "ctDNADetected",
                                       "survivalStatus", "survivalTime", "VAFg0p001"]
    synergy_meta = pd.DataFrame(fp_samples[columns], columns=columns)

    # fp_samples = pd.read_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fp_samples.csv", header=None)
    # survival_data = pd.read_csv("/Users/alexandra/PhD/PearlStudy/survivalInfo.csv", sep=",")
    # fp_samples.columns=["SampleID_FP"]
    # fp_samples["Patient"]=fp_samples["SampleID_FP"].str[5:15]
    #
    # meta_survival_data = pd.merge(fp_samples,survival_data,how="left",on="Patient")
    # meta_survival_data["SampleID"]=meta_survival_data["SampleID_FP"].str.replace("_fragment_size_summary.csv","")
    #
    # meta_survival_data = meta_survival_data[["SampleID", "status", "BPFS"]]
    synergy_meta.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/SynergyStudyMetaData.csv", sep="\t",
                       index=False, header=True)

    print("Done")

# samples_synergy = pd.read_csv("/Users/alexandra/PhD/SynergyStudy/synergy_samples.csv", sep="\t", header=None)
# samples_synergy[0] = samples_synergy[0].replace({"_nn_data.bed":""}, regex=True)
# samples_synergy = samples_synergy.rename(columns={0: "Sample"})
# samples_synergy.to_csv("/Users/alexandra/PhD/SynergyStudy/synergy_samples.csv", sep="\t",
#                        index=False, header=True)
# meta_data_synergy()


def gather_all_meta():
    synergy = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/SynergyStudyMetaData.csv", sep="\t")
    neorhea = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/NeoRheaStudyMetaData.csv", sep="\t")
    pearl = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/PearlStudyMetaData_.csv", sep="\t")
    healthy = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/healthyMetaData.csv", sep="\t")
    batch = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/batches.csv", sep="\t")
    healthy_common = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/patients_with_common_cleaned.csv", sep=",")
    depth_df = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/mosdepth.csv", sep=" ", header=None)
    depth = (
        depth_df
            .copy()
            .rename(columns={depth_df.columns[0]: 'raw_id', depth_df.columns[1]: 'depth'})
    )
    depth['PatientId'] = depth['raw_id'].str.replace(
        r'\.mosdepth\.summary\.txt$', '', regex=True
    )
    depth = depth[['PatientId', 'depth']]

    healthy_validation = healthy_common[healthy_common["common"]=="no"]["PatientId"]
    synergy["study"] = "Synergy"

    neorhea["study"] = "Neorhea"
    pearl["study"] = "Pearl"
    healthy["study"] = "Healthy"
    pearl = pearl[pearl["Label"] == "Cancer"]
    neorhea = neorhea[neorhea["Label"] == "Cancer"]

    median_status_synergy = np.median(synergy[["survivalTime"]].values)
    synergy["MedianStatus"] = np.where(((synergy["survivalTime"] <= median_status_synergy) & (synergy["survivalStatus"] ==1)), 1, 0)

    median_status_pearl = np.median(pearl[["survivalTime"]].values)
    pearl["MedianStatus"] = np.where(((pearl["survivalTime"] <= median_status_pearl) & (pearl["survivalStatus"] ==1)), 1, 0)

    median_status_neorhea = np.nanmedian(neorhea[["survivalTime"]].values)
    neorhea["MedianStatus"] = np.where(((neorhea["survivalTime"] <= median_status_neorhea) & (neorhea["survivalStatus"] ==1)) , 1, 0)

    healthy["MedianStatus"] = 0

    synergy["ctDNADetected"] = np.where(synergy["ctDNADetected"] == False, 0, synergy["ctDNADetected"])
    synergy["ctDNADetected"] = np.where(synergy["ctDNADetected"] == True, 1, synergy["ctDNADetected"])

    pearl["ctDNADetected"] = np.where(pearl["ctDNADetected"]==False, 0, pearl["ctDNADetected"])
    pearl["ctDNADetected"] = np.where(pearl["ctDNADetected"] == True, 1, pearl["ctDNADetected"])
    pearl = pearl[pearl["Label"] == "Cancer"]

    neorhea["ctDNADetected"] = np.where(neorhea["ctDNADetected"] == False, 0, neorhea["ctDNADetected"])
    neorhea["ctDNADetected"] = np.where(neorhea["ctDNADetected"] == True, 1, neorhea["ctDNADetected"])
    neorhea = neorhea[neorhea["Label"] == "Cancer"]

    all_data = pd.concat([neorhea, pearl, synergy, healthy], axis=0)

    all_data["Label"] = np.where(all_data["Label"]=="Cancer", 1, 0)
    all_data = pd.merge(all_data, batch, on="PatientId", how="left")
    all_data["batch"] = np.where(all_data["PatientId"]=="Genome-IJB-HP-30-xx_S18", "A00154:1106:H3L22DRX2", all_data["batch"])
    all_data["Process"] = np.where(all_data["study"] == "Pearl", "Validation", "Training")

    all_data["Process"] = np.where(all_data["PatientId"].isin(healthy_validation), "Validation", all_data["Process"])

    all_data = all_data.merge(depth, on='PatientId', how='inner')
    lowcov = ["NIPT-P0104-ArmB-SNRFR0085-W1-Screen_S37", "NIPT-P0154-ArmA-SNRFR0005-W3_S44",
              "NIPT-P0156-ArmA-SNRFR0046-W1-Screen_S46", "NIPT-P0156-ArmA-SNRFR0046-W13-EOC_S48",
              "NIPT-P0156-ArmA-SNRFR0046-W3_S47"]  #

    # metadata = metadata[(metadata["study"]!="Neorhea") ]
    all_data = all_data[~all_data["PatientId"].isin(lowcov)]
    all_data.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData_.csv", sep="\t",
                        index=False, header=True)

gather_all_meta()


def gather_all_meta_labels():
    synergy = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/SynergyStudyMetaData.csv", sep="\t")
    neorhea = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/NeoRheaStudyMetaData.csv", sep="\t")
    pearl = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/PearlStudyMetaData.csv", sep="\t")
    healthy = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/healthyMetaData.csv", sep="\t")
    batch = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/batches.csv", sep="\t")
    synergy["study"] = "Synergy"
    neorhea["study"] = "Neorhea"
    pearl["study"] = "Pearl"
    healthy["study"] = "Healthy"
    pearl = pearl[pearl["Label"] == "Cancer"]
    neorhea = neorhea[neorhea["Label"] == "Cancer"]
    median_status_synergy = np.median(synergy["BCFS"].values)
    synergy["MedianStatus"] = np.where(synergy["BCFS"]<= median_status_synergy)

    synergy_ = synergy
    neorhea_ = neorhea
    pearl_ = pearl

    synergy_["ctDNADetected"] = np.where(synergy_["ctDNADetected"] == False, "NotDetected", synergy_["ctDNADetected"] )
    synergy_["ctDNADetected"] = np.where(synergy_["ctDNADetected"] == True, "Detected", synergy_["ctDNADetected"])

    pearl_["ctDNADetected"] = np.where(pearl_["ctDNADetected"] == False, "NotDetected", pearl_["ctDNADetected"])
    pearl_["ctDNADetected"] = np.where(pearl_["ctDNADetected"] == True, "Detected", pearl_["ctDNADetected"])

    neorhea_["ctDNADetected"] = np.where(neorhea_["ctDNADetected"] == False, "NotDetected", neorhea_["ctDNADetected"])
    neorhea_["ctDNADetected"] = np.where(neorhea_["ctDNADetected"] == True, "Detected", neorhea_["ctDNADetected"])
    healthy_ = healthy
    healthy_["ctDNADetected"] = "NotDetected"
    all_data_labels = pd.concat([neorhea_, pearl_, synergy_, healthy_], axis=0)
    all_data_labels = pd.merge(all_data_labels, batch, on="PatientId", how="left")
    all_data_labels["batch"] = np.where(all_data_labels["PatientId"] == "Genome-IJB-HP-30-xx_S18", "A00154:1106:H3L22DRX2",
                                 all_data_labels["batch"])
    all_data_labels.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaDataLabels.csv", sep="\t",
                    index=False, header=True)

