# 🧬 ALFAssay
**A Python-based machine learning pipeline for circulating tumor DNA (ctDNA) quantification in metastatic breast cancer.**

![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![Status](https://img.shields.io/badge/status-active-success.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

ALFAssay leverages fragmentomic features derived from cell-free DNA (cfDNA) sequencing data to model and predict ctDNA tumor fractions. It is designed for large-scale cfDNA analysis and extracts highly specific predictive features from GC-corrected BAM files.

---

## 🧠 Overview

The ALFAssay pipeline is composed of four main stages:
1. **Preprocessing** – Extract fragmentomic signals from BAM files.
2. **Feature Engineering** – Build model-ready feature matrices.
3. **Model Training** – Train and evaluate neural network models.
4. **Prediction** – Apply trained models to new samples.

---

## ⚙️ Installation

This project uses [uv](https://github.com/astral-sh/uv) for lightning-fast dependency management and isolated virtual environments. 

### 1. Install `uv`
If you do not have `uv` installed, install it globally:
```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

### 2. Set up the Environment

**For Windows, Linux, and Intel Macs:**
Simply clone the repository and sync the environment:
```bash
uv sync
```

**⚠️ For Apple Silicon Macs (M1/M2/M3/M4):**
Because this project relies on legacy bioinformatics C-extensions (`sorted_nearest`, `ncls`) that do not support native Apple Silicon, you must build the environment using Intel Emulation (Rosetta) before syncing:
```bash
# 1. Force an emulated Intel Python environment
uv venv --python "3.10-x86_64" --seed

# 2. Sync dependencies
uv sync
```

---

## 🚀 Pipeline Usage

### Step 1: Pre-processing (BAM to Features)
**Prerequisites:** 
- Paired FASTQ files must be aligned, and GC correction must be applied using [GC Paragon](https://github.com/BGSpiegl/GCparagon) (using the GC tag).
- Place the GC-corrected BAM files (and their indexes) inside `ALFAssay/FragmentsManipulations/data/input`.

**Command:**
```bash
uv run -m FragmentsManipulations.runFragmentSizeCountMultiprocess \
  -b sample_gc_corrected.bam -o ./out/sample \
  -sfs 30 -sfe 700 -gv hg38 -t 22 \
  -ws 100000 -wf True -mp 60
```

**Outputs:**
Generated in `ALFAssay/FragmentsManipulations/data/output`:
1. `*_fragment_size_summary.csv`: Overall fragment count.
2. `*_fragment_size_summary_window.bed`: Fragment counts on window size.
3. `*_fragment_size_expanded.bed`: Detailed fragment data (chromosome, start, end, GC content, GC Paragon weight, 4/6 bp motif ends).
4. `*_nn_data.bed`: Fragment summaries per bin.

**Feature details inside `*_nn_data.bed`:**
| Feature Name | Description | Size Range (bp) | Type |
| --- | --- | --- | --- |
| `ultraShortFrag` | Number of ultra-short fragments | `< 90` | Count |
| `shortFrag` | Number of short fragments | `90 ≤ size < 151` | Count |
| `longFrag` | Number of long fragments | `151 ≤ size < 221` | Count |
| `ultraLongFrag` | Number of ultra-long fragments | `> 220` | Count |
| `ultraShortFragGCPara` | Sum of GC-weighted ultra-short fragments | `< 90` | Weighted sum |
| `shortFragGCPara`| Sum of GC-weighted short fragments | `90 ≤ size < 151` | Weighted sum |
| `longFragGCPara` | Sum of GC-weighted long fragments | `151 ≤ size < 221`| Weighted sum |
| `ultraLongFragGCPara` | Sum of GC-weighted ultra-long fragments | `> 220` | Weighted sum |

---

### Step 2: Feature Engineering
**Prerequisites:**
1. **Fragment Data Files:** Ensure the `*_nn_data.bed` files generated in Step 1 are available in `ALFAssay/FragmentsManipulations/data/output/`.
2. **Metadata File:** A patient-level metadata file must exist with the following columns:

| Column | Description |
| --- | --- |
| `PatientId` | Unique identifier matching the fragment file names. |
| `Label` | Clinical label (e.g., cancer / healthy). |
| `ctDNADetected` | Binary indicator of circulating tumor DNA presence (`ichorTF > 0.03`). |
| `VAF` | Variant Allele Frequency (can be ichorTF tumor fraction). |
| `ichorTF` | Tumor fraction estimated via [ichorCNA](https://github.com/broadinstitute/ichorCNA). |
| `study` | Study or cohort identifier. |

**Command:**
```bash
uv run -m FeatureEngineering.gather_features
```

---

### Step 3: Model Training
**Prerequisites:**
1. Metadata file (e.g., `refdata/AllMetaData.csv`).
2. Engineered features matrix (e.g., `FeatureEngineering/data/NN_5MB_features.csv`).

**Command:**
```bash
uv run -m NNModel.main
```

**Outputs:**
- `NNModel/results/results_cv_epoch_<epoch>_model.csv`: Cross-validation predictions.
- `NNModel/results/results_validation_ichorTF_model.csv`: External validation predictions.
- `NNModel/models/ALFAssay_model.pkl`: Final trained model.
- `NNModel/plots/Loss function value over the epochs.png`: Plot of training vs. test loss.

---

### Step 4: Prediction
Use the `predict_tf.py` script to generate tumor fraction predictions from precomputed feature files using a previously trained model. *(This step does not perform any training).*

**Command:**
```bash
uv run -m NNModel.predict_tf FeatureEngineering/data/NN_5MB_pivot.csv
```

**Optional Arguments:**
```bash
# Provide custom paths for the model and output destination
uv run -m NNModel.predict_tf FeatureEngineering/data/NN_5MB_pivot.csv \
  --model NNModel/models/ALFAssay_model.pkl \
  --output NNModel/results/predictions.tsv
```

**Output format (`predictions.tsv`):**
| Column | Description |
| --- | --- |
| `PatientId` | Patient/sample identifier |
| `PredictedValue` | Predicted tumor fraction value |
| `TrueValue` | True `ichorTF` value |

---

## 🏁 End-to-End Example

To run the entire pipeline end-to-end on provided example data:

```bash
uv run -m FragmentsManipulations.example_run
uv run -m FeatureEngineering.gather_features
uv run -m NNModel.main
uv run -m NNModel.predict_tf FeatureEngineering/data/NN_5MB_pivot.csv
```

---

## 💡 Notes

- **GC correction is mandatory** for accurate predictions.
- Patient IDs must strictly match across all BAM files, BED files, and Metadata sheets.
- ALFAssay is highly optimized for High-Performance Computing (HPC) environments and supports SLURM execution.

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).