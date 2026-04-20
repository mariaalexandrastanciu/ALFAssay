#!/bin/bash
#
#SBATCH --job-name=FSC
#SBATCH --output=job_logs/fragmentSizeCountingMulti%A_%a.txt
#SBATCH --array=0-97
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --time=03:00:00

#SBATCH --mem-per-cpu=3000

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export TMP=$LOCALSCRATCH
export TEMP=$LOCALSCRATCH
export TMPDIR=$LOCALSCRATCH


date

echo "Task ID: $SLURM_ARRAY_TASK_ID"

ml releases/2021b Python/3.9.6-GCCcore-11.2.0
source ~/cfDNAEnv/bin/activate

cd ALFAssay

study="healthy"
analysis="NN5Mb"
genome_v="hg38"
study_dir=/globalscratch/ulb/bctr/astanciu/$study

mapping_quality=60
size_fragment_start=30
size_fragment_end=700
window_size=5000000
by_window_flag=True
gc_correction_type=GCParagon
filename_key=_fragment_size_summary_window.csv
filename_key_simple_summary=_fragment_size_summary.csv
bamDir=${study_dir}/GCParagonBams
bamFiles=($(find $bamDir -type f -name "*bam"))
bamFile=${bamFiles[$SLURM_ARRAY_TASK_ID]}
id=$(basename $bamFile .markDup.GCtagged.bam)
mkdir -p ${study_dir}/FragmentationPatterns/
mkdir -p ${study_dir}/FragmentationPatterns/${analysis}
outDir=${study_dir}/FragmentationPatterns/${analysis}

echo $id
if [ -f $outDir/${id}_fragment_size_summary.csv ]; then
        echo "file $id exists"
else

python -Wignore  -m FragmentsManipulations.runFragmentSizeCountMultiprocess -b $bamFile -o $outDir/$id  -sfs $size_fragment_start \
        -sfe $size_fragment_end -gv $genome_v -t22 -ws $window_size -wf $by_window_flag -mp $mapping_quality -gcCorrection $gc_correction_type


echo "Fragment size counts calculated at: $outDir/$id.csv "

gzip $outDir/${id}_fragment_size_expanded.bed


fi