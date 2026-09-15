#!/bin/bash
#SBATCH --job-name=hsci_mci_ad
#SBATCH --output=hsci_mci_ad_%j.log
#SBATCH --error=hsci_mci_ad_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/blue/nessie/mostafarezaali/400M_PRISM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mostafarezaali@ufl.edu
set -euo pipefail
DATA_DIR=${DATA_DIR:-/blue/nessie/mostafarezaali/400M_PRISM}
CODE_DIR=${CODE_DIR:-$DATA_DIR/HSCI_Alzheimer}
PY=${PY:-/blue/nessie/mostafarezaali/.conda/envs/torch_b200/bin/python}
: "${PATIENT_CSV:?Set PATIENT_CSV to the patient CSV filename or full path}"
unset PYTHONHOME PYTHONPATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$DATA_DIR"
"$PY" -u "$CODE_DIR/exposures.py" --data-dir "$DATA_DIR" --patient-csv "$PATIENT_CSV" \
  --id-col "${ID_COL:-ID}" --zip-col "${ZIP_COL:-zip5}" \
  --mci-col "${MCI_COL:-MCI_DATE}" --ad-col "${AD_COL:-AD_DATE}" \
  --date-format "${DATE_FORMAT:-%Y-%m-%d}" --pcts "${PCTS:-90,95}" \
  --workers "${SLURM_CPUS_PER_TASK:-64}"
