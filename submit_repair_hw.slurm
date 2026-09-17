#!/bin/bash
#SBATCH --job-name=repair_hsci_hw
#SBATCH --output=repair_hsci_hw_%j.log
#SBATCH --error=repair_hsci_hw_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --chdir=/blue/nessie/mostafarezaali/400M_PRISM
set -euo pipefail
DATA_DIR=${DATA_DIR:-/blue/nessie/mostafarezaali/400M_PRISM}
CODE_DIR=${CODE_DIR:-$DATA_DIR/Alzheimer_Study/HSCI_Alzheimer}
INPUT_DIR=${INPUT_DIR:-$DATA_DIR/HSCI_Alzheimer_outputs}
OUTPUT_DIR=${OUTPUT_DIR:-$DATA_DIR/HSCI_Alzheimer_outputs_HW_fixed}
PY=${PY:-/blue/nessie/mostafarezaali/.conda/envs/torch_b200/bin/python}
unset PYTHONHOME PYTHONPATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PY" -u "$CODE_DIR/repair_hw_columns.py" --data-dir "$DATA_DIR" \
  --input-dir "$INPUT_DIR" --output-dir "$OUTPUT_DIR" --pcts "${PCTS:-90,95}"
