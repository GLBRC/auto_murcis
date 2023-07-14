#!/bin/bash
source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
unset $PYTHONPATH
conda activate /home/glbrc.org/user/.conda/envs/murcis
/home/glbrc.org/user/scripts/auto_murcis/murcs_script.py -f bamfiles -t gene_spaces_sequence-Nicole.txt
conda deactivate
