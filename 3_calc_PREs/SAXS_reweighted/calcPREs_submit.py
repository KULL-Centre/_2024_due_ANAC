import os
import subprocess
import time
from jinja2 import Template
from analyse import *

proteins = initProteins()
proteins.to_pickle('proteins.pkl')

submission = Template("""#!/bin/bash
#SBATCH --job-name=PRE{{name}}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=2000
#SBATCH -t 20:00:00
#SBATCH -o PRE_{{name}}.out
#SBATCH -e PRE_{{name}}.err
#SBATCH --partition=sbinlab_ib

#source /groups/sbinlab/giulio/.bashrc
#conda activate hoomd

python=/sbinlab/giulio/miniconda3/envs/hoomd/bin/python3.8

start=$(date +%s.%N)

$python calcPREs.py --name {{name}}

duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`

echo $execution_time""")

name='ANAC046'
with open('PRE_{:s}.sh'.format(name), 'w') as submit:
    submit.write(submission.render(name=name))
subprocess.run(['sbatch','PRE_{:s}.sh'.format(name)])
time.sleep(20)
