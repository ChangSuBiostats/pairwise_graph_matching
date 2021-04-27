#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=5:00
#SBATCH --mem-per-cpu=5G

cd ..
module load MATLAB
data_path="/gpfs/ysm/project/zf59/cs/brain_connectivity/precision/data"
matlab -nodisplay -nodesktop -r "try, match_to_a_subj(10,3e-4,'yeo','"$data_path"'),catch me, fprintf('%s / %s\n',me.identifier,me.message), exit(1), end, exit(0)"
