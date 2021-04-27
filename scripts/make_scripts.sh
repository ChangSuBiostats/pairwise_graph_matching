#!/bin/bash

# make a list of jobs, 
# where each job is for matching subject i (i=1,...,997) to all other 997 subjects,
# and save the relevant results.

data_path="/gpfs/ysm/project/zf59/cs/brain_connectivity/precision/data"

echo '' > pairwise_match_joblist.txt
for i_ref in {1..997}
do
	command="try, match_to_a_subj("$i_ref",3e-4,'yeo','"$data_path"'),catch me, fprintf('%s / %s\n',me.identifier,me.message), exit(1), end, exit(0)"
	echo "cd ..; module load MATLAB; matlab -nodisplay -nodesktop -r \""$command"\"" >> pairwise_match_joblist.txt
done
