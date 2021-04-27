# time matching one subject to all other 997 subjects
#
# Chang Su, c.su@yale.edu

module load MATLAB

cd ..

# set a reference subject
i_ref=1
data_path="/gpfs/ysm/project/zf59/cs/brain_connectivity/precision/data"

# reference:
# 1. make MATLAB take arguments from command line:
#	https://stackoverflow.com/questions/8981168/running-a-matlab-program-with-arguments
# 2. run matlab functions with arguments in a bash script:
#	https://stackoverflow.com/questions/28305451/running-matlab-m-file-with-arguments-in-a-bash-script

#time matlab -nodisplay -nodesktop -r 'try, match_to_a_subj("$i_ref",3e-4,'yeo','"$data_path"'),catch me, fprintf('%s / %s\n',me.identifier,me.message), exit(1), end, exit(0)'

time matlab -nodisplay -nodesktop -r "try, match_to_a_subj("$i_ref",3e-4,'yeo','"$data_path"'),catch me, fprintf('%s / %s\n',me.identifier,me.message), exit(1), end, exit(0)"

