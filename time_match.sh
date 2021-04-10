# time matching one subject to all other 997 subjects
#
# Chang Su, c.su@yale.edu

# set a reference subject
i_ref=1

# reference:
# 1. make MATLAB take arguments from command line:
#	https://stackoverflow.com/questions/8981168/running-a-matlab-program-with-arguments
# 2. run matlab functions with arguments in a bash script:
#	https://stackoverflow.com/questions/28305451/running-matlab-m-file-with-arguments-in-a-bash-script

time matlab -nodisplay -nodesktop -r "try, match_to_a_subj("$i_ref",3e-4,'two_region','/Users/chang/Documents/research/brain_connectivity/data/precision/'),catch me, fprintf('%s / %s\n',me.identifier,me.message), exit(1), end, exit(0)"

