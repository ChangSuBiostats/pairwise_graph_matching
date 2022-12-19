# Pairwise Graph Matching
Codes for matching all HCP precision FC maps to one another using Hungarian algorithm.

## Reference
**A graph-matching based metric of functional connectome distance between pairs of individuals varies with their ages, cognitive performances and familial relationships**.
Hussain Bukhari, Chang Su, Elvisha Dhamala, Zijin Gu, Keith Jamison, Amy Kuceyeski.
bioRxiv 2022.10.03.510660; doi: https://doi.org/10.1101/2022.10.03.510660

## Functions
1. Graph matching algorithm: iterative_procedure.m
2. Matching one subject to all test subjects in HCP data: match_to_a_subj.m
3. Others
    * Carefully designed penalty matrices: generate_penalty_matrix.m
    * Sample codes for analyzing the results: sample_results.m
    * Hierarchical clustering to generate functional based clusters: hclust_on_precision_FC.m

## Notes
'Graph_matching_algorithm.pdf' is a short note on the algorithm we implemented.

## Pre-computed data
In output/, we have three pre-computed penalty matrices, labels for left and right hemishpere, hclust results for cc400, etc.

## MISC
1. Helper functions for visualization: visualization_tools
3. Bash script used for computing all permutation matrices: /scripts
