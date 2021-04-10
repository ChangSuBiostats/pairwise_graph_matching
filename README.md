# Pairwise Graph Matching
Codes for matching all HCP precision FC maps to one anotherusing Hungarian algorithm 

## Main modules

1. Graph matching algorithm
  iterative_procedure.m
2. Matching one subject to all others in HCP data
  match_to_a_subj.m
3. Others
  3.1 Hierarchical clustering to generate functional based clusters
    hclust_on_precision_FC.m
  3.2 Carefully designed penalty matrices
    generate_penalty_matrix.m
    
## Pre-computed data

In output/, we have two pre-computed penalty matrices, left & right labels and hclust results for cc400.
