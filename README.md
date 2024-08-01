# GI-Score
Code used to calculate genetic interaction (GI) scores for Double-CRISPRi

# Usage for Double-CRISPRi
For each replicate:

1) Make a matrix of RF using the count files

```rep1_temp <- make_a_matrix_of_gammas(read.csv(rep1_files[comparisons[i,1]], stringsAsFactors = FALSE), read.csv(rep1_files[comparisons[i,2]], stringsAsFactors = FALSE), generations = comparisons[i,3])```
   
2) Use the matrix of RFs to calculate GI scores for both parts

```rep1_vals_temp1 <- get_expected2(rep1_temp[,1:316],use_additive = TRUE, use_all_rows = TRUE, use_all_cols = TRUE, correlation_mad_threshold = 5,row_sgRNA = matched_sgRNAs_row, col_sgRNA = matched_sgRNAs[1:316])```
```rep1_vals_temp2 <- get_expected2(rep1_temp[,317:1310],use_additive = TRUE, use_all_rows = TRUE, use_all_cols = TRUE, correlation_mad_threshold = 5, row_sgRNA = matched_sgRNAs_row, col_sgRNA = matched_sgRNAs[317:1310])```

3) Join the independently calculated chunks

```rep1_vals[[i]] <- join_get_expected2(rep1_vals_temp1, rep1_vals_temp2)```

4) Do all the downstream analysis

