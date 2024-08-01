#
Rcpp::sourceCpp("all_sums.cpp")

#this one does some better handling of negative gammas
make_a_matrix_of_gammas <- function(start_counts, stop_counts, generations=10, min_start = 100, symmetrical = FALSE){
  
  ############
  #normalize and add pseudocounts
  to_na <- which(start_counts[,3] < min_start)

  start_counts[to_na,3] <- NA
  stop_counts[to_na,3] <- NA
  
  start_counts[,3] <- start_counts[,3] + 1
  stop_counts[,3] <- stop_counts[,3] + 1
  
  #get the relevant categories
  unique_gene1 <- unique(start_counts[,1])
  unique_gene2 <- unique(start_counts[,2])
  
  if (symmetrical){
    both <- unique(c(unique_gene1, unique_gene2))
    unique_gene1 <- both
    unique_gene2 <- both
  }
  
  control_rows <- grep("control",unique_gene1)
  control_cols <- grep("control",unique_gene2)
  
  ############
  #calculate gammas
  gammas <- log(stop_counts[,3]/start_counts[,3],2)/generations
  all_control_gammas <- gammas[intersect(which(start_counts[,1] %in% unique_gene1[control_rows]), 
                                         which(start_counts[,2] %in% unique_gene2[control_cols]))]
  #if there are controls, use them. otherwise use the median
  if (!is.na(median(all_control_gammas, na.rm = TRUE))){
    gammas <- gammas - median(all_control_gammas, na.rm = TRUE) + 1
  } else {
    gammas <- gammas - median(gammas, na.rm = TRUE) + 1
  }

  ############
  #initialize the matrix
  results <- matrix(NA, length(unique_gene1), length(unique_gene2))  
  rownames(results) <- unique_gene1
  colnames(results) <- unique_gene2
  
  #put in the results
  for (i in seq_along(unique_gene1)){
    of_int <- which(start_counts[,1] == unique_gene1[i])
    if (length(of_int) > 0){
    for (j in seq_along(unique_gene2)){
      of_int2 <- of_int[which(start_counts[of_int,2] == unique_gene2[j])]
      if (length(of_int2) > 0){results[i,j] <- gammas[of_int2]}
    }
    }
    if (i%%10 == 0){print(i)}
  }
  
  return(results)
  
}
  
#
get_expected2 <- function(results, use_all_rows = FALSE, use_all_cols = FALSE, use_additive = FALSE, 
                          correlation_mad_threshold = NA, row_sgRNA = NA, col_sgRNA = NA){
  
 
  ############
  #initialize individual sgRNA 
  single_row_gamma <- vector("list", dim(results)[1])
  single_col_gamma <- vector("list", dim(results)[2])
  fraction_not_na_rows <- rep(NA, length(single_row_gamma))
  fraction_not_na_cols <- rep(NA, length(single_col_gamma))
  
  #get control cols
  control_rows <- grep("control",rownames(results))
  control_cols <- grep("control",colnames(results))
  
  #single sgRNAs
  #if we have control cols, use them
  if (length(control_cols) > 0 & !use_all_cols){
    number_total_rows <- length(control_cols)
    for (i in seq_along(single_row_gamma)){
      single_row_gamma[[i]] <- results[i,control_cols]
      to_keep <- !is.na(single_row_gamma[[i]])
      fraction_not_na_rows[i] <- sum(to_keep, na.rm = TRUE)/number_total_rows
      single_row_gamma[[i]] <- single_row_gamma[[i]][to_keep]}
  }
  
  #if not, lets get a set of largely control columns  
  if (length(control_cols) == 0 | use_all_cols){
    col_medians <- apply(results,2,median, na.rm = TRUE)
    temp_control_cols <- which(abs(col_medians-median(col_medians, na.rm = TRUE)) < mad(col_medians, na.rm = TRUE))
    number_total_rows <- length(temp_control_cols)
    for (i in seq_along(single_row_gamma)){
      single_row_gamma[[i]] <- results[i,temp_control_cols]
      to_keep <- !is.na(single_row_gamma[[i]])
      fraction_not_na_rows[i] <- sum(to_keep, na.rm = TRUE)/number_total_rows
      single_row_gamma[[i]] <- single_row_gamma[[i]][to_keep]}
  }
  
  #if we have control rows, else
  if (length(control_rows) > 0 & !use_all_rows){
    number_total_cols <- length(control_rows)
    for (j in seq_along(single_col_gamma)){
      single_col_gamma[[j]] <- results[control_rows,j]
      to_keep <- !is.na(single_col_gamma[[j]])
      fraction_not_na_cols[j] <- sum(to_keep, na.rm = TRUE)/number_total_cols
      single_col_gamma[[j]] <- single_col_gamma[[j]][to_keep]}
  }
  #if not, lets get a set of largely control columns  
  if (length(control_rows) == 0 | use_all_rows){
    row_medians <- apply(results,1,median, na.rm = TRUE)
    temp_control_rows <- which(abs(row_medians-median(row_medians, na.rm = TRUE)) < mad(row_medians, na.rm = TRUE))
    number_total_cols <- length(temp_control_rows)
    for (j in seq_along(single_col_gamma)){
      single_col_gamma[[j]] <- results[temp_control_rows,j]
      to_keep <- !is.na(single_col_gamma[[j]])
      fraction_not_na_cols[j] <- sum(to_keep, na.rm = TRUE)/number_total_cols
      single_col_gamma[[j]] <- single_col_gamma[[j]][to_keep]}
  }
  
  ############
  #zscores
  expected <- results
  mads <- results
  zscores <- matrix(NA, dim(results)[1], dim(results)[2])
  rownames(zscores) <- rownames(results)
  colnames(zscores) <- colnames(results)
  
  #
  for (i in seq_along(single_row_gamma)){
      
      if (fraction_not_na_rows[i] > 0.5){
      temp_vec <- t(single_row_gamma[[i]])
      
      for (j in seq_along(single_col_gamma)){ 
        
        if (fraction_not_na_cols[j] > 0.5){
          
          #multiplicative fitness
          if (!use_additive){
            matrix_of_options <- single_col_gamma[[j]] %*% temp_vec
            temp_median <- Rfast::med(matrix_of_options)
            temp_mad <- Rfast::med(abs(matrix_of_options - temp_median))*1.4826
            zscores[i,j] <- (results[i,j] - temp_median)/temp_mad
            expected[i,j] <- temp_median
            mads[i,j] <- temp_mad
          }
          
          #additive fitness
          if (use_additive){
            matrix_of_options <- all_sums(single_col_gamma[[j]], single_row_gamma[[i]])-1
            temp_median <- Rfast::med(matrix_of_options)
            temp_mad <- Rfast::med(abs(matrix_of_options - temp_median))*1.4826
            zscores[i,j] <- (results[i,j] - temp_median)/temp_mad
            expected[i,j] <- temp_median
            mads[i,j] <- temp_mad
          }
        }
      }
    }
    if (i %% 10 == 0){print(i)}
  }
  
  #
  uncorrected_zscores <- zscores
  
  ########
  #bad sgRNAs by their correlation to sequence
  if (!is.na(correlation_mad_threshold[1]) & !is.na(row_sgRNA[1]) & !is.na(col_sgRNA[1])){
    
    #calculate the row and column correlations for 
    #initialize
    row_corr <- rep(NA, length(row_sgRNA))
    col_corr <- rep(NA, length(col_sgRNA))
    
    row_start <- substr(row_sgRNA,1,2)
    col_start <- substr(col_sgRNA,1,2)
    
    #row
    for (i in seq_along(row_corr)){
      if (sum(is.na(zscores[i,])) < length(col_start)*0.5){
        row_corr[i] <- summary(lm(zscores[i,] ~ col_start))$r.squared
      }
    }
    
    #col
    for (i in seq_along(col_corr)){
      if (sum(is.na(zscores[,i])) < length(row_start)*0.5){
        col_corr[i] <- summary(lm(zscores[,i] ~ row_start))$r.squared
      }
    }
    
    rows_to_remove <- which(row_corr > median(row_corr, na.rm = TRUE)+
                              correlation_mad_threshold*mad(row_corr, na.rm = TRUE))
    cols_to_remove <- which(col_corr > median(col_corr, na.rm = TRUE)+
                              correlation_mad_threshold*mad(col_corr, na.rm = TRUE))
    
    #now, NA out the relevant stuff in results, then ?? refeed the whole thing to the pipeline again!
    results[rows_to_remove,] <- NA
    results[,cols_to_remove] <- NA
    
    #
    zscores[rows_to_remove,] <- NA
    zscores[,cols_to_remove] <- NA
    
  } else {
    rows_to_remove <- NA
    cols_to_remove <- NA
  }
  
  #return
  mads[which(mads == results, arr.ind = TRUE)] <- NA
  expected[which(expected == results, arr.ind = TRUE)] <- NA
  return(list(results, single_row_gamma, single_col_gamma, expected, zscores, rows_to_remove, cols_to_remove, mads,uncorrected_zscores))
  
}
