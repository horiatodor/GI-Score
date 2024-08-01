#join_get_expected2

#
join_get_expected2 <- function(object1, object2){
  
  #quit if objects arent the samel length
  if (length(object1) != length(object2)){return(NULL)}
  
  #initialize
  new_list <- vector("list", length(object1))
  
  #loop to join (cbind or c)
  for (a in seq_along(object1)){
    if (!is.null(dim(object1[[a]]))){
      new_list[[a]] <- cbind(object1[[a]], object2[[a]])
    } else {
      new_list[[a]] <- c(object1[[a]], object2[[a]])
    }
  }
  
  #return
  return(new_list)
}