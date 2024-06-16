#'@title convertToExpression.
#'@description
#'The function takes a data frame of rna expression counts and transforms it to 
#'a likely used format required as input for packages like edger and Deseq2.
#'@param df the data frame do be converted.
#'@param col_sp_names the collum with the sample names ids.
#'@param col_gene_names the collum with genes names ids.
#'@param col_counts the collum with the counts per gene.
#'@return A data frame called df_modified. it contains the data with only sample names, gene names and counts
#'in the expression format.
#' Para utilizar: convertToExpression()
#'@author Bruno Rodrigo Assunção(2024).
 

convertToExpression=function(df=expn, col_sp_names="barcode",                                    
                             col_gene_names="gene_id", 
                             col_counts="normalized_read_count"){ 
  require(dplyr)
  # sample names col 1, gene names col 2 e raw counts col 3 
  tempa <- df 
  tempc <- unique(tempa[[col_gene_names]])      
  tempg <- NULL   
  for (i in 1:length(tempc)) {
    tempd=tempa[tempa[,col_gene_names] %in% tempc[i], c(col_sp_names, col_counts)]     
    tempd=tempd[order(tempd[[col_sp_names]]),] 
    
    #order genes     
    tempe <- data.frame(t(tempd))     
    #temp uniques rows of final data.frame       
    tempf <- lapply(tempe[-1,-1], as.numeric)
    tempf <- as.data.frame(t(unlist(tempf)))
    colnames(tempf) <- t(tempd)[1,-1] #collumn name == samples names
    rownames(tempf) <- tempc[i] #rownames == gene(i)
    assign(paste0("df_",i), tempf)
    #tempg=rbind(tempg, tempf)   
  }   
  
  df_modified=NULL   
  for (e in 1:length(tempc)) {    
    #df_modified=rbind(df_modified, get(paste0("df_",e)))     
    # Modified to include dfs thas were duplicated and then different number of collumns    
    df_modified <- dplyr::bind_rows(df_modified, distinct(get(paste0("df_",e))))   
  } 
  return(df_modified) 
}
