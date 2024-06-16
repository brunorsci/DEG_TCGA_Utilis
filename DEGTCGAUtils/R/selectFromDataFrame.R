#'@title selectFromDataFrame
#'@description
#'This function is used to takes annotation of another data frame by a shared id/ barcode.
#'It was used in the DEG-TCGA-WORKFLOW (Bruno R. Assuncao, 2024. <brunor> RPubs).
#'Only for that propose and that use. That need to be modified if for other purpose.

selectFromDataFrame=function(xdata){
  selected=NULL
  for (i in xdata) {
    s = info_amostras_artigo %>%
      filter(tissue == 'BRCA') %>%
      filter(barcode==i) 
    selected<-rbind(selected, s)
  }
  return(selected)
}
