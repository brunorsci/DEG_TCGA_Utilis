#'@title PlotDEG
#'@description
#'This function plot Edger DEG results in a plot that takes log Fold Changes 
#'for each genes with different colors for better interpretation.
#'@param XData the results as showed by edger packages.
#'@param cutoff_up the cutoff for up regulated genes.
#'@param cutoff_down the cutoff for down regulated genes.
#'@param Title the phrase to be exibited as title.
#'@importFrom grDevices rainbow
#'@importFrom stats filter
#' Para utilizar: PlotDEG()
#'@author Bruno Rodrigo Assunção(2024)

PlotDEG<-function(XData=Edge_LumAvsLumB,cutoff_up=1, cutoff_down=-1, Title="Differentially Expressed Genes"){
  
  # cutoff 
  cutoff1 <- as.numeric(cutoff_up) 
  cutoff2 <- as.numeric(cutoff_down) 
  
  # Unique color for genes 
  cores_unicas <- rainbow(length(rownames(XData))) 
  XData <- XData %>%   
    mutate(color = cores_unicas[match(rownames(XData), 
                                      unique(rownames(XData)))]) 
  # Creating the plot
  grafico <- ggplot(data = XData, aes(x = logCPM, y = logFC)) +   
    geom_point(aes(color = as.factor(rownames(XData))), size = 3) +   
    geom_text(hjust = 1, vjust = 1, size = 1, aes(label = "")) +  # Não exibir o nome do gene no gráfico      
    labs(x = "log CPM", y = "Log Fold Change") +   
    ggtitle(paste(Title)) +   
    scale_color_manual(values = cores_unicas, name = "Genes", 
                       guide = guide_legend(ncol = 2)) +  
    # Color of the legend
    theme_minimal() +   geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black")# Add lines in x
  
  # To plot
  plot(grafico)
}