# DEG_TCGA_Utilis
To process TCGA expression gene counts data to  DEG analysis by Edger and DESeq2 packages.

---
title: "DEG-TCGA Utils"
author: "Bruno Rodrigo Assunção"
date: "2024-05-02"
output:
  pdf_document: default
  html_document: default
---

Expressão Diferencial de Genes em Câncer de Mama: workflow bioinformático na linguagem de programação R

Esse é um workflow escrito na linguagem de programação R para análise de genes diferencialmente expressos de glicosiltransferases em Transcriptoma de Cancer de Mama obtidos a partir do repositorio ICGC Data Portal.

Dependências:

library("DESeq2") 
library("dplyr") 
library("edgeR")
library("ggplot2") 
library(tidyverse)
library(DEGTCGAUtils) # Essa última foi criada para esse fim

FUN

00 Entrada de dados

rna <- read.delim("D:/brca_eua_atual/exp_seq.BRCA-US.tsv") 
info_amostras_artigo <- read.csv("D:/brca_eua_atual/info_amostras_artigo.csv")

Selecionamos apenas alguns genes de interesse. Os genes foram selecionados da forma
à seguir à partir da lista abaixo.

rna <- as.data.frame(rna) 
rna <- rna[rna$gene_id %in% c("ALG10B","ALG6",  "ALG8",  "A4GALT","A4GNT",  "ALG14"), ] 

Anotação de informações de subtipos moleculares e histológicos

De modo a obter informação confiável e validada a respeito do subtipo histológico/molecular das amostras obtidas, utilizamos as informações de um artigo cientifíco que utilizou esse conjunto experimental.

Anotar dados a partir do artigo

Função para seleção de anotação a partir do artigo:

# Outro modo de fazer a mesma coisa 
library(tidyverse) 
selectFromDataFrame <- function(xdata){ 
  selected=NULL   
  for (i in rna$submitted_sample_id) {     
    s = info_amostras_artigo %>%     
    filter(tissue == 'BRCA') %>%     
    filter(barcode==i)      
    selected <- rbind(selected, s)   
    }   
return(selected) 
}

Filtragem e anotação dos subtipos moleculares

Anotamos as informações dos subtipos molecuares à partir da tabela informações das amostras. Isso não é necessário caso essa informação esteja disponível no conjunto de dados.

selected=NULL
  for (i in rna$submitted_sample_id) {
    s = info_amostras_artigo %>%
    filter(tissue == 'BRCA') %>%
    filter(barcode==i) 
    selected<-rbind(selected, s)
  }
  
colnames(rna) <-c(colnames(rna)[1:4],"barcode", colnames(rna)[6:22])
exp=merge(rna, distinct(selected), by="barcode")# distinct paralinhas unicas

Modificando o nome das colunas e combinando os resultados da busca acima:

colnames(rna) <- c(colnames(rna)[1:4],"barcode", colnames(rna)[6:22]) 
exp <- merge(rna, distinct(selected), by="barcode") # distinct para linhas unicas

Converter para formato de expressão:

Modificar a tabela de expressão de modo que as colunas serão linhas e as linhas colunas no novo objeto. Antes, vamos observar como os dados estão organizados nessa tabela. Apenas as primeiras linhas e algumas colunas de interesse.

Os dados de contagens (expressão dos genes/rnas mensageiros) e os nomes das amostras estão dispostos de forma inversa as requeridas pelas bibliotecas de análise de DEGs, como DESeq2 e EdgeR, onde os genes estão nas colunas e as amostras nas linhas. Basta olhar os exemplos desses programas acima. Por esse motivo, criamos uma função para fazer essa "transposição".

Esta função espera um data frame com apenas três colunas na seguinte ordem: Samples (ids) na primeira coluna, genes (nomes) na segunda coluna e a contagem dos genes na terceira coluna.

Por isso, vamos selecionar apenas essas do objeto rna.

# Contagens normalizadas 
expn <- select(.data = exp, c("barcode", "gene_id", "normalized_read_count", "subtype")) 
# Contagens brutas (não normalizadas) 
expr <- select(.data = exp, c("barcode", "gene_id", "raw_read_count", "subtype"))

library(DEGTCGAUtils)

# Para a tabela com contagens normalizadas 
expn_new <- convertdf_to_expression(df = as.data.frame(expn), col_sp_names = "barcode", col_gene_names = "gene_id",  col_counts="normalized_read_count") 

# Para a tabela com contagens brutas 
expr_new <- convertdf_to_expression(df = expr, col_sp_names = "barcode", col_gene_names = "gene_id", col_counts="raw_read_count")

Vamos observar como ficou organizado os dados destas tabelas após a modificação.

{r}
expn_new[1:5,1:5]

Note que cada gene está distribuído no lugar das colunas e nas linhas estão distribuídas as contagens de cada gene respectivo. A tabela expn contém as contagens normalizadas, como obtida do banco de dados e preferidas ao se trabalhar com a biblioteca EDGER. Enquanto que a tabela expr contém os dados com contagens brutas, as quais são esperadas quando se trabalha com a biblioteca DESeq2.

3 EDGER

Análise de Genes Diferencialmente Expressos Utilizando a biblioteca edgeR.

{r include=FALSE}
library(edgeR)

Visualizar os dados

expn_new[1:5,1:3]

3.1 Para estimativa de disperção comum e dispersão tagwise em única corrida

Aqui não utilizamos a filtragem para valores baixos, devido aos dados estarem normalizados e assim serem valores muito pequenos. Em caso de utilização de contagens brutas (raw counts), como esperado pelo pacote DESeq2, utilizamos um valor de corte (como definido pelas linhas mutadas abaixo).

DGE <- DGEList(counts=expr_new, group=selecao_amostras_artigo$subtype) 
DGE$samples 
group=selecao_amostras_artigo$subtype 
keep <- filterByExpr(y = DGE, group=group,  min_count = 1) 
DGE <- DGE[keep, , keep.lib.sizes=T] 
DGE <- normLibSizes(DGE) 
DGE <- estimateCommonDisp(DGE)

3.2 Gerado as comparacoes

Edge_LumAvsLumB <- exactTest(DGE, pair=c("LumA", "LumB")) 
Edge_LumAvsLumB <- as.data.frame(topTags(lumAvsLumB))  
Edge_LumAvsBasal <- exactTest(DGE, pair=c("LumA", "Basal")) 
Edge_LumAvsBasal <- as.data.frame(topTags(lumAvsBasal))  
Edge_LumAvsHer2 <- exactTest(DGE, pair=c("LumA", "Her2")) 
Edge_LumAvsHer2 <- as.data.frame(topTags(lumAvsHer2))  
Edge_LumBvsBasal <- exactTest(DGE, pair=c("LumB", "Basal")) 
Edge_LumBvsBasal <- as.data.frame(topTags(LumBvsBasal))  
Edge_LumBvsHer2 <- exactTest(DGE, pair=c("LumB", "Her2")) 
Edge_LumBvsHer2 <- as.data.frame(topTags(LumBvsHer2)) 
Edge_BasalvsHer2 <- exactTest(DGE, pair=c("Basal", "Her2")) 
Edge_BasalvsHer2 <- as.data.frame(topTags(BasalvsHer2))

3.3 Definindo up ou Down de acordo ao cutoff

comparations<- c("Edge_LumAvsLumB","Edge_LumAvsBasal","Edge_LumAvsHer2","Edge_LumBvsBasal", "Edge_LumBvsHer2","Edge_BasalvsHer2")

# Loop for
for (i in comparations) {
  assign(paste(i),
         get(i) %>%   mutate(color = case_when(logFC > cutoff1 ~ "Up regulated", 
                           logFC < cutoff2 ~ "Down Regulated",  
                           TRUE ~ "Not Differentially Expressed"))) 
}

# TO add a color collumn

#for (i in comparations) {
#  i <- get(i) %>%   
#    mutate(color = case_when(logFC > cutoff1 ~ "Up regulated", 
#                             logFC < cutoff2 ~ "Down Regulated",  
#                             TRUE ~ "Not Differentially Expressed"))
#}

3.4 Exportação dos resultados.

write.table(x = lumAvsLumB, file = "Edge_lumAvsLumB.txt", sep = "\t")  
write.table(x = lumAvsBasal, file = "Edge_lumAvsBasal.txt", sep = "\t")  
write.table(x = lumAvsHer2, file = "Edge_lumAvsHer2.txt", sep = "\t")  
write.table(x = LumBvsBasal, file = "Edge_LumBvsBasal.txt", sep = "\t")  
write.table(x = LumBvsHer2, file = "Edge_LumBvsHer2.txt", sep = "\t")  
write.table(x = BasalvsHer2, file = "Edge_BasalvsHer2.txt", sep = "\t")

3.4 Visualização Gráfica

PlotMD do pacote edger

plotMD(Edge_LumAvsLumB[,-5])

Plot Master: plot modficado.

PlotDEG(XData = Edge_BasalvsHer2, Title = "Differentially Expressed Genes: Basal vs Her2")

DESEQ2

Análise de Genes Diferencialmente Expressos Utilizando a biblioteca DESeq2.


write.table(x = expr_new, file = "counts_exprnew.txt", sep = "\t")  
write.table(x = lumAvsLumB, file = "lumAvsLumB.txt", sep = "\t")  

dds <- DESeq2::DESeqDataSetFromMatrix(countData=expr_new, 
                              colData=selecao_amostras_artigo[,c(2,4)], 
                              design=~subtype)

# To exclude low counts

#keep <- rowSums(counts(dds) >= 10)
#dds <- dds[keep,]

# Run deseq2
dds <- DESeq2::DESeq(dds)

DESeq2::results(dds, alpha = 0.05, lfcThreshold =1)

# Gerando as comparações

assign(paste0("res_DESeq2_LumA_LumB"), DESeq2::results(dds, contrast=c("subtype","LumA","LumB") , alpha = 0.05, lfcThreshold =1))
assign(paste0("res_DESeq2_LumA_Basal"), DESeq2::results(dds, contrast=c("subtype","LumA","Basal") , alpha = 0.05, lfcThreshold =1))
assign(paste0("res_DESeq2_LumA_Her2"), DESeq2::results(dds, contrast=c("subtype","LumA","Her2") , alpha = 0.05, lfcThreshold =1))
assign(paste0("res_DESeq2_LumB_Basal"), DESeq2::results(dds, contrast=c("subtype","LumB","Basal") , alpha = 0.05, lfcThreshold =1))
assign(paste0("res_DESeq2_LumB_Her2"), DESeq2::results(dds, contrast=c("subtype","LumB","Her2") , alpha = 0.05, lfcThreshold =1))
assign(paste0("res_DESeq2_Basal_Her2"), DESeq2::results(dds, contrast=c("subtype","Basal","Her2") , alpha = 0.05, lfcThreshold =1))

Gerar tabelas para Genes Up and Down regulated

# Visualize 'DEGs' Up and Down regulated
for (i in grep("res_DESeq2", ls(), value = TRUE)) {
  assign(paste0("DEG_UP_",i), as.data.frame(subset(get(i),get(i)$log2FoldChange > 1 & get(i)$padj < 0.05)))
  assign(paste0("DEG_DOWN_",i),as.data.frame(subset(get(i),get(i)$log2FoldChange < 1 & get(i)$padj < 0.05)))
}

Exportar a tabela de resultados

# Down regulated
write.table(x = DEG_DOWN_res_DESeq2_Basal_Her2, file = "DEG_DOWN_res_DESeq2_Basal_Her2.txt", sep = "\t")  
write.table(x = DEG_DOWN_res_DESeq2_LumA_Basal, file = "DEG_DOWN_res_DESeq2_LumA_Basal.txt", sep = "\t")  
write.table(x = DEG_DOWN_res_DESeq2_LumA_Her2, file = "DEG_DOWN_res_DESeq2_LumA_Her2", sep = "\t")  
write.table(x = DEG_DOWN_res_DESeq2_LumA_LumB, file = "DEG_DOWN_res_DESeq2_LumA_LumB.txt", sep = "\t")  
write.table(x = DEG_DOWN_res_DESeq2_LumB_Basal, file = "DEG_DOWN_res_DESeq2_LumB_Basal.txt", sep = "\t")  
write.table(x = DEG_DOWN_res_DESeq2_LumB_Her2, file = "DEG_DOWN_res_DESeq2_LumB_Her2.txt", sep = "\t")

# Up regulated
write.table(x = DEG_UP_res_DESeq2_Basal_Her2, file = "DEG_DOWN_res_DESeq2_Basal_Her2.txt", sep = "\t")  
write.table(x = DEG_UP_res_DESeq2_LumA_Basal, file = "DEG_DOWN_res_DESeq2_LumA_Basal.txt", sep = "\t")  
write.table(x = DEG_UP_res_DESeq2_LumA_Her2, file = "DEG_DOWN_res_DESeq2_LumA_Her2", sep = "\t")  
write.table(x = DEG_UP_res_DESeq2_LumA_LumB, file = "DEG_DOWN_res_DESeq2_LumA_LumB.txt", sep = "\t")  
write.table(x = DEG_UP_res_DESeq2_LumB_Basal, file = "DEG_DOWN_res_DESeq2_LumB_Basal.txt", sep = "\t")  
write.table(x = DEG_UP_res_DESeq2_LumB_Her2, file = "DEG_DOWN_res_DESeq2_LumB_Her2.txt", sep = "\t")

Agora é avaliar os dados e comparar os resultados gerados pelos dois programas.

