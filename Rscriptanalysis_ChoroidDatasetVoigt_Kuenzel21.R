library(plotly)
library(ggrepel)
library(gtools)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gplots)
library(useful) 
library(Matrix)
library(Seurat)
library(Matrix.utils)
library(future)
library(xlsx)
library(MixGHD)

options(future.globals.maxSize = 4000 * 1024^2)

#########################read in and prepare the count matrices of all batches (7 donors with periphery and macula each)
matrixhumanchorrpe1 <- read.table(".../GSM4037981_macula_donor_1_expression.tsv")
matrixhumanchorrpe2 <- read.table(".../GSM4037982_macula_donor_2_expression.tsv")
matrixhumanchorrpe3 <- read.table(".../GSM4037983_macula_donor_3_expression.tsv")
matrixhumanchorrpe4 <- read.table(".../GSM4037984_peripheral_donor_1_expression.tsv")
matrixhumanchorrpe5 <- read.table(".../GSM4037985_peripheral_donor_2_expression.tsv")
matrixhumanchorrpe6 <- read.table(".../GSM4037986_peripheral_donor_3_expression.tsv")
matrixhumanchorrpe7 <- read.table(".../GSM4037987_macula_donor_4_enriched_expression.tsv")
matrixhumanchorrpe8 <- read.table(".../GSM4037988_macula_donor_5_enriched_expression.tsv")
matrixhumanchorrpe9 <- read.table(".../GSM4037989_macula_donor_6_enriched_expression.tsv")
matrixhumanchorrpe10 <- read.table(".../GSM4037990_macula_donor_7_enriched_expression.tsv")
matrixhumanchorrpe11 <- read.table(".../GSM4037991_peripheral_donor_4_enriched_expression.tsv")
matrixhumanchorrpe12 <- read.table(".../GSM4037992_peripheral_donor_5_enriched_expression.tsv")
matrixhumanchorrpe13 <- read.table(".../GSM4037993_peripheral_donor_6_enriched_expression.tsv")
matrixhumanchorrpe14 <- read.table(".../GSM4037994_peripheral_donor_7_enriched_expression.tsv")

processmatrixlist<-list(matrixhumanchorrpe1, matrixhumanchorrpe2, matrixhumanchorrpe3, matrixhumanchorrpe4, matrixhumanchorrpe5, matrixhumanchorrpe6, matrixhumanchorrpe7, matrixhumanchorrpe8, matrixhumanchorrpe9, matrixhumanchorrpe10, matrixhumanchorrpe11, matrixhumanchorrpe12, matrixhumanchorrpe13, matrixhumanchorrpe14)

Metaprocessmatrixlist<-processmatrixlist

for(i in 1:14){
  
  Metaprocessmatrixlist[[i]]<-processmatrixlist[[i]][,1:3]
  
}

Metaprocessmatrixlist

minohnenull <- function(x){
  x <- x[x!=0]
  return(min(x))
}


getcountmatrix<-function(matrix){
  
  matrix<-matrix[,4:dim(matrix)[2]]
  matrix<-exp(1)^matrix
  matrix<-matrix-1
  
  minima<-apply(matrix,1,minohnenull)
  countmatrix<-matrix/minima
  
  return(countmatrix)
  
}



countmatrixhumanchorrpe1 <- getcountmatrix(matrixhumanchorrpe1)
countmatrixhumanchorrpe2 <- getcountmatrix(matrixhumanchorrpe2)
countmatrixhumanchorrpe3 <- getcountmatrix(matrixhumanchorrpe3)
countmatrixhumanchorrpe4 <- getcountmatrix(matrixhumanchorrpe4)
countmatrixhumanchorrpe5 <- getcountmatrix(matrixhumanchorrpe5)
countmatrixhumanchorrpe6 <- getcountmatrix(matrixhumanchorrpe6)
countmatrixhumanchorrpe7 <- getcountmatrix(matrixhumanchorrpe7)
countmatrixhumanchorrpe8 <- getcountmatrix(matrixhumanchorrpe8)
countmatrixhumanchorrpe9 <- getcountmatrix(matrixhumanchorrpe9)
countmatrixhumanchorrpe10 <- getcountmatrix(matrixhumanchorrpe10)
countmatrixhumanchorrpe11 <- getcountmatrix(matrixhumanchorrpe11)
countmatrixhumanchorrpe12 <- getcountmatrix(matrixhumanchorrpe12)
countmatrixhumanchorrpe13 <- getcountmatrix(matrixhumanchorrpe13)
countmatrixhumanchorrpe14 <- getcountmatrix(matrixhumanchorrpe14)


countmatriceslist<-list(countmatrixhumanchorrpe1,countmatrixhumanchorrpe2,countmatrixhumanchorrpe3,countmatrixhumanchorrpe4,countmatrixhumanchorrpe5,countmatrixhumanchorrpe6,countmatrixhumanchorrpe7,countmatrixhumanchorrpe8,countmatrixhumanchorrpe9,countmatrixhumanchorrpe10,countmatrixhumanchorrpe11,countmatrixhumanchorrpe12,countmatrixhumanchorrpe13,countmatrixhumanchorrpe14)

for(i in 1:length(countmatriceslist)){
  row.names(countmatriceslist[[i]])<-processmatrixlist[[i]][,1]
}

for(i in 1:length(countmatriceslist)){
  countmatriceslist[[i]]<-t(countmatriceslist[[i]])
}

for(i in 1:length(countmatriceslist)){
  countmatriceslist[[i]]<-CreateSeuratObject(countmatriceslist[[i]], min.cells=3, min.features=200)
}

for (i in 1:length(countmatriceslist)) {
  countmatriceslist[[i]] <- SCTransform(countmatriceslist[[i]], verbose = FALSE)
}

chorrpe.features <- SelectIntegrationFeatures(object.list = countmatriceslist, nfeatures = 3000)
countmatriceslist <- PrepSCTIntegration(object.list = countmatriceslist, anchor.features = chorrpe.features, 
                                        verbose = FALSE)

for(i in 1:14){
  countmatriceslist[[i]][["percent.mt"]]<-PercentageFeatureSet(countmatriceslist[[i]], pattern = "^MT-", assay="RNA")
}

for(i in 1:14){
  countmatriceslist[[i]] <- subset(countmatriceslist[[i]], subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 5)
}


chorrpe.anchors <- FindIntegrationAnchors(object.list = countmatriceslist, normalization.method = "SCT", 
                                          anchor.features = chorrpe.features, verbose = FALSE)
chorrpe.integrated <- IntegrateData(anchorset = chorrpe.anchors, normalization.method = "SCT", 
                                    verbose = FALSE)

chorrpe.integrated <- RunPCA(chorrpe.integrated, verbose = FALSE)
chorrpe.integrated <- RunUMAP(chorrpe.integrated, dims = 1:30)
chorrpe.integrated <- RunTSNE(chorrpe.integrated, dims = 1:30, check_duplicates = FALSE)

### attach metadata
Metaprocessmatrixlistnew<-rbind(Metaprocessmatrixlist[[1]], Metaprocessmatrixlist[[2]], Metaprocessmatrixlist[[3]], Metaprocessmatrixlist[[4]], Metaprocessmatrixlist[[5]], Metaprocessmatrixlist[[6]], Metaprocessmatrixlist[[7]], Metaprocessmatrixlist[[8]], Metaprocessmatrixlist[[9]], Metaprocessmatrixlist[[10]], Metaprocessmatrixlist[[11]], Metaprocessmatrixlist[[12]], Metaprocessmatrixlist[[13]], Metaprocessmatrixlist[[14]])
row.names(Metaprocessmatrixlistnew)<-Metaprocessmatrixlistnew[,1]
Metaprocessmatrixlistnew<-Metaprocessmatrixlistnew[,-1]

chorrpe.integrated$cellname<-chorrpe.integrated@active.ident
chorrpe.integrated$cellname<-as.character(chorrpe.integrated$cellname)
chorrpe.integrated$cellname<-colnames(chorrpe.integrated[["RNA"]]@counts)

chorrpe.integrated$experiment<-chorrpe.integrated$cellname
chorrpe.integrated$voigtcluster<-chorrpe.integrated$cellname
chorrpe.integrated$perimac<-chorrpe.integrated$cellname
chorrpe.integrated$donor<-chorrpe.integrated$cellname

for(i in 1:length(colnames(chorrpe.integrated[["RNA"]]@counts))){
  chorrpe.integrated$voigtcluster[i]<-as.character(Metaprocessmatrixlistnew[colnames(chorrpe.integrated[["RNA"]]@counts)[i],1])
  chorrpe.integrated$perimac[i]<-as.character(Metaprocessmatrixlistnew[colnames(chorrpe.integrated[["RNA"]]@counts)[i],2])
  chorrpe.integrated$donor[i]<-as.character(Metaprocessmatrixlistnew[colnames(chorrpe.integrated[["RNA"]]@counts)[i],2])
}

chorrpe.integrated$experiment<-chorrpe.integrated$donor

chorrpe.integrated$experiment<-as.factor(chorrpe.integrated$experiment)

levels(chorrpe.integrated$experiment)<-c("exp1", "exp2", "exp3", "exp4", "exp5", "exp6", "exp7", "exp8", "exp9", "exp10", "exp11", "exp12","exp13", "exp14")

chorrpe.integrated$voigtcluster<-as.factor(chorrpe.integrated$voigtcluster)

chorrpe.integrated$perimac<-as.factor(chorrpe.integrated$perimac)
levels(chorrpe.integrated$perimac)[1:7]<-"Macula"
levels(chorrpe.integrated$perimac)[2:8]<-"Periphery"

chorrpe.integrated$donor<-as.factor(chorrpe.integrated$donor)
levels(chorrpe.integrated$donor)[1:7]<-c("pat1", "pat2", "pat3", "pat4", "pat5", "pat6", "pat7")
levels(chorrpe.integrated$donor)[8:14]<-c("pat1", "pat2", "pat3", "pat4", "pat5", "pat6", "pat7")



chorrpe.integrated <- FindNeighbors(chorrpe.integrated, dims = 1:15)
chorrpe.integrated <- FindClusters(chorrpe.integrated, resolution = 4)


DefaultAssay(chorrpe.integrated) <- "RNA"
chorrpe.integrated <- NormalizeData(chorrpe.integrated)
plan("multiprocess", workers = 5)
chorrpe.integrated.markers <- FindAllMarkers(chorrpe.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

seuratannnotations<-c(
  "Macroph.", #
  "Endot.", #
  "Endot.", # 
  "Fibrobl.", #
  "Endot.", #
  "delete", #
  "Endot.", #
  "Endot.", #
  "T_NKCell", #
  "Endot.", #
  "Endot.", #
  "delete", #
  "n.d.", #
  "Macroph.", #
  "Fibrobl.", #
  "Endot.", #
  "Endot.", #
  "n.d.", #
  "Endot.", #
  "Endot.", #
  "Macroph.", #
  "Endot.", #
  "T_NKCell", #
  "Macroph.", #
  "Macroph.", #
  "BCell", #
  "Endot.", #
  "Mural", #
  "Fibrobl.", #
  "Endot.", #
  "Fibrobl.",
  "Fibrobl.", #
  "Fibrobl.", #
  "Fibrobl.", #
  "BCell", #
  "Endot.", #
  "Schwann", #
  "Macroph.", #
  "Mural", #
  "Melanoc.", #
  "T_NKCell", #
  "n.d.", #
  "Schwann", #
  "delete", #
  "Endot.", #
  "delete", #
  "RPE", #
  "Fibrobl.", #
  "Macroph.", #
  "n.d.", #
  "Endot.", # 
  "Mural", # 
  "delete", #
  "Macroph.", #
  "Melanoc.", #
  "Melanoc.", #
  "Fibrobl.", # 
  "Endot.", #
  "Endot.", #
  "Endot.", #
  "T_NKCell", #
  "Melanoc.", #
  "delete", #
  "T_NKCell", #
  "Endot.", # 
  "delete", #
  "n.d.", #
  "delete", #
  "RPE" #
)

chorrpe.integrated$seuratannnotations <- chorrpe.integrated$seurat_clusters
levels(chorrpe.integrated$seuratannnotations) <- seuratannnotations

chorrpe.integrated@active.ident<-chorrpe.integrated$seuratannnotations
chorrpe.integrated <-subset(chorrpe.integrated, ident="delete", invert = TRUE)




#Visualization:
##functions DimPlot(), FeaturePlot(), Heatmap() were used with relevant metadata presentation.
##for ARI indexing we used the ARI() function with relevant metadata factors as input variables.






