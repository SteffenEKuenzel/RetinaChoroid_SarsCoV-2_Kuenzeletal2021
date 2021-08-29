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
library(MixGHD)
options(future.globals.maxSize = 8000 * 1024^2)

#read in the annotation file from MENONetal2019
fileMENONannot <- read.table("###/sample_annotationfile.tsv")

#get individ. cell names
CELLNAMES <- fileMENONannot[,1]
CELLNAMES <- CELLNAMES[-1]

#read in and prepare countmatrix from MENONetal2019 and create Seurat object (totalMENON_SEURAT) with given parameters
fileMENON_MATRIX <- readMM("###/countmatrix.mtx")
fileMENON_GENENAMES <- read.csv("###/gene_names.txt", header=FALSE)
fileMENON_GENENAMES <- fileMENON_GENENAMES[,1]
row.names(fileMENON_MATRIX)<-fileMENON_GENENAMES
colnames(fileMENON_MATRIX)<-CELLNAMES
totalMENON_SEURAT <- CreateSeuratObject(fileMENON_MATRIX, min.cells = 3, min.features = 100, project = "HUMAN_RETINA_COVID19")


#extract location(maculavsperiphery) and number of retina (donor1,2,3) and batch. put it in own objects for metadata annotation
RETandLOC<- fileMENONannot[,2]
RETandLOC<- RETandLOC[-1]

Retina<-RETandLOC
levels(Retina)<-c("Retina1","Retina2","Retina3","Retina1","Retina2","Retina3","Retina3")
Retina<-as.character(Retina)
names(Retina)<-CELLNAMES

LOC <- RETandLOC
levels(LOC)<-c("Macula","Macula","Macula","Periphery","Periphery","Periphery","Periphery")
LOC<-as.character(LOC)
names(LOC)<-CELLNAMES

Batch <- RETandLOC
Batch <- as.character(Batch)
names(Batch)<-CELLNAMES

totalMENON_SEURAT$LOC<-totalMENON_SEURAT@active.ident
totalMENON_SEURAT$LOC<-as.character(totalMENON_SEURAT$LOC)
totalMENON_SEURAT$LOC<-"NULL"
totalMENON_SEURAT$Retina<-totalMENON_SEURAT@active.ident
totalMENON_SEURAT$Retina<-"NULL"
totalMENON_SEURAT$Retina<-as.character(totalMENON_SEURAT$Retina)
totalMENON_SEURAT$Batch<-totalMENON_SEURAT@active.ident
totalMENON_SEURAT$Batch<-"NULL"
totalMENON_SEURAT$Batch<-as.character(totalMENON_SEURAT$Batch)
for(i in 1:length(totalMENON_SEURAT$LOC)){
  totalMENON_SEURAT$LOC[i]<- LOC[names(totalMENON_SEURAT$LOC[i])]
}
totalMENON_SEURAT$LOC<-as.factor(totalMENON_SEURAT$LOC)
for(i in 1:length(totalMENON_SEURAT$Retina)){
  totalMENON_SEURAT$Retina[i]<- Retina[names(totalMENON_SEURAT$Retina[i])]
}
totalMENON_SEURAT$Retina<-as.factor(totalMENON_SEURAT$Retina)
for(i in 1:length(totalMENON_SEURAT$Batch)){
  totalMENON_SEURAT$Batch[i]<- Batch[names(totalMENON_SEURAT$Batch[i])]
}
totalMENON_SEURAT$Batch<-as.factor(totalMENON_SEURAT$Batch)

# ORIGPUBLIC means the cell type annotation of the original publication, which will later be used for ARI indexing
preclusterORIGPUBLIC <- fileMENONannot[,c(1,49)]
row.names(preclusterORIGPUBLIC)<-fileMENONannot[,1]
preclusterORIGPUBLIC[,2]<-as.character(preclusterORIGPUBLIC[,2])

totalMENON_SEURAT$precluster<-totalMENON_SEURAT@active.ident
totalMENON_SEURAT$precluster<-"None"

for(i in 1:length(totalMENON_SEURAT$precluster)){
  totalMENON_SEURAT$precluster[i]<- preclusterORIGPUBLIC[names(totalMENON_SEURAT$precluster[i]),2]
}
totalMENON_SEURAT$precluster<-as.factor(totalMENON_SEURAT$precluster)


# compare Seurat workflow. extract cells that have a high percentage of mt rna, then follow seurat sctransform workflow
totalMENON_SEURAT[["percent.mt"]] <- PercentageFeatureSet(totalMENON_SEURAT, pattern = "MT")
totalMENON_SEURAT_reduced <- subset(totalMENON_SEURAT, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 10)
totalMENON_SEURAT_reduced <- SplitObject(totalMENON_SEURAT_reduced, split.by="Batch")
for(i in 1:length(totalMENON_SEURAT_reduced)){totalMENON_SEURAT_reduced[[i]] <- SCTransform(totalMENON_SEURAT_reduced[[i]],verbose=FALSE)} 
totalMENON_SEURAT_reduced.features <- SelectIntegrationFeatures(object.list = totalMENON_SEURAT_reduced, nfeatures = 3000)
totalMENON_SEURAT_reduced <- PrepSCTIntegration(object.list = totalMENON_SEURAT_reduced, anchor.features = totalMENON_SEURAT_reduced.features, verbose = FALSE)
totalMENON_SEURAT_reduced.anchors <- FindIntegrationAnchors(object.list = totalMENON_SEURAT_reduced, normalization.method = "SCT", anchor.features = totalMENON_SEURAT_reduced.features, verbose = FALSE)
totalMENON_SEURAT_reduced.INTEGRATED <- IntegrateData(totalMENON_SEURAT_reduced.anchors, normalization.method = "SCT",  verbose = FALSE)
totalMENON_SEURAT_reduced.INTEGRATED <- RunPCA(totalMENON_SEURAT_reduced.INTEGRATED, verbose = FALSE)
totalMENON_SEURAT_reduced.INTEGRATED <- RunUMAP(totalMENON_SEURAT_reduced.INTEGRATED, dims = 1:50)
totalMENON_SEURAT_reduced.INTEGRATED <- RunTSNE(totalMENON_SEURAT_reduced.INTEGRATED)

totalMENON_SEURAT_reduced.INTEGRATED <- FindNeighbors(totalMENON_SEURAT_reduced.INTEGRATED, dims = 1:10)
totalMENON_SEURAT_reduced.INTEGRATED <- FindClusters(totalMENON_SEURAT_reduced.INTEGRATED, resolution = 1.8)


DefaultAssay(totalMENON_SEURAT_reduced.INTEGRATED) <- "RNA"
totalMENON_SEURAT_reduced.INTEGRATED <- NormalizeData(totalMENON_SEURAT_reduced.INTEGRATED)
plan("multiprocess", workers = 4)
totalMENON_SEURAT_reduced.INTEGRATED.markers <- FindAllMarkers(totalMENON_SEURAT_reduced.INTEGRATED, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


#celltype annotation was performed by consulting the cell type markers as indicated in suppl. figure. Suppl. plots were generated by the DotPlot() function

#cell type annotation and removal of unclear small subclusters (mixed cell types and doublets)
annotation_clusterlevels<-
  c(
    "Rod PR",
    "Rod PR",
    "Rod PR",
    "Rod PR",
    "Rod PR",
    "Ganglion C.",
    "Fibroblast",
    "Rod PR",
    "Fibroblast",
    "Bipolar C.",
    "Rod PR",
    "Rod PR",
    "Fibroblast",
    "Macroglia",
    "Rod PR",
    "Bipolar C.",
    "Bipolar C.",
    "Macroglia",
    "Macroglia",
    "Fibroblast",
    "Rod PR",
    "Macroglia",
    "Macroglia",
    "Rod PR",
    "Amacrine C.",
    "Remove",
    "Cone PR",
    "Amacrine C.",
    "Macroglia",
    "Horizontal C.",
    "Bipolar C.",
    "Vasc. Cells + Microgl.",
    "Remove",
    "Rod PR",
    "Bipolar C.")


totalMENON_SEURAT_reduced.INTEGRATED <-subset(totalMENON_SEURAT_reduced.INTEGRATED, ident=c("25","32"), invert = TRUE)


totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated <- totalMENON_SEURAT_reduced.INTEGRATED$seurat_clusters
levels(totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated) <- annotation_clusterlevels
totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated_neworder <-totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated


totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated<-as.factor(totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated)



#in order to check whether a cell expresses at least one copy of a certain gene, we used code like this (example is for FURIN)
totalMENON_SEURAT_reduced.INTEGRATED$positiveforFURIN<-totalMENON_SEURAT_reduced.INTEGRATED$Batch
totalMENON_SEURAT_reduced.INTEGRATED$positiveforFURIN<-as.factor(totalMENON_SEURAT_reduced.INTEGRATED$positiveforFURIN)
totalMENON_SEURAT_reduced.INTEGRATED$positiveforFURIN<-"NULL"

for(i in 1:length(totalMENON_SEURAT_reduced.INTEGRATED$positiveforFURIN)){
  if(GetAssayData(object = totalMENON_SEURAT_reduced.INTEGRATED[["RNA"]], slot = "counts")["FURIN",i]>0){
    totalMENON_SEURAT_reduced.INTEGRATED$positiveforFURIN[i]<-"POS"
  }
}


#exclude the "Remove" level
totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated<-as.character(totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated)
totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated<-as.factor(totalMENON_SEURAT_reduced.INTEGRATED$seurat_cluster_annotated)




####the code from here can similarly be used for the choroid dataset (adaption of the seuratobject-name)


#in order to get CPM values (per cell type per batch), we used the following code (seuratobj is a general object, that can later be used for different seratobjects, e.g. the choroid dataset)
seuratobj<-totalMENON_SEURAT_reduced.INTEGRATED
seuratobj$batch<-seuratobj$Batch
seuratobj$annot<-seuratobj$seurat_cluster_annotated

seuratobj$batch<-as.factor(seuratobj$batch)
seuratobj$annot<-as.factor(seuratobj$annot)
seuratobj@active.ident<-seuratobj$batch

# creation of list that: each object of the list is a seurat object that contains cells of one batch
myls <- vector("list", length = length(levels(seuratobj$batch)))
for(i in 1:length(levels(seuratobj$batch))){
  myls[[i]]<-subset(seuratobj, ident=levels(seuratobj$batch)[[i]])
}

# creation of a matrix that will contain the number of cells per cell type (cols) and per batch (rows)
genematrix<-matrix(0, nrow = length(levels(seuratobj$batch)), ncol = length(levels(seuratobj$annot)))
row.names(genematrix)<-c(levels(seuratobj$batch))
colnames(genematrix)<-c(levels(seuratobj$annot))
for(i in 1:length(levels(seuratobj$batch))){
  genematrix[i,]<-as.vector(table(myls[[i]]$annot))
}

#if min(genematrix) is equal zero, the following code has to be modified
min(genematrix)  

#creation of a list of length of the number of batches. each object in the list has the length of the number of cell types.
mylsfinal<-myls
for(i in 1:length(levels(seuratobj$batch))){
  mylsfinal[[i]]<-vector("list", length = length(levels(seuratobj$annot)))
}

#the list mylsfinal should have the length of the number of batches. each object of this list is again a list has the length of the number of annotated cell types - so each of the smaller objects contains one cell type of a certain batch
for(i in 1:length(levels(seuratobj$batch))){
  for(k in 1:length(levels(seuratobj$annot))){
    myls[[i]]@active.ident<- myls[[i]]$annot
    mylsfinal[[i]][[k]]<-subset(myls[[i]], ident=(levels(seuratobj$annot)[[k]]))
  }
}

#create a genematrix that will give out the values for analysis
genematrix<-matrix(0, nrow = length(levels(seuratobj$batch)), ncol = length(levels(seuratobj$annot)))
row.names(genematrix)<-c(levels(seuratobj$batch))
colnames(genematrix)<-c(levels(seuratobj$annot))


#to check the CPM for a certain gene, just give in the given gene name in here (FURIN as example)
genename<-"FURIN"
for(i in 1:length(levels(seuratobj$batch))){
  for(k in 1:length(levels(seuratobj$annot))){
    genematrix[i,k]<- sum(GetAssayData(object = mylsfinal[[i]][[k]][["RNA"]], slot = "counts")[genename,]) / sum(GetAssayData(object = mylsfinal[[i]][[k]][["RNA"]], slot = "counts")[,]) *1000000
  } 
}


#to get a matrix in which you can see which celltypes in which batch are positive for a certain gene, use this code
#$positiveforFURIN has to be created before (see above)
positivecellsmatrix<-matrix(0, nrow = length(levels(seuratobj$batch)), ncol = length(levels(seuratobj$annot)))
row.names(positivecellsmatrix)<-c(levels(seuratobj$batch))
colnames(positivecellsmatrix)<-c(levels(seuratobj$annot))

for(i in 1:length(levels(seuratobj$batch))){
  for(k in 1:length(levels(seuratobj$annot))){
    positivecellsmatrix[i,k]<- 100*table(mylsfinal[[i]][[k]]$positiveforFURIN)[2]/(table(mylsfinal[[i]][[k]]$positiveforFURIN)[1]+table(mylsfinal[[i]][[k]]$positiveforFURIN)[2])
  } 
}

#Visualization:
##functions DimPlot(), FeaturePlot(), Heatmap() were used with relevant metadata presentation.
##for ARI indexing we used the ARI() function with relevant metadata factors as input variables.










