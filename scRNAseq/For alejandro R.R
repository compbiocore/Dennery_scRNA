library(Seurat)
library(cowplot)
library(ggplot2)
set.seed(61)

#must make rference to identify gene names on OSCAR before hand
#raw reads from Genewiz were downloaded on OSCAR
#cellranger command were used to filter the data and assign gene names to each gene 
#this is an overnight run!
#as well as transcription and number of cells 

#load data into Rstudio
AirD7.data <- Read10X(data.dir ="./AirD7out/filtered_feature_bc_matrix/")
O2D7.data <- Read10X(data.dir ="./O2D7out/filtered_feature_bc_matrix/")
SD7.data <- Read10X(data.dir ="./SD7out/filtered_feature_bc_matrix/")

#changing the data into a Seurat readable Matrix
#ignore cells with less than 200 genes as default
#at least 3 cells per cluster maybe
AirD7.data <- CreateSeuratObject(counts = AirD7.data,  project = "AirD7",min.cells = 3, min.features = 200)
O2D7.data <- CreateSeuratObject(counts = O2D7.data,  project = "O2D7",min.cells = 3, min.features = 200)
SD7.data <- CreateSeuratObject(counts = SD7.data,  project = "SD7",min.cells = 3, min.features = 200)

#name of each condition for the data. Metadata is information about the dataseet.
#Con, paste AirD7 on the condition AirD7
AirD7.data@meta.data$Con <- paste("AirD7")
O2D7.data@meta.data$Con <- paste("O2D7")
SD7.data@meta.data$Con <- paste("SD7")


#Test nCount and nFeatures:
#nCount is a quality check - need to define parameters for next steps.
#number of transcript reads per cell for each dataset
e <-VlnPlot(AirD7.data, features = c("nCount_RNA"), pt.size = 0.5, y.max = 50000) + labs(title = "AirD7") +theme(axis.title.x=element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),legend.position = 'none')   
f <-VlnPlot(O2D7.data, features = c("nCount_RNA"), pt.size = 0.5, y.max = 50000) + labs(title = "O2D7") +theme(axis.title.x=element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),legend.position = 'none')
g <-VlnPlot(SD7.data, features = c("nCount_RNA"), pt.size = 0.5, y.max = 50000) + labs(title = "SD7") +theme(axis.title.x=element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),legend.position = 'none')
plot_grid(e,f,g, nrow = 1, ncol = 3)

#nFeature is another quality check
#difference is the end picture represents the number of genes detected per cell
#We will only take the parameters that represent most of the cells
g <-VlnPlot(AirD7.data, features = c("nFeature_RNA"), pt.size = 0.5, y.max = 10000) + labs(title = "AirD7") +theme(axis.title.x=element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),legend.position = 'none')   
h <-VlnPlot(O2D7.data, features = c("nFeature_RNA"), pt.size = 0.5, y.max = 10000) + labs(title = "O2D7") +theme(axis.title.x=element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),legend.position = 'none')
i <-VlnPlot(SD7.data, features = c("nFeature_RNA"), pt.size = 0.5, y.max = 10000) + labs(title = "SD7") +theme(axis.title.x=element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),legend.position = 'none')
plot_grid(g,h,i, nrow = 1, ncol = 3)


#remove cells with less than 200 genes per cells
#similarly cells with more than 7000 (arbitrary cutoff at graph, most cells <7000)
AirD7.data <- subset(AirD7.data, subset = nFeature_RNA > 200 & nFeature_RNA < 7000)
O2D7.data <- subset(O2D7.data, subset = nFeature_RNA > 200 & nFeature_RNA < 7000)
SD7.data <- subset(SD7.data, subset = nFeature_RNA > 200 & nFeature_RNA < 7000)

#Normalize each dataset separately before comparing
#In order to have expresssion level per dataset 
AirD7.data <- NormalizeData(AirD7.data, verbose = FALSE)
O2D7.data <- NormalizeData(O2D7.data, verbose = FALSE)
SD7.data <- NormalizeData(SD7.data, verbose = FALSE)

#look for the transcript variation in each dataset
AirD7.data <- FindVariableFeatures(AirD7.data, selection.method = "vst", nfeatures = 7000, verbose = FALSE)
O2D7.data <- FindVariableFeatures(O2D7.data, selection.method = "vst", nfeatures = 7000, verbose = FALSE)
SD7.data <- FindVariableFeatures(SD7.data, selection.method = "vst", nfeatures = 7000, verbose = FALSE)

#integrate the data together with FindIntegrationAnchors making a list for each of the datasets
#smaller dimension the more stringent
Anchors <- FindIntegrationAnchors(object.list = list(AirD7.data, O2D7.data, SD7.data), dims = 1:20)
#RNA expression and integration dataset. To make clusters you need integrated.
#For futureplot you go back to the RNA
CombinedAleD7 <- IntegrateData(anchorset = Anchors, dims = 1:20)
DefaultAssay(CombinedAleD7) <- "integrated"

#add Scale data vars to regress (not in the initial script)
#vars.to.regress is a stats function
#looks at number of genes and transcripts per cell and normalize the data in the integrated dataset.
#normalize the expression to have a similar comparison accross the datasets
CombinedAleD7 <- ScaleData(CombinedAleD7, vars.to.regress = c("nFeature_RNA","nCount_RNA"))
#principal component analysis (PCA) - 3D graph x,y,z used to make tisney plot defining which cells are closer to each order
CombinedAleD7 <- RunPCA(CombinedAleD7, npcs = 30, verbose = FALSE)
plot<-DimPlot(CombinedAleD7, reduction = "pca", pt.size = 0.5)
#vaguely gives the elbow - look at background after the elbow. At 20 it's less stringent than at 15.
#it's a bit of hard to tell
ElbowPlot(CombinedAleD7)

DimHeatmap(CombinedAleD7, dims = 1:20, cells = 500, balanced = TRUE)

# t-SNE and Clustering
CombinedAleD7 <- RunTSNE(CombinedAleD7, dims = 1:15)
CombinedAleD7 <- FindNeighbors(CombinedAleD7, dims = 1:15)
CombinedAleD7 <- FindClusters(CombinedAleD7, resolution = 0.5)
DimPlot(CombinedAleD7, reduction = "tsne", group.by = "Con", pt.size = 0.5)
DimPlot(CombinedAleD7, reduction = "tsne", split.by = "Con",label = TRUE)
#look up in seurat how to chang the name of each clusters
DimPlot(CombinedAleD7, reduction = "tsne", split.by = "Con")

#Features plot for specific genes
DefaultAssay(CombinedAleD7) <- "RNA"
FeaturePlot(CombinedAleD7, features = c("Hopx"))
FeaturePlot(CombinedAleD7, features = c("Adgre1", "Vim", "Vegfa", "Fn1", "Acta2", "Epcam"), split.by = "Con", col=c("grey", "red"))
FeaturePlot(CombinedAleD7, features = c("Adar"), split.by = "Con", col=c("grey", "red"))
FeaturePlot(CombinedAleD7, features = c("Hopx"), split.by = "Con")

#create combined vlnplot
#vilin plot for specific genes
#Gene ezpression per cluster 
plot <-VlnPlot(CombinedAleD7, features = c("Hmox1"), split.by = "Con",  pt.size = 0, combine = FALSE)
CombinePlots(plot, ncol = 1)

#this step is important. Save everything previously done until here. Just open the RDS file and you will have all the info here.
#run libraries and then open the RDS file
saveRDS(CombinedAleD7,file=paste("CombinedAleD7Dim15_res0.5.rds",sep=""))


#Subset per condition:
#subset air separately for exmaple
Air <- subset(CombinedAleD7, subset= Con=='AirD7')
O2 <- subset(CombinedAleD7, subset= Con=='O2D7')
S <- subset(CombinedAleD7, subset= Con=='SD7')

#only look at air
VlnPlot(Air, features = c("Hmox1"), pt.size = 0.5) + 
  theme(axis.text.x = element_text(angle=0),legend.position = 'none') + labs(title = "Air Hmox1 Expression")

#only look at O2
VlnPlot(O2, features = c("Hmox1"), pt.size = 0.5) + 
  theme(axis.text.x = element_text(angle=0),legend.position = 'none') + labs(title = "O2 Hmox1 Expression")

VlnPlot(S, features = c("Hmox1"), pt.size = 0.5) + 
  theme(axis.text.x = element_text(angle=0),legend.position = 'none') + labs(title = "S Hmox1 Expression")



#compare one to one at each time separetly
#Stephany said to use this one: Change the number 1 by 1

Air7 <- subset(CombinedAleD7, subset= Con=='AirD7')

#excel file of conserved markers
#ident.1 in cluster 1
con.markers <- FindConservedMarkers(Air7, ident.1 = 1,logfc.threshold = 0.25, grouping.var = "Con", verbose = FALSE)
write.table(con.markers,file=paste("./ConservednewMarkersD7/","CONSERVEDAircluster1.csv",sep=""),quote=F,sep="\t")

con.markers <- FindConservedMarkers(CombinedAleD7, ident.1 = 1,logfc.threshold = 0.25, grouping.var = "Con", verbose = FALSE)
write.table(con.markers,file=paste("./ConservednewMarkersD7/","CONSERVEDALEcluster12.csv",sep=""),quote=F,sep="\t")
#find markers find more unique markerz. Doesn't care if it is present in other conditions
#find conserved markers

#Stephany Foster is the best: Differential gene expression per cluster per condition
output <- "./DifferentialMarkersD7onAirandS/"
CombinedAleD7$celltype.condition <- paste(Idents(CombinedAleD7), CombinedAleD7$Con, sep="_")
CombinedAleD7$celltype <- Idents(CombinedAleD7)
Idents(CombinedAleD7) <- "celltype.condition"

for (i in 0:22){ #or however many clusters you have, we have 22
  try({
    ident1 <- paste0(i,"_AirD7")
    ident2 <- paste0(i,"_SD7")
    condition.diffgenes <- FindMarkers(CombinedAleD7, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.5)
    write.table(condition.diffgenes, file=paste0(output,i,".txt"))
  })
}

#Maybe try this (also work)
#end of analysis, this is optional
Idents(CombinedAle) <- CombinedAle$seurat_clusters
DefaultAssay(CombinedAle) <- "RNA"
CombinedAle.cluster0 <- FindMarkers(CombinedAle, ident.1 = "AirD7", ident.2 = "O2D7", verbose = TRUE, group.by="Con", subset.ident = "0", logfc.threshold = 0.25)
write.table(CombinedAle.cluster0, file=paste("DifferentialMarkersD7onAirandO2/Cluster0.txt",sep=""),quote=F,sep="\t")


#give number of cells per cluster per condition:
sweet <- print((table(CombinedAleD7$celltype.condition)))
write.table(sweet, file= "CombinedD7cells_per_cluster.csv")



#give number of cells per condition:
print((table(CombinedAleD7$Con)))

saveRDS(CombinedAleD7,file=paste("CombinedAleD7B_res0.5.rds",sep=""))






