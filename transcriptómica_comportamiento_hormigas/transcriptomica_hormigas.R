library(DESeq2)

Counts = read.delim("Raw_counts.txt", sep="\t", header=TRUE) #importo la tabla de counts
# Transformo en una matriz la tabla de counts porque así lo necesita el Deseq2
rnames <- Counts [,1] # Columna 1 son los nombres
Counts_matrix <-data.matrix(Counts[,2:ncol(Counts)]) #Transformar en matrix ap columna2
rownames(Counts_matrix) <- rnames #Asignar nombre a la matrix

#Cargo diseño experimental
experiment = read.delim("experimental_design.txt", sep="\t", header=TRUE)

dds <- DESeqDataSetFromMatrix(countData = Counts_matrix,
                                      colData = experiment,
                                      design= ~Treatment) #design es el nombre de la segunda columna del archivo experimento

#Filtrado de genes poco expresados y/o muy variables
library(HTSFilter)

dds_filtered <- HTSFilter(dds, s.len=100, plot=TRUE)$filteredData #s.len=n modelizaciones
dim(dds_filtered) #Pasamos de 13084 a 9116

#Generamos los valores transformados que son útiles para hacer heatmaps de genes diferencialmente expresados
rlog<- rlogTransformation(dds_filtered, blind = TRUE, fitType ="parametric")
tab<-data.frame(assay(rlog))
write.table(tab, file="RLOG_Transformed_values.txt")

#Exploramos los datos
library("pheatmap")
library("RColorBrewer")

sampleDists <- dist(t(assay(rlog)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rlog$sample)
rownames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

tiff("Heatmap_RLD.tif",           
     width = 8*100,        # 5 x 100 pixels
     height = 8*100,
     res = 100,        # 100 pixels per inch
     pointsize = 8) 
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

#PCA plot 

library("ggpubr")
tiff("PCA_plotRLD_sintags.tif",           
     width = 6*200,        # 5 x 300 pixels
     height = 6*200,
     res = 100,            # 300 pixels per inch
     pointsize = 8)        #  font size
plotPCA(rlog, intgroup = c("Treatment"))
dev.off()

#Para generar plot con tags de las condiciones correr y salvar a partir de hacer zoom
p<-plotPCA(rlog, intgroup = c("Treatment")) + geom_text(aes(label=name), size=3, check_overlap = TRUE, vjust= -1, hjust=1)

#Analisis de expresion diferencial
dds_analysis <- DESeq(dds_filtered)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

resultsNames(dds_analysis)
#[1]"Intercept"                "Treatment_Treated_vs_Control"            

res <- results(dds_analysis, contrast=c("Treatment","Treated","Control"), independentFiltering = FALSE) #Fold change positivo  indica aumento en macho respecto a hembra
resLFC <- lfcShrink(dds_analysis, coef = 2, res = res, type = c("apeglm"))
write.table(resLFC, file="Deseq2_statistics.txt")

summary(resLFC)

#out of 9116 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 17, 0.19%
#LFC < 0 (down)     : 35, 0.38%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


tiff("MAplot_dds.tif",           
     width = 8*300,        # 5 x 300 pixels
     height = 8*300,
     res = 300,            # 300 pixels per inch
     pointsize = 11)
plotMA(resLFC, ylim=c(-4,4), main = "Male vs Female", cex = 0.5) 
abline(h=c(-1,1), col="blue") 
dev.off()

===
#Solo considerando 3 x 3, C1,2 y C5 del control y T3,T4 y T5 del tratado
  
library(DESeq2)

Counts = read.delim("Raw_counts3x3.txt", sep="\t", header=TRUE) #importo la tabla de counts
# Transformo en una matriz la tabla de counts porque así lo necesita el Deseq2
rnames <- Counts [,1] # Columna 1 son los nombres
Counts_matrix <-data.matrix(Counts[,2:ncol(Counts)]) #Transformar en matrix ap columna2
rownames(Counts_matrix) <- rnames #Asignar nombre a la matrix

#Cargo diseño experimental
experiment = read.delim("experimental_design3x3.txt", sep="\t", header=TRUE)

dds <- DESeqDataSetFromMatrix(countData = Counts_matrix,
                              colData = experiment,
                              design= ~Treatment) #design es el nombre de la segunda columna del archivo experimento

#Filtrado de genes poco expresados y/o muy variables
library(HTSFilter)

dds_filtered <- HTSFilter(dds, s.len=100, plot=TRUE)$filteredData #s.len=n modelizaciones
dim(dds_filtered) #Pasamos de 13084 a 8921

#Generamos los valores transformados que son útiles para hacer heatmaps de genes diferencialmente expresados
rlog<- rlogTransformation(dds_filtered, blind = TRUE, fitType ="parametric")
tab<-data.frame(assay(rlog))
write.table(tab, file="RLOG_Transformed_values.txt")

#Exploramos los datos
library("pheatmap")
library("RColorBrewer")

sampleDists <- dist(t(assay(rlog)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rlog$sample)
rownames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

tiff("Heatmap_RLD3x3.tif",           
     width = 8*100,        # 5 x 100 pixels
     height = 8*100,
     res = 100,        # 100 pixels per inch
     pointsize = 8) 
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

#No se pq no corre a la primera y hay que correrlo de vuelta para que se genere el resultado
#PCA plot 

library("ggpubr")
tiff("PCA_plotRLD_sintags3X3.tif",           
     width = 6*200,        # 5 x 300 pixels
     height = 6*200,
     res = 100,            # 300 pixels per inch
     pointsize = 8)        #  font size
plotPCA(rlog, intgroup = c("Treatment"))
dev.off()

#Para generar plot con tags de las condiciones correr y salvar a partir de hacer zoom
p<-plotPCA(rlog, intgroup = c("Treatment")) + geom_text(aes(label=name), size=3, check_overlap = TRUE, vjust= -0.5, hjust=1)

#Analisis de expresion diferencial
dds_analysis <- DESeq(dds_filtered)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

resultsNames(dds_analysis)
#[1]"Intercept"                "Treatment_Treated_vs_Control"            

res <- results(dds_analysis, contrast=c("Treatment","Treated","Control"), independentFiltering = FALSE) #Fold change positivo  indica aumento en macho respecto a hembra
resLFC <- lfcShrink(dds_analysis, coef = 2, res = res, type = c("apeglm"))
write.table(resLFC, file="Deseq2_statistics3x3.txt")

summary(resLFC)

#out of 8921 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 743, 8.3%
#LFC < 0 (down)     : 1404, 16%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 0)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


tiff("MAplot_3x3.tif",           
     width = 8*300,        # 5 x 300 pixels
     height = 8*300,
     res = 300,            # 300 pixels per inch
     pointsize = 11)
plotMA(resLFC, ylim=c(-4,4), main = "1% sucrose diet vs 20% sucrose diet", cex = 0.5) 
abline(h=c(-1,1), col="blue") 
dev.off()

====
#https://towardsdatascience.com/pheatmap-draws-pretty-heatmaps-483dab9a3cc
library(pheatmap)
library(gplots)
library(grid)

Rlog_DEGs= read.delim("Rlog_DE_genes.txt", sep="\t", header=TRUE)
rnames <- Rlog_DEGs [,1] # Columna 1 son lons nombres
Rlog_DEGs <-data.matrix(Rlog_DEGs[,2:ncol(Rlog_DEGs)]) #Transformar en matrix ap columna2
rownames(Rlog_DEGs) <- rnames #Asignar nombre a la matrix
tiff("Heatmap_DEGs.tif",   
     width = 3*300,        # 5 x 300 pixels
     height = 3*300,
     res = 300,            # 300 pixels per inch
     pointsize = 10)
pheatmap(Rlog_DEGs, center = TRUE, clustering_distance_rows = "euclidean", clustering_method ="ward.D2", scale ="row", cluster_cols = F, cutree_rows = 2, border_color = NA, gaps_col = c(5,10), fontsize_row = 3,
               fontsize_col = 3)
dev.off()

===
  #Analisis de enriquecimiento de funcion
  
  #Installar ermineJ
  #install.packages("remotes")
  #remotes::install_github("PavlidisLab/ermineR")
  
  #Instalar rjava seguir instrucciones de https://www.r-bloggers.com/2018/02/installing-rjava-on-ubuntu/
  
library(ermineR)
library(rJava)

#Lista GO de Hymenopteramine https://hymenoptera.elsiklab.missouri.edu/hgd-go-annotation
#Descargar archivo y extraer columnas con Gene ID y GO terms
#Si un gen tiene varios GO terms, cada uno va en una fila

goAnnots<-read.table("GO_term_list.txt", sep="\t", header=T) #primera fila nombre de los genes, segunda fila GO
goAnnots<-as.data.frame(goAnnots)
head(goAnnots)
summary(goAnnots)

#names              Goterm         
#Length:159145      Length:159145     
#Class :character   Class :character  
#Mode  :character   Mode  :character

#Armo objeto lista de GOterms por genID
GOList <- split(goAnnots[,2],goAnnots$names)
GOList
class(GOList)

#Formateo para q lo lea el paquete
annotationList <- makeAnnotation(GOList,return = T)
annotationList
summary(annotationList)

#IDs            GeneSymbols         GeneNames           GOTerms         
#Length:9716        Length:9716        Length:9716        Length:9716       
#Class :character   Class :character   Class :character   Class :character  
#Mode  :character   Mode  :character   Mode  :character   Mode  :character   


#Del archivo DEA_results_postlfcShrink saco los valores de LogFc y p-value ajustado para los genes analizados

scoreListFC = read.table("ID_logFc_list.txt", sep="\t", header=T, row.names = 1) # Valores absolutos
scoreListFDR = read.table("ID_padjust_list.txt", sep="\t", header=T, row.names = 1)


#Foldchange como scores:
#Utiliza el valor absoluto del LogFC y por lo tanto no distingue sub o sobreexpresion
#A diferencia del ORA no usa un thershold para seleccionar genes
#En aspects selecciono el tipo de clasificacion de los
#GOterm (procesos biologicos (B), funcion molecular(M) o tipo de componenete celulares(C))

#logTrans y bigIsBetter seteados segun score

#Si la conexión de internet es lenta puede que dé un error x lo que usar
library(httr)
options(timeout = max(10000, getOption("timeout")))

#Procesos biol
ermineR(annotation = annotationList, expression = NULL, aspects = "B",
        scores = scoreListFC, scoreColumn = 1, logTrans = F, bigIsBetter = T,
        test = "GSR", pAdjust = "FDR", iterations = 200000, geneReplicates = "mean", stats = "mean",
        return = F, minClassSize = 10, maxClassSize = 200,
        output = "GSR_LogFc_PB.txt")
#No dio nada

#Funcion molecular
ermineR(annotation = annotationList, expression = NULL, aspects = "M",
        scores = scoreListFC, scoreColumn = 1, logTrans = F, bigIsBetter = T,
        test = "GSR", pAdjust = "FDR", iterations = 200000, geneReplicates = "mean", stats = "mean",
        return = F, minClassSize = 10, maxClassSize = 200,
        output = "GSR_LogFc_FM.txt")

#No dio nada

#Loc celular - Da error y no corre
ermineR(annotation = GOList, expression = NULL, aspects = "C",
        scores = scoreListFC, scoreColumn = 1, logTrans = F, bigIsBetter = T,
        test = "GSR", pAdjust = "FDR", iterations = 200000, geneReplicates = "mean", stats = "mean",
        return = F, minClassSize = 10, maxClassSize = 200,
        output = "GSR_LogFc_LC.txt")

#Con adjusted p-value como score y por eso cambia logTrans a T y bigIsbetter a F

#Procesos biol
ermineR(annotation = annotationList, expression = NULL, aspects = "B",
        scores = scoreListFDR, scoreColumn = 1, logTrans = T, bigIsBetter = F,
        test = "GSR", pAdjust = "FDR", iterations = 200000, geneReplicates = "mean", stats = "mean",
        return = F, minClassSize = 10, maxClassSize = 200,
        output = "GSR_pvalue_PB.txt")
#No dio ningún termino enriquecido

#Fun molecular
ermineR(annotation = annotationList, expression = NULL, aspects = "M",
        scores = scoreListFDR, scoreColumn = 1, logTrans = T, bigIsBetter = F,
        test = "GSR", pAdjust = "FDR", iterations = 200000, geneReplicates = "mean", stats = "mean",
        return = F, minClassSize = 10, maxClassSize = 200,
        output = "GSR_pvalue_FM.txt")

#Solo un término enriquecido