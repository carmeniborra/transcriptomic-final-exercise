# Cargamos las librerias
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("EnhancedVolcano")

# Cargar matriz de conteos y metadatos
counts <- read.delim("rawcounts.tsv", row.names = 1)
metadata <- read.delim("metadata.tsv", row.names = 1)

# Comprobar que coiciden los nombres de columnas en counts con filas en metadata
all(colnames(counts) == rownames(metadata)) #debe ser TRUE

# Filtrar muestras de 24 horas
metadata_24h <- metadata[metadata$time == "24h", ]
counts_24h <- counts[, rownames(metadata_24h)]


# Crear el objeto DESeqDataSet
dds_24h <- DESeqDataSetFromMatrix(countData = counts_24h,
                                  colData = metadata_24h,
                                  design = ~ agent)

## Prefiltrado. Vamos a eliminar genes
## con tan poquitas cuentas que no merece la pena conservar.
## Así ahorramos memoria, aunque no ganamos nada a efectos estadísticos.
keep <- rowSums(counts(dds_24h) >= 10) >= 4 ## conservar genes que tengan al menos 10 cuentas en al menos 1/3 de las muestras
dds_24h <- dds_24h[keep, ]

## Antes de hacer la DGE, hagamos un análisis exploratorio
## la función VST normaliza y estabiliza la varianza de las counts. Ideal para
## clustering o PCA.

## La diferencia entre VST y la función DESeq radica en que la expresión diferencial
## no usa las cuentas normalizadas (y estabilizadas) per sé, sino que los factores de 
## normalización y dispersión se incluyen en un modelo BN con las cuentas crudas.
## Para visualización, clustering etc. la función vst transforma las counts y estabiliza
## la varianza APARTE para que podamos trabajar con cuentas ya "listas".

vsd_24h <- vst(dds_24h, blind = TRUE) # no se tendrá en cuenta el diseño experimental en la transformación, para exploración.
plotPCA(vsd_24h, intgroup = "agent") # mostrará si las muestras se agrupan según el tratamiento.
plotPCA(vsd_24h, intgroup = "patient") #mostrará si las muestras se agrupan por paciente

## Calculamos las distancias a partir de las cuentas normalizadas y variance-stabilized (vst)
sampleDists <- dist(t(assay(vsd_24h)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_24h$agent, vsd$patient, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main = "Distancia entre muestras (DPN vs Control, 24h)")

## Nos preparamos para la DGE haciendo el modelo lineal a partir de la BN
## La función DESeq realiza todos los pasos de DESeq2 desde estimar los size factors
## hasta controla la dispersión
dds2_24h <- DESeq(dds_24h, test = "Wald")

## Vamos a verificar cómo ha quedado la estimación de la dispersión 
plotDispEsts(dds2_24h)

## Exploremos cómo quedan los cambios de fold entre condiciones con respecto
## a las cuentas normalizadas
plotMA(dds2_24h)

## Obtengamos nuestra lista de genes DEG
my_results <- results(object = dds2_24h,
                      contrast = c("agent", "DPN", "Control"),
                      alpha = 0.05,
                      pAdjustMethod = "BH",
                      tidy = TRUE)

##Anotamos los genes para mostrar los símbolos para facilitar el estudio biológico
# Anotar los genes con símbolos humanos
genesID <- mygene::queryMany(res_DPN$row, 
                             scopes = "ensembl.gene", 
                             fields = "symbol", 
                             species = "human")

# Eliminar duplicados para asegurar que cada gen tenga una anotación
genesID <- genesID[!duplicated(genesID$query), ]

# Fusionar la información anotada con los resultados
res_DPN$row <- ifelse(is.na(genesID$symbol),
                      genesID$query,
                      genesID$symbol)

# Convertir factores si es necesario
metadata$Treatment <- factor(metadata$Treatment, levels = c("Control", "OHT", "DPN"))
metadata$Time <- factor(metadata$Time, levels = c("0h", "24h"))

# Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Treatment)

## Suele ser buena idea establecer un corte a priori de log fold
res_DPN_thresh <- results(object = dds2_24h,
                          contrast = c("agent", "DPN", "Control"),
                          lfcThreshold = 1,
                          alpha = 0.05,
                          pAdjustMethod = "BH",
                          tidy = TRUE)
##Anotamos los genes para mostrar los símbolos para facilitar el estudio biológico
library(mygene)

genesID_threshold <- mygene::queryMany(res_DPN_thresh$row,
                                       scopes = "ensembl.gene",
                                       fields = "symbol",
                                       species = "human")

genesID_threshold <- genesID_threshold[!duplicated(genesID_threshold$query), ]

res_DPN_thresh$row <- ifelse(is.na(genesID_threshold$symbol),
                             genesID_threshold$query,
                             genesID_threshold$symbol)
write.csv(res_DPN_thresh, "DEG_DPN_vs_Control_24h_thresh_annotated.csv", row.names = FALSE)

## Heatmap de los genes TOP DGE por p-valor ajustado
mat <- assay(vsd_24h)[head(order(res_DPN_thresh$padj), 30), ]
rownames(mat) <- res_DPN_thresh$row[head(order(res_DPN_thresh$padj), 30)]
pheatmap(mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         main = "TOP 30 genes DPN vs Control (24h)")

# Creamos el Volcano plot

EnhancedVolcano(res_DPN_thresh,
                lab = res_DPN_thresh$row,
                x = "log2FoldChange",
                y = "padj",
                title = "DEG DPN vs Control (24h)",
                FCcutoff = 1,
                pCutoff = 0.05,
                subtitle = NULL,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                labSize = 6.0)

