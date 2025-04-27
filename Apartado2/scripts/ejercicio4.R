library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library(EnhancedVolcano)
library(dplyr)
#Cargamos la matriz de cuentas y los metadatos
counts <- read.delim("rawcounts.tsv", row.names = 1)
metadata <- read.delim("metadata.tsv", row.names = 1)  

#creamos una variable que combine agent y time
metadata_at <-
  metadata |>
  mutate(group=paste(agent, time, sep="_"))

#Para crear el objeto DESeq necesitamos que todas las variables sean factores, transformamos group en factor
metadata_at <-
  metadata_at |>
  mutate(group = factor(group))

# Creamos el objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata_at,
                              design = ~ group)

## Prefiltrado. Vamos a eliminar genes
## con tan poquitas cuentas que no merece la pena conservar.
## Así ahorramos memoria, aunque no ganamos nada a efectos estadísticos.
keep <- rowSums(counts(dds)) >= 10 ## Seleccionar genes con más de 10 cuentas en todos los samples
dds <- dds[keep, ]

## Nos preparamos para la DGE haciendo el modelo lineal a partir de la BN
## La función DESeq realiza todos los pasos de DESeq2 desde estimar los size factors
## hasta controla la dispersión
dds2 <- DESeq(dds, test = "Wald")

## Vamos a verificar cómo ha quedado la estimación de la dispersión 
plotDispEsts(dds2)

## Exploremos cómo quedan los cambios de fold entre condiciones con respecto
## a las cuentas normalizadas
plotMA(dds2)

## Obtengamos nuestra lista de genes DEG para cada grupo
res_oht <- results(dds2,
                   contrast = c("group", "Control_24h", "OHT_24h"),
                   alpha = 0.05,
                   pAdjustMethod = "BH",
                   tidy = TRUE)

res_dpn <- results(dds2,
                   contrast = c("group", "Control_24h", "DPN_24h"),
                   alpha = 0.05,
                   pAdjustMethod = "BH",
                   tidy = TRUE)

#Anotamos los genes para mostrar los símbolos para facilitar el estudio biológico
genesID_oht <-mygene::queryMany(my_results$row, scopes="ensembl.gene", fields="symbol", species="human")
genesID_oht <- genesID_oht[!duplicated(genesID_oht$query),]  
res_oht$row <- ifelse(is.na(genesID_oht$symbol),genesID_oht$query,genesID_oht$symbol) 


genesID_dpn <-mygene::queryMany(my_results$row, scopes="ensembl.gene", fields="symbol", species="human")
genesID_dpn <- genesID_dpn[!duplicated(genesID_dpn$query),]  
res_dpn$row <- ifelse(is.na(genesID_dpn$symbol),genesID_dpn$query,genesID_dpn$symbol) 

## la función VST normaliza y estabiliza la varianza de las counts, necesario previo al heatmap
vsd <- vst(dds, blind = TRUE)
rownames(vsd) <- ifelse(is.na(genesID_oht$symbol), genesID_oht$query, genesID_oht$symbol)

res_sig_oht <- res_oht %>% filter(padj < 0.05)

# Heatmap de los genes TOP DGE por p-valor ajustado
mat_oht<- assay(vsd)[head(order(res_oht$padj), 30), ]  
pheatmap(mat_oht)


# Filtramos el metadata para tener solo las muestras 'Control_24h' y 'OHT_24h'
metadata_24_CO <- metadata_at %>%
  filter(group %in% c("Control_24h", "OHT_24h"))

# Filtramos la matriz de expresión VST (vsd) para que solo contenga las muestras de 'Control_24h' y 'OHT_24h'
vsd_24_CO <- vsd[, rownames(metadata_24_CO)]

# Extraemos los 30 genes más significativos según padj
top_genes <- head(order(res_oht$padj), 30)

# Creamos la matriz de expresión con los top genes 
mat_oht <- assay(vsd_24_CO)[top_genes, ]

# Creamos el heatmap con la anotación usando 'group'
pheatmap(mat_oht,
         annotation_col = metadata_filtered[, "group", drop = FALSE], # Usamos 'group' para anotación
         fontsize_row = 8,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")

# Creamos el Volcano plot
EnhancedVolcano(res_oht,
                lab = res_oht$row,
                x = "log2FoldChange",
                y = "padj",
                title = "DEG OHT vs control",
                FCcutoff = 1,
                pCutoff = 0.05,
                subtitle = NULL,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                labSize = 6.0)