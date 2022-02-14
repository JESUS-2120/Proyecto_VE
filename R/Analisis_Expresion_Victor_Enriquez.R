library("recount3")
library("edgeR")
library("variancePartition")
library("limma")
library("pheatmap")
library("RColorBrewer")

## Revisamos los proyectos con datos de humano en recount3
human_projects <- available_projects()

## Exploramos los proyectos disponibles en recount3
proj_info_interactive <- interactiveDisplayBase::display(human_projects)

## Verificamos que solo seleccionamos un renglón.
stopifnot(nrow(proj_info_interactive) == 1)

## Creamos el objeto RSE
rse_gene_SRP051765 <- create_rse(proj_info_interactive)

## class: RangedSummarizedExperiment
## dim: 63856 36
## metadata(8): time_created recount3_version ... annotation recount3_url
## assays(1): raw_counts
## rownames(63856): ENSG00000278704.1 ENSG00000277400.1 ... ENSG00000182484.15_PAR_Y
## ENSG00000227159.8_PAR_Y
## rowData names(10): source type ... havana_gene tag
## colnames(36): SRR7268163 SRR7268164 ... SRR7268197 SRR7268198
## colData names(175): rail_id external_id ... recount_pred.curated.cell_line BigWigURL

## Convertimos las cuentas por nucleotido a cuentas por lectura
assay(rse_gene_SRP051765,"counts") <- compute_read_counts(rse_gene_SRP051765)

## class: RangedSummarizedExperiment
## dim: 63856 36
## metadata(8): time_created recount3_version ...
## annotation recount3_url
## assays(2): raw_counts counts
## rownames(63856): ENSG00000278704.1
## ENSG00000277400.1 ...
## ENSG00000182484.15_PAR_Y
## ENSG00000227159.8_PAR_Y
## rowData names(10): source type ... havana_gene
## tag
## colnames(36): SRR7268163 SRR7268164 ...
## SRR7268197 SRR7268198
## colData names(175): rail_id external_id ...
## recount_pred.curated.cell_line BigWigURL

## Tras analizar los sample_attributes pudimos observar que no hay ningun
##problema con los datos
rse_gene_SRP051765$sra.sample_attributes[sample(1:41,size = 3)]

## Expandimos la informacion
rse_gene_SRP051765 <- expand_sra_attributes(rse_gene_SRP051765)

colData(rse_gene_SRP051765)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP051765)))
]

## Pasamos de character a factor
rse_gene_SRP051765$sra_attribute.drug_treatment <- factor(rse_gene_SRP051765$sra_attribute.drug_treatment)

rse_gene_SRP051765$sra_attribute.origin <- factor(rse_gene_SRP051765$sra_attribute.origin)

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP051765)[
  ,
  grepl("^sra_attribute.[drug_treatment|source_name|origin]", colnames(colData(rse_gene_SRP051765)))
]))

## Buscamos diferencias entre las muestras generadas in vitro y las muestras obtenidad de un modelo de xenoinjerto
rse_gene_SRP051765$origen <- factor(ifelse(rse_gene_SRP051765$sra_attribute.origin == "cells in vitro", "in_vitro", "xenoinjerto"))

table(rse_gene_SRP051765$origen)

rse_gene_SRP051765$assigned_gene_prop <- rse_gene_SRP051765$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP051765$recount_qc.gene_fc_count_all.total

with(colData(rse_gene_SRP051765), tapply(assigned_gene_prop, origen, summary))

## Analizamos los datos para ver si eliminaremos muestras de mala calidad
hist(rse_gene_SRP051765$assigned_gene_prop)

table(rse_gene_SRP051765$assigned_gene_prop > 0.4)

table(rse_gene_SRP051765$assigned_gene_prop > 0.45)

## Normalizamos los datos
dge <- DGEList(
  counts = assay(rse_gene_SRP051765, "counts"),
  genes = rowData(rse_gene_SRP051765)
)

dge <- calcNormFactors(dge)

## Definimos nuestro modelo
datos <- data.frame(
 "treatment" = rse_gene_SRP051765$sra_attribute.drug_treatment,
 "origin" = rse_gene_SRP051765$sra_attribute.origin,
 "line" = rse_gene_SRP051765$sra_attribute.cell_line,
 "type" = rse_gene_SRP051765$sra_attribute.cell_type
 )

design <- model.matrix( ~ origin, datos)

vobjGenes <- voom(dge, design)

formula <- ~ (1|line) + (1|type) + (1|treatment) + (1|origin)

var_exp<- fitExtractVarPartModel( vobjGenes, formula, datos )

plotVarPart(var_exp)

modelo <- model.matrix(~ sra_attribute.origin + sra_attribute.drug_treatment,
                    data = colData(rse_gene_SRP051765)
)

## Comenzamos con el analisis de expresion
vGene <- voom(dge, modelo, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = c(3,4,5,6,7),
  number = nrow(rse_gene_SRP051765),
  sort.by = "none"
)

## contamos los genes diferencialmente expresados entre los distintos tratamientos
table(de_results$adj.P.Val < 0.05)

## Agregamos los nombres de los genes
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

nombres_heatmap <- rownames(exprs_heatmap)

rownames(exprs_heatmap) <- rowRanges(rse_gene_SRP051765)$gene_name[
  match(rownames(exprs_heatmap), rowRanges(rse_gene_SRP051765)$gene_id)
]

## Creamos una tabla con información de las muestras
## y con nombres de columnas mas faciles de leer
df <- as.data.frame(colData(rse_gene_SRP051765)[, c("sra_attribute.origin", "sra_attribute.drug_treatment")])
colnames(df) <- c("Origin", "Treatment")

## Hagamos un heatmap
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df
)

## Convertimos los  los tratamientos a colores
col.Treatment <- df$Treatment

levels(col.Treatment) <- brewer.pal(nlevels(col.Treatment), "Set1")

col.Treatment <- as.character(col.Treatment)

## MDS por tratamiento
plotMDS(vGene$E, labels = df$Treatment, col = col.Treatment)

## Hacemos una grafica para el gen mas expresado diferencialmente
i <- which.min(de_results$P.Value)

df_temp <- data.frame(
  Expression = vGene$E[i,],
  Tratamiento = rse_gene_SRP051765$sra_attribute.drug_treatment,
  Origen = rse_gene_SRP051765$sra_attribute.origin
)

ggplot(df_temp, aes(y = Expression, x = Origen, fill = Tratamiento)) +
  geom_boxplot() +
  theme_dark(base_size = 20)
