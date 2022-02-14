library("recount3")
library("edgeR")
library("variancePartition")

## Revisamos los proyectos con datos de humano en recount3
human_projects <- available_projects()

## Exploramos los proyectos disponibles en recount3
proj_info <- subset(
  human_projects,
  project == "SRP051765" & project_type == "data_sources"
)
## Crea un objetio de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP051765 <- create_rse(proj_info)

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
rse_gene_SRP051765$sra_attribute.source_name <- factor(rse_gene_SRP051765$sra_attribute.source_name)

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
dataf <- data.frame(
  "cell line" = rse_gene_SRP051765$sra_attribute.cell_line,
  "cell type" = rse_gene_SRP051765$sra_attribute.cell_type,
  "treatment" = rse_gene_SRP051765$sra_attribute.drug_treatment,
  "origin" = rse_gene_SRP051765$sra_attribute.origin,
  "name" = rse_gene_SRP051765$sra_attribute.source_name
)

design <- model.matrix( ~ treatment, dataf)

vobjGenes <- voom(dge, design )

formula <- ~ (1|cell.line) + (1|name) + (1|treatment) + (1|origin)

varPart <- fitExtractVarPartModel( vobjGenes, formula, dataf )
