path <- "C:/Users/Martin/Documents/Master/Thesis/CCLE_data/"

################### expression data#############
library(CePa)
expression.data <- 
  read.gct(
    paste0(path,"Expression_Files/CCLE_Expression_Entrez_2012-09-29.gct"))

# convert Affymetrics to HUGO
library(biomaRt)

#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# ensemble.hugo.ids <- getBM(attributes=c("hgnc_symbol",'affy_hg_u133_plus_2'),
#                            filter = 'affy_hg_u133_plus_2',
#                            values=rownames(expression.data),
#                            mart=mart)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
                host = "https://sep2019.archive.ensembl.org",
                path = "/biomart/martservice", archive = FALSE)

gene.ids <- getBM(
  attributes=c('entrezgene_id','ensembl_gene_id','external_gene_name'),
  filters = "entrezgene_id",
  values = sub("\\_.*", "",rownames(expression.data)),
  mart = mart)






# takes very long
expression.data.res.1 <- read.csv(paste0(
  path,"Expression_Files/CCLE_Expression_2012-09-29.res"),
  header = T,sep = "\t",row.names = NULL)

expression.data.res.2 <- read.csv(paste0(
  path,"Expression_Files/CCLE_Expression_Entrez_2012-10-18.res"),
  header = T,sep = "\t",row.names = NULL)

library(cmapR)
expression.data.gctx <- parse_gctx(paste0(
  path,"Expression_Files/CCLE_expression_CN_muts_GENEE_2010-04-16.gctx"))

library(rtrim)
expression.data.tdf <- read_tdf(paste0(
  path,"Expression_Files/CCLE_Expression_Entrez_IGV_2012-09-29.tdf"))


drug.data <- 
  read.csv(paste0(path,
   "Pharmacological_Profiling_Files/CCLE_NP24.2009_Drug_data_2015.02.24.csv"))

profiling.data <- 
  read.csv(paste0(path,
   "Pharmacological_Profiling_Files/CCLE_NP24.2009_profiling_2012.02.20.csv"))

drug.sensitivity.data <- read.csv(paste0(path,
"Pharmacological_Profiling_Files/Drug_sensitivity_(PRISM_Repurposing_Primary_Screen).csv"))

library("readxl")
gnf.data <- read_excel(paste0(path,
             "Pharmacological_Profiling_Files/CCLE_GNF_data_090613.xls"))

mut.data <- 
  read.csv(paste0(path,"Mutation_Files/CCLE_mutations.csv"))

copy.numbers = read.table(paste0(path,
              "Mutation_Files/CCLE_copynumber_byGene_2013-12-03.txt"),sep='\t',
              header = T)



