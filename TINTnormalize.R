library(tidyverse)
library(oligo)
library(clariomdhumantranscriptcluster.db)
library(affycoretools)

params <- list(datadir = '../../data/TINT2026', intermediate = '../../data/TINT2026/intermediate', esetdir = '../../data/TINT2026/esets')

## Read combined metadata file which include name of dataset, celfile names, run (batch) information and clinical data
meta <- read_delim(file.path(params$datadir, "sdrf.tsv"), "\t")
meta <- meta |> transmute(sampleid = `Source Name`, 
  dataset = `Characteristics[collection method]`, 
  patientid = `Characteristics[individual]`,
  ISUPgr = `Characteristics[isup group]`,
  age = `Characteristics[age]`,
  run = `Characteristics[experimental run]`,
  sampletype = `Characteristics[sample type]`,
  disease = `Characteristics[disease]`) |> 
  mutate(dataset2 = ifelse(dataset=="biopsy", ifelse(sampletype=="Bx-tumor", "Bxtumor", "Bxbenign"), dataset),
          celfile = file.path(params$datadir, "CELS", paste0(sampleid, ".CEL")))

## Read cel files, RMA normalize and annotate
readnormannot <- function(pheno) {
  rawData <- read.celfiles(pheno$celfile, verbose=FALSE)
  colnames(rawData) <- pheno$sampleid
  pData(rawData) <- pheno |> dplyr::select(where(~ !(all(is.na(.)) | all(. == "")))) |> mutate(rn=sampleid) |> column_to_rownames("rn")
  eset <- oligo::rma(rawData, target="core")
  eset <- annotateEset(eset, clariomdhumantranscriptcluster.db, columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME","GENETYPE"))
  return(eset)
}

##Normalize separately per dataset
esets <- list()
esets$Bxbenign <- readnormannot(meta |> filter(dataset2=="Bxbenign"))
esets$TINT <- readnormannot(meta |> filter(dataset2==TINT))
##Normalize both biopsy datasets together
esets[["biopsy"]] <- readnormannot(meta |> filter(dataset=="biopsy"))
### Normalize radical prostatectomy and biopsies together
esets[["TINTbiopsy"]] <- readnormannot(meta)

for (nm in names(esets)){
  saveRDS(esets[[nm]], file.path(params$esetdir, paste0("eset", nm, ".Rds")))
}
