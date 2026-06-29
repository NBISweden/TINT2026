library(tidyverse)
library(oligo)
library(clariomdhumantranscriptcluster.db)
library(affycoretools)

params <- list(intermediate = '../../data/intermediate', esetdir='../../data/esets')

## Read combined metadata file which include name of dataset, celfile names, run (batch) information and clinical data
meta <- readRDS(file.path(params$intermediate, "meta.Rds"))
meta <- meta |> filter(dataset %in% c("TINT", "bombiopsi"))

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
esets <- meta |> group_by(dataset2) |> summarise(n=n(), eset=list(readnormannot(pick(everything()))))
nms <- esets$dataset2
esets <- esets$eset
names(esets) <- nms

for (nm in names(esets)){
  saveRDS(esets[[nm]], file.path(params$esetdir, paste0("eset", nm, ".Rds")))
}

##Normalize both bombiopsi datasets together
esets[["bombiopsi"]] <- readnormannot(meta |> filter(dataset=="bombiopsi"))
saveRDS(esets$bombiopsi, file.path(params$esetdir, paste0("eset", "bombiopsi", ".Rds")))

### Normalize TINT and bombiopsi together
esets[["TINTbombiopsi"]] <- readnormannot(meta)
saveRDS(esets$TINTbombiopsi, file.path(params$esetdir, paste0("eset", "TINTbombiopsi", ".Rds")))


