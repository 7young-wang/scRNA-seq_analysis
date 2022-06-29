# scRNA-seq_analysis

This is a log describing how to analyze scRNA-seq data from GSE158142. My main purpose is to 
1. identify DE genes (with p-value)
2. identify clusterMarkers
3. construct WGCNA gene network if I can

I mainly follow the [Seurat official tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow-1) accompanied by [UCD and UCSF Workshop](https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART6.html). But in this log I will also specify what each function is doing (based on my understanding). 

## Step 1: Prepare and Merge the Data

The data I am working at includes 9 replicates of 5dpf zebrafish. So the first step is to merge them into one giant Seurat object and treat the merged object as one individual fish. I started doing this with `IntegrateData` workflow but then found out that is mainly to compare different groups rather than integrating replicates (but it's possible I use it for temperal analysis or KO vs. WT RNA-seq analysis).

Load the Seurat Library and create an empty list.

```
library(Seurat)
seurat_5dpf_list <- NULL
```

Create a pipeline loading and converting all the raw seurat data I have from the GSE quary into Seurat objects.

```
import_seurat_data <- function(filename) {
  filename <- paste("./GSE158142_RAW/", filename, sep = "")
  expression_matrix <- ReadMtx(
    mtx = paste(filename, "-matrix.mtx.gz", sep = ""), features = paste(filename, "-genes.tsv.gz", sep = ""),
    cells = paste(filename, "-barcodes.tsv.gz", sep = "")
  )
  seurat_5dpf_list <<- append(seurat_5dpf_list, CreateSeuratObject(counts = expression_matrix)) #<<- is global assignment
}
#make a list of all 5dpf file names
all5dpf <- c("GSM4793214_zBr5dpf2_S3", "GSM4793215_zBr5dpf3_S4",
             "GSM4793216_comb-zBr5dpf4_S1", "GSM4793218_comb-zBr5dpf6_S3", 
             "GSM4793217_comb-zBr5dpf5_S2", "GSM4793219_comb-zBr5dpf7_S4",
             "GSM4793220_comb-zBr5dpf8_S5", "GSM4793221_comb-zBr5dpf9_S6")
one5dfp <- ReadMtx(
  mtx = paste("./GSE158142_RAW/GSM4793213_zBr5dpf1_S7", "-matrix.mtx.gz", sep = ""), features = paste("./GSE158142_RAW/GSM4793213_zBr5dpf1_S7", "-genes.tsv.gz", sep = ""),
  cells = paste("./GSE158142_RAW/GSM4793213_zBr5dpf1_S7", "-barcodes.tsv.gz", sep = "")
)
one5dfp <- CreateSeuratObject(counts = one5dfp)
lapply(X = all5dpf, FUN = import_seurat_data)
```

Merge the single Seurat object with the Seurat object list (9 objects in total). Give each object the cell ids though they are not used downstream.

```
comb_5dpf <- merge(one5dfp, y = seurat_5dpf_list, add.cell.ids = c("1_S7", "2_S3", "3_S4", "4_S1", "6_S3", "5_S2", "7_S4", "8_S5", "9_S6"))
```

Now the preparation step is down. 

## Step 2: Filter the Cells





