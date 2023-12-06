library(here)
library(biomaRt)
library(rentrez)
library(tidyverse)

library(AnnotationDbi)
library(org.Mm.eg.db)
library(xml2)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")


mouse.mapping = as.list(org.Mm.egALIAS2EG)

gene.entrez = mouse.mapping %>% unlist() %>% enframe()

unique.entrez = unique(gene.entrez$value)

# df <- AnnotationDbi::select(org.Mm.eg.db, head(gene.names), columns(org.Mm.eg.db), "SYMBOL")


# filename = here("notebooks/all_BM_cells_scNym_DEGs.xlsx")
# 
# sheet.names = readxl::excel_sheets(filename)
# 
# gene.list = readxl::read_excel(filename, sheet = sheet.names[1])
# 
# gene.names = gene.list$names
# 
# genes.info = entrez_summary(db = "gene", id = unique.entrez[1:400], use_history = TRUE)

all.genes.info = list()
for (i in 1:(length(unique.entrez) %/% 300 + 1) ) {
    cat(i, "\n")
    interval = seq((i - 1) * 300,  min(i * 300, length(unique.entrez)) )
    all.genes.info[[i]] = entrez_summary(db = "gene", id = unique.entrez[interval],
                                         use_history = TRUE)
}

saveRDS(object = all.genes.info, file = "mouse.entrez.info.RDS")

