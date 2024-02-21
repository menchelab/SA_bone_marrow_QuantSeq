#' @param cell.type

get_celltype_DEseq2 = function(ct, counts, smp, min.gene.counts = 3) {
    
    # ct = "Eosinophil"
    
    smp.tbl.subset = smp %>% filter(cell.type == ct)
    
    gene.data = counts[, 1:8]
    
    common.samples = intersect(smp.tbl.subset$Sample_Name, colnames(counts))
    smp.tbl.subset = smp.tbl.subset %>% filter(Sample_Name %in% common.samples)
    counts = counts[, common.samples ]
    
    rownames(counts) = gene.data$gene_id
    
    counts = counts[rowSums(counts) >= min.gene.counts, ]
    
    se = SummarizedExperiment(assays = as.matrix(counts), 
                              rowData = gene.data %>% filter(gene_id %in% rownames(counts)),
                              colData = smp.tbl.subset)
    
    dds = DESeqDataSet(se, design= ~ treatment)
    
    dds <- DESeq(dds)    
    
    resnames = resultsNames(dds) # lists the coefficients
    
    # res <- results(dds, name=resnames[2])
    
    # or to shrink log fold changes association with condition:
    res <- lfcShrink(dds, coef = resnames[2], type="apeglm")

    
    return(list (dds = dds, res = res ))
}


read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, function(y) {
        y = as.data.frame(y)
        # y[[1]] = NULL
        return(y)
    })
    names(x) <- sheets
    x
}

write_excel_allsheets <- function(file.list, filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    write.xlsx(file.list, file = filename)
}

add_mouse_geneinfo = function(filename, suffix = "gene_info", entrez.info.file = NULL) {
    
    if (!exists("entrez.gene.condensed")) {
        
        if (is.null(entrez.info.file)) {
            entrez.info.file = here("scripts/mouse.entrez.info.RDS")
        }
        entrez.gene.all = readRDS(entrez.info.file)
        
        entrez.gene.condensed = do.call(
            rbind,
            lapply(
                entrez.gene.all, 
                function(gene.bunch) {
                    out = sapply(gene.bunch, 
                                 function(x) return(x[c("uid", "name", "description", 
                                                        "summary", "otheraliases", "otherdesignations")])) %>% 
                        t() %>% as.data.frame()
                    return(out)
                } ) ) 
        
        for(cn in colnames(entrez.gene.condensed)) {
            entrez.gene.condensed[[cn]] = unlist(entrez.gene.condensed[[cn]])
        }
    }
    
    outfile = paste0(tools::file_path_sans_ext(filename), "_", suffix, ".xlsx")
    
    file.sheets.list = read_excel_allsheets(filename)
    
    sheet.names = readxl::excel_sheets(filename)
    
    for(sheet.name in sheet.names) {
        if (!exists("entrez.gene.condensed")) {
        
        if (is.null(entrez.info.file)) {
            entrez.info.file = here("scripts/mouse.entrez.info.RDS")
        }
        entrez.gene.all = readRDS(entrez.info.file)
        
        entrez.gene.condensed = do.call(
            rbind,
            lapply(
                entrez.gene.all, 
                function(gene.bunch) {
                    out = sapply(gene.bunch, 
                                 function(x) return(x[c("uid", "name", "description", 
                                                        "summary", "otheraliases", "otherdesignations")])) %>% 
                        t() %>% as.data.frame()
                    return(out)
                } ) ) 
        
        for(cn in colnames(entrez.gene.condensed)) {
            entrez.gene.condensed[[cn]] = unlist(entrez.gene.condensed[[cn]])
        }
    }
    file.sheets.list[[sheet.name]] = left_join(
            dplyr::rename(file.sheets.list[[sheet.name]], name = gene_name), 
            entrez.gene.condensed, by = "name")
    }
    
    write_excel_allsheets(file.sheets.list, filename = outfile)
}



library(dplyr)
# mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

## https://www.biostars.org/p/9567892/#9568018

convert_mouse_to_human <- function(gene_list) { 
    output = c()
    
    for(gene in gene_list) {
        class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
        if( !identical(class_key, integer(0)) ) {
            human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
            for(human_gene in human_genes) {
                output = rbind(c(gene, human_gene), output)
            }
        }
    }
    return (output)
}

convert_human_to_mouse <- function(gene_list) {
    output = c()
    
    if (!exists("mouse_human_genes")) {
        mouse_human_genes <<- read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")    
    }
    
    for(gene in gene_list) {
        class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
        if( !identical(class_key, integer(0)) ) {
            human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
            for(human_gene in human_genes) {
                output = rbind(c(gene, human_gene), output)
            }
        }
    }
    return (output)
}
