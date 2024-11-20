
# Pathway analysis functions 
#
#
# UPDATE in 2023-11-21 : modify functions 
# - set gene.annotation entrez to character
# - modify setentrezToSymbol , addGeneFoldChange with tidyverse left_join
# 
# UPDATE in 2023-10-19 : using mouse/human pathway both .... 
# 
# 


library(glue)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(data.table)


# pathway.version <- 'E:/2022_SysPharm_HSY/13_Scripts/PathwayAnalysis/Pathways_202310'

message(glue('Pathway version : {pathway.version}'))

pwfiles <- list.files(path = pathway.version,pattern = "term.*\\.RDS$")
pwfiles.path <- list.files(path = pathway.version,pattern = "term.*\\.RDS$", full.names = T)
pwnames <- gsub(pwfiles, replacement = '',pattern = '\\.RDS$')
pwlist <- list()
for(i in 1:length(pwfiles)){
    pwlist[[pwnames[i]]] <- readRDS(file = pwfiles.path[i])
}

pwlist.human <- pwlist[grepl(names(pwlist), pattern = '^human') ]
pwlist.mouse <- pwlist[grepl(names(pwlist), pattern = '^mouse') ]

gene.pw.annotation <- readRDS(file.path(pathway.version, 'gene.annotation.RDS'))
gene.pw.annotation$entrez <- as.character(gene.pw.annotation$entrez)

gene.pw.annotation.human <- gene.pw.annotation %>% dplyr::filter(tax_id == '9606')
gene.pw.annotation.mouse <- gene.pw.annotation %>% dplyr::filter(tax_id == '10090')

all.entrez.id <- as.character(unique(gene.pw.annotation$entrez))


getAllenrichPathways <- function(symbols=NULL, 
                                 entrez=NULL, 
                                 pwlist = pwlist.human,
                                 all.entrez.id = all.entrez.id,
                                 dds.all.entrez.id = NULL,
                                 pvalcutoff = 0.2, 
                                 minCount = 2,
                                 GO.Evidence.IEA = TRUE,
                                 addGeneFoldChange=TRUE, 
                                 fc_data = NULL,
                                 minGSSize = 2,
                                 maxGSSize = Inf,
                                 gene.pw.annotation = gene.pw.annotation,
                                 simplifyByTopGenes = FALSE, 
                                 GO_CC_MF = TRUE,
                                 only.go.kegg = TRUE, 
                                 simplify = TRUE, 
                                 make2df = TRUE ){
    
    pathways <- list()
    
    if(!exists('all.entrez.id')){
        message('...no all entrez id'); return('')}
    
    if(!is.null(dds.all.entrez.id)){
        message("intersect dds all entrez id & base all entrez id ")
        all.entrez.id <- intersect(dds.all.entrez.id, all.entrez.id)
    }
    
    
    if(addGeneFoldChange == TRUE){
        if(is.null(fc_data)){ message("...add fold change data ..."); return('') }}
    
    message(glue('Enrichment analysis with \n {paste0(names(pwlist), collapse = " ")}'))
    
    if(is.null(entrez)){
        message("Convert Symbol To Entrez ID...")
        entrez <- gene.pw.annotation$entrez[gene.pw.annotation$SYMBOL %in% symbols]
        }

    names(pwlist) <- gsub(x = names(pwlist), pattern = '^(human|mouse)\\.', replacement = '')
    
    if(!GO.Evidence.IEA){
        go.pways <- grep(names(pwlist) , pattern = 'GO.*\\.term2gene')
        pwlist[go.pways] <- sapply(pwlist[go.pways], 
                                   function(x){x <- x[x$Evidence != 'IEA',]},
                                   simplify = F, USE.NAMES = F)
    }
    
    
    cat('GO BP enrichment analysis...\n')
    
    try(expr = {
        
        pathways[['GO BP']] <- enricher(gene = entrez, 
                                        TERM2GENE = pwlist$GOBP.term2gene, 
                                        TERM2NAME = pwlist$GOBP.term2name,
                                        pvalueCutoff = pvalcutoff,
                                        qvalueCutoff = 1,
                                        minGSSize = minGSSize,
                                        maxGSSize = maxGSSize,
                                        pAdjustMethod = "BH", 
                                        universe = all.entrez.id
        )
        
        pathways[['GO BP']]@result <- pathways[['GO BP']]@result[
                            which(pathways[['GO BP']]@result$p.adjust < pvalcutoff & 
                                      pathways[['GO BP']]@result$Count >= minCount), ] 

        
        if(simplify){
            # pathways[['GO BP']] <- simplify(pathways[['GO BP']])
        }
        
    })
    
    if(GO_CC_MF == TRUE){
        
        cat('GO CC enrichment analysis...\n')
        try(expr = {
            
            pathways[['GO CC']] <- enricher(gene = entrez, 
                                            TERM2GENE = pwlist$GOCC.term2gene, 
                                            TERM2NAME = pwlist$GOCC.term2name,
                                            pvalueCutoff = pvalcutoff,
                                            qvalueCutoff = 1,
                                            minGSSize = minGSSize,
                                            maxGSSize = maxGSSize,
                                            pAdjustMethod = "BH",
                                            universe = all.entrez.id
            )
            pathways[['GO CC']]@result <- pathways[['GO CC']]@result[
                which(pathways[['GO CC']]@result$p.adjust < pvalcutoff & 
                          pathways[['GO CC']]@result$Count >= minCount), ] 

            
            if(simplify){
                #   pathways[['GO CC']] <- simplify(pathways[['GO CC']])
            }
            
        })
        
        cat('GO MF enrichment analysis...\n')
        try(expr = {
            
            pathways[['GO MF']] <- enricher(gene = entrez, 
                                            TERM2GENE = pwlist$GOMF.term2gene, 
                                            TERM2NAME = pwlist$GOMF.term2name,
                                            pvalueCutoff = pvalcutoff,
                                            qvalueCutoff = 1,
                                            minGSSize = minGSSize,
                                            maxGSSize = maxGSSize,
                                            pAdjustMethod = "BH",
                                            universe = all.entrez.id
                                            )
            pathways[['GO MF']]@result <- pathways[['GO MF']]@result[
                which(pathways[['GO MF']]@result$p.adjust < pvalcutoff & 
                          pathways[['GO MF']]@result$Count >= minCount), ] 

            
            if(simplify){
                # pathways[['GO MF']] <- simplify(pathways[['GO MF']])
            }
            
        })   
        
    }
    
    cat('KEGG enrichment analysis...\n')
    try(expr = {
        pathways[['KEGG']] <- enricher(gene = entrez, 
                                       TERM2GENE = pwlist$KEGG.term2gene, 
                                       TERM2NAME = pwlist$KEGG.term2name,
                                       pvalueCutoff = pvalcutoff,
                                       qvalueCutoff = 1,
                                       minGSSize = minGSSize,
                                       maxGSSize = maxGSSize,
                                       pAdjustMethod = "BH",
                                       universe = all.entrez.id
                                       )
        pathways[['KEGG']]@result <- pathways[['KEGG']]@result[
            which(pathways[['KEGG']]@result$p.adjust < pvalcutoff & 
                      pathways[['KEGG']]@result$Count >= minCount), ] 

        
    })
    
    
    if(only.go.kegg ==FALSE){
        try(expr = {
            cat('KEGG Module enrichment analysis...\n')
            pathways[['KEGG.M']] <- enricher(gene = entrez, 
                                             TERM2GENE = pwlist$MKEGG.term2gene, 
                                             TERM2NAME = pwlist$MKEGG.term2name,
                                             pvalueCutoff = pvalcutoff,
                                             qvalueCutoff = 1,
                                             minGSSize = minGSSize,
                                             maxGSSize = maxGSSize,
                                             pAdjustMethod = "BH",
                                             universe = all.entrez.id
                                            )
            pathways[['KEGG.M']]@result <- pathways[['KEGG.M']]@result[
                which(pathways[['KEGG.M']]@result$p.adjust < pvalcutoff & 
                          pathways[['KEGG.M']]@result$Count >= minCount), ] 

        })
        
        
        try(expr = {
            cat('WikiPathway enrichment analysis...\n')
            pathways[['WikiPathways']] <- enricher(gene = entrez, 
                                                   TERM2GENE = pwlist$wp.term2gene, 
                                                   TERM2NAME = pwlist$wp.term2name,
                                                   pvalueCutoff = pvalcutoff,
                                                   qvalueCutoff = 1,
                                                   minGSSize = minGSSize,
                                                   maxGSSize = maxGSSize,
                                                   pAdjustMethod = "BH",
                                                   universe = all.entrez.id
                                            )
            pathways[['WikiPathways']]@result <- pathways[['WikiPathways']]@result[
                which(pathways[['WikiPathways']]@result$p.adjust < pvalcutoff & 
                          pathways[['WikiPathways']]@result$Count >= minCount), ] 

        })
        
        
        try(expr = {
            cat('Reactome enrichment analysis...\n')
            pathways[['Reactome']] <- enricher(gene = entrez, 
                                               TERM2GENE = pwlist$reactome.term2gene, 
                                               TERM2NAME = pwlist$reactome.term2name,
                                               pvalueCutoff = pvalcutoff,
                                               qvalueCutoff = 1,
                                               minGSSize = minGSSize,
                                               maxGSSize = maxGSSize,
                                               pAdjustMethod = "BH",
                                               universe = all.entrez.id
                                            )
            pathways[['Reactome']]@result <- pathways[['Reactome']]@result[
                which(pathways[['Reactome']]@result$p.adjust < pvalcutoff & 
                          pathways[['Reactome']]@result$Count >= minCount), ] 

        })
        
        
        try(expr = {
            cat("MsigDB Hallmark enrichment analysis...\n")
            pathways[['MsigHall']] <- enricher(gene = entrez, 
                                               TERM2GENE = pwlist$hallmark.term2gene, 
                                               pvalueCutoff = pvalcutoff,
                                               qvalueCutoff = 1,
                                               minGSSize = minGSSize,
                                               maxGSSize = maxGSSize,
                                               pAdjustMethod = "BH",
                                               universe = all.entrez.id
                                            )
            pathways[['MsigHall']]@result <- pathways[['MsigHall']]@result[
                which(pathways[['MsigHall']]@result$p.adjust < pvalcutoff & 
                          pathways[['MsigHall']]@result$Count >= minCount), ] 

        })
    }
    
    modules <- names(pathways[!sapply(pathways, is.null)])
    
    if(addGeneFoldChange == TRUE){
        
        for(i in 1:length(modules)){
            try(expr = {
                pathways[[modules[i]]] <- addGeneFoldChange(pathways[[modules[i]]], fc_data = fc_data)
                message(glue('{modules[i]} add gene fold change'))
            })
        }
    }else{
        
        for(i in 1:length(modules)){
            try(expr = {
                pathways[[modules[i]]] <- setentrezToSymbol(setentrezToSymbol[[modules[i]]], 
                                                            gene.annotation = gene.pw.annotation)
                message(glue('{modules[i]} set Entrez ID to SYMBOL'))
            })
        
    }
    }

    
    if(simplifyByTopGenes == TRUE){
        
        for(i in 1:length(modules)){
            try(expr = {
                pathways[[modules[i]]] <- simplifyByTopGenes(pathways[[modules[i]]])
                message(glue('{modules[i]} simplified by Top Genes'))
            })
        }
    }
    
    print(cat(glue('{paste0(modules, collapse =", ")} loaded\\n')))
    
    if(make2df == TRUE){
        
        message(glue('Pathways To DataFrame'))
        pathways <- pathways_to_df(pathways)
    
    }
    
    return(pathways)
    
}


setentrezToSymbol <- function(enrichresult, gene.annotation = gene.annotation, gsea = FALSE){
    
    x <- enrichresult@result
    if(is.null(x)){
        message('no enrich result...')
        break
    }
    message(' set entrez ID to gene symbol using gene annotation file')
    
    
    x <- x %>% separate_rows(geneID, sep = "/") %>% 
        left_join(gene.annotation %>% dplyr::select(entrez, SYMBOL) %>% distinct(), 
                  by = c('geneID'='entrez')) %>% 
        mutate(geneID = SYMBOL) %>% 
        dplyr::select(-SYMBOL) %>% 
        group_by(ID) %>% 
        mutate(geneID = paste0(geneID, collapse = "/")) %>% distinct()
    
    enrichresult@result <- x
    
    return(enrichresult)
    
}


splitGeneSet <- function(gs){
    gs <- stringr::str_split(string = gs, pattern = '/')[[1]]
    return(gs)
}

addGeneFoldChange <- function(enrichresult, fc_data = top_genes_filt, top_genes = 10){
    
    # add gene fold change data to GO Enrichment analysis result. 
    # x is enrichGo@result object output of clusterprofiler::enrichGO
    # fc_data : dataframe contains fold change data of genes. columns :  Symbol, logFC
    # top_genes : make a new column that only contains N number of top fold change genes in each GO term. 
    
    x <- enrichresult@result
    if(is.null(x)){
        message('no enrich result...')
        break
    }
    
    x <- x %>% separate_rows(geneID, sep = "/") %>% 
        left_join(fc_data %>% dplyr::select(gene_id, SYMBOL,log2FoldChange) %>% distinct(gene_id, .keep_all = T), 
                  by = c('geneID'='gene_id')) %>% 
        group_by(ID) %>% 
        dplyr::arrange(desc(abs(log2FoldChange))) %>% 
        ungroup() %>% 
        dplyr::mutate(geneID_FC = paste0(SYMBOL, '(',round(2^(log2FoldChange),2),')')) %>% 
        dplyr::select(-log2FoldChange) %>% 
        group_by(ID) %>% 
        dplyr::mutate(geneID_FC = paste0(geneID_FC, collapse = '/'), 
                      geneID = paste0(geneID, collapse = '/'),
                      SYMBOL = paste0(SYMBOL, collapse = '/')) %>% 
        distinct() %>% 
        arrange(pvalue)
    
    
    if(!is.null(top_genes)){
        
        pattern <- paste0("^((?:[^/]+/){", top_genes - 1, "}[^/]+).*$")
        x <- x %>% mutate(geneID_FC_top_genes =gsub(pattern,'\\1',x = geneID_FC))
        
    }
    
    
    enrichresult@result <- x
    
    return(enrichresult)
}


simplifyByTopGenes <- function(enrichresult){
    
    # add gene fold change data to GO Enrichment analysis result. 
    # x is enrichGo@result object output of clusterprofiler::enrichGO
    # fc_data : dataframe contains fold change data of genes. columns :  Symbol, logFC
    # top_genes : make a new column that only contains N number of top fold change genes in each GO term. 
    
    x <- data.frame(enrichresult)
    if(is.null(x)){
        message('no enrich result...')
        break
    }
    x <- x %>% distinct(geneID_FC_top_genes, .keep_all = TRUE)
    enrichresult@result <- x
    return(enrichresult)
    
}

pathways_to_df <- function(enrichresult){
    
    resultframe_list <- list()
    
    for(i in 1:length(enrichresult)){
        
        if(!is.null(enrichresult[[i]])){
            df <- enrichresult[[i]]@result
            
            if(nrow(df) !=0)
            {
                df$geneset_count <- as.numeric(gsub(x=df$BgRatio,pattern = "/.*",replacement = ""))
                df$bg_count <- as.numeric(gsub(x=df$BgRatio,pattern = "^.*/",replacement = ""))
                df$DEG_count <- as.numeric(gsub(x=df$GeneRatio, pattern = "^.*/", replacement = ""))
                
                df$GeneRatio <- df$GeneRatio
                df$GeneRatio_calc <- sapply(df$GeneRatio, function(txt) eval(parse(text=txt)))
                
                df$Description <- gsub(x = df$Description, pattern = "_", replacement = " ") 
                df$Description <- gsub(x = df$Description, pattern = " - Mus musculus \\(house mouse\\)", replacement = "") 
                df$Database <- names(enrichresult)[i]
                df <- df %>% dplyr::select(Database, ID, Description, Count, geneset_count, 
                                           DEG_count, bg_count, GeneRatio, GeneRatio_calc,  everything())
                resultframe_list[[names(enrichresult[i])]] <- df
            }
        }
        
    }
    
    resultframe_df <- rbindlist(resultframe_list, fill = TRUE, use.names = TRUE) %>% as_tibble()
    
    return(resultframe_df)
    
}

getAllGSEA <- function(geneList, 
                       pvalcutoff = 0.2, 
                       minCount = 2,
                       GO.Evidence.IEA = TRUE,
                       GO_CC_MF = TRUE,
                       only.go.kegg = TRUE, 
                       addGeneFoldChange = TRUE, 
                       fc_data = NULL,
                       make2df = FALSE,
                       pwlist = pwlist.human,
                       minGSSize = 2,
                       maxGSSize = Inf,
                       gene.pw.annotation = gene.pw.annotation,
                       simplifyByTopGenes = FALSE
){
    
    #### geneList :
    #### also should modify addgenefoldchange ...
    
    pathways <- list()
    
    
    gene.pw.annotation$entrez <- as.character(gene.pw.annotation$entrez)
    
    if(addGeneFoldChange == TRUE){
        if(is.null(fc_data)){ message("...add fold change data ..."); return('') }}
    
    message(glue('Enrichment analysis with \n {paste0(names(pwlist), collapse = " ")}'))
    
    # if(is.null(entrez)){
    #     message("Convert Symbol To Entrez ID...")
    #     entrez <- gene.annotation$entrez[gene.annotation$SYMBOL %in% symbols]
    # }
    
    names(pwlist) <- gsub(x = names(pwlist), pattern = '^(human|mouse)\\.', replacement = '')
    
    if(!GO.Evidence.IEA){
        go.pways <- grep(names(pwlist) , pattern = 'GO.*\\.term2gene')
        pwlist[go.pways] <- sapply(pwlist[go.pways], 
                                   function(x){x <- x[x$Evidence != 'IEA',]},
                                   simplify = F, USE.NAMES = F)
    }
    
    
    if(addGeneFoldChange == TRUE){
        if(is.null(fc_data)){ message("...add fold change data ..."); return('') }}
    
    message(glue('GSEA Enrichment analysis with \n {paste0(names(pwlist), collapse = " ")}'))
    
    
    cat('GO BP GSEA...\n')
    
    try(expr = {
        
        pathways[['GO BP']] <- clusterProfiler::GSEA(geneList = geneList,
                                    TERM2GENE = pwlist$GOBP.term2gene, 
                                    TERM2NAME = pwlist$GOBP.term2name,
                                    pvalueCutoff = pvalcutoff, 
                                    minGSSize = minGSSize, 
                                    maxGSSize = maxGSSize, 
                                    pAdjustMethod = 'BH', 
                                    by = 'fgsea',
                                    BPPARAM=SnowParam(workers = 20))

    })
    
    if(GO_CC_MF == TRUE){
        
        cat('GO CC enrichment analysis...\n')
        try(expr = {
            
            pathways[['GO CC']] <- GSEA(geneList = geneList,
                                        TERM2GENE = pwlist$GOCC.term2gene,
                                        TERM2NAME = pwlist$GOCC.term2name,
                                        pvalueCutoff = pvalcutoff, 
                                        minGSSize = minGSSize, 
                                        maxGSSize = maxGSSize, 
                                        pAdjustMethod = 'BH', 
                                        by = 'fgsea',
                                        BPPARAM=SnowParam(workers = 20))
            
        })
        
        cat('GO MF enrichment analysis...\n')
        try(expr = {
            
            pathways[['GO MF']] <- GSEA(geneList = geneList,
                                        TERM2GENE = pwlist$GOMF.term2gene,
                                        TERM2NAME = pwlist$GOMF.term2name,
                                        pvalueCutoff = pvalcutoff, 
                                        minGSSize = minGSSize, 
                                        maxGSSize = maxGSSize, 
                                        pAdjustMethod = 'BH', 
                                        by = 'fgsea',
                                        BPPARAM=SnowParam(workers = 20))
            
            
        })   
        
    }
    
    cat('KEGG enrichment analysis...\n')
    try(expr = {
        pathways[['KEGG']] <- GSEA(geneList = geneList,
                                   TERM2GENE = pwlist$KEGG.term2gene,
                                   TERM2NAME = pwlist$KEGG.term2name,
                                   pvalueCutoff = pvalcutoff, 
                                   minGSSize = minGSSize, 
                                   maxGSSize = maxGSSize, 
                                   pAdjustMethod = 'BH', 
                                   by = 'fgsea',
                                   BPPARAM=SnowParam(workers = 20))

    })
    
    
    if(only.go.kegg ==FALSE){
        try(expr = {
            cat('KEGG Module enrichment analysis...\n')
            pathways[['KEGG.M']] <- GSEA(geneList = geneList,
                                         TERM2GENE = pwlist$MKEGG.term2gene,
                                         TERM2NAME = pwlist$MKEGG.term2name,
                                         pvalueCutoff = pvalcutoff, 
                                         minGSSize = minGSSize, 
                                         maxGSSize = maxGSSize, 
                                         pAdjustMethod = 'BH', 
                                         by = 'fgsea',
                                         BPPARAM=SnowParam(workers = 20))

        })
        
        
        try(expr = {
            cat('WikiPathway enrichment analysis...\n')
            pathways[['WikiPathways']] <- GSEA(geneList = geneList,
                                               TERM2GENE = pwlist$wp.term2gene,
                                               TERM2NAME = pwlist$wp.term2name,
                                               pvalueCutoff = pvalcutoff, 
                                               minGSSize = minGSSize, 
                                               maxGSSize = maxGSSize, 
                                               pAdjustMethod = 'BH', 
                                               by = 'fgsea',
                                               BPPARAM=SnowParam(workers = 20))

        })
        
        
        try(expr = {
            cat('Reactome enrichment analysis...\n')
            pathways[['Reactome']] <- GSEA(geneList = geneList,
                                           TERM2GENE = pwlist$reactome.term2gene,
                                           TERM2NAME = pwlist$reactome.term2gene,
                                           pvalueCutoff = pvalcutoff, 
                                           minGSSize = minGSSize, 
                                           maxGSSize = maxGSSize, 
                                           pAdjustMethod = 'BH', 
                                           by = 'fgsea',
                                           BPPARAM=SnowParam(workers = 20))

        }
        )
        
        
        try(expr = {
            cat("MsigDB Hallmark enrichment analysis...\n")
            pathways[['MsigHall']] <- GSEA(geneList = geneList,
                                           TERM2GENE = pwlist$hallmark.term2gene,
                                           TERM2NAME = pwlist$hallmark.term2name,
                                           pvalueCutoff = pvalcutoff, 
                                           minGSSize = minGSSize, 
                                           maxGSSize = maxGSSize, 
                                           pAdjustMethod = 'BH', 
                                           by = 'fgsea',
                                           BPPARAM=SnowParam(workers = 20))

        })
    }
    
    modules <- names(pathways[!sapply(pathways, is.null)])
    
    if(addGeneFoldChange == TRUE){
        
        for(i in 1:length(modules)){
            try(expr = {
                pathways[[modules[i]]] <- addGeneFoldChange_GSEA(pathways[[modules[i]]], fc_data = fc_data)
                message(glue('{modules[i]} add gene fold change'))
            })
        }
    }else{
        
        try(expr = {
            pathways[[modules[i]]] <- setentrezToSymbol_GSEA(pathways[[modules[i]]], 
                                                             gene.annotation = gene.pw.annotation)
            message(glue('{modules[i]} set Entrez id to SYMBOL'))
        })
        
    }
    

    cat(glue('{paste0(modules, collapse =",")} loaded'), sep = '\n')
    
    if(make2df == TRUE){
        pathways <- gsea_to_df(pathways)
    }
    
    return(pathways)
    
}



addGeneFoldChange_GSEA <- function(enrichresult, fc_data = top_genes_filt, top_genes = 10){
    # add gene fold change data to GSEA result. 
    # x is enrichGo@result object output of clusterprofiler::gseGO
    # fc_data : dataframe contains fold change data of genes. columns : Symbol, logFC
    # top_genes : make a new column that only contains N number of top fold change genes in each GO term. 
    
    x <- enrichresult@result
    if(is.null(x)){
        message('no enrich result...')
        break
    }
    
    x$geneID <- x$core_enrichment
    x <- x %>% separate_rows(geneID, sep = "/") %>% 
        left_join(fc_data %>% dplyr::select(gene_id, SYMBOL,log2FoldChange) %>% distinct(gene_id, .keep_all = T), 
                  by = c('geneID'='gene_id')) %>% 
        group_by(ID) %>% 
        dplyr::arrange(desc(abs(log2FoldChange))) %>% 
        ungroup() %>% 
        dplyr::mutate(geneID_FC = paste0(SYMBOL, '(',round(2^(log2FoldChange),2),')')) %>% 
        dplyr::select(-log2FoldChange) %>% 
        group_by(ID) %>% 
        dplyr::mutate(geneID_FC = paste0(geneID_FC, collapse = '/'), 
                      geneID = paste0(geneID, collapse = '/'),
                      SYMBOL = paste0(SYMBOL, collapse = '/')) %>% 
        distinct() %>% 
        arrange(pvalue)
    
    
    # x$core_enrichment_FC <- 1
    # attributes(enrichresult)$top.genes <- list()
    # 
    # for(i in 1:length(x$core_enrichment)){
    #     
    #     gs <- stringr::str_split(string = x$core_enrichment[i], pattern = '/')[[1]]
    #     
    #     fcs <- fc_data[fc_data$SYMBOL %in% gs, c("log2FoldChange",'SYMBOL'), drop = FALSE] %>% as.data.frame() %>% 
    #         dplyr::arrange(desc(abs(log2FoldChange))) %>% distinct(SYMBOL, .keep_all = TRUE)
    #     
    #     g_fc <- paste0(fcs$SYMBOL, '(',round(2^(fcs$log2FoldChange),2),')')
    #     
    #     attributes(enrichresult)$top.genes[[x$ID[i]]] <- fcs$SYMBOL
    #     
    #     if(!is.null(top_genes)){
    #         x$geneID_FC_top_genes[i] <- paste0(g_fc[1:min(length(g_fc),top_genes)], collapse = '/')
    #     }
    #     
    #     x$core_enrichment_FC[i] <- paste0(g_fc, collapse = '/')
    #     x$core_enrichment[i] <- paste0(fcs$SYMBOL, collapse = '/')
    # }
    enrichresult@result <- x
    
    return(enrichresult)
}

setentrezToSymbol_GSEA <- function(enrichresult, gene.annotation = gene.pw.annotation){
    
    x <- enrichresult@result
    if(is.null(x)){
        message('no enrich result...')
        break
    }
    message(' set entrez ID to gene symbol using gene annotation file')
    
    x$geneID <- x$core_enrichment
    x <- x %>% separate_rows(geneID, sep = "/") %>% 
        left_join(gene.annotation %>% dplyr::select(entrez, SYMBOL) %>% distinct(), 
                  by = c('geneID'='entrez')) %>% 
        mutate(geneID = SYMBOL) %>% 
        dplyr::select(-SYMBOL) %>% 
        group_by(ID) %>% 
        mutate(geneID = paste0(geneID, collapse = "/")) %>% distinct()
    
    enrichresult@result <- x
    
    return(enrichresult)
    
}

gsea_to_df <- function(enrichresult){
    
    resultframe_list <- list()
    
    for(i in 1:length(enrichresult)){
        
        if(!is.null(enrichresult[[i]])){
            df <- enrichresult[[i]]@result
            
            if(nrow(df) !=0)
            {
                df$Description <- gsub(x = df$Description, pattern = "_", replacement = " ") 
                df$Description <- gsub(x = df$Description, pattern = " - Mus musculus \\(house mouse\\)", replacement = "") 
                df$Database <- names(enrichresult)[i]
                df <- df %>% dplyr::select(Database, ID, Description, everything())
                resultframe_list[[names(enrichresult[i])]] <- df
            }
        }
        
    }
    
    resultframe_df <- rbindlist(resultframe_list, fill = TRUE, use.names = TRUE) %>% as_tibble()
    
    return(resultframe_df)
    
}


enrichBarPlot <- function(enrich.df,
                          database = 'GO BP',
                          title = NULL,
                          annotation = TRUE,
                          alpha = 0.3, 
                          description_width = 30,
                          annotation_width = 80,
                          showCategory = 12, 
                          sorted = 'p.adjust', 
                          desc= FALSE, 
                          subset = NULL
){
    df <- enrich.df %>% dplyr::filter(Database == all_of(database))
    
    df$Description <- paste0(df$Count,'/',df$geneset_count, ' ',df$Description)
    total_deg <- unique(df$DEG_count)[1]
    xlab <- paste0('GeneRatio (Total DEG:', total_deg, ')')
    
    if(desc == TRUE){
        df <- df %>% arrange(desc(.[[sorted]])) %>% 
            slice_head(n = showCategory) %>% 
            mutate(Description = factor(Description, levels = Description))
        
    }else{
        df <- df %>% arrange(.[[sorted]]) %>% 
            slice_head(n = showCategory) %>% 
            mutate(Description = factor(Description, levels = Description))
        
    }
    if(!is.null(subset)){
        df <- df[subset,]
    }
    
    # Draw Gene ratio vs Enrich Term Plot 
    
    p <- ggplot(df, aes(x = GeneRatio, y = Description)) +
        geom_col(fill = 'darkblue', width = 0.8, alpha = alpha) +
        scale_y_discrete(label = function(x) stringr::str_wrap(x, width = description_width), limits = rev) +
        scale_x_continuous(expand = expansion(mult = c(0,0.1))) +
        xlab(xlab) +
        theme_classic() + 
        theme(axis.text.y = element_text(face='bold',size=11, color='black'), 
              axis.title.y = element_blank(),
              legend.position = 'none') +
        ggtitle(title)
    
    if(annotation == TRUE){
        tryCatch(expr={ 
            p <- p + 
                geom_text(aes(label=stringr::str_wrap(geneID_FC_top_genes, width = annotation_width), x= 0), 
                          hjust = -0.01, size=4, fontface='bold') 
        })
    }
    p
    return(p)
    
    
}


goBarPlot <- function(enrichresult, 
                      title='GO BP results', 
                      annotation = TRUE, alpha = 0.3, 
                      description_width = 30,
                      annotation_width = 80,
                      showCategory = 12, sorted = 'p.adjust', desc= FALSE, subset = NULL){
    
    df <- enrichresult@result
    df$geneset_count <- gsub(x=df$BgRatio,pattern = "/.*",replacement = "")
    df$GeneRatio <- sapply(df$GeneRatio, function(txt) eval(parse(text=txt)))
    df$Description <- gsub(x = df$Description, pattern = "_", replacement = " ") 
    df$Description <- paste0(df$Count,'/',df$geneset_count, ' ',df$Description)
    
    total_deg <- length(enrichresult@gene)
    xlab <- paste0('GeneRatio (Total DEG:', total_deg, ')')
    
    if(desc == TRUE){
        df <- df %>% arrange(desc(.[[sorted]])) %>% 
            slice_head(n = showCategory) %>% 
            mutate(Description = factor(Description, levels = Description))
        
    }else{
        df <- df %>% arrange(.[[sorted]]) %>% 
            slice_head(n = showCategory) %>% 
            mutate(Description = factor(Description, levels = Description))
        
    }
    
    if(!is.null(subset)){
        df <- df[subset,]
    }
    
    
    # Draw Gene ratio vs Enrich Term Plot 
    
    p <- ggplot(df, aes(x = GeneRatio, y = Description)) +
        geom_col(fill = 'darkblue', width = 0.8, alpha = alpha) +
        scale_y_discrete(label = function(x) stringr::str_wrap(x, width = description_width), limits = rev) +
        scale_x_continuous(expand = expansion(mult = c(0,0.1))) +
        xlab(xlab) +
        theme_classic() + 
        theme(axis.text.y = element_text(face='bold',size=11, color='black'), 
              axis.title.y = element_blank(),
              legend.position = 'none') +
        ggtitle(title)
    
    if(annotation == TRUE){
        tryCatch(expr={ 
            p <- p + 
                geom_text(aes(label=stringr::str_wrap(geneID_FC_top_genes, width = annotation_width), x= 0), 
                          hjust = -0.01, size=4, fontface='bold') 
        })
    }
    p
    return(p)
}


cluster.path.go.plot <- function(cluster.all, subcluster, 
                                 title ="HFD KO DOWN C1", 
                                 fc_data = res.hfd.ko.hfd.wt.table, pvalcutoff = 0.05, qvalcutoff = 0.1,
                                 only.go.kegg = TRUE) {
    
    
    subcluster.symbols <- cluster.all %>% filter(cluster %in% subcluster) %>% pull(SYMBOL) %>% unique()
    
    path <- getAllenrichPathways(subcluster.symbols, 
                                 qvalcutoff, 
                                 pvalcutoff, 
                                 addGeneFoldChange=TRUE, 
                                 fc_data = fc_data,
                                 simplifyByTopGenes = TRUE, 
                                 only.go.kegg = only.go.kegg )
    print(goBarPlot(path[['GO BP']], title = title, sorted = 'GeneRatio', desc = TRUE))
    print(clusterPlot2(cluster.all %>% filter(cluster %in% subcluster)))
    return(path)
} 


GOtoGeneSet <- function(goid, goname,
                        godb = pwlist.mouse$mouse.GOBP.term2gene, 
                        gonamedb = pwlist.mouse$mouse.GOBP.term2name, 
                        toSYMBOL = TRUE, 
                        gene.annotation = gene.annotation ){
    
    
    if(is.null(goid)){
        goid <-  gonamedb %>% dplyr::filter(name == goname) %>% pull(term) %>% unique()
    }
    
    geneset <- godb %>% left_join(gonamedb) %>% filter(term == goid) %>% distinct()
    
    if(toSYMBOL){
        
        geneset <- geneset %>% dplyr::rename(entrez = gene) %>% 
            left_join(gene.annotation, by = 'entrez') %>% distinct()
        
    }
    
    return(geneset)
    
}


get_intersect_pathway <- function(l, anti = FALSE){
    
    # l ; list of pathways df object by getAllenrichPathways-> get_intersect_pathway 
    # l[[i]] ; each of pathway eg. GO BP 
    intersect_pw <- list()
    pathways <- names(l[[1]])
    
    for(i in 1:length(pathways)){
        pathway_i <- pathways[i]
        
        pw_per_samples <- list()
        
        for(j in 1:length(l)){
            pw_per_samples[[names(l)[j]]] <- l[[j]][[i]]
            
        }
        
        intersect_pw[[pathway_i]] <- pw_per_samples
        
    }
    
    for(i in 1:length(intersect_pw)){
        
        if(anti){
            intersect_pw[[i]] <- purrr::reduce(intersect_pw[[i]], .f = anti_join, 
                                               by = c('ID','Description')
            )
            
        }else{
            colns <- colnames(intersect_pw[[i]][[1]])[-c(1,2)] # make new col names without ID,Description
            colns_w_suffix <- paste(rep(colns, each = length(names(l))), names(l), sep = "_")
            intersect_pw[[i]] <- purrr::reduce(intersect_pw[[i]], .f = inner_join, 
                                               by = c('ID','Description'), 
                                               suffix = paste0('_',names(l))
            )
            intersect_pw[[i]]  <- intersect_pw[[i]] %>% dplyr::select(ID,Description, colns_w_suffix)
        }
        
    }
    return(intersect_pw)
    
}

