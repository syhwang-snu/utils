library(msigdbr)

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)



getGOenrichPathways <- function(symbols, qvalcutoff = 0.5, pvalueCutoff = 0.5){
    
    entrez <- mapIds(org.Mm.eg.db,keys = symbols, column = 'ENTREZID', keytype = 'SYMBOL')
    pathways <- list()
    cat('GO BP enrichment analysis...\n')
    pathways[['GO BP']] <- enrichGO(gene = symbols, 
                                    OrgDb = org.Mm.eg.db, 
                                    keyType = 'SYMBOL', 
                                    ont = "BP", 
                                    universe = normalized_counts$SYMBOL, 
                                    qvalueCutoff = qvalcutoff,pvalueCutoff =  pvalueCutoff,
                                    pAdjustMethod = "none" )
    pathways[['GO BP']] <- simplify(pathways[['GO BP']])
    
    return(pathways)
    
}


getAllenrichPathways <- function(symbols, qvalcutoff = 0.05, pvalueCutoff = 0.05, 
                                 addGeneFoldChange=TRUE, fc_data = NULL,
                                 simplifyByTopGenes = TRUE, only.go.kegg = TRUE ){

    entrez <- mapIds(org.Mm.eg.db,keys = symbols, column = 'ENTREZID', keytype = 'SYMBOL')
    pathways <- list()
    
    
    if(addGeneFoldChange == TRUE){
        if(is.null(fc_data)){ 
            message("...add fold change data ...") 
            return('')
            }
        
        }
    
    cat('GO BP enrichment analysis...\n')
    try(expr = {
    
    pathways[['GO BP']] <- enrichGO(gene = symbols, 
                                                    OrgDb = org.Mm.eg.db, 
                                                    keyType = 'SYMBOL', 
                                                    ont = "BP", 
                                                    universe = normalized_counts$SYMBOL, 
                                                    qvalueCutoff = qvalcutoff,pvalueCutoff =  pvalueCutoff,
                                                    pAdjustMethod = "BH" )
    pathways[['GO BP']] <- simplify(pathways[['GO BP']])
    })
    
    
    cat('KEGG enrichment analysis...\n')
    
    try(expr = {
    pathways[['KEGG']] <- enrichKEGG(gene = entrez, 
                                                     organism = "mmu",
                                                     universe = all.entrez.id,
                                                     qvalueCutoff = 0.05 )
    pathways[['KEGG']] <- setReadable(OrgDb = org.Mm.eg.db, x = pathways[['KEGG']], keyType = 'ENTREZID')
    })
    cat('KEGG Module enrichment analysis...\n')
    
    try(expr = {
    pathways[['KEGG.M']] <- enrichMKEGG(gene = entrez, 
                                                        organism = "mmu",
                                                        universe = all.entrez.id,
                                                        qvalueCutoff = 0.2)
    pathways[['KEGG.M']]  <- setReadable(OrgDb = org.Mm.eg.db, x = pathways[['KEGG.M']] , keyType = 'ENTREZID')
    })
    
    cat('WikiPathway enrichment analysis...\n')
    try(expr = {
    pathways[['WikiPathways']] <- enrichWP(gene = entrez,
             organism = "Mus musculus",
             universe = all.entrez.id,
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.1 )
    pathways[['WikiPathways']] <- setReadable(OrgDb = org.Mm.eg.db, x = pathways[['WikiPathways']], keyType = 'ENTREZID')
    })
    cat('Reactome enrichment analysis...\n')
    
    try(expr = {
    pathways[['Reactome']] <- enrichPathway(gene = entrez, 
                                      organism = "mouse",
                                      universe = all.entrez.id,
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.1 )
    pathways[['Reactome']] <- setReadable(OrgDb = org.Mm.eg.db, x = pathways[['Reactome']], keyType = 'ENTREZID')
    })
    
    cat("MsigDB Hallmark enrichment analysis...\n")
    try(expr = {
    pathways[['MsigHall']] <- enricher(gene = entrez,
                                       universe = all.entrez.id,
                                        pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.1, 
                                       TERM2GENE = m_t2g
                                       
                                       )
    pathways[['MsigHall']] <- setReadable(OrgDb = org.Mm.eg.db, x = pathways[['MsigHall']], keyType = 'ENTREZID')
    })
    
    
    modules <- names(pathways[!sapply(pathways, is.null)])
    if(addGeneFoldChange == TRUE){
        
        for(i in 1:length(modules)){
            try(expr = {
            pathways[[modules[i]]] <- addGeneFoldChange(pathways[[modules[i]]], fc_data = fc_data)
            message(glue('{modules[i]} add gene fold change'))
            })
        }
    }
    
    if(simplifyByTopGenes == TRUE){
        
        for(i in 1:length(modules)){
            try(expr = {
            pathways[[modules[i]]] <- simplifyByTopGenes(pathways[[modules[i]]])
            message(glue('{modules[i]} simplified '))
            })
        }
        
    }
    cat(glue('{paste0(modules, collapse =",")} loaded'))
    return(pathways)
        
    }
    
    


addGeneFoldChange <- function(enrichresult, fc_data = top_genes_filt, top_genes = 10){
    # add gene fold change data to GO Enrichment analysis result. 
    # x is enrichGo@result object output of clusterprofiler::enrichGO
    # fc_data : dataframe contains fold change data of genes. columns :  Symbol, logFC
    # top_genes : make a new column that only contains N number of top fold change genes in each GO term. 
    x <- data.frame(enrichresult)
    if(is.null(x)){
        message('no enrich result...')
        break
    }
    x$geneID_FC <- 1
    attributes(enrichresult)$top.genes <- list()
    for(i in 1:length(x$geneID)){
        gs <- stringr::str_split(string = x$geneID[i], pattern = '/')[[1]]
        fcs <- fc_data[fc_data$SYMBOL %in% gs, c("log2FoldChange",'SYMBOL'), drop = FALSE] %>% as.data.frame() %>% 
            dplyr::arrange(desc(abs(log2FoldChange))) %>% distinct(SYMBOL, .keep_all = TRUE)
        g_fc <- paste0(fcs$SYMBOL, '(',round(fcs$log2FoldChange,2),')')
        
        attributes(enrichresult)$top.genes[[x$ID[i]]] <- fcs$SYMBOL
        
        if(!is.null(top_genes)){
            x$geneID_FC_top_genes[i] <- paste0(g_fc[1:min(length(g_fc),top_genes)], collapse = '/')
        }
        
        x$geneID_FC[i] <- paste0(g_fc, collapse = '/')
        
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


goBarPlot <- function(enrichresult, title='GO BP results', 
                      annotation = TRUE, alpha = 0.3, 
                      description_width = 30,
                      annotation_width = 80,
                      showCategory = 12){
    
    df <- enrichresult@result
    df$GeneRatio <- sapply(df$GeneRatio, function(txt) eval(parse(text=txt)))
    df <- df %>% arrange(p.adjust) %>% slice_head(n = showCategory) %>% mutate(Description = factor(Description, levels = Description))
    
    p <- ggplot(df, aes(x = GeneRatio, y = Description)) +
        geom_col(fill = 'darkblue', width = 0.8, alpha = alpha) +
        scale_y_discrete(label = function(x) stringr::str_wrap(x, width = description_width), limits = rev) +
        scale_x_continuous(expand = expansion(mult = c(0,0.1))) +
        theme_classic() + 
        theme(axis.text.y = element_text(face='bold',size=11, color='black'), 
              axis.title.y = element_blank(),
              legend.position = 'none') +
        ggtitle(title)

    if(annotation == TRUE){
        p <- p + 
            geom_text(aes(label=stringr::str_wrap(geneID_FC_top_genes, width = annotation_width), x= 0), 
                      hjust = -0.01, size=4, fontface='bold') 
    }
    p
    return(p)
}


cluster.path.go.plot <- function(cluster.all, subcluster, title ="HFD KO DOWN C1", fc_data = res.hfd.ko.hfd.wt.table) {
    
    path <-  getGOenrichPathways(cluster.all %>% filter(cluster %in% subcluster) %>% pull(SYMBOL) %>% unique(), qvalcutoff = 0.7)
    path[['GO BP']] <- addGeneFoldChange(path[['GO BP']], fc_data = fc_data)
    print(goBarPlot(simplifyByTopGenes( path[['GO BP']]), title = title))
    print(clusterPlot2(cluster.all %>% filter(cluster %in% subcluster)))
    
} 


GOtoGeneSet <- function(genesetName, enrichresult, only.topgene = TRUE){
    
    goid <- enrichresult@result %>% filter(Description %in% genesetName) %>% pull(ID)
    print(goid)
    geneset <- enrichresult@geneSets[names(enrichresult@geneSets) == goid][[1]]
    
    if(only.topgene == TRUE){
        top_genes <- attributes(enrichresult)$top.genes[names(attributes(enrichresult)$top.genes) == goid ][[1]]
        geneset <- top_genes[top_genes %in% geneset]
        
    }
    
    return(geneset)
    
}


