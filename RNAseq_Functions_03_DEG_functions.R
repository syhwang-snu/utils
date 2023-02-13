

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#
# RNAseq functions 3 - DEG Functions #####
#

getresGeneTable <- function(dds= dds.all.clean, name=NULL, contrast=NULL, compare.groups, alpha = 0.1){
    
    res <- NULL
    if(!is.null(name)){
        res <- results(dds, 
                       name = name,
                       alpha = alpha,
                       filter = rowMeans(normal.count.mtx.all[,compare.groups])
                       
        )
    }
    if(!is.null(contrast)){
        res <- results(dds, 
                       contrast = contrast,
                       alpha = alpha,
                       filter = rowMeans(normal.count.mtx.all[,compare.groups])
                       
        ) 
        
    }
    
    
    res.df <- res %>% data.frame() %>% rownames_to_column(var="gene_id") %>% left_join(id2symbol.g29)
    res.df <- res.df %>% as_tibble() %>% left_join(normalized.counts.group) 
    
    res.df$maxExp <- rowMaxs(
        as.matrix(res.df[,compare.groups ]))
    
    res.df$maxMedian <- res.df %>% 
        left_join(getMaxMedian(compare.groups)) %>% 
        pull(maxMedian)
    
    res.df$maxMin <- res.df %>% 
        left_join(getMaxMin(compare.groups)) %>% 
        pull(maxMin)
    
    res.df$wilcox.p <- res.df %>% 
        left_join(getWilcox(compare.groups)) %>% 
        pull(wilcox.p)
    
    summary(res)
    return(res.df)
    
    
}




group.dge <- function(dds=dds.all.clean, 
                      group1, group2, 
                      all.res.tables = NULL,
                      padj_val = 0.1, 
                      log2FC = 1, 
                      maxExp = 10, 
                      maxMedian = 20, 
                      maxMin = 10, 
                      wilcox.p = 0.4, 
                      minMaxDiff = 5,
                      expression_cutoff=NULL ){
    
    if(is.null(all.res.tables)){
        res <- getresGeneTable(dds= dds, contrast = c('group',group1, group2), compare.groups = c(group1, group2))
    }else{
        
        res <- all.res.tables[[paste0(group1, "_vs_", group2)]]
    }
    
    sig <- res[res$padj < padj_val & 
                   abs(res$log2FoldChange) > log2FC & 
                   res$maxExp > maxExp &
                   res$maxMedian > maxMedian &
                   res$maxMin > maxMin & 
                   res$wilcox.p < wilcox.p  &
                   (res$maxMin * minMaxDiff > res$maxExp)
        , ] %>% 
        arrange(padj)
    if(!is.null(expression_cutoff)){
        passed.gene_id <- getExpressionCutoffGenes(select_groups = c(group1, group2), cutoff = expression_cutoff)
        sig <- sig %>% dplyr::filter(gene_id %in% passed.gene_id)
        print(nrow(sig))
    }

    return(sig)
    
}

getAllSigGeneTables <- function(combination_table = all.combination, 
                                dds = dds.all.clean,
                                all.res.tables = all.res.tables,
                                log2FoldChange = 1,
                                expression_cutoff = NULL,
                                padj = 0.1, 
                                maxExp = 10, 
                                maxMedian = 20, 
                                maxMin = 10, 
                                wilcox.p = 0.4, 
                                minMaxDiff = 5){
    # combination_table : group1, group2, compare ( group1_vs_group2)
    
    all.combination.sig.tables <- list()

    
    for(i in 1:nrow(combination_table)){
        cat(paste0(i, '/',nrow(combination_table),' ',all.combination$compare[i]), sep = '\n')
        
        suppressMessages(expr = { 
            all.combination.sig.tables[[combination_table$compare[i]]] <- 
                group.dge(dds = dds,
                          group1 = combination_table$group1[i], 
                          group2 = combination_table$group2[i],
                          all.res.tables = all.res.tables,
                          log2FC =   log2FoldChange,
                          padj_val =   padj,                                 
                          maxExp = maxExp, 
                          maxMedian = maxMedian, 
                          maxMin = maxMin, 
                          wilcox.p = wilcox.p, 
                          minMaxDiff = minMaxDiff,
                          expression_cutoff = expression_cutoff)
        })
    }

    
    return(all.combination.sig.tables)
    
}


getAllResGeneTables <- function(combination_table = all.combination,
                                dds = dds.all.clean,
                                expression_cutoff = NULL){
    # combination_table : group1, group2, compare ( group1_vs_group2)
    
    all.combination.res.tables <- list()
    
    for(i in 1:nrow(combination_table)){
        cat(paste0(i, '/',nrow(combination_table),' ',all.combination$compare[i]), sep = '\n')
        
        suppressMessages(expr = { 
            all.combination.res.tables[[combination_table$compare[i]]] <- 
                getresGeneTable(dds = dds, 
                                contrast = c("group", combination_table$group1[i], combination_table$group2[i]), 
                                compare.groups = c(combination_table$group1[i],combination_table$group2[i] ), alpha = 0.1)
        })
        
    }
    
    return(all.combination.res.tables)
}



group.specific.venn <- function(combination = all.combination, 
                                combination.tables.list = all.combination.sig.tables ,
                                group, 
                                others,
                                main_compare,
                                direction = 'UP'){
    
    
    if(direction == 'UP'){
        
        group1 <- combination %>% filter(group1 == group & group2 %in% others)
        group2 <- combination %>% filter(group2 == group & group1 %in% others)
        
        up.plus.compares <- group1$compare
        up.minus.compares <- group2$compare
        
        up.lists <- append(
            lapply(combination.tables.list[up.plus.compares], function(x){x %>% dplyr::filter(log2FoldChange > 0)}) ,
            lapply(combination.tables.list[up.minus.compares], function(x){x %>% dplyr::filter(log2FoldChange < 0)})
        )
        
        up.p <- ggVennDiagram(x = lapply(up.lists, `[[`, 'gene_id') , 
                              category.names = paste0(group ,' vs ' , c(group1$group2, group2$group1)))   + 
            scale_x_continuous(expand = expansion(0.2)) + 
            scale_fill_material('blue') + 
            theme(legend.position = "none") + 
            ggtitle(glue("{group} - UP Regulated genes"))
        
        
        id.gene <- as.character(paste0(1:length(others), collapse = ''))
        
        up.venn <- Venn(lapply(up.lists, `[[`, 'gene_id'))
        up.venn <- process_data(up.venn)
        up.venn.gene <- up.venn@region %>% dplyr::filter(id == id.gene) %>% pull(item) %>% .[[1]]
        
        print(up.p)
        venn.gene.table <- combination.tables.list[[main_compare]] %>% filter(gene_id %in% up.venn.gene)
        
    }else{
        
        group1 <- combination %>% filter(group1 == group & group2 %in% others)
        group2 <- combination %>% filter(group2 == group & group1 %in% others)
        
        dn.plus.compares <- group1$compare
        dn.minus.compares <- group2$compare
        
        dn.lists <- append(
            lapply(combination.tables.list[dn.plus.compares], function(x){x %>% dplyr::filter(log2FoldChange < 0)}) ,
            lapply(combination.tables.list[dn.minus.compares], function(x){x %>% dplyr::filter(log2FoldChange > 0)})
        )
        
        dn.p <- ggVennDiagram(x = lapply(dn.lists, `[[`, 'gene_id') , 
                              category.names = paste0(group ,' vs ' , c(group1$group2, group2$group1)))   + 
            scale_x_continuous(expand = expansion(0.2)) + 
            scale_fill_material('red') + 
            theme(legend.position = "none") + 
            ggtitle(glue("{group} - DOWN Regulated genes"))
        
        
        id.gene <- as.character(paste0(1:length(others), collapse = ''))
        
        dn.venn <- Venn(lapply(dn.lists, `[[`, 'gene_id'))
        dn.venn <- process_data(dn.venn)
        dn.venn.gene <- dn.venn@region %>% dplyr::filter(id == id.gene) %>% pull(item) %>% .[[1]]
        
        print(dn.p)
        venn.gene.table <- combination.tables.list[[main_compare]] %>% filter(gene_id %in% dn.venn.gene)
        
        
        
        
    }
    
    
    return(venn.gene.table)
    
}


getExpressionCutoffGenes <- function(select_groups, cutoff = 0.5){
    
    sample_ids <- meta.data %>% filter(group %in% select_groups) %>% dplyr::select(sample_id, group)
    
    normalized.counts.cut <- normalized.counts.unclean %>% 
        select(gene_id, SYMBOL, all_of(sample_ids$sample_id)) %>% 
        group_by(gene_id,SYMBOL) %>% 
        pivot_longer(cols = -c(SYMBOL,gene_id), names_to = 'sample_id',values_to = 'counts') %>% 
        left_join(sample_ids[,c("sample_id", "group")]) %>% ungroup() %>% group_by(group, gene_id, SYMBOL) %>% 
        summarise(median_counts = median(counts)) %>% ungroup() %>% group_by(group) %>% 
        mutate( median_cut = ifelse( median_counts > quantile(median_counts, probs = cutoff), TRUE, FALSE)) %>% 
        select(-median_counts) %>% 
        ungroup() %>% filter(median_cut == TRUE) %>% distinct(gene_id) %>% pull(gene_id)
    
    return(normalized.counts.cut)
    
    
}


# DEG Statistics 

dge.stat <- function(x, compare){
    
    if(is.character(x)) {x <- get(x) }
    
    
    up <- nrow(x %>% filter( log2FoldChange > 0))
    down <- nrow(x %>% filter( log2FoldChange < 0))
    
    return(data.frame(
        Compare = compare, 
        ALL = up + down,
        UP = up, DOWN = down))
}



dge.stat.cutoff <- function(x, cutoff = c('Antrum_Ctrl_GF','Antrum_Ctrl_SPF')){
    
    if(is.character(x)) {x <- get(x) }
    
    cutoff_genes <- normalized.counts.median.cut %>% dplyr::select(gene_id, {{cutoff}}) %>% 
        filter(.[[cutoff[1]]] == TRUE | .[[cutoff[2]]] == TRUE) %>% distinct(gene_id) %>% pull(gene_id)
    
    up <- nrow(x %>% filter( log2FoldChange > 0 & gene_id %in% cutoff_genes))
    down <- nrow(x %>% filter( log2FoldChange < 0 & gene_id %in% cutoff_genes))
    
    return(data.frame(
        Compare = paste0(cutoff[1], ' VS ',cutoff[2]), 
        UP = up, DOWN = down))
}

dge.cut <- function(x, cutoff = c('Antrum_Ctrl_GF','Antrum_Ctrl_SPF')){
    
    if(is.character(x)) {x <- get(x) }
    
    cutoff_genes <- normalized.counts.median.cut %>% dplyr::select(gene_id, {{cutoff}}) %>% 
        filter(.[[cutoff[1]]] == TRUE | .[[cutoff[2]]] == TRUE) %>% distinct(gene_id) %>% pull(gene_id)
    
    df <- x %>% dplyr::filter(gene_id %in% cutoff_genes)
    
    return(df)
}

VennToTable <- function(lst, table, id.gene = as.character(paste0(1:length(lst), collapse = ''))){
    
    
    up.venn <- Venn(lapply(lst, `[[`, 'gene_id'))
    up.venn <- process_data(up.venn)
    up.venn.gene <- up.venn@region %>% dplyr::filter(id %in% id.gene) %>% pull(item)
    up.venn.gene <- unlist(up.venn.gene)
    venn.gene.table <- table %>% filter(gene_id %in% up.venn.gene)
    
    return(venn.gene.table)
}


ggVennFromTable <- function(list, 
                            category.names =  paste0(" vs ",names(list) ),
                            title = "Venn", 
                            expand.x = 0.2, 
                            col = 'blue'){
    
    ggVennDiagram(x = lapply(list, `[[`, 'gene_id') , 
                  category.names = category.names)   + 
        scale_x_continuous(expand = expansion(expand.x)) + 
        ggsci::scale_fill_material(col) + 
        theme(legend.position = "none") + 
        ggtitle(title)
    
}


getMaxMedian <- function(compare.groups){
    compare.groups.samples <- meta.data %>% filter(group %in% compare.groups) %>% select(sample_id, group)
    maxMedian <- normalized.counts.unclean %>% pivot_longer(cols = -c('gene_id','SYMBOL'), 
                                                            names_to = 'sample_id', 
                                                            values_to = 'normalized_counts') %>% 
        inner_join(compare.groups.samples) %>%  
        group_by(gene_id,SYMBOL, group) %>% 
        summarise(median_counts = median(normalized_counts)) %>% 
        summarise(maxMedian = max(median_counts)) %>%
        ungroup() %>% select(gene_id, maxMedian)
    
    return(maxMedian)
} 

getMaxMin <- function(compare.groups){
    compare.groups.samples <- meta.data %>% filter(group %in% compare.groups) %>% select(sample_id, group)
    maxMin <- normalized.counts.unclean %>% pivot_longer(cols = -c('gene_id','SYMBOL'), 
                                                         names_to = 'sample_id', 
                                                         values_to = 'normalized_counts') %>% 
        inner_join(compare.groups.samples) %>%  
        group_by(gene_id,SYMBOL, group) %>% 
        summarise(min_counts = min(normalized_counts)) %>% 
        summarise(maxMin = max(min_counts)) %>%
        ungroup() %>% select(gene_id, maxMin)
    
    return(maxMin)
} 


getWilcox <- function(compare.groups){
    
    # Wilcox p calculation using all samples ( include pancreatic dirty samples )
    
    compare.groups.samples <- meta.data %>% 
        filter(group %in% compare.groups) %>% 
        select(sample_id, group)
    wilcox.p <- normalized.counts.unclean %>% pivot_longer(cols = -c('gene_id','SYMBOL'), 
                                                           names_to = 'sample_id', 
                                                           values_to = 'normalized_counts') %>% 
        inner_join(compare.groups.samples) %>%  
        group_by(gene_id,SYMBOL) %>% 
        summarise(wilcox.p = wilcox.test(normalized_counts ~ group, paired = FALSE)$p.value) %>% 
        ungroup() %>% select(gene_id, wilcox.p) 
    
    return(wilcox.p)
} 




