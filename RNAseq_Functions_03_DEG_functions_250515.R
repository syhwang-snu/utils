

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#
# RNAseq functions 3 - DEG Functions #####
# 
# 2025-05-15
# Expression Cutoff 먼저 적용 후 DEseq pvalue 및 padjust value 계산
# 
# 2025-04-07
# Add volcanoplot
# 
# 
# 
# 2025-01-05 fix 
# 2024-12-30 
# - fix getWilcox Error :  cannot use 'paired' in formula method

library(ggVennDiagram)
library(ggvenn)
library(ggrepel)

getresGeneTable <- function(dds= dds, name=NULL, contrast=NULL, 
                            compare.groups, 
                            alpha = 0.1, 
                            lfcThreshold = 0,
                            normalized.counts.group){
    
    res <- NULL
    if(!is.null(name)){
        res <- results(dds, 
                       name = name,
                       alpha = alpha, 
                       lfcThreshold = lfcThreshold )
                       
        
    }
    if(!is.null(contrast)){
        print(contrast)
        res <- results(dds, 
                       contrast = contrast,
                       alpha = alpha,
                       lfcThreshold = lfcThreshold)
        
    }
    
    
    res.df <- res %>% data.frame() %>% rownames_to_column(var="gene_id") %>% left_join(gene.annotation)
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

getAllResGeneTables <- function(combination_table = compare.combination,
                                normalized.counts.group = normalized.counts.group,
                                dds = dds, 
                                lfcThreshold= 0,
                                alpha = 0.1){
    # combination_table : group1, group2, compare ( group1_vs_group2)
    
    all.combination.res.tables <- list()
    
    for(i in 1:nrow(combination_table)){
        cat(paste0(combination_table$compare[i], ' ', i, '/',nrow(combination_table),' ',combination_table$compare[i]), sep = '\n')
        
        suppressMessages(expr = { 
            all.combination.res.tables[[combination_table$compare[i]]] <- 
                getresGeneTable(dds = dds, 
                                contrast = c("group", combination_table$group1[i], combination_table$group2[i]), 
                                compare.groups = c(combination_table$group1[i],combination_table$group2[i] ), 
                                alpha = alpha, 
                                lfcThreshold= lfcThreshold,
                                normalized.counts.group = normalized.counts.group)
        })
        
    }
    
    return(all.combination.res.tables)
}

ResTablesToUpDown <- function(list){
    
    sigList <- list()
    
    for(i in 1:length(list)){
        
        t <- list[[i]]
        up_name <- paste0(names(list)[i], "_UP")
        dn_name <- paste0(names(list)[i], "_DN")
        
        sigList[[up_name]] <- t %>% dplyr::filter(log2FoldChange > 0)
        sigList[[dn_name]] <- t %>% dplyr::filter(log2FoldChange < 0) 
        
        
    }
    
    return(sigList)
}


# Get DEG 
getGroupDEGs <- function(dds,
                         group1, group2,
                         all.res.tables = NULL,
                         lfcThreshold = 0,
                         padj_val = NULL,
                         pvalue = NULL,
                         log2FC = NULL,
                         maxExp = NULL,
                         maxMedian = NULL,
                         maxMin = NULL,
                         wilcox.p = NULL,
                         minMaxDiff = NULL,
                         expression_cutoff = NULL,
                         normalized.counts = NULL) {
    

    # 1. Expression cutoff 먼저 적용

    if (!is.null(expression_cutoff)) {
        passed.gene_id <- getExpressionCutoffGenes(
            select_groups = c(group1, group2),
            cutoff = expression_cutoff,
            normalized.counts = normalized.counts
        )
        
        # DESeqDataSet에서 필터링
        dds <- dds[rownames(dds) %in% passed.gene_id, ]
        
        # 2. 필터링된 dds로 DE 분석
        
        res <- getresGeneTable(
            dds = dds,
            contrast = c("group", group1, group2),
            compare.groups = c(group1, group2),
            lfcThreshold = lfcThreshold,  
            normalized.counts.group = normalized.counts
        )
        
    }else{
        if(!is.null(all.res.tables))
        {
            res <- all.res.tables[[paste0(group1, "_vs_", group2)]]
        }else{
            errorCondition(message = 'Give all.res.tables')
        }

    }
    

    # 3. 나머지 필터 조건들 적용

    filter_idx <- res$gene_id > 0 &
        (if (!is.null(pvalue))     { res$pvalue < pvalue }       else TRUE) &
        (if (!is.null(padj_val))   { res$padj < padj_val }       else TRUE) &
        (if (!is.null(log2FC))     { abs(res$log2FoldChange) > log2FC } else TRUE) &
        (if (!is.null(maxExp))     { res$maxExp > maxExp }       else TRUE) &
        (if (!is.null(maxMedian))  { res$maxMedian > maxMedian } else TRUE) &
        (if (!is.null(maxMin))     { res$maxMin > maxMin }       else TRUE) &
        (if (!is.null(wilcox.p))   { res$wilcox.p < wilcox.p }   else TRUE) &
        (if (!is.null(minMaxDiff)) { res$maxMin * minMaxDiff > res$maxExp } else TRUE)
    
    filter_idx[is.na(filter_idx)] <- FALSE
    sig <- res[filter_idx, ] %>% arrange(padj)
    
    print(nrow(sig))
    return(sig)
}


getExpressionCutoffGenes <- function(select_groups, cutoff = 0.5, normalized.counts = normalized.counts){
    
    sample_ids <- meta.data %>% filter(group %in% select_groups) %>% dplyr::select(sample_id, group)
    
    normalized.counts.cut <- normalized.counts %>% 
        dplyr::select(gene_id, all_of(sample_ids$sample_id)) %>% 
        group_by(gene_id) %>% 
        pivot_longer(cols = -c(gene_id), names_to = 'sample_id',values_to = 'counts') %>% 
        left_join(sample_ids[,c("sample_id", "group")]) %>% ungroup() %>% group_by(group, gene_id) %>% 
        summarise(median_counts = median(counts)) %>% ungroup() %>% group_by(group) %>% 
        mutate( median_cut = ifelse( median_counts >= quantile(median_counts, probs = cutoff), TRUE, FALSE)) %>% 
        dplyr::select(-median_counts) %>% 
        ungroup() %>% filter(median_cut == TRUE) %>% distinct(gene_id) %>% pull(gene_id)
    
    return(normalized.counts.cut)
    
}


getAllSigGeneTables <- function(combination_table = all.combination, 
                                dds = dds,
                                all.res.tables = all.res.tables,
                                lfcThreshold = 0,
                                log2FoldChange = NULL,
                                expression_cutoff = NULL,
                                pvalue = NULL,
                                padj = NULL, 
                                maxExp = NULL, 
                                maxMedian = NULL, 
                                maxMin = NULL, 
                                wilcox.p = NULL, 
                                minMaxDiff = NULL, 
                                normalized.counts = normalized.counts){
    # combination_table : group1, group2, compare ( group1_vs_group2)
    
    all.combination.sig.tables <- list()

    
    for(i in 1:nrow(combination_table)){
        cat(paste0(i, '/',nrow(combination_table),' ',combination_table$compare[i]), sep = '\n')
        
        suppressMessages(expr = { 
            all.combination.sig.tables[[combination_table$compare[i]]] <- 
                getGroupDEGs(dds = dds,
                          group1 = combination_table$group1[i], 
                          group2 = combination_table$group2[i],
                          all.res.tables = all.res.tables,
                          lfcThreshold = lfcThreshold,
                          log2FC =   log2FoldChange,
                          pvalue = pvalue, 
                          padj_val =   padj,                                 
                          maxExp = maxExp, 
                          maxMedian = maxMedian, 
                          maxMin = maxMin, 
                          wilcox.p = wilcox.p, 
                          minMaxDiff = minMaxDiff,
                          expression_cutoff = expression_cutoff, 
                          normalized.counts = normalized.counts)
        })
    }

    
    return(all.combination.sig.tables)
    
}


# DEG Statistics 

# dge.sta%>% <- function(x, compare){
#     
#     if(is.character(x)) {x <- get(x) }
#     
#     
#     up <- nrow(x %>% filter( log2FoldChange > 0))
#     down <- nrow(x %>% filter( log2FoldChange < 0))
#     
#     return(data.frame(
#         Compare = compare, 
#         ALL = up + down,
#         UP = up, DOWN = down))
# }


deg.updn.list.stats <- function(lst, combination.table = compare.combination){
    
    
    stats <- sapply(lst, nrow)
    compare.stat.df <- tibble(  compare_group = names(lst), 
                                stats = stats)
    compare.stat.df <- compare.stat.df %>% mutate(compare = gsub(pattern = '_UP$|_DN$', replacement = '', x = compare_group), 
                                                  direction = ifelse(grepl(pattern = 'UP$', x= compare_group), 'UP','DN'))
    compare.stat.df <- compare.stat.df %>% dplyr::select(compare, direction, stats) %>% 
        pivot_wider(id_cols = compare, names_from = 'direction', values_from = 'stats')
   
    return(compare.stat.df)
    
}

#### Venn Diagram Compare ####


group.specific.venn <- function(combination = compare.combination, 
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




VennToTable <- function(lst, table, 
                        id.gene = as.character(paste0(1:length(lst), collapse = '/')),
                        id_type = 'gene_id'){
    
    # 2024-12-30 changed function
    # 2025-04-07 changed : id_type
    
    up.venn <- Venn(lapply(lst, `[[`, id_type))
    up.venn <- process_data(up.venn)
    # up.venn.gene <- up.venn@region %>% dplyr::filter(id %in% id.gene) %>% pull(item)
    up.venn.gene <- up.venn$regionData %>% dplyr::filter(id %in% id.gene) %>% pull(item)
    up.venn.gene <- unlist(up.venn.gene)
    venn.gene.table <- table %>% filter(.[[id_type]] %in% up.venn.gene) %>% arrange(pvalue)
    
    return(venn.gene.table)
}


ggVennFromTable <- function(lst, 
                            category.names =  names(lst),
                            title = "Venn", 
                            expand.x = 0.2, 
                            col = NULL){
    
    # ggVennDiagram(x = lapply(lst, `[[`, 'gene_id') , 
    #               category.names = category.names, set_color = "black")   + 
    #     scale_x_continuous(expand = expansion(expand.x)) + 
    #     scale_color_manual(values = rep('black', length(lst))) +
    #     ggsci::scale_fill_material(col) + 
    #     theme(legend.position = "none") + 
    #     ggtitle(title)
    
    if(is.null(col)){ col <- rep('transparent', length(lst))}
    
    ggvenn(data = lapply(lst, `[[`, 'gene_id'), show_percentage = F, text_size = 8,
           set_name_size = 5,auto_scale=FALSE, fill_color = col) +
        theme(
            plot.background = element_rect(fill = 'transparent', colour = 'transparent'),
            panel.background = element_rect(fill = 'transparent', colour = 'transparent'),
            legend.position = 'none'
        )
    
}

getMaxMedian <- function(compare.groups){
    compare.groups.samples <- meta.data %>% filter(group %in% compare.groups) %>% dplyr::select(sample_id, group)
    maxMedian <- normalized.counts %>% pivot_longer(cols = -c('gene_id','SYMBOL'), 
                                                            names_to = 'sample_id', 
                                                            values_to = 'normalized_counts') %>% 
        inner_join(compare.groups.samples) %>%  
        group_by(gene_id,SYMBOL, group) %>% 
        summarise(median_counts = median(normalized_counts)) %>% 
        summarise(maxMedian = max(median_counts)) %>%
        ungroup() %>% dplyr::select(gene_id, maxMedian)
    
    return(maxMedian)
} 

getMaxMin <- function(compare.groups){
    compare.groups.samples <- meta.data %>% filter(group %in% compare.groups) %>% dplyr::select(sample_id, group)
    maxMin <- normalized.counts %>% pivot_longer(cols = -c('gene_id','SYMBOL'), 
                                                         names_to = 'sample_id', 
                                                         values_to = 'normalized_counts') %>% 
        inner_join(compare.groups.samples) %>%  
        group_by(gene_id,SYMBOL, group) %>% 
        summarise(min_counts = min(normalized_counts)) %>% 
        summarise(maxMin = max(min_counts)) %>%
        ungroup() %>% dplyr::select(gene_id, maxMin)
    
    return(maxMin)
} 


getWilcox <- function(compare.groups){
    
    # Wilcox p calculation using all samples
    
    compare.groups.samples <- meta.data %>% 
        dplyr::filter(group %in% compare.groups) %>% 
        dplyr::select(sample_id, group)
    wilcox.p <- normalized.counts %>% pivot_longer(cols = -c('gene_id','SYMBOL'), 
                                                           names_to = 'sample_id', 
                                                           values_to = 'normalized.counts') %>% 
        inner_join(compare.groups.samples) %>%  
        dplyr::group_by(gene_id,SYMBOL) %>% 
        dplyr::summarise(wilcox.p = wilcox.test(formula = normalized.counts ~ group)$p.value) %>% 
        ungroup() %>% dplyr::select(gene_id, wilcox.p) 
    
    return(wilcox.p)
} 



restoSigNcuts <- function(list, Ncutoff = 300){
    
    sigList = list()
    
    for(i in 1:length(list)){
        
        t <- list[[i]]
        up_name <- paste0(names(list)[i], "_UP")
        dn_name <- paste0(names(list)[i], "_DN")
        pass_geneid = normalized.counts.q1.cut %>%  
            filter(.[[all.combination$group1[i]]] == TRUE | .[[all.combination$group2[i]]] == TRUE) %>% pull(gene_id)
        sigList[[up_name]] <- t %>% filter(log2FoldChange > 0 & gene_id %in% pass_geneid) %>% slice_head(n= Ncutoff)
        sigList[[dn_name]] <- t %>% filter(log2FoldChange < 0 & gene_id %in% pass_geneid) %>% slice_head(n= Ncutoff)
        
        
    }
    
    return(sigList)
}


ResTablesToUpDown <- function(list){
    
    sigList = list()
    
    for(i in 1:length(list)){
        
        t <- list[[i]]
        up_name <- paste0(names(list)[i], "_UP")
        dn_name <- paste0(names(list)[i], "_DN")
        
        sigList[[up_name]] <- t %>% filter(log2FoldChange > 0)
        sigList[[dn_name]] <- t %>% filter(log2FoldChange < 0) 
        
        
    }
    
    return(sigList)
}






