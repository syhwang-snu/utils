

#### 2024-05-30 Microarrya DEG functions ####


deg.updn.list.stats <- function(lst){
    
    
    stats <- sapply(lst, nrow)
    compare.stat.df <- tibble(  compare_group = names(lst), 
                                stats = stats)
    compare.stat.df <- compare.stat.df %>% mutate(compare = gsub(pattern = '_UP$|_DN$', replacement = '', x = compare_group), 
                                                  direction = ifelse(grepl(pattern = 'UP$', x= compare_group), 'UP','DN'))
    compare.stat.df <- compare.stat.df %>% dplyr::select(compare, direction, stats) %>% 
        pivot_wider(id_cols = compare, names_from = 'direction', values_from = 'stats')
    
    return(compare.stat.df)
    
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


getMaxMedian <- function(compare.groups){
    
    compare.groups.samples <- meta.data %>% filter(group %in% compare.groups) %>% dplyr::select(sample_id, group)
    
    maxMedian <- normalized.counts %>% dplyr::select(-gene_id, - SYMBOL, -gene_name) %>% 
        pivot_longer(cols = -c('probe_id'), 
            names_to = 'sample_id', 
            values_to = 'normalized_counts') %>% 
        inner_join(compare.groups.samples) %>%  
        group_by(probe_id, group) %>% 
        summarise(median_counts = median(normalized_counts)) %>% 
        summarise(maxMedian = max(median_counts)) %>%
        ungroup() %>% dplyr::select(probe_id, maxMedian)
    
    return(maxMedian)
} 

getMaxMin <- function(compare.groups){
    
    compare.groups.samples <- meta.data %>% filter(group %in% compare.groups) %>% dplyr::select(sample_id, group)
    maxMin <- normalized.counts %>% dplyr::select(-gene_id, - SYMBOL, -gene_name) %>% 
        pivot_longer(cols = -c('probe_id'), 
                 names_to = 'sample_id', 
                 values_to = 'normalized_counts') %>% 
        inner_join(compare.groups.samples) %>%  
        group_by(probe_id, group) %>% 
        summarise(min_counts = min(normalized_counts)) %>% 
        summarise(maxMin = max(min_counts)) %>%
        ungroup() %>% dplyr::select(probe_id, maxMin)
    
    return(maxMin)
} 


getWilcox <- function(compare.groups){
    
    # Wilcox p calculation using all samples
    
    compare.groups.samples <- meta.data %>% 
        dplyr::filter(group %in% compare.groups) %>% 
        dplyr::select(sample_id, group)
    wilcox.p <- normalized.counts %>% dplyr::select(-gene_id, - SYMBOL, -gene_name) %>% 
        pivot_longer(cols = -c('probe_id'), 
                               names_to = 'sample_id', 
                               values_to = 'normalized.counts') %>% 
        inner_join(compare.groups.samples) %>%  
        dplyr::group_by(probe_id) %>% 
        dplyr::summarise(wilcox.p = wilcox.test(formula = normalized.counts ~ group, paired = FALSE)$p.value) %>% 
        ungroup() %>% dplyr::select(probe_id, wilcox.p) 
    
    return(wilcox.p)
} 



VennToTable <- function(lst, table, id.gene = as.character(paste0(1:length(lst), collapse = ''))){
    
    
    up.venn <- Venn(lapply(lst, `[[`, 'probe_id'))
    up.venn <- process_data(up.venn)
    up.venn.gene <- up.venn@region %>% dplyr::filter(id %in% id.gene) %>% pull(item)
    up.venn.gene <- unlist(up.venn.gene)
    venn.gene.table <- table %>% filter(probe_id %in% up.venn.gene)
    
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

genePlot <- function(count = normalized.counts.excel, 
                     symbol=NULL,
                     pb_id=NULL, 
                     wrap_plot=F,
                     color = 'group',
                     meta.data = meta.data,
                     x.angle=90, 
                     ncol= 1,
                     force_zero = TRUE,
                     stats = FALSE, 
                     title.size = 10, x.size = 10, y.size = 10, subtitle.size = 8,
                     saveplot = T,savename = NULL, savedir = './gene', save.x = 7, save.y = 9, save.format = 'png'){
    
    # Microarray count plots
    
    if(is.null(pb_id)){
        
        pb_ids <-  count %>% dplyr::filter(SYMBOL == symbol) %>%
            mutate(Amean = rowMeans(.[,unique(meta.data$group)])) %>%
            arrange(desc(Amean)) %>%
            pull(probe_id) %>% unique()
        
    }else{ 
        pb_ids <- pb_id
        symbol <- count[count$probe_id == pb_id, 'SYMBOL']
    }
    
    if(length(pb_ids) == 0){ print('No probe_ids'); return()}
    p <- list()
    
    for(i in 1:length(pb_ids)){
        
        pb_id <- pb_ids[i]
        
        geneCounts <- count %>% dplyr::filter(probe_id == pb_id) %>%
            dplyr::select(probe_id, SYMBOL, gene_name, all_of(meta.data$sample_id)) %>%
            pivot_longer(cols = meta.data$sample_id, names_to = 'sample_id', values_to = 'exp') %>%
            left_join(meta.data)
        
        gene_name <- geneCounts$gene_name[1]
        
        title <- paste0(symbol, "\n(", pb_id, ")" )
        subtitle <- gene_name
        
        cols = c('black','blue', 'red')
        seed <- as.numeric(charToRaw(title)) %>% sum()
        set.seed(seed)
        
        p[[pb_id]] <- ggbarplot(geneCounts,
                                x = "group",
                                y = "exp",
                                add = "mean_se",
                                title = title,
                                subtitle = subtitle,
                                color = "black",
                                fill = color,
                                alpha = 0.4,
                                palette = cols,
                                ylab = "Normalized Intensity") +
            geom_jitter(aes(fill = group), height = 0, width = 0.2, shape = 21) +
            scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
            #       scale_x_discrete(label = label) +
            theme(title = element_text(size=title.size, face='bold'),
                  plot.subtitle = element_text(size = subtitle.size, face='bold'),
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(size = x.size, face='bold', 
                                             angle = x.angle, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size = y.size, face='bold'),
                  axis.title.y = element_text(size = y.size, face='bold'),
                  plot.background = element_rect(fill = 'transparent', colour = 'transparent'),
                  panel.background = element_rect(fill = 'transparent', colour = 'transparent'),
                  legend.position = 'none')
        
        
        if(force_zero == TRUE){
            
            p[[pb_id]] <- p[[pb_id]] + expand_limits( y = 0)
            
        }
        
        if(stats == TRUE)
        {
            anova.p <- compare_means(count ~ group, geneCounts, method = "kruskal.test")
            anova.subtitle <- paste0(anova.p$method, ' p: ', anova.p$p.format)
            p[[pb_id]] <- p[[pb_id]] + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + # Add pairwise comparisons p-value
                labs(subtitle = anova.subtitle)  # Add global p-value
            
        }
        
        if(saveplot){
            if(!dir.exists(savedir)){
                dir.create(savedir)
            }
            
            savename <- paste0(symbol, '_', pb_id)
            ggsave( paste0(savename, '.', save.format), 
                    egg::set_panel_size(p[[pb_id]], width=unit(4, "cm"), height=unit(5, "cm")), 
                    width = save.x, height = save.y, units = 'cm', dpi = 300, 
                    path = savedir, device = save.format, bg = 'transparent')
            
        }
        
    }
    
    if(wrap_plot){
        p <- wrap_plots(p, byrow = T, ncol = ncol)
    }
    
    p
    
    return(p)
    
}
