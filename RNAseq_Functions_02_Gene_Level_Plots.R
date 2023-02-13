library(patchwork)
library(cowplot)
library(ggpubr)
library(ggh4x)


genePlot <- function(symbol=NULL,gene_id=NULL,dds=dds, stats = TRUE, 
                     group = NULL, 
                     title = NULL, 
                     my_comparisons = list( c("Antrum_Ctrl_SPF", "Antrum_Ctrl_GF")), 
                     force_zero = TRUE, x.angle = 90)
{
    
    if(is.null(gene_id)){
        g2id <- id2symbol.g29$gene_id[id2symbol.g29$SYMBOL == symbol][1]
        
        geneCounts <- plotCounts(dds = dds, gene = g2id, intgroup = "group", normalized = TRUE,
                                 returnData = TRUE, transform =FALSE)
        
    }else{
        g2id <- gene_id
        geneCounts <- plotCounts(dds = dds, gene = g2id, intgroup = "group", normalized = TRUE,
                                 returnData = TRUE, transform =FALSE)
        symbol <- id2symbol.g29$SYMBOL[id2symbol.g29$gene_id == gene_id][1]
    }
    
    
    gene_name <- symbol2genename$gene_name[symbol2genename$SYMBOL == symbol]
    if(is.null(title)){
        title <- paste0(symbol, " Expression")
        subtitle <- gene_name
    }
    
    if(!is.null(group)){
        geneCounts$group <- factor(geneCounts$group, levels =  group)
    }
    
    p <- ggboxplot(geneCounts, x = "group", y = "count", 
                   add = "jitter", 
                   title = title, subtitle = subtitle,
                   ylab = "Normalized Counts", 
                   width = 0.7, 
                   lwd = 1.2, 
                   fatten = 3) + 
        theme(title = element_text(size=20, face='bold'), 
              plot.subtitle = element_text(size=15, face='bold'),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 18, face='bold', angle = x.angle), 
              axis.text.y = element_text(size = 18, face='bold'), 
              axis.title.y = element_text(size = 18, face='bold')) 
    
    
    if(force_zero == TRUE){
        
        p <- p + expand_limits( y = 0) 
        
    }
    
    
    if(stats == TRUE)
    {
        anova.p <- compare_means(count ~ group, geneCounts, method = "kruskal.test")
        anova.subtitle <- paste0(anova.p$method, ' p: ', anova.p$p.format)
        p <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + # Add pairwise comparisons p-value
            labs(subtitle = anova.subtitle)  # Add global p-value
        
    }
    
    
    p
    
}


ExpCompPlot <- function(table = res.hfd.ko.chow.ko.table, x= 'CHOW_KO', y = 'CDAHFD_KO', sig.gene.id=NULL, label=NULL, title = 'HFD W2KO vs CHOW W2KO'){
    
    
    table$sig.gene <- ifelse(table$gene_id %in% sig.gene.id,1,0)
    table$sig.gene <- factor(table$sig.gene, levels = c(0,1))
    table$label <- ''
    
    if(!is.null(label)){
        table$label <- ifelse(table$SYMBOL %in% label,table$SYMBOL,'')
        
    }
    
    p <- ggplot(table, aes(x= .data[[x]], y=.data[[y]], label = SYMBOL)) +
        geom_point(aes(color = sig.gene)) + 
        scale_color_manual(values = c('grey','red')) + 
        geom_label_repel(max.overlaps = 1000, aes(label= label)) + 
        scale_x_log10(labels = scales::label_number()) + 
        scale_y_log10(labels = scales::label_number()) +
        stat_regline_equation(formula = y ~ x, aes(label = ..rr.label..), 
                              label.x.npc = 0.1, label.y.npc = 0.99, 
                              size= 4, fontface='bold')+
        coord_cartesian(xlim = c(0.1, 100000), ylim=c(0.1, 100000)) +
        theme_bw() + theme(legend.position = 'none') +
        labs(title = title)
    
    
    p
    
    
}


clusterPlot <- function(gene_ids, h= 4, cluster_k= 4, group=NULL, interaction = NULL, ncol=NULL, 
                        normalized_counts.group=normalized_counts.group ){
    
    
    normal.count.mtx.all <- normalized_counts.group %>% 
        dplyr::select(-SYMBOL) %>%
        column_to_rownames('gene_id') %>% 
        as.matrix()
    normal.count.mtx <- normal.count.mtx.all[gene_ids,group]
    normal.count.mtx <- normal.count.mtx %>% t() %>% scale() %>% t()
    gene_dist <- dist(normal.count.mtx)
    normal.count.mtx.df <- as.data.frame(normal.count.mtx) %>% rownames_to_column('gene_id')
    
    
    gene_hclust <- hclust(gene_dist, method = "ward.D2")
    # The default `plot()` function can be used to produce a simple dendrogram
    #p1 <- plot(gene_hclust, labels = FALSE)
    #p1 <- p1 + abline(h = 4, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
    
    gene_cluster <- cutree(gene_hclust, k = cluster_k) %>% 
        # turn the named vector into a tibble
        enframe() %>% 
        # rename some of the columns
        rename(gene_id = name, cluster = value)
    
    trans_cts_cluster <- normal.count.mtx.df %>% 
        inner_join(gene_cluster, by = "gene_id") %>% group_by(cluster) %>% mutate(gene_n=n()) %>% 
        pivot_longer(cols = c(-gene_id,-cluster,-gene_n), names_to = 'group',values_to = 'cts') %>% 
        ungroup() %>% mutate(cluster_label = paste0(cluster, " (",gene_n," genes)")) %>% left_join(id2symbol.g29)
    
    trans_cts_cluster$group <- factor(trans_cts_cluster$group, levels = group)
    
    
    p<- trans_cts_cluster %>% 
        ggplot(aes(group, cts)) +
        geom_line(aes(group = gene_id), alpha = 0.3) +
        geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1.5, 
                  aes(group = 1)) + 
        geom_hline(yintercept = 0) +
        facet_wrap(vars(cluster_label)) +
        ylab("z-expression") + 
        theme_bw() + 
        theme(text = element_text(face='bold', size = 15), 
              axis.title.x = element_blank(),axis.text.x = element_text(angle = 90))
    
    
    print(p)
    return(trans_cts_cluster)
    
}


clusterPlot.only <- function(trans_cts_cluster){
    
    p<- trans_cts_cluster %>% 
        ggplot(aes(group, cts)) +
        geom_line(aes(group = gene_id), alpha = 0.3) +
        geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1.5, 
                  aes(group = 1)) + 
        geom_hline(yintercept = 0) +
        facet_wrap(vars(cluster_label)) +
        ylab("z-expression") + 
        theme_bw() + 
        theme(text = element_text(face='bold', size = 15), 
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90))
    return(p)
    
    
}

clusterPlot2 <-function(cluster_df){
    
    cluster_df %>% ggplot(aes(group, cts)) +
        geom_line(aes(group = gene_id), alpha = 0.3) +
        geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1.5, 
                  aes(group = 1)) + 
        geom_hline(yintercept = 0) +
        facet_wrap(vars(cluster_label)) +
        ylab("z-expression")+ 
        theme_bw() + 
        theme(text = element_text(face='bold', size = 15), 
              axis.title.x = element_blank(), 
              axis.text.x = element_text(angle = 90))
}


saveGenePlot <- function(symbols=NULL, gene_ids = NULL, path='./Figure/', width =18 , height = 13, device = 'png', dpi = 300){
    
    if(!is.null(gene_ids)){
        stopifnot(length(symbols) == length(gene_ids))
    }
    
    l <- list()
    for(i in 1:length(symbols)){
        print(symbols[i])
        if(!is.null(gene_ids)){
            l[[i]] <- genePlot(gene_id = gene_ids[i], dds = dds, stats = FALSE)
            
        }else{ l[[i]] <- genePlot(symbols[i], dds = dds, stats = FALSE)}
    }
    
    for(i in 1:length(l)){
        print(i)
        ggsave(filename = glue('{symbols[i]}.png'), 
               plot = l[[i]],
               device = device, units = "cm", width = width, height = height, path = path, dpi = dpi)
        
    }
    
    
}



