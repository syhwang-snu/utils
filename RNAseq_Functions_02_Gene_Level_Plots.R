library(patchwork)
library(cowplot)
library(ggpubr)
library(ggh4x)

genePlot <- function(symbol, stats = TRUE, 
                     group = c("CHOW_WT",   "CDAHFD_WT" ,"CHOW_KO"   ,"CDAHFD_KO"), 
                     title = NULL, 
                     my_comparisons = list( c("CHOW_WT", "CHOW_KO"), 
                                            c("CDAHFD_WT", "CDAHFD_KO"),
                                            c("CHOW_KO","CDAHFD_KO"), 
                                            c("CHOW_WT", "CDAHFD_WT")), force_zero = TRUE)
    {

    g2id <- id2symbol.g29$gene_id[id2symbol.g29$SYMBOL == symbol]
    
    geneCounts <- plotCounts(dds, gene = g2id, intgroup = "group", normalized = TRUE,
                             returnData = TRUE)
    if(is.null(title)){
        title <- paste0(symbol," Expression")
    }
    
    if(!is.null(group)){
        geneCounts$group <- factor(geneCounts$group, levels =  group)
    }
    
    p <- ggboxplot(geneCounts, x = "group", y = "count", 
                   add = "jitter", 
                   title = title, 
                   ylab = "Normalized Counts") + 
        theme(axis.title = element_text(size=15), 
              axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9, face='bold')) 
    
    
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


clusterPlot <- function(gene_ids, h= 4, cluster_k= 4, group= c('CHOW_WT','CDAHFD_WT','CHOW_KO','CDAHFD_KO'), interaction = NULL, ncol=NULL){
    
    
    normal.count.mtx <- normal.count.mtx.all[gene_ids, ]
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
    
    if(!is.null(interaction)){
        
        trans_cts_cluster$isKO <- ifelse(trans_cts_cluster$group %in% c('CHOW_WT','CDAHFD_WT'),'WT','KO')
        trans_cts_cluster$isKO <- factor(trans_cts_cluster$isKO, levels = c('WT','KO'))
        trans_cts_cluster$org.group <- trans_cts_cluster$group 
        levels(trans_cts_cluster$group) <- c('CHOW','HFD','CHOW','HFD')
        
        p<- trans_cts_cluster %>% 
            ggplot(aes(group, cts)) +
            geom_line(aes(group = gene_id), alpha = 0.3) +
            geom_line(stat = "summary", fun = "median", colour = "brown", linewidth = 1.5, 
                      aes(group = 1)) + 
            geom_hline(yintercept = 0) +
            facet_nested_wrap(vars(cluster_label,isKO), dir = "h", ncol = ncol) +
            ylab("z-expression") + scale_x_discrete(expand = expansion(mult = 0.2)) + 
            theme_bw() + 
            theme(text = element_text(face='bold', size = 15), 
                  axis.title.x = element_blank())
        
    }
    
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
              axis.title.x = element_blank())
    
    
    print(p)
    return(trans_cts_cluster)
    
}

