# 2022-11-04 scRNAseq Utils by HSY
# BarPlotTheme - change barplot theme to highlight letters
# GetLegendOnly - get legend of plot using cowplot
# FeatureColGreyRed - new scale of featureplot color in grey to red 
# BarPlot2 - modified BarPlot from ann081993(Using FetchData() from Seurat)
# addGeneFoldChange - for gene ontology analysis using clusterprofiler, goBarPlot
# goBarPlot
# goBarPlot2 - more concise representation of goBarPlot ... now modifying

library(ggplot2)
library(ggpubr)
library(ggsci)
library(viridis)
library(cowplot)
library(patchwork)
library(dplyr)
library(reshape2)
library(glue)

BarPlotTheme <- function(){
    return(theme(plot.title = element_text(size = 20, face = 'bold'),
                 axis.text = element_text(size = 15,face='bold'), 
                 axis.title=element_text(size = 15,face='bold'), 
                 legend.title = element_text(size=15, face='bold'),
                 legend.text = element_text(size=13)) & NoAxesTitle(no_text_x = FALSE))
}

GetLegendOnly <- function(x, legend = 'horizon'){
    legend <- cowplot::get_legend(x + theme(legend.position = "top")) 
    grid.newpage()
    grid.draw(legend)
    
    
}

FeatureColGreyRed <- function() {
    return(scale_color_gradientn(colors = c("gray", "gray", "red","darkred")))
}

BarPlot2 <- function(object, features = g, ncol = NULL, cols = NULL, error = "mean_se",
                    group.by = NULL, split.by = NULL, slot = "data", size = NULL, nolegend = FALSE) {
    found <- features %in% c(rownames(object), colnames(object[[]]))
    if(any(!found)){cat(paste0("The following requested features were not found: ",
                                   paste0(features[!found], collapse = ", "), "\n"))}

    features <- features[found] # feature is in gene or metadata 
    if(is.null(group.by)) { group.by = 'ident'}
    g_ex <- FetchData(object = object, slot = slot, vars = c(features, group.by,split.by))
    
    od <- order(object@reductions$pca@cell.embeddings[, "PC_1"])
    
    ncell <- nrow(g_ex)
    nfeat <- ncol(g_ex)-1
    if(!is.null(split.by)){nfeat <- nfeat - 1}
    if(is.null(ncol)) { ncol <- ceiling(sqrt(nfeat))}
    
    
    df <- g_ex %>% tibble::rownames_to_column('cell') %>% pivot_longer(., cols = -c(cell, {{group.by}}, {{split.by}}), names_to = 'gene', values_to = 'value') %>% dplyr::rename(split = {{split.by}}, ident = {{group.by}})
        
    df <- df[with(df, order(gene, ident)), ]
    df <- df %>% group_by(ident, gene) %>% mutate(med = quantile(value)[4])
    df <- df %>% group_by(ident, gene) %>% mutate(med = median(value, na.rm = TRUE))

    
    plist <- list()
    n = 1
    for(g in features) {
        df_subset <- df[df$gene == g, ]
        fill <- ifelse(is.null(split.by), "ident", "split")
        h2 <- ggbarplot(df_subset, x = "ident", y = "value", fill = fill, size = size,
                        palette = cols, color = "black", position = position_dodge(0.9), 
                        add = error, error.plot = "upper_errorbar", add.params = list(width = 0.3, size = size)) +
            ggtitle(g) + FeatureTitle() +
            theme(panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1), # element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank()) +
            scale_y_continuous(expand = expansion(mult = c(0,0.1))) + rotate_x_text(45)
        if(nolegend){h2 <- h2 & NoLegend()}
        #p <- wrap_plots(list(h1, h2), widths = c(3,1))
        plist[[n]] <- h2; n = n + 1
    }
    return(wrap_plots(plist, ncol = ncol))
}


addGeneFoldChange <- function(x, fc_data = top_genes_filt, top_genes = FALSE){
    
    # add gene fold change data to GO Enrichment analysis result. 
    # x is enrichGo@result object output of clusterprofiler::enrichGO
    # fc_data : dataframe contains fold change data of genes. columns :  gene, avg_log2FC
    # top_genes : make a new column that only contains N number of top fold change genes in each GO term. 
    
    x$geneID_FC <- 1
    for(i in 1:length(x$geneID)){
        gs <- stringr::str_split(string = x$geneID[i], pattern = '/')[[1]]
        fcs <- fc_data[fc_data$gene %in% gs, c("avg_log2FC",'gene'), drop = FALSE] %>% as.data.frame() %>% 
            dplyr::arrange(desc(abs(avg_log2FC))) %>% distinct(gene, .keep_all = TRUE)
        g_fc <- paste0(fcs$gene, '(',round(fcs$avg_log2FC,2),')')
        
        if(top_genes){
            x$geneID_FC_top_genes[i] <- paste0(g_fc[1:min(length(g_fc),top_genes)], collapse = '/')
        }
        
        x$geneID_FC[i] <- paste0(g_fc, collapse = '/')
        
    }
    return(x)
}


goBarPlot <- function(x, plot.title='',pvalue = 'p.adjust' ,showCategory = 12, bar.col='grey80', gene.text.size=3,gene.width=70,
                      layout.width=c(1,0.5,3), include.genelist = TRUE, go.text.size=13, go.text.width=30
                      ){
    
    category <- min(nrow(x), showCategory)
    df <- x %>% arrange(p.adjust) %>% 
        slice_min(p.adjust, n = category, with_ties = FALSE)  %>% 
        mutate(Description = factor(Description, level = rev(Description)), 
               log10padj = log10(pvalue),
               position = rep(1,category))
    # p1 : Count Log Plot 
    p1 <- df %>% 
        ggplot(aes(x = Description ,y= Count)) +
        geom_col(fill = bar.col, color = 'black', size = 0.7, width = 0.8) +
        scale_x_discrete(labels =function(x){ stringr::str_wrap(x, width = go.text.width)}) +
        coord_flip() +
        theme_void() + 
        theme(
            axis.title.x = element_text(size = 12, face = "bold") ,
            axis.title.y = element_blank() ,
            axis.text = element_text(face = "bold", size = 12),
            axis.text.x = element_text(face = "bold", size = 11),
            axis.text.y = element_text(hjust = 0.95, size=go.text.size), # go.text.size
            legend.position = 'none',
            plot.background = element_rect(fill = 'transparent', color = NA),
            panel.background = element_rect(fill = 'transparent',color = 'black'),
            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm") ,
            legend.background = element_rect(fill='transparent',color = NA),
            legend.box.background = element_rect(fill='transparent', color = NA)
        ) 
    
    # p2 : Geom Text Description of log P value 
    p2 <- df %>% 
        ggplot(aes(x = Description )) +
        geom_text(aes(y = position, label = round(log10padj,2)), hjust = 0.5,vjust = -0, size = 4.5, fontface = 'bold') + coord_flip() + theme_void()+ geom_vline(xintercept = seq(0.5,(category + 1.5), 1)) +
        labs(y = '', x= NULL) + ggtitle('log10(P-Value)') + theme(text = element_text(face = 'bold'), 
                                                                  plot.title = element_text(hjust = 0.5))
    # p3 : Gene list
    p3 <- df %>% 
        ggplot(aes(x = Description )) +
        geom_text(aes(y = position, label = stringr::str_wrap(geneID_FC_top_genes, width = gene.width)), size=gene.text.size) + 
        geom_vline(xintercept = seq(0.5,(category + 1.5), 1)) +
        coord_flip() + theme_void()+
        labs(y = '', x= NULL) + ggtitle('Genes') + theme(text = element_text(face = 'bold'), 
                                                         plot.title = element_text(hjust = 0.5))
    if(include.genelist){
    p <- p1 + p2 + p3 + plot_layout(widths = layout.width) + 
        patchwork::plot_annotation(title = plot.title, theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.4),
                                                                     plot.background = element_rect(fill = 'transparent', color = NA),
                                                                     panel.background = element_rect(fill = 'transparent',color = NA),
                                                                     legend.background = element_rect(fill='transparent', color = NA),
                                                                     legend.box.background = element_rect(fill='transparent', color = NA)
        ))}
    else{p <- p1 + p2 + plot_layout(widths = layout.width) + 
        patchwork::plot_annotation(title = plot.title, theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.4),
                                                                     plot.background = element_rect(fill = 'transparent', color = NA),
                                                                     panel.background = element_rect(fill = 'transparent',color = NA),
                                                                     legend.background = element_rect(fill='transparent', color = NA),
                                                                     legend.box.background = element_rect(fill='transparent', color = NA)
        ))}
    return(p)
    
}

goBarPlot2 <- function(x, plot.title='',pvalue = 'p.adjust' ,showCategory = 12, bar.col='grey80', gene.text.size=3,gene.width=70,
                      layout.width=c(1,0.1), include.genelist = TRUE, go.text.size=13, go.text.width=30
){
    
    category <- min(nrow(x), showCategory)
    df <- x %>% arrange(p.adjust) %>% 
        slice_min(p.adjust, n = category, with_ties = FALSE)  %>% 
        mutate(Description = factor(Description, level = rev(Description)), 
               log10padj = log10(pvalue),
               position = rep(1,category))
    # p1 : Count Log Plot 
    p1 <- df %>% 
        ggplot(aes(x = Description ,y= Count)) +
        geom_col(fill = bar.col, color = 'black', size = 0.7, width = 0.8) +
        scale_x_discrete(labels =function(x){ stringr::str_wrap(x, width = go.text.width)}) +
        geom_text(aes(label = stringr::str_wrap(geneID_FC_top_genes, width = gene.width, exdent = 0)), 
                  size=gene.text.size, position=position_stack(vjust = 0), hjust = -0.01) + 
        geom_text(aes(label = round(log10padj,2)), position=position_stack(vjust = 1), hjust = 0.01) + 
        
        coord_flip() + 
        scale_y_continuous(expand = c(0,1))+
        theme_void() + 
        theme(
            axis.title.x = element_text(size = 12, face = "bold") ,
            axis.title.y = element_blank() ,
            axis.text = element_text(face = "bold", size = 12),
            axis.text.x = element_text(face = "bold", size = 11),
            axis.text.y = element_text(hjust = 0.95, size=go.text.size), # go.text.size
            legend.position = 'none',
            plot.background = element_rect(fill = 'transparent', color = NA),
            panel.background = element_rect(fill = 'transparent',color = 'black'),
            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm") ,
            legend.background = element_rect(fill='transparent',color = NA),
            legend.box.background = element_rect(fill='transparent', color = NA)
        ) 
    
    # p2 : Geom Text Description of log P value 
    p2 <- df %>% 
        ggplot(aes(x = Description )) +
        geom_text(aes(y = position, label = round(log10padj,2)), hjust = 0.5,vjust = -0, size = 4.5, fontface = 'bold') + coord_flip() + theme_void()+ geom_vline(xintercept = seq(0.5,(category + 1.5), 1)) +
        labs(y = '', x= NULL) + ggtitle('log10(P-Value)') + theme(text = element_text(face = 'bold'), 
                                                                  plot.title = element_text(hjust = 0.5))

    if(include.genelist){
        p <- p1 + p2 + plot_layout(widths = layout.width) + 
            patchwork::plot_annotation(title = plot.title, theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.4),
                                                                         plot.background = element_rect(fill = 'transparent', color = NA),
                                                                         panel.background = element_rect(fill = 'transparent',color = NA),
                                                                         legend.background = element_rect(fill='transparent', color = NA),
                                                                         legend.box.background = element_rect(fill='transparent', color = NA)
            ))}
    else{p <- p1 + p2 + plot_layout(widths = layout.width) + 
        patchwork::plot_annotation(title = plot.title, theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.4),
                                                                     plot.background = element_rect(fill = 'transparent', color = NA),
                                                                     panel.background = element_rect(fill = 'transparent',color = NA),
                                                                     legend.background = element_rect(fill='transparent', color = NA),
                                                                     legend.box.background = element_rect(fill='transparent', color = NA)
        ))}
    return(p)
    
}



