#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# 2022-06-29 ASV Plot Utils HSY
# utils for ASV processing & visualization. 
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(ggtranscript)
library(rtracklayer)





gtf2exon <- function(ens_gtf){
    # make gtf to contain only exons. 
    
    ens_gtf_exons <- as.data.frame(ens_gtf) %>% dplyr::filter(type == "exon")
    ens_gtf_exons <- ens_gtf_exons %>% 
        dplyr::select(
            seqnames,
            start,
            end,
            strand,
            type,
            gene_name,
            transcript_name,
            transcript_id,
            transcript_type
        )
    return(ens_gtf_exons)
    
}


geneTranscriptPlot <- function(g = 'Atp6v0a1' , x = topGenes_table, gtf_exon_table = ens_gtf_exons, CPM_columns = c('Ctrl','PA'), rescaled = TRUE){
    # x : gene & log Fold Change Information table : adapted from edgeR toptag 
    # gtf_exon_table : gtf imported by rtracklayer that contains only exon. as.data.frame() %>% filter() ...
    # 
    
    g_df <- gtf_exon_table %>% dplyr::filter(gene_name == g)
    g_df <- g_df %>% left_join(x[,c('tx_name','logFC',CPM_columns)], by=c('transcript_id'='tx_name'))
    g_df_rescaled <- shorten_gaps(
        g_df, 
        to_intron(g_df, "transcript_name"), 
        group_var = "transcript_name"
    )
    
    p1 <- g_df_rescaled %>%
        dplyr::filter(type == "exon") %>%
        ggplot(aes(
            xstart = start,
            xend = end,
            y = transcript_name
        )) +
        geom_range(
            aes(fill = logFC) 
        ) + scale_fill_gradientn(colors=c("Blue","white","Red"), limits = c(-3,3),
                                 na.value = "white", 
                                 oob = scales::squish)+
        geom_intron(
            data = g_df_rescaled %>% dplyr::filter(type == "intron"), 
            arrow.min.intron.length = 200
        ) + theme_bw() + theme(axis.text.y = element_blank(), 
                               axis.title.y = element_blank(),
                               plot.background = element_rect(fill = 'transparent', color = NA),
                               panel.background = element_rect(fill = 'transparent',color = NA),
                               legend.background = element_rect(fill='transparent',color = NA),
                               legend.box.background = element_rect(fill='transparent',color = NA))
    if(rescaled == FALSE){
        # no rescale for introns. 
        g_df_rescaled <- g_df
        p1 <- g_df_rescaled %>%
            dplyr::filter(type == "exon") %>%
            ggplot(aes(
                xstart = start,
                xend = end,
                y = transcript_id
            )) +
            geom_range(
                aes(fill = logFC) 
            ) + scale_fill_gradientn(colors=c("Blue","white","Red"), limits = c(-3,3),
                                     na.value = "white", 
                                     oob = scales::squish)+
            geom_intron(
                data = to_intron(g_df_rescaled, "transcript_id"),
                aes(strand = strand)
            ) +
            theme_bw() + theme(axis.text.y = element_blank(), 
                               axis.title.y = element_blank())
        
    }
    
    g_tx_names <- unique(g_df[,c('transcript_id','transcript_name')]) %>% dplyr::rename(tx_name = transcript_id)
    g_tx_names$transcript_name <- factor(g_tx_names$transcript_name)
    g_cpms <- cpm_y_all %>% dplyr::select(Symbol, tx_name, transcript_name, Ctrl, PA) %>% 
        dplyr::filter(Symbol == g) %>% 
        group_by(Symbol) %>% 
        pivot_longer(cols = c('Ctrl','PA'), names_to = "Sample",values_to = "CPM")  %>% 
        mutate(transcript_name = factor(transcript_name, levels = rev(levels(g_tx_names$transcript_name))))
    
    
    p2 <- g_cpms  %>% 
        ggplot(aes(x = Sample ,y= CPM),fill='Black') +
        geom_col() + 
        coord_flip() +
        theme_light()  + scale_y_reverse(lim= c(max(10, max(g_cpms$CPM)),0)) +
        facet_grid(rows = vars(transcript_name), 
                   scales = "free_y", # Let the x axis vary across facets.
                   space = "free_y",  # Let the width of facets vary and force all bars to have the same width.
                   switch = "y", ) +
        theme(strip.placement = "outside",                      # Place facet labels outside x axis labels.
              strip.background = element_blank(), # Make facet label background white.
              axis.title.x = element_text(size = 12, face = "bold") ,
              axis.title.y = element_blank() ,
              axis.text = element_text(face = "bold", size = 9),
              axis.text.x = element_text(face = "bold", size = 11),
              strip.text = element_text(color = 'black', face = "bold", size = 12),
              strip.text.y.left = element_text(angle = 0), 
              legend.position = 'none',
              plot.background = element_rect(fill = 'transparent', color = NA),
              panel.background = element_rect(fill = 'transparent',color = NA),
              legend.background = element_rect(fill='transparent',color = NA),
              legend.box.background = element_rect(fill='transparent', color = NA)
        ) 
    p <- p2 + p1 + plot_layout( widths = c(1,2)) + 
        patchwork::plot_annotation(title = g, theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.4),
                                                            plot.background = element_rect(fill = 'transparent', color = NA),
                                                            panel.background = element_rect(fill = 'transparent',color = NA),
                                                            legend.background = element_rect(fill='transparent', color = NA),
                                                            legend.box.background = element_rect(fill='transparent', color = NA)
        ))
    
    return(p)
    
}



addGeneFoldChange <- function(x, fc_data = top_genes_filt, top_genes = NULL){
    # add gene fold change data to GO Enrichment analysis result. 
    # x is enrichGo@result object output of clusterprofiler::enrichGO
    # fc_data : dataframe contains fold change data of genes. columns :  Symbol, logFC
    # top_genes : make a new column that only contains N number of top fold change genes in each GO term. 
    
    x$geneID_FC <- 1
    for(i in 1:length(x$geneID)){
        gs <- stringr::str_split(string = x$geneID[i], pattern = '/')[[1]]
        fcs <- fc_data[fc_data$Symbol %in% gs, c("logFC",'Symbol'), drop = FALSE] %>% as.data.frame() %>% 
            dplyr::arrange(desc(abs(logFC))) %>% distinct(Symbol, .keep_all = TRUE)
        g_fc <- paste0(fcs$Symbol, '(',round(fcs$logFC,2),')')
        if(!is.null(top_genes)){
            x$geneID_FC_top_genes[i] <- paste0(g_fc[1:min(length(g_fc),top_genes)], collapse = '/')
        }
        x$geneID_FC[i] <- paste0(g_fc, collapse = '/')
        
    }
    return(x)
}


print(cat(" !!!
      gtf2exon
      geneTranscriptPlot
      addGeneFoldChange
      
      Loaded
      !!!
      "))



