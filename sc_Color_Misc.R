# 2022-11-04 scRNAseq Utils by HSY
# BarPlotTheme - change barplot theme to highlight letters
# GetLegendOnly - get legend of plot using cowplot
# FeatureColGreyRed - new scale of featureplot color in grey to red 

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
