#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# 2022-06-29 HSY 
# ggplot graph to ppt 
#
#



library(ggpubr)


graph2ppt <- function(l, path='figs.pptx', left = 0, top = 1, width = 14, height = 4.95, temp_path = 'template.pptx'){
    # width / height (inch) * 0.394 -> (cm) 
    
    width <- width * 0.394
    height <- height * 0.394
    
    dmls <- list()
    for(i in 1:length(l)){
        try({
            img <- l[[i]]
            dml <- rvg::dml(ggobj = img, bg = 'transparent')
            dmls[[i]] <- dml
        })
    }
    
    if (!file.exists(path)) {
        out <- officer::read_pptx(path = temp_path)
    }
    # if file exist, append slides to existing file ----
    else {
        out <- officer::read_pptx(path)
    }
    
    free_loc <- officer::ph_location(
        left = left, top = top, 
        width = width, height = height)
    
    for(i in 1:length(dmls)){
        
        try({
            dml <- dmls[[i]]
            out <- out %>% 
                officer::add_slide() %>% 
                officer::ph_with(value = dml, location = free_loc) })
    }
    print(out,target = path)
}


