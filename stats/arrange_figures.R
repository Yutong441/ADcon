# from https://github.com/Yutong441/TBdev/blob/main/R/arrange_figures.R
get_highlight_font <- function (fontface, fontsize, font_fam){
        return (grid::gpar (fontface=fontface, fontsize=fontsize, fontfamily=font_fam))
}

add_grob_text <- function (grob_plot, content, ...){
        highlight_font <- get_highlight_font (...)
        return (gridExtra::arrangeGrob (grob_plot, 
                             left=grid::textGrob (content, 
                                            x = unit(1, "npc"),
                                            y = unit(0.95, "npc"),
                                            gp = highlight_font ) ))
}

#" Replicate a vector x to a target length target_n
expand_length <- function (x, target_n){
        discrep <- length (x) - target_n
        if (discrep < 0 ){
                print ("not enough values, replicating the last item to fill")
                new_x <- c(x, rep (x[length(x)], -discrep) )
        }else if (discrep > 0){ 
                print (paste("too many values, only the first", target_n, "will be used"))
                new_x <- x [1:target_n]
        }else{new_x <- x
        }
        return (new_x)
}

push_plot_viewport <- function (plot_ob, text_args=NULL){
        if ( "ggplot" %in% class (plot_ob) ){
                #grid.draw (grid.grabExpr ( print (grob_list[[i]])))
                # prevent text from being cutoff from the borders
                gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot_ob))
                gt$layout$clip[gt$layout$name == "panel"] <- "off"
                grid::grid.draw(gt)
        }else if ("gtable" %in% class (plot_ob)){
                grid::grid.draw(plot_ob)
        }else if ("data.frame" %in% class (plot_ob) ) {
                colnames (plot_ob) <- gsub ("\\.", " ", colnames (plot_ob))
                gridExtra::grid.table (plot_ob, rows=NULL)
        }else if ("character" %in% class (plot_ob) & !is.null (text_args)){
                custom_text_str (plot_ob, text_args)
        }else if ("array" %in% class (plot_ob)){
                grid::grid.raster (plot_ob)
        }
}

#" Same function as `grid.arrange`, with the functionality of integrating
#" ComplexHeatmap and adding panel labels
#" 
#" @param grob_list a list of plots. They can be ggplots, ComplexHeatmap,
#" data.frame, character, image array
#" @param save_path where the plots are saved
#" @param grid_layout, same as the argument in `grid.arrange`
#" @param margin_width the width of the margin in inches. The default is 0.79in
#" (2cm). Same unit in `margin_height`, `page_width`, `page_height`
#" @param plot_width the width of each subplot. 
#" @param plot_height same concept as plot_width
#" @importFrom grid pushViewport viewport grid.layout grid.rect grid.draw
#" popViewport grid.text grid.grabExpr 
#" @export
arrange_plots <- function (grob_list, save_path, grid_layout, panel_label=NULL,
                           margin_width=0.79, margin_height=0.79,
                           #panel_spacing=0.1975, 
                           panel_spacing=0.01, 
                           plot_width=9, plot_height=9, 
                           aes_param=list (highlight_font = list (fontface="bold", 
                                        fontsize=40), font_fam = "Arial") 
                           ){
        
        highlight_font <- get_highlight_font (aes_param$highlight_font$fontface, 
                                              aes_param$highlight_font$fontsize, 
                                              aes_param$font_fam)
        plot_width <- expand_length (plot_width, ncol (grid_layout)) # a vector
        plot_height <- expand_length (plot_height, nrow (grid_layout))
        # calculate the width and height for all the plots
        page_width <-  sum (plot_width) # a single scalar
        page_height <- sum (plot_height)
        # calculate the margin for the entire page
        width_pro <- 1 - margin_width/page_width
        height_pro <- 1 - margin_height/page_height

        # default panel labels: A, B, C, D ... 
        if (is.null (panel_label)){ panel_label <- LETTERS [1:length (grob_list)]}
        # cairo_pdf seems to offer better support with UTF-8 encoding
        grDevices::cairo_pdf (save_path, width=page_width, height=page_height)
        grid::grid.newpage ()
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow (grid_layout), ncol(grid_layout)) ,
                              width=width_pro, height=height_pro))

        # start plotting
        for (i in 1:length (grob_list)){
                print (paste ("arranging figure", i))
                row_pos <- unique (which (grid_layout == i, arr.ind=T) [, "row"])
                col_pos <- unique (which (grid_layout == i, arr.ind=T) [, "col"])

                # calculate the plotting dimensions
                # calculate the margin for each subplot
                subplot_width <- sum (plot_width [col_pos])
                subplot_height <-sum (plot_height[row_pos])
                width_pro_panel <- 1 - 2*panel_spacing/subplot_width
                height_pro_panel <- 1 - 2*panel_spacing/subplot_height

                grid::pushViewport (grid::viewport (layout.pos.col=col_pos,layout.pos.row=row_pos))
                # in order to increase spacing between panels, it is necessary
                # to draw an empty rectangle
                grid::grid.rect(gp=grid::gpar(lty=0)) #lty=0 means blank line
                grid::pushViewport (grid::viewport (layout.pos.col=col_pos,layout.pos.row=row_pos, 
                                        width=width_pro_panel, height=height_pro_panel))

                # need different functions for ggplot, table and ComplexHeatmap
                # objects
                push_plot_viewport (grob_list[[i]])
                grid::grid.text (panel_label[i], x=unit (0.01, "npc"), y=unit (0.99, "npc"), 
                           gp=highlight_font)
                print ("finished")
                grid::popViewport (2)
        }
        grDevices::dev.off ()
}

#" Save a list of plots individually
#"
#" @description This function and its arguments are based upon `arrange_plots`.
#" Instead of saving all the plots in a single pdf, this function saves the
#" plots individually in a directory. The usage is exactly the same as
#" `arrange_plots`.
#" @export
save_indiv_plots <- function (grob_list, save_dir, grid_layout, plot_width=9,
                              plot_height=9, panel_label=NULL,
                              aes_param=list (highlight_font = list (fontface="bold", 
                                        fontsize=40), font_fam = "Arial") 
                     ) {
    plot_width <- expand_length(plot_width, ncol(grid_layout)) # a vector
    plot_height <- expand_length(plot_height, nrow(grid_layout))
    highlight_font <- get_highlight_font(aes_param$highlight_font$fontface,
                                         aes_param$highlight_font$fontsize,
                                         aes_param$font_fam)
    if (!dir.exists (save_dir) ){dir.create (save_dir)}
    if (is.null (panel_label)){ panel_label <- LETTERS [1:length (grob_list)]}
    for (i in 1:length (grob_list)){
        print(paste("arranging figure", i))
        row_pos <- unique(which(grid_layout == i, arr.ind=T) [, "row"])
        col_pos <- unique(which(grid_layout == i, arr.ind=T) [, "col"])
        page_width <- sum(plot_width [col_pos])
        page_height <- sum(plot_height[row_pos])

        save_path <- paste(save_dir, "/panel", panel_label[i], ".pdf", sep="")
        grDevices::cairo_pdf(save_path, width=page_width, height=page_height)
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 1)))
        # move the plot labels slightly inwards horizontally
        # (unit(0.05, "npc")) because individual plots do not have as
        # much padding around them
        push_plot_viewport(grob_list[[i]])
        #grid::grid.text (panel_label[i], x=grid::unit (0.05, "npc"), 
        #                 y=grid::unit (0.95, "npc"), gp=highlight_font)
        grDevices::dev.off()
        }
}
