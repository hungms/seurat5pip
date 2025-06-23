#' Plot cell cycle
#'
#' @param obj Seurat object
#' @param split.by Metadata column to split the object by
#' @param output_dir Directory to save the output
#' @return List of plots
#' @export
plot_cc <- function(obj, var = "Phase", split.by = NULL, output_dir = NULL){

    # validate split.by
    split.by <- validate_split.by(split.by, obj)

    # plot
    metadata <- obj@meta.data

    if(is.factor(metadata[[var]])){
        levels <- levels(metadata[[var]])}
    else{
        levels <- unique(metadata[[var]])}

    # get palette
    pal <- get_palette("viridis", n = length(levels))
    names(pal) <- levels
    
    plot <- plot_percent(metadata, group.by = split.by, var = var, pal = pal, legend.ncol = 1, filename = "plot_cc_seurat.png", output_dir = output_dir)

    # return plot
    return(plot)
}
