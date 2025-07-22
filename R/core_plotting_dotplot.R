#' plot_dotplot
#' 
#' @param object Seurat object
#' @param features Features to plot
#' @param assay Assay to use
#' @param cols Colors to use
#' @param col.min Minimum value for color scale
#' @param col.max Maximum value for color scale
#' @param dot.min Minimum value for dot size
#' @param dot.scale Dot scale
#' @param idents Identities to plot
#' @param group.by Group by
#' @param split.by Split by
#' @param cluster.idents Cluster identities
#' @param scale Scale
#' @param scale.by Scale by
#' @param scale.min Minimum value for scale
#' @param scale.max Maximum value for scale
#' @param plot.title Plot title
#' 
#' @return A ggplot object
#' @export
#' @importFrom Seurat FetchData

plot_dotplot <- function(
  object,
  features,
  assay = NULL,
  cols = get_palette("Reds", n = 7)[c(1,7)],
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  base_size = 16,
  dot.scale = 14,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA,
  plot.title = " ",
  diffexp = NULL
) {

  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  
  cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]), idents = idents))
  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  
  # Add split information to data.features
  if (!is.null(x = split.by)) {
    splits <- FetchData(object = object, vars = split.by)[cells, split.by]
    data.features$split <- splits
  }
  
  # Create data.plot with split information
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      feature_cols <- colnames(data.features)[colnames(data.features) %in% features]
      # For each split value in this identity group
      if (!is.null(x = split.by)) {
        split_vals <- unique(data.features$split[data.features$id == ident])
        # Create separate entries for each split value
        result <- lapply(split_vals, function(split_val) {
          data.use <- data.features[data.features$id == ident & data.features$split == split_val, feature_cols, drop = FALSE]
          if (nrow(data.use) == 0) return(NULL)  # Skip if no data for this combination
          avg.exp <- apply(
            X = data.use,
            MARGIN = 2,
            FUN = function(x) {
              return(mean(x = expm1(x = x)))
            }
          )
          pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
          list(avg.exp = avg.exp, pct.exp = pct.exp, split = split_val, id = ident)
        })
        result <- result[!sapply(result, is.null)]  # Remove NULL entries
        return(result)
      } else {
        # Original behavior for no split
        data.use <- data.features[data.features$id == ident, feature_cols, drop = FALSE]
        avg.exp <- apply(
          X = data.use,
          MARGIN = 2,
          FUN = function(x) {
            return(mean(x = expm1(x = x)))
          }
        )
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
        return(list(list(avg.exp = avg.exp, pct.exp = pct.exp, split = NULL, id = ident)))
      }
    }
  )
  
  # Flatten the nested list
  data.plot <- unlist(data.plot, recursive = FALSE)
  names(x = data.plot) <- make.unique(sapply(data.plot, function(x) x$id))
  
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = function(x) unlist(x[c("avg.exp", "pct.exp")]))
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  
  # Convert data.plot to data frame with split information
  data.plot <- lapply(
    X = data.plot,
    FUN = function(x) {
      data.use <- as.data.frame(x[c("avg.exp", "pct.exp")])
      data.use$features.plot <- names(x$avg.exp)
      data.use$id <- x$id
      if (!is.null(x = split.by)) {
        data.use$split <- x$split
      }
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = log1p(data.use))
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  data.plot$id <- factor(x = data.plot$id, levels = rev(id.levels))

  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp.scaled')) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme_prism(border = T, base_size = 16) +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.5, linetype = "dashed"),
      axis.title.y = element_blank(), 
      legend.title = element_text(), 
      strip.background = element_rect(),
      axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(size = guide_legend(title = 'Percent\nExpressed')) +
    labs(
      title = plot.title,
      x = '',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    scale_color_gradient(low = cols[1], high = cols[2], breaks = c(-1, 0, 1), labels = c("-1", "0", "1")) +
    guides(color = guide_colorbar(
        title = 'Average\nExpression', 
        order = 1, 
        title.position = "top",
        direction = "horizontal",
        frame.colour = "black",
        ticks.colour = "black",
        barwidth = 7,
        barheight = 1.5))

  ## add diffexp
  if(!is.null(diffexp)){
    columns_required <- c("features.plot", "id", "p_val_adj")
    if(!is.null(split.by)){
      columns_required <- c(columns_required, "split")}
    stopifnot(all(columns_required %in% colnames(diffexp)))

    diffexp$p_signif <- ifelse(diffexp$p_val_adj < 0.05, "padj < 0.05", "ns")
    diffexp$p_signif <- factor(diffexp$p_signif, levels = c("ns", "padj < 0.05"))

    diffexp <- diffexp %>%
      filter(features.plot %in% features) %>%
      merge(., data.plot, by = c("features.plot", "id"), all.x = T)
    
    plot <- plot + 
      geom_point(data = diffexp, mapping = aes_string(size = 'pct.exp', color = 'avg.exp.scaled', stroke = 'p_signif'), shape = 21) #+
      #scale_stroke(values = c("ns" = 0, "padj < 0.05" = 1))
    }
  
  # Handle faceting
  if (!is.null(x = feature.groups) || !is.null(x = split.by)) {
    # Determine facet formula based on what's available
    if (!is.null(x = feature.groups) && !is.null(x = split.by)) {
      # Both feature groups and split by are present
      facet_formula <- as.formula("split ~ feature.groups")
    } else if (!is.null(x = feature.groups)) {
      # Only feature groups
      facet_formula <- as.formula("~ feature.groups")
    } else {
      # Only split by
      facet_formula <- as.formula("~ split")
    }
    
    plot <- plot + facet_grid(
      facets = facet_formula,
      scales = "free",
      space = "free",
      switch = NULL
    ) + theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_rect(fill = "#f7f7f7"),
      strip.text.x = element_text(size = 16, margin = margin(0.3,0,0.3,0, "cm")),
      strip.text.y = element_text(size = 16, margin = margin(0,0.3,0,0.3, "cm"))
    )
  }
  
  return(plot)
}