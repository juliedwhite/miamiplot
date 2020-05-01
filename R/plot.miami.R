#' Miami plot
#'
#' Create a Miami plot ggplot2 object.
#'
#' @param data A data.frame object. Required.
#' @param split.by A character vector. The name of the column to use for
#'   splitting into top and bottom sections of the Miami plot. Required.
#' @param split.at A character or numeric vector. If numeric, the top plot will
#'   contain your results where the values in the `split.by` column are
#'   >= `split.at`. If character, top plot will contain your results where
#'   the values in the `split.by` column are equal to `split.at`. Required.
#' @param chr The name of the column containing your chromosome information.
#'   Defaults to "chr"
#' @param pos The name of the column containing your position information.
#'   Defaults to "pos"
#' @param p The name of the column containing your p-value information.
#'   Defaults to "P"
#' @param genomewideline Should we draw a genome-wide significance line? Set to
#'   FALSE if you do not want this line. Defaults to 5e-8, or supply your own number.
#' @param suggestiveline Should we draw a suggestive significance line? Set to
#'   FALSE if you do not want this line. Defaults to 1e-5, or supply your own number.
#' @param hits.label Either the name of the column(s), max. 2, to use for
#'   automatically labeling n hits, determined using `top.n.hits`, or a
#'   character vector of probes/genes/SNPs to label. If supplying a list of
#'   genes, it is helpful to also supply `top.n.hits` to limit the labels
#'   to the top N results, though this isn't necessary.
#'   Defaults to FALSE: labels aren't displayed.
#' @param top.n.hits How many of the top hits do you want to label? Defaults to
#'   FALSE: the number of labels aren't filtered.
#' @export
#' @return a \code{ggplot2} object
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#' @keywords manhattan miami plot SNP genetics
#'

ggmiami <- function(
  data,
  split.by,
  split.at,
  chr="chr",
  pos="pos",
  p = "P",
  genomewideline = 5e-8,
  suggestiveline = 1e-5,
  hits.label = FALSE,
  top.n.hits = FALSE
  ) {

  # Check the required input
  check.miami.input(data = data, split.by = split.by, split.at = split.at)

  # Identify column index for chromosome information.
  chr.indx <- find.col.indx(data = data, col = chr)
  chr.name <- chr

  # Identify column index for position information
  pos.indx <- find.col.indx(data = data, col = pos)
  pos.name <- pos

  # Identify column index for p-value information
  p.indx <- find.col.indx(data = data, col = p)
  p.name <- p

  # Test if chromosome, position, pvalue columns are numeric or integer.
  check.miami.numeric(data = data, col.indx = chr.indx)
  check.miami.numeric(data = data, col.indx = pos.indx)
  check.miami.numeric(data = data, col.indx = p.indx)

  # To make the plot, we need to know the cumulative position of each probe/snp
    # across the whole genome.
  data <- data %>%
    #Group by chromosome
    group_by(!!sym(chr.name)) %>%
    # Compute chromosome size
    summarise(chrlength = max(!!sym(pos.name))) %>%
    # Calculate cumulative position of each chromosome
    mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>%
    # Remove chr length column
    select(-chrlength) %>%
    # Temporarily add the cumulative length of each chromosome to the initial
      # dataset
    left_join(data, ., by=chr.name) %>%
    # Sort by chr then position
    arrange(!!sym(chr.name), !!sym(pos.name)) %>%
    # Add the position to the cumulative chromosome length to get the position
      # of this probe relative to all other probes
    mutate(rel.pos = !!sym(pos.name) + cumulativechrlength) %>%
    # Remove cumulative chr length column
    select(-cumulativechrlength)

  # However, we don't want to print the cumulative position on the x-axis, so we
    # need to provide the chromosome position labels relative to the entire genome.
  axis.data = data %>%
    #Group by chromosome
    group_by(!!sym(chr.name)) %>%
    #Find the center of the chromosome
    summarize(chr.center=(max(rel.pos) + min(rel.pos))/2)

  # To create a symmetric looking plot, calculate the maximum p-value
  maxp <- ceiling(max(-log10(data[,p.indx])))

  # Depending on what the user has input for split.by and split.at, make top
    # and bottom data.
  if(is.numeric(split.at)){
    # Create top data using numeric info.
    top.data <- data %>%
      filter(!!sym(split.by) >= split.at)
    # Create bottom data using numeric info
    bot.data <- data %>%
      filter(!!sym(split.by) < split.at)
  } else if(is.character(split.at)){
    # Create top data using character specified by user
    top.data <- data %>%
      filter(!!sym(split.by) == split.at)
    # Create bottom data using character not specified by user
    bot.data <- data %>%
      filter(!!sym(split.by) != split.at)
  }

  # Create base top plot
  top.plot <- ggplot(data = top.data, aes(x=rel.pos, y=-log10(!!sym(p.name)))) +
    geom_point(aes(color=as.factor(!!sym(chr.name))), size=0.25) +
    scale_color_manual(values = rep(c("black", "grey"), nrow(axis.data))) +
    scale_x_continuous(labels = pull(axis.data[,1]),
                       breaks = axis.data$chr.center,
                       expand = expansion(mult=0.01)) +
    scale_y_continuous(limits = c(0, maxp),
                       expand = expansion(mult = c(0.02,0))) +
    labs(x = "", y = expression('-log'[10]*'(P)')) +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank(),
          plot.margin = margin(b=0))

  # Create base bottom plot
  bottom.plot <- ggplot(data = bot.data, aes(x=rel.pos, y=-log10(!!sym(p.name)))) +
    geom_point(aes(color=as.factor(!!sym(chr.name))), size=0.25) +
    scale_color_manual(values = rep(c("black", "grey"), nrow(axis.data))) +
    scale_x_continuous(breaks = axis.data$chr.center, position = "top",
                       expand = expansion(mult=0.01)) +
    scale_y_reverse(limits = c(maxp, 0), expand = expansion(mult = c(0,0.02))) +
    labs(x = "", y = expression('-log'[10]*'(P)')) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_blank(),
          axis.title.x = element_blank(), plot.margin = margin(t=0))

  # If the user has requested a suggetive line, add:
  if(!suggestiveline){
    top.plot <- top.plot +
      geom_hline(yintercept = -log10(suggestiveline), color = "blue",
                 linetype = "dashed", size = 0.3)
    bottom.plot <- bottom.plot +
      geom_hline(yintercept = -log10(suggestiveline),color = "blue",
                 linetype = "dashed", size = 0.3)
  }

  # If the user has requested a genome-wide line, add:
  if(!genomewideline){
    top.plot <- top.plot +
      geom_hline(yintercept = -log10(genomewideline), color = "red",
                 linetype = "solid", size = 0.3)
    bottom.plot <- bottom.plot +
      geom_hline(yintercept = -log10(genomewideline), color = "red",
                 linetype = "solid", size = 0.3)
  }

  # If the user requests labels:
  if(!hits.label){
    # Create the labels for the top and bottom plot
    top.labels <- make.labels(data = top.data, p.name = p.name,
                              hits.label = hits.label, top.n.hits = top.n.hits)
    bot.labels <- make.labels(data = bot.data, p.name = p.name,
                              hits.label = hits.label, top.n.hits = top.n.hits)

    # Add to plots
    top.plot <- top.plot +
      geom_label_repel(data = top.labels, aes(label = label), size = 2,
                       segment.size = 0.2, point.padding = 0.2,
                       ylim = c(maxp/2, NA),
                       min.segment.length = 0, force = 2, box.padding = 0.5)

    bottom.plot <- bottom.plot +
      geom_label_repel(data = bot.labels, aes(label = label), size = 2,
                       segment.size = 0.2, point.padding = 0.2,
                       ylim = c(NA, -(maxp/2)),
                       min.segment.length = 0, force = 2, box.padding = 0.5)
  }

  # Put the two together
  p <- top.plot + bottom.plot + plot_layout(ncol=1)

  return(p)

}


