#' Create a Miami plot ggplot2 object.
#'
#' @param data A data.frame object. Required.
#' @param split.by A character vector. The name of the column to use for
#'   splitting into top and bottom sections of the Miami plot. Required.
#' @param split.at A character or numeric vector. If numeric, the top plot will
#'   contain your results where the values in the \code{split.by} column are
#'   >= \code{split.at}. If character, top plot will contain your results where
#'   the values in the \code{split.by} column are equal to \code{split.at}. Required.
#' @param chr The name of the column containing your chromosome information.
#'   Defaults to "chr"
#' @param pos The name of the column containing your position information.
#'   Defaults to "pos"
#' @param p The name of the column containing your p-value information.
#'   Defaults to "p"
#' @param chrcolors Either a vector of two colors to alternate across chromosomes
#'   or a vector of colors to use for coloring chromosomes, with length equal to
#'   the number of chromosomes being plotted. Defaults to alternating black and
#'   grey.
#' @param genomewideline Should we draw a genome-wide significance line? Set to
#'   NULL if you do not want this line. Defaults to 5e-8, or supply your own number.
#' @param suggestiveline Should we draw a suggestive significance line? Set to
#'   NULL if you do not want this line. Defaults to 1e-5, or supply your own number.
#' @param hits.label.col Either the name of the column(s), max = 2, to use for
#'   automatically labeling n hits, determined using \code{top.n.hits}, or the
#'   column where the values you provide in \code{hits.label} can be found.
#'   Defaults to NULL: labels aren't displayed.
#' @param hits.label A user-specified character vector of probes/genes/SNPs
#'   to label. Defaults to NULL.
#' @param top.n.hits How many of the top hits do you want to label? Defaults to
#'   5. Set to NULL to turn off this filtering.
#' @export
#' @return a \code{ggplot2} object
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#' @keywords manhattan miami plot SNP genetics
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 aes

ggmiami <- function(
  data,
  split.by,
  split.at,
  chr="chr",
  pos="pos",
  p = "p",
  chrcolors = c("black", "grey"),
  genomewideline = 5e-8,
  suggestiveline = 1e-5,
  hits.label.col = NULL,
  hits.label = NULL,
  top.n.hits = 5) {

  # Prepare the data
  plot.data <- prep_miami_data(data = data, split.by = split.by,
                               split.at = split.at, chr = chr, pos = pos, p = p)

  # Prepare the colors
  if(length(chrcolors) == 2){
    chrcolors <- rep(chrcolors, length.out = nrow(plot.data$axis))
  } else if (length(chrcolors) == nrow(plot.data$axis)){
    chrcolors <- chrcolors
  } else {
    stop("The number of colors specified in {chrcolors} does not match the
         number of chromosomes to be displayed.")
  }

  # Create base top plot.
  top.plot <- ggplot2::ggplot(data = plot.data$top,
                              aes(x=.data$rel.pos, y=.data$loggedp)) +
    ggplot2::geom_point(aes(color=as.factor(.data$chr)), size=0.25) +
    ggplot2::scale_color_manual(values = chrcolors) +
    ggplot2::scale_x_continuous(labels = plot.data$axis$chr,
                                breaks = plot.data$axis$chr.center,
                                expand = ggplot2::expansion(mult=0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, plot.data$maxp),
                                expand = ggplot2::expansion(mult = c(0.02,0))) +
    ggplot2::labs(x = "", y = expression('-log'[10]*'(P)')) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(b=0))

  # Create base bottom plot
  bottom.plot <- ggplot2::ggplot(data = plot.data$bottom,
                                 aes(x=.data$rel.pos, y=.data$loggedp)) +
    ggplot2::geom_point(aes(color=as.factor(.data$chr)), size=0.25) +
    ggplot2::scale_color_manual(values = chrcolors) +
    ggplot2::scale_x_continuous(breaks = plot.data$axis$chr.center,
                                position = "top",
                                expand = ggplot2::expansion(mult=0.01)) +
    ggplot2::scale_y_reverse(limits = c(plot.data$maxp, 0),
                             expand = ggplot2::expansion(mult = c(0,0.02))) +
    ggplot2::labs(x = "", y = expression('-log'[10]*'(P)')) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t=0))

  # If the user has requested a suggetive line, add:
  if(!is.null(suggestiveline)){
    top.plot <- top.plot +
      ggplot2::geom_hline(yintercept = -log10(suggestiveline), color = "blue",
                          linetype = "solid", size = 0.3)
    bottom.plot <- bottom.plot +
      ggplot2::geom_hline(yintercept = -log10(suggestiveline),color = "blue",
                          linetype = "solid", size = 0.3)
  }

  # If the user has requested a genome-wide line, add:
  if(!is.null(genomewideline)){
    top.plot <- top.plot +
      ggplot2::geom_hline(yintercept = -log10(genomewideline), color = "red",
                          linetype = "dashed", size = 0.3)
    bottom.plot <- bottom.plot +
      ggplot2::geom_hline(yintercept = -log10(genomewideline), color = "red",
                          linetype = "dashed", size = 0.3)
  }

  # If the user requests labels:
  if(!is.null(hits.label)){
    # Create the labels for the top and bottom plot
    top.labels <- make_miami_labels(data = plot.data$top,
                                    hits.label.col = hits.label.col,
                                    hits.label = hits.label,
                                    top.n.hits = top.n.hits)

    bot.labels <- make_miami_labels(data = plot.data$bottom,
                                    hits.label.col = hits.label.col,
                                    hits.label = hits.label,
                                    top.n.hits = top.n.hits)

    # Add to plots
    top.plot <- top.plot +
      ggrepel::geom_label_repel(data = top.labels, aes(label = label), size = 2,
                                segment.size = 0.2, point.padding = 0.2,
                                ylim = c(plot.data$maxp/2, NA),
                                min.segment.length = 0,
                                force = 2, box.padding = 0.5)

    bottom.plot <- bottom.plot +
      ggrepel::geom_label_repel(data = bot.labels, aes(label = label), size = 2,
                                segment.size = 0.2, point.padding = 0.2,
                                ylim = c(NA, -(plot.data$maxp/2)),
                                min.segment.length = 0,
                                force = 2, box.padding = 0.5)
  }

  # Put the two together
  p <- top.plot + bottom.plot + patchwork::plot_layout(ncol=1)

  return(p)

}


