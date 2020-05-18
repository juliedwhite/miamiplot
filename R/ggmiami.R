#' Create a Miami plot ggplot2 object.
#'
#' @param data A data.frame object. Required.
#' @param split_by A character vector. The name of the column to use for
#'   splitting into top and bottom sections of the Miami plot. Required.
#' @param split_at A character or numeric vector. If numeric, the top plot will
#'   contain your results where the values in the \code{split_by} column are
#'   >= \code{split_at}. If character, top plot will contain your results where
#'   the values in the \code{split_by} column are equal to \code{split_at}.
#'   Required.
#' @param chr The name of the column containing your chromosome information.
#'   Defaults to "chr"
#' @param pos The name of the column containing your position information.
#'   Defaults to "pos"
#' @param p The name of the column containing your p-value information.
#'   Defaults to "p"
#' @param chr_colors Either a vector of two colors to alternate across
#'   chromosomes or a vector of colors to use for coloring chromosomes, with
#'   length equal to the number of chromosomes being plotted. Defaults to
#'   alternating black and grey.
#' @param genome_line Should we draw a genome-wide significance line? Set to
#'   NULL if you do not want this line. Defaults to 5e-8, or supply your own
#'   number.
#' @param suggestive_line Should we draw a suggestive significance line? Set to
#'   NULL if you do not want this line. Defaults to 1e-5, or supply your own
#'   number.
#' @param hits_labe_lcol Either the name of the column(s), max = 2, to use for
#'   automatically labeling n hits, determined using \code{top_n_hits}, or the
#'   column where the values you provide in \code{hits_label} can be found.
#'   Defaults to NULL: labels aren't displayed.
#' @param hits_label A user-specified character vector of probes/genes/SNPs
#'   to label. Defaults to NULL.
#' @param top_n_hits How many of the top hits do you want to label? Defaults to
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
  split_by,
  split_at,
  chr = "chr",
  pos = "pos",
  p  =  "p",
  chr_colors = c("black", "grey"),
  genome_line = 5e-8,
  suggestive_line = 1e-5,
  hits_label_col = NULL,
  hits_label = NULL,
  top_n_hits = 5) {

  # Prepare the data
  plot_data <- prep_miami_data(data = data, split_by = split_by,
                               split_at = split_at, chr = chr, pos = pos, p = p)

  # Prepare the colors
  if (length(chr_colors) == 2) {
    chr_colors <- rep(chr_colors, length.out = nrow(plot_data$axis))
  } else if (length(chr_colors) == nrow(plot_data$axis)) {
    chr_colors <- chr_colors
  } else {
    stop("The number of colors specified in {chr_colors} does not match the
         number of chromosomes to be displayed.")
  }

  # Create base top plot.
  top_plot <- ggplot2::ggplot(data = plot_data$top,
                              aes(x = .data$rel_pos, y = .data$loggedp)) +
    ggplot2::geom_point(aes(color = as.factor(.data$chr)), size = 0.25) +
    ggplot2::scale_color_manual(values = chr_colors) +
    ggplot2::scale_x_continuous(labels = plot_data$axis$chr,
                                breaks = plot_data$axis$chr_center,
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, plot_data$maxp),
                                expand =
                                  ggplot2::expansion(mult = c(0.02, 0))) +
    ggplot2::labs(x = "", y = expression("-log" [10]* "(P)")) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(b = 0))

  # Create base bottom plot
  bottom_plot <- ggplot2::ggplot(data = plot_data$bottom,
                                 aes(x = .data$rel_pos, y = .data$loggedp)) +
    ggplot2::geom_point(aes(color = as.factor(.data$chr)), size = 0.25) +
    ggplot2::scale_color_manual(values = chr_colors) +
    ggplot2::scale_x_continuous(breaks = plot_data$axis$chr_center,
                                position = "top",
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_reverse(limits = c(plot_data$maxp, 0),
                             expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = "", y = expression("-log" [10]* "(P)")) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 0))

  # If the user has requested a suggetive line, add:
  if (!is.null(suggestive_line)) {
    top_plot <- top_plot +
      ggplot2::geom_hline(yintercept = -log10(suggestive_line), color = "blue",
                          linetype = "solid", size = 0.3)
    bottom_plot <- bottom_plot +
      ggplot2::geom_hline(yintercept = -log10(suggestive_line), color = "blue",
                          linetype = "solid", size = 0.3)
  }

  # If the user has requested a genome-wide line, add:
  if (!is.null(genome_line)) {
    top_plot <- top_plot +
      ggplot2::geom_hline(yintercept = -log10(genome_line), color = "red",
                          linetype = "dashed", size = 0.3)
    bottom_plot <- bottom_plot +
      ggplot2::geom_hline(yintercept = -log10(genome_line), color = "red",
                          linetype = "dashed", size = 0.3)
  }

  # If the user requests labels:
  if (!is.null(hits_label_col)) {
    # Create the labels for the top and bottom plot
    top_labels <- make_miami_labels(data = plot_data$top,
                                    hits_label_col = hits_label_col,
                                    hits_label = hits_label,
                                    top_n_hits = top_n_hits)

    bot_labels <- make_miami_labels(data = plot_data$bottom,
                                    hits_label_col = hits_label_col,
                                    hits_label = hits_label,
                                    top_n_hits = top_n_hits)

    # Add to plots
    top_plot <- top_plot +
      ggrepel::geom_label_repel(data = top_labels, aes(label = label), size = 2,
                                segment.size = 0.2, point.padding = 0.2,
                                ylim = c(plot_data$maxp / 2, NA),
                                min.segment.length = 0,
                                force = 2, box.padding = 0.5)

    bottom_plot <- bottom_plot +
      ggrepel::geom_label_repel(data = bot_labels, aes(label = label), size = 2,
                                segment.size = 0.2, point.padding = 0.2,
                                ylim = c(NA, -(plot_data$maxp / 2)),
                                min.segment.length = 0,
                                force = 2, box.padding = 0.5)
  }

  # Put the two together
  p <- top_plot + bottom_plot + patchwork::plot_layout(ncol = 1)

  return(p)

}
