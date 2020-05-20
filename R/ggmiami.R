#' Create a Miami plot ggplot2 object.
#'
#' @param data A data.frame object. Required.
#' @param split_by A character vector. The name of the column to use for
#'   splitting into upper and lower sections of the Miami plot. Required.
#' @param split_at A character or numeric vector. If numeric, the upper plot
#'   will contain your results where the values in the \code{split_by} column
#'   are >= \code{split_at}. If character, the upper plot will contain your
#'   results where the values in the \code{split_by} column are equal to
#'   \code{split_at}. Required.
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
#' @param upper_ylab What would you like the y-axis title for the upper plot to
#'   be? Defaults to "-log10(p)". Text added here will result in a two-line axis
#'   with your text on top and "-log10(p)" below.
#' @param lower_ylab Same as \code{upper_ylab}, but for the lower plot.
#' @param genome_line Should we draw a genome-wide significance line? Set to
#'   NULL if you do not want this line. Defaults to 5e-8, or supply your own
#'   number.
#' @param genome_line_color What color should your genome-wide significance line
#'   be? Defaults to red.
#' @param suggestive_line Should we draw a suggestive significance line? Set to
#'   NULL if you do not want this line. Defaults to 1e-5, or supply your own
#'   number.
#' @param suggestive_line_color What color should your suggestive significance
#'   line be? Defaults to blue.
#' @param hits_label_col Either the name of the column(s), max = 2, to use for
#'   automatically labeling n hits, defined by \code{top_n_hits}, or the
#'   column where the values you provide in \code{hits_label} can be found.
#'   Defaults to NULL: labels aren't displayed.
#' @param hits_label A user-specified character vector of probes/genes/SNPs
#'   to label. Defaults to NULL.
#' @param top_n_hits How many of the top hits do you want to label? Defaults to
#'   5. Set to NULL to turn off this filtering.
#' @examples
#'   If you want to put SNPs with positive beta values in the upper plot and
#'   netative beta values in the lower plot:
#'   ggmiami(data = df, split_by = "beta", split_at = 0)
#'
#'   If you want to put results from study A in the upper plot and study B
#'   in the lower plot:
#'   ggmiami(data = df, split_by = "study", split_at = "A")
#'
#'   If you want to add a genome wide line at 9e-8, instead of 5e-8:
#'   ggmiami(data = df, split_by = "study", split_at = "A", genome_line = 9e-8)
#'
#'   If you want to label the top 5 hits with the SNP and gene name:
#'   ggmiami(data = df, split_by = "study", split_at = "A",
#'           hits_label_col = c("SNP", "Gene"))
#'
#'   If you want to label the top 10 hits from the AHRR, TBX15, and RARA genes:
#'   ggmiami(data = df, split_by = "study", split_at = "A",
#'           hits_label_col = "Gene", hits_label = c("AHRR", "TBX15", "RARA"),
#'           top_n_hits = 10)
#'
#'   When labeling, the items in hits_label must all come from a single column.
#'   Specifying multiple columns will return an error.
#'
#'   \dontrun{
#'   ggmiami(data = df, split_by = "study", split_at = "A",
#'           hits_label_col = c("chr", Gene"),
#'           hits_label = c("1", AHRR", "2", "TBX15", "RARA"),
#'           top_n_hits = 10)
#'   }
#'
#' @export
#' @return a \code{ggplot2} object
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#' @keywords manhattan miami plot SNP genetics
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 aes
#'

ggmiami <- function(
  data,
  split_by,
  split_at,
  chr = "chr",
  pos = "pos",
  p  =  "p",
  chr_colors = c("black", "grey"),
  upper_ylab = "-log10(p)",
  lower_ylab = "-log10(p)",
  genome_line = 5e-8,
  genome_line_color = "red",
  suggestive_line = 1e-5,
  suggestive_line_color = "blue",
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

  # Prepare the axis titles
  if (upper_ylab == "-log10(p)") {
    upper_ylab <- expression("-log" [10]* "(p)")
  } else {
    upper_ylab <- bquote(atop(.(upper_ylab), "-log" [10]* "(p)"))
    # upper_ylab <- bquote(atop(NA, atop(textstyle(.(upper_ylab)),
                                 # textstyle("-log" [10]* "(p)"))))

  }

  # Prepare the axis titles
  if (lower_ylab == "-log10(p)") {
    lower_ylab <- expression("-log" [10]* "(p)")
  } else {
    lower_ylab <- bquote(atop(.(lower_ylab), "-log" [10]* "(p)"))
    # lower_ylab <- bquote(atop(atop(textstyle(.(lower_ylab)),
                                    # textstyle("-log" [10]* "(p)")), NA))
  }

  # When ggtext is published on CRAN, this is a better solution for axis text.
  # lower_ylab <- paste0(lower_ylab, "<br>-log<sub>10</sub>(p)")
  # Then in ggplot2 call: axis.title.y = ggtext::element_markdown(),

  # Create base upper plot.
  upper_plot <- ggplot2::ggplot(data = plot_data$upper,
                              aes(x = .data$rel_pos, y = .data$logged_p)) +
    ggplot2::geom_point(aes(color = as.factor(.data$chr)), size = 0.25) +
    ggplot2::scale_color_manual(values = chr_colors) +
    ggplot2::scale_x_continuous(labels = plot_data$axis$chr,
                                breaks = plot_data$axis$chr_center,
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, plot_data$maxp),
                                expand =
                                  ggplot2::expansion(mult = c(0.02, 0))) +
    ggplot2::labs(x = "", y = upper_ylab) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(b = 0, l = 10))

  # Create base lower plot
  lower_plot <- ggplot2::ggplot(data = plot_data$lower,
                                 aes(x = .data$rel_pos, y = .data$logged_p)) +
    ggplot2::geom_point(aes(color = as.factor(.data$chr)), size = 0.25) +
    ggplot2::scale_color_manual(values = chr_colors) +
    ggplot2::scale_x_continuous(breaks = plot_data$axis$chr_center,
                                position = "top",
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_reverse(limits = c(plot_data$maxp, 0),
                             expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = "", y = lower_ylab) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 0, l = 10))

  # If the user has requested a suggetive line, add:
  if (!is.null(suggestive_line)) {
    upper_plot <- upper_plot +
      ggplot2::geom_hline(yintercept = -log10(suggestive_line),
                          color = suggestive_line_color,
                          linetype = "solid", size = 0.4)
    lower_plot <- lower_plot +
      ggplot2::geom_hline(yintercept = -log10(suggestive_line),
                          color = suggestive_line_color,
                          linetype = "solid", size = 0.4)
  }

  # If the user has requested a genome-wide line, add:
  if (!is.null(genome_line)) {
    upper_plot <- upper_plot +
      ggplot2::geom_hline(yintercept = -log10(genome_line),
                          color = genome_line_color,
                          linetype = "dashed", size = 0.4)
    lower_plot <- lower_plot +
      ggplot2::geom_hline(yintercept = -log10(genome_line),
                          color = genome_line_color,
                          linetype = "dashed", size = 0.4)
  }

  # If the user requests labels:
  if (!is.null(hits_label_col)) {
    # Create the labels for the upper and lower plot
    upper_labels <- make_miami_labels(data = plot_data$upper,
                                    hits_label_col = hits_label_col,
                                    hits_label = hits_label,
                                    top_n_hits = top_n_hits)

    lower_labels <- make_miami_labels(data = plot_data$lower,
                                    hits_label_col = hits_label_col,
                                    hits_label = hits_label,
                                    top_n_hits = top_n_hits)

    # Add to plots
    upper_plot <- upper_plot +
      ggrepel::geom_label_repel(data = upper_labels, aes(label = label),
                                size = 2, segment.size = 0.2,
                                point.padding = 0.2,
                                ylim = c(plot_data$maxp / 2, NA),
                                min.segment.length = 0,
                                force = 2, box.padding = 0.5)

    lower_plot <- lower_plot +
      ggrepel::geom_label_repel(data = lower_labels, aes(label = label),
                                size = 2, segment.size = 0.2,
                                point.padding = 0.2,
                                ylim = c(NA, -(plot_data$maxp / 2)),
                                min.segment.length = 0,
                                force = 2, box.padding = 0.5)
  }

  # Put the two together
  p <- upper_plot + lower_plot + patchwork::plot_layout(ncol = 1)

  return(p)

}
