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
#'   to label. Defaults to NULL. If you specify this, you must also specify
#'   hits_labels_col so that we know where to look for the items.
#' @param top_n_hits How many of the top hits do you want to label? Defaults to
#'   5. Set to NULL to turn off this filtering.
#' @param upper_labels_df A dataframe containing the relative position, logged
#'   p-value, and label to use for labelling the upper plot. Column names should
#'   be c("rel_pos", "logged_p", "label"). Defaults to NULL.
#' @param lower_labels_df Same as \code{upper_labels_df} but for the lower plot.
#' @param upper_highlight A vector of SNPs or gene names you would like to
#'   highlight in the upper plot. Defaults to NULL. If you specify this, you
#'   must also specify upper_highlight_col so that we know where to find the
#'   items.
#' @param upper_highlight_col The column where the values you provide in
#'   \code{upper_highlight} can be found. Defaults to NULL.
#' @param upper_highlight_color The color that you would like the highlighted
#'   points in the upper plot to be. Defaults to "green."
#' @param lower_highlight Same as \code{upper_highlight} but for the lower plot.
#' @param lower_highlight_col Same as \code{upper_highlight_col} but for the
#'   lower plot.
#' @param lower_highlight_color Same as \code{lower_highlight_color} but for the
#'   lower plot.
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
  top_n_hits = 5,
  upper_labels_df = NULL,
  lower_labels_df = NULL,
  upper_highlight = NULL,
  upper_highlight_col = NULL,
  upper_highlight_color = "green",
  lower_highlight = NULL,
  lower_highlight_col = NULL,
  lower_highlight_color = "green") {

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
  # Then in ggplot2 call: axis.title.y = ggtext::element_markdown() and un-set
  # the l margin.

  # Create base upper plot.
  upper_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = plot_data$upper,
                        aes(x = .data$rel_pos, y = .data$logged_p,
                            color = as.factor(.data$chr)), size = 0.25) +
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
  lower_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = plot_data$lower,
                        aes(x = .data$rel_pos, y = .data$logged_p,
                            color = as.factor(.data$chr)), size = 0.25) +
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

  # If they try to specify both hits_label_col and dataframe(s), return an
  # error message.
  if (all(!is.null(hits_label_col), any(!is.null(upper_labels_df),
                                        !is.null(lower_labels_df)))) {
    stop("You have specified both hits_label_col and a *_labels_df. This
         package does not know how to use both information simultaneously.
         Please only use one method for labelling: either hits_label_col (with
         or without hits_label), or *_labels_df.")
  }

  # If the user requests labels based on the hits_label_col
  if (all(!is.null(hits_label_col), is.null(upper_labels_df),
          is.null(lower_labels_df))) {
    # Create the labels for the upper and lower plot
    upper_labels_df <- make_miami_labels(data = plot_data$upper,
                                         hits_label_col = hits_label_col,
                                         hits_label = hits_label,
                                         top_n_hits = top_n_hits)

    lower_labels_df <- make_miami_labels(data = plot_data$lower,
                                         hits_label_col = hits_label_col,
                                         hits_label = hits_label,
                                         top_n_hits = top_n_hits)

    # Add to the plots
    upper_plot <- upper_plot +
      ggrepel::geom_label_repel(data = upper_labels_df,
                                aes(x = .data$rel_pos,
                                    y = .data$logged_p,
                                    label = .data$label),
                                size = 2, segment.size = 0.2,
                                point.padding = 0.3,
                                ylim = c(plot_data$maxp / 2, NA),
                                min.segment.length = 0, force = 2,
                                box.padding = 0.5)

    lower_plot <- lower_plot +
      ggrepel::geom_label_repel(data = lower_labels_df,
                                aes(x = .data$rel_pos,
                                    y = .data$logged_p,
                                    label = .data$label),
                                size = 2, segment.size = 0.2,
                                point.padding = 0.3,
                                ylim = c(NA, -(plot_data$maxp / 2)),
                                min.segment.length = 0, force = 2,
                                box.padding = 0.5)
  }

  # If the user requests labels based on a dataframe for the upper plot
  if (all(is.null(hits_label_col), !is.null(upper_labels_df))) {
      # Make sure that the column names in upper_labels_df are correct
      checkmate::assertNames(colnames(upper_labels_df),
                             identical.to = c("rel_pos", "logged_p", "label"))

      # Add to plot
      upper_plot <- upper_plot +
        ggrepel::geom_label_repel(data = upper_labels_df,
                                  aes(x = .data$rel_pos,
                                      y = .data$logged_p,
                                      label = .data$label),
                                  size = 2, segment.size = 0.2,
                                  point.padding = 0.3,
                                  ylim = c(plot_data$maxp / 2, NA),
                                  min.segment.length = 0, force = 2,
                                  box.padding = 0.5)
  }

  # If the user requests labels based on a dataframe for the lower plot
  if (all(is.null(hits_label_col), !is.null(lower_labels_df))) {
      # Make sure that the column names in lower_labels_df are correct
      checkmate::assertNames(colnames(lower_labels_df),
                             identical.to = c("rel_pos", "logged_p", "label"))

      # Add to plot
      lower_plot <- lower_plot +
        ggrepel::geom_label_repel(data = lower_labels_df,
                                  aes(x = .data$rel_pos,
                                      y = .data$logged_p,
                                      label = .data$label),
                                  size = 2, segment.size = 0.2,
                                  point.padding = 0.3,
                                  ylim = c(NA, -(plot_data$maxp / 2)),
                                  min.segment.length = 0, force = 2,
                                  box.padding = 0.5)
  }

  # Highlight upper snps
  if (all(!is.null(upper_highlight), !is.null(upper_highlight_col))) {
    # Make highlighting df
    upper_highlight_df <- highlight_miami(data = plot_data$upper,
                                          highlight = upper_highlight,
                                          highlight_col = upper_highlight_col,
                                          highlight_color =
                                            upper_highlight_color)

    # Add to plot
    upper_plot <- upper_plot +
      ggplot2::geom_point(data = upper_highlight_df,
                          aes(x = .data$rel_pos, y = .data$logged_p),
                          color = upper_highlight_df$color,
                          size = 0.25)
  }

  # Highlight lower snps
  if (all(!is.null(lower_highlight), !is.null(lower_highlight_col))) {
    # Make highlighting df
    lower_highlight_df <- highlight_miami(data = plot_data$lower,
                                          highlight = lower_highlight,
                                          highlight_col = lower_highlight_col,
                                          highlight_color =
                                            lower_highlight_color)

    # Add to plot
    lower_plot <- lower_plot +
      ggplot2::geom_point(data = lower_highlight_df,
                          aes(x = .data$rel_pos, y = .data$logged_p),
                          color = lower_highlight_df$color,
                          size = 0.25)
  }

  # Put the two together
  p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)

  return(p)

}
