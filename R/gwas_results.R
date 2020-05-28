#' Simulated data for 60,000 SNPs from two studies
#'
#' A simulated dataset containing the rsid, chromosome, position, beta value,
#'   standard error, t-statistic, p-value, and study name for 60,000 SNPs total.
#'   Each study has 30,000 SNPs spread across 22 chromosomes.
#'
#' @format A data.frame with 60,000 rows and 8 variables:
#' \describe{
#'   \item{rsid}{Variant name}
#'   \item{chr}{Chromosome}
#'   \item{pos}{Base pair position}
#'   \item{beta}{Beta value from simulated GWAS}
#'   \item{se}{Standard error from simulated GWAS}
#'   \item{tstat}{T-statistic from simulated GWAS}
#'   \item{pval}{P-value from simulated GWAS}
#'   \item{study}{Either study A or study B, to signify two GWAS of the same
#'                phenotypes that have been combined for plotting purposes.}
#' }
#'
"gwas_results"
