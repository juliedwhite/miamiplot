# MiamiPlot
 R function to create a ggplot2 based miami plot

## Installation
Install directly from github using devtools
```coffee
library(devtools)
install_github("juliedwhite/miamiplot")
```
Load the package each time you use it:
```coffee
library(miamiplot)
```

## Usage
See the vignette for more examples:
```coffee
vignette("miamiplot")
```

Take a look at the built-in data:
```coffee
head(gwas_results)
```

Basic miami plot using built-in data, with positive beta values in the upper
plot and negative beta values in the lower plot.
```coffee
ggmiami(data = gwas_results[which(gwas_results$study == "A"),], 
        split_by = "beta", split_at = 0, p = "pval", 
        upper_ylab = "Positive beta values",
        lower_ylab = "Negative beta values")
```
