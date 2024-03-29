# pgainsim
An R-package to assess the mode of inheritance for quantitative trait loci in GWAS.
The corresponding publication has been published in Bioinformatics:

Nora Scherer, Peggy Sekula, Peter Pfaffelhuber, Pascal Schlosser, pgainsim: an R-package to assess the mode of inheritance for quantitative trait loci in GWAS, Bioinformatics, Volume 37, Issue 18, 15 September 2021, Pages 3061–3063, https://doi.org/10.1093/bioinformatics/btab150

# Introduction

After performing genome-wide association studies under an additive, recessive and
dominant model for a quantitative trait the pgainsim package allows 
to determine study specific critical values when one of the three models
is more informative than the other ones for a quantitative trait locus.
First, study specific p-gain values are simulated. Second, based on the
simulated values quantiles of the empirical density of the p-gain are computed.
Finally, the quantiles are exported or interpolated and the critical values are obtained.

# Installation
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("rmarkdown", quietly = TRUE))
    install.packages("rmarkdown")
devtools::install_github("genepi-freiburg/pgainsim",build_vignettes = TRUE)
```

# Example
```R
browseVignettes("pgainsim")
```

# Contact
If you have any issues using the package then please get in touch with Pascal Schlosser (pascal.schlosser at uniklinik-freiburg.de).
Bug reports etc are most welcome, we want the package to be easy to use for everyone!
