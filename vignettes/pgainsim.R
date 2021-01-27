## ----setup, cache = F, echo = FALSE-------------------------------------------
knitr::opts_chunk$set(error = TRUE)

## ----install_package----------------------------------------------------------
if (!requireNamespace("devtools", quietly = TRUE))
   install.packages("devtools")
if (!requireNamespace("pgainsim", quietly = TRUE))
  devtools::install_github("genepi-freiburg/pgainsim",build_vignettes = TRUE)

## ----simulation---------------------------------------------------------------
require(pgainsim)
sim_dat <- p_gain_simulation(pgain_types=c("add","rec"), AFs=c(0.3,0.5,0.7,0.9), n=10000L, cores=1L)

## ----density_plot-------------------------------------------------------------
p_gain_density_plot(pgain_type="rec", sim_data=sim_dat, xlim=c(0,5), print_pdf=FALSE)

## ----quantiles----------------------------------------------------------------
quantile <- p_gain_quantiles(n_tests=5L, sim_data=sim_dat)

## ----quantiles_fit------------------------------------------------------------
list_fits <- p_gain_quantile_fit(pgain_quantile=quantile, test_number=50L, print_pdf=FALSE)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
warnings()

