## ----setup, cache = F, echo = FALSE-------------------------------------------
knitr::opts_chunk$set(error = TRUE)

## ----install_package----------------------------------------------------------
devtools::install_github("genepi-freiburg/pgainsim")

## ----simulation---------------------------------------------------------------
library(pgainsim)
sim_dat <- p_gain_simulation(pgain_types=c("add","rec"),AFs=c(0.3,0.7),n=10000L, cores=4L)

## ----density_plot-------------------------------------------------------------
p_gain_density_plot(pgain_type="rec", sim_data=sim_dat, xlim=c(0,5))

## ----quantiles----------------------------------------------------------------
quantile <- p_gain_quantiles(n_tests=5L,sim_data=sim_dat)

## ----quantiles_fit------------------------------------------------------------
list_fits <- p_gain_quantile_fit(pgain_quantile=quantile,test_number=50L)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
warnings()

