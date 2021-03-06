% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgainsim.R
\name{p_gain_quantiles}
\alias{p_gain_quantiles}
\title{Computation of p-gain-quantiles for numbers of tests based on the result of function p_gain_simulation.}
\usage{
p_gain_quantiles(sim_data, pgain_types = names(sim_data), n_tests)
}
\arguments{
\item{sim_data}{list of data frames. For every type of p-gain there is a data frame (in the same order as they are listed in pgain_types). Each column describes simulated p-gain-values for the allele frequency, which is in the column name. Output of function p_gain_simulation.}

\item{pgain_types}{Character vector of the types of p-gains, which were simulated in the function p_gain_simulation and from which quantiles should be computed. ("add" means additive p-gain, "rec" means recessive p-gain, "dom" means dominant p-gain) (default = names(sim_data)).}

\item{n_tests}{Integer. The number of tests for which the p-gain-quantile should be computed. It depends on the available number of datapoints.}
}
\value{
list_pgains_quant list of data frames. For every type of p-gain there is a data frame (in the same order as they are listed in pgain_types). Columns describe p-gain-quantiles for different allele frequenciey (numeric values) and rows discribe number of tests. For the additive case simulated p-gain values of allele frequencies AF and 1-AF are combined before computing the quantiles.
}
\description{
Computation of p-gain-quantiles for numbers of tests based on the result of function p_gain_simulation.
}
\examples{
sim_dat <- p_gain_simulation(pgain_types=c("add","rec"),AFs=c(0.3,0.7),n=10000L, cores=4L)
quantile <- p_gain_quantiles(n_tests=5L,sim_data=sim_dat)

}
