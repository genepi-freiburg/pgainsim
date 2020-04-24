# werden am Ende entfernt(werden wo anders als notwendig angegeben) sind nur noch hier für testzwecke
library(minpack.lm)
library(parallel)

#' Random draw of recessive and dominante pgains under no association.
#'
#'		## hier kann man einen Beschreibung einfügen, hab jetzt erstmal die Kommentare aus dem Code gesammelt. Könntest du einen zusammenhängenden englsichen Text daraus machen?
#' wenn wir einen data frame mit allen AF als output haben, haben wir bei sehr großer Anzahl der random draws bzw. vielen SNPs per trait und vielen AF potentiell wieder ein Größen-Problem
#' Simulation von Genotypen-Wahrscheinlichkeiten eines SNP mit bestimmter AF im HW-Equilibrium
#'
#' @name p_gain_simulation
#' @param AFs    Numeric vector of assumed allele frequencies.
#' @param n    Integer. The number of random draws.
#' @param snps_per_trait    Integer. The number of single nucleotide polymorphisms to be simulated per random draw of the trait (default = 1L). snps_per_trait can be increased for efficient simulation. By increasing snps_per_trait you are reducing the number of independent draws of the trait.
#' @param n_study    Integer. The number of samples per simulation / study size (default = 1000L).
#' @param cores    Integer. Amount of CPU cores used (<=1 : sequential)
#' @return pgain_AF data frame, der für jede Allelfrequenz (in erster Spalte) pgain_rec-Datenpunkte enthält 	### <- Hier return values beschreiben. habe do.call raus genommen dh wir geben eine Liste zurück und die folgenden Funktionen müssen angepasst werden###########
#'
#' @examples
#' sim_data <- p_gain_simulation(AFs=c(0.1,0.5),n=10000L,snps_per_trait=1L,n_study=1000L,cores=2L)
#'
#' @export
p_gain_simulation<-function(AFs,n=100000L,snps_per_trait=1L,n_study=1000L,cores=1L)
{

library(parallel)

if (is.null(AFs) || !is.numeric(AFs))
stop("pgainsim: Error: AFs must be a numeric.")

if (is.null(n) || !is.integer(n))
stop("pgainsim: Error: n must be an integer.")

if (is.null(snps_per_trait) || !is.integer(snps_per_trait))
stop("pgainsim: Error: snps_per_trait must be an integer.")

if (is.null(n_study) || !is.integer(n_study))
stop("pgainsim: Error: n_study must be an integer.")

if (is.null(cores) || !is.integer(cores))
stop("pgainsim: Error: cores must be an integer.")

n_traits <- ceiling(n / snps_per_trait)
if(n==(n_traits*snps_per_trait)){
print(paste0("Generating ",n," random draws by drawing ",n_traits," standardnormal trait(s) with ",snps_per_trait," SNP(s) each."))
} else{
warning("The number of random draws is not a multiple of the number of SNPs per trait. Increasing the number of random draws to match.")
n <- n_traits * snps_per_trait
print(paste0("Generating ",n," random draws by drawing ",n_traits," standardnormal trait(s) with ",snps_per_trait," SNP(s) each."))
}

pgain_AF <- vector("list",length(AFs))

for(k in 1:length(AFs)){

AF <- AFs[k]
pvals <- mclapply(X=1:n_traits,FUN=function(j){
trait <- rnorm(n_study)
pval1<-pval2<-pval3<-rep(NA, snps_per_trait)
for(i in 1:snps_per_trait)
{
snpAA <- rbinom(n=n_study,size=1, prob=(1-AF)^2)
snpAB <- rbinom(n=n_study,size=1,prob=2*AF*(1-AF))
snpBB <- rbinom(n=n_study,size=1,prob=AF^2)
correction <- rowSums(cbind(snpAA,snpAB,snpBB))
snpAA[correction==0] <- 1
snpAB[correction==0] <- 1
snpBB[correction==0] <- 1
correction <- rowSums(cbind(snpAA,snpAB,snpBB))
snpAA <- snpAA/correction
snpAB <- snpAB/correction
snpBB <- snpBB/correction

add <- snpAB + 2* snpBB
dom <- snpAB+snpBB

pval1[i]<-summary(lm(trait~dom))$coefficients[2,4]
pval2[i]<-summary(lm(trait~snpBB))$coefficients[2,4]
pval3[i]<-summary(lm(trait~add))$coefficients[2,4]
}

return(cbind(pval1,pval2,pval3))
},mc.cores=cores)

pvals <- do.call("rbind",pvals)


test2_rec<-cbind(apply(pvals[,c(1,3)],1,min), pvals[,2])
pgain_rec<-test2_rec[,1]/test2_rec[,2]

pgain_AF[[k]] <- data.frame(AF=AF,pgain_rec=pgain_rec)
}
#pgain_AF <- do.call(rbind,pgain_AF)  
invisible(pgain_AF)
}




##2. Funktion p.gain.quantile: Quantile der pgain-Dichte bestimmen für #Tests=1 bis für #Tests=numb_tests

#Input der Funktion p_gain_quantile: sim_data, n_tests;
#pgain_simulation ist ein data frame, der für jede Allelfrequenz pgain_rec-Datenpunkte enthält (Output der Funktion p_gain_simulation)
#n_tests ist die Anzahl der Tests, bis zu der die zugehörigen Quantile der p-gain-Dichte erechnet werden sollen, z.B. n_tests=100
#Output der Funktion p_gain_quantile: data frame, der für jede Allelfrequenz in den Spalten Quantile der p-gain-Dichte für #tests=1 bis #tests=numb_tests enthält




### plural für den Funktionsnamen?

#' Computation of pgain-quantiles for numbers of tests based on the result of function p_gain_simulation.
#'
#' @name p_gain_quantile
#' @param n_tests    Integer. The number of tests for which the pgain-quantile should be computed. It depends on the available number of datapoints.
#' @param sim_data    data frame. Columns describe allele frequency (first column) and corresponding pgain (second column). Output of function p_gain_simulation.
#' @return invisible null
#'
#' @examples
#' sim_data <- p_gain_simulation(AFs=c(0.1,0.5),n=10000L,snps_per_trait=1L,n_study=1000L,cores=2L)   # Beispiele müssen self contained sein
#' pgain_quantile <- p_gain_quantile(n_tests=50L,sim_data)
#'
#' @export
p_gain_quantile<-function(n_tests, sim_data)
{

if (is.null(n_tests) || !is.numeric(n_tests))
stop("pgainsim: Error: n_tests must be an integer.")

if (is.null(sim_data) || !is.data.frame(sim_data) || !ncol(sim_data)==2 || !is.numeric(sim_data[,1]) || !is.numeric(sim_data[,2]))
stop("pgainsim: Error: sim_data must be a data frame with two columns with numeric values.")

AFs <- unique(sim_data[,1])

n_datapoints <- nrow(sim_data[sim_data[,1]==AFs[1],])


if (0.05/n_tests*n_datapoints<10)  #####hier weiß ich nicht, welche Zahl als Grenze sinnvoll wäre, damit die Quantile nicht so verzerrt sind nahe bei n_tests. 				Ich auch nicht :-) probierst du mal ein bisschen rum? 10, 100, 1000?
warning(paste0("n_tests (=",ntests,") for which the p-gain quantile should be computed is too large compared to the available number of datapoints. Reduce n_tests."))

number_tests <- c(1:n_tests)

alphas <- 0.05/number_tests


list_pgain_quantile <- vector("list",length=length(AFs))



for(i in 1:length(AFs)){

AF <- AFs[i]

res_tmp <- sim_data[sim_data[,1]==AF,]


pgain.order_rec<-res_tmp$pgain_rec[order(res_tmp$pgain_rec, decreasing=TRUE)]


list_pgain_quantile[[i]] <- as.data.frame(pgain.order_rec[floor(length(pgain.order_rec)*alphas)])

}


pgain_quantiles <- do.call(cbind,list_pgain_quantile)


colnames(pgain_quantiles) <- as.character(AFs)

invisible(pgain_quantiles)

}




##3. Funktion p_gain_quantile_fit: Fitten der Datenpunkte #tests-Quantil der p-gain-Dichte mittels log-linearer Funktion und dazugehöriger Plot und Auswertung des fits für bestimmte Anzahl von Tests
#Anmerkungen: die ersten 5 Datenpunkte werden ignoriert, damit der fit durch log(1)=0 nicht verzerrt wird
#			  falls es sehr viele Datenpunkte #tests-Quantil gibt (d.h. n_tests aus der vorigen Funktion p_gain_quantile sehr groß ist), so ist ein Ausdünnen der Datenpunkte sinnvoll, wofür die Variable n_data_ff benutzt werden kann

#Input der Funktion p_gain_quantile_fit: pgain_quantile, n_data_ff, start_vec, test_number
#pgain_quantile ist ein data frame, der für jede Allelfrequenz Quantile der p-gain-Dichte für #tests=1 bis #tests=n_tests als Vektor enthält (Output der Funktion p_gain_quantile)
#start_vec ist numeric vector, der die Startwerte für den log-linearen fit enthält
	#Problem: es wäre geschickt, wenn man für jede Allelfrequenz einzeln die Startwerte für den log-linearen fit wählen könnte, aber dann wird Funktion unübersichtlich
#Output der Funktion p_gain_quantile_fit: zum einen Liste, die für jede Allelfrequenz den log-linearen Fit der Datenpunkte #tests-Quantil der p-gain-Dichte enthält und zum anderen den Plot von Datenpunkten und Fit für jede Allelfrequenz erstellt und Auswertung des fits für test_number viele Tests



#' Computation of log-linear fit of the pgain-quantiles (dependent on the number of tests) and evaluation for a determined number of tests
#'
#' @name p_gain_quantile_fit
#' @param pgain_quantile    data frame. Columns describe pgain-quantiles for different allele frequencies (numeric values) and rows discribe number of tests. Output of function p_gain_quantile.
#' @param n_data_ff    Integer. Number of quantile datapoints that should be used for the fit. n_data_ff is a divider of the number of available datapoints (default = nrow(pgain_quantile)).
#' @param start_vec    Numeric vector of starting estimates for the log-linear fit of length 3 (default = c(1,0.5,1.2)).
#' @param test_number    Integer. Number of tests for which the p-gain threshold should be determined.
#' @return invisible null and plot of log-linear fit of the quantiles and approximated quantile for test_number many tests (wie schreibt man das?)
#'
#' @examples			### beim Testen des Beispiels (aktuell mit eingeschobenem sim_data <- do.call(rbind,sim_data)) wirft er mir beim Fitten noch NaN warnings aus.
#' sim_data <- p_gain_simulation(AFs=c(0.1,0.5),n=10000L,snps_per_trait=1L,n_study=1000L,cores=2L)   # Beispiele müssen self contained sein
#' pgain_quantile <- p_gain_quantile(n_tests=50L,sim_data)
#' fits <- p_gain_quantile_fit(pgain_quantile,n_data_ff=nrow(pgain_quantile),start_vec=c(1,1,1.2),test_number=200L)
#'
#' @export
p_gain_quantile_fit<-function(pgain_quantile,n_data_ff=nrow(pgain_quantile),start_vec = c(1,0.5,1.2),test_number)
{

if (is.null(pgain_quantile) || !is.data.frame(pgain_quantile) || !sum(!(sapply(pgain_quantile,class))=="numeric")==0)
stop("pgainsim: Error: pgain_quantile must be a data frame with numeric values in each column.")

if (is.null(n_data_ff) || !is.integer(n_data_ff) || !nrow(pgain_quantile)%%n_data_ff==0)
stop("pgainsim: Error: n_data_ff must be an integer and a divider of the number of available datapoints.")

if (is.null(start_vec) || !is.numeric(start_vec) || !length(start_vec)==3)
stop("pgainsim: Error: start_vec must be a numeric vector of length 3.")

if (is.null(test_number) || !is.integer(test_number))
stop("pgainsim: Error: test_number must be an integer.")


AFs <- as.numeric(colnames(pgain_quantile))

number_tests <- c(1:nrow(pgain_quantile))

fits <- vector("list",length=length(AFs))

for(i in 1:length(AFs)){

AF <- AFs[i]

a <- nrow(pgain_quantile)/n_data_ff

quantile_red <- pgain_quantile[,i][seq(6,nrow(pgain_quantile),a)]
number_tests_red <- number_tests[seq(6,nrow(pgain_quantile),a)]


fits[[i]] <- nlsLM(formula=quantile_red~log(a+b*number_tests_red,base=d),start=list(a=start_vec[1],b=start_vec[2],d=start_vec[3]), control=nls.control(maxiter = 1000))
pdf(paste0("Plot_#tests_vs_pgain_quantile_",AF,"_log-linear_fit.pdf"))
plot(number_tests_red,quantile_red, main=paste0("#tests vs pgain-quantile ", AF," log-linear fit"))
lines(number_tests_red, predict(fits[[i]]), col="red", type="l")
dev.off()

print_text <- paste0("p-gain-threshold for ",test_number," tests, allele frequency ",AF)
print(print_text)
print(log(coef(fits[[i]])[1] + coef(fits[[i]])[2]*test_number,base=coef(fits[[i]])[3]))


}

invisible(fits)

}