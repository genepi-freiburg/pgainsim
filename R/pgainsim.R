#####p-gain simulation functions for R-package

##1. Funktion p.gain.simulation: Datensimulation für 100k Datenpunkte (d.h. pgain-Werte) pro Allelfrequenz

#Input der Funktion p.gain.simulation: AFs
#AFs ist Vektor bestehend aus den gewünschten Allelfrequenzen z.B AFs <- c(seq(0.1,0.9,0.1))
#Output der Funktion p.gain.simulation: Liste, in der für jede Allelfrequenz 100k pgain_rec-Datenpunkte enthalten sind: jeder Listeneintrag enthält data frame mit Spalten Allelfrequenz AF und rezessivem p-gain pgain_rec

#' Random draw of recessive and dominante pgains under no association.
#'
#' @name p_gain_simulation
#' @param AFs    Numeric vector of assumed allele frequencies.
#' @param n    Integer. The number of random draws.
#' @param snps_per_trait    Integer. The number of single nucleotide polymorphisms to be simulated per random draw of the trait (default = 1L). snps_per_trait can be increased for efficient simulation. By increasing snps_per_trait you are reducing the number of independent draws of the trait.
#' @param n_study    Integer. The number of samples per simulation / study size (default = 1000L).
#' @param cores    Integer. Amount of CPU cores used (<=1 : sequential)
#' @return invisible null
#'
#' @examples
#' sim_data <- p_gain_simulation(AFs=c(0.1,0.5),n=10000L,snps_per_trait=1L,n_study=1000L,cores=2L)
#'
#' @export
p_gain_simulation<-function(AFs,n=100000L,snps_per_trait=1L,n_study=1000L,cores=1L)
{
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
#pval.l <- list(pval1.l,pval2.l,pval3.l)
#pval1.l<-pval2.l<-pval3.l<-NULL

pvals <- mclapply(X=1:n_traits,FUN=function(j){
#traitnkonzentration wird als normalverteilt simuliert
trait <- rnorm(n_study)
pval1<-pval2<-pval3<-rep(NA, snps_per_trait)
for(i in 1:snps_per_trait)
{
#Simulation von Genotypen-Wahrscheinlichkeiten eines SNP mit bestimmter AF im HW-Equilibrium
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
pgain_AF <- do.call(rbind,pgain_AF)

invisible(pgain_AF)
}


pgain_simulation <- p.gain.simulation(AFs)


##2. Funktion p.gain.quantile: Quantile der pgain-Dichte bestimmen für #Tests=1 bis für #Tests=numb_tests

#Input der Funktion p.gain.quantile: pgain_simulation und numb_tests;
#pgain_simulation ist eine Liste, die für jede Allelfrequenz pgain_rec-Datenpunkte enthält (Output der Funktion p.gain.simulation)
#numb_tests ist die Anzahl der Tests, bis zu der die zugehörigen Quantile der p-gain-Dichte erechnet werden sollen, z.B. numb_tests=200
#Output der Funktion p.gain.quantile: Liste, in der für jede Allelfrequenz Quantile der p-gain-Dichte für #tests=1 bis #tests=numb_tests enthalten sind, also Vektor der Länge numb_tests

p.gain.quantile<-function(pgain_simulation,numb_tests)
{

number_tests <- c(1:numb_tests)

alphas <- 0.05/number_tests


list_pgain_quantile <- vector("list",length=length(AFs))

for(i in 1:length(AFs)){

AF <- AFs[i]

res_tmp <- pgain_simulation[[i]]


pgain.order_rec<-res_tmp$pgain_rec[order(res_tmp$pgain_rec, decreasing=TRUE)]


list_pgain_quantile[[i]] <- pgain.order_rec[floor(length(pgain.order_rec)*alphas)]

}

invisible(list_pgain_quantile)

}


pgain_quantile <- p.gain.quantile(pgain_simulation,numb_tests)



##3. Funktion p.gain.quantile.fit: Fitten der Datenpunkte #tests-Quantil der p-gain-Dichte mittels log-linearer Funktion und dazugehöriger Plot
#Anmerkungen: die ersten 5 Datenpunkte werden ignoriert, damit der fit durch log(1)=0 nicht verzerrt wird
#			  falls es sehr viele Datenpunkte #tests-Quantil gibt (d.h. numb_tests aus der vorigen Funktion p.gain.quantile sehr groß ist (dafür müsste die p.gain.simulation-Funktion jedoch mehrmals angewendet und die Outputs zusammengeklebt worden sein)), so ist ein Ausdünnen der Datenpunkte sinnvoll

#Input der Funktionp.gain.quantile.fit: pgain_quantile, start_a, start_b, start_d
#pgain_quantile ist eine Liste, die für jede Allelfrequenz Quantile der p-gain-Dichte für #tests=1 bis #tests=numb_tests als Vektor enthält (Output der Funktion p.gain.quantile)
#start_a, start_b, start_d sind Zahlen als Startwerte für den log-linearen fit, z.B. start_a <- 1, start_b <- 0.5, start_d <- 1.2
	#Problem: es wäre geschickt, wenn man für jede Allelfrequenz einzeln die Startwerte für den log-linearen fit wählen könnte, aber dann wird Funktion unübersichtlich
#Output der Funktion p.gain.quantile.fit: zum einen Liste, die für jede Allelfrequenz den log-linearen Fit der Datenpunkte #tests-Quantil der p-gain-Dichte enthält und zum anderen den Plot von Datenpunkten und Fit für jede Allelfrequenz erstellt

p.gain.quantile.fit<-function(pgain_quantile,start_a,start_b,start_d)
{

library(minpack.lm)

number_tests <- c(1:length(pgain_quantile[[1]]))

fits <- vector("list",length=length(AFs))

for(i in 1:length(AFs)){

AF <- AFs[i]

quantile_red <- pgain_quantile[[i]][5:length(number_tests)]
number_tests_red <- number_tests[5:length(number_tests)]

fits[[i]] <- nlsLM(formula=quantile_red~log(a+b*number_tests_red,base=d),start=list(a=start_a,b=start_b,d=start_d), control=nls.control(maxiter = 1000))
pdf(paste0("Plot_#tests_vs_pgain_quantile_",AF,"_log-linear_fit.pdf"))
plot(number_tests_red,quantile_red, main=paste0("#tests vs pgain-quantile ", AF," log-linear fit"))
lines(number_tests_red, predict(fits[[i]]), col="red", type="l")
dev.off()


}

invisible(fits)

}



fits <- p.gain.quantile.fit(pgain_quantile,start_a,start_b,start_d)




##4. Funktion p.gain.threshold: Quantil der p-gain-Dichte für beliebige Anzahl von Tests (test_number) mittels log-linearem fit bestimmen

#Input der Funktion p.gain.threshold: fits, test_number
#fits ist eine Liste, die für jede Allelfrequenz den log-linearen fit enthält (Output der Funktion p.gain.quantile.fit)
#test_number ist eine Zahl (die Anzahl von Tests, für die man den p-gain-threshold bestimmen will), z.B. test_number <- 2.46e+04
#Output der Funktion p.gain.threshold: approximierte Quantile der p-gain-Dichte für die bestimmte Anzahl test_number von Tests für jede Allelfrequenz


p.gain.threshold<-function(fits,test_number)
{


for(i in 1:length(AFs)){

AF <- AFs[i]

print_text <- paste0("p-gain-threshold for ",test_number," tests, allele frequency ",AF)
print(print_text)
print(log(coef(fits[[i]])[1] + coef(fits[[i]])[2]*test_number,base=coef(fits[[i]])[3]))

}

}



p.gain.threshold(fits,test_number)

