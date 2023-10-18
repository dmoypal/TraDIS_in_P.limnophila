#Install the necessary packages

#install.packages("MASS")
#install.packages("fitdistrplus")
#install.packages("survival")
#install.packages("dplyr")

library(survival)
library(MASS)
library(fitdistrplus)
library(dplyr)
library(stringr)

dir="/Users/davidmoyanopalazuelo/Library/CloudStorage/Dropbox/proyectos/prueba_que_los_scripts_de_github_funcionan_bien/"
setwd(dir)

#load the data
data<-read.delim(paste0(dir,"with_trans.tnseq.tsv"),header=TRUE,sep="\t")

#get only gene and pseudogene
filtrado <- filter(data, type %in% c("gene", "pseudogene")) %>% 
  mutate(Name = str_match(info, "ID=([^;]*)")[,2]) %>%
  mutate(all_div_length = count_all / (end - start)) %>%
  dplyr::select(Name,all_div_length)

filtrado$Name <- str_replace(filtrado$Name, "gene-", "")

all_plim <- all(str_detect(filtrado$Name, "Plim"))#check all have the same locus_tag

insertion_index <- filtrado$all_div_length
histrogram <- hist(insertion_index, breaks ="FD",prob=TRUE, main="Distribution")
nclass.FD(insertion_index) ## calculamos el bind width en base al numero de agrupaciones que veremos 


#Define two cut-offs
cutoff1 <- 0.013
cutoff2 <- 0.022

fit_exp <- fitdistr(insertion_index[insertion_index < cutoff1],densfun = "exponential")
fit_gamma <- fitdistr(insertion_index[insertion_index > cutoff2 ],densfun = "gamma")
hist(insertion_index, pch=20, breaks=36, prob=TRUE)
abline(v=cutoff1, col="green", lwd=2)
abline(v=cutoff2, col="orange", lwd=2)

#draw the distribution and fits
curve(dgamma(x, shape = fit_gamma$estimate[1], rate = fit_gamma$estimate[2]), col="blue", lwd=2, add=T)
curve(dexp(x, rate = fit_exp$estimate[1]) ,col="red", lwd=2, add=T)
legend("topright", legend=c("Gamma", "Exponential", "Cut-off 1", "Cut-off 2"), col=c("blue", "red", "green", "orange"), lwd=2)

prob_exp <- dexp(insertion_index, rate = fit_exp$estimate) ## Probabilities of the genes that are under exponential distribution
prob_gamma <- dgamma(insertion_index, shape = fit_gamma$estimate[1], rate = fit_gamma$estimate[2]) # the same but with gamma dist
log_ratio <- log2(prob_exp/prob_gamma)

essential_exp_gamma<- sum(log_ratio > log2(12)) 
essential_exp_gamma_array <- log_ratio > log2(12)
non_essential_exp_gamma <- sum(log_ratio < -log2(12)) 
non_essential_exp_gamma_array <- log_ratio < -log2(12)
unclear_exp_gamma <- sum(log_ratio > -log2(12) & log_ratio < log2(12)) 
unclear_exp_gamma_array <- log_ratio > -log2(12) & log_ratio < log2(12)
  
#add essentiality label. 
filtrado[,"clase"] <- NA
for(i in seq(1,length(essential_exp_gamma_array))){
  if (essential_exp_gamma_array[i]==TRUE){
    filtrado[ i, "clase"] <- "essential"
  }
}
for(i in seq(1,length(non_essential_exp_gamma_array))){
  if (non_essential_exp_gamma_array[i]==TRUE){
    filtrado[ i, "clase"] <- "non-essential"
  }
}
for(i in seq(1,length(unclear_exp_gamma_array))){
  if (unclear_exp_gamma_array[i]==TRUE){
    filtrado[ i, "clase"] <- "unclear"
  }
}
  
#write the essentiality file
write.table(filtrado, paste0(dir,"essentiality.tsv"), sep="\t", row.names=FALSE, quote=FALSE)