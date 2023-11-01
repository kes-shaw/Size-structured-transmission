#Size structured transmission figure 3
#susceptibility, exposure, and transmission rate changing with body size
library (tidyverse)
library(doBy)
library(ggplot2)
library(MASS) #for the covariance matrix to find the good pars
library(GGally) #to visualize how the good pars vary together
#library(usethis) 
#usethis::edit_r_environ()

#panel A: exposure
#Need to pull from the 95%CI of good pars, similar to field data analysis
# We can approximate the 95% CI on parameter SETs by drawing the parameters from the multivariate normal, with the means specified by the parameter estimates and the variance-covariance matrix approximated by the information matrix

# You can use the fitted information matrix from the model to approximate the var-cov matrix for pars
#here "m_size" is from the R script "size strucutred transmission parameters"
par_sets = data.frame(mvrnorm(10000, mu=coef(m_size_batch), Sigma = vcov(m_size_batch)))
##need to change to batched one!!!!!
#you can visualize how the different parameter estimates co-vary
pairs(par_sets)

#now calculate the neg log likelihoods of these parameter sets that are drawn from the covariance matrix
NLL_for_par_sets = mapply(NLL_experiment_size_batches,
                          sigma_0_A = par_sets$sigma_0_A, sigma_0_B = par_sets$sigma_0_B, sigma_0_C = par_sets$sigma_0_C, sigma_0_D = par_sets$sigma_0_D, sigma_0_E = par_sets$sigma_0_E,epsilon_0 = par_sets$epsilon_0,
                          epsilon_sz=par_sets$epsilon_sz, sigma_sz=par_sets$sigma_sz)


par_sets[,"NLL"] = NLL_for_par_sets

# Now This cuts these distributions to the 95% for parameter sets and lets  you visualize it
ggpairs(subset(par_sets, NLL <= min(par_sets$NLL + 2))) #+2 is the 95% CI for the NLLs; only positive because an NLL is maximally negative, duh
#ok save the best pars sets!
good_pars = subset(par_sets, NLL <= min(par_sets$NLL + 2))
head(good_pars)
nrow(good_pars)
##decided not to do this and limit them down the line
#now take the min, mean, and max of those epsilons and use mapply to loop over many
#best_pars = rbind(apply(good_pars,2,mean),apply(good_pars,2,min),apply(good_pars,2,max))
#rownames(best_pars) <- c("pars_mean", "pars_min", "pars_max")

#now write a function to apply this epsilon across a ton of body sizes
eps_L = function(eps_0, eps_sz, Length){
  eps_L = eps_0*exp(eps_sz*Length)
  eps_L
}


#ok now collect all the epsilon_0 and epsilon_size values from the good parameters
#eps_0 = good_pars[,2]
#eps_sz = good_pars[,4]
#ok put them in the same dataframe
#eps_groups = as_data_frame(cbind(eps_0,eps_sz))
#colnames(eps_groups)=c("eps_0", "eps_sz")
#head(eps_groups)
# ^^^^^^^
#decided not to do the above because will just use the good_pars for epsilon and sigma and call directly from that

#replicate all the good pars as many times as there are snail lengths you want
pars_repeat = bind_rows(replicate(191, good_pars, simplify = FALSE))
nrow(pars_repeat)

#make the sizes you want
Sizes = seq(from=1, to=20, by=.1)
#repeat them as many times as there are parameter sets
#will need to change  this each time you run good_pars because the number that fall within in the 
#95%CI will vary with each draw!!
Length=rep(Sizes,each=1103980)
length(Length)
#bring them all together
pars_w_L=cbind(pars_repeat,Length)
#head(good_pars)
colnames(pars_w_L) = c("sig_0","eps_0", "sig_sz","eps_sz", "NLL","Length")

#in mapply, list the function for it to use and then each argument in the function is a column name of 
#the dataframe it should scroll over. MoreArgs are any arguements  that don't change in each iteration and
#must be a list
eps_predictions = mapply(FUN=eps_L, eps_0 = pars_w_L$eps_0, eps_sz = pars_w_L$eps_sz, Length = pars_w_L$Length)
#length(predictions) == dim(pars_w_Length)[1]#should be the same length and return "true"

#put them back together in the same dataframe
results = cbind(pars_w_L, eps_predictions)
tail(results)
#get the mean, min, and max which represent the 95%CI
epsilon_preds_clean = summaryBy(eps_predictions~Length, data= results, FUN=c(mean,min,max))
ggplot(epsilon_preds_clean, aes(x=Length, y=eps_predictions.mean))+
  geom_point()+
  ylim(0,25)+
  geom_errorbar(aes(ymin=eps_predictions.min, ymax=eps_predictions.max), width=.05,
                position=position_dodge(.9), colour = "black")+
  labs(x = "Length", y = "Epsilon")+
  #geom_line(y=0.1686205, lty=6)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 13), breaks = pretty_breaks())+
  #geom_ribbon(aes(ymin = 0.1475, 
                 # ymax=0.1918), colour="light grey", fill="light grey",alpha=0.5)+
  theme_classic(base_size=16)

ggplot(epsilon_preds_clean, aes(x=Length, y=eps_predictions.mean))+
  geom_line(color = "black")+
  ylim(0,25)+
  labs(x = "Length", y = "Epsilon")+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 13), breaks = pretty_breaks())+
  geom_ribbon(aes(ymin = eps_predictions.min, 
   ymax = eps_predictions.max), colour="light grey", fill="light grey",alpha=0.5)+
  theme_classic(base_size=16)

##ok try the same thing but with susceptibility
suscep_L = function(suscep_0, suscep_sz, Length){
  suscep_L = suscep_0*exp(suscep_sz*Length)
  suscep_L
}

#we already did all the needed legwork above to get pars_w_L so we  don't need to repeat that
head(pars_w_L)
suscep_predictions = mapply(FUN=suscep_L, suscep_0 = pars_w_L$sig_0, suscep_sz = pars_w_L$sig_sz, Length = pars_w_L$Length)
#length(predictions) == dim(pars_w_Length)[1]#should be the same length and return "true"

#put them back together in the same dataframe
results = cbind(pars_w_L, suscep_predictions)
tail(results)
#get the mean, min, and max which represent the 95%CI
suscep_preds_clean = summaryBy(suscep_predictions~Length, data= results, FUN=c(mean,min,max))
head(suscep_preds_clean)
ggplot(suscep_preds_clean, aes(x=Length, y=suscep_predictions.mean))+
  geom_line()+
  #ylim(0,25)+
  #geom_errorbar(aes(ymin=suscep_predictions.min, ymax=suscep_predictions.max), width=.05,
                #position=position_dodge(.9), colour = "black")+
  labs(x = "Length", y = "Sigma")+
  #geom_line(y=0.1686205, lty=6)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 13), breaks = pretty_breaks())+
  geom_ribbon(aes(ymin = suscep_predictions.min, 
  ymax=suscep_predictions.max), colour="light grey", fill="blue",alpha=0.5)+
  theme_classic(base_size=16)


###add in ganesh estimates to sigma graph
#setwd("~/Desktop/Emory/Civitello Lab/Data/Size Structure 2021/analysis")
#size_grad = read.csv("ganesh_snail_size.csv")
#get rid of false precision of there being so many sig figs to the size measurements (an artifact of imageJ)
size_grad[,"Size"]=round(size_grad[,"Size"], digits = 1)
head(size_grad)
s1<-subset(size_grad, Size<=2) 
#s1
m1<-mle2(Infected ~ dbinom(prob = 1-exp(-sigma_ganesh* s1[,"Miracidia"]),size=1),data=s1,start=list(sigma_ganesh=0.6), 
         control=list(parscale=c(sigma_ganesh=0.6)))
m1

s2<-subset(size_grad,Size>2 & Size<=4)
#s2
m2<-mle2(Infected ~ dbinom(prob = 1-exp(-sigma_ganesh* s2[,"Miracidia"]),size=1),data=s2,start=list(sigma_ganesh=0.6), 
         control=list(parscale=c(sigma_ganesh=0.6)))
m2

s3<-subset(size_grad,Size>4 & Size<=6)
#s3
m3<-mle2(Infected ~ dbinom(prob = 1-exp(-sigma_ganesh* s3[,"Miracidia"]),size=1),data=s3,start=list(sigma_ganesh=0.6), 
         control=list(parscale=c(sigma_ganesh=0.6)))
m3

s4<-subset(size_grad,Size>6 & Size<=8)
#s4
m4<-mle2(Infected ~ dbinom(prob = 1-exp(-sigma_ganesh* s4[,"Miracidia"]),size=1),data=s4,start=list(sigma_ganesh=0.6), 
         control=list(parscale=c(sigma_ganesh=0.6)))
m4

s5<-subset(size_grad,Size>8 & Size<=10)
#s5
m5<-mle2(Infected ~ dbinom(prob = 1-exp(-sigma_ganesh* s5[,"Miracidia"]),size=1),data=s5,start=list(sigma_ganesh=0.6), 
         control=list(parscale=c(sigma_ganesh=0.6)))
m5

s6<-subset(size_grad,Size>10 & Size<=12)
#s6
m6<-mle2(Infected ~ dbinom(prob = 1-exp(-sigma_ganesh* s6[,"Miracidia"]),size=1),data=s6,start=list(sigma_ganesh=0.6), 
         control=list(parscale=c(sigma_ganesh=0.6)))
m6

s7<-subset(size_grad,Size>12 & Size<=14)
#s7
m7<-mle2(Infected ~ dbinom(prob = 1-exp(-sigma_ganesh* s7[,"Miracidia"]),size=1),data=s7,start=list(sigma_ganesh=0.6), 
         control=list(parscale=c(sigma_ganesh=0.6)))
m7

ganesh_sigs = c(coef(m1), coef(m2), coef(m3), coef(m4), coef(m5), coef(m6), coef(m7))
size_grps = c(1.5,3,5,7,9,11,13)
ganesh_sig_est=as_tibble(cbind(ganesh_sigs, size_grps))
colnames(ganesh_sig_est)=c("Sigma", "Size")

###add ganeshs to the plot
a = ggplot(suscep_preds_clean, aes(x=Length, y=suscep_predictions.mean))+
  geom_line()+
  #ylim(0,25)+
  #geom_errorbar(aes(ymin=suscep_predictions.min, ymax=suscep_predictions.max), width=.05,
  #position=position_dodge(.9), colour = "black")+
  labs(x = "Length", y = "Sigma")+
  #geom_line(y=0.1686205, lty=6)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 13), breaks = pretty_breaks())+
  geom_ribbon(aes(ymin = suscep_predictions.min, 
                  ymax=suscep_predictions.max), colour="light grey", fill="light grey",alpha=0.5)+
  theme_classic(base_size=16)
a+  geom_point(data=ganesh_sig_est, mapping=aes(x=size_grps, y=ganesh_sigs), color = "red")

#####ok now combine them for beta
beta_L = function(eps_0, eps_sz,suscep_0, suscep_sz, Length){
  suscep_L = suscep_0*exp(suscep_sz*Length)
  eps_L = eps_0*exp(eps_sz*Length)
  beta = eps_L*suscep_L
  beta
}
head(pars_w_L)
beta_predictions = mapply(FUN=beta_L, eps_0=pars_w_L$eps_0, eps_sz=pars_w_L$eps_sz, suscep_0 = pars_w_L$sig_0, suscep_sz = pars_w_L$sig_sz, Length = pars_w_L$Length)

#put them back together in the same dataframe
results = cbind(pars_w_L, beta_predictions)
tail(results)
#get the mean, min, and max which represent the 95%CI
beta_preds_clean = summaryBy(beta_predictions~Length, data= results, FUN=c(mean,min,max))
head(beta_preds_clean)
ggplot(beta_preds_clean, aes(x=Length, y=beta_predictions.mean))+
  geom_line()+
  #ylim(0,25)+
  #geom_errorbar(aes(ymin=beta_predictions.min, ymax=beta_predictions.max), width=.05,
               # position=position_dodge(.9), colour = "black")+
  labs(x = "Length", y = "beta")+
  #geom_line(y=0.1686205, lty=6)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 13), breaks = pretty_breaks())+
  geom_ribbon(aes(ymin = beta_predictions.min, 
   ymax=beta_predictions.max), colour="light grey", fill="#952EA0",alpha=0.5)+
  theme_classic(base_size=16)

