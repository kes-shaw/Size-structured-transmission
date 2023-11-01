#Size structured transmission supplementary fig 1
#linear model params to show bad fits
#most of code slightly modified from R file "figure 1 tank level prevalence with model preds"

library(ggplot2)
library(tidyverse)
library(doBy)
library(viridis)
library(forcats)
library(gridExtra)
library(ggpubr)
#load dataset of the pop structures over a range of parasite levels
setwd("~/Desktop/Emory/Civitello Lab/My papers/Size structured transmission/figures")
dataset=as_tibble(read.csv("size class predictions_pop structures_csv.csv"))
##fully size dependent
prev_size = function(t=1, sigma_0, epsilon_0, sigma_sz, epsilon_sz, P0, N_S, N_M, N_L) {
  #make size dependent things, just hardcode in the average length of each size class (so 2.5 for small snails)
  epsilon_S = epsilon_0*exp(epsilon_sz*3.02)
  epsilon_M = epsilon_0*exp(epsilon_sz*7.57)
  epsilon_L =epsilon_0*exp(epsilon_sz*13.7)
  sigma_S = sigma_0*exp(sigma_sz*3.02)
  sigma_M = sigma_0*exp(sigma_sz*7.57)
  sigma_L =sigma_0*exp(sigma_sz*13.7)
  prev_S = 1- exp(-sigma_S*(epsilon_S/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_M = 1- exp(-sigma_M*(epsilon_M/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_L = 1- exp(-sigma_L*(epsilon_L/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_tank=((prev_S*N_S)+(prev_M*N_M)+(prev_L*N_L))/(N_S+N_M+N_L)
  prev_tank
}
# check it works
prev_size(t=1, sigma_0=0.7, epsilon_0=0.2, sigma_sz=-0.22, epsilon_sz=0.19, P0=250/18, N_S=12, N_M=3, N_L=3)
#now run it over the dataset with many levels of P0
#get coefficients (parameters) from the file "size strucutred transmission parameters_DJC_April 25"

coef(m_size_batch_linear)
#for sigma_0, take average of the different sigma estimates by batch
#sigma_0=mean(c( 1.7368561,  1.7423139,  1.8910352,  1.8945988,  0.7706458))
#sigma_0=1.60709
#epsilon_0=1.0000000
#sigma_sz=-0.2373604
#epsilon_sz= 0.0100000 
#vol=15
size_preds_lin = mapply(FUN=prev_size, P0 =dataset[,"P_Level"]/18, N_S=dataset[,"N_S"]/15, N_M=dataset[,"N_M"]/15, N_L=dataset[,"N_L"]/15,
                    MoreArgs = list(t=1, sigma_0 = 1.60709, epsilon_0 = 1.0000000, sigma_sz = -0.2373604, epsilon_sz =  0.0100000 ))


size_lin = cbind(dataset,size_preds_lin)
colnames(size_lin)=c("Structure", "P0", "N_S", "N_M", "N_L", "Prevalence")
#clean experimental data to get rid of infectivity control
prev_no_ctlr[,"Structure"] = factor(prev_no_ctlr[,"Structure"], levels= c("Uniform Small", "Small Skewed", "Equal", "Uniform Medium", "Uniform Large"))

###fully size dependent with linear stuff
size_plot_lin = ggplot(data=prev_no_ctlr, mapping=aes(x=P_Level*18, y=Prevalence.mean, colour = factor(Structure)))+
  geom_point(size = 3)+
  ylim(0,1)+
  geom_errorbar(data=prev_no_ctlr,mapping=aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(name="", type = c("#0375B4","#007849","black","#EB6E80", "#590494"))+
  scale_x_continuous(breaks=c(0,50,100,150,200,250))+
  labs(x = "Parasite density, L-1", y = "Prevalence ± SE", title = "D. Linear")+
  theme_classic(base_size=14)

q= size_plot_lin+geom_line(data=size_lin, aes(x=P0, y=Prevalence, color=Structure), lwd=1)
q

###now with exponential as the first panel
size_plot = ggplot(data=prev_no_ctlr, mapping=aes(x=P_Level*18, y=Prevalence.mean, colour = factor(Structure)))+
  geom_point(size = 3)+
  ylim(0,1)+
  geom_errorbar(data=prev_no_ctlr,mapping=aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(guide="none", type = c("#0375B4","#007849","black","#EB6E80", "#590494"))+
  scale_x_continuous(breaks=c(0,50,100,150,200,250))+
  labs(x = "", y = "Prevalence ± SE", title = "C. Exponential")+
  theme_classic(base_size=14)

p= size_plot+geom_line(data=size, aes(x=P0, y=Prevalence, color=Structure), lwd=1)
p

models = ggarrange(p,q, nrow=2, legend="bottom")

########
#ok panel A for showing the relationships
#postiive exponential = exposure
crv=as.tibble(exp(seq(0,5,.01)))

?curve
exp_size
coef(m_size_batch)
curve(expr=(.25*exp(.183*x)),from = 0, to = 25, n=500)
crv=crv%>%
  mutate(dog=seq(0,20,.04))

exp_exp=ggplot(crv, aes(x=dog, y=value/100))+
  geom_line(size=2)+
  ylim(0,1)+
  labs(x="", y="Exposure (ε)", title="A.Exponential")+
  theme_classic(base_size=14)

#neg = suscep
crvn=as.tibble(exp(-(seq(0,5,.01))))

crvn=crvn%>%
  mutate(dog=seq(0,20,.04))

sus_exp=ggplot(crvn, aes(x=dog, y=value))+
  geom_line(size=2)+
  labs(x="", y="Susceptibility (σ)", title=" ")+
  theme_classic(base_size=14)

#linear
lin=as.tibble(seq(0,5,.01))

lin=lin%>%
  mutate(dog=seq(0,20,.04))

lin_exp=ggplot(lin, aes(x=dog, y=value/5))+
  geom_line(size=2)+
  ylim(0,1)+
  labs(x="Diameter (mm)", y="Exposure (ε)", title="B. Linear")+
  theme_classic(base_size=14)

#linear susceptibility
lin_neg=as.tibble(10-(seq(5,10,.01)))

lin_neg=lin_neg%>%
  mutate(dog=seq(0,20,.04))

lin_sus=ggplot(lin_neg, aes(x=dog, y=(value/5)))+
  geom_line(size=2)+
  ylim(0,1)+
  labs(x="Diameter (mm)", y="Susceptibility(σ)", title="")+
  theme_classic(base_size=14)
left = ggarrange(exp_exp, sus_exp, lin_exp, lin_sus, heights=c(3,3,3,3))
ggarrange(left, models)
#####################################
#exposure only
prev_exp= function(t=1, sigma_0, epsilon_0, epsilon_sz, P0, N_S, N_M, N_L) {
  #make size dependent things, just hardcode in the average length of each size class (so 2.5 for small snails)
  epsilon_S = epsilon_0*exp(epsilon_sz*3.02)
  epsilon_M = epsilon_0*exp(epsilon_sz*7.57)
  epsilon_L =epsilon_0*exp(epsilon_sz*13.5)
  sigma_S = sigma_0
  sigma_M = sigma_0
  sigma_L =sigma_0
  prev_S = 1- exp(-sigma_S*(epsilon_S/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_M = 1- exp(-sigma_M*(epsilon_M/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_L = 1- exp(-sigma_L*(epsilon_L/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_tank=((prev_S*N_S)+(prev_M*N_M)+(prev_L*N_L))/(N_S+N_M+N_L)
  prev_tank
}

#now run it over the dataset with many levels of P0
#get coefficients (parameters) from the file "size strucutred transmission parameters"
coef(m_exp_only_batch_lin)
#sigma_0=mean(c(0.07552609, 0.08095818, 0.22989395, 0.23338486, 0.22011562 ))
#sigma_0=0.1679757
#epsilon_0= 0.19000000
#epsilon_sz=0.19000000
#vol=15
exp_preds_lin = mapply(FUN=prev_exp, P0 =dataset[,"P_Level"]/18, N_S=dataset[,"N_S"]/15, N_M=dataset[,"N_M"]/15, N_L=dataset[,"N_L"]/15,
                   MoreArgs = list(t=1, sigma_0 = 0.1679757, epsilon_0 =0.19000000 , epsilon_sz = 0.19000000))


exp_lin = cbind(dataset,exp_preds_lin)
colnames(exp_lin)=c("Structure", "N_S", "N_M", "N_L", "P0","Prevalence")

exp_plot_lin=ggplot(data=prev_no_ctlr, mapping=aes(x=P_Level*18, y=Prevalence.mean, colour = factor(Structure)))+
  geom_point(size=3)+
  ylim(0,1)+
  geom_errorbar(data=prev_no_ctlr,mapping=aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(guide="none", type = c("#952EA0","black","#2335be","#3CBB75FF", "#FDE725FF"))+
  scale_x_continuous(breaks=c(0,50,100,150,200,250))+
  labs(x = "", y = "", title = "B. Size-dependent exposure")+
  theme_classic(base_size=16)

r=exp_plot_lin+geom_line(data=exp_lin, aes(x=P0, y=Prevalence, color=Structure), lwd =1)

###############################
#suscep only
prev_sus = function(t=1, sigma_0, epsilon_0, sigma_sz, P0, N_S, N_M, N_L) {
  #make size dependent things, just hardcode in the average length of each size class (so 2.5 for small snails)
  epsilon_S = epsilon_0
  epsilon_M = epsilon_0
  epsilon_L =epsilon_0
  sigma_S = sigma_0*exp(sigma_sz*2.5)
  sigma_M = sigma_0*exp(sigma_sz*7)
  sigma_L =sigma_0*exp(sigma_sz*13.7)
  prev_S = 1- exp(-sigma_S*(epsilon_S/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_M = 1- exp(-sigma_M*(epsilon_M/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_L = 1- exp(-sigma_L*(epsilon_L/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_tank=((prev_S*N_S)+(prev_M*N_M)+(prev_L*N_L))/(N_S+N_M+N_L)
  prev_tank
}

coef(m_sus_only_batch_lin)
#sigma_0=mean(c(1.7366188,  1.7416310,  1.8894367,  1.8947735,  0.7710212))
#sigma_0=1.606696
#epsilon_0=0.0100000
#sigma_sz=-0.2372832
#vol=15
sus_preds_lin = mapply(FUN=prev_sus, P0 =dataset[,"P_Level"]/18, N_S=dataset[,"N_S"]/15, N_M=dataset[,"N_M"]/15, N_L=dataset[,"N_L"]/15,
                   MoreArgs = list(t=1, sigma_0 =1.606696, epsilon_0 = 0.0100000, sigma_sz = -0.2372832))
sus_lin = cbind(dataset,sus_preds_lin)
colnames(sus_lin)=c("Structure", "N_S", "N_M", "N_L", "P0","Prevalence")

sus_plot_lin=ggplot(data=prev_no_ctlr, mapping=aes(x=P_Level*18, y=Prevalence.mean, colour = factor(Structure)))+
  geom_point(size=3)+
  ylim(0,1)+
  geom_errorbar(data=prev_no_ctlr,mapping=aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(guide="none", type = c("#952EA0","black","#2335be","#3CBB75FF", "#FDE725FF"))+
  scale_x_continuous(breaks=c(0,50,100,150,200,250))+
  labs(x = "Parasite density, L-1", y = "Prevalence ± SE", title = "C. Size-dependent susceptibility")+
  theme_classic(base_size=16)

s=sus_plot_lin+geom_line(data=sus_lin, aes(x=P0, y=Prevalence, color=Structure), lwd=1)
