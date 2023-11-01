#Figure 2
#prevalence within size classes
#size-dependent model for main paper, the rest of them for the supplement
library(tidyverse)
library(doBy)
library(ggplot2)
##stuck on line 48, trying to subset my toy dataset to just the small snail structures
##then need to run through mapply with best prediction values and see if it works
#then think about adding in 95%CI

setwd("~/Desktop/Emory/Civitello Lab/Data/Size Structure 2021")
Data<-read.csv("combined prevalence data.csv")

####smalls first
Data = Data %>%
  mutate (Prevalence = NumInfec/(NumInfec+NumUninfec))
head(Data)
prev_size = summaryBy(Prevalence ~ Structure+P_Level+Size, data=Data, FUN=c(length,mean,sd))
prev_small = prev_size[prev_size$Size == "Small" & prev_size$Structure != "Control Medium",]
head(prev_small)
##plot smalls just to look
ggplot(prev_small, aes(P_Level, Prevalence.mean, colour = factor(Structure)))+
  geom_point(size = 3)+
  ylim(0,1)+
  geom_errorbar(aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(name = "Size Structure", type = c("black", "#f76f73", "#027fdc"))+
  labs(title = "A", y = " ", x = " ")+
  theme_classic(base_size=12)

##ok now chug through the model to get model predictions
prev_S = function(t=1, sigma_0, epsilon_0, sigma_sz, epsilon_sz, P0, N_S, N_M, N_L) {
  #make size dependent things, just hardcode in the average length of each size class (so 2.5 for small snails)
  epsilon_S = epsilon_0*exp(epsilon_sz*2.5)
  epsilon_M = epsilon_0*exp(epsilon_sz*7)
  epsilon_L =epsilon_0*exp(epsilon_sz*13.5)
  sigma_S = sigma_0*exp(sigma_sz*2.5)
  sigma_M = sigma_0*exp(sigma_sz*7)
  sigma_L =sigma_0*exp(sigma_sz*13.5)
  prev_S = 1- exp(-sigma_S*(epsilon_S/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_S
}
#check to make sure it works
prev_S(t=1, sigma_0=.2, epsilon_0=1, sigma_sz=.05, epsilon_sz = 0.001, P0=144, N_S=6, N_M=6, N_L=6, D_S=6, D_M=6, D_L=6, I_S=1, I_M=2, I_L=3)
#load dataset of the pop structures over a range of parasite levels
setwd("~/Desktop/Emory/Civitello Lab/My papers/Size structured transmission/figures")
dataset=as_tibble(read.csv("size class predictions_pop structures_csv.csv"))

#ok now just look at ones with small snails

dataset = dataset%>%
  filter(Structure!="Uniform Medium")%>%
  filter(Structure!="Uniform Large")

#first just use the best param estimates from the model
coef(m_size_batch)
sigma_0=0.6085019
epsilon_0=0.2498385
sigma_sz=-0.2139004
epsilon_sz= 0.1835933
vol=15
head(dataset)
#prev_S(t=1, sigma_0=.2, epsilon_0=0.5, sigma_sz=.05, epsilon_sz = 0.001, P0=144, N_S=6/15, N_M=6/15, N_L=6/15)

small_preds = mapply(FUN=prev_S, P0 =dataset[,"P_Level"]/18, N_S=dataset[,"N_S"]/15, N_M=dataset[,"N_M"]/15, N_L=dataset[,"N_L"]/15,
                     MoreArgs = list(t=1, sigma_0 = 0.6895559, epsilon_0 = 0.1904044, sigma_sz = -0.2282718, epsilon_sz = 0.1950911))

small_predictions = cbind(dataset,small_preds)
colnames(small_predictions)=c("Structure", "N_S", "N_M", "N_L", "P_Level", "Prevalence")
##Med
dataset=as_tibble(read.csv("size class predictions_pop structures_csv.csv"))
dataset = dataset%>%
  filter(Structure!="Uniform Small")%>%
  filter(Structure!="Uniform Large")

prev_M = function(t=1, sigma_0, epsilon_0, sigma_sz, epsilon_sz, P0, N_S, N_M, N_L) {
  #make size dependent things, just hardcode in the average length of each size class (so 2.5 for small snails)
  epsilon_S = epsilon_0*exp(epsilon_sz*2.5)
  epsilon_M = epsilon_0*exp(epsilon_sz*7)
  epsilon_L =epsilon_0*exp(epsilon_sz*13.5)
  sigma_S = sigma_0*exp(sigma_sz*2.5)
  sigma_M = sigma_0*exp(sigma_sz*7)
  sigma_L =sigma_0*exp(sigma_sz*13.5)
  prev_M = 1- exp(-sigma_M*(epsilon_M/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_M
}

med_preds = mapply(FUN=prev_M, P0 =dataset[,"P_Level"]/18, N_S=dataset[,"N_S"]/15, N_M=dataset[,"N_M"]/15, N_L=dataset[,"N_L"]/15,
                     MoreArgs = list(t=1, sigma_0 = 0.6895559, epsilon_0 = 0.1904044, sigma_sz = -0.2282718, epsilon_sz = 0.1950911))


med_predictions = cbind(dataset, med_preds)
colnames(med_predictions)=c("Structure", "N_S", "N_M", "N_L", "P_Level", "Prevalence")
#Large
dataset=as_tibble(read.csv("size class predictions_pop structures_csv.csv"))
dataset = dataset%>%
  filter(Structure!="Uniform Small")%>%
  filter(Structure!="Uniform Medium")

prev_L = function(t=1, sigma_0, epsilon_0, sigma_sz, epsilon_sz, P0, N_S, N_M, N_L) {
  #make size dependent things, just hardcode in the average length of each size class (so 2.5 for small snails)
  epsilon_S = epsilon_0*exp(epsilon_sz*2.5)
  epsilon_M = epsilon_0*exp(epsilon_sz*7)
  epsilon_L =epsilon_0*exp(epsilon_sz*13.5)
  sigma_S = sigma_0*exp(sigma_sz*2.5)
  sigma_M = sigma_0*exp(sigma_sz*7)
  sigma_L =sigma_0*exp(sigma_sz*13.5)
  prev_L = 1- exp(-sigma_L*(epsilon_L/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_L
}

lg_preds = mapply(FUN=prev_L, P0 =dataset[,"P_Level"]/18, N_S=dataset[,"N_S"]/15, N_M=dataset[,"N_M"]/15, N_L=dataset[,"N_L"]/15,
                   MoreArgs = list(t=1, sigma_0 = 0.6895559, epsilon_0 = 0.1904044, sigma_sz = -0.2282718, epsilon_sz = 0.1950911))


lg_predictions = cbind(dataset, lg_preds)
colnames(lg_predictions)=c("Structure", "N_S", "N_M", "N_L", "P_Level", "Prevalence")

##ok so now plot the predictions with the expeirmental data

#experimental data from here: 
setwd("~/Desktop/Emory/Civitello Lab/Data/Size Structure 2021")
Data<-read.csv("combined prevalence data.csv")
Data = Data %>%
  mutate (Prevalence = NumInfec/(NumInfec+NumUninfec))
prev_size = summaryBy(Prevalence ~ Structure+P_Level+Size, data=Data, FUN=c(length,mean,sd))
prev_small = prev_size[prev_size$Size == "Small",]
prev_no_ctlr = subset(prev_size, Structure != "Control Medium")
prev_med = prev_no_ctlr[prev_no_ctlr$Size == "Medium",]
prev_lg = prev_size[1:3 & prev_size$Size == "Large",]

a=ggplot(data=small_predictions, aes(x=P_Level, y=Prevalence, color=Structure))+
  geom_point(data=small_predictions)+
  #ylim(0,25)+
  labs(x = "Parasite Level", y = "Prevalence")+
  scale_x_continuous(breaks=c(0,50,100,150,200,250))+
  #geom_line(y=0.1686205, lty=6)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 13), breaks = pretty_breaks())+
  #geom_ribbon(aes(ymin = 0.1475, 
  # ymax=0.1918), colour="light grey", fill="light grey",alpha=0.5)+
  theme_classic(base_size=16)

sm=ggplot(data=prev_small, mapping=aes(x=P_Level*18, y=Prevalence.mean, colour = factor(Structure)))+
  geom_point(size=3)+
  ylim(0,1)+
  geom_errorbar(data=prev_small,mapping=aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(name = "Size Structure", type = c("black","black", "#c51b8a"))+
  labs(x = "", y = "Prevalence", title = "A.")+
  scale_x_continuous(breaks=c(0,50,100,150,200,250))+
  theme_classic(base_size=16)

a=sm+geom_line(data=small_predictions, aes(x=P_Level, y=Prevalence, color=Structure))

###plot med
med=ggplot(data=prev_med, mapping=aes(x=P_Level*18, y=Prevalence.mean, colour = factor(Structure)))+
  geom_point(size = 3)+
  ylim(0,1)+
  geom_errorbar(data=prev_med,mapping=aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(name = "", type = c("#fb6a4a", "#de2d26", "#a50f15"))+
  labs(x = "", y = "Prevalence", title = "B.")+
  scale_x_continuous(breaks=c(0,50,100,150,200,250))+
  theme_classic(base_size=16)

b=med+geom_line(data=med_predictions, aes(x=P_Level, y=Prevalence, color=Structure))


###plot lg
lg=ggplot(prev_lg, mapping=aes(x=P_Level*18, y=Prevalence.mean, colour = factor(Structure)))+
  geom_point(size = 3)+
  ylim(0,1)+
  geom_errorbar(data=prev_lg,mapping=aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(name = "", type = c("#9e9ac8", "#8856a7", "#810f7c"))+
  labs(x = "Parasite Level", y = "Prevalence", title = "C.")+
  scale_x_continuous(breaks=c(0,50,100,150,200,250))+
  theme_classic(base_size=16)

c=lg+geom_line(data=lg_predictions, aes(x=P_Level, y=Prevalence))
library(ggpubr)
ggarrange(a,b,c, nrow=3, widths=c(.1,.1,.1), heights = c(3,3,3))

###replotting large to just have experimental data
##at lowest parasite level 2 points overlap; jittering options all ugly
##for EEID talk going to change parasite level in the dataframe to be slightly higher in order
##to achieve an x jitter for this one point
##for equal, parasites = 2 per snail going to make it 38 parasites total (instead of the actual 36)
setwd("~/Desktop/Emory/Civitello Lab/Data/Size Structure 2021/analysis")
Data_jitter = read.csv("lg snail plot prevalence data.csv") 

Data_jitter = Data_jitter %>%
  mutate (Prevalence = NumInfec/(NumInfec+NumUninfec))
prev_size_jitter = summaryBy(Prevalence ~ Structure+P_Level+Size, data=Data_jitter, FUN=c(length,mean,sd))
prev_lg_jitter = prev_size_jitter[prev_size_jitter$Size == "Large",]

lg_jitter_plot=ggplot(prev_lg_jitter, mapping=aes(x=P_Level*18, y=Prevalence.mean, colour = factor(Structure)))+
  geom_point(size = 3)+
  ylim(0,1)+
  geom_errorbar(data=prev_lg_jitter,mapping=aes(ymin=Prevalence.mean-(Prevalence.sd/4), ymax=Prevalence.mean+(Prevalence.sd/4)), width=.05,
                position=position_dodge(.9), colour = "black")+
  scale_colour_discrete(name = "Size Structure", type = c("#9e9ac8", "#8856a7", "#810f7c"))+
  labs(x = "Parasite Level", y = "Prevalence", title = "")+
  theme_classic(base_size=16)
