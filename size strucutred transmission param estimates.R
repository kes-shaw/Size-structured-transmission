#Size structured param estimation
#Feb 2022
library(tidyverse)
library(bbmle)
library(scales)

setwd("~/Desktop/Emory/Civitello Lab/Data/Size Structure 2021/analysis")
dataset = read.csv("Cleaned Prevalence Data V2.csv")
controls = read.csv("Controls.csv")
size_grad = read.csv("ganesh_snail_size.csv")


####Current functionality###
# we can generate prevalence estimates for the experiment (name the funciton here)
# we can calculate neg log likelihoods for becky's experiment (name the function here)
#we can incorporate the data from becky's experiment, becky's controls, ganesh's susceptibiltiy experiment
#fits size dependent models and compares their AIC scores

###Pending tasks###
#deal with Ganesh data having batches, figure out what is true. 
#way later: explore different more complicate sigma/epsilon functions?

#See kelsey/becky lab notebook pg. 74 for the prevalence equation written out

#start with baseline sigma and epsilon and change later
#t is always 1
#going to write a transmission function then use mapply
#name all the different arguements in the function
#N_i is the number of snails of Small, Med, Large class at exposure. D=diagnosed; I= infected

#make prevalence function and make sure the prevalences make sense
prev_all = function(t=1, sigma_0, epsilon_0, P0, N_S, N_M, N_L, D_S, D_M, D_L, I_S, I_M, I_L) {
  prev_S = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  prev_M = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  prev_L = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  c(prev_S, prev_M, prev_L)
}

prev_all(sigma_0 = .2, epsilon_0 = 10, P0 = 8, N_S = 6, N_M = 6, N_L = 6, D_S = 0, D_M = 0, D_L = 0)
#if it worked these should all be zero

##ok now mess around with it more
#these prevlaneces are model estimates so depend on the experimental conditions and the parameters we give. not related to the data we have yet
#aka #predictions
#so now we need likelihood scores to compare predictions to the data so we can choose good parameters!
#Neg log likes time babyyyyyy
#3 predictions for each size class, even in tanks where not all sizes are present. It will always be scored with the prediction versus the data
    #when there's not data then you get a perfect score

NLL_S = -(dbinom(x = I_S, prob = prev_S, size=D_S, log = TRUE))
NLL_M = -(dbinom(x = I_M, prob = prev_M, size=D_M, log = TRUE))
NLL_L = -(dbinom(x = I_L, prob = prev_L, size=D_L, log = TRUE))
#?dbinom; x is the observation; prob is the probability based on our model; size is the number diagnosed



##now take all these NLL's and put them in your original function(won't work outside of it)
?
NLL_rep <- function(t=1, sigma_0, epsilon_0, P0, N_S, N_M, N_L, D_S, D_M, D_L, I_S, I_M, I_L){
  sigma_0 = 1/(1+exp(-sigma_0)) #logit transform makes any number constrained to between zero and 1; important for later mle function
  epsilon_0 = exp(epsilon_0) #constrains it to be positive
  prev_S = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  prev_M = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  prev_L = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  c(prev_S, prev_M, prev_L)
  NLL_S = -(dbinom(x = I_S, prob = prev_S, size=D_S, log = TRUE))
  NLL_M = -(dbinom(x = I_M, prob = prev_M, size=D_M, log = TRUE))
  NLL_L = -(dbinom(x = I_L, prob = prev_L, size=D_L, log = TRUE))
  #Total NLL is what get returned
  NLL_S+NLL_M+NLL_L
}

##above doesnt work for some weird ass reason so below is from dave but lacks comments

NLL_rep = function(t=1, sigma_0, epsilon_0, P0, N_S, N_M, N_L, D_S, D_M, D_L, I_S, I_M, I_L) {
  
  #transformations
  
  #sigma_0 = 1/(1+exp(-sigma_0))
  
  #epsilon_0 = exp(epsilon_0)
  
  
  
  prev_S = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  
  prev_M = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  
  prev_L = 1- exp(-sigma_0*(epsilon_0/(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))*P0*(1-exp(-t*(epsilon_0*N_S+epsilon_0*N_M+epsilon_0*N_L))))
  
  NLL_S = -(dbinom(x=I_S, prob = prev_S, size = D_S, log = T))
  
  NLL_M = -(dbinom(x=I_M, prob = prev_M, size = D_M, log = T))
  
  NLL_L = -(dbinom(x=I_L, prob = prev_L, size = D_L, log = T))
  
  NLL_S + NLL_M + NLL_L
}

#so the above looks at one container/replicate and based on how you set it up it'll look at the 3 size classes and score based on how well the data matches the predictions
NLL_rep(t=1, sigma_0=.2, epsilon_0=1, P0=144, N_S=6, N_M=6, N_L=6, D_S=6, D_M=6, D_L=6, I_S=1, I_M=2, I_L=3)

#so to do with all the replicates, use the magic of mapply. apply function similar to for loop but faster
#make new function using old function 
#give mapply the name of the function you want to cal; the rest of the arguments are anything that's a column in your dataset that is part of the function you call
#then you do moreargs and you give it a list of the other arguments
NLL_experiment= function(t=1, sigma_0, epsilon_0, df, control_df){
  NLLs = mapply(NLL_rep, P0=df[,"P_Level"]*18/vol, N_S=df[,"Start_Small"], N_M=df[,"Start_Med"], 
                N_L=df[,"Start_Large"], D_S=df[,"D_Small"], D_M=df[,"D_Med"],
                D_L = df[,"D_Large"], I_S = df[,"Infec_Small"], I_M = df[,"Infec_Med"], I_L = df[,"Infec_Large"],
                MoreArgs = list(t=t, sigma_0 = sigma_0, epsilon_0 = epsilon_0))
sum(NLLs)
  }
#dataset = read.csv("Cleaned Prevalence Data V2.csv")
NLL_experiment(sigma_0 = .2, epsilon_0 = 1, df = df, control_df=control_df)

###real nll expeirment
#right now epsilon is snails per bucket; want it as snails per L so divide by 15
NLL_experiment= function(t=1, sigma_0, epsilon_0, df=dataset, control = controls, vol = 15){
  NLLs = mapply(NLL_rep, P0=dataset[,"P_Level"]*18/vol, N_S=dataset[,"Start_Small"]/vol, N_M=dataset[,"Start_Med"]/vol, 
                N_L=dataset[,"Start_Large"]/vol, D_S=dataset[,"D_Small"], D_M=dataset[,"D_Med"],
                D_L = dataset[,"D_Large"], I_S = dataset[,"Infec_Small"], I_M = dataset[,"Infec_Med"], I_L = dataset[,"Infec_Large"],
                MoreArgs = list(t=t, sigma_0 = sigma_0, epsilon_0 = epsilon_0))
  
  #now evaluate the NLL for Becky control data to account for different infectivity of batches
  #now predict the sigma for the individuals of the body size used in the controls; will change from sigma_0 later, this is placeholder
  sigma_control = sigma_0
  #can use the density function for the binomial distribution to get the log likelihood score. We want to know the likelihood of each observation (each snail a replicate here)
  control_NLLs = -sum(dbinom(x=control[,"Infected"], prob = 1-exp(-sigma_0* control[,"P_Level"]), size = 1, log = T))
  #change it so that if a bad prediction is made, give an NA as a score; if the prediction isn't nonsense, sum the scores
  #function "any" is very useful
  if(any(!is.finite(NLLs)|!is.finite(control_NLLs))){return(NA)}else{sum(NLLs, control_NLLs)}
  sum(NLLs)
}
NLL_experiment(t=1, sigma_0 = .2, epsilon_0 = 1, df = df, control=controls)
#sigma is the chance a single miricidium infects a snail; epsilon is the amount of water in which a miricidium can find every snail in a day

#returns the score for those two parameters, given this model
#so now given this (size-ignoring, not good) model, what are the best param values?
#make sure bbmle is loaded up top
#we wrote our neg log likelhiood function so use that; start is the place to start looking
m_null = mle2(minuslogl = NLL_experiment, start = list(epsilon_0 = .1, sigma_0 = .2))
coef(m_null)
#sigma is the chance a single miricidium infects a snail; epsilon is the amount of water in which a miricidium can find every snail in a day
#can profile the output to look at the confidence intervals and the shape of the curve of the scores of the params it tried
#p_null= profile(m_null)
#confint(p_null)
#plot(p_null)

#need the logit transformation for sigma and exponentiate epsilon so epsilon is always positive and sigma is between zero and 1 in order to get reasonable predictions


##ok now need control data to pin down sigma to better estimate epsilon
controls = read.csv("Controls.csv")
head(controls)
#get rid of weird NAs
controls = controls[, sapply(controls, function(i) !all(is.na(i)))]
#weird columns from choosing the wrong type of csv to save it as (not the UHT 8 thing)
head(controls)
m_controls_null=mle2(Infected~dbinom(prob = 1-exp(-sigma*P_Level), size = 1), data = controls, start = list(sigma = .1))

##Dave code: there is an impact of batch on sigma so we should incorporate that going forward

#add in controls
NLL_experiment= function(t=1, sigma_0, epsilon_0, df=dataset, control = controls, vol = 15){
  NLLs = mapply(NLL_rep, P0=dataset[,"P_Level"]*18/vol, N_S=dataset[,"Start_Small"]/vol, N_M=dataset[,"Start_Med"]/vol, 
                N_L=dataset[,"Start_Large"]/vol, D_S=dataset[,"D_Small"], D_M=dataset[,"D_Med"],
                D_L = dataset[,"D_Large"], I_S = dataset[,"Infec_Small"], I_M = dataset[,"Infec_Med"], I_L = dataset[,"Infec_Large"],
                MoreArgs = list(t=t, sigma_0 = sigma_0, epsilon_0 = epsilon_0))
  
  #now evaluate the NLL for Becky control data to account for different infectivity of batches
  #now predict the sigma for the individuals of the body size used in the controls; will change from sigma_0 later, this is placeholder
  sigma_control = sigma_0
  #can use the density function for the binomial distribution to get the log likelihood score. We want to know the likelihood of each observation (each snail a replicate here)
  #the probability so simple here because make the assumption that every parasite found the host in these well plates
  control_NLLs = -sum(dbinom(x=control[,"Infected"], prob = 1-exp(-sigma_control* control[,"P_Level"]), size = 1, log = T))
  #change it so that if a bad prediction is made, give an NA as a score; if the prediction isn't nonsense, sum the scores
  #function "any" is very useful
  if(any(!is.finite(NLLs))|!is.finite(control_NLLs)){return(NA)}else{sum(c(NLLs, control_NLLs))}
}
NLL_experiment(t=1, sigma_0 = .2, epsilon_0 = 1, df = df, control=controls)

m_control = mle2(minuslogl = NLL_experiment, start = list(epsilon_0 = .1, sigma_0 = .2))
m_control

#ok now add in susceptibility data from ganesh
#called size_grad, top of code
head(size_grad)
#get rid of false precision of there being so many sig figs to the size measurements (an artifact of imageJ)
size_grad[,"Size"]=round(size_grad[,"Size"], digits = 1)
head(size_grad)
#NLL code same as controls because parasites all got in

NLL_experiment= function(t=1, sigma_0, epsilon_0, df=dataset, control = controls, size_df = size_grad, vol = 15){
  NLLs = mapply(NLL_rep, P0=dataset[,"P_Level"]*18/vol, N_S=dataset[,"Start_Small"]/vol, N_M=dataset[,"Start_Med"]/vol, 
                N_L=dataset[,"Start_Large"]/vol, D_S=dataset[,"D_Small"], D_M=dataset[,"D_Med"],
                D_L = dataset[,"D_Large"], I_S = dataset[,"Infec_Small"], I_M = dataset[,"Infec_Med"], I_L = dataset[,"Infec_Large"],
                MoreArgs = list(t=t, sigma_0 = sigma_0, epsilon_0 = epsilon_0))
  
  #now evaluate the NLL for Becky control data to account for different infectivity of batches
  #now predict the sigma for the individuals of the body size used in the controls; will change from sigma_0 later, this is placeholder
  sigma_control = sigma_0
  #can use the density function for the binomial distribution to get the log likelihood score. We want to know the likelihood of each observation (each snail a replicate here)
  #the probability so simple here because make the assumption that every parasite found the host in these well plates
  control_NLLs = -sum(dbinom(x=control[,"Infected"], prob = 1-exp(-sigma_control* control[,"P_Level"]), size = 1, log = T))
  
  #ok, now add in susceptibility data from ganesh
  sigma_ganesh = sigma_0
  size_gradNLL = -sum(dbinom(x=size_grad[,"Infected"], prob = 1-exp(-sigma_ganesh* size_grad[,"Miracidia"]), size = 1, log = T))
  #change it so that if a bad prediction is made, give an NA as a score; if the prediction isn't nonsense, sum the scores
  #function "any" is very useful
  if(any(!is.finite(NLLs))|!is.finite(control_NLLs)|!is.finite(size_gradNLL)){return(NA)}else{sum(c(NLLs, control_NLLs, size_gradNLL))}
}
NLL_experiment(t=1, sigma_0 = .2, epsilon_0 = 1, df = df, control=controls)

m_with_suscep = mle2(minuslogl = NLL_experiment, start = list(epsilon_0 = 1, sigma_0 = .2))
m_with_suscep

##woohoo now we have a null model for all our data!!!
#ok, so now need a size deoendent model. start with the below which is very simple and smooth, it will work since we have few distinct size classes. we can go back and make it something more complicated later
#look at pg.74 lab notebook for answers
#sigma(l) = sigma0*e^sigma_L*L
#take old NLLrep and instead of sigma0/epislon0 for eevrything, need those AND size dependent ones (sigma_L and epsilon_L)
NLL_rep_size = function(t=1, sigma_0, epsilon_0, sigma_sz, epsilon_sz, P0, N_S, N_M, N_L, D_S, D_M, D_L, I_S, I_M, I_L) {
  #make size dependent things, just hardcode in the average length of each size class (so 2.5 for small snails)
  epsilon_S = epsilon_0*exp(epsilon_sz*2.5)
  epsilon_M = epsilon_0*exp(epsilon_sz*7)
  epsilon_L =epsilon_0*exp(epsilon_sz*13.5)
  sigma_S = sigma_0*exp(sigma_sz*2.5)
  sigma_M = sigma_0*exp(sigma_sz*7)
  sigma_L =sigma_0*exp(sigma_sz*13.5)
  prev_S = 1- exp(-sigma_S*(epsilon_S/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_M = 1- exp(-sigma_M*(epsilon_M/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  prev_L = 1- exp(-sigma_L*(epsilon_L/(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))*P0*(1-exp(-t*(epsilon_S*N_S+epsilon_M*N_M+epsilon_L*N_L))))
  
  NLL_S = -(dbinom(x=I_S, prob = prev_S, size = D_S, log = T))
  NLL_M = -(dbinom(x=I_M, prob = prev_M, size = D_M, log = T))
  NLL_L = -(dbinom(x=I_L, prob = prev_L, size = D_L, log = T))
  NLL_S + NLL_M + NLL_L
}
#try it and see if it works
NLL_rep_size(t=1, sigma_0=.2, epsilon_0=1, sigma_sz=.05, epsilon_sz = 0.001, P0=144, N_S=6, N_M=6, N_L=6, D_S=6, D_M=6, D_L=6, I_S=1, I_M=2, I_L=3)

NLL_experiment_size= function(t=1, sigma_0, epsilon_0, sigma_sz, epsilon_sz, df=dataset, control = controls, size_df = size_grad, vol = 15){
  NLLs = mapply(NLL_rep_size, P0=dataset[,"P_Level"]*18/vol, N_S=dataset[,"Start_Small"]/vol, N_M=dataset[,"Start_Med"]/vol, 
                N_L=dataset[,"Start_Large"]/vol, D_S=dataset[,"D_Small"], D_M=dataset[,"D_Med"],
                D_L = dataset[,"D_Large"], I_S = dataset[,"Infec_Small"], I_M = dataset[,"Infec_Med"], I_L = dataset[,"Infec_Large"],
                #this tells mapply that these things don't change throughout the datasheet
                MoreArgs = list(t=t, sigma_0 = sigma_0, epsilon_0 = epsilon_0, sigma_sz = sigma_sz, epsilon_sz = epsilon_sz))
  
  #now evaluate the NLL for Becky control data to account for different infectivity of batches
  #now predict the sigma for the individuals of the body size used in the controls; will change from sigma_0 later, this is placeholder
  sigma_control = sigma_0*exp(sigma_sz*7)
  #can use the density function for the binomial distribution to get the log likelihood score. We want to know the likelihood of each observation (each snail a replicate here)
  #the probability so simple here because make the assumption that every parasite found the host in these well plates
  control_NLLs = -sum(dbinom(x=control[,"Infected"], prob = 1-exp(-sigma_control* control[,"P_Level"]), size = 1, log = T))
  
  #ok, now add in susceptibility data from ganesh
  #code in the size dependent suscept for ganeshs data because every snail was a different size so a slightly more complicaed calculation
  size_gradNLL = -sum(dbinom(x=size_grad[,"Infected"], prob = 1-exp((-sigma_0*exp(sigma_sz*size_df[,"Size"])) * size_df[,"Miracidia"]), size = 1, log = T))
  #change it so that if a bad prediction is made, give an NA as a score; if the prediction isn't nonsense, sum the scores
  #function "any" is very useful
  if(any(!is.finite(NLLs))|!is.finite(control_NLLs)|!is.finite(size_gradNLL)){return(NA)}else{sum(c(NLLs, control_NLLs, size_gradNLL))}
}

NLL_experiment_size(t=1, sigma_0 = .2, epsilon_0 = 1, df = df, sigma_sz=0, epsilon_sz=0, control=controls, size_df = size_grad)

m_size = mle2(minuslogl = NLL_experiment_size, start = list(epsilon_0 = 1, sigma_0 = .2, sigma_sz=0, epsilon_sz=0))
m_size
#woot, have these fits now with the size dependence
#do these results make sense? sigma going down with increasing size, epsilon going up...yep, makes sense!

#now we can compare the log likes with the null!
#how to score a model? model deviance is 2*NLL to say how far the model is from the data. also need to know the complexity, which is equal to the number of free parameters it has
#good ole Akaike
#so it's -2*NLL + 2*#parameters
#aikake weights is an estimate for the probability that a given mdoel is the best model for all the models you consider
#the difference in weights is in favor or against the models
AICtab(m_null, m_size, weights = T, base = T)
#lower score is better; the difference between the (dAIC) being bigger than 10ish is good

#ok, what if you just needed exposure or just needed susceptibility, to see if you can still get good fits
#can hold a param constant in the mle to just collapse these models real easily
m_suscep_only = mle2(minuslogl = NLL_experiment_size, start = list(epsilon_0 = 1, sigma_0 = .2, sigma_sz=0, epsilon_sz=0), fixed = list(epsilon_sz=0))
m_exp_only = mle2(minuslogl = NLL_experiment_size, start = list(epsilon_0 = 1, sigma_0 = .2, sigma_sz=0, epsilon_sz=0), fixed = list(sigma_sz=0))
#can also just do the null this way
m_null = mle2(minuslogl = NLL_experiment_size, start = list(epsilon_0 = 1, sigma_0 = .2, sigma_sz=0, epsilon_sz=0), fixed = list(sigma_sz=0, epsilon_sz=0))

AICtab(m_null, m_size, m_suscep_only, m_exp_only, weights = T, base = T)

##ok so for batch effects, we can apply the function to estimate different sigma_0 but same other param values by *batch*
#can write a wrapper around the NLL function to chop the dataset into the batches, run the function on each batch to give different sigma_0 to estimate but keep other parameters the same
#then you add up 6 numbers at the end!

##BATCH IT UP
# we will hardcode in 6 batches and batch effect only on sigma_0, so not super generalizable or easy to change
#can't tell if ganesh's data has 1 or 2 batches!! will treat as one for now
NLL_experiment_size_batch= function(t=1, sigma_0_A, sigma_0_B, sigma_0_C, sigma_0_D, sigma_0_E, epsilon_0, sigma_sz, epsilon_sz, df=dataset, control = controls, size_df = size_grad, vol = 15){
  df_A=subset(dataset, Exposure.Date =="3/25/2021")
  df_B=subset(dataset, Exposure.Date =="9/24/2021")
  df_C=subset(dataset, Exposure.Date =="10/19/2021")
  df_D=subset(dataset, Exposure.Date =="11/9/2021")
  control_A=subset(control, Exposure.Date =="3/25/2021")
  control_B=subset(control, Exposure.Date =="9/24/2021")
  control_C=subset(control, Exposure.Date =="10/19/2021")
  control_D=subset(control, Exposure.Date =="11/9/2021")
  
  #NLLs = mapply(NLL_rep_size, P0=dataset[,"P_Level"]*18/vol, N_S=dataset[,"Start_Small"]/vol, N_M=dataset[,"Start_Med"]/vol, 
                #N_L=dataset[,"Start_Large"]/vol, D_S=dataset[,"D_Small"], D_M=dataset[,"D_Med"],
                #D_L = dataset[,"D_Large"], I_S = dataset[,"Infec_Small"], I_M = dataset[,"Infec_Med"], I_L = dataset[,"Infec_Large"],
                #this tells mapply that these things don't change throughout the datasheet
                #MoreArgs = list(t=t, sigma_0 = sigma_0, epsilon_0 = epsilon_0, sigma_sz = sigma_sz, epsilon_sz = epsilon_sz))
  
                #run the NLL's calculation for each of the 4 batches
              NLLs_A = mapply(NLL_rep_size, P0=df_A[,"P_Level"]*18/vol, N_S=df_A[,"Start_Small"]/vol, N_M=df_A[,"Start_Med"]/vol, 
                              N_L=df_A[,"Start_Large"]/vol, D_S=df_A[,"D_Small"], D_M=df_A[,"D_Med"],
                              D_L = df_A[,"D_Large"], I_S = df_A[,"Infec_Small"], I_M = df_A[,"Infec_Med"], I_L = df_A[,"Infec_Large"],
                              #this tells mapply that these things don't change throughout the datasheet
                MoreArgs = list(t=t, sigma_0_A = sigma_0_A, epsilon_0 = epsilon_0, sigma_sz = sigma_sz, epsilon_sz = epsilon_sz ))
  
                NLLs_B = mapply(NLL_rep_size, P0=df_B[,"P_Level"]*18/vol, N_S=df_B[,"Start_Small"]/vol, N_M=df_B[,"Start_Med"]/vol, 
                                N_L=df_B[,"Start_Large"]/vol, D_S=df_B[,"D_Small"], D_M=df_B[,"D_Med"],
                                D_L = df_B[,"D_Large"], I_S = df_B[,"Infec_Small"], I_M = df_B[,"Infec_Med"], I_L = df_B[,"Infec_Large"],
                                #this tells mapply that these things don't change throughout the datasheet
                                MoreArgs = list(t=t, sigma_0_B = sigma_0_B, epsilon_0 = epsilon_0, sigma_sz = sigma_sz, epsilon_sz = epsilon_sz ))
                
                NLLs_C = mapply(NLL_rep_size, P0=df_C[,"P_Level"]*18/vol, N_S=df_C[,"Start_Small"]/vol, N_M=df_C[,"Start_Med"]/vol, 
                                N_L=df_C[,"Start_Large"]/vol, D_S=df_C[,"D_Small"], D_M=df_C[,"D_Med"],
                                D_L = df_C[,"D_Large"], I_S = df_C[,"Infec_Small"], I_M = df_C[,"Infec_Med"], I_L = df_C[,"Infec_Large"],
                                #this tells mapply that these things don't change throughout the datasheet
                                MoreArgs = list(t=t, sigma_0_C = sigma_0_C, epsilon_0 = epsilon_0, sigma_sz = sigma_sz, epsilon_sz = epsilon_sz ))
                
                NLLs_D = mapply(NLL_rep_size, P0=df_D[,"P_Level"]*18/vol, N_S=df_D[,"Start_Small"]/vol, N_M=df_D[,"Start_Med"]/vol, 
                                N_L=df_D[,"Start_Large"]/vol, D_S=df_D[,"D_Small"], D_M=df_D[,"D_Med"],
                                D_L = df_D[,"D_Large"], I_S = df_D[,"Infec_Small"], I_M = df_D[,"Infec_Med"], I_L = df_D[,"Infec_Large"],
                                #this tells mapply that these things don't change throughout the datasheet
                                MoreArgs = list(t=t, sigma_0_D = sigma_0_D, epsilon_0 = epsilon_0, sigma_sz = sigma_sz, epsilon_sz = epsilon_sz ))
                
  #now evaluate the NLL for Becky control data to account for different infectivity of batches
  #now predict the sigma for the individuals of the body size used in the controls; will change from sigma_0 later, this is placeholder
  #sigma_control = sigma_0*exp(sigma_sz*7)
  #can use the density function for the binomial distribution to get the log likelihood score. We want to know the likelihood of each observation (each snail a replicate here)
  #the probability so simple here because make the assumption that every parasite found the host in these well plates
  #control_NLLs = -sum(dbinom(x=control[,"Infected"], prob = 1-exp(-sigma_control* control[,"P_Level"]), size = 1, log = T))
  
  sigma_control_A = sigma_0_A*exp(sigma_sz*7)
  control_NLLs_A = -sum(dbinom(x=control_A[,"Infected"], prob = 1-exp(-sigma_control_A* control_A[,"P_Level"]), size = 1, log = T))
  
  sigma_control_B = sigma_0_B*exp(sigma_sz*7)
  control_NLLs_B = -sum(dbinom(x=control_B[,"Infected"], prob = 1-exp(-sigma_control_B* control_B[,"P_Level"]), size = 1, log = T))
  
  sigma_control_C = sigma_0_C*exp(sigma_sz*7)
  control_NLLs_C = -sum(dbinom(x=control_C[,"Infected"], prob = 1-exp(-sigma_control_C* control_C[,"P_Level"]), size = 1, log = T))
  
  sigma_control_D = sigma_0_D*exp(sigma_sz*7)
  control_NLLs_D = -sum(dbinom(x=control_D[,"Infected"], prob = 1-exp(-sigma_control_D* control_D[,"P_Level"]), size = 1, log = T))
  
  #ok, now add in susceptibility data from ganesh
  #code in the size dependent suscept for ganeshs data because every snail was a different size so a slightly more complicaed calculation
  #just said sigma_E in the probability since ganesh's all one batch
  size_gradNLL = -sum(dbinom(x=size_grad[,"Infected"], prob = 1-exp((-sigma_0_E*exp(sigma_sz*size_df[,"Size"])) * size_df[,"Miracidia"]), size = 1, log = T))
  #change it so that if a bad prediction is made, give an NA as a score; if the prediction isn't nonsense, sum the scores
  #function "any" is very useful
  NLLs = c(NLLs_A, NLLs_B, NLLs_C, NLLs_D)
  control_NLLs=c(control_NLLs_A,control_NLLs_B, control_NLLs_C, control_NLLs_D)
  #if(any(!is.finite(NLLs))|!is.finite(control_NLLs)|!is.finite(size_gradNLL)){return(NA)}else{sum(c(NLLs, control_NLLs, size_gradNLL))}

  }

#check that if theyre all the same it works ok
NLL_experiment_size_batch(t=1, sigma_0_A = .2, sigma_0_B = .2, sigma_0_C = .2, sigma_0_D = .2, sigma_0_E= .2, epsilon_0 = 1, sigma_sz=0, epsilon_sz=0, control=controls, df = dataset, size_df = size_grad)
head(dataset)

####make predictions across many parasite levels to draw lines on a graph


