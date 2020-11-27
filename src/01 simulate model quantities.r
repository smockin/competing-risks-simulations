
# clear environment

rm(list=ls(all.names = T))

library(data.table)
library(magrittr)
library(survival)
library(cmprsk)
library(miscF)
# loading data from the repository

load("~/data/cleaned data.rda")


# Sample sizes of simulated datasets
Ns<-c(500,1000, 2000, 3000, 4000,  5000
)

# Proportion(p) of subjects (with covariates equal to zero) who experience the 
# event of interest when t -> infinity.
pToUSe<- seq(from=0.1, to=0.9, by=0.1)

# Function to simulate data based on; 
# 1.) P value (pToUSe) and
# 2.) Datasets of various follow-up period(0 to 14)
#The function also determines whether lambda values have met the
#threshold (p.anyevent AND event.prop) are approaximately equal


# datasets to genereta per scenario
datasetsToGenerate=1000 

all.dataOut.store<<- data.table()
load("cache/model rate parameter-converged censoring rates.rda")
rateParametersToUse<- converged.rates.store
rm(converged.rates.store)
load("cache/model data generating values.rda")

losModelstouse<- LosModels


# function to generate datasets for each scenario 



get.Model.Quantities<- function(){
  lapply(Ns,
         function(N){
           lapply(pToUSe,
                  function(p){
                    lapply(unique(losModelstouse$LOS),
                           function(los.n){ 
                             lapply(1:datasetsToGenerate
                                    , function(iter){
                                      # We generate  binary covariates with prevalences equal to those of the
                                      # binary covariates in the CIN dataset.
                                      
                                      set.seed(iter)
                                      
                                      x1 <- (rbinom(N
                                                    ,1
                                                    ,  unique(losModelstouse[LOS==los.n, P1])
                                      ))
                                      
                                      x2 <- (rbinom(N
                                                    ,1
                                                    ,  unique(losModelstouse[LOS==los.n, P2])
                                      ))
                                      
                                      
                                      x3 <- (rbinom(N
                                                    ,1
                                                    ,  unique(losModelstouse[LOS==los.n, P3])
                                      ))
                                      
                                      
                                      # bind variable to a datasets
                                      
                                      X <- cbind(x1,x2,x3
                                                
                                      )
                                      
                                      # Extract coefficients for main primary event fit in CIN data
                                      B=losModelstouse[LOS==los.n & Model==1,.(isFemale
                                                                               ,underFive
                                                                               ,hasDangersigns)] %>%
                                        unclass() %>% 
                                        c() %>%
                                        unlist() 
                                      
                                      # linear predictor for main event
                                      XB <- X %*% B
                                      
                                      # Extract regressions coefficients for competing event 1 fit in CIN data
                                      
                                      B.other1= losModelstouse[LOS==los.n & Model==2,.(isFemale
                                                                                       ,underFive
                                                                                       ,hasDangersigns)] %>%
                                        unclass() %>% 
                                        c() %>%
                                        unlist() 
                                      
                                      # linear predictor for competing event 1
                                      XB.other1 <- X %*% B.other1
                                      
                                      # Extract regressions coefficients for competing event 2 fit in CIN data
                                      
                                      B.other2= losModelstouse[LOS==los.n & Model==3,.(isFemale
                                                                                       ,underFive
                                                                                       ,hasDangersigns)] %>%
                                        unclass() %>% 
                                        c() %>%
                                        unlist() 
                                      
                                      # linear predictor for competing event 2
                                      XB.other2 <- X %*% B.other2
                                      
                                      # Extract regressions coefficients for comepting event 2 fit in CIN data
                                      
                                      B.other3= losModelstouse[LOS==los.n & Model==4,.(isFemale
                                                                                       ,underFive
                                                                                       ,hasDangersigns)] %>%
                                        unclass() %>% 
                                        c() %>%
                                        unlist() 
                                      
                                      # linear predictor for competing event 3
                                      XB.other3 <- X %*% B.other3
                                      
                                      # The regression coefficients are obtained from the 4 models fitted earlier
                                      
                                      p1 <- 1 - (1-p)^(exp(XB))
                                      p2 <- 1 - (1-p)^(exp(XB.other1))
                                      p3 <- 1 - (1-p)^(exp(XB.other2))
                                      p4 <- 1 - (1-p)^(exp(XB.other3))
                                      
                                      
                                      # probailities of various events
                                      eventProps<- cbind(p1, p2,p3, p4)
                                      
                                      # align probabilities to be used in the multnormial experiment
                                      
                                      eventProps.m<- apply(eventProps,
                                                           1L,
                                                           function(props){
                                                             return (props/sum(props))
                                                           })
                                      
                                      # transpose to matrix: will appear as columns for all 4 possible events
                                      
                                      eventProps.m<- eventProps.m %>% t()
                                      
                                      # Determine the type of event that occurred in a multnormial experiment
                                      # based on matrix of probabilities, per subject it will return 1 for event 1, 2 for event 2, 3 for event  3, or 4 for event 4 ).
                                      
                                      set.seed(iter)
                                      
                                      event.type<- rMultinom(p=eventProps.m) 
                                      
                                      # Determine whether a type 1 event occurred.
                                      
                                      event.type1<- ifelse(event.type==1, 1, 0)
                                      
                                      # Determine whether other events occurred.
                                      
                                      event.type2<- ifelse(event.type==2, 1, 0)
                                      event.type3<- ifelse(event.type==3, 1, 0)
                                      event.type4<- ifelse(event.type==4, 1, 0)
                                      
                                      # The distribution of the competing events follow an Exponential distribution
                                      # with defined rate parameter (see page 144).
                                      set.seed(iter)
                                      
                                      T2 <- rexp(N,exp(XB.other1))
                                      
                                      set.seed(iter)
                                      
                                      T3 <- rexp(N,exp(XB.other2))
                                      
                                      set.seed(iter)
                                      
                                      T4 <- rexp(N,exp(XB.other3))
                                      
                                      # For T1 distribution, we have a formulas for the CDF. We inverted this in 
                                      # Mathematica, and then evaluate it at u ~ U(0,1)
                                      
                                      set.seed(iter)
                                      
                                      u <- runif(N)
                                      T1 <- -log(-(1-((-u + 1/(1-(1-p)^exp(XB))) * (1-(1-p)^exp(XB)) )^(1/exp(XB)) - p)/p)
                                      
                                      T.event <- (T1*event.type1) + (T2*event.type2) + (T3*event.type3) + (T4*event.type4)
                                      
                                      # Time at which subjects are censored.
                                      
                                      # load the appropriate rate parameter(simulated)
                                      
                                      searched.rates.parameter<- rateParametersToUse
                                      setDT(searched.rates.parameter)
                                      lambda<-  searched.rates.parameter[LOS==los.n & P==p, Rate]
                                      
                                      set.seed(iter)
                                      T.censor <- rexp(N,lambda)
                                      
                                      event.type <- ifelse(T.event <= T.censor,event.type,0)
                                      T.event <- ifelse(T.event <= T.censor,T.event,T.censor)
                                      
                                      
                                      # check if the proportion of deaths match that of the true data
                                      
                                      # confidence interval about the proportion of simulated primary events
                                      propCI=prop.test(n=length(event.type),x = sum(event.type==1))
                                      propCI=as.character(propCI$conf.int)
                                      propCI.L=propCI[1]
                                      propCI.U=propCI[2]
                                      
                                      
                                      # combine all these summaries into one data object to be returned
                                      toreturn<- data.table::data.table(
                                        iter=iter,
                                        LOS=los.n,
                                        N=N,
                                        P=p
                                      )
                                      
                                      empirical.ds<- copy(dsToUse)
                                      #main event
                                      true.anyEvent.prop= (empirical.ds[los<=los.n
                                                                        ,event.type=="Died"] %>% table() %>% prop.table())["TRUE"]
                                     # probability of any event
                                       sim.anyEvent.prop= sum(event.type==1)/length(event.type)
                                      
                                      toreturn[, true.anyEvent.prop:=true.anyEvent.prop]
                                      toreturn[, sim.anyEvent.prop:=sim.anyEvent.prop]
                                      toreturn[, sim.anyEvent.LCI:=propCI.L]
                                      toreturn[, sim.anyEvent.UCI:=propCI.U]
                                      
                                      all.dataOut.store<<- rbindlist(list(all.dataOut.store,toreturn))
                                      save(all.dataOut.store, file="cache/simulated model estimates.rda")
                                      
                                      cat(paste0("\nProgress= LOS: ",los.n, " P: ",p ," Iteration: ", iter," N: ", N))
                                      return(toreturn)
                                      
                                    })
                           })
                  })
         })
}

all.dataOut.store<- data.table()

quantities.data<- get.Model.Quantities()






