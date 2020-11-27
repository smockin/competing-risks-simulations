

# ___Part 2 ___ Getting censoring rate parameter for each scenario

# Sample sizes of simulated datasets
Ns<-c(500,1000, 2000, 3000, 4000,  5000)

# Proportion(p) of subjects (with covariates equal to zero) who experience the 
# event of interest when t -> infinity.
pToUSe<- seq(from=0.1, to=0.9, by=0.1)

# Function to simulate data based on; 
# 1.) P value (pToUSe) and
# 2.) Datasets of various follow-up period(0 to 14)
#The function also determines whether lambda values have met the
#threshold (p.anyevent AND event.prop) are approaximately equal



# datasets to genereta per scenario
#datasetsToGenerate=100

N.large <- 1000000

# function to generate datasets for each scenario 

get.censoring.rate<- function(){
  cat(" Start of simulations")
  pToUSe<- seq(from=0.1, to=0.9, by=0.1)
  N.large <- 1000000 #population to sample from 
  lapply(pToUSe,
         function(p){
           load('cache/model data generating values.rda')
           
           lapply(unique(LosModels$LOS),
                  function(los){ 
                    event.prob<- unique(LosModels[LOS==los, prop.any.event])
                    
                    # Get probability of experiencing any event in the simulated data
                    
                    get.any.event.prop<- function(lambda, iter) { 
                      
                      
                      LosModels<- get("LosModels",envir = parent.frame())
                      los<- get("los",envir = parent.frame())
                      p<- get("p",envir = parent.frame())
                      N.large<- get("N.large",envir = parent.frame())
                      
                      N<- N.large
                      
                      # We generate  binary covariates with prevalences equal to those of the
                      # binary covariates in the CIN dataset.
                      
                      set.seed(iter)
                      
                      x1 <- (rbinom(N
                                    ,1
                                    ,  unique(LosModels[LOS==los, P1])
                      ))
                      
                      x2 <- (rbinom(N
                                    ,1
                                    ,  unique(LosModels[LOS==los, P2])
                      ))
                      
                      x3 <- (rbinom(N
                                    ,1
                                    ,  unique(LosModels[LOS==los, P3])
                      ))
                      
                      
                      # bind variable to a datasets
                      
                      X <- cbind(x1,x2,x3
                      )
                      
                      # Extract coefficients for main primary event fit in CIN data
                      B=LosModels[LOS==los & Model==1,.(isFemale
                                                        ,underFive
                                                        ,hasDangersigns)] %>%
                        unclass() %>% 
                        c() %>%
                        unlist() 
                      
                      # linear predictor for main event
                      XB <- X %*% B
                      
                      # Extract regressions coefficients for competing event 1 fit in CIN data
                      
                      B.other1= LosModels[LOS==los & Model==2,.(isFemale
                                                                ,underFive
                                                                ,hasDangersigns)] %>%
                        unclass() %>% 
                        c() %>%
                        unlist() 
                      
                      # linear predictor for competing event 1
                      XB.other1 <- X %*% B.other1
                      
                      # Extract regressions coefficients for competing event 2 fit in CIN data
                      
                      B.other2= LosModels[LOS==los & Model==3,.(isFemale
                                                                ,underFive
                                                                ,hasDangersigns)] %>%
                        unclass() %>% 
                        c() %>%
                        unlist() 
                      
                      # linear predictor for competing event 2
                      XB.other2 <- X %*% B.other2
                      
                      # Extract regressions coefficients for comepting event 2 fit in CIN data
                      
                      B.other3= LosModels[LOS==los & Model==4,.(isFemale
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
                      
                      # probabilities of various events
                      
                      eventProps<- cbind(p1, p2,p3, p4)
                      
                      # align probabilities to be used in the multnormial experiment
                      
                      eventProps.m<- apply(eventProps,
                                           1L,
                                           function(props){
                                             return (props/sum(props))
                                           })
                      
                      # transpose to matrix: will appear as columns for all 4 possible events
                      
                      eventProps.m<- eventProps.m %>% t()
                      
                      # Determine the type of event that occurred in an multnormial experiment
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
                      # with defined rate parameter (see page 144 of Bayesman book).
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
                      set.seed(iter)
                      T.censor <- rexp(N,lambda)
                      
                      
                      event.type <- ifelse(T.event <= T.censor,event.type,0)
                      p.anyevent <- length(event.type[event.type != 0])/N.large
                      
                      p.anyevent
                    }
                    
                    #__END of part of censoring function ___
                    
                    #___________The Bisection approach determine the censoring rate for each scenario: 
                     #the difference btn probability of primary event and that of any event should be negligible
                    low.rate <- 0.01
                    high.rate <- 100
                    mid.rate <<- (low.rate + high.rate)/2
                    lambda.store[1]<<-mid.rate
                    p.event.store[1] <<- 1
                    iter<<-1
                    
                    while (
                      abs(event.prob - p.event.store[length(p.event.store)]) > 0.0001
                    ) {
                      isConverged<- try(
                        all.equal(lambda.store[iter-2] ,
                                  lambda.store[iter-1L])
                        ,silent = T)
                      if(
                        # rate has converged
                        any(isTRUE(isConverged)) & iter>5L
                      ){
                        cat("\nNo new rates. Exiting ...")
                        break
                      }
                      set.seed(iter)
                      p.anyEvent<<- get.any.event.prop(lambda=mid.rate, iter = iter)
                      p.event.store[iter] <- p.anyEvent
                      lambda.store[iter] <<-mid.rate
                      rates.store<-data.table(LOS=los, P=p,Rate=mid.rate,
                                              Iteration=iter, 
                                              Event.prop=event.prob,
                                              AnyEvent.prop=p.event.store[iter]
                      )
                      if (p.anyEvent < event.prob){
                        high.rate <- mid.rate
                      } else{
                        low.rate <- mid.rate
                      }
                      
                      mid.rate <<- (low.rate + high.rate)/2
                      all.rates.store<<- rbindlist(list(all.rates.store,rates.store),use.names = T, fill = T)
                      save(all.rates.store, file="cache/model rate parameter-All censoring rates.rda")
                      cat(paste0("\nProgress= LOS: ",los, " P: ",p ," Iteration: ", iter))
                      
                      iter <<- iter + 1
                      cat("\n...")
                    }
                    
                    
                    toreturn<-data.table(LOS=los, P=p,Rate=mid.rate, Iteration=iter,Event.prop=event.prob,AnyEvent.prop=p.anyEvent)
                    converged.rates.store<<- rbindlist(list(converged.rates.store,toreturn),use.names = T, fill = T)
                    save(converged.rates.store, file="cache/model rate parameter-converged censoring rates.rda")
                    return(toreturn)
                  })
         })
}


all.rates.store<- data.table()

converged.rates.store<- data.table()

p.event.store<-NULL
lambda.store<-NULL
iter<- 1
rate.quantities<- get.censoring.rate()
