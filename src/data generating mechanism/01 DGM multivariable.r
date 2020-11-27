


#Required libs

library(data.table)
library(magrittr)
library(survival)
library(cmprsk)
library(miscF)
# loading data from the repository

load("~/cleaned data.rda")

#'
# ___Part 1 ___ Getting True Data Generating values in different scenarios
LosModels<-NULL
LosModels<<- split(dsToUse, dsToUse$los) %>%
  lapply(function(data.los){ 
    setDT(dsToUse)
    los.num<- unique(data.los$los)
    dataset<- copy(dsToUse[los<=los.num, .SD])
    data.los<-copy(dataset)
    # Create variable denoting occurrence of in-hospital deaths(primary event).
    dataset$died <- ifelse(dataset$event.type=="Died" ,1,0)
    
    n.events<-sum(dataset$died, na.rm=T)
    # The proportion of subjects for whom ANY event is observed.
    prop.any.event<-n.events/nrow(data.los)
    
    # Create variable denoting occurrence of competing event 1. 
    dataset$discharged <- ifelse(dataset$event.type == "Discharged",2,0)
    
    # Create variable denoting occurrence of competing event 2.
    dataset$refer <- ifelse(dataset$event.type == "Refer",3,0)
    
    # Create variable denoting occurrence of competing event 3. 
    dataset$discharged.a.advc <- ifelse(dataset$event.type == "Discharged against advice",4,0)
    
    
    # Create event indicator for any type of event observed(NB absconded is censored).
    dataset$event <- ifelse(dataset$event.type %in%
                              c("Discharged",
                                "Died",
                                "Refer",
                                "Discharged against advice")
                            ,1,0)
    
    
    ################################################################################
    # Fit Fine-Gray competing risk regression model to the emprical sample(CIN data).
    # This will be used to set the regression coefficients in the simulations.
    ################################################################################
    suppressMessages(attach(dataset))
    
    cat(paste0("\nFitting models for LOS= ",los.num, " ..."))
    
    # primary event: Died
    crr1 <- crr(ftime=los,fstatus=died,
                cov1=cbind(isFemale
                           ,underFive
                           ,hasDangersigns
                           
                ),
                failcode=1,cencode=0)
    cat("\nModel 1 fitted")
    
    # competing event: Discharged
    crr2 <- crr(ftime=los,fstatus=discharged,
                cov1=cbind(isFemale
                           ,underFive
                           ,hasDangersigns 
                ),
                failcode=2,cencode=0)
    cat("\nModel 2 fitted")
    
    # competing event: referrals
    crr3 <- crr(ftime=los,fstatus=refer,
                cov1=cbind(isFemale
                           ,underFive
                           ,hasDangersigns
                ),
                failcode=3,cencode=0)
    
    cat("\n Model 3 fitted")
    
    # competing event: discharged against advice
    crr4 <- crr(ftime=los,fstatus=discharged.a.advc,
                cov1=cbind(isFemale
                           ,underFive
                           ,hasDangersigns
                ),
                failcode=4,cencode=0)
    
    cat("\nModel 4 fitted")
    
    #Prevalence to be used for Plasmode-type sinulations
    p1<-mean(isFemale==T, na.rm = T)
    p2<-mean(underFive==T, na.rm = T)
    p3<-mean(hasDangersigns==T, na.rm = T)
    
    #_______Extract Model estimates_________
    
    #model 1
    model1=data.table(Model=1)
    B <- crr1$coef
    cnames<- names(B)
    model1coefs<- as.numeric(B)
    model1[, (cnames):=as.list(model1coefs)]
    
    #model 2
    B.other1 <- crr2$coef
    
    model2=data.table(Model=2)
    
    cnames2<- names(B.other1)
    model2coefs<- as.numeric(B.other1)
    model2[, (cnames2):=as.list(model2coefs)]
    
    # model 3
    B.other2 <- crr3$coef
    
    model3=data.table(Model=3)
    
    cnames3<- names(B.other2)
    model3coefs<- as.numeric(B.other2)
    model3[, (cnames3):=as.list(model3coefs)]
    
    # model 4
    B.other3 <- crr4$coef
    
    model4=data.table(Model=4)
    
    cnames4<- names(B.other3)
    model4coefs<- as.numeric(B.other3)
    model4[, (cnames4):=as.list(model4coefs)]
    
    # bind estimates into one
    models.coefs<- rbindlist(list(model1, model2, model3, model4))
    models.coefs[, LOS:=los.num]
    models.coefs[, prop.any.event:=prop.any.event]
    models.coefs[, P1:=p1]
    models.coefs[, P2:=p2]
    models.coefs[, P3:=p3]
    
    suppressWarnings(detach(dataset))
    return(models.coefs)
  })->tobind
tobind<- tobind %>% do.call(rbind,.)
LosModels<- tobind

