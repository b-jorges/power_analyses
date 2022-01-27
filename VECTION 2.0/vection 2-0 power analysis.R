require(lme4)
require(dplyr)
require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


SimulateData = function(nReps = 1,nParticipants = 6){
  
  StartVelocity = 0
  Acceleration = 0.8
  TargetDistances = c(6,10,14,18,22,26,30,34,38,42)
  Sessions = c("BDC1","EarlyFlight1","EarlyFlight2","EarlyFlight3","LateFlight","EarlyPostFlight1","EarlyPostFlight2","EarlyPostFlight3","LatePostFlight")
  Jitter2 = c("No Jitter","Vertical Jitter", "Horizontal Jitter")
  
  #make a dataframe with all conditions
  Dataframe = expand.grid(
    StartVelocity = StartVelocity,
    Accleration = Acceleration,
    TargetDistance = TargetDistances,
    Session = Sessions,
    Jitter2 = Jitter2,
    Rep = 1:nReps,
    ID = paste0("s",1:nParticipants))
  
  #we expect both types of jitter to have the same effect, so we won't distinguish between them for the sake of the power analysis
  Dataframe = Dataframe %>% 
    mutate(Jitter = case_when(Jitter2 == "No Jitter" ~ "1No Jitter",
                              Jitter2 == "Vertical Jitter" ~ "Jitter",
                              Jitter2 == "Horizontal Jitter" ~ "Jitter")) %>% 
    #the baseline k is 1
    mutate(k_session = case_when(
      Session == "BDC1" ~ 1,
      Session == "EarlyFlight1" ~ 1,
      Session == "EarlyFlight2" ~ 1,
      Session == "EarlyFlight3" ~ 1,
      Session == "LateFlight" ~ 1,
      #we add the expected effect upon return to earth, with some variability in the effect. The earlier, the stronger the effect
      Session == "EarlyPostFlight1" ~ 1 - abs(rnorm(length(ID),0.15,0.15*0.15)),
      Session == "EarlyPostFlight2" ~ 1 - abs(rnorm(length(ID),0.125,0.125*0.15)),
      Session == "EarlyPostFlight3" ~ 1 - abs(rnorm(length(ID),0.1,0.1*0.15)), #this is roughly the VECTION 1.0 time after launch, therefore we use the effect observed in VECTION 1.0
      Session == "LatePostFlight" ~ 1),
      alpha_session = case_when(
        #baseline alpha is 0
        Session == "BDC1" ~ 0,
        #again, stronger effect earlier after launch, again with some variability in the strength of the effect
        Session == "EarlyFlight1" ~ rnorm(length(ID),0.015,0.015*0.15),
        Session == "EarlyFlight2" ~ rnorm(length(ID),0.0125,0.0125*0.15),
        Session == "EarlyFlight3" ~ rnorm(length(ID),0.01,0.01*0.15), #this is roughly the VECTION 1.0 time after launch, therefore we use the effect observed in VECTION 1.0
        Session == "LateFlight" ~ 0,
        Session == "EarlyPostFlight1" ~ 0,
        Session == "EarlyPostFlight2" ~ 0,
        Session == "EarlyPostFlight3" ~ 0,
        Session == "LatePostFlight" ~ 0),
      #Jitter to remedy the effect we observe (fully), but again with some variability 
      #in the way it offsets the biases in response to microgravity exposure
      JitterEffect_k = case_when(
        Session == "BDC1" ~ 0,
        Session == "EarlyFlight1" ~ 0,
        Session == "EarlyFlight2" ~ 0,
        Session == "EarlyFlight3" ~ 0,
        Session == "LateFlight" ~ 0,
        Session == "EarlyPostFlight1" & Jitter == "Jitter" ~ rnorm(length(ID),0.15,0.15*0.15),
        Session == "EarlyPostFlight2" & Jitter == "Jitter" ~ rnorm(length(ID),0.125,0.125*0.15),
        Session == "EarlyPostFlight3" & Jitter == "Jitter" ~ rnorm(length(ID),0.1,0.1*0.15),
        Session == "EarlyPostFlight1" & Jitter == "1No Jitter" ~ 0,
        Session == "EarlyPostFlight2" & Jitter == "1No Jitter" ~ 0,
        Session == "EarlyPostFlight3" & Jitter == "1No Jitter" ~ 0,
        Session == "LatePostFlight" ~ 0),
      JitterEffect_alpha = case_when(
        Session == "BDC1" ~ 0,
        Session == "EarlyFlight1" & Jitter == "Jitter" ~ rnorm(length(ID),-0.015,0.015*0.15),
        Session == "EarlyFlight2" & Jitter == "Jitter" ~ rnorm(length(ID),-0.0125,0.0125*0.15),
        Session == "EarlyFlight3" & Jitter == "Jitter" ~ rnorm(length(ID),-0.01,0.01*0.15),
        Session == "EarlyFlight1" & Jitter == "1No Jitter" ~ 0,
        Session == "EarlyFlight2" & Jitter == "1No Jitter" ~ 0,
        Session == "EarlyFlight3" & Jitter == "1No Jitter" ~ 0,
        Session == "LateFlight" ~ 0,
        Session == "EarlyPostFlight1" ~ 0,
        Session == "EarlyPostFlight2" ~ 0,
        Session == "EarlyPostFlight3" ~ 0,
        Session == "LatePostFlight" ~ 0))
  
  #each person has a different baseline k and alpha, this is from the VECTION 1.0 data
  Dataframe = Dataframe %>%
    group_by(ID) %>%
    mutate(k_PerPerson = rnorm(1,1,0.2),
           alpha_PerPerson = abs(rnorm(1,0.002,0.003)))
  
  #we put together all the values drawn above and use the Lappe model to predict values for each combination of k and alpha per trial
  Dataframe = Dataframe %>% 
    mutate(SessionEffect = k_PerPerson*k_session,
           JitterEffect = k_PerPerson*JitterEffect_k,
           k = rnorm(length(ID),k_PerPerson*k_session+k_PerPerson*JitterEffect_k,0.2),
           alpha = abs(rnorm(length(ID),alpha_PerPerson + alpha_session + JitterEffect_alpha,abs(alpha_PerPerson + alpha_session + JitterEffect_alpha)*0.2)),
           PredictionLappeCountingUp_Obj = (k/alpha)*(1-exp(-TargetDistance*alpha)), #Adjust Target
           PredictionLappeCountingDown_Obj = (log(TargetDistance + k/alpha) - log(k/alpha))/alpha)
  
  #crude outlier analysis
  Dataframe %>% filter(PredictionLappeCountingDown_Obj < 2*TargetDistance &
                         PredictionLappeCountingUp_Obj < 2*TargetDistance)
}

###########fitting the data#########
GetLappeParameters = function(Dataframe){
  #in this function we fit the simulated data to the Lappe model
  #first we get the mean simulated travel distance per condition
  DataframeLappe_Simulated = Dataframe %>%
    group_by(Jitter,TargetDistance,ID,Session) %>%
    mutate(MTT_MeanPerDistance = mean(PredictionLappeCountingUp_Obj),
           AT_MeanPerDistance = mean(PredictionLappeCountingDown_Obj)) %>%
    slice(1) %>%
    select(Jitter,TargetDistance,ID,Session,
           MTT_MeanPerDistance,
           AT_MeanPerDistance,
           Session,
           PredictionLappeCountingUp_Obj,
           PredictionLappeCountingDown_Obj)
  
  #then we define an error function that computes the RMSE between the simulated data and fitted lappe values
  errfn_Lappe = function(p,
                         VectorRealDistances,
                         VectorResponseMTT,
                         VectorResponseAT){
    
    if(p[2] > 0){
      RMSE_MTT = mean((VectorResponseMTT-lappeMTT(VectorRealDistances, p[1], p[2]))^2)^0.5 #p(1)= gain, p(2)= decay
      RMSE_AT = mean((VectorResponseAT-lappeAT(VectorRealDistances, p[1], p[2]))^2)^0.5 #p(1)= gain, p(2)= decay
      # print(RMSE_MTT)
      # print(RMSE_AT)
      mean(c(RMSE_MTT,RMSE_AT))    
    }
    else{
      10000
    }
    
  }
  
  initalParam_Lappe <- c(1, 0.01)
  
  #this is the lappe model for the move to target task
  lappeMTT = function(d0, gain, alpha){
    (gain/alpha)*(1-exp(-d0*alpha))
  }
  
  #this is the lappe model for the adjust target task
  lappeAT = function(d0, gain, alpha){
    (log(d0 + gain/alpha) - log(gain/alpha))/alpha
  }
  
  #here we use the error function defined above to compute the lappe gain, lappe decay and RMSE for each condition, session, and participant
  DataframeLappe_Simulated = DataframeLappe_Simulated[complete.cases(DataframeLappe_Simulated),] %>% 
    group_by(Jitter,Session,ID) %>%
    mutate(LappeGain = optim(initalParam_Lappe, 
                             errfn_Lappe, 
                             VectorRealDistances = TargetDistance, 
                             VectorResponseMTT = MTT_MeanPerDistance,
                             VectorResponseAT = AT_MeanPerDistance)$par[1],
           LappeDecay = optim(initalParam_Lappe, 
                              errfn_Lappe, 
                              VectorRealDistances = TargetDistance, 
                              VectorResponseMTT = MTT_MeanPerDistance,
                              VectorResponseAT = MTT_MeanPerDistance)$par[2],
           RMSE = optim(initalParam_Lappe, 
                        errfn_Lappe, 
                        VectorRealDistances = TargetDistance, 
                        VectorResponseMTT = MTT_MeanPerDistance,
                        VectorResponseAT = MTT_MeanPerDistance)$value,
           
           Prediction_MTT_Lappe = lappeMTT(TargetDistance,LappeGain,LappeDecay),
           Prediction_AT_Lappe = lappeAT(TargetDistance,LappeGain,LappeDecay)) %>%
    group_by(Jitter,Session,ID,TargetDistance) %>% 
    slice(1)
  
  DataframeLappe_Simulated
}

#then we use the functions built above to simulate 200 datasets and perform the planned analysis over each of these simulated datasets
#we save the p value for each contrast we are interested in, and the fraction of cases where p < 0.05 out of these 200 simulated datasets is the power
Powerframe_Vection_2.0 = data.frame()
nIterationsLower = 1
nIterationsUpper = 200
time_init_overall = Sys.time()

#this takes several hours to simulate. You can instead load the results of the simulation with the following line:
load(file=paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SavedVariables/Powerframe_Vection_2.0.RData"))

#if you want to run the simulations yourself, uncomment the following for loop
# for (iteration in nIterationsLower:nIterationsUpper){
# 
#   time_init = Sys.time()
# 
#   for (j in 1:3){
#     for (k in c(6,8,10)){
# 
#       Dataframe = SimulateData(nReps = j, nParticipants = k)
#       DataframeLappe_Simulated = GetLappeParameters(Dataframe)
# 
#       ModelGain = lmer(LappeGain ~ Session*Jitter + (Session + Jitter|ID),
#                        data = DataframeLappe_Simulated %>%
#                          group_by(Jitter,Session,ID) %>%
#                          slice(1))
#       summary(ModelGain)
#       ModelDecay = lmer(LappeDecay ~ Session*Jitter + (Session + Jitter|ID),
#                         data = DataframeLappe_Simulated %>%
#                           group_by(Jitter,Session,ID) %>%
#                           slice(1))
#       summary(ModelDecay)
# 
#       Powerframe_Vection_2.0 = rbind(Powerframe_Vection_2.0,
#                                      data.frame(Label = rownames(summary(ModelGain)$coef),
#                                                 pvalue = summary(ModelGain)$coef[,"Pr(>|t|)"],
#                                                 coefficient = summary(ModelGain)$coef[,"Estimate"],
#                                                 DependentVariable = "LappeGain",
#                                                 nParticipants = k,
#                                                 nReps = j,
#                                                 nIteration = iteration),
#                                      data.frame(Label = rownames(summary(ModelDecay)$coef),
#                                                 pvalue = summary(ModelDecay)$coef[,"Pr(>|t|)"],
#                                                 coefficient = summary(ModelDecay)$coef[,"Estimate"],
#                                                 DependentVariable = "LappeDecay",
#                                                 nParticipants = k,
#                                                 nReps = j,
#                                                 nIteration = iteration))
# 
#       save(Powerframe_Vection_2.0, file = paste0(dirname(rstudioapi::getSourceEditorContext()$path),
#                                                  "/SavedVariables/Powerframe_Vection_2.0.RData"))
#       load(file=paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/SavedVariables/Powerframe_Vection_2.0.RData"))
# 
#     }
#   }
# 
#   print(paste0("Iteration #", iteration, " has passed. It took ", round(difftime(Sys.time(), time_init, units='mins'),2), " minute(s). That means, the whole thing should take about ",
#                round((nIterationsUpper-nIterationsLower-1)*difftime(Sys.time(), time_init, units='mins'),2), " minute(s). ", round(difftime(Sys.time(), time_init_overall, units='mins'),2),
#                " minutes have already passed."))
# }



#Just add some prettier labels:
Power = Powerframe_Vection_2.0 %>% 
  group_by(DependentVariable,nParticipants,nReps,Label) %>% 
  mutate(Power = mean(pvalue < 0.05),
         Labels2 = case_when(
           Label == "SessionEarlyPostFlight1" ~ "Early Post Flight 1",
           Label == "SessionEarlyPostFlight2"~ "Early Post Flight 2",
           Label == "SessionEarlyPostFlight3"~ "Early Post Flight 3",
           Label == "SessionEarlyPostFlight1:JitterJitter"~ "Early Post Flight 1 * Jitter",
           Label == "SessionEarlyPostFlight2:JitterJitter"~ "Early Post Flight 2 * Jitter",
           Label == "SessionEarlyPostFlight3:JitterJitter"~ "Early Post Flight 3 * Jitter",
           Label == "SessionEarlyFlight1" ~ "Early Flight 1",
           Label == "SessionEarlyFlight2"~ "Early Flight 2",
           Label == "SessionEarlyFlight3"~ "Early Flight 3",
           Label == "SessionEarlyFlight1:JitterJitter"~ "Early Flight 1 * Jitter",
           Label == "SessionEarlyFlight2:JitterJitter"~ "Early Flight 2 * Jitter",
           Label == "SessionEarlyFlight3:JitterJitter"~ "Early Flight 3 * Jitter"))  %>% 
  slice(1)

#plot for the power for the lappe gain
LappeGain = ggplot(Power %>% filter(DependentVariable == "LappeGain" & 
                          Label %in% c("SessionEarlyPostFlight1",
                                       "SessionEarlyPostFlight2",
                                       "SessionEarlyPostFlight3",
                                       "SessionEarlyPostFlight1:JitterJitter",
                                       "SessionEarlyPostFlight2:JitterJitter",
                                       "SessionEarlyPostFlight3:JitterJitter")),
       aes(nReps,Power,color = as.factor(nParticipants))) +
  geom_line(size = 1.5) +
  facet_wrap(.~Labels2) +
  scale_color_discrete(name = "", 
                       labels = c("6 Participants","8 Participants","10 Participants")) +
  scale_x_continuous(breaks = c(1,2,3), labels = c("1\n       (10 min)","2\n(20 min)","3\n(30 min)       ")) +
  xlab("") +
  geom_hline(yintercept = 0.8, linetype = 2) +
  geom_hline(yintercept = 0.9, linetype = 4) +
  geom_hline(yintercept = 0.95, linetype = 3) +
  ggtitle("A. Power for Gain Parameter")

#plot for the power for the lappe decay
LappeDecay = ggplot(Power %>% filter(DependentVariable == "LappeDecay" & 
                          Label %in% c("SessionEarlyFlight1",
                                       "SessionEarlyFlight2",
                                       "SessionEarlyFlight3",
                                       "SessionEarlyFlight1:JitterJitter",
                                       "SessionEarlyFlight2:JitterJitter",
                                       "SessionEarlyFlight3:JitterJitter")),
       aes(nReps,Power,color = as.factor(nParticipants))) +
  geom_line(size = 1.5) +
  facet_wrap(.~Labels2) +
  scale_color_discrete(name = "", 
                       labels = c("6 Participants","8 Participants","10 Participants")) +
  scale_x_continuous(breaks = c(1,2,3), labels = c("1\n       (10 min)","2\n(20 min)","3\n(30 min)       ")) +
  xlab("Repetitions per Condition\n(Crew Time per Session)") +
  geom_hline(yintercept = 0.8, linetype = 2) +
  geom_hline(yintercept = 0.9, linetype = 4) +
  geom_hline(yintercept = 0.95, linetype = 3) +
  ggtitle("B. Power for Decay Parameter")

#put both plots together and save them
plot_grid(LappeGain,LappeDecay,nrow = 2)
ggsave("Power Lappe VECTION 2.0.jpg", w = 10, h = 8)