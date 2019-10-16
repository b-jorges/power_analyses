library(lme4)

#####function for 1 random variable, two fixed + interaction, one fixed effect is between subjects
sim_between <- function(bFix1, levelsFix1, bFix2, levelsFix2, bInteraction, Intercept, sdSubject, sdError, nParticipants, nTrials, nBlocks) {
  Subject <- rep(1:(2*nParticipants), each=nTrials*nBlocks) #Vector with subject IDs
  Fix1 <- rep(levelsFix1, each=nParticipants*nTrials*nBlocks) #Vector with IDs for first fixed effect
  Fix2 = rep(levelsFix2,nTrials*nParticipants*2) #Vector with conditions per block
  
  # random effects per subject
  S.re <- rnorm(nParticipants*nTrials*nBlocks*2, 0, sdSubject)
  
  #inter-subject error that is not correlated with effect
  eps <- rnorm(nParticipants*2*nTrials*nBlocks, 0, sdError)
  
  # put it all together
  Response <- Intercept + bFix1*(Fix1==unique(levelsFix1)[1]) + bFix2*(Fix2==unique(levelsFix2)[2]) + 
    bInteraction*(Fix1==unique(levelsFix1)[1])*(Fix2==unique(levelsFix2)[2]) +
    S.re[Subject] + eps
  
  # put into a data frame
  mydata <- data.frame(Subject = paste('s',Subject, sep=''), 
                       Fix1=Fix1, Fix2=Fix2,
                       Response = Response,
                       Block=rep((1:nBlocks),nTrials*nParticipants*2)) #add block numbers
  
  # Is there a relevant interaction between Fixed Effect 1 and Fixed Effect 2?
  fit1 <- lmer(Response ~ (Fix1*Fix2) + (1|Subject), data=mydata, REML = FALSE)
  fit2 <- lmer(Response ~ Fix1 + Fix2 + (1|Subject), data=mydata, REML = FALSE)
  anova(fit2,fit1)$`Pr(>Chisq)`[2] #Model comparison is preferrable to getting p-values directly from model
}


#####Task 2: power for effect of gender in space
Power_Gender_Task2=c()
for (i in c(3:6)){ #We simulate effect sizes for 3 to 6 subjects in each group
  out1 <- replicate(nIterations, {
    sim_between(bFix1=0, levelsFix1=c("M","F"), bFix2=0, levelsFix2=c("no","no","yes","yes","no","no","no","no"), 
                bInteraction=0.05, Intercept = 1, sdSubject=0.3, nParticipants=i, sdError=0.17, nTrials=30, nBlocks = 8)})
  #bInteraction = difference between men and women in space; estimated from current data
  #sdSubject = inter subject variability; estimated from current data among astronauts across all conditions
  #sdError = intra-subject variability: mean standard deviation across all astronauts
  #we have 30 trials in this task, and 8 blocks. LevelsFix2 referring to whether participants are in space or not.
  hist(out1) #Outputs distribution of p values for each subject
  Power_Gender_Task2 = c(Power_Gender_Task2,mean(out1 < 0.05)) #Percentage of p below 0.05
  print(i)
  i}
Power_Gender_Task2


#####Task 3: power for effect of gender in space
#analogous to Power_Gender_Task2
Power_Gender_Task3=c()
for (i in c(3:6)){
  out1 <- replicate(nIterations, {
    sim_between(bFix1=0, levelsFix1=c("M","F"), bFix2=0, levelsFix2=c("no","no","yes","yes","no","no","no","no"), 
                bInteraction=0.2, Intercept = 1, sdSubject=0.20, nParticipants=i, 
                sdError=0.23, nTrials=3, nBlocks = 8)})
  hist(out1)
  Power_Gender_Task3 = c(Power_Gender_Task3,mean(out1 < 0.05))
  print(i)
  i}
Power_Gender_Task3