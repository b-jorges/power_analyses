library(lme4)

#####function for 1 random variable, two fixed + interaction, one fixed effect is between subjects
sim_between <- function(bFix1, levelsFix1, bFix2, levelsFix2, bInteraction, Intercept, sdSubject, sdError, nParticipants, nTrials, nBlocks) {
  Subject <- rep(1:(2*nParticipants), each=nTrials*nBlocks) #Vector with subject IDs
  Fix1 <- rep(levelsFix1, each=nParticipants*nTrials*nBlocks) #Vector with IDs for first fixed effect
  
  #  assume frequency is constant accross word, random from 1-100
  Fix2 = rep(levelsFix2,nTrials*nParticipants*2) 
  
  # random effects per subject
  S.re <- rnorm(nParticipants*nTrials*nBlocks*2, 0, sdSubject)
  
  # epsilons----------------------------------------------------------------------------------------
  eps <- rnorm(nParticipants*2*nTrials*nBlocks, 0, sdError)
  
  # put it all together
  Response <- Intercept + bFix1*(Fix1==unique(levelsFix1)[1]) + bFix2*(Fix2==unique(levelsFix2)[2]) + 
    bInteraction*(Fix1==unique(levelsFix1)[1])*(Fix2==unique(levelsFix2)[2]) +
    S.re[Subject] + eps
  
  # put into a data frame
  mydata <- data.frame(Subject = paste('s',Subject, sep=''), 
                       Fix1=Fix1, Fix2=Fix2,
                       Response = Response,
                       Block=rep((1:nBlocks),nTrials*nParticipants*2))
  
  # analyze looking at interaction term with LR test
  fit1 <- lmer(Response ~ (Fix1*Fix2) + (1|Subject), data=mydata, REML = FALSE)
  fit2 <- lmer(Response ~ Fix1 + Fix2 + (1|Subject), data=mydata, REML = FALSE)
  summary(fit1)
  anova(fit2,fit1)$`Pr(>Chisq)`[2]
}

#####function for 1 random variable, two fixed + interaction, both fixed effects are within subjects
sim_within <- function(bFix1, levelsFix1, bFix2, levelsFix2, bInteraction, Intercept, sdSubject, sdError, nParticipants, nTrials, nBlocks) {
  Subject <- rep( 1:(nParticipants), each=nTrials*nBlocks) #Vector with subject IDs
  Fix1 <- rep(levelsFix1, each=nParticipants*nTrials*nBlocks/2) #Vector with IDs for first fixed effect
  
  #  assume frequency is constant accross word, random from 1-100
  Fix2 = rep(levelsFix2,nTrials*nParticipants) 
  
  # random effects per subject
  S.re <- rnorm(nParticipants*nTrials*nBlocks, 0, sdSubject)
  
  # epsilons----------------------------------------------------------------------------------------
  eps <- rnorm(nParticipants*nTrials*nBlocks, 0, sdError)
  
  # put it all together
  Response <- Intercept + bFix1*(Fix1==unique(levelsFix1)[1]) + bFix2*(Fix2==unique(levelsFix2)[2]) + 
    bInteraction*(Fix1==unique(levelsFix1)[1])*(Fix2==unique(levelsFix2)[2]) +
    S.re[Subject] + eps
  
  # put into a data frame
  mydata <- data.frame(Subject = paste('s',Subject, sep=''), 
                       Fix1=Fix1, Fix2=Fix2,
                       Response = Response,
                       Block=rep((1:nBlocks),nTrials*nParticipants))
  
  # analyze looking at interaction term with LR test
  fit1 <- lmer(Response ~ (Fix1*Fix2) + (1|Subject), data=mydata, REML = FALSE)
  fit2 <- lmer(Response ~ Fix1 + Fix2 + (1|Subject), data=mydata, REML = FALSE)
  
  anova(fit2,fit1)$`Pr(>Chisq)`[2]
}

#####function for 1 random variable, one fixed
sim_Main_Task2 <- function(bFix, levelsFix, Intercept, sdSubject, sdError, nParticipants, nTrials, nBlocks) {
  Subject <- rep(1:(2*nParticipants), each=nTrials*nBlocks) #Vector with subject IDs
  Fix <- rep(levelsFix, nParticipants*nTrials) #Vector with IDs for first fixed effect
  
  # random effects per subject
  S.re <- rnorm(nParticipants*nTrials*nBlocks*2, 0, sdSubject)
  
  # epsilons----------------------------------------------------------------------------------------
  eps <- rnorm(nParticipants*2*nTrials*nBlocks, 0, sdError)
  
  # put it all together
  Response <- Intercept + bFix*(Fix==unique(levelsFix)[2]) +
    S.re[Subject] + eps
  
  # put into a data frame
  mydata <- data.frame(Subject = paste('s',Subject, sep=''), 
                       Fix=Fix, 
                       Response = Response,
                       Block=rep((1:nBlocks),nTrials*nParticipants*2))
  
  # analyze looking at interaction term with LR test
  fit1 <- lmer(Response ~ Fix + (1|Subject), data=mydata, REML = FALSE)
  fit2 <- lmer(Response ~ + (1|Subject), data=mydata, REML = FALSE)
  
  anova(fit2,fit1)$`Pr(>Chisq)`[2]
}

#####how precise do we want to simulate power? 500 simulations per effect seems reasonable
nIterations = 500

#####Task 2: power for effect of gender in space
Power_Gender_Task2=c()
for (i in c(3:6)){
  out1 <- replicate(nIterations, {
    sim_between(bFix1=0, levelsFix1=c("M","F"), bFix2=0, levelsFix2=c("no","no","yes","yes","no","no","no","no"), 
                bInteraction=0.075, Intercept = 1, sdSubject=0.27, nParticipants=i, sdError=0.19*1.5, nTrials=30, nBlocks = 8)})
  hist(out1)
  Power_Gender_Task2 = c(Power_Gender_Task2,mean(out1 < 0.05))
  print(i)
  i}
Power_Gender_Task2

#####Task 3: power for effect of gender in space
Power_Gender_Task3=c()
for (i in c(3:6)){
  out1 <- replicate(nIterations, {
    sim_between(bFix1=0, levelsFix1=c("M","F"), bFix2=0, levelsFix2=c("no","no","yes","no","no","no","no","no"), 
                bInteraction=0.15, Intercept = 1, sdSubject=0.36, nParticipants=i, 
                sdError=0.12, nTrials=3, nBlocks = 8)})
  hist(out1)
  Power_Gender_Task3 = c(Power_Gender_Task3,mean(out1 < 0.05))
  print(i)
  i}
Power_Gender_Task3
