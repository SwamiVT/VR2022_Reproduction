## Data analysis script accompanying:
## Veríssimo, J., Verhaeghen, P., Goldman, N., Weinstein, M., & Ullman, M. T. (2021)
## Evidence that aging yields improvements as well as declines across attention and executive functions

## João Veríssimo
## Email: jlverissimo@edu.ulisboa.pt
## 30 November 2020

## This script fits all mixed-effects models reported in the paper
## The already-fit models can also be loaded as such:
# load("~/models/ANT-models.RData")

## Data & functions ####
setwd("C:/Users/vipla/Documents/R/Working Drirectory")
ant <- read.table("Verissimo et al/data/ANT-data.csv", sep=",", header=T, fileEncoding="utf-8")

## Functions
# Centering
c.  <- function(x) {scale(as.numeric(x), scale=F)}
# p-values
p.lmer <- function(model){
  p.from.t <- function (tval, obs, nf){2*(1-pt(abs(tval), obs-nf))}
  msum <- summary(model)
  mts <- msum$coefficients[, "t value"]
  mdp <- msum$devcomp$dims["N"]
  mfe <- nrow(msum$coefficients)
  mps <- as.matrix(round(p.from.t(mts, mdp, mfe), 4))
  colnames(mps) <- "p value"
  return(mps)
}

## Linear analysis (main model) ####
## See Results section, Figure 1, and Supplementary Table 2

## Subset (excluded: participants less than 75% accuracy, incorrect responses, double cue trials)
nrow(ant.main <- droplevels(ant[ant$AccuracyAtLeast75==1 & ant$Cue!="DOUBLE" & ant$Accuracy==1,]))
length(unique(ant.main$Participant_ID))
## Contrasts (for Sex, Flanker, and Cue)
ant.main$Sex <- factor(ant.main$Sex, levels=c("Female", "Male"))
(contrasts(ant.main$Sex) <- matrix(dimnames=list(levels(ant.main$Sex), c(".ME")), c(-0.5, 0.5), nrow=2))
ant.main$Flanker <- factor(ant.main$Flanker, levels=c("CONGRUENT", "INCONGRUENT"))
(contrasts(ant.main$Flanker) <- matrix(dimnames=list(levels(ant.main$Flanker), c(".EXECUTIVE")), c(-0.5, 0.5), nrow=2))
ant.main$Cue <- factor(ant.main$Cue, levels=c("CENTRAL", "NOCUE", "SPATIAL"))
(contrasts(ant.main$Cue) <- matrix(dimnames=list(levels(ant.main$Cue), c(".ALERTING", ".ORIENTING")), c(-1/3, 2/3, -1/3, 1/3, 1/3, -2/3), nrow=3))

## Model
library(lme4)
m <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
            Cue * Flanker +
            Cue * c.(Age) +
            Flanker * c.(Age) +
            Sex + c.(Education) + c.(Trial),
          ant.main)
print(summary(m), cor=F)
cbind(summary(m)$coefficients, p.lmer(m))
confint(m, parm="beta_")

## Linear analysis (follow-ups) ####

## Follow-ups at different ages (see Results section)
## Baseline: Age=58
m.age.min <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                    Cue * Flanker +
                    Cue * I(Age-58) +
                    Flanker * I(Age-58) +
                    Sex + c.(Education) + c.(Trial),
                  ant.main)
print(summary(m.age.min), cor=F)
cbind(summary(m.age.min)$coefficients, p.lmer(m.age.min))
confint(m.age.min, parm=c("Cue.ALERTING", "Cue.ORIENTING", "Flanker.EXECUTIVE"))

## Baseline: Age=98
m.age.max <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                    Cue * Flanker +
                    Cue * I(Age-98) +
                    Flanker * I(Age-98) +
                    Sex + c.(Education) + c.(Trial),
                  ant.main)
print(summary(m.age.max), cor=F)
cbind(summary(m.age.max)$coefficients, p.lmer(m.age.max))
confint(m.age.max, parm=c("Cue.ALERTING", "Cue.ORIENTING", "Flanker.EXECUTIVE"))

## Baseline: Age=90
m.age.90 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                   Cue * Flanker +
                   Cue * I(Age-90) +
                   Flanker * I(Age-90) +
                   Sex + c.(Education) + c.(Trial),
                 ant.main)
print(summary(m.age.90), cor=F)
cbind(summary(m.age.90)$coefficients, p.lmer(m.age.90))
confint(m.age.90, parm=c("Cue.ALERTING", "Cue.ORIENTING", "Flanker.EXECUTIVE"))

## Baseline: Age=74 (smallest alerting effect, in absolute value of the estimate)
m.age.74 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                   Cue * Flanker +
                   Cue * I(Age-74) +
                   Flanker * I(Age-74) +
                   Sex + c.(Education) + c.(Trial),
                 ant.main)
print(summary(m.age.74), cor=F)
cbind(summary(m.age.74)$coefficients, p.lmer(m.age.74))
confint(m.age.74, parm="Cue.ALERTING")

## Follow-ups for different flankers (see caption of Supplementary Table 2 and caption of Supplementary Figure 1)
## Baseline: Flanker INCONGRUENT
ant.main$Flanker <- factor(ant.main$Flanker, c("INCONGRUENT", "CONGRUENT"))
contrasts(ant.main$Flanker)
m.incong <- update(m)
print(summary(m.incong), cor=F)
cbind(summary(m.incong)$coefficients, p.lmer(m.incong))
confint(m.incong, parm=c("c.(Age)", "Cue.ALERTING", "Cue.ORIENTING"))

## Baseline: Flanker CONGRUENT
ant.main$Flanker <- factor(ant.main$Flanker, c("CONGRUENT", "INCONGRUENT"))
contrasts(ant.main$Flanker)
m.cong <- update(m)
print(summary(m.cong), cor=F)
cbind(summary(m.cong)$coefficients, p.lmer(m.cong))
confint(m.cong, parm=c("c.(Age)", "Cue.ALERTING", "Cue.ORIENTING"))

## Revert contrasts for Flanker
ant.main$Flanker <- factor(ant.main$Flanker, levels=c("CONGRUENT", "INCONGRUENT"))
(contrasts(ant.main$Flanker) <- matrix(dimnames=list(levels(ant.main$Flanker), c(".EXECUTIVE")), c(-0.5, 0.5), nrow=2))

## Follow-ups for different cues (see caption of Supplementary Table 2 and caption of Supplementary Figure 1)
## Baseline: Cue CENTRAL
ant.main$Cue <- factor(ant.main$Cue, c("CENTRAL", "NOCUE", "SPATIAL"))
contrasts(ant.main$Cue)
m.central <- update(m)
print(summary(m.central), cor=F)
cbind(summary(m.central)$coefficients, p.lmer(m.central))
confint(m.central, parm=c("c.(Age)", "Flanker.EXECUTIVE"))

## Baseline: Cue NOCUE
ant.main$Cue <- factor(ant.main$Cue, c("NOCUE", "CENTRAL", "SPATIAL"))
contrasts(ant.main$Cue)
m.nocue <- update(m)
print(summary(m.nocue), cor=F)
cbind(summary(m.nocue)$coefficients, p.lmer(m.nocue))
confint(m.nocue, parm=c("c.(Age)", "Flanker.EXECUTIVE"))

## Baseline: Cue SPATIAL
ant.main$Cue <- factor(ant.main$Cue, c("SPATIAL", "CENTRAL", "NOCUE"))
contrasts(ant.main$Cue)
m.spatial <- update(m)
print(summary(m.spatial), cor=F)
cbind(summary(m.spatial)$coefficients, p.lmer(m.spatial))
confint(m.spatial, parm=c("c.(Age)", "Flanker.EXECUTIVE"))

## Revert contrasts for Cue
ant.main$Cue <- factor(ant.main$Cue, levels=c("CENTRAL", "NOCUE", "SPATIAL"))
(contrasts(ant.main$Cue) <- matrix(dimnames=list(levels(ant.main$Cue), c(".ALERTING", ".ORIENTING")), c(-1/3, 2/3, -1/3, 1/3, 1/3, -2/3), nrow=3))

## Linear analysis (sensitivity analyses) ####
## See Supplementary Table 4

## 1. Participant mean log-RT also as covariate
## Compute participant mean log-RT
ant.main$Participant.RT.log <- ave(log(ant.main$RT), list(ant.main$Participant_ID), FUN=mean)

## Model (is singular, because random intercept is now 0)
m.rt <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
               Cue * Flanker +
               Cue * c.(Age) +
               Flanker * c.(Age) +
               Sex + c.(Education) + c.(Trial) + c.(Participant.RT.log),
             ant.main)
print(summary(m.rt), cor=F)
cbind(summary(m.rt)$coefficients, p.lmer(m.rt))
confint(m.rt, parm=c("Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 2. Raw (unlogged) RTs as dependent measure
m.raw <- lmer(RT ~ (1+Flanker|Participant_ID) +
                Cue * Flanker +
                Cue * c.(Age) +
                Flanker * c.(Age) +
                Sex + c.(Education) + c.(Trial),
              ant.main)
print(summary(m.raw), cor=F)
cbind(summary(m.raw)$coefficients, p.lmer(m.raw))
confint(m.raw, parm=c("Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 3. No covariates
## (does not converge with Flanker slope)
m.nc <- lmer(log(RT) ~ (1|Participant_ID) +
               Cue * Flanker +
               Cue * c.(Age) +
               Flanker * c.(Age),
             ant.main)
print(summary(m.nc), cor=F)
cbind(summary(m.nc)$coefficients, p.lmer(m.nc))
confint(m.nc, parm=c("Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 4. Visual acuity also as covariate
## Data points and n of participants without valid visual acuity
nrow(ant.main[!is.na(ant.main$Visual_Acuity),])
length(unique(ant.main[!is.na(ant.main$Visual_Acuity),]$Participant_ID))

## Model
m.vis <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                Cue * Flanker +
                Cue * c.(Age) +
                Flanker * c.(Age) +
                Sex + c.(Education) + c.(Trial) + c.(Visual_Acuity),
              ant.main[!is.na(ant.main$Visual_Acuity),])
print(summary(m.vis), cor=F)
cbind(summary(m.vis)$coefficients, p.lmer(m.vis))
confint(m.vis, parm=c("c.(Visual_Acuity)", "Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 5.1. Only responses faster than 1600ms
## Data points below 1600ms
nrow(ant.main[ant.main$RT<1600,])

## Model
m.1600 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                 Cue * Flanker +
                 Cue * c.(Age) +
                 Flanker * c.(Age) +
                 Sex + c.(Education) + c.(Trial),
               ant.main[ant.main$RT<1600,])
print(summary(m.1600), cor=F)
cbind(summary(m.1600)$coefficients, p.lmer(m.1600))
confint(m.1600, parm=c("Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 5.2. Only responses faster than 1500ms
## Data points below 1500ms
nrow(ant.main[ant.main$RT<1500,])

## Model
m.1500 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                 Cue * Flanker +
                 Cue * c.(Age) +
                 Flanker * c.(Age) +
                 Sex + c.(Education) + c.(Trial),
               ant.main[ant.main$RT<1500,])
print(summary(m.1500), cor=F)
cbind(summary(m.1500)$coefficients, p.lmer(m.1500))
confint(m.1500, parm=c("Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 6. Both central and double cues included
## Subset (excluded: participants less than 75% accuracy, incorrect responses - double cue trials are kept)
## Data points and n of participants with central and double cues
nrow(ant.double <- droplevels(ant[ant$AccuracyAtLeast75==1 & ant$Accuracy==1,]))
length(unique(ant.double$Participant_ID))
## Recode DOUBLE are CENTRAL
levels(ant.double$Cue)[levels(ant.double$Cue)=="DOUBLE"] <- "CENTRAL"
## Contrasts (for Sex, Flanker, and Cue)
ant.double$Sex <- factor(ant.double$Sex, levels=c("Female", "Male"))
(contrasts(ant.double$Sex) <- matrix(dimnames=list(levels(ant.double$Sex), c(".ME")), c(-0.5, 0.5), nrow=2))
ant.double$Flanker <- factor(ant.double$Flanker, levels=c("CONGRUENT", "INCONGRUENT"))
(contrasts(ant.double$Flanker) <- matrix(dimnames=list(levels(ant.double$Flanker), c(".EXECUTIVE")), c(-0.5, 0.5), nrow=2))
ant.double$Cue <- factor(ant.double$Cue, levels=c("CENTRAL", "NOCUE", "SPATIAL"))
(contrasts(ant.double$Cue) <- matrix(dimnames=list(levels(ant.double$Cue), c(".ALERTING", ".ORIENTING")), c(-1/3, 2/3, -1/3, 1/3, 1/3, -2/3), nrow=3))

## Model
## Alerting and executive can be estimated for both double and central cues, but orienting should not be estimated in this model
m.double <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                   Cue * Flanker +
                   Cue * c.(Age) +
                   Flanker * c.(Age) +
                   Sex + c.(Education) + c.(Trial),
                 ant.double)
print(summary(m.double), cor=F)
cbind(summary(m.double)$coefficients, p.lmer(m.double))
confint(m.double, parm=c("Cue.ALERTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 7. Lenient participant exclusion: at least 50% accuracy
## Subset (excluded: incorrect responses, double cue trials - participants with less than 75% accuracy, but >=50%, are kept)
## Data points and n of participants with at least 50% accuracy
nrow(ant.acc50 <- droplevels(ant[ant$Cue!="DOUBLE" & ant$Accuracy==1,]))
length(unique(ant.acc50$Participant_ID))
## Contrasts (for Sex, Flanker, and Cue)
ant.acc50$Sex <- factor(ant.acc50$Sex, levels=c("Female", "Male"))
(contrasts(ant.acc50$Sex) <- matrix(dimnames=list(levels(ant.acc50$Sex), c(".ME")), c(-0.5, 0.5), nrow=2))
ant.acc50$Flanker <- factor(ant.acc50$Flanker, levels=c("CONGRUENT", "INCONGRUENT"))
(contrasts(ant.acc50$Flanker) <- matrix(dimnames=list(levels(ant.acc50$Flanker), c(".EXECUTIVE")), c(-0.5, 0.5), nrow=2))
ant.acc50$Cue <- factor(ant.acc50$Cue, levels=c("CENTRAL", "NOCUE", "SPATIAL"))
(contrasts(ant.acc50$Cue) <- matrix(dimnames=list(levels(ant.acc50$Cue), c(".ALERTING", ".ORIENTING")), c(-1/3, 2/3, -1/3, 1/3, 1/3, -2/3), nrow=3))

## Model
m.acc50 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                  Cue * Flanker +
                  Cue * c.(Age) +
                  Flanker * c.(Age) +
                  Sex + c.(Education) + c.(Trial),
                ant.acc50)
print(summary(m.acc50), cor=F)
cbind(summary(m.acc50)$coefficients, p.lmer(m.acc50))
confint(m.acc50, parm=c("Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 8. Stringent participant exclusion: at least 90% accuracy
## Data points and n of participants with at least 90% accuracy
nrow(ant.main[ant.main$AccuracyAtLeast90==1,])
length(unique(ant.main[ant.main$AccuracyAtLeast90==1,]$Participant_ID))

## Model
m.acc90 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                  Cue * Flanker +
                  Cue * c.(Age) +
                  Flanker * c.(Age) +
                  Sex + c.(Education) + c.(Trial),
                ant.main[ant.main$AccuracyAtLeast90==1,])
print(summary(m.acc90), cor=F)
cbind(summary(m.acc90)$coefficients, p.lmer(m.acc90))
confint(m.acc90, parm=c("Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## 9. Extremely old participants excluded
## Data points and n of participants with age<90
nrow(ant.main[ant.main$Age<90,])
length(unique(ant.main[ant.main$Age<90,]$Participant_ID))

## Model
m.neo <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                Cue * Flanker +
                Cue * c.(Age) +
                Flanker * c.(Age) +
                Sex + c.(Education) + c.(Trial),
              ant.main[ant.main$Age<90,])
print(summary(m.neo), cor=F)
cbind(summary(m.neo)$coefficients, p.lmer(m.neo))
confint(m.neo, parm=c("Cue.ALERTING:c.(Age)", "Cue.ORIENTING:c.(Age)", "Flanker.EXECUTIVE:c.(Age)"))

## Nonlinear analysis: Quadratic model ####
## See Results section and Figure 2a

m.quad <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                 Cue * Flanker +
                 Cue * poly(c.(Age), 2) +
                 Flanker * poly(c.(Age), 2) +
                 Sex + c.(Education) + c.(Trial),
               ant.main)
print(summary(m.quad), cor=F)
cbind(summary(m.quad)$coefficients, p.lmer(m.quad))
confint(m.quad, parm=c("Cue.ALERTING:poly(c.(Age), 2)2", "Cue.ORIENTING:poly(c.(Age), 2)2", "Flanker.EXECUTIVE:poly(c.(Age), 2)2"))

## Nonlinear analysis: Breakpoint model ####

## Discovery of best breakpoint
## See Supplementary Figure 4
## (no random slope due to non-convergence; REML=F more appropriate for model comparisons)
models.fit <- NULL
for(bp in seq(min(ant.main$Age), max(ant.main$Age), 1)){
  ant.main$adj.Age <- ant.main$Age - bp
  ant.main$post.bp <- as.numeric(ant.main$adj.Age>0)
  m.current <- lmer(log(RT) ~ (1|Participant_ID) +
                      Cue * Flanker +
                      Cue * (adj.Age + adj.Age:post.bp) +
                      Flanker * (adj.Age + adj.Age:post.bp) +
                      Sex + c.(Education) + c.(Trial),
                    ant.main, REML=F)
  m.current.fit <- data.frame(Breakpoint=bp, Deviance=deviance(m.current, REML=F), AIC=AIC(m.current))
  models.fit <- rbind(models.fit, m.current.fit)}
## Best breakpoint
models.fit[models.fit$AIC==min(models.fit$AIC),]
bp <- models.fit[models.fit$AIC==min(models.fit$AIC),]$Breakpoint

## Make indicator variable
ant.main$adj.Age <- ant.main$Age - bp
ant.main$post.bp <- as.numeric(ant.main$adj.Age>0)

## Model
## See Results section and Figure 2b
m.break <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                  Cue * Flanker +
                  Cue * (adj.Age + adj.Age:post.bp) +
                  Flanker * (adj.Age + adj.Age:post.bp) +
                  Sex + c.(Education) + c.(Trial),
                ant.main)
print(summary(m.break), cor=F)
cbind(summary(m.break)$coefficients, p.lmer(m.break))
confint(m.break, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp",
                        "Flanker.EXECUTIVE:adj.Age"))

## Nonlinear analysis: Breakpoint model (follow-ups) ####
## See Results section

## Follow-ups at different ages
## Baseline: Age=58
m.break.age.min <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                          Cue * Flanker +
                          Cue * (I(Age-58)+ adj.Age:post.bp) +
                          Flanker * (I(Age-58) + adj.Age:post.bp) +
                          Sex + c.(Education) + c.(Trial),
                        ant.main)
print(summary(m.break.age.min), cor=F)
cbind(summary(m.break.age.min)$coefficients, p.lmer(m.break.age.min))
confint(m.break.age.min, parm=c("Flanker.EXECUTIVE"))

## Baseline: Age=76
m.break.age.76 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                         Cue * Flanker +
                         Cue * (I(Age-76)+ adj.Age:post.bp) +
                         Flanker * (I(Age-76) + adj.Age:post.bp) +
                         Sex + c.(Education) + c.(Trial),
                       ant.main)
print(summary(m.break.age.76), cor=F)
cbind(summary(m.break.age.76)$coefficients, p.lmer(m.break.age.76))
confint(m.break.age.76, parm=c("Flanker.EXECUTIVE"))

## Follow-up for age slope greater than 76
## Make reversed indicator variable
ant.main$adj.Age <- ant.main$Age - bp
ant.main$post.bp <- as.numeric(ant.main$adj.Age<=0)
## Model
m.break.rev <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                      Cue * Flanker +
                      Cue * (adj.Age + adj.Age:post.bp) +
                      Flanker * (adj.Age + adj.Age:post.bp) +
                      Sex + c.(Education) + c.(Trial),
                    ant.main)
print(summary(m.break.rev), cor=F)
cbind(summary(m.break.rev)$coefficients, p.lmer(m.break.rev))
confint(m.break.rev, parm=c("Flanker.EXECUTIVE:adj.Age"))

## Nonlinear analysis: Breakpoint model (sensitivity analyses) ####
## See Supplementary Table 5

## Make indicator variable
ant.main$adj.Age <- ant.main$Age - bp
ant.main$post.bp <- as.numeric(ant.main$adj.Age>0)

## 1. Participant mean log-RT also as covariate
## Model (is singular, because random intercept is now 0)
m.break.rt <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                     Cue * Flanker +
                     Cue * (adj.Age + adj.Age:post.bp) +
                     Flanker * (adj.Age + adj.Age:post.bp) +
                     Sex + c.(Education) + c.(Trial) + c.(Participant.RT.log),
                   ant.main)
print(summary(m.break.rt), cor=F)
cbind(summary(m.break.rt)$coefficients, p.lmer(m.break.rt))
confint(m.break.rt, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 2. Raw (unlogged) RTs as dependent measure
m.break.raw <- lmer(RT ~ (1+Flanker|Participant_ID) +
                      Cue * Flanker +
                      Cue * (adj.Age + adj.Age:post.bp) +
                      Flanker * (adj.Age + adj.Age:post.bp) +
                      Sex + c.(Education) + c.(Trial),
                    ant.main)
print(summary(m.break.raw), cor=F)
cbind(summary(m.break.raw)$coefficients, p.lmer(m.break.raw))
confint(m.break.raw, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 3. No covariates
m.break.nc <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                     Cue * Flanker +
                     Cue * (adj.Age + adj.Age:post.bp) +
                     Flanker * (adj.Age + adj.Age:post.bp),
                   ant.main)
print(summary(m.break.nc), cor=F)
cbind(summary(m.break.nc)$coefficients, p.lmer(m.break.nc))
confint(m.break.nc, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 4. Visual acuity also as covariate
m.break.vis <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                      Cue * Flanker +
                      Cue * (adj.Age + adj.Age:post.bp) +
                      Flanker * (adj.Age + adj.Age:post.bp) +
                      Sex + c.(Education) + c.(Trial) + c.(Visual_Acuity),
                    ant.main[!is.na(ant.main$Visual_Acuity),])
print(summary(m.break.vis), cor=F)
cbind(summary(m.break.vis)$coefficients, p.lmer(m.break.vis))
confint(m.break.vis, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 5.1. Only responses faster than 1600ms
## (does not converge with Flanker slope)
m.break.1600 <- lmer(log(RT) ~ (1|Participant_ID) +
                       Cue * Flanker +
                       Cue * (adj.Age + adj.Age:post.bp) +
                       Flanker * (adj.Age + adj.Age:post.bp) +
                       Sex + c.(Education) + c.(Trial),
                     ant.main[ant.main$RT<1600,])
print(summary(m.break.1600), cor=F)
cbind(summary(m.break.1600)$coefficients, p.lmer(m.break.1600))
confint(m.break.1600, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 5.2. Only responses faster than 1500ms
m.break.1500 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                       Cue * Flanker +
                       Cue * (adj.Age + adj.Age:post.bp) +
                       Flanker * (adj.Age + adj.Age:post.bp) +
                       Sex + c.(Education) + c.(Trial),
                     ant.main[ant.main$RT<1500,])
print(summary(m.break.1500), cor=F)
cbind(summary(m.break.1500)$coefficients, p.lmer(m.break.1500))
confint(m.break.1500, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 6. Both central and double cues included
## Subset (excluded: participants less than 75% accuracy, incorrect responses - double cue trials are kept)
## Make indicator variable
ant.double$adj.Age <- ant.double$Age - bp
ant.double$post.bp <- as.numeric(ant.double$adj.Age>0)

## Model
## Alerting and executive can be estimated for both double and central cues, but orienting should not be estimated in this model
m.break.double <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                         Cue * Flanker +
                         Cue * (adj.Age + adj.Age:post.bp) +
                         Flanker * (adj.Age + adj.Age:post.bp) +
                         Sex + c.(Education) + c.(Trial),
                       ant.double)
print(summary(m.break.double), cor=F)
cbind(summary(m.break.double)$coefficients, p.lmer(m.break.double))
confint(m.break.double, parm=c("Cue.ALERTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 7. Lenient participant exclusion: at least 50% accuracy
## Subset (excluded: incorrect responses, double cue trials - participants with less than 75% accuracy, but >=50%, are kept)
## Make indicator variable
ant.acc50$adj.Age <- ant.acc50$Age - bp
ant.acc50$post.bp <- as.numeric(ant.acc50$adj.Age>0)

## Model
m.break.acc50 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                        Cue * Flanker +
                        Cue * (adj.Age + adj.Age:post.bp) +
                        Flanker * (adj.Age + adj.Age:post.bp) +
                        Sex + c.(Education) + c.(Trial),
                      ant.acc50)
print(summary(m.break.acc50), cor=F)
cbind(summary(m.break.acc50)$coefficients, p.lmer(m.break.acc50))
confint(m.break.acc50, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 8. Stringent participant exclusion: at least 90% accuracy
m.break.acc90 <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                        Cue * Flanker +
                        Cue * (adj.Age + adj.Age:post.bp) +
                        Flanker * (adj.Age + adj.Age:post.bp) +
                        Sex + c.(Education) + c.(Trial),
                      ant.main[ant.main$AccuracyAtLeast90==1,])
print(summary(m.break.acc90), cor=F)
cbind(summary(m.break.acc90)$coefficients, p.lmer(m.break.acc90))
confint(m.break.acc90, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## 9. Extremely old participants excluded
m.break.neo <- lmer(log(RT) ~ (1+Flanker|Participant_ID) +
                      Cue * Flanker +
                      Cue * (adj.Age + adj.Age:post.bp) +
                      Flanker * (adj.Age + adj.Age:post.bp) +
                      Sex + c.(Education) + c.(Trial),
                    ant.main[ant.main$Age<90,])
print(summary(m.break.neo), cor=F)
cbind(summary(m.break.neo)$coefficients, p.lmer(m.break.neo))
confint(m.break.neo, parm=c("Cue.ALERTING:adj.Age:post.bp", "Cue.ORIENTING:adj.Age:post.bp", "Flanker.EXECUTIVE:adj.Age:post.bp"))

## Accuracy analysis (main model) ####
## See Supplementary Table 6

## Subset (excluded: participants less than 75% accuracy, double cue trials - incorrect responses are kept)
nrow(ant.incorrect <- droplevels(ant[ant$AccuracyAtLeast75==1 & ant$Cue!="DOUBLE",]))
length(unique(ant.incorrect$Participant_ID))
## Contrasts (for Sex, Flanker, and Cue)
ant.incorrect$Sex <- factor(ant.incorrect$Sex, levels=c("Female", "Male"))
(contrasts(ant.incorrect$Sex) <- matrix(dimnames=list(levels(ant.incorrect$Sex), c(".ME")), c(-0.5, 0.5), nrow=2))
ant.incorrect$Flanker <- factor(ant.incorrect$Flanker, levels=c("CONGRUENT", "INCONGRUENT"))
(contrasts(ant.incorrect$Flanker) <- matrix(dimnames=list(levels(ant.incorrect$Flanker), c(".EXECUTIVE")), c(-0.5, 0.5), nrow=2))
ant.incorrect$Cue <- factor(ant.incorrect$Cue, levels=c("CENTRAL", "NOCUE", "SPATIAL"))
(contrasts(ant.incorrect$Cue) <- matrix(dimnames=list(levels(ant.incorrect$Cue), c(".ALERTING", ".ORIENTING")), c(-1/3, 2/3, -1/3, 1/3, 1/3, -2/3), nrow=3))

## Model
## (does not converge with Flanker slope)
library(blme)
m.accuracy <- bglmer(Accuracy ~ (1|Participant_ID) +
                       Cue * Flanker +
                       Cue * c.(Age) +
                       Flanker * c.(Age) +
                       Sex + c.(Education) + c.(Trial),
                     ant.incorrect,
                     family="binomial",
                     fixef.prior=normal(c(3, rep(0.4, 12))),
                     control=glmerControl(optimizer="bobyqa"))
print(summary(m.accuracy), cor=F)
confint(m.accuracy, parm="beta_", method="Wald")

## Accuracy analysis (follow-ups) ####
## See caption of Supplementary Table 6

## Follow-ups at different ages (see Results section)
## Baseline: Age=58
m.accuracy.age.min <- bglmer(Accuracy ~ (1|Participant_ID) +
                               Cue * Flanker +
                               Cue * I(Age-58) +
                               Flanker * I(Age-58) +
                               Sex + c.(Education) + c.(Trial),
                             ant.incorrect, family="binomial",
                             fixef.prior=normal(c(3, rep(0.4, 12))),
                             control=glmerControl(optimizer="bobyqa"))
print(summary(m.accuracy.age.min), cor=F)
confint(m.accuracy.age.min, parm="beta_", method="Wald")

## Baseline: Age=98
m.accuracy.age.max <- bglmer(Accuracy ~ (1|Participant_ID) +
                               Cue * Flanker +
                               Cue * I(Age-98) +
                               Flanker * I(Age-98) +
                               Sex + c.(Education) + c.(Trial),
                             ant.incorrect, family="binomial",
                             fixef.prior=normal(c(3, rep(0.4, 12))),
                             control=glmerControl(optimizer="bobyqa"))
print(summary(m.accuracy.age.max), cor=F)
confint(m.accuracy.age.max, parm="beta_", method="Wald")

## Baseline: Age=90
m.accuracy.age.90 <- bglmer(Accuracy ~ (1|Participant_ID) +
                              Cue * Flanker +
                              Cue * I(Age-90) +
                              Flanker * I(Age-90) +
                              Sex + c.(Education) + c.(Trial),
                            ant.incorrect, family="binomial",
                            fixef.prior=normal(c(3, rep(0.4, 12))),
                            control=glmerControl(optimizer="bobyqa"))
print(summary(m.accuracy.age.90), cor=F)
confint(m.accuracy.age.90, parm="beta_", method="Wald")
