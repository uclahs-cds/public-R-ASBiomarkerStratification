
library(BoutrosLab.plotting.survival);
library(BoutrosLab.statistics.survival)


source('data.frame.confusion.matrix.R');

a <- data.frame(
  biodb$BiopsyUpgraded,
  biodb$PHI,
  biodb$DaysBxToLastReview,
  biodb$DaysDxToUpgrade
  );
colnames(a) <- c('actual', 'predicted', 'DaysBxToLastReview', 'time');

for (i in 1:dim(a)[1]){
   if(is.na(a$time[i])){
      a$time[i] <- a$DaysBxToLastReview[i]
     }

  }

threshold <- 44.70000;

a <- a[complete.cases(a),];

b <- data.frame.confusion.matrix(a[,c('actual', 'predicted')], threshold)

c <- BoutrosLab.plotting.survival::create.km.plot(
  survival.object = Surv(a$time, a$actual), 
  patient.groups = as.factor(b$predicted),
  statistical.method = "logrank",
  #ph.assumption.check = "warning.and.plot"
  xaxis.cex = 1.4,
  yaxis.cex = 1.4,
  xlab.label = 'Time since diagnosis (Days)',
  xlab.cex = 1.5,
  ylab.label = 'Estimated survival probability',
  ylab.cex = 1.5,
  key.groups.corner = c(0,0), 
  key.groups.x.pos = 0, 
  key.groups.y.pos = 0.01, 
  key.groups.cex = 1.3,
  key.stats.corner = c(1,0),
  key.stats.x.pos = 1, 
  key.stats.y.pos = 0.9,
  risk.label.pos = -860,
  risktable.fontsize = 15,
  left.padding = 7
);

## # fit.coxmodel
#tmp_model <- coxph(
#                Surv(time_to_bcr, bcr) ~ meth,
#                data = data_df
#                )

d <-BoutrosLab.statistics.survival::fit.coxmodel(
  groups = as.factor(b$predicted),
  survobj = Surv(a$time, a$actual),
  #other.data = data.frame(gender),
  #stratification.factor = numbered.age,
  #stratification.value = 60,
  return.cox.model = TRUE
  )
