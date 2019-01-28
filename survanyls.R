data(GBSG2,package = 'TH.data')
library('survival')
library('survminer')

# Count censored and uncensored data
num_cens <- table(GBSG2$cens)
num_cens

# Create barplot of censored and uncensored data
barplot(num_cens)

# Use help() to look at cens
help(GBSG2, package = "TH.data")

# Create Surv-Object
sobj <- Surv(GBSG2$time, GBSG2$cens)

# Look at 10 first elements
sobj[1:10]

# Look at summary
summary(sobj)

# Look at structure
str(sobj)

# Kaplan-Meier estimate
GBSG2
km <- survfit(Surv(time,cens)~1, data = GBSG2)

# plot of the Kaplan-Meier estimate
ggsurvplot(km)

# add the risk table to plot
ggsurvplot(km, risk.table = TRUE)

# add a line showing the median survival time
ggsurvplot(km, risk.table = TRUE, surv.median.line = "hv")

# Weibull model
wb <- survreg(Surv(time, cens) ~ 1, data = GBSG2)

# Compute the median survival from the model
predict(wb, type = "quantile", p = 0.5 , newdata = data.frame(1))

# 70 Percent of patients survive beyond time point...
predict(wb, type ='quantile', p = 1-0.7, newdata = data.frame(1))


# Retrieve survival curve from model probabilities 
surv <- seq(.99, .01, by = -.01)

# Get time for each probability
t <- predict(wb, type = 'quantile', p =1-surv , newdata = data.frame(1))

# Create data frame with the information
surv_wb <- data.frame(time = t, surv = surv)

# Look at first few lines of the result
head(surv_wb)

surv_wb <- data.frame(time =t, surv = surv, 
                      upper = NA, lower = NA, std.err = NA)

# Plot
ggsurvplot_df(fit = surv_wb, surv.geom = geom_line)

# Weibull model
wbmod <- survreg(Surv(time,cens) ~ horTh, data =GBSG2)

# Retrieve survival curve from model
surv <- seq(.99, .01, by = -.01)
t_yes <- predict(wbmod, type = "quantile", p = 1-surv,
                 newdata = data.frame(horTh='yes'))

# Take a look at survival curve
str(t_yes)

# Lung Cancer data
 data('lung')
# Look at the data set
str(lung)

# Estimate a Weibull model
wbmod <- survreg(Surv(time, status) ~ sex, data = lung)
coef(wbmod)
#WE can clearly see that, The sexfemale coefficient 
#is positive which means women tend to survive longer.

##Weibull visualization data=GBSG2
# 1.Weibull model
wbmod <- survreg(Surv(time,cens) ~ horTh + tsize, data = GBSG2)

# 2.Imaginary patients
newdat <- expand.grid(
  horTh = levels(GBSG2$horTh),
  tsize = quantile(GBSG2$tsize, probs = c(0.25, 0.50, 0.75)))
newdat

# 3.Compute survival curves
surv <- seq(.99, .01, by = -.01)
t <- predict(wbmod, type = 'quantile', p = 1-surv,newdat = newdat)

# How many rows and columns does t have?
dim(t)
# 4.Use cbind() to combine the information in newdat with t
surv_wbmod_wide <- cbind(newdat, t)

# 4.Use melt() to bring the data.frame to long format
library('reshape2')
surv_wbmod <- melt(
  surv_wbmod_wide, id.vars = c("horTh", "tsize"), 
  variable.name = "surv_id", value.name = "time")

# 4.Use surv_wbmod$surv_id to add the correct survival probabilities surv
surv_wbmod$surv <- surv[as.numeric(surv_wbmod$surv_id)]
# Add columns upper, lower, std.err, and strata to the data.frame
surv_wbmod[, c("upper", "lower", "std.err", "strata")] <- NA

# Take a look at the structure of the object
str(surv_wbmod)

#5. Plot the survival curves
ggsurvplot_df(surv_wbmod, surv.geom = geom_line,
              linetype = 'horTh', color = 'tsize', legend.title = NULL)

# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2)

# Log-Normal model
lnmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2, dist = "lognormal")

# Newdata
newdat <- data.frame(horTh = levels(GBSG2$horTh))

# Surv
surv <- seq(.99, .01, by = -.01)

# Survival curve from Weibull model and log-normal model
wbt <- predict(wbmod, type = "quantile", p = 1-surv, newdata = newdat)
lnt <- predict(lnmod, type = 'quantile', p = 1-surv, newdata = newdat)
surv_wide = cbind(wbt,lnt)
# Melt the data.frame into long format.
surv_long <- melt(surv_wide, id.vars = c("horTh", "dist"),
                  variable.name = "surv_id", value.name = "time")

# Add column for the survival probabilities
##isme ek error h
surv_long$surv <- surv[as.numeric(surv_long$surv_id)]

# Add columns upper, lower, std.err, and strata contianing NA values
surv_long[, c("upper", "lower", "std.err", "strata")] <- NA

# Plot the survival curves
ggsurvplot_df(surv_long, surv.geom = geom_line,linetype = 'horTh',
              color = 'dist', legend.title = NULL)


# Cox model
cxmod <- coxph(Surv(time, cens) ~ horTh + tsize, data = GBSG2)

# Imaginary patients
newdat <- expand.grid(
  horTh = levels(GBSG2$horTh),
  tsize = quantile(GBSG2$tsize, probs = c(0.25, 0.5, 0.75)))
rownames(newdat) <- letters[1:6]

# Inspect newdat
newdat

# Compute survival curves
cxsf <- survfit(cxmod, data = GBSG2, newdata = newdat, conf.type = "none")
# Look at first 6 rows of cxsf$surv and time points
head(cxsf$surv)
head(cxsf$time)

# Compute data.frame needed for plotting
surv_cxmod0 <- surv_summary(cxsf)

# Get a character vector of patient letters (patient IDs)
pid <- as.character(surv_cxmod0$strata)

# Multiple of the rows in newdat so that it fits with surv_cxmod0
m_newdat <- newdat[pid, ]

# Add patient info to data.frame
surv_cxmod <- cbind(surv_cxmod0, m_newdat)

# Plot
ggsurvplot_df(surv_cxmod, linetype = 'horTh', color = 'tsize',
              legend.title = NULL, censor = FALSE)

##one more cox
# Compute Cox model and survival curves
cxmod <- coxph(Surv(time, status) ~ performance, data = lung)
new_lung <- data.frame(performance = c(60, 70, 80, 90))
cxsf <- survfit(cxmod, data = lung, newdata = new_lung, conf.type = "none")

# Use the summary of cxsf to take a vector of patient IDs
surv_cxmod0 <- surv_summary(cxsf)
pid <- as.character(surv_cxmod0$strata)

# Duplicate rows in newdat to fit with surv_cxmod0 and add them in
m_newdat <- new_lung[pid, , drop = FALSE]
surv_cxmod <- cbind(surv_cxmod0, m_newdat)

# Plot
ggsurvplot_df(surv_cxmod, color = 'performance', legend.title = NULL, censor = FALSE)

##mix
# Compute Kaplan-Meier curve
km <- survfit(Surv(time, status) ~ 1, data = lung)

# Compute Cox model
cxmod <- coxph(Surv(time, status) ~ performance, data = lung)

# Compute Cox model survival curves
new_lung <- data.frame(performance = c(60, 70, 80, 90))
cxsf <- survfit(cxmod, data = lung, newdata = new_lung, conf.type = "none")

# Plot Kaplan-Meier curve
ggsurvplot(km, conf.int = FALSE)

# Plot Cox model survival curves
ggsurvplot(cxsf, censor = FALSE)




