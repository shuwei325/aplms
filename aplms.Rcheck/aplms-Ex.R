pkgname <- "aplms"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "aplms-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('aplms')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AIC.aplms")
### * AIC.aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: AIC.aplms
### Title: Akaike information criterion
### Aliases: AIC.aplms

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
AIC(model)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("AIC.aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("BIC.aplms")
### * BIC.aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BIC.aplms
### Title: Bayesian information criterion
### Aliases: BIC.aplms

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
BIC(model)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BIC.aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("aplms")
### * aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aplms
### Title: Fitting Additive partial linear models with symmetric errors
### Aliases: aplms print.aplms

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
plot(model)
summary(model)
print(model)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("aplms.diag.plot")
### * aplms.diag.plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aplms.diag.plot
### Title: Diagnostic Plots for additive partial linear models with
###   symmetric errors
### Aliases: aplms.diag.plot
### Keywords: Additive Residuals errors linear models partial symmetric
###   with

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
aplms.diag.plot(model, which = 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aplms.diag.plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coef.aplms")
### * coef.aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coef.aplms
### Title: Extract the coefficients of the fitted APLMS model
### Aliases: coef.aplms

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
coef(model)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coef.aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("family.elliptical")
### * family.elliptical

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: family.elliptical
### Title: Family Objects for Elliptical Models
### Aliases: family.elliptical Cauchy Cnormal Gstudent Glogis Normal
###   Powerexp Student LogisI LogisII GNormal
### Keywords: Elliptical distributions

### ** Examples

Normal()
Powerexp(k=0.1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("family.elliptical", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hospitalization")
### * hospitalization

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hospitalization
### Title: Respiratory diseases hospitalization Dataset
### Aliases: hospitalization
### Keywords: datasets

### ** Examples

data(hospitalization)
head(hospitalization)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hospitalization", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("influence.aplms")
### * influence.aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: influence.aplms
### Title: local influence analysis of the object 'aplms()'
### Aliases: influence.aplms

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
influence(model, perturbation = c("case-weight"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("influence.aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("influenceplot.aplms")
### * influenceplot.aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: influenceplot.aplms
### Title: Local influence plots of the object 'aplms()'
### Aliases: influenceplot.aplms

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
influenceplot.aplms(model, perturbation = c("case-weight"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("influenceplot.aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.aplms")
### * plot.aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.aplms
### Title: Default APLMS plotting
### Aliases: plot.aplms
### Keywords: Additive Residuals errors linear partial symmetric with xs

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
plot(model)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("residuals.aplms")
### * residuals.aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: residuals.aplms
### Title: Extract Residuals for APLMS fits
### Aliases: residuals.aplms
### Keywords: Additive Residuals errors linear models partial symmetric
###   with

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
residuals(model)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("residuals.aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.aplms")
### * summary.aplms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.aplms
### Title: Print method for "aplms" class
### Aliases: summary.aplms

### ** Examples

data(temperature)
temperature.df = data.frame(temperature,time=1:length(temperature))
model<-aplms::aplms(temperature ~ 1,
                   npc=c("time"), basis=c("cr"),Knot=c(60),
                   data=temperature.df,family=Powerexp(k=0.3),p=1,
                   control = list(tol = 0.001,
                                  algorithm1 = c("P-GAM"),
                                  algorithm2 = c("BFGS"),
                                  Maxiter1 = 20,
                                  Maxiter2 = 25),
                   lam=c(10))
summary(model)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.aplms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("temperature")
### * temperature

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: temperature
### Title: Global Annual Mean Surface Air Temperature Change
### Aliases: temperature
### Keywords: datasets

### ** Examples

data(temperature)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("temperature", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
