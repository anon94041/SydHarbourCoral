# IPM for Sydney Harbour population of coral Plesiastrea versipora
# accompanies
# How does a widespread reef coral maintain populations in isolated environments? 
# Kristin Precoda, Andrew H. Baird, Alisha Madsen, Toni Mizerek, Brigitte Sommer, Sheena N. Su, Joshua S. Madin

library("xlsx")  # for read.csv
library("Matrix")  # for Matrix()
library("fields")  # for tim.colors()

source("R/functions.R")

#####
# read in various types of data

# growth data from Dec. 2011 - Aug. 2012
growth.data.a <- read.table("data/growth-Dec2011-Aug2012.txt", header=T)
# growth data from Dec. 2012 - Dec. 2013
growth.data.s <- read.table("data/growth-Dec2012-Dec2013.txt", header=T)
names(growth.data.s) <- c("x", "y", "survival")

# size structure data
size.struct <- read.table("data/sizestructure-Dec2011.txt", header=T)

# data on sex & colony sizes; discard ones where size wasn't measured
sex.size.1 <- na.omit(unique(read.table("data/sex-size-Madsenetal2014.txt", header=T)))

# data on sex, reproductivity, & colony size derived from Withers 2000, Fairlight, chapter 3
# Withers 2000 grouped colonies into size classes, so sizes of each
# colony aren't known.  Represent all colonies with the size at the 
# midpoint of their class
sex.size.w <- read.csv("data/sex-size-Withers2000-Fairlight.csv", header=T)
sex.size.w2 <- cbind(sex.size.w, apply(sex.size.w[, 1:2], 1, mean))
names(sex.size.w2)[ncol(sex.size.w2)] <- "size"

# survival data
# assume the probability of seeing a coral in a later year if it's there is
# close to 1 (analogous to saying mark/recapture would have a recapture rate
# of essentially 1)

# survival data from Withers 2000, table 3.4.  Will also add in this
# study's survival data from growth data files (which were read in above)
# set the size to the midpoint of the min and max of each size class
survival.data.w <- read.csv("data/survival-size-Withers2000-Fairlight.csv")
survival.data.w2 <- cbind(survival.data.w,
                          apply(survival.data.w[, 1:2], 1, mean))
names(survival.data.w2)[ncol(survival.data.w2)] <- "sizeRaw"

# end: read in various types of data

#####
# process/transform data

# growth & survival data set
# create single growth & survival data set with both 2011-2 and 2012-3 sizes
# start by multiplicatively extending 2011-2 data to 365 days
# there are no recruits that appeared in 2012 in this data, so no special cases in multiplying
tmp <- data.frame(cbind(growth.data.a$x,
             growth.data.a$x*(growth.data.a$y/growth.data.a$x)^(365/266),
             growth.data.a$survival))
names(tmp) <- c("x", "y", "survival")
# then bind annualized 2011-2 data together with the 2012-3 data
growth.data.r <- rbind(growth.data.s, tmp)
names(growth.data.r) <- c("sizeRaw", "sizeNextRaw", "survival")

x1 <- growth.data.r[complete.cases(growth.data.r), ]$sizeRaw
y1 <- growth.data.r[complete.cases(growth.data.r), ]$sizeNextRaw
# note: bcNLL suggests a transform for both predictor and response
# variable, which in the growth model are both size (Rees, Childs, 
# & Ellner 2014): power suggested is 0.1691593 which is approx 1/6;
# 1/6 is closer to interpretable
optimize(bcNLL, c(0.001, 1), x=x1, y=y1)$minimum

transformtype <- 3

# apply size transformation to growth data
growth.data <- cbind(growth.data.r,
                     transformSize(growth.data.r$sizeRaw, transformtype),
                     transformSize(growth.data.r$sizeNextRaw, transformtype))
names(growth.data) <- c("sizeRaw", "sizeNextRaw", "survival", "size", "sizeNext")

# sex & colony sizes
# examined residuals of logistic regression of sex on colony size, with
# various transforms of size; residuals looked best with untransformed size.
# store transformed size data in sex.size for use in plotting; not used in analyses
sex.size <- cbind(sex.size.1, transformSize(sex.size.1$colony_size, transformtype))
names(sex.size) <- c("date", "sex", "sizeRaw", "size")
sex.size$sex <- -1*(sex.size$sex - 2)      # remap so 1=F, 0=M

# reproductivity & colony sizes
# 0 = not reproductive, 1 = reproductive
# start with Withers 2000 data
repro.w <- rep(0, nrow(sex.size.w2))
repro.w[sex.size.w2$status!="not reproductive"] <- 1
# only have size data for 80 of Madsen's 140 colonies examined for reproductivity, but know 106/140 were reproductive
# Use bootstrap to get sizes for (106-80) reproductive colonies and for (140-106) nonreproductive colonies
set.seed(1234)
bootstrap.repro <- sample(sex.size$sizeRaw, 106-80, replace=T)
set.seed(2345)
bootstrap.nonrepro <- sample(sex.size$sizeRaw, 140-106, replace=T)
repro.size <- data.frame(rbind(cbind(repro.w, sex.size.w2$size),
                               cbind(rep(1, 80), sex.size$sizeRaw),
                               cbind(rep(1, 106-80), bootstrap.repro),
                               cbind(rep(0, 140-106), bootstrap.nonrepro)))
names(repro.size) <- c("reproductive", "size")
# examined residuals of glm(reproductive ~ size, family=binomial, data=repro.size)
# without transformation and with transformations; considerably closer to
# normality using the same 1/6th root transform as in the growth model, and
# also makes modelling a bit easier to have the same transform and not a
# new one
repro.size <- cbind(repro.size, transformSize(repro.size$size, transformtype))
names(repro.size) <- c("reproductive", "sizeRaw", "size")

# survival & colony sizes
# examined residuals of glm(survival ~ size, family=binomial, data=surv.data) with
# and without transformations; residuals are closer to normality using the
# 1/6th root transform (although the residuals for the nonsurvivor points
# are still outliers)
survival.data.w3 <- cbind(survival.data.w2,
                          transformSize(survival.data.w2$sizeRaw, transformtype))
names(survival.data.w3)[ncol(survival.data.w3)] <- "size"

surv.data <- data.frame(rbind(cbind(survival.data.w3$size, survival.data.w3$survival),
                              cbind(growth.data$size, growth.data$survival)))
names(surv.data) <- c("size", "survival")

# end: process/transform data

#####
# build component models

# growth: use lm() to predict next size from current size
# histogram of residuals looks plausibly normal, though both tails are a
# little heavy; this is more evident on the qqplot of the residuals
grow.mod <- lm(sizeNext ~ size, data=growth.data, na.action=na.exclude)
summary(grow.mod)

# logistic regression of sex on size
sex.mod <- glm(sex ~ sizeRaw, family=binomial, data=sex.size)
summary(sex.mod)

# logistic regression to get probability of reproductivity on size
repro.mod <- glm(reproductive ~ size, family=binomial, data=repro.size)
summary(repro.mod)

# logistic regression of survival on size
surv.mod <- glm(survival ~ size, family=binomial,
                data=growth.data[!is.na(growth.data$size), ],
                na.action=na.exclude)
# force the intercept to be 0 because it wasn't significantly different from 0 
# (and thus surv.mod$coeff[1] will be the coefficient for size)
surv.mod <- glm(survival ~ 0 + size, family=binomial,
                data=surv.data[!is.na(surv.data$size), ],
                na.action=na.exclude)
summary(surv.mod)

# end: build component models

#---------------------------------------------
# setup for matrix representing IPM
# use 0.35 cm^2 (= 35 mm^2) as min size and as recruit size; 
# Withers reports about 2/3 of recruits are < 35mm^2 and nearly all the rest are < 50mm^2
minsize <- transformSize(0.35, transformtype)
recruitsize <- transformSize(.35, transformtype)
# max size: greater of (1600 cm^2, 30% larger than largest actual observed size)
maxsize <- transformSize(max(1600, size.struct[, 2]*1.3), transformtype)

sizeindepmort <- 0.01   # size-independent mortality, as a percentage
minreliablesize <- 20   # assume minimum reliably detected size in photos

# parameters for survival and growth models
s.params <- c(surv.mod$coeff[1], sizeindepmort)
g.params <- c(grow.mod$coeff[1], grow.mod$coeff[2], summary(grow.mod)$sigma^2)

# estimate best Q via log likelihood: use density at nearest mesh point, 
# for each observed data point
# also, use large matrix (800x800)

set.seed(1234)
mat.dim <- 800            # matrix dimension (mat.dim x mat.dim)
# y is the set of sizes used as mesh points
y <- minsize + c(0:(mat.dim-1)) * (maxsize - minsize) / (mat.dim-1)
qseqlength <- 601         # how many Q values to calculate LL for
qseq <- seq(100e-4, 700e-4, length.out=qseqlength)  # calculate LL over a range of Q values
ll <- rep(0, qseqlength)
# observed sizes and the closest mesh points to the observed sizes
ssd.obs <- transformSize(size.struct[size.struct$area_T > minreliablesize, 2], transformtype)
closestys <- rep(0, length(ssd.obs))
for (i in 1:length(ssd.obs))
  closestys[i] <- which.min(abs(y-ssd.obs[i]))
for (j in 1:qseqlength) {
  # reproduction parameters: includes Q value to get LL for
  r.params <- c(sex.mod$coeff[1], sex.mod$coeff[2],
                repro.mod$coeff[1], repro.mod$coeff[2],
                qseq[j], recruitsize, transformtype)
  
  M <- bigmatrix(y, s.params, r.params, g.params)
  mat <- as.matrix(Matrix(M$G) %*% Matrix(M$S) + Matrix(M$R))
  # in M$G, the second index is the current size; first index is next size
  Pvers.eigen <- get.eigen.stuff(mat)
  ssd <- as.vector(Pvers.eigen[[2]])   # stable size distribution
  normalizer <- sum(ssd[y > transformSize(minreliablesize, transformtype)])
  dens <- ssd[closestys] / normalizer   # normalize to be a probability density
  ll[j] <- sum(log(dens))
}  
qseq[which.max(ll)]

# LR-based 95% CI: from
# http://support.sas.com/documentation/cdl/en/etsug/60372/HTML/default/viewer.htm#etsug_model_sect052.htm
thresh <- qchisq(.975, 1)
# in CI should be new q's such that 2*abs(LL(new q) - LL(best q) <= thresh
LLbest <- max(ll)                          # q=347 (ll=-1186.585) with mat.dim=800
qbest <- qseq[which.max(ll)]
# q is in recruits/cm^2; q*1e4 is recruits/m^2
# note: below indices of qseq hold for qseq from 100 to 700 in 601 steps
qcilower <- min(qseq[2*(abs(ll - LLbest)) <= thresh])  # q=207; qseq[108], ll[108]=-1188.51
qciupper <- max(qseq[2*(abs(ll - LLbest)) <= thresh])  # q=600; qseq[501], ll[501]=-1188.539
plot(qseq*1e4, ll, type="l", xlab=paste("Q: recruits per m^2 of reproductive female colonies"), 
     ylab="Log likelihood")
abline(h=LLbest-thresh/2, lty=2)
abline(v=qseq[which.max(ll)]*1e4, lty=3)
abline(v=qcilower*1e4, lty=3)
abline(v=qciupper*1e4, lty=3)

# having estimated best Q using log likelihood, calculate IPM
# Use a smaller matrix, because a little faster
set.seed(1234)
mat.dim <- 300     # matrix dimension (mat.dim x mat.dim)
# don't bother with boundary points, just set up the mesh points,
# which will include the min and max sizes.  Mesh points will be equally
# spaced on the transformed scale (if using a transform)
y <- minsize + c(0:(mat.dim-1))*(maxsize - minsize)/(mat.dim-1)

q <- qbest
r.params <- c(sex.mod$coeff[1], sex.mod$coeff[2],
              repro.mod$coeff[1], repro.mod$coeff[2],
              q, recruitsize, transformtype)

M <- bigmatrix(y, s.params, r.params, g.params)
mat <- M$G %*% M$S + M$R
# in M$G, the second index is the current size; first index is next size
Pvers.eigen <- get.eigen.stuff(mat)
lambda <- Pvers.eigen[[1]]
ssd <- as.vector(Pvers.eigen[[2]])      # Stable size distribution
rv <- as.vector(Pvers.eigen[[3]]/Pvers.eigen[[3]][1]) # Reproductive value


# number of recruits per m2 if the population has the size structure
# of the ssd
sum((qbest*1e4 * Pr.female(y, tr=transformtype) * Pr.repro(y) * ssd)[-1])
# 49.0
# 95% CI for number of recruits per m^2 if the population has the size structure of the ssd
sum((qcilower*1e4 * Pr.female(y, tr=transformtype) * Pr.repro(y) * ssd)[-1])
# 29.3 is lower edge of confidence interval
sum((qciupper*1e4 * Pr.female(y, tr=transformtype) * Pr.repro(y) * ssd)[-1])
# 84.9 is upper edge of CI
# mean recruits/m^2 over all sizes, unweighted by frequency of sizes: 61.9
mean((qbest*1e4 * Pr.female(y, tr=transformtype) * Pr.repro(y))[-1])
# area under the curve of # of recruits/m^2 given ssd and
# reproductivity and sex (but excluding recruit size)
sum((qbest * Pr.female(y, tr=transformtype) * Pr.repro(y) * ssd * y^6)[-1])
# 0.2535515
# max number of recruits is at y[159], i.e. .0031 recruits
# look for window of sizes around y[159] that produces 50% of recruits:
# about 50% of recruits are produced between sizes y[137] and y[180]
sum((qbest * Pr.female(y, tr=transformtype) * Pr.repro(y) * ssd * y^6)[137:180]) /
  sum((qbest * Pr.female(y, tr=transformtype) * Pr.repro(y) * ssd * y^6)[-1])
y[137]^6    # 67 cm2
y[180]^6    # 184 cm2

# calculate eviction probability (colonies getting too big for matrix)
# look at sizes up to 10^8 cm^2, which is effectively Inf
mat.dim <- 1000
infsize <- transformSize(10^8, transformtype)
y.bigger <- minsize + c(0:(mat.dim-1))*(infsize - minsize)/(mat.dim-1)
M.bigger <- bigmatrix(y.bigger, s.params, r.params, g.params)
mat.bigger <- as.matrix(Matrix(M.bigger$G) %*% Matrix(M.bigger$S) + Matrix(M.bigger$R))
# in M$G, the second index is the current size; first index is next size
Pvers.eigen.bigger <- get.eigen.stuff(mat.bigger)
ssd.bigger <- as.vector(Pvers.eigen.bigger[[2]]) # Stable size distribution
sum(ssd.bigger[y.bigger>maxsize])
# probability of eviction: about 1.1e-6
# end of calculating eviction probability

# elasticity for each component model param individually, by adjusting each 
# by +/-1%, and seeing what its effect on population lambda is
mat.dim <- 300
nparams <- 10       # total number of component model params
ssd.sn <- matrix(0, nrow=nparams*2, ncol=mat.dim)
lam.sn <- rep(0, nparams*2)
r.params <- c(sex.mod$coeff[1], sex.mod$coeff[2],
              repro.mod$coeff[1], repro.mod$coeff[2],
              qbest, recruitsize, transformtype)

all.params.orig <- c(g.params, s.params, r.params)
for (i in 1:nparams) {         # tweak each param in turn
  for (j in 1:2) {             # by 1% down and then 1% up
    all.params <- all.params.orig
    all.params[i] <- all.params[i]*(1+2*(j-1.5)/100)
    M.el <- bigmatrix(y, all.params[4:5], all.params[6:12], all.params[1:3])  
    mat.el <- as.matrix(Matrix(M.el$G) %*% Matrix(M.el$S) + Matrix(M.el$R))
    Pvers.eigen.el <- get.eigen.stuff(mat.el)
    lam.sn[2*(i-1)+j] <- Pvers.eigen.el[[1]]
    ssd.sn[2*(i-1)+j, ] <- as.vector(Pvers.eigen.el[[2]]) # Stable size distribution
  }
}
# for plotting purposes, order the params by how big a difference they make to lambda
param.order <- order(abs(lam.sn[seq(2, 20, by=2)] - lam.sn[seq(1, 19, by=2)]), decreasing=T)

# end elasticity for each param

# elasticity & sensitivity for transitions from any state to another
# based on Merow et al. 2014: but Merow et al. have 
# sens <- outer(rv, ssd)/v.dot.w
# elas <- matrix(as.vector(sens)*as.vector(mat)/lambda, nrow=mat.dim)
# which produce a matrix that I think has flipped x/y.
# rows of elas would be size at time t, cols would be size at t+1 (and likewise for sens)
v.dot.w <- sum(ssd * rv)
sens <- outer(ssd, rv)/v.dot.w
elas <- matrix(as.vector(mat) * as.vector(sens)/lambda, nrow=mat.dim)
# now rows of elas are size at time t+1, cols are size at t

# end elasticity & sensitivity for transitions from any state to another


##### figures #######

alltcks <- c(.35, 1, 10, 20, 50, 100, 250, 500, 1000, 1500)
tcks <- alltcks[c(2, 3, 5:10)]    # ticks that will get labelled
# plots of growth, survival, pr(f), pr(reproductive) functions
xlim <- c(.5, 3.5)
modelcolor <- "black"

postscript("figs/fig1.eps", horizontal=F, onefile=F, paper="special",
           width=8, height=8, pointsize=13)
par(mfrow=c(2, 2))
par(mar=c(3.5, 4.5, 2, .2)+.1)
par(mgp=c(3, 1, 0))
plot(growth.data$size, growth.data$sizeNext, xlim=xlim, ylim=xlim,
     asp=1, xlab="",
     ylab=expression(paste("Size (cm"^2, ") at time t+1")), axes=F)
axis(1, at=tcks^(1/6), labels=paste(tcks), las=2)
axis(2, at=tcks^(1/6), labels=paste(tcks), las=2)
abline(a=grow.mod$coeff[1], b=grow.mod$coeff[2], col=modelcolor, lwd=3)
abline(0, 1, col="black", lty=2, lwd=2)
legend(x="bottomright", legend=c("Growth model", "45 degrees"), 
       col=c(modelcolor, "black"), lty=c(1, 2), lwd=c(3, 2), bty="n")
legend(.0008^(1/6), 3800^(1/6), legend=c("(a)"), bty="n", xpd=NA)
box()

plot(y, surv.f(y, m=.01), ylab="Pr(survival)",
     xlab="", xlim=xlim, ylim=c(-.05, 1.05), type="l",
     lwd=2, col=modelcolor, axes=F)
axis(1, at=tcks^(1/6), labels=paste(tcks), las=2)
axis(2)
par(new=T)
set.seed(1234)
plot(surv.data$size, surv.data$surv+rnorm(nrow(surv.data), 0, .026),
     xlab="", ylab="", xlim=xlim, ylim=c(-.05, 1.05), cex=.9, axes=F)
legend(.0008^(1/6), 1.22, legend=c("(b)"), bty="n", xpd=NA)
box()

par(mar=c(4.5, 4.5, 1, .2)+.1)
plot(y, Pr.repro(y), ylab="Pr(reproductive)",
     xlab=expression(paste("Size (cm"^2, ") at time t")), 
     xlim=xlim, ylim=c(-.05, 1.05), type="l", lwd=2, col=modelcolor, axes=F)
axis(1, at=tcks^(1/6), labels=paste(tcks), las=2)
axis(2)
par(new=T)
set.seed(1234)
plot(repro.size$size, repro.size$repro+rnorm(nrow(repro.size), 0, .026),
     xlim=xlim, ylim=c(-.05, 1.05), cex=.9,
     xlab="", ylab="", axes=F)
legend(.0008^(1/6), 1.22, legend=c("(c)"), bty="n", xpd=NA)
box()

plot(y, Pr.female(y, tr=transformtype), ylab="Pr(female)",
     xlab=expression(paste("Size (cm"^2, ") at time t")),
     xlim=xlim, ylim=c(-.05, 1.05), type="l", lwd=2, col=modelcolor, axes=F)
axis(1, at=tcks^(1/6), labels=paste(tcks), las=2)
axis(2)
par(new=T)
set.seed(1234)
plot(sex.size$size, sex.size$sex+rnorm(nrow(sex.size), 0, .026),
     xlim=xlim, ylim=c(-.05, 1.05), cex=.9,
     xlab="", ylab="", axes=F)
legend(.0008^(1/6), 1.22, legend=c("(d)"), bty="n", xpd=NA)
box()
dev.off()
par(mgp=c(3, 1, 0))
par(mar=c(5, 4, 4, 2)+.1)
# end: plots of growth, survival, pr(f), pr(reproductive) functions

# figure showing Q
# note: below indices of qseq hold for qseq from 100 to 700 in 601 steps
# (used mat.dim=800 in calculating this)

postscript("figs/fig2.eps", horizontal=F, onefile=F, paper="special",
           width=8, height=8, pointsize=12)
par(mfrow=c(2, 2))
par(mar=c(5, 3.7, 2, 4)+.1)
plot(y[1:230], (Pr.female(y, tr=transformtype)*Pr.repro(y))[1:230],
     xlab=expression(paste("Size (cm"^2, ")")), ylab="",
     type="l", lwd=2, cex.lab=1.1, axes=F)
axis(1, at=tcks^(1/6), labels=paste(tcks), las=2)
axis(2)
legend(.08^(1/6), .54, legend=c("(a)"), bty="n", cex=1.2, xpd=NA)
mtext("Probability of reproductive female", 2, line=2.7, cex=.95)
par(new=T)
plot(y[1:230], (qbest*Pr.female(y, tr=transformtype)*Pr.repro(y))[1:230],
     xlab="", ylab="", axes=F,
     type="n", lwd=2)
axis(4, at=seq(0, .0150, by=.0050))
abline(v=50^(1/6), lty=3)
abline(v=100^(1/6), lty=3)
abline(v=250^(1/6), lty=3)
mtext(expression(paste("Recruits per cm"^2)), 4, line=2.7,
      cex=.95)
box()

par(mar=c(5, 4.5, 2, 2)+.1)
plot(qseq, ll, type="l",
     xlab=expression(paste("Q: recruits per cm"^2, " of reproductive female colonies")),
     ylab="Log likelihood", cex.lab=1.1, lwd=2)
abline(h=LLbest-thresh/2, lty=2)
abline(v=qseq[which.max(ll)], lty=3)
abline(v=qseq[108], lty=3)
abline(v=qseq[501], lty=3)
legend(.006, -1183, legend=c("(b)"), bty="n", cex=1.2, xpd=NA)

par(mar=c(5, 3.7, 2, 4)+.1)
plot(y[1:230], (qbest*Pr.female(y, tr=transformtype)*Pr.repro(y)*ssd*y^6)[1:230],
     xlab=expression(paste("Size (cm"^2, ")")), ylab="",
     type="l", lwd=2, cex.lab=1.1, axes=F)
axis(1, at=tcks^(1/6), labels=paste(tcks), las=2)
abline(v=50^(1/6), lty=3)
abline(v=100^(1/6), lty=3)
abline(v=250^(1/6), lty=3)
mtext("Relative recruits under", 2, line=2.7, cex=.95)
mtext("stable size distribution", 2, line=1.4, cex=.95)
legend(.08^(1/6), .0037, legend=c("(c)"), bty="n", cex=1.2, xpd=NA)
box()
# end figure showing Q

# plot of recruits per colony
par(mar=c(5, 4.5, 2, 2)+.1)
plot(y[1:230], (qbest*1e4*Pr.female(y, tr=transformtype)*Pr.repro(y)*y^6/1e4)[1:230],
     xlab=expression(paste("Size (cm"^2, ")")),
     ylab="Recruits per colony", axes=F, cex.lab=1.1,
     type="l", lwd=2)
axis(1, at=tcks^(1/6), labels=paste(tcks), las=2)
axis(2)
abline(v=50^(1/6), lty=3)
abline(v=100^(1/6), lty=3)
abline(v=250^(1/6), lty=3)
box()
legend(.08^(1/6), 2.83, legend=c("(d)"), bty="n", cex=1.2, xpd=NA)
dev.off()
# end plot of recruits per colony
par(mar=c(5, 4, 4, 2)+.1)

# figure showing stable size distribution & empirical size structure
# with CI around ssd based on CI for Q
# ssd and confidence interval around it (based on CI around Q)
r.params.cil <- c(sex.mod$coeff[1], sex.mod$coeff[2],
              repro.mod$coeff[1], repro.mod$coeff[2],
              qcilower, recruitsize, transformtype)
M.cil <- bigmatrix(y, s.params, r.params.cil, g.params)
mat.cil <- as.matrix(Matrix(M.cil$G) %*% Matrix(M.cil$S) + Matrix(M.cil$R))
Pvers.eigen.cil <- get.eigen.stuff(mat.cil)
lambda.cil <- Pvers.eigen.cil[[1]]
ssd.cil <- as.vector(Pvers.eigen.cil[[2]]) 
r.params.ciu <- c(sex.mod$coeff[1], sex.mod$coeff[2],
                  repro.mod$coeff[1], repro.mod$coeff[2],
                  qciupper, recruitsize, transformtype)

M.ciu <- bigmatrix(y, s.params, r.params.ciu, g.params)
mat.ciu <- as.matrix(Matrix(M.ciu$G) %*% Matrix(M.ciu$S) + Matrix(M.ciu$R))
Pvers.eigen.ciu <- get.eigen.stuff(mat.ciu)
lambda.ciu <- Pvers.eigen.ciu[[1]]
ssd.ciu <- as.vector(Pvers.eigen.ciu[[2]])

finalfigure <- F
if (finalfigure) 
  tiff("figs/fig3.tif", pointsize=48, width=2000, height=2000, antialias="none", compression="lzw")
if (!finalfigure) 
  png(file="figs/fig3.png", pointsize=62, width=2500, height=2500, antialias="none")
par(mfrow=c(2, 2))
par(mar=c(5, 4.5, 1.5, 3.5)+.1)
plot(M$meshpts, ssd, type='l', lwd=2, lty=4,
     xlab=expression(paste("Size (cm"^2, ")")),
     ylab="IPM stable size distribution",
     xlim=c(minsize, 3.5),
     ylim=range(ssd[-1]), axes=F)
axis(1, at=alltcks^(1/6), labels=paste(alltcks), las=2)
axis(2)
box()
polygon(c(M.cil$meshpts, rev(M.ciu$meshpts)),
        c(ssd.cil, rev(ssd.ciu)), col=rgb(0, 0, 0, .2), border=NA)
par(new=T)
hist(transformSize(size.struct[, 2], transformtype), breaks=20,
     xlim=c(minsize, 3.5), axes=F, ylab="", xlab="", main="", freq=F,
     ylim=c(min(ssd[-1]), 1.11428571*2.023406))
axis(4, at=seq(0, 2, by=.2))
mtext("Empirical size structure", 4, line=2.2, cex=0.85)
abline(v=transformSize(minreliablesize, transformtype), col="black",
       lwd=2, lty=3)
legend(.08^(1/6), 2.7, legend=c("(a)"), bty="n", cex=1.2, xpd=NA)
par(mar=c(5, 4, 4, 2)+.1)
# end: figure showing stable size distribution & empirical size structure

# plot kernel matrix (K in literature; mat here)
par(mar=c(5, 4.7, 1, 2.5)+.1)
par(mar=c(5, 4.7, 1.5, 2)+.1)
broaden.recruits <- c(rep(1, 3), 2:length(y))
y.broaden.recruits <- c(seq(y[1]*.9, y[1]*.99, length=2), y)
image(y.broaden.recruits, y.broaden.recruits,
      t(log(mat))[broaden.recruits, broaden.recruits],
      xlab=expression(paste("Size (cm"^2, ") at time t")),
      ylab=expression(paste("Size (cm"^2, ") at time t+1")),
      col=tim.colors(),
      axes=F, asp=1)
axis(1, at=alltcks^(1/6), labels=paste(alltcks), las=2)
axis(2, at=alltcks^(1/6), labels=paste(alltcks), las=2)
lines(x=c(y[1], y.broaden.recruits[length(y.broaden.recruits)]),
      y=rep(y[1], 2), col="white", lwd=1)
lines(y=c(y[1], y.broaden.recruits[length(y.broaden.recruits)]),
      x=rep(y[1], 2), col="white", lwd=1)
abline(a=0, b=1, col="white", lty=3, lwd=1.5)
legend(.08^(1/6), 3100^(1/6), legend=c("(b)"), bty="n", cex=1.2, xpd=NA)
# end plot kernel matrix (K in literature; mat here)

# plot elasticity & sensitivity: with recruits emphasized and diagonal line
# the axes are just at the edges of the plots if the 
# aspect ratio of the plot is just right, and away from the 
# plots otherwise
par(mar=c(5, 4.7, 1.5, 2)+.1)
broaden.recruits <- c(rep(1, 2), 2:length(y))
y.broaden.recruits <- c(seq(y[1]*.9, y[1]*.99, length=2), y)
# rows of sens are size at time t+1, cols at time t
image(y.broaden.recruits, y.broaden.recruits,
      t(log(sens))[broaden.recruits, broaden.recruits],
      xlab=expression(paste("Size (cm"^2, ") at time t")),
      ylab=expression(paste("Size (cm"^2, ") at time t+1")),
      col=tim.colors(),
      axes=F, asp=1)
axis(1, at=alltcks^(1/6), labels=paste(alltcks), las=2)
axis(2, at=alltcks^(1/6), labels=paste(alltcks), las=2)
lines(x=c(y[1], y.broaden.recruits[length(y.broaden.recruits)]),
      y=rep(y[1], 2), col="white", lwd=1)
lines(y=c(y[1], y.broaden.recruits[length(y.broaden.recruits)]),
      x=rep(y[1], 2), col="white", lwd=1)
abline(a=0, b=1, col="black", lty=3, lwd=1.5)
legend(.08^(1/6), 3100^(1/6), legend=c("(c)"), bty="n", cex=1.2, xpd=NA)

image(y.broaden.recruits, y.broaden.recruits,
           t(log(elas))[broaden.recruits, broaden.recruits],
      xlab=expression(paste("Size (cm"^2, ") at time t")),
      ylab=expression(paste("Size (cm"^2, ") at time t+1")),
      col=tim.colors(),
      axes=F, asp=1)
axis(1, at=alltcks^(1/6), labels=paste(alltcks), las=2)
axis(2, at=alltcks^(1/6), labels=paste(alltcks), las=2)
par(new=T)
lines(x=c(y[1], y.broaden.recruits[length(y.broaden.recruits)]),
      y=rep(y[1], 2), col="white", lwd=1)
lines(y=c(y[1], y.broaden.recruits[length(y.broaden.recruits)]),
      x=rep(y[1], 2), col="white", lwd=1)
abline(a=0, b=1, col="white", lty=3, lwd=1.5)
legend(.08^(1/6), 3100^(1/6), legend=c("(d)"), bty="n", cex=1.2, xpd=NA)
dev.off()
# end plot elasticity & sensitivity: with recruits emphasized

# plot elasticity results for each param
postscript("figs/fig2A.eps", horizontal=F, onefile=F, paper="special",
           width=6, height=3.5, pointsize=10)
par(mfrow=c(1, 2))
par(mar=c(9.5, 4, 1, 1)+.1)
x <- rep(1:nparams, each=2)
lam.sn.order <- rep(0, nparams*2)
lam.sn.order[(1:nparams)*2] <- param.order*2
lam.sn.order[(1:nparams)*2-1] <- param.order*2-1
plot(x, lam.sn[lam.sn.order], type="n", xlab="", ylab=expression(paste(lambda)),
     xaxt="n")
for (i in 1:nparams) {
  arrows(x[2*i-1], lam.sn[2*param.order[i]-1],
         x[2*i-1], lam.sn[2*param.order[i]],
         length=.05, angle=90, code=3)
}
param.labels <- c("Growth intercept", "Growth slope", "Growth variance",
                  "Survival intercept", "Size-indep. mortality",
                  "Sex intercept", "Sex slope", "Reproduction intercept",
                  "Reproduction slope", "Q")
axis(1, labels=param.labels[param.order], at=1:nparams, padj=0, las=3)
legend(.56, 1.189, legend=c("(a)"), bty="n", cex=1.2, xpd=NA)
par(mar=c(5, 4, 1, 1)+.1)

plot(y[-1], ssd.sn[1, -1], ylim=c(0, max(ssd.sn[, -1])),
     xlab=expression(paste("Size (cm"^2, ")")),   
     ylab="IPM stable size distribution", type="n", axes=F)
axis(1, at=alltcks^(1/6), labels=paste(alltcks), las=2)
axis(2)
for (i in 1:(nparams*2)) {
  par(new=T)
  if (i==3)
    plot(y[-1], ssd.sn[i, -1], type="l", axes=F, xlab="", ylab="",
         lty=2, lwd=2)
  if (i==4)
    plot(y[-1], ssd.sn[i, -1], type="l", axes=F, xlab="", ylab="",
         lty=3, lwd=2)
  if (i < 3 || i > 4)
    plot(y[-1], ssd.sn[i, -1], type="l", axes=F, xlab="", ylab="")
}
legend(.08^(1/6), .00955, legend=c("(b)"), bty="n", cex=1.2, xpd=NA)
box()
dev.off()
par(mar=c(5, 4, 4, 2)+.1)
# end plot elasticity results for each param

### end figures
