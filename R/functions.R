# functions in support of IPM model for Plesiastrea versipora

transformSize <- function(x, tr=1, inverse=0) {
  inv <- rep(inverse, length(x))
  if (tr==1)
    return(ifelse(inv, 10^x, log10(x)))
  if (tr==2)
    return(ifelse(inv, x^4, x^0.25))
  if (tr==3)
    return(ifelse(inv, x^6, x^(1/6)))
  if (tr==4)
    return(1/x)
  if (tr==0)
    return(x)
  
  stop("transformSize: Unrecognized value for tr")
}

# Survival function S(z, t) 
# z is raw (size in cm^2) or transformed
#surv.f <- function(z, a=surv.mod$coeff[1], b=surv.mod$coeff[2], m=0) {
surv.f <- function(z, b=surv.mod$coeff[1], m=0) {
  # m is size-independent mortality
  
  if (m < 0 | m > 1) { 
    warning("problem in surv.f: m not in [0,1]; setting m to 0 or 1"); 
    m <- ifelse(m < 0, 0, 1) 
  }
  u <- exp(b*z)
  return((1-m)*u/(1+u))
}

# size-dependent probability of being female Pr.female(z)
# z is size
Pr.female <- function(z, a=sex.mod$coeff[1], b=sex.mod$coeff[2], tr) {
  # undo the transform to get back to units of cm^2
  z.untr <- transformSize(z, tr, inverse=1)
  
  u <- exp(a + b*z.untr)
  return(u/(1+u))
}

# size-dependent probability of a colony being reproductive
Pr.repro <- function(z, a=repro.mod$coeff[1], b=repro.mod$coeff[2]) {
  u <- exp(a + b*z)
  return(u/(1+u))
}

# growth function G(zz|z), i.e. the probability of a new size zz given
# the prior size z (in cm2), or zz|z in transformed units.  Modelling zz 
# as a Gaussian distribution given z,
# because the residuals of the growth model are kinda sorta normal (reasonably
# symmetric, anyway, although both tails are too heavy)
# note the output of grow.f will need to be scaled so that G(zz|z) over
# all zz for a given z sum up to 1
grow.f <- function(z, zz, mean.intercept=grow.mod$coeff[1],
                   mean.slope=grow.mod$coeff[2],
                   var.intercept=summary(grow.mod)$sigma^2) {
  mean.zz <- mean.intercept + z * mean.slope
  return(dnorm(zz, mean.zz, sqrt(var.intercept)))
}

# recruitment function: this describes the number of recruits arriving
# in the population; implies it's a closed system
# from Edmunds et al. 2014: "Closed system modeling can be justified if 
# inter-connected populations all tend to occupy similar habitats and 
# environmental changes operate at scales larger than the meta-population, 
# and therefore affect all populations similarly. In this case, a closed 
# meta-population model will provide an approximation for local dynamics"
# recruitsize is assumed max size of a recruit in its first year, in cm2
# q is recruits per square cm of colony area
recruit.f <- function(z, y, q, recruitsize, tr) {
  # z is potential parent; y is potential recruit
  # recruitsize is the max size of recruits that are produced
  
  # undo the transform to get back to units of m^2
  z.untr <- transformSize(z, tr, inverse=1)
  
  #  n <- sum(abs(y-recruitsize) < .03)
  
  r <- q * z.untr
  r[y > recruitsize] <- 0
  #  r[y > recruitsize || z <= recruitsize] <- 0
  #  r[abs(y-recruitsize)>=.03] <- 0
  
  return(r)
}

fecund.f <- function(z, y, f.intercept=sex.mod$coeff[1],
                     f.slope=sex.mod$coeff[2],
                     r.intercept=repro.mod$coeff[1],
                     r.slope=repro.mod$coeff[2],
                     q, recruitsize, tr) {
  f <- Pr.female(z, f.intercept, f.slope, tr=tr) * 
    Pr.repro(z, r.intercept, r.slope) * 
    recruit.f(z, y, q, recruitsize, tr)
  return(f)
}

############## The 'big matrix' M of size n x n
#### Thanks to Stephen P. Ellner for code -- taken & modified from Coulson 2012,
# who took it from SPE's website
bigmatrix<-function(y, s.params, r.params, g.params) {  
  
  # create S, R, and G matrices
  S <- (diag(surv.f(y, s.params[1], s.params[2])))
  R <- (t(outer(y, y, fecund.f, r.params[1], r.params[2], r.params[3],
                r.params[4], r.params[5], r.params[6], r.params[7])))
  G <- (t(outer(y, y, grow.f, g.params[1], g.params[2], g.params[3])))
  
  # scale G so columns sum to 1
  G <- G/matrix(as.vector(apply(G, 2, sum)), nrow=length(y), ncol=length(y), byrow=TRUE)
  return(list(S=S, R=R, G=G, meshpts=y))
}

# For large matrices it is faster to calculate the dominant eigenvalue and 
# eigenvectors associated with it via iteration rather than using eigen.
get.eigen.stuff <- function(mat){ 
  sz <- dim(mat)[1]
  t.now <- runif(sz)
  t.now <- t.now/sum(t.now)
  t.next <- mat%*%t.now
  t.next <- t.next/sum(t.next)
  i <- 0
  while (sum(abs(t.next-t.now))>0.0000001){
    i <- i+1
    t.now <- t.next
    t.next <- mat%*%t.now
    lambda <- sum(t.next)/sum(t.now)
    t.next <- t.next/sum(t.next)
  }
  r.now <- runif(sz)
  r.now <- r.now/sum(r.now)
  r.next <- r.now%*%mat
  r.next <- r.next/sum(r.next)
  while (sum(abs(r.next-r.now))>0.0000001){
    r.now <- r.next
    r.next <- r.now%*%mat
    r.next <- r.next/sum(r.next)
  }
  return(list(lambda, t.next, r.next))
}

# Box-Cox code from Rees, Childs, Ellner 2014, supplement 5
# for the growth model, both the predictor and response variables are
# size, so want to apply the same transformation to both
bcNLL <- function(lambda, x, y) {
  xl <- x^lambda
  yl <- y^lambda
  fit <- lm(yl ~ xl)
  s2hat <- mean(fit$residuals^2)
  return(0.5*length(x)*log(s2hat/lambda^2) - (lambda-1)*sum(log(y)))
}

# code for the image.scale function below is from http://menugget.blogspot.de/2013/12/new-version-of-imagescale-function.html
# Modified slightly to handle log probabilities (sets zlim max to <= 0),
# to not draw a box around the scale, and to use only a very fine axis line

# This function creates a color scale for use with the image()
# function. Input parameters should be consistent with those
# used in the corresponding image plot. The "axis.pos" argument
# defines the side of the axis. The "add.axis" argument defines
# whether the axis is added (default: TRUE) or not (FALSE).
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, logprob=F, ...) {
  if (!missing(breaks)) {
    if (length(breaks) != (length(col)+1)) {stop("must have one more break than colour")}
  }
  if (missing(breaks) & !missing(zlim)) {
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if (missing(breaks) & missing(zlim)) {
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3) # adds a bit to the range in both directions
    if (logprob) zlim[2] <- min(0, zlim[2])
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for (i in seq(poly)) {
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if (axis.pos %in% c(1, 3)) {ylim<-c(0, 1); xlim<-range(breaks)}
  if (axis.pos %in% c(2, 4)) {ylim<-range(breaks); xlim<-c(0, 1)}
  plot(1, 1, t="n", ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for (i in seq(poly)) {
    if (axis.pos %in% c(1, 3)) {
      polygon(poly[[i]], c(0, 0, 1, 1), col=col[i], border=NA)
    }
    if (axis.pos %in% c(2, 4)) {
      polygon(c(0, 0, 1, 1), poly[[i]], col=col[i], border=NA)
    }
  }
  if (add.axis) {axis(axis.pos, lwd=0, lwd.ticks=1)}
}
