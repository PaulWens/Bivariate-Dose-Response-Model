# FIT BIVARIATE DOSE-RESPONSE MODEL USED IN CURE ET AL. 2025
# v13. Fits the models of the published paper.

library(rjags) # This package seems to need R v4.0 or below (start RStudio holding down Control button)
library(car) # To plot the estimated bivariate normal distribution.
library(mixtools) # To plot ellipse
library(truncnorm)
library(matlab)
library(plot3D)
library(gdata)
library(coda)
load.module("dic")

# Function to draw ellipse for bivariate normal data
ellipse_bvn <- function(bvn, alpha, col){
  Xbar <- apply(bvn,2,mean)
  S <- cov(bvn)
  ellipse(Xbar, S, alpha = alpha, col=col)
}

# Select dose metric & severity level
dose="SPL" # Choose SPL, SEL or SELcum
severity="medium" # Choose low: 1-3, medium: 1-6, or high: 1-9

# Prepare data
o <- read.csv('Table_dose_response_modelling_v3.csv') 
if(severity=="low"){o = o[o$Severity==1,]}
if(severity=="medium"){o = o[o$Severity==2,]}
if(severity=="high"){o = o[o$Severity==3,]}
whale.id = o$Whale
lower.bound = c(60,-50)
upper.bound = c(200,0)  
measurement.precision = c(2.5,0.25)^(-2)
beta.precision = c(30,10)^(-2)
if(dose=="SPL"){y = cbind(o$SPLmax, -o$Rangemin/1000)} # onset thresholds
if(dose=="SEL"){y = cbind(o$SELmax, -o$Rangemin/1000)} # onset thresholds
if(dose=="SELcum"){y = cbind(o$SELcum, -o$Rangemin/1000)} # onset thresholds
n.whales = max(whale.id)
n.sessions = length(y[,1])
L = cbind(rep(upper.bound[1],n.sessions), rep(upper.bound[2],n.sessions)) # right censoring lower bounds
L[o$Censored==1,] = y[o$Censored==1,]
y[o$Censored==1,] = NA
I.censored = cbind(o$Censored,o$Censored)  # Right censored y/n

# Plot data with 95% confidence ellipse
par(mfrow=c(1,1))
plot(y, cex=2, xlim=c(60, 200), ylim=c(-25, 0), col="black", 
     xlab=expression(paste("SELmax", " (dB re 1 ", muPa^2," s)")), ylab="DISTmin (km)")
ellipse_bvn(y[!is.na(y[,1]),], .05, "black")
points(L[L[,1]<200,],col="black",pch=2)


# JAGS code of model
model_string <- "
# with trick that uses correlated part of residual in x[,1] in mean function of x[,2]

model{

###priors###
# Population mean thresholds
mu[1]~dunif(lower.bound[1],upper.bound[1]) # RL
mu[2]~dunif(lower.bound[2],upper.bound[2]) # Distance

# Between-whale variation in threshold (sd)
phi[1]~dunif(0,30) # RL
phi[2]~dunif(0,10) # Distance      # PRIOR ASSUMPTION OF NO CORRELATION
inv_phi2 <- phi^(-2) # precision

# Within-whale between-session variation in threshold (sd)
sigma[1]~dunif(0,30) # RL         
sigma[2]~dunif(0,10) # Distance    
inv_sigma2 <- sigma^(-2) # precision

# Parameters for correlation between individual threshold RL and Distance
# for between-whale variation
# lambda = rho * sdy / sdx
lam ~ dunif(-0.333333,0.333333)
var1[1] <- phi[1]^2
var2[1] <- phi[2]^2 + lam^2 * var1[1]
rho[1] <- lam * var1[1] / sqrt(var1[1]*var2[1]) # correlation coefficient
phi.tot <- sqrt(var2[1]) # sd of the correlated and uncorrelated parts combined

# Precision matrix for observation model
iSigma[1,1] <- measurement.precision[1] # RL
iSigma[1,2] <- 0
iSigma[2,1] <- 0
iSigma[2,2] <- measurement.precision[2] # Distance

###process model###
# Between whales
for(i in 1:n.whales){
  mu_i[i,1]~dnorm(mu[1], inv_phi2[1]) T(lower.bound[1],upper.bound[1]) # RL
  mu_i[i,2]~dnorm(mu[2] + lam * (mu_i[i,1] - mu[1]), inv_phi2[2]) T(lower.bound[2],upper.bound[2]) # Distance 
}

# Between sessions within whale
for(j in 1:n.sessions){
  t[j,1]~dnorm(mu_i[whale.id[j],1], inv_sigma2[1]) T(lower.bound[1],upper.bound[1]) # RL
  t[j,2]~dnorm(mu_i[whale.id[j],2], inv_sigma2[2]) T(lower.bound[2],upper.bound[2]) # Distance
}

####observation model###
for(j in 1:n.sessions){
  I.censored[j,1]~dinterval(t[j,1],L[j,1]) # RL
  I.censored[j,2]~dinterval(t[j,2],L[j,2]) # Distance
  y[j,1:2]~dmnorm(t[j,], iSigma[,])
}
}
"
data = list(y=y,n.whales=n.whales,n.sessions=n.sessions,lower.bound=lower.bound,upper.bound=upper.bound, 
            L=L,measurement.precision=measurement.precision,whale.id=whale.id,I.censored=I.censored)
#data = list(y=c(), n.whales=0, n.sessions=0, lower.bound=lower.bound, upper.bound=upper.bound, 
#            L=c(), measurement.precision=measurement.precision, whale.id=whale.id,I.censored=matrix(0,n.sessions,2)) # to plot prior distributions

# Use midpoints as starting values (not sensitive to starting values)
rl.mid = lower.bound[1]+(upper.bound[1]-lower.bound[1])/2
r.mid = lower.bound[2]+(upper.bound[2]-lower.bound[2])/2
t.init = cbind(rep(rl.mid,n.sessions), rep(r.mid,n.sessions)) # non-censored
L.mid = cbind(L[,1] + (upper.bound[1] - L[,1])/2, L[,2] + (upper.bound[2] - L[,2])/2)
t.init[o$Censored==1,] = L.mid[o$Censored==1,] # censored
inits = list(t = t.init, mu = c(rl.mid,r.mid), phi = c(15,5), sigma = c(15,5))

set.seed(1)
m1 <- jags.model(textConnection(model_string), data = data, inits = inits, # bivariate model
                         n.adapt = 500, n.chains = 3, quiet = F)

update(m1, 5000)
n.iter = 100000 # Use 20k for quick run, 100k with thin=10 for final
thin1 = 10
samples.all <- coda.samples(m1, c("mu","sigma","rho","phi","phi.tot","lam","mu_i"), #, "deviance" 
                            n.iter = n.iter, thin=thin1)

# Backtransform Range by flipping sign
samples.all2 = samples.all
samples.all2[[1]][,"mu[2]"] = -samples.all[[1]][,"mu[2]"] 
samples.all2[[2]][,"mu[2]"] = -samples.all[[2]][,"mu[2]"] 
samples.all2[[3]][,"mu[2]"] = -samples.all[[3]][,"mu[2]"] 
samples.all2[[1]][,"rho"] = -samples.all[[1]][,"rho"] 
samples.all2[[2]][,"rho"] = -samples.all[[2]][,"rho"] 
samples.all2[[3]][,"rho"] = -samples.all[[3]][,"rho"] 
samples.all2[[1]][,"lam"] = -samples.all[[1]][,"lam"] 
samples.all2[[2]][,"lam"] = -samples.all[[2]][,"lam"] 
samples.all2[[3]][,"lam"] = -samples.all[[3]][,"lam"] 


# Assess model fit
stats = summary(samples.all2)
print(stats)
gelman.diag(samples.all2)

#plot(samples.all2)
#traceplot(samples.all2)
#densplot(samples.all2)
#gelman.plot(samples.all2)
par(mar=c(4.1,4.4,2,2), lwd=2, las=1, cex=1.1)
par(mfrow=c(3,3))
gelman.plot(samples.all2[,c(2,3,1,28:30,32,33,31)])
traceplot(samples.all2[,c(2,3,1,28:30,32,33,31)])
densplot(samples.all2[,c(2,3,1,28:30,32,33,31)])
lines(x=c(-0.39,-0.39),y=c(0,8),col="red", lty=2)

# Compare correlation coefficient of raw data to posterior of rho
samples.all.rho <- rbind(as.matrix(samples.all2[[1]][,"rho"]),as.matrix(samples.all2[[2]][,"rho"]),as.matrix(samples.all2[[3]][,"rho"]))
quantile(samples.all.rho,c(0.5,0.025,0.975)) # median and 95% credible interval
HPDinterval(as.mcmc(samples.all.rho), prob = 0.95) # 95% HPD interval
Fn.ecdf <- ecdf(samples.all.rho)
#Fn.ecdf <- ecdf(samples.all2[[1]][,"rho"]) # no different between chains for 2 decimals precision at n.iter=100k
#Fn.ecdf <- ecdf(samples.all2[[2]][,"rho"])
#Fn.ecdf <- ecdf(samples.all2[[3]][,"rho"])
Fn.ecdf(-0.41) # quantile for SPLmax-DISTmin
Fn.ecdf(-0.39) # quantile for SELmax-DISTmin


#### Plot dose-response functions
# note: range not backtransformed here
thin2 = 10 # use at least 1k iter after thinning for plotting
samples = samples.all[1][seq(1,n.iter/thin1, by=thin2),] # Keep one chain and thin

# dose vector [RL,Range]
Ld = 100 # length of dose vectors (100 is nice, 50 acceptable)
dose1<-cbind(seq(lower.bound[1],upper.bound[1],length=Ld), seq(lower.bound[2],upper.bound[2],length=Ld)) # dose for 1D DR functions
dose2<-meshgrid(seq(lower.bound[1],upper.bound[1],length=Ld), seq(lower.bound[2],upper.bound[2],length=Ld)) # 2D grid of LdxLd evenly spaced points
dose3<-cbind(as.vector(dose2$x), as.vector(dose2$y)) # input for distribution function

# summary statistics
est<-list()
est$mu<-cbind(quantile(samples[[1]][,"mu[1]"],c(0.025,0.25,0.5,0.75,0.975)), quantile(samples[[1]][,"mu[2]"],c(0.025,0.25,0.5,0.75,0.975)))
est$phi<-cbind(quantile(samples[[1]][,"phi[1]"],c(0.025,0.25,0.5,0.75,0.975)), quantile(samples[[1]][,"phi.tot"],c(0.025,0.25,0.5,0.75,0.975)))
est$sigma<-cbind(quantile(samples[[1]][,"sigma[1]"],c(0.025,0.25,0.5,0.75,0.975)), quantile(samples[[1]][,"sigma[2]"],c(0.025,0.25,0.5,0.75,0.975)))

# posterior samples
est$sds<-cbind(samples[[1]][,"sigma[1]"], samples[[1]][,"sigma[2]"]) 
#est$rhos<-samples[[1]][,"rho[2]"]

# Dose-response surface for each sample
Li = n.iter/thin1/thin2 # number of iterations for plotting
responses <- matrix(0,Li*n.whales,Ld*Ld)
for(i in 1:Li){ 
  for(j in 1:n.whales){
    mu <-cbind(samples[[1]][,paste0("mu_i[",j,",1]")][i], samples[[1]][,paste0("mu_i[",j,",2]")][i])
    p1 <- ptruncnorm(dose2$x,lower.bound[1],upper.bound[1],mean=mu[1],sd=est$sds[i,1])
    p2 <- ptruncnorm(dose2$y,lower.bound[2],upper.bound[2],mean=mu[2],sd=est$sds[i,2]) 
    responses[j+n.whales*(i-1),] <- p1 * p2
  }
  print(paste(i,"/",Li,"samples processed"))
}
response<-apply(responses, 2, median) # posterior median surface
qtile<-apply(responses, 2, quantile, probs=c(0.005,0.025,0.25,0.5,0.75,0.975,0.995))
qtile05<-apply(responses, 2, quantile, probs=0.05)
qtile95<-apply(responses, 2, quantile, probs=0.95)

# 3D surface of median posterior probability
Z <- reshape(as.array(response),Ld,Ld)

# Other posterior probabilities
Z.upper <- reshape(as.array(qtile[6,]),Ld,Ld) # 0.975
Z.lower <- reshape(as.array(qtile[2,]),Ld,Ld) # 0.025
Z.upper95 <- reshape(as.array(qtile95),Ld,Ld) # 0.95
Z.lower05 <- reshape(as.array(qtile05),Ld,Ld) # 0.05

if(dose=="SPL"){
  RLlabel = "SPLmax (dB re 1 uPa)"
  RLlabel2 = expression(paste("SPLmax", " (dB re 1 ", mu, "Pa)"))}
if(dose=="SEL"){
  RLlabel = "SELmax (dB re 1 muPa^2 s)"
  RLlabel2 = expression(paste("SELmax", " (dB re 1 ", mu, Pa^2," s)"))}
if(dose=="SELcum"){
  RLlabel = "SELcum (dB re 1 uPa^2 s)"
  RLlabel2 = expression(paste("SELmax (dB re 1 ",mu,Pa^2," s)"))}

# Values for MS
round(approx(Z[Ld,], dose2$x[Ld,], 0.5)$y) # SEL at p50 and 0 m
round(approx(Z[94,], dose2$x[94,], 0.5)$y) # SEL at p50 and ~3 km
approx(Z.upper[,Ld], dose2$y[,Ld], 0.05)$y # distance at p05 and 200 dB (max no-response in 0.975)


################# MAKE PLOT #####################
# use 1000x825
layout(matrix(c(1,1,1,1,2,3), 3, 2, byrow = TRUE), heights=c(1,1,1.5)) #layout.show(3)
par(mar=c(4.1,4.4,2,2), lwd=2, las=2, cex=1.1)
#par(mfrow=c(1,1))
irl = which(dose2$x[1,]>=90) # limit rl axis
irange = which(dose2$y[,1]>=-11.5) # limit range axis
persp3D(x=dose2$x[irange,irl], y=-dose2$y[irange,irl], z=Z[irange,irl], shade = 0.1, phi = 20, theta = -45-90, 
        ticktype = "detailed", bty = "b2", lphi = 30, ltheta = 80, expand = 0.5, 
        clab = "", zlab = "P(response)", xlab=RLlabel, ylab="DISTmin (km)", main="a)")

contour2D(x=dose1[,1], y=-dose1[,2], z=t(Z), lwd = 2, colkey = F, 
          clab = "p(response)", xlab=RLlabel2, ylab="DISTmin (km)", main="b)", col='black',
          levels = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
          labcex=0.8, xlim=c(60,200), ylim=c(0,18))
contour2D(x=dose1[,1], y=-dose1[,2], z=t(Z), lwd = 2, colkey = F, 
          clab = "", add=TRUE, col='black',
          levels = c(0.01,0.025,0.05), lty=3,
          labcex=0.8, xlim=c(60,200), ylim=c(0,18))

contour2D(x=dose1[,1], y=-dose1[,2], z=t(Z.upper), lwd = 2, colkey = F, 
          clab = "p(response)", xlab=RLlabel2, ylab="DISTmin (km)", main="c)", col='black',
          levels = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
          labcex=0.8, xlim=c(60,200), ylim=c(0,18))
contour2D(x=dose1[,1], y=-dose1[,2], z=t(Z.upper), lwd = 2, colkey = F, 
          clab = "", add=TRUE, col='black',
          levels = c(0.01,0.025,0.05), lty=3,
          labcex=0.8, xlim=c(60,200), ylim=c(0,18))



# # Export 2D function
predsurf<-Z
colnames(predsurf) = dose1[,1]
rownames(predsurf) = -dose1[,2]
write.csv(predsurf, file = "2D_function_median.csv")
predsurf<-Z.lower
colnames(predsurf) = dose1[,1]
rownames(predsurf) = -dose1[,2]
write.csv(predsurf, file = "2D_function_p025.csv")
predsurf<-Z.upper
colnames(predsurf) = dose1[,1]
rownames(predsurf) = -dose1[,2]
write.csv(predsurf, file = "2D_function_p975.csv")
predsurf<-Z.lower05
colnames(predsurf) = dose1[,1]
rownames(predsurf) = -dose1[,2]
write.csv(predsurf, file = "2D_function_p050.csv")
predsurf<-Z.upper95
colnames(predsurf) = dose1[,1]
rownames(predsurf) = -dose1[,2]
write.csv(predsurf, file = "2D_function_p950.csv")

