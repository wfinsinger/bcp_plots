### Example 1 ####
# Change-point analysis of a sequence with a change in mean

#### Load libraries and R functions ####
library(changepoint)
library(bcp)
library(strucchange)

source("./R/bcp_plot_function.r")
source("./R/cpt_v12.r")



#### Create data ####
# Create a continuous dataset with a change in mean of "Counts" (e.g. counted charcoal pieces)
CmI <- seq(0, 200, 1)  # continuous sequence of sample depths
AgeI <- seq(0, 2000, 10) # continuous sequence of sample ages
prox.dat <- data.frame(CmI, AgeI)

sedArI <- diff(prox.dat$CmI) / diff(prox.dat$AgeI)  # calculate sediment-accumulation rates (here constant)
prox.dat$sedArI <- c(sedArI, NA)                    # add an NA in the last cell
prox.dat$countI <- c(rpois(100, 15), rpois(101, 3)) # draw counts from poisson distribution with given means
prox.dat$volI <- 1
prox.dat$ConcI <- prox.dat$countI / prox.dat$volI # calculate concentrations
prox.dat$ArI <- prox.dat$ConcI * prox.dat$sedArI  # calculate proxy-accumulation rate
prox <- prox.dat[-nrow(prox.dat), ] # delete last row because sedArI=NA in the last cell



#### Performs change-point analyses ####
# Change-point analysis of charcoal-accumulation rate with method="BinSeg"
ansmeanvar <- cpt.meanvar(prox$ArI, method="BinSeg")
plot(ansmeanvar, cpt.width=3)
print(ansmeanvar)

# Performs test to investigate influence of sediment-accumulation rate on proxy-accumulation rate
# NB: there is no influence because sediment-accumulation rates are constant in this example
cpt1 <- proxy.cpt(serie=prox, Name="Example_cpt_analysis", bootstrap=F)



### Performs bayesian change-point analyses ####
# Define dataset and series
bcp.data <- prox
series <- prox$ArI

# Perform bayesian change-point analysis (BCP)
bcp.0 <- bcp(series, return.mcmc=T)
#plot(bcp.0, main="Change Points: z.ages")

# Perform  breakpoint analysis (BP)
bp.0 <- breakpoints(series ~ 1, h = 2)$breakpoints


# Plot results
pdf("./cpt_output/Example_bcp-bp_analyses.pdf")
plot.bcp.proxy(params=bcp.data, series=series, bcp.run=bcp.0, bp.run=bp.0,
               age.scale="calBP", yr.min=0, yr.max=2000, yr.steps=100,
               title="bcp")
dev.off()

# Extract Posterior Probabilities from BCP analysis
bcp.0.postprob <- data.frame(bcp.data$AgeI, bcp.0$posterior.prob)
colnames(bcp.0.postprob) <- c("AgeI", "Posterior Probs")

interval.prob(bcp.0, 90, 110) # Estimates the probability of at least one change point
#                                  in the specified interval of sequential observations.
#                                 thus, 1000 cal yr BP

# Extract Breakpoints from BP analyses
bp.0.changes <- data.frame(bp.0, bcp.data$CmI[bp.0], bcp.data$AgeI[bp.0])
colnames(bp.0.changes) <- c("Position", "CmI", "AgeI")
