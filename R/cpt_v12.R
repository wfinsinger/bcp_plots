#----------------------------------------------------------------------------------------#
#  Determines zone boundaries for single influx records using the change-point analysis  #
#  as described in Finsinger et al. (2016).                                              #
#----------------------------------------------------------------------------------------#
#
# The function uses the output from the 'pretreatment_full.R' script (pretreat.full() function), thus 7 columns:
#    Column 1: CmI            =   top sample depth, interpolated
#    Column 2: AgeI           =   top sample age (as cal yrs BP), interpolated
#    Column 3: sedArI         =   sediment accumulation rate as from age-depth model (as cm year-1), interpolated
#    Column 4: countI         =   raw data (counts, areas, other), interpolated
#    Column 5: volI           =   sample volume, interpolated
#    Column 6: ConcI          =   concentration (as pieces cm-3), interpolated
#    Column 7: ArI            =   influx (or accumulation rate) value (as pieces cm-2 yr-1), interpolated

# The function also requires following additional parameters
#  (which can be left at default values)
#    Name      =   Site name
#    bootstrap =   if FALSE the random dataset is generated with the runif() function
#    q         =   number of random datasets (Default=1000) generated to determine
#    n.Q       =   the maximum number of change points, Default=10
#    n.screen  =   a change point in the random datasets is validated if it occurs in more than
#                   n.screen datasets. By default n.screen = q * 0.025 (thus with q=1000 this equals
#                   2.5% chance of occurrence)
#   output.dir =   the name of the output directory (Default = ./cpt_output)
#
#
# Suggested citation: Finsinger W., Magyari E.K., Fevre J., Orban I., Pal I., Vincze I., Hubay K,
#                     Birks H.H., Braun M., Toth M.  (2016) – Holocene fire regimes near the treeline
#                     in the Retezat Mts. (Southern Carpathians). Quaternary International.
#                     doi: 10.1016/j.quaint.2016.04.029. In press
#

# ------ Defines FUNCTION ------------
proxy.cpt = function(serie, Name, bootstrap=F, q=1000, n.Q=10, n.screen=q*0.025, penalty="4*log(n)",
                     output.dir=file.path(".","cpt_output"))
{
  
# ------------ SETUP ------------ #
  
  # -------- Load required libraries ----------- #
  require(changepoint)

  # -------- Create output directory
  if (dir.exists(output.dir) == T) {
    cat("NB: Output directory already exists. Files are overwritten or added.\n\n")
  } else {
  dir.create(output.dir)
  }
  
  
# ---------- Defines default parameters
#   bootstrap = F
#   q = 1000
#   n.Q = 10
#   n.screen = q*0.025
#   penalty = "4*log(n)"
  
# ---------- Settings giving relatively robust results with Lia and Brazi CHARc and 'runif' function
  cpt.test=cpt.meanvar
  meth.cpt="BinSeg"
  pen="Manual"
  pen.val= penalty
  t.stat="Normal"

# --------- First sets some values useful in the following
A <- serie
ConcI = A[ ,6]
sedArI = A[ ,3]
A.length = length(ConcI)
A.min = min(ConcI)
A.max = max(ConcI)


# --------- Then starts the Change-point analysis

## CPT analysis with bootstrapped Concentration Aset ####

if(bootstrap) {
  
  cpts.rand.constant = matrix() # stores cpts from rand series
  cpts.rand.model = matrix() # stores cpts from rand series with age-depth model from A
  
  ## Loop cpt analyses through rand Asets
  for (i in 1:q) {
    Conc.rand = sample(ConcI, size=A.length, replace=T)
    
    ### Merge rand Aset with A
    rand = cbind(A, Conc.rand)
    # plot(rand[,4], rand[,6], type="l", main=paste(Name, "- randised Conc"), xlab="A Conc", ylab="rand Conc")
    
    ## Calculate randAR with linear age-depth model and with records' age-depth model
    
    sed.AR.median = median(rand$sedArI)
    
    rand$ArI.rand.constant = rand$Conc.rand * sed.AR.median
    rand$ArI.rand.model = rand$Conc.rand * sedArI
    
    #plot(rand$age, rand$ArI.rand.constant, type="l")
    #lines(rand$age, rand$ArI.rand.model, col="red")
    #lines(rand$age, rand$ArI, col="blue")
    
    ### Run change-point analysis with rand AR Asets
    #cpt.AR.rand = cpt.mean(rand$ArI.rand, method=meth.cpt)
    cpt.AR.rand.constant = cpt.test(rand$ArI.rand.constant, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
    cpts.rand.constant = c(cpts.rand.constant, cpts(cpt.AR.rand.constant))
    
    #cpt.AR.rand = cpt.mean(rand$ArI.A.rand, method=meth.cpt)
    cpt.AR.rand.model = cpt.test(rand$ArI.rand.model, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
    cpts.rand.model = c(cpts.rand.model, cpts(cpt.AR.rand.model))
  }
  
  ### Run cpt analysis with true CHAR record
  cpt.AR = cpt.test(rand$ArI, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
  cpts(cpt.AR)
  
  # screens cpts from rand series
  cpts.rand.screen = data.frame(tabulate(cpts.rand.constant), 1.2*max(rand$ArI))
  colnames(cpts.rand.screen) = c("cpt.count", "X1")
  cpts.rand.screen$id = seq(1:(length(cpts.rand.screen[,1])))
  rand.screen = cpts.rand.screen[ which(cpts.rand.screen[,1]>n.screen), ]
  dummy = matrix(data=NA, nrow=1, ncol=3, dimnames=list("1", c("cpt.count", "X1", "id"))) # adds a 0 line to plot in case of no cpt occurrence
  rand.screen = rbind(rand.screen, dummy)
  
  cpts.A.rand.screen = data.frame(tabulate(cpts.rand.model), 1.2*max(rand$ArI))
  colnames(cpts.A.rand.screen) = c("cpt.count", "X1")
  cpts.A.rand.screen$id = seq(1:(length(cpts.A.rand.screen[,1])))
  A.rand.screen = cpts.A.rand.screen[ which(cpts.A.rand.screen[,1]>n.screen), ]
  dummy = matrix(data=NA, nrow=1, ncol=3, dimnames=list("1", c("cpt.count", "X1", "id"))) # adds a 0 line to plot in case of no cpt occurrence
  A.rand.screen = rbind(A.rand.screen, dummy)
  
  
  ## Diagnostic plot with true CHAR record and cpts
  pdf(file.path(output.dir, paste(Name,".pdf", sep="")))
  par(mfrow=c(2,1))
  par(mar = c(0.5,5,0.5,1))
  par(oma = c(3,1,2,1), cex=0.8)
  x.lim.cpt = rev(range(c(1, A.length)))
  x.lim.deptime = rev(range(c(min(rand$age), max(rand$age))))
  y.lim=c(0, 1.4*max(rand$ArI))
  plot(cpt.AR, type="s", ylab="CHAR", xaxt="n", xlim=x.lim.cpt, main=Name)
  abline(v=c(cpts(cpt.AR)), lwd=2, col="blue")
  par(new=T)
  plot(rand.screen$id, rand.screen$X1, ylab="", yaxt="n", xaxt="n", pch=1, xlim=x.lim.cpt, ylim=y.lim)
  par(new=T)
  plot(A.rand.screen$id, A.rand.screen$X1, ylab="", yaxt="n", xaxt="n", pch=1, xlim=x.lim.cpt, ylim=y.lim)
  plot(rand$age, rand$sedArI, type="l", ylab="Sediment-accumulation\nrate (cm/yr)", xlim=x.lim.deptime)
  mtext(("Age (cal yrs BP)"), side = 1, line = 2.5, cex=0.9)
  dev.off()
  
  } else {  
  
  
  ## CPT Analysis with random CharConcentration generated with the 'runif' function ####

  # First creates a large record of random charcoal concentration A
  rand.A = data.frame(seq(1:(100*A.length)), runif((100*A.length), A.min, A.max))
  colnames(rand.A) [2] = "Conc.rand"
  
  ## Loop cpt analyses through rand Asets
  cpts.rand.constant = matrix() # makes space to store cpts from rand series with constant SAR
  cpts.rand.model = matrix() # makes space to store cpts from rand series with age-depth model from A
  
  for (i in 1:q) {
    # Samples the rand.A with replacement
    Conc.rand = sample(rand.A$Conc.rand, size=A.length, replace=T)
    #plot(rand.A[,1], rand.A[,2], type="l", main="randomised_charA", xlab="Age", ylab="rand CHAR#")
    
    ### Merge rand Aset with A
    rand = cbind(A, Conc.rand)
    # plot(rand[,2], rand[,6], type="l", main=paste(Name, "- randomised Conc"), xlab="A Conc", ylab="rand Conc")
    
    ## Calculate randAR with linear age-depth model and with records' age-depth model
    sed.AR.median = median(sedArI)
    
    rand$ArI.rand.constant = rand$Conc.rand * sed.AR.median # random CHAR with constant SAR
    rand$ArI.rand.model = rand$Conc.rand * sedArI      # random CHAR with modelled SAR
    
#     par(mfrow=c(2,1))
#     plot(rand$age, rand$ArI.rand.constant, type="l")
#     lines(rand$age, rand$ArI.rand.model, col="red")
#     lines(rand$age, rand$ArI, col="blue")
#     plot(rand$age, rand$sedArI, type="l")
    
    ### Run change-point analysis with rand CHAR Asets
    cpt.AR.rand.constant = cpt.test(rand$ArI.rand.constant, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
    cpts(cpt.AR.rand.constant)
    cpts.rand.constant = c(cpts.rand.constant, cpts(cpt.AR.rand.constant))
    
    cpt.AR.rand.model = cpt.test(rand$ArI.rand.model, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
    cpts.rand.model = c(cpts.rand.model, cpts(cpt.AR.rand.model))
  }
  
  ### Run cpt analysis with true CHAR record
  cpt.AR = cpt.test(rand$ArI, method=meth.cpt, Q=n.Q, test.stat=t.stat, penalty=pen, pen.value=pen.val)
  cpts(cpt.AR)
  
  # screens cpts from both rand series
  cpts.rand.screen = data.frame(tabulate(cpts.rand.constant), 1.2*max(rand$ArI))
  colnames(cpts.rand.screen) = c("cpt.count", "Y1")
  cpts.rand.screen$id = seq(1:(length(cpts.rand.screen[,1])))
  rand.screen = cpts.rand.screen[ which(cpts.rand.screen[,1]>n.screen), ]
  dummy = matrix(data=NA, nrow=1, ncol=3, dimnames=list("1", c("cpt.count", "Y1", "id"))) # adds a 0 line to plot in case of no cpt occurrence
  rand.screen.const = rbind(rand.screen, dummy)
  
  cpts.A.rand.screen = data.frame(tabulate(cpts.rand.model), 1.2*max(rand$ArI))
  colnames(cpts.A.rand.screen) = c("cpt.count", "Y1")
  cpts.A.rand.screen$id = seq(1:(length(cpts.A.rand.screen[,1])))
  A.rand.screen = cpts.A.rand.screen[ which(cpts.A.rand.screen[,1]>n.screen), ]
  dummy = matrix(data=NA, nrow=1, ncol=3, dimnames=list("1", c("cpt.count", "Y1", "id"))) # adds a 0 line to plot in case of no cpt occurrence
  rand.screen.model = rbind(A.rand.screen, dummy)
  
  
  ## Diagnostic plot with true CHAR record and cpts
  pdf(file.path(output.dir, paste(Name,".pdf", sep="")))
  par(mfrow=c(3,1))
  par(mar = c(0.5,5,0.5,1))
  par(oma = c(5,1,1,1), cex=0.7)
  x.lim.cpt = rev(range(c(1, A.length)))
  x.lim.age = rev(range(c(min(rand$AgeI), max(rand$AgeI))))
  y.lim=c(0, 1.4*max(rand$ArI))
  plot(A$AgeI, A$ConcI, type="s", ylab="Concentration\ninterpolated", xaxt="n", xlim=x.lim.age)
  plot(cpt.AR, type="s", ylab="Accumulation rate\ninterpolated", xaxt="n", xlim=x.lim.cpt)
  #plot(A$AgeI, A$ArI, type="s")
  abline(v=c(cpts(cpt.AR)), lwd=2, col="blue")
  par(new=T)
  plot(rand.screen.const$id, rand.screen.const$Y1, ylab="", yaxt="n", xaxt="n", pch=1, xlim=x.lim.cpt, ylim=y.lim, col="red")
  par(new=T)
  plot(rand.screen.model$id, rand.screen.model$Y1, ylab="", yaxt="n", xaxt="n", pch=1, xlim=x.lim.cpt, ylim=y.lim, col="red")
  plot(rand$AgeI, rand$sedArI, xlim=x.lim.age, type="l", ylab="Sediment-accumulation\nrate (cm/yr)", xaxt="n")
  at = seq(from=0, to=max(x.lim.age), by=1000)
  axis(side=1, at=at, las=2, hadj=0.9)
  mtext(("Age (cal yrs BP)"), side = 1, line = 3.5, cex=0.8)
  mtext(paste(Name), side = 3, line = 29.5, cex=0.8)
  dev.off()
  }
  
  # Extract depths and ages of change-points detected in the interpolated proxy series
 Positions_cpts <- cpt.AR@cpts
 Depths_cpts <- rand$CmI[cpt.AR@cpts]
 Ages_cpts <- rand$AgeI[cpt.AR@cpts]

  # Return output
  output <- structure(list(Positions_cpts=Positions_cpts, Depths_cpts=Depths_cpts,
                    Ages_cpts=Ages_cpts))
  class(output) <- "ProxyCPT"
  return(output)
 }


# Welcome
cat("Hi there, welcome to Change-point analysis with influx records (v.1.2)\n")
cat(" \n")
cat("The function requires one input file with 7 columns:\n")
cat("  Column 1: CmI            =   top sample depth, interpolated\n")
cat("  Column 2: AgeI           =   top sample age (as cal yrs BP), interpolated\n")
cat("  Column 3: sedArI         =   sediment accumulation rate as from age-depth model (as cm year-1), interpolated\n")
cat("  Column 4: countI         =   raw A, interpolated\n")
cat("  Column 5: volI           =   sample volume, interpolated\n")
cat("  Column 6: ConcI          =   concentration (as pieces cm-3), interpolated\n")
cat("  Column 7: ArI            =   influx (or accumulation rate) value (as pieces cm-2 yr-1), interpolated\n")