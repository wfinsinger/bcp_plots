plot.bcp.proxy <- function(params=NULL, series=NULL, bcp.run=NULL, bp.run=NULL,
                           age.scale=NULL, yr.min=NULL, yr.max=NULL, yr.steps=NULL, title=NULL) {
  
  #age.scale <- "CE"
  # yr.min <- 2000    # youngest x-axis tick mark label
  # yr.max <- -5500   # oldest x-axis tick mark label
  # yr.steps <- 500   # interval for x-axis tick marks
  # 
  # age.scale <- "calBP"
  # yr.min <- 0
  # yr.max <- 7500
  # yr.steps <- 500
  # 
  # age.scale <- NULL
  # yr.min <- NULL
  # yr.max <- NULL
  # yr.steps <- NULL
  
  
  # Set variables used for plotting
  if (is.null(age.scale) == T) {
    xax <- seq(from=min(params$AgeI), to=max(params$AgeI), length=length(series))
    x.lim <- c(max(xax), min(xax))
    x.lab <- axisTicks(usr=x.lim, log=F)
    x.title <- "Age (cal yr BP)"
  }
  
  if (is.null(age.scale) == F) {
    if (age.scale == "calBP") {
      xax <- seq(from=min(params$AgeI), to=max(params$AgeI), length=length(series))
      x.lim <- c(max(xax), min(xax))
      x.lab <- seq(from=yr.min, to=yr.max, by=yr.steps)
      x.title <- "Age (cal yr BP)"
    }
    if (age.scale == "CE") {
      xax <- 1950 - (seq(from=min(params$AgeI), to=max(params$AgeI), length=length(series)))
      x.lim <- c(min(xax), max(xax))
      x.lab <- seq(from=yr.max, to=yr.min, by=500)
      x.title <- "Age (BCE/CE)"
    }
  }
  
  rho <- rep(0, length(series))
  rho[bp.run] <- 1
  b.num <- 1 + c(0, cumsum(rho[1:(length(rho)-1)]))
  bp.mean <- unlist(lapply(split(series, b.num), mean))
  bp.ri <- rep(0, length(series))
  for (i in 1:length(bp.ri)) {
    bp.ri[i] <- bp.mean[b.num[i]]
  }
  op <- par(mfrow=c(2,1), col.lab="black", col.main="black")
  op2 <- par(mar=c(0,4,4,2), cex.axis=0.75)
  plot(x=xax, y=bcp.run$data[,2], type="l", lwd=0.5, xlim=x.lim, axes=F,
       xlab="", ylab="Posterior Mean", main=title)
  lines(x=xax, bcp.run$posterior.mean, lwd=2)
  lines(x=xax, bp.ri, col="blue")
  axis(1, at=x.lab, line=-0.1, labels=F)

  axis(2)
  par(op2)
  op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
  plot(x=xax, y=bcp.run$posterior.prob, axes=F, type="l", ylim=c(0,1), xlim=x.lim,
       xlab=x.title, ylab="Posterior Probability", lwd=1.5)
  axis(1, at=x.lab)

  for (i in 1:length(bp.ri)) {
    abline(v=xax[bp.run[i]], col="blue")
  }
  axis(2, yaxp=c(0, 0.9, 3))
  par(op3)
  par(op)
}
