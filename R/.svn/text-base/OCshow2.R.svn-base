
"OCshow2" <- 
  function (x, ..., global = TRUE, percentage = TRUE, top = 0.1, ngenes = 200,what=c("FDR","TD","pTD"),
            legend, lty, col, main, xlab, ylab, labout=FALSE,cex=1,bg="white",shift, ylim=NULL, xlim=NULL,
            padj=NA, adj=NA)
{
  what <- match.arg(what)
  args = list(x, ...)
  cl = sapply(args, function(x) class(x)[1])
  if (length(setdiff(cl, c("fdr1d.result", "fdr2d.result", 
                           "FDR.result"))) > 0) 
    stop("arguments must have class fdr.result or FDR.result")
  if (!global & "FDR.result" %in% cl) {
    warning("argument of class FDR.result - global=FALSE will be ignored")
    global = TRUE
  }
  k = length(args)
  if (missing(lty))
    if(what %in% c("TD","pTD"))
      lty = 2:(k+1)
    else
      lty = 1:k
  if (missing(col))
    if(what == "TD" & !labout)
      col = rep(par("fg"), k+1)
    else
      col = rep(par("fg"), k)

  if (missing(xlab)) 
    xlab = if (percentage) 
      "Proportion DE"
    else "Number of genes declared DE"

  if (missing(ylab)) 
    ylab = if (global)
      switch(what,"FDR"="FDR","TD"="Estimated number of true DE genes",
             "Estimated proportion of true DE genes")
    else "fdr"

  xtract = function(x) {
    n = nrow(x)
    is.FDR = class(x)[1] == "FDR.result"
    Fdr = if (is.FDR) 
      x$FDR
    else x$fdr
    ord = order(Fdr)
    Fdr = Fdr[ord]
    cnt = 1:n

    
    prc = cnt/n
    ndx = prc <= top

    if (global & !is.FDR) 
      Fdr = cumsum(Fdr)/(1:length(Fdr))

    Fdr <- switch(what,"TD" = cnt * (1 - Fdr),"pTD" = 1-Fdr, Fdr)

    
    if(!percentage)
      {
        Fdr = Fdr[1:ngenes]
        cnt = cnt[1:ngenes]
      }
    else
      {
        Fdr = Fdr[ndx]
        cnt = cnt[ndx]
      }
    prc = prc[ndx]
    
    list(count = cnt, perc = prc, Fdr = Fdr, maxFDR = max(Fdr), minFDR = min(Fdr))
  }
  pp = lapply(args, xtract)
  xr = range(unlist(sapply(pp, function(x) range(if (percentage) 
    x$perc
  else x$count))))
  yr = range(unlist(sapply(pp, function(x) range(x$Fdr))))

  if(labout & !missing(legend))
    {
      NCHAR <- strwidth(legend[which.max(nchar(legend))],"inches",cex=cex)
      par(omi=c(0,0,0,NCHAR),las=1,bg=bg)
    }
  
  plot(xr, yr, type = "n", xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim)
  for (i in 1:k) {
    ll = pp[[i]]
    lines(if (percentage) 
          ll$perc
    else ll$count, ll$Fdr, lty = lty[i], col = col[i])

    
  }

  if(labout)
    {
      if(what == "pTD")
          wAt <- sapply(pp,function(x) x[["minFDR"]])
      else
          wAt <- sapply(pp,function(x) x[["maxFDR"]])
      if(!missing(shift))
        {
          sshift <- do.call("rbind",shift)
          ordwAt <- sshift[,1]
          wAt[ordwAt] <- wAt[ordwAt] + sshift[,2]
        }
      mtext(legend,cex=cex,col=col,
            at = wAt
            ,side = 4, padj = padj, adj = adj)
      
    }


  
  if(what == "TD")
    abline(0,1)
  if (!missing(legend))
    {
      if(!labout)
        {
        if(what=="TD")
          {
            legend = c("Line of Identity",legend)
            lty <- c(1,lty)
            col <- c("black",col)
          }
        legend("topleft", legend = legend, lty = lty, col = col)
      }
      
    }
  if (!missing(main)) 
    title(main)
  
}
