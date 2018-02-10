"TDshow" <-
function (x, ...,  ngenes = 200, DEG = "DEG.pset",FDR="fdr.local",
            legend, lty, col, main, xlab, ylab, labout=FALSE,cex=1,shift,bg="white",
            out.sens)
{

  args = list(x, ...)
     
  k = length(args)
  if(!missing(out.sens))
    out <- vector("list",k)

  if (missing(lty))
      lty = 2:(k+1)
  else
    if(length(lty) < k)
      lty <- rep(lty,ceiling(k/length(lty)))
  
  if (missing(col))
    col = rep(par("fg"), k)

  if (missing(xlab)) 
    xlab = "Number Declared DE"

  if (missing(ylab)) 
    ylab = "Number of True DEG"
  
  xtract = function(x) {
    n = nrow(x)
    ord = order(x[[FDR]])
    TDE <- x[[DEG]]
    cTDE <- cumsum(TDE[ord])
    cTDE <- cTDE[1:ngenes]
    cnt = 1:n
    cnt = cnt[1:ngenes]

    list(count = cnt,cTDE = cTDE, maxTDE = max(cTDE))
  }
  
  pp = lapply(args, xtract)
  xr = yr = c(0,ngenes)

  if(labout & !missing(legend))
    {
      NCHAR <- strwidth(legend[which.max(nchar(legend))],"inches",cex=cex)
      opp <- par(omi=c(0,0,0,NCHAR),las=1,bg=bg)
    }
  
  plot(xr, yr, type = "n", xlab = xlab, ylab = ylab)
  for (i in 1:k) {
    ll = pp[[i]]
    lines(ll$count, ll$cTDE, lty = lty[i], col = col[i])
    if(!missing(out.sens))
      {
        if(out.sens[1]=="all")
          out[[i]] <- cbind(DEG=ll$count,sens=ll$cTDE/ll$count)
        else
          out[[i]] <- cbind(DEG=ll$count,sens=ll$cTDE/ll$count)[ll$count %in% out.sens,]
        }
    
  }
  abline(0,1)

  if(labout)
    {
      wAt <- sapply(pp,function(x) x[["maxTDE"]])

      if(!missing(shift))
        {
          sshift <- do.call("rbind",shift)
          ordwAt <- order(wAt,decreasing=T)[sshift[,1]]
          wAt[ordwAt] <- wAt[ordwAt] + sshift[,2]
        }
      
      mtext(legend,cex=cex,col=col,
            at = wAt
            ,side=4,padj=-.5,adj=-.1)
      
    }

  
  if (!missing(legend))
    {
      if(!labout)
        {
          legend = c("Line of Identity",legend)
          lty <- c(1,lty)
          col <- c("black",col)
          legend("topleft", legend = legend, lty = lty, col = col)
        }
      if(!missing(out.sens))
        names(out) <- legend
    }
      
  if (!missing(main)) 
    title(main)

  if(!missing(out.sens))
    invisible(out)
  
}

