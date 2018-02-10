## Generic function

normalize.lvs <- function(object,ref.fun=c("median","mean"),lvs.id,use.loess=FALSE,...)
  {
    
    if(missing(lvs.id))
      stop("Must provide ids for LVS genes")
    
    if(!is.logical(lvs.id) && !is.numeric(lvs.id))
      stop("lvs.id must be a vector either logical or numeric")
    
    nk <- ncol(object)
    
    ref.fun <- match.arg(ref.fun)
    ref.data <- switch(ref.fun,"median"=rowMedians(object),rowMeans(object))
    
    out <- NULL
    for (i in 1:nk)
      {
        if(use.loess)
          {
            sm <- loess(object[lvs.id,i]~ref.data[lvs.id])
            a = approx(sm$fit, sm$x, xout=object[,i], rule = 2)$y
          }
        else
          {
            sm = smooth.spline(y=object[lvs.id,i] , x= ref.data[lvs.id])
            a = approx(sm$y, sm$x, xout=object[,i], rule = 2)$y
          }
        out = cbind(out,a)
      }
    
    colnames(out) <- colnames(object)
    rownames(out) <- rownames(object)
    
    return(out)
  }


normalize.ExpressionSet.lvs <- function(eset,ref.fun=c("median","mean"),lvs.id,
                                        use.loess=FALSE,...)
  {
    
    if(missing(lvs.id))
      stop("Must provide ids for LVS genes")
    
    if(!is.logical(lvs.id) && !is.numeric(lvs.id))
      stop("lvs.id must be a vector either logical or numeric")
    
    x <- exprs(eset)
    nk <- ncol(x)
    
    ref.fun <- match.arg(ref.fun)
    ref.data <- switch(ref.fun,"median"=rowMedians(x),rowMeans(x))
    
    out <- NULL
    for (i in 1:nk)
      {
        if(use.loess)
          {
            sm <- loess(x[lvs.id,i]~ref.data[lvs.id])
            a = approx(sm$fit, sm$x, xout=x[,i], rule = 2)$y
          }
        else
          {
            sm = smooth.spline(y=x[lvs.id,i] , x= ref.data[lvs.id])
            a = approx(sm$y, sm$x, xout=x[,i], rule = 2)$y
          }
        out = cbind(out,a)
      }
    
    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    
    exprs(eset) <- out

    return(eset)
  }



normalize.AffyBatch.lvs <- function(aBatch,lvs.id,ref.fun=c("median","mean"),
                                    use.loess=FALSE,
                                    bgcorrect.method="mas", pmcorrect.method = "mas",
                                    summary.method="mas",...)
  {
    
    if(missing(lvs.id))
      stop("Must provide ids for LVS genes")
    
    if(!is.logical(lvs.id) && !is.numeric(lvs.id))
      stop("lvs.id must be a vector either logical or numeric")

    eset <- expresso(aBatch, normalize=FALSE,
                         bgcorrect.method=bgcorrect.method,
                         pmcorrect.method=pmcorrect.method,
                         summary.method=summary.method)

    if(summary.method == "mas")
      x <- log2(exprs(eset))
    else
      x <- exprs(eset)
    
    nk <- ncol(x)
    
    ref.fun <- match.arg(ref.fun)
    ref.data <- switch(ref.fun,"median"=rowMedians(x),rowMeans(x))
    
    out <- NULL
    for (i in 1:nk)
      {
        if(use.loess)
          {
            sm <- loess(x[lvs.id,i]~ref.data[lvs.id])
            a = approx(sm$fit, sm$x, xout=x[,i], rule = 2)$y
          }
        else
          {
            sm = smooth.spline(y=x[lvs.id,i] , x= ref.data[lvs.id])
            a = approx(sm$y, sm$x, xout=x[,i], rule = 2)$y
          }
        out = cbind(out,a)
      }
    
    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    
    exprs(eset) <- out

    return(eset)
  }
