## --- Wrapper funtion for C code interface

RA.fit <- function(X,y,nr,p,ln,start.p,n.b)
  {

    ## fit the model we need only betas, residuals and weights
    fit <- .C("lvs_rlm_fit_R",as.double(X),as.double(y),as.integer(nr),
              as.integer(p), double(p), double(ln),double(ln))[5:7]

    ## get se_estimates & residSE
    vfit <- unlist(.C("lvs_rlm_compute_se_R",as.double(X),as.double(y),
                      as.integer(nr), as.integer(p),as.double(fit[[1]]),as.double(fit[[2]]),
                      as.double(fit[[3]]),double(1),
                      double(1),as.integer(start.p),
                      as.integer(n.b))[c(8,9)])

    ## verify which is needed and then get only 2
    ##names(vfit) <- c("residSD","stddev","Array")

    return(vfit)
    
  }



## --- Create design matrix for fitting RLM models

make.design <- function(nsamples,nprobes)
  {
    Xprob <- matrix(rep(t(contr.helmert(nprobes)),nsamples),byrow=T,ncol=nprobes-1)
    Xsmp <- matrix(rep(contr.helmert(nsamples),each=nprobes),ncol=nsamples-1)
    ans <- cbind(rep(1,nrow(Xprob)),Xprob,Xsmp)
    ans
  }


## --- higer level function to create RA object

## TODO: add gcrma method for background correction

## abatch: object of class "AffyBatch"
## Value: returns a "numeric" matrix with 2 columns


setGeneric("compute.RA",function(abatch, verbose = FALSE, bg.RA = c("none","imm","rma"),
                                 do.parallel=TRUE,cores=NULL,use.snow=FALSE,
                                 type.cluster=c("MPI","PVM"))
           standardGeneric("compute.RA")
           ,package="FLUSH.LVS.bundle")

setMethod("compute.RA","AffyBatch",
          function (abatch,verbose=TRUE, bg.RA=c("none","imm","rma"),
                    do.parallel=TRUE,cores=NULL,use.snow=FALSE,
                    type.cluster=c("MPI","PVM"))
          {


            if(do.parallel || use.snow)
              {
                type.cluster <- match.arg(type.cluster)
                
                if(.Platform$OS.type == "unix")
                  {
                    if(use.snow)
                      {
                        if(is.null(cores))
                          cores <- 2
                        basicClusterInit(clusterNumberNodes=cores,typeCluster=type.cluster)
                      }
                    else
                      {
                        require(multicore)
                        if(!is.null(cores))
                          mc.cores <- cores
                        else
                          mc.cores <- getOption("cores")
                      }
                  }
                else
                  {
                    require(snow)
                    if(is.null(cores))
                      cores <- 2
                    
                    basicClusterInit(clusterNumberNodes=cores,typeCluster=type.cluster)
                  }
              }
            
            
            ProbeNames <- probeNames(abatch)
            GeneNames <- featureNames(abatch)
            nsamples <- nrow(pData(abatch))
            
            bgcorrect <- match.arg(bg.RA)
            
            ## modified from "rma.background.correct"
            
            bg.rma <- function (object, copy = FALSE)
              {

                x <- pm(object)
                rows <- dim(x)[1]
                cols <- dim(x)[2]

                if (is.integer(x)) {
                  x <- matrix(as.double(x), rows, cols)
                  copy <- FALSE
                }

                ans <- .Call("R_rma_bg_correct", x , copy, PACKAGE = "preprocessCore")

                return(ans)
              }
            
            
            
            ## Actual fitting function
            
            fit.fun <- function(ind)
              {
                mat <- PM[ind,]
                X <- make.design(nsamples=ncol(mat),nprobes=nrow(mat))
                start.p <- nrow(mat) + 1
                out <- RA.fit(X = X, y = mat, nr = nrow(X), p = ncol(X), ln = length(mat), 
                              start.p=start.p, n.b=n.b)
                out
                
              }
            
            PM <- switch(bgcorrect,
                         "none" = log2(pm(abatch)),
                         "imm" = log2(matrix(.C("IdealMM_correct", as.double(pm(abatch)), 
                           as.double(mm(abatch)), as.integer(length(ProbeNames)), 
                           as.integer(nsamples), as.character(ProbeNames))[[1]], 
                           ncol = nsamples)),
                         "rma" = log2(bg.rma(abatch)),
                         ##"gcrma" = log2(pm(bg.adjust.gcrma(abatch,type="affinities"))),
                         stop("Background correction method not implemented")
                         )
            
            
            n.b <- nsamples - 1
                        
            ind <- .Internal(split(seq_along(ProbeNames), factor(ProbeNames)))
            
            if(do.parallel)
              {
                if(use.snow)
                  {
                    out.lst <- parLapply(cl,ind,fit.fun)
                    stopCluster(cl)
                    
                  }
                else
                  out.lst <- mclapply(ind,fit.fun,mc.cores=mc.cores)
                
              }
            else
              out.lst <- lapply(ind,fit.fun)
            
            
            ## this would take time so better just split the index
            ## gene.list <- .local.split(PM,factor(ProbeNames))
            ## ans <- mclapply(gene.list,fit.fun) 
            
            out <- do.call("rbind",out.lst[GeneNames])
            colnames(out) <- c("logResSD","sqrtArrays")
            
            out[,"logResSD"] <- log(out[,"logResSD"])
            out[,"sqrtArrays"] <- sqrt(out[,"sqrtArrays"])
            out <- as.data.frame(out)
            
            ans <- new("RA",data = out, ProbeID = GeneNames, cdfName = cdfName(abatch))
            
            
            return(ans)
          }
          
          )
