## Time-stamp: <20-02-2008 06:40:59 ste PC69>

################################################################
################################################################
## This code defines several class and methods for the FLUSH  ##
## algorithm.                                                 ##
##                                                            ##
## Derived from FLUSH package, version 1.1.0                  ##
##                                                            ##
## Author: Stefano Calza                                      ##
## First version: 12 Feb 2008                                 ##
##                                                            ##
################################################################
################################################################


## --- Class definition

setClass("RA",representation(data = "data.frame", ProbeID = "character", keep = "matrix",
                             fitted = "matrix",
                             fdr = "numeric", cdfName = "character"),
         prototype = prototype(data = data.frame(), ProbeID = character(), keep = matrix(),
           fitted = matrix(),
           fdr = numeric(), cdfName = character())
         ,package="FLUSH.LVS.bundle")

setClass("FLUSH",representation(exprs = "matrix", ProbeID = "character", RA = "data.frame",
                                fdr = "numeric",
                                onlyExprs = "logical", keep = "logical", cdfName = "character"),
         prototype = prototype(exprs = matrix(), ProbeID = character(),RA = data.frame(),
           fdr = numeric(), onlyExprs = logical(), keep = logical(), cdfName = character())
         ,package="FLUSH.LVS.bundle")


## --- Generics definition

setGeneric("Flush",function(object,RA, proportion = 0.6, delta = NULL, lambda = NULL,
                               df = 10, check = TRUE, ...)
           standardGeneric("Flush")
           ,package="FLUSH.LVS.bundle")

setGeneric("FlushSet",function(object,proportion = 0.6, delta = NULL, lambda = NULL,
                               df = 10,...)
           standardGeneric("FlushSet")
           ,package="FLUSH.LVS.bundle")


setGeneric("Flushit",function(object,RA, proportion=0.6, df=10, delta=NULL, lambda=NULL,
                              check = TRUE, ...)
           standardGeneric("Flushit")
           ,package="FLUSH.LVS.bundle")



## Deprecated: this is kept only for compatibility with FLUSH
setGeneric("fitRA",function(object, verbose = FALSE, bg.RA = c("none","imm","rma"))
           standardGeneric("fitRA")
           ,package="FLUSH.LVS.bundle")

setGeneric("RAplot",function(object,...)
           standardGeneric("RAplot")
           ,package="FLUSH.LVS.bundle")



## --- Methods definition

## For functions evolutions, see package FLUSH_1.1.0

## Deprecated: this is kept only for compatibility with FLUSH
setMethod("fitRA","AffyBatch",
          function (object, verbose = FALSE, bg.RA = c("none","imm","rma"))
          {
            .Deprecated("compute.RA",package="FLUSH.LVS.bundle",
                        msg="fitRA is deprecated and kept only for backward compatibility. Use compute.RA instead")
            ans <- compute.RA(object ,verbose = verbose, bg.RA = bg.RA)


            return(ans)
            
          })

## setMethod("Flush","ExpressionSet",
##           function(object, RA,  check = TRUE, onlyExprs = TRUE, which = 1)
##           {
##             if(!length(RA@keep))
##               stop("Need first use FlushSet on RA object")
            
##             if(which > ncol(RA@keep))
##               stop(paste("Arg which must be [0,",ncol(RA@keep),"]",sep=""))
               
##             keep <- RA@keep[,which]
               
##             if(check)
##               if(!all.equal(featureNames(object),RA@ProbeID))
##                 if(!all(feaureNames(object) %in% RA@ProbeID))
##                   stop("Object gene names don't correspond to RA Probes ID. Please check!")
##                 else
##                   {
##                     keep <- keep[match(featureNames(object),RA@ProbeID)]
##                     RA@data <- RA@data[match(featureNames(object),RA@ProbeID)]
##                   }

##             out <- exprs(object)[keep,]
##             probeID <- featureNames(object)[keep]
            
##             if(onlyExprs)
##               {
                
##                 ans <- new("ExpressionSet",exprs=out,phenoData = phenoData(object),
##                            annotation=annotation(object))
##                 description(ans) <- description(object)
##               }
##             else
##               {
##                 ares <- data.frame(RA@data[,c("logResSD","sqrtArrays")],fitted=RA@fitted)

##                 ans <- new("FLUSH",exprs = out, ProbeID = probeID, RA = ares,
##                            onlyExprs = onlyExprs,
##                            keep = keep, cdfName = RA@cdfName)
##               }

            
##             return(ans)
            
##           })


setMethod("Flush","ExpressionSet",
          function(object, RA,  proportion=0.6, delta=NULL, lambda=NULL,df=10, check = TRUE,
                   onlyExprs = TRUE, which = 1,...)
          {

            ## If RA doesn't come from FlushSet
            if(!length(RA@keep) > 1)
              {
                onlyExprs = TRUE
                RA <- FlushSet(RA, proportion=proportion, df=df, delta=delta, lambda=lambda,...)
              }

            if(which > ncol(RA@keep))
                  stop(paste("Arg which must be [0,",ncol(RA@keep),"]",sep=""))
                
            keep <- RA@keep[,which]
            
            if(check)
              if(!all.equal(featureNames(object),RA@ProbeID))
                if(!all(feaureNames(object) %in% RA@ProbeID))
                  stop("Object gene names don't correspond to RA Probes ID. Please check!")
                else
                  {
                    keep <- keep[match(featureNames(object),RA@ProbeID)]
                    RA@data <- RA@data[match(featureNames(object),RA@ProbeID)]
                  }

            out <- exprs(object)[keep,]
            probeID <- featureNames(object)[keep]
            
            if(onlyExprs)
              {
                
                ans <- new("ExpressionSet",exprs=out,phenoData = phenoData(object),
                           annotation=annotation(object))
                description(ans) <- description(object)
              }
            else
              {
                ares <- data.frame(RA@data[,c("logResSD","sqrtArrays")],fitted=RA@fitted)

                ans <- new("FLUSH",exprs = out, ProbeID = probeID, RA = ares,
                           onlyExprs = onlyExprs,
                           keep = keep, cdfName = RA@cdfName)
              }

            
            return(ans)
            
          })



setMethod("Flush","matrix",
          function(object, RA,  proportion=0.6, delta=NULL, lambda=NULL,df=10, check = FALSE,
                   which = 1,...)
          {
            

            ## If RA doesn't come from FlushSet
            if(!length(RA@keep))
              {
                onlyExprs = TRUE
                RA <- FlushSet(RA, proportion=proportion, df=df, delta=delta, lambda=lambda,...)
              }
            
            if(which > ncol(RA@keep))
              stop(paste("Arg which must be [0,",ncol(RA@keep),"]",sep=""))
            
            keep <- RA@keep[,which]
            
            
            if(check)
              warning("No check defined for object of class matrix")
            
            out <- object[keep,]
            return(out)
            
          })



## CHANGES: 29Aug06 specify better the args delta and lambda
##          29Aug06 allow to use more than 1 value for proportion
##          1May07 allow to specify the df of the B-spline

setMethod("FlushSet","RA",
          function(object, proportion = 0.6, delta = NULL, lambda = NULL, df = 10,...)
          {
            require(quantreg, quietly = TRUE)
            require(splines, quietly = TRUE)


            
            if(!is.null(delta) && is.null(lambda))
              lambda <- delta
            
            if(!is.null(lambda) && is.null(delta))
              delta <- lambda
              
            weights <- numeric()
            
            tau <- proportion


            
            logResSD <- object@data[["logResSD"]]
            Y <- object@data[["sqrtArrays"]]

            if(!is.null(delta))
              weights = (rank(logResSD)/length(logResSD)+delta)/lambda


            if(length(weights))
              fitted.values <- fitted(rq(Y ~ bs(logResSD,df=df),weights=weights,tau = tau))
            else
              fitted.values <- fitted(rq(Y ~ bs(logResSD,df=df),tau = tau))

            if(length(tau) > 1)
              {
                keep <- matrix(Y,ncol=ncol(fitted.values),nrow=nrow(fitted.values),
                               byrow=FALSE) > fitted.values
              }
            else
              {
                keep <- Y > fitted.values
                dim(keep) <- c(length(Y),1)
                dim(fitted.values) <- c(length(fitted.values),1)
              }
            ##keep <- drop(Y > fitted.values)

            ans <- object
            ans@keep <- keep
            ans@fitted <- fitted.values
            #ans@fitted <- drop(fitted.values)

            ans
            
            })



## setMethod("Flushit","ExpressionSet",
##           function(object,RA, proportion=0.6, df=10, delta=NULL, lambda=NULL,
##                               check = TRUE,...)
##           {


##             if(length(proportion) > 1)
##               {
##               warning("Only one value possible for \"proportion\" argument in Flushit. Set to first value:", proportion[1])
##               proportion <- proportion[1]
##             }
            
##             RA <- FlushSet(RA, proportion=proportion, df=df, delta=delta, lambda=lambda,...)
            
##             keep <- RA@keep
               
##             if(check)
##               if(!all.equal(featureNames(object),RA@ProbeID))
##                 if(!all(feaureNames(object) %in% RA@ProbeID))
##                   stop("Object gene names don't correspond to RA Probes ID. Please check!")
##                 else
##                   {
##                     keep <- keep[match(featureNames(object),RA@ProbeID)]
##                     RA@data <- RA@data[match(featureNames(object),RA@ProbeID)]
##                   }

##             out <- exprs(object)[keep,]
##             probeID <- featureNames(object)[keep]
            
##             ans <- new("ExpressionSet",exprs=out,phenoData = phenoData(object),
##                        annotation=annotation(object))
##             description(ans) <- description(object)

            
##             return(ans)
            
##           })
