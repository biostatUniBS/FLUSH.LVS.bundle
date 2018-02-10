setGeneric("lvs.fit",function(object, proportion = 0.6, DF=10,nonpar=FALSE,...)
           standardGeneric("lvs.fit")
           ,package="FLUSH.LVS.bundle")






setMethod("lvs.fit","RA",
          function(object, proportion = 0.6, DF=10,nonpar=FALSE)
          {
            tau <- proportion           
            logResSD <- object@data[,"logResSD"]
            Y <- object@data[,"sqrtArrays"]
            fitted.values <- fitted(rq(Y ~ bs(logResSD,df=DF),tau = tau))
            keep <- Y < fitted.values
            keep
          }
          )

## For compatibility with LVS package
setMethod("lvs.fit","data.frame",
          function(object, proportion = 0.6, DF=10,nonpar=FALSE,log.it = TRUE, sqrt.it = log.it,
                   resid.name=ifelse(log.it,"residSD","logResSD"),
                   array.name=ifelse(sqrt.it,"Array","sqrtArrays"))
          {
            if(!all(c("residSD","Array") %in% names(object)))
              stop("object of class data.frame must contain columns \"residSD\" and \"Array\"")
            tau <- proportion
            
            if(log.it)
              logResSD <- log(object[,resid.name])
            else
              logResSD <- object[,resid.name]

            if(sqrt.it)
              Y <- sqrt(object[,array.name])
            else
              Y <- object[,array.name]
            
            fitted.values <- fitted(rq(Y ~ bs(logResSD,df=DF),tau = tau))
            keep <- Y < fitted.values
            keep
          }
          )



## lvs.fit <- function(object, proportion = 0.6, DF=10,nonpar=FALSE ,...)
##   {
##     tau <- proportion           
##     logResSD <- log(object[,"residSD"])
##     Y <- sqrt(object[,"Array"])
## ##     if(nonpar)
## ##       fitted.values <- fitted(rqss(Y ~ qss(logResSD),tau = tau))
## ##    else
##     fitted.values <- fitted(rq(Y ~ bs(logResSD,df=DF),tau = tau))
##     keep <- Y < fitted.values
##     keep
##   }




## Does all the job from AffyBatch --> Normalized expression

setGeneric("lvs",function(aBatch, proportion=0.6, DF=10, bg.RA=c("none","imm","rma"),
                          bgcorrect.method="mas", pmcorrect.method = "mas",
                          summary.method="mas")
           standardGeneric("lvs")
           ,package="FLUSH.LVS.bundle")


setMethod("lvs","AffyBatch",
          function(aBatch, proportion=0.6, DF = 10, bg.RA=c("none","imm","rma"),
                   bgcorrect.method="mas", pmcorrect.method = "mas",
                   summary.method="mas")
          {
            data.RA <- compute.RA(aBatch,bg.RA=bg.RA)
            lvs.id <- lvs.fit(data.RA, DF=DF, proportion=proportion)
            ##         cat("Computing expression values using ",summary.method," default settings\n")
            ##         cat("------------------------\n")
            lvs.prep <- expresso(aBatch, normalize=FALSE,
                                 bgcorrect.method=bgcorrect.method,
                                 pmcorrect.method=pmcorrect.method,
                                 summary.method=summary.method)
            cat("LVS Normalization...\n")
            normalize(lvs.prep,lvs.id=lvs.id,method="lvs")
          }
          )
