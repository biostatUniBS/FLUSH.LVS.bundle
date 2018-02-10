## Time-stamp: <28-05-2007 14:33:33 ste PC69>

########################################################
## This code implements some plotting functions
##
## Author: Stefano Calza
## First version: 28 May 2007
##
## Note: radical library redesign from version 0.1.0
##
#######################################################



#############
## CHANGES ##
#############

## 29Aug06: corrected typo in getSYMBOL (ProbesID -> ProbeID)
## 10Spet06: added cex argument for added points
## 28May07: redesign library. No changes in this functions
##          logResSE -> logResSD

setMethod("RAplot","RA",
          function(object,cutfdr,breaks=5,colorize=FALSE,
                   add,add.col,fg="grey",where.legend="topright",legend.cex=1,
                   reverse=FALSE,col.seq,identify=FALSE,id.symbols=identify,add.cex=.3,
                   heat, heat.all = FALSE, nrgcols=13, add.rq = FALSE, FDR, nameResid,
                   nameArrays, add.pch = 16, ...)
          {

            
            sqrtArrays <- object@data[[ifelse(missing(nameArrays),"sqrtArrays",nameArrays)]]
            logResSD <- object@data[[ifelse(missing(nameResid),"logResSD",nameResid)]]
            
            fdr <- object@fdr
            
            if(!missing(heat))
              {
                
                rowObj <- object@ProbeID
                
                if(is(heat,"exprSet"))
                  {
                    wGn <- match(rowObj,geneNames(heat))
                    MuExpr <- esApply(heat[wGn,],1,mean)
                  }
                else
                  if(is(heat,"matrix"))
                    {
                      wGn <- match(rowObj,rownames(heat))
                      MuExpr <- apply(heat[wGn,],1,mean)
                    }
                  else
                    if(is(heat,"vector"))
                      MuExpr <- heat
                    else
                      stop("heat must be an object of class exprSet, matrix or vector")
                
                qntExpr <- quantile(MuExpr,seq(0,1,1/nrgcols))
                CMuExpr <- cut(MuExpr,qntExpr,right=T,include=T)
                levExpr <- 1:nlevels(CMuExpr)
                ColRGB <- rev(RGBColVec(nrgcols=nrgcols))
                Colr <- ColRGB[as.numeric(CMuExpr)]
                
                levExpr <- c(levExpr,max(levExpr)+1)

                par.orig <- par(no.readonly=TRUE)

                on.exit(par(par.orig))
        
                mar.orig <- par.orig$mar

                w <- (2 + mar.orig[2]) * par("csi") * 2.54
                ##layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
                layout(matrix(c(2,2,0,1), nc = 2), widths = c(1, lcm(w)))
                mar <- mar.orig
                mar[4] <- mar[2]
                ##mar[2] <- 1
                mar[2] <- 2
                par(mar = mar, las=1)
                plot.new()
                plot.window(xlim = c(0,1), ylim = range(levExpr), xaxs = "i", 
                            yaxs = "i")
                rect(0, levExpr[-length(levExpr)], 1, levExpr[-1], col = ColRGB)
                axis(4, at=levExpr, label=round(qntExpr,2))
                
                par(mar=par.orig$mar,las=par.orig$las)

                if(heat.all)
                  plot(y = sqrtArrays, x = logResSD,cex=.2,ylab="sqrt Array Effect",
                       xlab="log(residual Std Dev)",pch=16,col=Colr)
                else
                  plot(y = sqrtArrays, x = logResSD,cex=.2,ylab="sqrt Array Effect",
                       xlab="log(residual Std Dev)",pch=16,col=fg)
              }
            else
              plot(y = sqrtArrays, x = logResSD,cex=.2,ylab="sqrt Array Effect",
                   xlab="log(residual Std Dev)",pch=16,col=fg)
            
            
            if(!missing(cutfdr))
              {
                if(!length(object@fdr) & missing(FDR))
                  stop("Either slot fdr or argument FDR non missing needed.")

                fdr <- object@fdr

                if(!missing(FDR))
                  fdr <- FDR
                
                if(length(breaks) > 5 || (length(breaks)==1 && breaks > 5))
                  stop("Max 5 fdr classes")
                
                
                colTab <- data.frame(cex=c(2.5,2,1.5,1,0.5),pch=c(rep(16,5)),
                                     col = if(!missing(col.seq)) col.seq else
                                     I(c("red","blue","purple","darkgreen","skyblue")))
                
                if(!colorize)
                  colTab[,"col"] <- "black"
                
                if(length(breaks) > 1 & breaks[length(breaks)] > cutfdr)
                  breaks <- breaks[breaks <= cutfdr]
                
                objFDR <- cbind(object@data[fdr < cutfdr,c("sqrtArrays","logResSD")],
                                fdr[fdr < cutfdr] )
                
                
                BREAK <- hist(objFDR[["fdr"]],breaks=breaks,plot=FALSE)
                cutsb <- BREAK$breaks
                cnt <- BREAK$counts
                
                clFDR <- cut(objFDR[["fdr"]],cutsb,right=T,include=T)
                ## If not needed not use the cex=2.5 but start from 2
                if(nlevels(clFDR) <= 4)
                  colTab <- colTab[-1,]
                
                nclFDR <- as.numeric(clFDR)
                
                if(reverse)
                  {
                    nclFDR <- nlevels(clFDR) - nclFDR + 1
                    gCol <- nlevels(clFDR):1
                  }
                else
                  gCol <- 1:nlevels(clFDR)
                
                pCex <- colTab[nclFDR,"cex"]
                pPch <- colTab[nclFDR,"pch"]
                
                if(missing(heat))
                  Colfdr <- colTab[nclFDR,"col"]
                else
                  Colfdr <- Colr[fdr <= cutfdr]
                
                
                
                colPoints <- function(obj,Col)
                  {
                    points(y = obj[,"sqrtArrays"],x = obj[,"logResSD"],col=Col,cex=pCex,pch=pPch)
                  }
                
                colPoints(obj=objFDR,Col=Colfdr)
                
                legend(where.legend,legend=paste(levels(clFDR)," (",cnt,")",sep=""),
                       cex=legend.cex,
                       pch=colTab[gCol,"pch"],pt.cex=colTab[gCol,"cex"],col=colTab[gCol,"col"],
                       title="fdr (N genes)")
              }

            
            if(!missing(add))
              {
                if(is.list(add))
                  {
                    if(missing(add.col))
                      add.col <- as.list(rainbow(length(add)))

                    
                    if(length(add.col) == 1 && length(add) > 1)
                      add.col <- as.list(rep(unlist(add.col),length(add)))

                    if(length(add.pch) == 1)
                      add.pch <- as.list(rep(unlist(add.pch),length(add)))

                    if(length(add.cex) == 1)
                      add.cex <- as.list(rep(unlist(add.cex),length(add)))
                    
                    for(i in seq(along=add))
                      matpoints(x = logResSD[add[[i]]], y = sqrtArrays[add[[i]]],
                                col=add.col[[i]],pch=add.pch[[i]],cex=add.cex[[i]])
                  }
                else
                  {
                    if(missing(add.col))
                      add.col <- "red"
                    matpoints(x = logResSD, y = add,pch=add.pch,cex=add.cex,col=add.col)
                  }
              }
            
            if(add.rq)
              if(!length(object@fitted))
                stop("No fitted RQ values. Run FlushSet on RA object first")
              else
                {
                  matpoints(x = logResSD,y=object@fitted,pch=16,cex=.3,col="black")
                }
            
            if(identify)
              {
                cdfName <- tolower(object@cdfName)
                cdfName <- gsub("_","",cdfName)
                labels <- getSYMBOL(object@ProbeID,cdfName)
                
                idSymb <- identify(x = logResSD, y = sqrtArrays, labels=labels)
              }
          })



setMethod("RAplot","FLUSH",
          function(object,cutfdr,breaks=5,colorize=FALSE,
                   add,add.col,fg="grey",where.legend="topright",legend.cex=1,
                   reverse=FALSE,col.seq,identify=FALSE,id.symbols=identify,add.cex=.3,
                   heat=FALSE,nrgcols=13,heat.all, add.rq = FALSE, FDR, nameResid,
                   nameArrays, ...)
          {
            
            if(object@onlyExprs)
              stop("RAplot requires a FLUSH object with onlyExprs set to FALSE")
            
            
            sqrtArrays <- object@RA[[ifelse(missing(nameArrays),"sqrtArrays",nameArrays)]]
            logResSD <- object@RA[[ifelse(missing(nameResid),"logResSD",nameResid)]]

            fdr <- object@fdr

            nameArrays <- "sqrtArrays"
            nameResid <- "logResSD"


            if(heat)
              {
                heat <- object@exprs
                
                MuExpr <- apply(heat,1,mean)
                
                
                qntExpr <- quantile(MuExpr,seq(0,1,1/nrgcols))
                CMuExpr <- cut(MuExpr,qntExpr,right=T,include=T)
                levExpr <- 1:nlevels(CMuExpr)
                ColRGB <- rev(RGBColVec(nrgcols=nrgcols))
                Colr <- ColRGB[as.numeric(CMuExpr)]
                
                levExpr <- c(levExpr,max(levExpr)+1)

                par.orig <- par(no.readonly=TRUE)
                
                on.exit(par(par.orig))
                
                mar.orig <- par.orig$mar

                w <- (2 + mar.orig[2]) * par("csi") * 2.54
                ##layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
                layout(matrix(c(2,2,0,1), nc = 2), widths = c(1, lcm(w)))
                mar <- mar.orig
                mar[4] <- mar[2]
                ##mar[2] <- 1
                mar[2] <- 2
                par(mar = mar, las=1)
                plot.new()
                plot.window(xlim = c(0,1), ylim = range(levExpr), xaxs = "i", 
                            yaxs = "i")
                rect(0, levExpr[-length(levExpr)], 1, levExpr[-1], col = ColRGB)
                axis(4, at=levExpr, label=round(qntExpr,2))

                par(mar=par.orig$mar,las=par.orig$las)
                                
                plot(y = sqrtArrays, x = logResSD,cex=.2,ylab="sqrt Array Effect",
                     xlab="log(residual Std Dev)",pch=16,col=fg)
                points(y = sqrtArrays[object@keep], x = logResSD[object@keep],cex=.2,
                       pch=16,col=Colr)
              }
            else
              plot(y = sqrtArrays, x = logResSD,cex=.2,ylab="sqrt Array Effect",
                   xlab="log(residual Std Dev)",pch=16,col=fg)
            
            
            if(!missing(cutfdr))
              {
                if(!length(object@fdr) & missing(FDR))
                  stop("Either slot fdr or argument FDR non missing needed.")

                fdr <- object@fdr

                if(!missing(FDR))
                  fdr <- FDR
                
                if(length(breaks) > 5 || (length(breaks)==1 && breaks > 5))
                  stop("Max 5 fdr classes")
                
                
                colTab <- data.frame(cex=c(2.5,2,1.5,1,0.5),pch=c(rep(16,5)),
                                     col = if(!missing(col.seq)) col.seq else
                                     I(c("red","blue","purple","darkgreen","skyblue")))
                
                if(!colorize)
                  colTab[,"col"] <- "black"
                
                if(length(breaks) > 1 & breaks[length(breaks)] > cutfdr)
                  breaks <- breaks[breaks <= cutfdr]
                
                objFDR <- cbind(object@RA[fdr < cutfdr,c("sqrtArrays","logResSD")], fdr[fdr < cutfdr] )
                
                
                BREAK <- hist(objFDR[["fdr"]],breaks=breaks,plot=FALSE)
                cutsb <- BREAK$breaks
                cnt <- BREAK$counts
                
                clFDR <- cut(objFDR[["fdr"]],cutsb,right=TRUE,include=TRUE)
                ## If not needed not use the cex=2.5 but start from 2
                if(nlevels(clFDR) <= 4)
                  colTab <- colTab[-1,]
                
                nclFDR <- as.numeric(clFDR)
                
                if(reverse)
                  {
                    nclFDR <- nlevels(clFDR) - nclFDR + 1
                    gCol <- nlevels(clFDR):1
                  }
                else
                  gCol <- 1:nlevels(clFDR)
                
                pCex <- colTab[nclFDR,"cex"]
                pPch <- colTab[nclFDR,"pch"]
                
                if(!heat)
                  Colfdr <- colTab[nclFDR,"col"]
                else
                  Colfdr <- Colr[fdr <= cutfdr]
                
                
                
                colPoints <- function(obj,Col)
                  {
                    points(y = obj[,nameArrays],x = obj[,nameResid],col=Col,cex=pCex,pch=pPch)
                  }
                
                colPoints(obj=objFDR,Col=Colfdr)
                
                legend(where.legend,legend=paste(levels(clFDR)," (",cnt,")",sep=""),cex=legend.cex,
                       pch=colTab[gCol,"pch"],pt.cex=colTab[gCol,"cex"],col=colTab[gCol,"col"],
                       title="fdr (N genes)")
              }


                        if(!missing(add))
              {

                if(is.list(add))
                  {
                    if(missing(add.col))
                      add.col <- as.list(rainbow(length(add)))

                    if(length(add.col) == 1 && length(add) > 1)
                      add.col <- as.list(rep(unlist(add.col),length(add)))
                    
                    if(length(add.pch) == 1)
                      add.pch <- as.list(rep(unlist(add.pch),length(add)))
                    
                    if(length(add.cex) == 1)
                      add.cex <- as.list(rep(unlist(add.cex),length(add)))
                    
                    for(i in seq(along=add))
                      matpoints(x = logResSD[add[[i]]], y = sqrtArrays[add[[i]]],
                                col=add.col[[i]],pch=add.pch[[i]],cex=add.cex[[i]])
                  }
                else
                  {
                    if(missing(add.col))
                      add.col <- "red"
                    matpoints(x = logResSD, y = add,pch=add.pch,cex=add.cex,col=add.col)
                  }
              }
                        
            if(add.rq)
              if(!length(object@RA$fitted))
                stop("No fitted RQ values. Run FlushSet on RA object first")
              else
                {
                  matpoints(x = logResSD, y = object@RA$fitted,pch=16,cex=.3,col="black")
                }
                        
            if(identify)
              {
                labels <- getSYMBOL(object@ProbeID,object@cdfName)
                
                idSymb <- identify(x = logResSD, y = sqrtArrays ,labels=labels)
              }

          })
