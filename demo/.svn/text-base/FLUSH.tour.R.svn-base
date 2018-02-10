#####################################
##                                 ##
## DEMO for package FLUSH          ##
##                                 ##
## v3 20 Feb 2008  - Stefano Calza ##
##                                 ##
#####################################



load("FLUSH.RData")

choe.RA <- compute.RA(choe)


RAplot(choe.RA,add=list(FldChg < 1,FldChg==1,
                 FldChg > 1, FldChg >= 2, FldChg >= 3),
       add.col=list("yellow","green","blue","brown","red"))


choe.fSet <- FlushSet(choe.RA)

RAplot(choe.fSet,add=list(FldChg < 1,FldChg==1,
                     FldChg > 1, FldChg >= 2, FldChg >= 3),
         add.col=list("yellow","green","blue","brown","red"),add.rq=TRUE)


choe.fSet2 <- FlushSet(choe.RA,delta=.5)
RAplot(choe.fSet2,add.rq = TRUE)



choe.Flush <- Flush(choe.MAS5,choe.fSet, onlyExprs = FALSE)

RAplot(choe.Flush, heat = TRUE, add.rq = TRUE)

choe.exprs <- Flush(choe.MAS5,choe.RA,proportion=0.6)
dim(exprs(choe.exprs))
