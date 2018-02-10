## Modified from affyPLM 1.18.0 - Copyright Ben Bolstad 2008

## .initfunction <- function(where){

##   assign("normalize.ExpressionSet.methods",
##          "lvs",envir=as.environment(where))
## }



.First.lib <- function(libname, pkgname) {

  
  library.dynam("FLUSH.LVS.bundle",pkgname,libname)
  
  ## register lvs as a normalization method with the affyPLM package, if that is loaded:
  current.normmethods <- affy::normalize.AffyBatch.methods()
  
  upDate.normalize.AffyBatch.methods(c(current.normmethods,"lvs"))

  ## We require affyPLM for the moment, so no need to check for its presence
  ##  if ("package:affyPLM" %in% search())
  ##    {
  if(!"lvs" %in% .affyPLMInternalEnv[["normalize.ExpressionSet.methods"]])
    {
      assign("normalize.ExpressionSet.methods",
             c(.affyPLMInternalEnv[["normalize.ExpressionSet.methods"]],"lvs"),
             .affyPLMInternalEnv)
    }
  
  ##     }
  ##   else
  ##     {
  ##       .initfunction(match(paste("package:",pkgname,sep=""),search())) 
  ##     }
  
  .C("lvs_Lapack_Init",PACKAGE="FLUSH.LVS.bundle")

}



