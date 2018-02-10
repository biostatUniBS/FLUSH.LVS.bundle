basicClusterInit <- function(clusterNumberNodes = 2,
                             nameCluster = "cl",
                             typeCluster = c("MPI","PVM")) {
  ## Initialize a cluster
  ## Note that Rmpi already checks if another cluster is up and running,
  ## so we do not need to do that.

  require(snow)
  
  if(!(typeCluster %in% c("MPI", "PVM"))) stop("typeCluster needs to be PVM or MPI")

  if(typeCluster == "MPI") {
      require(Rmpi)
  }
  if(typeCluster == "PVM") {
      require(rpvm)
  }
  ## if(length(find(nameCluster)))
##       stop("\nThere is another object called ", nameCluster,".\n",
##            "This could mean that a cluster with that name already exists;\n",
##            "   in this case, please use the existing cluster \n",
##            "   ---you do not need to initialize the cluster, \n",
##            "   just pass its name as the parameter for 'nameCluster'---\n",
##            "   or stop that cluster and initialize a new one. \n",
##            "It could also mean that there is\n",
##            "   already an object with this name; either remove the object\n",
##            "   or use another name for the cluster.\n")

  assign(nameCluster,
         makeCluster(clusterNumberNodes,
                     type = typeCluster),
         env = .GlobalEnv)    

##   sprng.seed <- round(2^32 * runif(1))
##   print(paste("Using as seed for SPRNG", sprng.seed))
##   clusterSetupSPRNG(eval(parse(text = nameCluster)), seed = sprng.seed)
  clusterEvalQ(eval(parse(text = nameCluster)), library(FLUSH.LVS.bundle))
}
