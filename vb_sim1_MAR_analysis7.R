#===========================================================================
#
#
# Simultaneous imputation of the missing bacteriological test results in the 
# analysis model under the assumption that the missing bacteriological test 
# results are missing at random (MAR)
#===========================================================================

rm(list = ls())  #Clear the memory of any preloaded objects


#Point in the working directory
project.path = "/kyukon/scratch/gent/438/vsc43892/vb_sim_analysis"





#Read the data
load(file.path(project.path = getwd(),"sim_dat_Final1.RData"))



#Load the packages
# install.packages("doSNOW", repos = 'http://R-Forge.R-project.org')
dynamic.require <- function( package ) {
  
  if ( eval( parse( text = paste( 'require(', package, ')' ) ) ) ) {
    
    return( 'done' )
    
  }
  
  install.packages( package )
  
  return( eval( parse( text = paste( 'require(', package, ')' ) ) ) )
  
}

dynamic.require( 'doParallel' )

dynamic.require( 'rjags' )

dynamic.require( 'random' )




#Number of datasets
M=10
j=60  # will consider datasets c(1 : M)+j in sim_dat created in the data_for_sim_model1_MAR_4_groups.R code


#Source the dataset
source("data_for_sim1_model1_3_groups.R")

#Source the model
source("sim_model1_MAR_3_groups.R")





set.seed(79538)
jags.inits.1 <- function( ) {
  
  return( list( beta = rep(rnorm(0,10),3), alpha = rep(rnorm(0,10),3), a_prev = rep(rnorm(0,10),3),
                .RNG.name = 'lecuyer::RngStream',
                .RNG.seed = round( 1e+06 * runif( 1 ) ) ) )
  
}



# Multi-threaded run
# Find out how many cores are available (if you don't already know)
n = detectCores()
# n
# Find out how many cores are being used
# getDoParWorkers()
# n = M*3
cl <- makeCluster( n-1 ) # change to makeCluster( n ) to use n processes
registerDoParallel( cl )








set.seed( 74123 )

time.1 <- system.time( #getDoParWorkers( )
  vb.sim1.MAR.res61_70 <- foreach( m = 1:M, .inorder = FALSE, 
                                 .packages = c('R2jags'),#.packages = c( 'rjags','R2jags', 'random','coda','dclone' ), 
                                 .multicombine = TRUE) %dopar% {
                                   
          # load.module( 'lecuyer' )
          print(paste("Fitting model",m, sep = ""))
                                   
          result = jags.parallel(data = sim1_dat[[m]], inits = jags.inits.1,
           parameters.to.save = c("beta","alpha","a_prev",
                                  "se","sp","pi","pv"),
           model.file = 'sim_model1_MAR_3_groups.txt', n.chain = 3, 
           n.iter = 20000, n.burnin = 10000,
            n.thin = 10, DIC = T, jags.module = "dic", n.cluster = 3)
           return( result )
                                   
                                 }
)
print(time.1)



stopCluster(cl)


saveRDS(vb.sim1.MAR.res61_70, file = "/kyukon/data/gent/438/vsc43892/vb_sim_analysis/vb.sim1.MAR.res61_70.RData")


# print(vb.sim.MAR.res1[1], intervals = c(0.5, 0.025, 0.975))


