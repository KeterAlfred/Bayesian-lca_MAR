cat("model{
  ## likelihood -----------------------------------------------------------------------------------------
  
  
  # Freqobs[1:J] ~ dmulti(p[1:J], N)
  
  nv[1:G] ~ dmulti(pv[1:G],N)

 for(j in 1:N){
 

    for (k in 1:K) {



    y[j, (k+1)] ~ dbern(p[j, k])  #Also imputes missing data for Ultra & culture based on provided priors


    p[j, k] <- (cp1[j, k,1]^d[j] * cp0[j, k,1]^(1 - d[j]))^(1-y[j,8]-y[j,9])*
               (cp1[j, k,2]^d[j] * cp0[j, k,2]^(1 - d[j]))^y[j,8] *
               (cp1[j, k,3]^d[j] * cp0[j, k,3]^(1 - d[j]))^y[j,9] 
    }
  
  
    d[j] ~ dbern( pD[j]  )


    
    
    
    
#Modeling pr(y_j = + | y_j-1 ... y_1, D=+)
  for(g in 1:G){
    cp1[j,1,g]  =  phi( inprod( y[j,c(1:1)], beta[1:1,1,g] ) )
    cp1[j,2,g]  =  phi( inprod( y[j,c(1:1)], beta[1:1,2,g] ) )
    cp1[j,3,g]  =  phi( inprod( y[j,c(1,3)], beta[1:2,3,g] ) )
    cp1[j,4,g]  =  phi( inprod( y[j,c(1:1)], beta[1:1,4,g] ) )
    cp1[j,5,g]  =  phi( inprod( y[j,c(1,5)], beta[1:2,5,g] ) )

  }
    
    




#Modeling pr(y_j = + | y_j-1 ... y_1, D=-)
for( g in 1:G){
    cp0[j,1,g]  =  phi( inprod( y[j,c(1:1)], alpha[1:1,1,g] ) )
    cp0[j,2,g]  =  phi( inprod( y[j,c(1:2)], alpha[1:2,2,g] ) )
    cp0[j,3,g]  =  phi( inprod( y[j,c(1:3)], alpha[1:3,3,g] ) )
    cp0[j,4,g]  =  phi( inprod( y[j,c(1:1)], alpha[1:1,4,g] ) )
    cp0[j,5,g]  =  phi( inprod( y[j,c(1:1)], alpha[1:1,5,g] ) )
}

    
    
    


#prob of TB 
pD[j] = phi( a_prev[1] + a_prev[2]*y[j,8] + a_prev[3] *y[j,9] )


#By group
    pij[j,1] <- phi( a_prev[1] )
    pij[j,2] <- phi( a_prev[1] + a_prev[2] )
    pij[j,3] <- phi( a_prev[1] + a_prev[3] )


    }





  
#Probability of disease by group
  pi[1] = phi( a_prev[1] )
  pi[2] = phi( a_prev[1] + a_prev[2] )
  pi[3] = phi( a_prev[1] + a_prev[3] )
  

  
  
#Overall probability of disease
  pi[4] = pi[1]*pv[1] + pi[2]*pv[2] + pi[3]*pv[3]
  
  
  
 
# Estimation of sensitivity & specificity
# for( g in 1:G){
#  se[1,g] = mean( cp1[,1,g] )
#  se[2,g] = mean( cp1[,2,g] )
#  se[3,g] = mean( cp1[,3,g] ) * se[2,g] +
#            mean( cp1[,3,g] ) * (1 - se[2,g] )
#  se[4,g] = mean( cp1[,4,g] )
#  se[5,g] = mean( cp1[,5,g] ) * se[4,g] +
#            mean( cp1[,5,g] ) * (1 - se[4,g] )
# 
#  sp[1,g] = 1 - mean( cp0[,1,g] )
#  sp[2,g] = 1 - ( mean( cp0[,2,g] ) * ( 1 - sp[1,g] ) +
#                  mean( cp0[,2,g] ) * sp[1,g] )
#  sp[3,g] = 1 - ( mean( cp0[,3,g] ) * ( 1 - sp[2,g] ) +
#                  mean( cp0[,3,g] ) *  sp[2,g] )
#  sp[4,g] = 1 - mean( cp0[,4,g] )
#  sp[5,g] = 1 - mean( cp0[,5,g] )
# 
# }



 for(j in 1:dim_D){

for( g in 1:G){
    s[j,1,g]  =  phi( inprod( D[j,c(1)],           beta[1:1,1,g] ) )
    s[j,2,g]  =  phi( inprod( D[j,c(1)],           beta[1:1,2,g] ) )
    s[j,3,g]  =  phi( inprod( D[j,c(1,3)],         beta[1:2,3,g] ) )
    s[j,4,g]  =  phi( inprod( D[j,c(1)],           beta[1:1,4,g] ) )
    s[j,5,g]  =  phi( inprod( D[j,c(1,5)],         beta[1:2,5,g] ) )



      #Modeling pr(y_j = + | y_j-1 ... y_1, D=-)
    c[j,1,g]  =  phi( inprod( D[j,c(1)],             alpha[1:1,1,g] ) )
    c[j,2,g]  =  phi( inprod( D[j,c(1,2)],           alpha[1:2,2,g] ) )
    c[j,3,g]  =  phi( inprod( D[j,c(1,2,3)],         alpha[1:3,3,g] ) )
    c[j,4,g]  =  phi( inprod( D[j,c(1)],             alpha[1:1,4,g] ) )
    c[j,5,g]  =  phi( inprod( D[j,c(1)],             alpha[1:1,5,g] ) )

    
    
    
  for( k in 1:K){
   ps[j,k,g] = s[j,k,g]*D[j,(k+1)] + (1 - s[j,k,g])*(1 - D[j,(k+1)])
   pc[j,k,g] = c[j,k,g]*D[j,(k+1)] + (1 - c[j,k,g])*(1 - D[j,(k+1)])
    }
    
    jp1[j,1,g] = ps[j,1,g]
    jp1[j,2,g] = ps[j,2,g]
    jp1[j,3,g] = prod( ps[j,2:3,g] )
    jp1[j,4,g] = ps[j,4,g]
    jp1[j,5,g] = prod( ps[j,4:5,g] )

    
    
    jp0[j,1,g] = pc[j,1,g]
    jp0[j,2,g] = prod( pc[j,1:2,g] )
    jp0[j,3,g] = prod( pc[j,1:3,g] )
    jp0[j,4,g] = pc[j,4,g]
    jp0[j,5,g] = pc[j,5,g]


}



}

for( g in 1:G){

    se[1,g] = sum( jp1[1:1,1,g] )
    se[2,g] = sum( jp1[1:1,2,g] )
    se[3,g] = sum( jp1[c(1,3),3,g] )
    se[4,g] = sum( jp1[1:1,4,g] )
    se[5,g] = sum( jp1[c(1,9),5,g] )

    
    
    sp[1,g] = 1 - sum( jp0[1:1,1,g] )
    sp[2,g] = 1 - sum( jp0[1:2,2,g] )
    sp[3,g] = 1 - sum( jp0[1:4,3,g] )
    sp[4,g] = 1 - sum( jp0[1:1,4,g] )
    sp[5,g] = 1 - sum( jp0[1:1,5,g] )


}



#Overall estimates

for( k in 1:K){
    se[k,4] =  (se[k,1]*pi[1]*pv[1] + se[k,2]*pi[2]*pv[2] + 
                                      se[k,3]*pi[3]*pv[3] ) / pi[4]
                  
    sp[k,4] =  (sp[k,1]*(1-pi[1])*pv[1] + sp[k,2]*(1-pi[2])*pv[2] + 
                                          sp[k,3]*(1-pi[3])*pv[3] ) / (1-pi[4])
}





  ## priors -------------------------------------------------------------
  
# TRUE POSITIVE CASES
  beta[1,1,1] ~ dnorm(0, 5)
  
  beta[1,2,1] ~ dnorm(0, 5)
  
  beta[1,3,1] ~ dnorm(0, 5)
  beta[2,3,1] ~ dnorm(0, 5)
  
  beta[1,4,1] ~ dnorm(0, 5)
  
  beta[1,5,1] ~ dnorm(0, 5)
  beta[2,5,1] ~ dnorm(0, 5)

#-------------------------------------------------------------
  beta[1,1,2] ~ dnorm(0, 5)
  beta[1,2,2] ~ dnorm(0, 5)
  
  beta[1,3,2] ~ dnorm(0, 5)
  beta[2,3,2] ~ dnorm(0, 5)
  
  beta[1,4,2] = beta[1,4,1] #~ dnorm(0, 1)
  
  beta[1,5,2] = beta[1,5,1] #~ dnorm(0, 1)
  beta[2,5,2] = beta[2,5,1] #~ dnorm(0, 1)


#-------------------------------------------------------------
  beta[1,1,3] ~ dnorm(0, 5)
  
  beta[1,2,3] ~ dnorm(0, 5)
  
  beta[1,3,3] ~ dnorm(0, 5)
  beta[2,3,3] ~ dnorm(0, 5)
  
  beta[1,4,3] = beta[1,4,1] #~ dnorm(0, 1)
  
  beta[1,5,3] = beta[1,5,1] #~ dnorm(0, 1)
  beta[2,5,3] = beta[2,5,1] #~ dnorm(0, 1)


#-------------------------------------------------------------

  
  
    
#TRUE NEGATIVE CASES
  alpha[1,1,1] ~ dnorm(0, 1)
  
  alpha[1,2,1] ~ dnorm(0, 1)
  alpha[2,2,1] ~ dnorm(0, 1)
  
  alpha[1,3,1] ~ dnorm(0, 1)
  alpha[2,3,1] ~ dnorm(0, 1)
  alpha[3,3,1] ~ dnorm(0, 1)
  
  alpha[1,4,1] ~ dnorm(-3, 10)
  
  alpha[1,5,1] ~ dnorm(-3, 10)

#-------------------------------------------------------------
  alpha[1,1,2] ~ dnorm(0, 1)
  
  alpha[1,2,2] ~ dnorm(0, 1)
  alpha[2,2,2] ~ dnorm(0, 1)
  
  alpha[1,3,2] ~ dnorm(0, 1)
  alpha[2,3,2] ~ dnorm(0, 1)
  alpha[3,3,2] ~ dnorm(0, 1)
  
  alpha[1,4,2] = alpha[1,4,1] #~ dnorm(-3, 10)
  
  alpha[1,5,2] = alpha[1,5,1] #~ dnorm(-3, 10)

#-------------------------------------------------------------  
  alpha[1,1,3] ~ dnorm(0, 1)
  
  alpha[1,2,3] ~ dnorm(0, 1)
  alpha[2,2,3] ~ dnorm(0, 1)
  
  alpha[1,3,3] ~ dnorm(0, 1)
  alpha[2,3,3] ~ dnorm(0, 1)
  alpha[3,3,3] ~ dnorm(0, 1)
  
  alpha[1,4,3] = alpha[1,4,1] #~ dnorm(-3, 10)
  
  alpha[1,5,3] = alpha[1,5,1] #~ dnorm(-3, 10)

  #-------------------------------------------------------------

  

  

  # Prior for prevalence
  a_prev[1] ~ dnorm( -3, 5)
  a_prev[2] ~ dnorm(  0, 5)
  a_prev[3] ~ dnorm(  -1, 5)

  
  

  
  
#Prior for group probability distribution
  pv[1:G] ~ ddirich(c(1,1,1))

  
}", file = "sim_model1_MAR_3_groups.txt")
#================================================================================================================


















# 
# 
# cat("model{
#   ## likelihood -----------------------------------------------------------------------------------------
# 
# 
#  for(j in 1:N){
#  
# 
#     for (k in 1:K) {
# 
# 
# 
#     y[j, (k+1)] ~ dbern(p[j, k])  #Also imputes missing data for Ultra & culture based on provided priors
# 
#       p[j, k] <- (cp1[j, k]^d[j] * cp0[j, k]^(1 - d[j]))
# 
#     }
#   
#   
#     d[j] ~ dbern( pD[j]  )
# 
# 
#     
#     
#     
#     
# #Modeling pr(y_j = + | y_j-1 ... y_1, D=+)
#     cp1[j,1]  =  phi( inprod( y[j,c(1,8,9)],   beta[1:3,1] ) )
#     cp1[j,2]  =  phi( inprod( y[j,c(1,8,9)],   beta[1:3,2] ) )
#     cp1[j,3]  =  phi( inprod( y[j,c(1,3,8,9)], beta[1:4,3] ) )
#     cp1[j,4]  =  phi( inprod( y[j,c(1)],       beta[1:1,4] ) )  #MAR assumption
#     cp1[j,5]  =  phi( inprod( y[j,c(1,5)],     beta[1:2,5] ) )  #MAR assumption
#     
#     # cp1[j,4]  =  phi( inprod( y[j,c(1,8,9)],   beta[1:3,4] ) )
#     # cp1[j,5]  =  phi( inprod( y[j,c(1,5,8,9)], beta[1:4,5] ) )
# 
#  
#     cp1x[j,1,1]  =  phi( inprod( y[j,c(1)],   beta[1:1,1] ) )
#     cp1x[j,2,1]  =  phi( inprod( y[j,c(1)],   beta[1:1,2] ) )
#     cp1x[j,3,1]  =  phi( inprod( y[j,c(1,3)], beta[1:2,3] ) )
#     cp1x[j,4,1]  =  phi( inprod( y[j,c(1)],   beta[1:1,4] ) )
#     cp1x[j,5,1]  =  phi( inprod( y[j,c(1,5)], beta[1:2,5] ) )
#     
#     cp1x[j,1,2]  =  phi( inprod( y[j,c(1)],   beta[1:1,1] ) + beta[2,1] )
#     cp1x[j,2,2]  =  phi( inprod( y[j,c(1)],   beta[1:1,2] ) + beta[2,2] )
#     cp1x[j,3,2]  =  phi( inprod( y[j,c(1,3)], beta[1:2,3] ) + beta[3,3] )
#     cp1x[j,4,2]  =  phi( inprod( y[j,c(1)],   beta[1:1,4] ) ) #+ beta[2,4] )
#     cp1x[j,5,2]  =  phi( inprod( y[j,c(1,5)], beta[1:2,5] ) ) #+ beta[3,5] )
#     
#     cp1x[j,1,3]  =  phi( inprod( y[j,c(1)],   beta[1:1,1] ) + beta[3,1] )
#     cp1x[j,2,3]  =  phi( inprod( y[j,c(1)],   beta[1:1,2] ) + beta[3,2] )
#     cp1x[j,3,3]  =  phi( inprod( y[j,c(1,3)], beta[1:2,3] ) + beta[4,3] )
#     cp1x[j,4,3]  =  phi( inprod( y[j,c(1)],   beta[1:1,4] ) ) #+ beta[3,4] )
#     cp1x[j,5,3]  =  phi( inprod( y[j,c(1,5)], beta[1:2,5] ) ) #+ beta[4,5] )    
#     
# 
# 
# #Modeling pr(y_j = + | y_j-1 ... y_1, D=-)
#     cp0[j,1]  =  phi( inprod( y[j,c(1,8,9)],   alpha[1:3,1] ) )
#     cp0[j,2]  =  phi( inprod( y[j,c(1:2,8,9)], alpha[1:4,2] ) )
#     cp0[j,3]  =  phi( inprod( y[j,c(1:3,8,9)], alpha[1:5,3] ) )
#     cp0[j,4]  =  phi( inprod( y[j,c(1)],       alpha[1:1,4] ) )  #MAR assumption
#     cp0[j,5]  =  phi( inprod( y[j,c(1)],       alpha[1:1,5] ) )  #MAR assumption
#     
#     # cp0[j,4]  =  phi( inprod( y[j,c(1,8,9)],   alpha[1:3,4] ) )
#     # cp0[j,5]  =  phi( inprod( y[j,c(1,8,9)],   alpha[1:3,5] ) )
#     
#     cp0x[j,1,1]  =  phi( inprod( y[j,c(1)],   alpha[1:1,1] ) )
#     cp0x[j,2,1]  =  phi( inprod( y[j,c(1:2)], alpha[1:2,2] ) )
#     cp0x[j,3,1]  =  phi( inprod( y[j,c(1:3)], alpha[1:3,3] ) )
#     cp0x[j,4,1]  =  phi( inprod( y[j,c(1)],   alpha[1:1,4] ) )
#     cp0x[j,5,1]  =  phi( inprod( y[j,c(1)],   alpha[1:1,5] ) )
#     
#     cp0x[j,1,2]  =  phi( inprod( y[j,c(1)],   alpha[1:1,1] ) + alpha[2,1] )
#     cp0x[j,2,2]  =  phi( inprod( y[j,c(1:2)], alpha[1:2,2] ) + alpha[3,2] )
#     cp0x[j,3,2]  =  phi( inprod( y[j,c(1:3)], alpha[1:3,3] ) + alpha[4,3] )
#     cp0x[j,4,2]  =  phi( inprod( y[j,c(1)],   alpha[1:1,4] ) ) #+ alpha[2,4] )
#     cp0x[j,5,2]  =  phi( inprod( y[j,c(1)],   alpha[1:1,5] ) ) #+ alpha[2,5] )    
#     
#     cp0x[j,1,3]  =  phi( inprod( y[j,c(1)],   alpha[1:1,1] ) + alpha[3,1] )
#     cp0x[j,2,3]  =  phi( inprod( y[j,c(1:2)], alpha[1:2,2] ) + alpha[4,2] )
#     cp0x[j,3,3]  =  phi( inprod( y[j,c(1:3)], alpha[1:3,3] ) + alpha[5,3] )
#     cp0x[j,4,3]  =  phi( inprod( y[j,c(1)],   alpha[1:1,4] ) ) #+ alpha[3,4] )
#     cp0x[j,5,3]  =  phi( inprod( y[j,c(1)],   alpha[1:1,5] ) ) #+ alpha[3,5] )    
#     
# 
# 
# #prob of TB 
# pD[j] = phi( a_prev[1] + a_prev[2]*y[j,8] + a_prev[3] *y[j,9] )
# 
# 
# #By group
#     pij[j,1] <- phi( a_prev[1] )
#     pij[j,2] <- phi( a_prev[1] + a_prev[2] )
#     pij[j,3] <- phi( a_prev[1] + a_prev[3] )
# 
# 
#     }
# 
# 
# 
# 
# 
#   
# #Probability of disease by group
#   pi[1] = phi( a_prev[1] )
#   pi[2] = phi( a_prev[1] + a_prev[2] )
#   pi[3] = phi( a_prev[1] + a_prev[3] )
#   
# 
#   
#   
# #Overall probability of disease
#   pi[4] = mean( pD[] )
#   
#   
#   
#  
# # Estimation of sensitivity & specificity
# for( g in 1:G){
#  se[1,g] = mean( cp1x[,1,g] )
#  se[2,g] = mean( cp1x[,2,g] )
#  se[3,g] = mean( cp1x[,3,g] ) * se[2,g] +
#            mean( cp1x[,3,g] ) * (1 - se[2,g] )
#  se[4,g] = mean( cp1x[,4,g] )
#  se[5,g] = mean( cp1x[,5,g] ) * se[4,g] +
#            mean( cp1x[,5,g] ) * (1 - se[4,g] )
# 
#  sp[1,g] = 1 - mean( cp0x[,1,g] )
#  sp[2,g] = 1 - ( mean( cp0x[,2,g] ) * ( 1 - sp[1,g] ) +
#                  mean( cp0x[,2,g] ) * sp[1,g] )
#  sp[3,g] = 1 - ( mean( cp0x[,3,g] ) * ( 1 - sp[2,g] ) +
#                  mean( cp0x[,3,g] ) *  sp[2,g] )
#  sp[4,g] = 1 - mean( cp0x[,4,g] )
#  sp[5,g] = 1 - mean( cp0x[,5,g] )
# 
# }
#  
#  
#  #Overall
#  se[1,4] = mean( cp1[,1] )
#  se[2,4] = mean( cp1[,2] )
#  se[3,4] = mean( cp1[,3] ) * se[2,4] +
#            mean( cp1[,3] ) * (1 - se[2,4] )
#  se[4,4] = mean( cp1[,4] )
#  se[5,4] = mean( cp1[,5] ) * se[4,4] +
#            mean( cp1[,5] ) * (1 - se[4,4] )
# 
#  sp[1,4] = 1 - mean( cp0[,1] )
#  sp[2,4] = 1 - ( mean( cp0[,2] ) * ( 1 - sp[1,4] ) +
#                  mean( cp0[,2] ) * sp[1,4] )
#  sp[3,4] = 1 - ( mean( cp0[,3] ) * ( 1 - sp[2,4] ) +
#                  mean( cp0[,3] ) *  sp[2,4] )
#  sp[4,4] = 1 - mean( cp0[,4] )
#  sp[5,4] = 1 - mean( cp0[,5] )
#  
#  
#  
# # for( k in 1:K){
# #     se[k,4] =  (se[k,1]*pi[1]*pv[1] + se[k,2]*pi[2]*pv[2] + 
# #                                       se[k,3]*pi[3]*pv[3] ) / pi[4]
# #                   
# #     sp[k,4] =  (sp[k,1]*(1-pi[1])*pv[1] + sp[k,2]*(1-pi[2])*pv[2] + 
# #                                           sp[k,3]*(1-pi[3])*pv[3] ) / (1-pi[4])
# # }
# 
# 
# 
# 
# 
#   ## priors -------------------------------------------------------------
#   
# # TRUE POSITIVE CASES
#   beta[1,1] ~ dnorm(0, 5)
#   beta[2,1] ~ dnorm(0, 5)
#   beta[3,1] ~ dnorm(0, 5)
#   
#   
#   beta[1,2] ~ dnorm(0, 5)
#   beta[2,2] ~ dnorm(0, 5)
#   beta[3,2] ~ dnorm(0, 5)
#   
#   
#   beta[1,3] ~ dnorm(0, 5)
#   beta[2,3] ~ dnorm(0, 5)
#   beta[3,3] ~ dnorm(0, 5)
#   beta[4,3] ~ dnorm(0, 5)
#   
#   beta[1,4] ~ dnorm(0, 5)
#   # beta[2,4] ~ dnorm(0, 5)
#   # beta[3,4] ~ dnorm(0, 5)
#   
#   
#   beta[1,5] ~ dnorm(0, 5)
#   beta[2,5] ~ dnorm(0, 5)
#   # beta[3,5] ~ dnorm(0, 5)
#   # beta[4,5] ~ dnorm(0, 5)
#   
# 
# #-------------------------------------------------------------
#   # beta[1,1,2] ~ dnorm(0, 5)
#   # beta[1,2,2] ~ dnorm(0, 5)
#   # 
#   # beta[1,3,2] ~ dnorm(0, 5)
#   # beta[2,3,2] ~ dnorm(0, 5)
#   # 
#   # beta[1,4,2] = beta[1,4,1] #~ dnorm(0, 1)
#   # 
#   # beta[1,5,2] = beta[1,5,1] #~ dnorm(0, 1)
#   # beta[2,5,2] = beta[2,5,1] #~ dnorm(0, 1)
# 
# 
# #-------------------------------------------------------------
#   # beta[1,1,3] ~ dnorm(0, 5)
#   # 
#   # beta[1,2,3] ~ dnorm(0, 5)
#   # 
#   # beta[1,3,3] ~ dnorm(0, 5)
#   # beta[2,3,3] ~ dnorm(0, 5)
#   # 
#   # beta[1,4,3] = beta[1,4,1] #~ dnorm(0, 1)
#   # 
#   # beta[1,5,3] = beta[1,5,1] #~ dnorm(0, 1)
#   # beta[2,5,3] = beta[2,5,1] #~ dnorm(0, 1)
# 
# 
# #-------------------------------------------------------------
# 
#   
#   
#     
# #TRUE NEGATIVE CASES
#   alpha[1,1] ~ dnorm(0, .1)
#   alpha[2,1] ~ dnorm(0, .1)
#   alpha[3,1] ~ dnorm(0, .1)
#   
#   
#   alpha[1,2] ~ dnorm(0, .1)
#   alpha[2,2] ~ dnorm(0, .1)
#   alpha[3,2] ~ dnorm(0, .1)
#   alpha[4,2] ~ dnorm(0, .1)
#   
#   
#   alpha[1,3] ~ dnorm(0, .1)
#   alpha[2,3] ~ dnorm(0, .1)
#   alpha[3,3] ~ dnorm(0, .1)
#   alpha[4,3] ~ dnorm(0, .1)
#   alpha[5,3] ~ dnorm(0, .1)
#   
#   
#   alpha[1,4] ~ dnorm(-3, 10)
#   # alpha[2,4] ~ dnorm(0, 10)
#   # alpha[3,4] ~ dnorm(0, 10)
#   
#   
#   alpha[1,5] ~ dnorm(-3, 10)
#   # alpha[2,5] ~ dnorm(0, 10)
#   # alpha[3,5] ~ dnorm(0, 10)
# 
# 
# #-------------------------------------------------------------
#   # alpha[1,1,2] ~ dnorm(0, .1)
#   # 
#   # alpha[1,2,2] ~ dnorm(0, .1)
#   # alpha[2,2,2] ~ dnorm(0, .1)
#   # 
#   # alpha[1,3,2] ~ dnorm(0, .1)
#   # alpha[2,3,2] ~ dnorm(0, .1)
#   # alpha[3,3,2] ~ dnorm(0, .1)
#   # 
#   # alpha[1,4,2] = alpha[1,4,1] #~ dnorm(-3, 10)
#   # 
#   # alpha[1,5,2] = alpha[1,5,1] #~ dnorm(-3, 10)
# 
# #-------------------------------------------------------------  
#   # alpha[1,1,3] ~ dnorm(0, .1)
#   # 
#   # alpha[1,2,3] ~ dnorm(0, .1)
#   # alpha[2,2,3] ~ dnorm(0, .1)
#   # 
#   # alpha[1,3,3] ~ dnorm(0, .1)
#   # alpha[2,3,3] ~ dnorm(0, .1)
#   # alpha[3,3,3] ~ dnorm(0, .1)
#   # 
#   # alpha[1,4,3] = alpha[1,4,1] #~ dnorm(-3, 10)
#   # 
#   # alpha[1,5,3] = alpha[1,5,1] #~ dnorm(-3, 10)
# 
#   #-------------------------------------------------------------
# 
#   
# 
#   
# 
#   # Prior for prevalence
#   a_prev[1] ~ dnorm( -3, 10)
#   a_prev[2] ~ dnorm(  0, 10)
#   a_prev[3] ~ dnorm(  0, 10)
# 
#   
#   
# 
#   
#   
# #Prior for group probability distribution
#   # pv[1:G] ~ ddirich(c(1,1,1))
# 
#   
# }", file = "sim_model1_MAR_3_groups.txt")
# #================================================================================================================
# 





