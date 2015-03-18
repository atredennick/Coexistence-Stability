# Quick attempt for semi-discrete model coded
# by hand (i.e., not useing 'odeSolve')

updateDNR <- function(state_variables, parameters){
  D1 <- state_variables[1]
  D2 <- state_variables[2]
  N1 <- state_variables[3]
  N2 <- state_variables[4]
  R <- state_variables[5]
  mD <- parameters$mD
  r <- parameters$r
  K <- parameters$K
  K2 <- parameters$K2
  mN <- parameters$mN
  a <- parameters$a
  S <- parameters$S
  newD1 = D1 + (-(mD[1]*D1))
  newD2 = D2 + (-(mD[2]*D2))
  newN1 = N1 + (N1*(r[1]*exp(-K[1]*(exp(-K2[1]*R)))) - mN[1]*N1)
  newN2 = N2 + (N2*(r[2]*exp(-K[2]*(exp(-K2[2]*R)))) - mN[2]*N2)
  newR = max(0, R + (a*(S-R) - ((N1*(r[1]*exp(-K[1]*(exp(-K2[1]*R))))) + 
                           (N2*(r[2]*exp(-K[2]*(exp(-K2[2]*R))))))))
  return(c(newD1, newD2, newN1, newN2, newR))
}

# Test first model function

#TODO: make S a semi-normal function with increasing then decreasing values
#      over the 120 days.
S <- rlnorm(120,0,1)
S[50:120] <- 0
grow_days <- 120 #120 growing days

# Initial values
state_variables <- c(50, 50, 1, 1, 0)
save_states <- matrix(nrow=grow_days, ncol=length(state_variables))
save_states[1,] <- state_variables

# Run through one growing season
for(day in 2:grow_days){
  parameters <- list(mD = c(0.0001, 0.0001),
                     r = c(0.7, 0.3),
                     K = c(10,5),
                     K2 = c(1, 1.5),
                     mN = c(0.01, 0.01),
                     a = 0.5,
                     S = S[day])
  new_states <- updateDNR(state_variables, parameters)
  state_variables <- new_states
  save_states[day,] <- new_states
}

# Plot this to make sure it is right
par(mfrow=c(1,3))
matplot(c(1:grow_days), save_states[,c(3,4)], type="l",
        xlab="Day", ylab="Biomass")
plot(c(1:grow_days), save_states[,5], type="l",
     xlab="Day", ylab="Resource")

# plot the growth rate functions for each species
R <- seq(0,(4),0.01)
getR <- function(r, R, K, K2){
  out1 <- (r[1]*exp(-K[1]*(exp(-K2[1]*R))))
  out2 <- (r[2]*exp(-K[2]*(exp(-K2[2]*R)))) 
  return(cbind(out1, out2))
}
tmp<-getR(parameters$r, R, parameters$K, parameters$K2)
matplot(R, tmp, type="l", lty=c(1,1), 
        ylab="Instantaneous growth rate (r)", xlab="Resource density (R)")
