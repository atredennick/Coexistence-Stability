# Quick attempt for semi-discrete model coded
# by hand (i.e., not useing 'odeSolve')

updateDNR <- function(state_variables, parameters){
  D1 <- state_variables[1]
  D2 <- state_variables[2]
  N1 <- state_variables[3]
  N1 <- state_variables[4]
  newN = -(mD[1]*D1)
  dD2dt = -(mD[2]*D2)
  dN1dt = N1*(r[1]*exp(-K[1]*(exp(-K2[1]*R)))) - mN[1]*N1
  dN2dt = N2*(r[2]*exp(-K[2]*(exp(-K2[2]*R)))) - mN[2]*N2
  dRdt = a*(S-R) - ((N1*(r[1]*exp(-K[1]*(exp(-K2[1]*R))))) + (N2*(r[2]*exp(-K[2]*(exp(-K2[2]*R))))))
}

