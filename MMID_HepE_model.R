

library(deSolve)

# set up a function to solve and fit the model
HepE<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         # define variables
         P <- (S+E+I+R)
         beta<-(R0*(muo+tau)*(gamma+muo))/gamma
         lam <- beta*I/P
         
         # rate of change
         dS <- mui*P-muo*S-lam*S+omega*R
         dE <- -muo*E+lam*S-gamma*E
         dI <- -muo*I+gamma*E-tau*I
         dR <- -muo*R+tau*I-omega*R
         
         # return the rate of change
         list(c(dS, dE, dI, dR))
       }
  ) 
}  

 # model initial conditions
  initP<-209100
  initE<-1
  initI<-0
  initR<-0
  initS<-initP-initE-initI-initR
  state <- c(S = initS, E=initE, I = initI,R = initR)
  
  # model parameters
  parameters <- c(mui=(1/(50*52)),    # birth
                  muo=(1/(50*52)),    # death
                  R0= 4,              # basic reproduction number
                  omega=(1/(10*52)),  # rate of loss of immunity = 1/(average duration of immunity)
                  gamma=(1/2),        # rate of movement from latent to infectious stage = 1/(average latent period)
                  tau=(1/4),          # rate of recovery = 1/(average duration of infection)
                  report=1/7         # proportion of all infections that are reported
  )
  
 
  # define the number of weeks to run the model
  times <- seq(0, 20, by = 1)
  # solve the ODEs
  out <- ode(y = state, times = times, func = HepE, parms = parameters)
  
  pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]  
  plot(out[,"time"],pop,type='l',xlab='week')
  plot(out)
#map with real data
  cases <- read.csv('hepEdata_begin.csv')
  week_start <- 38
  week_stop <- 78
  dev.new()
  plot(cases$week,cases$new.case,pch=19,col='red',xlim=c(week_start,week_stop),ylim=c(0,120),main = "Outbreak cases: model versus data",xlab = "Time",ylab="new reported cases")
#work out incidence
  beta  <- (parameters["R0"]*(parameters["muo"]+parameters["tau"])*(parameters["gamma"]+parameters["muo"]))/parameters["gamma"]
  lam <- beta*out[,"I"]/pop
  inc <- parameters["report"]*lam*out[,"S"]
#plot model output onto data  
  lines(out[,"time"]+week_start,inc,lwd=3)


  


  
   

  
