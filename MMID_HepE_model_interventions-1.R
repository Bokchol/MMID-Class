#setwd("C:/...") # only needed if script is opened outside Rproj

library(deSolve)

cases <- read.csv('hepEdata_end.csv')

# set up a function to solve the model
HepE<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         # define intervention variables
         vac_rate <- vac_cov_i
         vaccinate <- (t>=(week_start+week_interv))*(t<week_start+week_interv+campaignweeks)*vac_rate
         lat_cov <- (t>=(week_start+week_interv))*lat_cov_i
         latrine <- (1-lat_eff*lat_cov)
         
         # define variables
         P <- (S+E+I+R+V)
         seas<-1+amp*cos(2*pi*(t-phi)/52)
         beta<-(R0*(muo+tau)*(gamma+muo))/gamma
         # include reduction in transmission due to latrines
         lam <- latrine*beta*seas*I/P
         
         
         # rate of change
         dS <- mui*P-muo*S-lam*S+omega*R-vaccinate*S+omega*V
         dE <- -muo*E+lam*S-gamma*E
         dI <- -muo*I+gamma*E-tau*I
         dR <- -muo*R+tau*I-omega*R
         dV <- -muo*V-omega*V+vaccinate*S
         
         # return the rate of change
         list(c(dS, dE, dI, dR, dV),new.inf=report*lam*S)
       }
  ) 
  
}

# MODEL INITIAL CONDITIONS
initP<-3000

initE<-1
initI<-0
initR<-0
initV<-0
initS<-initP-initE-initI-initR-initV
state <- c(S = initS, E=initE, I = initI,R = initR, V = initV)

#previously fitted parameter estimates
# estimate
#R0    6.076261 inital pop 2892.532196

#MODEL PARAMETERS
parameters <- c(mui=(1/(50*52)),    # birth
                muo=(1/(50*52)),    # death
                R0=6.1,               # basic reproduction number
                omega=(1/(10*52)),  # rate of loss of immunity = 1/(average duration of immunity)
                gamma=1/2,          # rate of movement from latent to infectious stage = 1/(average latent period)
                tau=1/4,            # rate of recovery = 1/(average duration of infection)
                report=1/7,         # proportion of all infections that are reported
                amp=0,              # relative amplitude of seasonal forcing
                phi=0,              # week of peak in seasonal forcing
                week_interv = 5,    # weeks after first case when the intervention starts
                lat_cov_i = 0.24,   # intervention coverage of latrines (baseline)
                lat_eff = 0.7,      # efficacy of latrines in reducing transmission
                vac_cov_i = 0,      # coverage of vaccine campaign
                vac_eff = 0,        # efficacy of vaccine 
                campaignweeks = 4   # time taken to reach target coverage
)

# Set the start and end time for the model simulation
week_start <- 38
week_stop <- 73
times <- seq(week_start, week_stop , by = 1)

#Numerical solving of ODEs
out <- ode(y = state, times = times, func = HepE, parms = parameters)
pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]+out[,"V"]
# define the actual week for the model output

inc <- out[,"new.inf"]
par(mfrow=c(2,1))

plot(cases$week,cases$new.case,pch=19,col='red',xlim=c(week_start,week_stop),ylim=c(0,120),main = "Outbreak cases: model versus data",xlab = "Time",ylab="new reported cases")
lines(out[,"time"],inc,lwd=3)

# The news report states that the final size of the epidemic was 1040
cumulativereports<- parameters["report"]*out[,"R"]
plot(tail(cases$week,n=1),1040,type='p',xlim=c(week_start,week_stop),ylim=c(0,2000),col='red',pch=19,main = "Cumulative case reports",xlab = "Time",ylab="cumulaive reported cases")
lines(out[,"time"],cumulativereports,lwd=3)     

