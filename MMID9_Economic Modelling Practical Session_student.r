#
#
# Mathematical Modelling of Infectious Diseases
###############################################################
## ECONOMIC MODELLING IN R PRACTICAL SESSION I  
###############################################################
#by Nantasit Luangasanati


# Extend Hepatitis E Model to include economic evaluation of interventions; 
# Improving hygiene
# vaccaination



library(deSolve)
rm(list=ls())
data <- read.csv('hepEdata_end.csv')

# Set the start and end time for the model simulation
week_start <- 38
week_stop <- 78
times <- seq(week_start, week_stop, by = (1/7))

# MODEL INITIAL CONDITIONS
initP<-3700 # population size
initE<-1 # Exposed
initI<-0 # Infectious
initR<-0 # Immune
initS<-initP-initE-initI-initR # Susceptible (non-immune)

state <- c(S = initS, E=initE, I = initI,R = initR)
# set up a function to solve the equations
HepE<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         # define variables
         P <- (S+E+I+R)
         seas<-1+amp*cos(2*pi*(t-phi)/52)
         beta<-R0*(muo+tau)*(gamma+muo)/gamma
         
         lat_cov <- (t>=(week_start+week_interv))*lat_cov_i
         latrine <- (1-lat_eff*lat_cov)
         
         # simple vaccination rate
         vac_rate <- vac_cov_i
         vaccinate <- (t>=(week_start+week_interv))*(t<week_start+week_interv+campaignweeks)*vac_rate
         
         lam <- latrine*beta*seas*I/P
         
         # rate of change
         dS <- mui*P-muo*S-lam*S+omega*R-vac_eff*vaccinate*S
         dE <- -muo*E+lam*S-gamma*E
         dI <- -muo*I+gamma*E-tau*I
         dR <- -muo*R+tau*I-omega*R+vac_eff*vaccinate*S
         
         # return the rate of change
         list(c(dS, dE, dI, dR))
       }
  ) 
  
}

##########################
base_par <- c(mui=(1/(50*52)),    # birth
                muo=(1/(50*52)),    # death
                R0=6.6,               # basic reproduction number
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



############################
#Intervention 1: Latrines
############################
#MODEL PARAMETERS
lat_par<-base_par
lat_par['lat_cov_i'] <-  0.8 # COMPLETE: intervention coverage of latrines (*3 times from baseline) 
lat_par['lat_eff'] <- 0.7 # COMPLETE: effectiveness of latrines
lat_par

############################
#Intervention2: Vaccine
############################
vac_par<-base_par
vac_par['vac_cov_i'] <-0.2  # COMPLETE: Vaccination coverage
vac_par['vac_eff'] <- 0.3 # COMPLETE: Vaccination effectiveness
vac_par


# Run the model
#Runs of 3 sets of parameters and compare the results of all options
out_base <- ode(y = state, times = times, func = HepE, parms = base_par)
out_lat <- ode(y = state, times = times, func = HepE, parms = lat_par)
out_vac <- ode(y = state, times = times, func = HepE, parms = vac_par)

##############################
#Outputs of all interventions
##############################
#Total cases from model prediction

Inc_base<-base_par["report"]*base_par["gamma"]*out_base[,"E"]
Inc_lat<-lat_par["report"]*lat_par["gamma"]*out_lat[,"E"]
Inc_vac<-vac_par["report"]*vac_par["gamma"]*out_vac[,"E"]

Total_case_base <- sum(Inc_base)  #the outbreak finish at day 281
Total_case_lat <- sum(Inc_lat)
Total_case_vac <- sum(Inc_vac)

Total_case_base
Total_case_lat
Total_case_vac


#Death Outcomes (number of death)
risk_death<- 0.05  #case fatality of hepatitis E

Total_death_base <- Total_case_base*risk_death      #Total cases from model prediction
Total_death_lat <- Total_case_lat*risk_death
Total_death_vac <-Total_case_vac*risk_death

Total_death_base
Total_death_lat
Total_death_vac

############################
#Economic Evaluation
############################
#Costs
C_latrine <-     200      # COMPLETE: Cost of latrine per unit ($US)
Num_household <- initP/5  # NCOMPLETE: umber of household in the population (average 5 people per household) 
C_vaccine<-   350         # COMPLETE: Cost of vaccine per person ($US)
C_medical <-  30          # COMPLETE: Average medical cost of hepatitis E per case ($US)

Total_cost_base <- C_medical*Total_case_base
Total_cost_lat <- C_latrine*Num_household*lat_par[['lat_cov_i']]+C_medical*Total_case_lat
Total_cost_vac <- C_vaccine*initP*vac_par[['vac_cov_i']]+C_medical*Total_case_vac

#Check total cost of each strategy
Total_cost_base
Total_cost_lat 
Total_cost_vac


#Health outcomes
YLD.case<- 0.037  # COMPLETE: Life year with disability per case (disability weight at 0.3 for 45 days, {0.3*45/365})
YLL.death<- 30    # COMPLETE: Life years lost per death


DALYsloss_base <-Total_death_base*YLL.death+Total_case_base*YLD.case
DALYsloss_lat <- Total_death_lat*YLL.death+Total_case_lat*YLD.case
DALYsloss_vac <- Total_death_vac*YLL.death+Total_case_vac*YLD.case

DALYsloss_base
DALYsloss_lat
DALYsloss_vac

############################
#Cost-Effectiveness analysis
############################
# No intervention vs Latrine
Inc_cost_lat <- Total_cost_lat-Total_cost_base  #Incremental cost due to latrines
Inc_cost_vac <- Total_cost_vac-Total_cost_base  #Incremental cost due to vaccination

DALY_averted_lat <- DALYsloss_base-DALYsloss_lat       #Incremental benefit provided by latrines (Number of case averted)
DALY_averted_vac <- DALYsloss_base-DALYsloss_vac       #Incremental benefit provided by vaccine (Number of case averted)

ICER_lat<-Inc_cost_lat/DALY_averted_lat
ICER_vac<-Inc_cost_vac/DALY_averted_vac

ICER_lat
ICER_vac

All_res<-matrix(NA, ncol=3, nrow=2)
All_res[1,]<-c(Inc_cost_lat, DALY_averted_lat, ICER_lat)
All_res[2,]<-c(Inc_cost_vac, DALY_averted_vac, ICER_vac)
colnames(All_res)<-c("Incremental Cost", "DALYs_averted", "ICER")
rownames(All_res)<-c("Latrines", "Vaccine")
All_res


#Plot to see results on the ICER plane
plot(All_res[,2], All_res[,1], xlim=c(-1,1500), ylim=c(-100,800000), ylab="Incremental costs ($USD)", xlab="DALY loss", pch=18, col=c("red","blue"), main="CEA of Interventions for HepE")
text(All_res[,2], All_res[,1],label=c("Latrines","Vaccine"), cex=0.8, pos=3)
abline(a=0, b=1000, lty=3, col='black')
abline(v = 0, h = 0)

