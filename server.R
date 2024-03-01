
library(deSolve)


week_start <- 38
week_stop <- 73
times <- seq(week_start, week_stop, by = (1/7))
parameters <- c(mui=(1/(50*52)),    
                muo=(1/(50*52)),    
                beta=1.76,              
                omega=(1/(10*52)),  
                gamma=1/2,          
                nui=1/4,            
                report=1/7,         
                amp=0,              
                phi=0,              
                nuitreat = 1/2, 
                gammatreat = 1/10,   
                week_interv = 5,
                lat_cov_i = 0.24,   
                lat_eff = 0.7,      
                vac_cov_i = 0.75,   
                prop_treat = 0.0,    
                infect_treat = 0.25, 
                campaignweeks = 4   
)

data<-read.csv(file="hepEdata_end.csv")

initP<-3700 
initE<-1 
initI<-0 
initR<-0 
initT<-0 
initS<-initP-initE-initI-initR+initT  

state <- c(S = initS, E=initE, I = initI,R = initR, TRT = initT)

HepE<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         P <- (S+E+I+R+TRT)
         seas<-1+amp*cos(2*pi*(t-phi)/52)
         R0<-beta*gamma/(muo+nui*(1-prop_treat)+prop_treat*nuitreat)*(gamma+muo)
         
         lat_cov <- (t>=(week_start+week_interv))*lat_cov_i
         latrine <- (1-lat_eff*lat_cov)
         lam <- latrine*beta*seas*(I/P+infect_treat*TRT/P)
         
         vac_rate <- (-log(1-vac_cov_i)/campaignweeks)
         vaccinate <- (t>=(week_start+week_interv))*(t<=week_start+week_interv+campaignweeks)*vac_rate
         
         dS <- mui*P-muo*S-lam*S+omega*R-vaccinate*S
         dE <- -muo*E+lam*S-gamma*E
         dI <- -muo*I - (1-prop_treat)*nui*I - prop_treat*nuitreat*I + gamma*E
         dR <- -muo*R - omega*R + vaccinate*S + gammatreat*TRT + (1-prop_treat)*nui*I
         dTRT <- -muo*TRT - gammatreat*TRT + prop_treat*nuitreat*I
         
         list(c(dS, dE, dI, dR, dTRT))
       }
  ) 
}




RunOde <- function(parms){

  parameters <- c(
  mui=(1/(50*52)),
  muo=(1/(50*52)),
  omega=(1/(10*52)),  
  gamma=1/2,
  nui=1/4,
  report=1/7,
  amp=0,
  phi=0,
  campaignweeks = 4,
  nuitreat = 1/2, 
  week_start = 38
, parms 

)

  
out <- ode(y = state, times = times, func = HepE, parms = parameters)



return(out)
}


shinyServer(

function(input, output, session) {


  parms <- reactive(c(
    beta= input$beta,
    lat_cov_i= input$lat_cov_i,
    lat_eff= input$lat_eff,
    vac_cov_i= input$vac_cov_i,
    prop_treat = input$prop_treat,    
    infect_treat = input$infect_treat,
    gammatreat = input$recovery_treat,
    week_interv = input$week_interv
    
  ))


  outode <- reactive(RunOde(parms()))


#plotting function
plotX <- function(){
  out <- outode()
  pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]+out[,"TRT"]
  inc <- parameters["report"]*parameters["gamma"]*out[,"E"]

  time<-out[,"time"]
  plot(time,inc,type='l',lwd=3,main = "Predicted HepE Incidence",xlab = "Time in weeks",ylab="New reported cases per week",ylim=c(0,85))
  points(data[,"week"],data[,"new.case"],pch=19,col='red')
}


  output$graphs <- renderPlot({
    plotX()
  })
  
  output$myImage <- renderImage({
    list(src = "model_interventions.bmp",height=500,width = 850)
  }, deleteFile = TRUE)

})






