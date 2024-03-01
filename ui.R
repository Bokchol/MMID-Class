



shinyUI(

  pageWithSidebar(

    # title
    headerPanel("Interventions"),

    # Sidebar with sliders
    sidebarPanel(
      
sliderInput("beta","Beta",min=1,max=4,step=0.1,value=1.75),
sliderInput("lat_cov_i","LatrineCoverage",min=0,max=1,step=0.1,value=0.0),
sliderInput("lat_eff","Latrineefficacy",min=0,max=1,step=0.1,value=0.0),
sliderInput("vac_cov_i","Vaccinecoverage",min=0,max=1,step=.05,value=0.0),
sliderInput("prop_treat","ProportionTreated",min=0,max=1,step=.05,value=0.0),
sliderInput("infect_treat","InfectiousnessTreated",min=0,max=1,step=.05,value=1.0),
sliderInput("recovery_treat","RecoveryTreated",min=0,max=1,step=.05,value=0.2),
sliderInput("week_interv","InteventionStart",min=1,max=30,step=1,value=5)


      #add the quote
      # , br(), hr(),
      # tags$p("Throughout my academic career, I'd given some pretty good talks. But being considered the best speaker in the computer science department is like being known as the tallest of the Seven Dwarfs._(Randy Pausch)")
    ),


    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("graphs"),
      plotOutput("myImage")
    )
  )
)



