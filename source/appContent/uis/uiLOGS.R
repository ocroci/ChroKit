tabLOGS<-tabItem(tabName = "LOGSblock",
  fluidRow(
    HTML("<b><h3>Logs: </h3></b>"),
    wellPanel(id = "logPanel",style = "overflow-y:scroll; max-height: 250px",
        htmlOutput("showlogs")
    )
  )

)