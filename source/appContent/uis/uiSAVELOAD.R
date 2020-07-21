

tabSAVELOAD<-tabItem (tabName = "SAVELOADblock",
  tabsetPanel(


    tabPanel("Save/Load",
      fluidRow(
        box(width=6,
          shinySaveButton("saveWork","Save","Save rds working environment...",filetype=list(rds="rds"))
        ),
        box(width=6,
          shinyFilesButton('loadenv', 'Load', 'Select rds file to load', FALSE)

        )  
      ) 
    )




  )         
)
