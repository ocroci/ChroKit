tabBAM <- tabItem(tabName = "BAMblock",

  fluidRow(
    #box to select 
    box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
      title=boxHelp(ID="msg_enrichmentFiles_importEnrichment",title="Import enrichment file"),



      radioButtons("loadEnrichmentsource",NULL,choices=c(
                                          "Choose one or more files from filesystem"="filesystem",
                                          "Manually type the path of a file"="path"
                                  ),selected="filesystem"),      
      uiOutput("loadEnrichmentsource"),

    ),


    uiOutput("show_boxremovebam"),

    uiOutput("show_boxrenamebam")


  )      


)

