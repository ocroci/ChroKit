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
      
      # HTML("Choose a file...:<br>"),
      # shinyFilesButton('fileBAM', 'Choose file', 'Please select an enrichment file', FALSE,multiple=TRUE),
      # HTML("<br><br><br>"),
      # HTML("...or select the path:"),
      # textInput("BAMfrompath",NULL,value=NULL,placeholder = "/path/to/BAM.bam or bigWig.bw"),
      # actionButton("confirmImportBAMfrompath", "Open file")




    ),
    box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
      title=boxHelp(ID="msg_enrichmentFiles_deleteEnrichment",title="Opened enrichment files"),

      HTML("<b>Select file references to be deleted:</b>"),
      wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 450px; background-color: #ffffff;",
        checkboxGroupInput("selectedBAMtoRemove",NULL,NULL)
      ),
      HTML("<br><br>"),
      actionButton("deleteBAM", "Delete")
    ),
    box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
      title=boxHelp(ID="msg_enrichmentFiles_renameEnrichment",title="Rename enrichment file"),

      selectInput("selectedBAMfiletoRename",label="Select enrichment(s):",NULL),
      HTML("<br><br>"),
      textInput("newfilenameBAMfile","New enrichment name:",placeholder="type new enrichment name here",value=""),
      HTML("<br><br>"),
      actionButton("renameBAMfile", "Rename")
    )               


  )      


)

