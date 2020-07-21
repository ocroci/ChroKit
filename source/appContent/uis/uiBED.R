####################################################################################
# TAB BED files
####################################################################################

tabBED <- tabItem(tabName = "BEDblock",
   

  #tabsetPanel(
    #tabPanel("New coordinate files", 
    fluidRow(  
      column(width=8,
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_coordinateFiles_chooseCoordinates",title="From file"),

          column(width=4,
            HTML("<h3>File to import</h3><br>"),
            HTML("<h4>Parameters:</h4>"),
            checkboxInput("readheader","Header",value=TRUE),
            numericInput(inputId = 'skiplines',label="Lines to skip:",min = 0, step = 1,value=0),
            HTML("<br>"),
            HTML('<hr size=3>'),
            HTML("<h4>Open file:</h4>"),
            HTML("<b>Select a file...:</b><br>"),
            shinyFilesButton('file', 'Choose file', 'Please select a file', FALSE),
            HTML("<br><br>"),
            HTML("<b>... or choose a file path:</b>"),
            textInput("BEDfrompath",NULL,value="",placeholder = "/path/to/BEDorGTF"),
            #checkboxInput("readheaderpath","Header?",value=TRUE),
            #numericInput(inputId = 'skiplines2',label="Lines to skip:",min = 0, step = 1,value=0),
            actionButton("confirmImportBEDfrompath", "Open file")   
          ),
          column(width=8,
            HTML("<h3>File preview</h3>"),
            HTML("<br>"),
            htmlOutput("showcurrentfile"),
            dataTableOutput("fileHead"),
            HTML("<br><br>"),
            #open button and cancel button
            fluidRow(
              column(4,uiOutput('openfilebutton')),
              column(2,uiOutput('cancelfilebutton'))
            )          
          )
        ),
    
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                      title=boxHelp(ID="msg_genelists_importGenelist",title="From genelist"),
          column(width=6,  
            HTML("<h3>Genes to import</h3>"),
            HTML("<br>"),  
            HTML("<b>Open a text file...:</b><br>"),
            shinyFilesButton('fileGENELISTS', 'Choose gene list', 'Please select a txt file', FALSE),
            HTML("<br><br>"),
            HTML("<b>...or select a path...:</b>"),
            textInput("GENELISTSfrompath",NULL,value=NULL,placeholder = "/path/to/geneList.txt"),
            #files must be clean and not have the header. Gene list name is automatically given by the file name
            actionButton("createGENELISTSfrompath", "Open gene list"),
            HTML("<br><br>"),
            HTML("<b>...or put IDs/symbols here:</b><br>"),
            #wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
            textAreaInput("pastedGENELISTS",NULL,value="",height=150),
            textInput("nameGENELISTS",NULL,placeholder="genelist name",value=""),
            actionButton("createGENELISTSfrompaste", "Import")  
          ),

          column(width=6,
            HTML("<h3>Parameters</h3>"),
            HTML("<br>"),
            HTML("<br>"),
            radioButtons("symbolORid","What kind of identifiers are you importing?",choices=c(
                                                      "ENTREZ IDs"="entrez",
                                                      "ENSEMBL IDs"="ensembl",
                                                      "Symbols"="symbol",
                                                      "RefSeq IDs"="refseq"
                                                            ),selected="symbol"),
            HTML("<br>"),
            HTML("<b>Max length for transcripts:</b>"),
            numericInput(inputId = 'thresholdTranscripts',label=NULL,min = 0, step = 100000,value=200000)       
          )

        ) 
      
      
      ),

      column(width=4,
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_deleteRois_deleteRois",title="Delete ROIs"),

          HTML("<b>ROI to delete:</b>"),
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 300px; background-color: #ffffff;",
            checkboxGroupInput("selectedCustomROItoRemove",NULL,NULL)
          ),
          # HTML("<br>"),
          actionButton("deleteROI", "Delete")
        ),
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_deleteRois_renameRois",title="Rename ROIs"),

          # HTML("<br>"),
          HTML("<b>ROI to rename:</b>"),
          selectInput("selectedCustomROItoRename",label=NULL,NULL),
          # HTML("<br>"),
          textInput("newfilenameROI","New ROI name:",placeholder="type new ROI name here",value=""),
          # HTML("<br>"),
          actionButton("renameROI", "Rename")
        ),
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_deleteRois_reorderRois",title="Reorder ROIs"),

          wellPanel(id = "logPanelROI",style = "overflow-y:scroll; max-height: 300px",
            fluidPage(
              uiOutput("dinamicROI")
            )
          ),
          #button to confirm to reorder
          actionButton("reorderROI","Reorder!")
        )               
      )


        
       
        
    )

        


          
)
