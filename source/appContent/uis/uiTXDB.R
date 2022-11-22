tabTXDB <- tabItem(tabName = "TXDBblock",
  
  #tabsetPanel(
    
    #tabPanel("Manage databases",


      fluidRow(
        #box to select 
        box(width=9,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_databases_extractAnnotatedElements",title="Select an assembly"),

          HTML("<b>Choose assembly to use from those available on the system:</b>"),
          fluidRow(
            column(8,
              selectInput('searchASSEMBLYforuse', label=NULL,choices=as.list(availASSEMBLIES))
            ),
            column(4,
              HTML("")
              
            )
          ),
          HTML("<br>"),
          list(HTML("<b>Select bp upstream and downstream TSS/TES:</b>"),htmlhelp("","help_txdb_windowupdownstream")),
          

          #here, display slider inputs for TSS/TESS up/downstream
          #to be improved
          HTML("<br><br><br>"),

          fluidRow(
            column(width=2,
              HTML("<h3>TSS regions</h3>")
            ),
            column(width=2,
              numericInput(inputId = 'absoluteFilterUpstreamTSS',label="Upstream",min = 0, max = 0, step = 50,value=300)
            ),
            column(width=2,
              numericInput(inputId = 'absoluteFilterDownstreamTSS',label="Downstream",min = 0, max = 0, step = 50,value=300)
            ),
            column(width=6)
          ),          

          HTML("<br><br><br><br>"),

          fluidRow(
            column(width=2,
              HTML("<h3>TES regions</h3>")
            ),
            column(width=2,
              numericInput(inputId = 'absoluteFilterUpstreamTES',label="Upstream",min = 0, max = 0, step = 50,value=300)
            ),
            column(width=2,
              numericInput(inputId = 'absoluteFilterDownstreamTES',label="Downstream",min = 0, max = 0, step = 50,value=300)
            ),
            column(width=6)
          ),

          # HTML("<h3><b>TES</b> base pairs upstream and downstream:</h3>"),
          # numericInput("TESdefinition",NULL,value=300),
          # HTML("<br><br>"),
          # checkboxGroupButtons(
          #   inputId = "includeOnly", label = "Transcripts for annotation must contain (click):", 
          #   choices = c("ENTREZ", "SYMBOL", "ENSEMBL", "RefSeq"), 
          #   justified = TRUE, status = "light",
          #   checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))
          # ),
          
          HTML("<br><br>"),
          actionButton("confirmASSEMBLYforuse", "Load assembly!")
          
        ),           

        box(width=3,collapsible = TRUE,status = "primary",solidHeader = TRUE,

          title=boxHelp(ID="msg_databases_downloadDatabases",title="... or download an assembly"),
          
          HTML("<b>Choose missing assembly to download from bioconductor:</b>"),
          # fluidRow(
          #   column(8,
              selectInput('searchASSEMBLYfordownload',label=NULL,choices=as.list(missingASSEMBLIES)),
            # ),
            # column(4,
              actionButton("confirmASSEMBLYfordownload", "Download database!")
          #   )
          # )
        )

      )



    #)

  #)
)
