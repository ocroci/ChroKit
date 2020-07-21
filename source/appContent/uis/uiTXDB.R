tabTXDB <- tabItem(tabName = "TXDBblock",
  
  #tabsetPanel(
    
    #tabPanel("Manage databases",


      fluidRow(
        #box to select 
        box(width=7,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_databases_extractAnnotatedElements",title="Extract annotated elements from database"),

          HTML("Choose assembly to use:"),
          fluidRow(
            column(8,
              selectInput('searchASSEMBLYforuse', label=NULL,choices=as.list(availASSEMBLIES))
            ),
            column(4,
              HTML("")
              
            )
          ),
          HTML("<br>"),
          HTML("Select bp upstream and downstream TSS/TES:"),
          

          #here, display slider inputs for TSS/TESS up/downstream
          #to be improved
          HTML("<br><br><br>"),

          # div(style="display: inline-block;vertical-align:top; width: 300px;",sliderInput("upstreamTSS", "upstream:",min=-5000,max=0,value=-300,step=1)),
          # div(style="display: inline-block;vertical-align:top; width: 30px;",HTML("<br>")),
          # div(style="display: inline-block;vertical-align:top; width: 50px;",HTML("<h3><b>TSS</b></h3>")),
          # div(style="display: inline-block;vertical-align:top; width: 30px;",HTML("<br>")),
          # div(style="display: inline-block;vertical-align:top; width: 300px;",sliderInput("downstreamTSS", "downstream:",min=0,max=5000,value=300,step=1)),
          fluidRow(
            column(width=2,HTML("<b>Upstream</b>")),
            column(width=8),
            column(width=2,HTML("<b>Downstream</b>"))
          ),
          fluidRow(
            column(width=2,
              numericInput(inputId = 'absoluteFilterUpstreamTSS',label=NULL,min = 0, max = 0, step = 50,value=0)
            ),
            column(width=8,
              HTML("")
            ),
            column(width=2,
              numericInput(inputId = 'absoluteFilterDownstreamTSS',label=NULL,min = 0, max = 0, step = 50,value=0) 
            )
          ),          
          div(style="display: inline-block;vertical-align:top; width: 300px;",noUiSliderInput(inputId="upstreamTSS", label=NULL,min=0,max=5000,value=300,step=50,orientation="horizontal",direction="rtl",tooltips=FALSE)),
          div(style="display: inline-block;vertical-align:top; width: 30px;",HTML("<br>")),
          div(style="display: inline-block;vertical-align:top; width: 50px;",HTML("<h3><b>TSS</b></h3>")),
          div(style="display: inline-block;vertical-align:top; width: 30px;",HTML("<br>")),
          div(style="display: inline-block;vertical-align:top; width: 300px;",noUiSliderInput(inputId="downstreamTSS", label=NULL,min=0,max=5000,value=300,step=50,orientation="horizontal",direction="ltr",tooltips=FALSE)),


          HTML("<br><br><br><br>"),
          # div(style="display: inline-block;vertical-align:top; width: 300px;",sliderInput("upstreamTES", "upstream:",min=-5000,max=0,value=-300,step=1)),
          # div(style="display: inline-block;vertical-align:top; width: 30px;",HTML("<br>")),
          # div(style="display: inline-block;vertical-align:top; width: 50px;",HTML("<h3><b>TES</b></h3>")),
          # div(style="display: inline-block;vertical-align:top; width: 30px;",HTML("<br>")),
          # div(style="display: inline-block;vertical-align:top; width: 300px;",sliderInput("downstreamTES", "downstream:",min=0,max=5000,value=300,step=1)),
          
          fluidRow(
            column(width=2,HTML("<b>Upstream</b>")),
            column(width=8),
            column(width=2,HTML("<b>Downstream</b>"))
          ),
          fluidRow(
            column(width=2,
              numericInput(inputId = 'absoluteFilterUpstreamTES',label=NULL,min = 0, max = 0, step = 50,value=0)
            ),
            column(width=8,
              HTML("")
            ),
            column(width=2,
              numericInput(inputId = 'absoluteFilterDownstreamTES',label=NULL,min = 0, max = 0, step = 50,value=0) 
            )
          ),
          div(style="display: inline-block;vertical-align:top; width: 300px;",noUiSliderInput(inputId="upstreamTES", label=NULL,min=0,max=5000,value=300,step=50,orientation="horizontal",direction="rtl",tooltips=FALSE)),
          div(style="display: inline-block;vertical-align:top; width: 30px;",HTML("<br>")),
          div(style="display: inline-block;vertical-align:top; width: 50px;",HTML("<h3><b>TES</b></h3>")),
          div(style="display: inline-block;vertical-align:top; width: 30px;",HTML("<br>")),
          div(style="display: inline-block;vertical-align:top; width: 300px;",noUiSliderInput(inputId="downstreamTES", label=NULL,min=0,max=5000,value=300,step=50,orientation="horizontal",direction="ltr",tooltips=FALSE)),


          # HTML("<h3><b>TSS</b> base pairs upstream and downstream:</h3>"),
          # numericInput("TSSdefinition",NULL,value=300),
            

          # HTML("<br>"),

          # HTML("<h3><b>TES</b> base pairs upstream and downstream:</h3>"),
          # numericInput("TESdefinition",NULL,value=300),
          HTML("<br><br>"),
          checkboxGroupButtons(
            inputId = "includeOnly", label = "Transcripts for annotation must contain (click):", 
            choices = c("ENTREZ", "SYMBOL", "ENSEMBL", "RefSeq"), 
            justified = TRUE, status = "light",
            checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))
          ),
          
          HTML("<br><br>"),
          actionButton("confirmASSEMBLYforuse", "Load assembly!")
          
        ),           

        box(width=5,collapsible = TRUE,status = "primary",solidHeader = TRUE,

          title=boxHelp(ID="msg_databases_downloadDatabases",title="Download databases"),
          
          HTML("Choose missing assembly to download from bioconductor:"),
          fluidRow(
            column(8,
              selectInput('searchASSEMBLYfordownload',label=NULL,choices=as.list(missingASSEMBLIES))
            ),
            column(4,
              actionButton("confirmASSEMBLYfordownload", "Download database!")
            )
          )
        )

      )



    #)

  #)
)
