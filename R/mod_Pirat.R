
#' @title ## To be customized ##
#' 
#' @description 
#' ## To be customized ##
#' 
#' @param id The id of the module
#' @param obj An list of two items
#' @param reset xxx
#' @param verbose A boolean (FALSE as default) which indicates whether to 
#' display more details on the process
#' 
#' @name mod_Pirat
#' 
#' @return A shiny app
#'
NULL

#' @rdname mod_Pirat
#' @return A shiny app
#' @export
#' 
mod_Pirat_ui <- function(id){
  
  ns <- NS(id)
  
  tagList(
    uiOutput(ns('extension_ui')),
    shiny::actionButton(ns("run"), "Run Pirat"),
    uiOutput(ns('valid_btn_ui')),
    plotOutput(ns('correlation_plot_UI'))
  )
}


#' @rdname mod_Pirat
#' @export
#' @return A shiny app
#' 
mod_Pirat_server <- function(id,
                             obj = reactive({NULL}),
                             reset = reactive({NULL}),
                             verbose = FALSE
                             ) {
  
  
  
  # Define default selected values for widgets
  # This is only for simple workflows
  # This list must contain one item per widget (defined in the ui() function
  # The name of each item is the same as in the ui minus the suffix '_ui')
  widgets.default.values <- list(
    extension = 'base'
   # widget2 = NULL,
   # widget3 = NULL
  )
  
  moduleServer(id,function(input, output, session) {
    ns <- session$ns
    
   maxval <- reactiveVal(0)
   
    dataOut <- reactiveValues(
      trigger = NULL,
      value = NULL,
      widgets = NULL
    )
    
    output$extension_ui <- renderUI({
      selectInput(ns('extension'), 'Algorithm',
                  choices = c('base' = 'base', 
                              '2 pg' = '2',
                              'S' = 'S',
                              'T' = 'T'),
                  selected = widgets.default.values$extension,
                  width = '150px')
    })
    
    
    GetNbPg <- reactive({
      5
    })
    
    shiny::observeEvent(input$run, {
      
      maxval(GetNbPg())
      shiny::withProgress(
        withCallingHandlers(
         # out <- long_run_op(num_iter=10),
           dataOut$value <- my_pipeline_llkimpute(obj(),
                                               extension = input$extension,
                                               verbose = verbose),
           dataOut$trigger <- as.numeric(Sys.time()),
           dataOut$widgets <- list(extension = input$extension),
           
           
          message = function(m){
            msg <- unlist(strsplit(m$message, split=' '))
            if(msg[1] == 'Peptide_group') {
              val <- as.numeric(msg[2])
              shiny::setProgress(value = val, message = m$message)
            }
          }
        ),
        max = maxval()
      )
      
      
      
    })
    
    
    output$correlation_plot_UI <- renderPlot({
      req(obj())
      plot_pep_correlations(pep.data = obj())
    })
     
    # observeEvent(input$valid_btn, {
    #   
    #   dataOut$value <- my_pipeline_llkimpute(bouyssie,  extension = input$extension)
    #   
    #   print('test')
    #   dataOut$trigger <- as.numeric(Sys.time())
    #   dataOut$widgets <- list(extension = input$extension)
    # })
    
    
    reactive({dataOut})
  })
  
}


