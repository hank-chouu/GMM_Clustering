#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)

header <- dashboardHeader(title = "GMM clustering")


sidebar <- dashboardSidebar(
  hr(), 
  sidebarMenu(id='sidebar',
              menuItem('Plot', tabName = 'plot', icon = icon("line-chart")),
              menuItem('Table import', tabName = 'table', icon = icon("table")),
              menuItem('About', tabName = 'about', selected = T, icon = icon("question")))
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName='table', 
      fluidRow(
        box(width = 12, status = "primary", solidHeader = TRUE, title= "Table", 
          fileInput('file1', "User data (.csv format)", accept = ".csv"),
          checkboxInput("use_example", "Or use example Data"),
          uiOutput("select_example"),
          tableOutput('table1'), 
          verbatimTextOutput('message1')
        )
      )
    ),
    tabItem(tabName = 'plot',
      fluidRow(
        column(width=12, 
          tabBox(width=NULL, 
            tabPanel('Clustering',
            fluidRow(column(width=4, 
                            numericInput('num_k', label = 'k components', value = NA)),
                     column(width=8, 
                            verbatimTextOutput('message2'), 
                            plotOutput('plot2'),
                            verbatimTextOutput('para1')))
            ),
            tabPanel('BIC',
            fluidRow(column(width=4, 
                            numericInput('bic_ub', label = 'Upper bound', value=6),
                            numericInput('bic_lb', label = 'Lower bound', value=2)
                            ),
                     column(width=8, 
                            verbatimTextOutput('message3'), 
                            plotOutput('plot3'),
                            verbatimTextOutput('message5')))
            ),
            tabPanel('AIC',
                     fluidRow(column(width=4, 
                                     numericInput('aic_ub', label = 'Upper bound', value=6),
                                     numericInput('aic_lb', label = 'Lower bound', value=2)
                     ),
                     column(width=8, 
                            verbatimTextOutput('message4'), 
                            plotOutput('plot4'),
                            verbatimTextOutput('message6')))
            )
                    
          
        )
      )
    )
  ), 
  tabItem(tabName = 'about',
          fluidPage(
            mainPanel(htmlOutput("showfile")
          )
  ))
))

shinyUI(
  dashboardPage(header, sidebar, body)
)

