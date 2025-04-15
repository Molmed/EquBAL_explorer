library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
fullobjconf = readRDS("fullobjconf.rds")
fullobjdef  = readRDS("fullobjdef.rds")



tcellsobjconf = readRDS("tcellsobjconf.rds")
tcellsobjdef  = readRDS("tcellsobjdef.rds")



macrophagesobjconf = readRDS("macrophagesobjconf.rds")
macrophagesobjdef  = readRDS("macrophagesobjdef.rds")



dendriticobjconf = readRDS("dendriticobjconf.rds")
dendriticobjdef  = readRDS("dendriticobjdef.rds")



neutrophilsobjconf = readRDS("neutrophilsobjconf.rds")
neutrophilsobjdef  = readRDS("neutrophilsobjdef.rds")



mastcellsobjconf = readRDS("mastcellsobjconf.rds")
mastcellsobjdef  = readRDS("mastcellsobjdef.rds")



### Start server code 
shinyUI(fluidPage( 
### HTML formatting of error messages 
 
tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))), 
list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))), 
 
   
### Page title 
titlePanel("EquBAL explorer"),  
navbarPage( 
  NULL,  
 navbarMenu("Full dataset",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("fullobja1drX", "X-axis:", choices = fullobjconf[dimred == TRUE]$UI, 
                           selected = fullobjdef$dimred[1]), 
            selectInput("fullobja1drY", "Y-axis:", choices = fullobjconf[dimred == TRUE]$UI, 
                        selected = fullobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("fullobja1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.fullobja1togL % 2 == 1", 
          selectInput("fullobja1sub1", "Cell information to subset:", 
                      choices = fullobjconf[grp == TRUE]$UI, 
                      selected = fullobjdef$grp1), 
          uiOutput("fullobja1sub1.ui"), 
          actionButton("fullobja1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("fullobja1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("fullobja1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.fullobja1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("fullobja1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("fullobja1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("fullobja1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("fullobja1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("fullobja1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("fullobja1inp1", "Cell information:", 
                           choices = fullobjconf$UI, 
                           selected = fullobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("fullobja1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.fullobja1tog1 % 2 == 1", 
              radioButtons("fullobja1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("fullobja1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("fullobja1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("fullobja1oup1.ui"))), 
        downloadButton("fullobja1oup1.pdf", "Download PDF"), 
        downloadButton("fullobja1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("fullobja1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("fullobja1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("fullobja1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.fullobja1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("fullobja1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("fullobja1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("fullobja1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("fullobja1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.fullobja1tog2 % 2 == 1", 
              radioButtons("fullobja1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("fullobja1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("fullobja1oup2.ui"))), 
        downloadButton("fullobja1oup2.pdf", "Download PDF"), 
        downloadButton("fullobja1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("fullobja1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("fullobja1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("fullobja2drX", "X-axis:", choices = fullobjconf[dimred == TRUE]$UI, 
                           selected = fullobjdef$dimred[1]), 
            selectInput("fullobja2drY", "Y-axis:", choices = fullobjconf[dimred == TRUE]$UI, 
                        selected = fullobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("fullobja2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.fullobja2togL % 2 == 1", 
          selectInput("fullobja2sub1", "Cell information to subset:", 
                      choices = fullobjconf[grp == TRUE]$UI, 
                      selected = fullobjdef$grp1), 
          uiOutput("fullobja2sub1.ui"), 
          actionButton("fullobja2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("fullobja2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("fullobja2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.fullobja2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("fullobja2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("fullobja2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("fullobja2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("fullobja2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("fullobja2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("fullobja2inp1", "Cell information:", 
                           choices = fullobjconf$UI, 
                           selected = fullobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("fullobja2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.fullobja2tog1 % 2 == 1", 
              radioButtons("fullobja2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("fullobja2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("fullobja2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("fullobja2oup1.ui"))), 
        downloadButton("fullobja2oup1.pdf", "Download PDF"), 
        downloadButton("fullobja2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("fullobja2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("fullobja2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("fullobja2inp2", "Cell information:", 
                           choices = fullobjconf$UI, 
                           selected = fullobjdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("fullobja2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.fullobja2tog2 % 2 == 1", 
              radioButtons("fullobja2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("fullobja2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("fullobja2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("fullobja2oup2.ui"))), 
        downloadButton("fullobja2oup2.pdf", "Download PDF"), 
        downloadButton("fullobja2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("fullobja2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("fullobja2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("fullobja3drX", "X-axis:", choices = fullobjconf[dimred == TRUE]$UI, 
                           selected = fullobjdef$dimred[1]), 
            selectInput("fullobja3drY", "Y-axis:", choices = fullobjconf[dimred == TRUE]$UI, 
                        selected = fullobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("fullobja3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.fullobja3togL % 2 == 1", 
          selectInput("fullobja3sub1", "Cell information to subset:", 
                      choices = fullobjconf[grp == TRUE]$UI, 
                      selected = fullobjdef$grp1), 
          uiOutput("fullobja3sub1.ui"), 
          actionButton("fullobja3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("fullobja3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("fullobja3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.fullobja3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("fullobja3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("fullobja3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("fullobja3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("fullobja3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("fullobja3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("fullobja3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("fullobja3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.fullobja3tog1 % 2 == 1", 
              radioButtons("fullobja3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("fullobja3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("fullobja3oup1.ui"))), 
        downloadButton("fullobja3oup1.pdf", "Download PDF"), 
        downloadButton("fullobja3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("fullobja3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("fullobja3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("fullobja3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("fullobja3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.fullobja3tog2 % 2 == 1", 
              radioButtons("fullobja3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("fullobja3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("fullobja3oup2.ui"))), 
        downloadButton("fullobja3oup2.pdf", "Download PDF"), 
        downloadButton("fullobja3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("fullobja3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("fullobja3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("fullobjb2drX", "X-axis:", choices = fullobjconf[dimred == TRUE]$UI, 
                           selected = fullobjdef$dimred[1]), 
           selectInput("fullobjb2drY", "Y-axis:", choices = fullobjconf[dimred == TRUE]$UI, 
                       selected = fullobjdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("fullobjb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.fullobjb2togL % 2 == 1", 
         selectInput("fullobjb2sub1", "Cell information to subset:", 
                     choices = fullobjconf[grp == TRUE]$UI, 
                    selected = fullobjdef$grp1), 
         uiOutput("fullobjb2sub1.ui"), 
         actionButton("fullobjb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("fullobjb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("fullobjb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.fullobjb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("fullobjb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("fullobjb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("fullobjb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("fullobjb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("fullobjb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("fullobjb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("fullobjb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("fullobjb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.fullobjb2tog1 % 2 == 1", 
         radioButtons("fullobjb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("fullobjb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("fullobjb2oup1.ui"), 
       downloadButton("fullobjb2oup1.pdf", "Download PDF"), 
       downloadButton("fullobjb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("fullobjb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("fullobjb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("fullobjb2oup2.ui"), 
       downloadButton("fullobjb2oup2.pdf", "Download PDF"), 
       downloadButton("fullobjb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("fullobjb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("fullobjc1inp1", "Cell information (X-axis):", 
                   choices = fullobjconf[grp == TRUE]$UI, 
                   selected = fullobjdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("fullobjc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("fullobjc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("fullobjc1pts", "Show data points", value = FALSE), 
       actionButton("fullobjc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.fullobjc1togL % 2 == 1", 
         selectInput("fullobjc1sub1", "Cell information to subset:", 
                     choices = fullobjconf[grp == TRUE]$UI, 
                     selected = fullobjdef$grp1), 
         uiOutput("fullobjc1sub1.ui"), 
         actionButton("fullobjc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("fullobjc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("fullobjc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.fullobjc1tog % 2 == 1", 
         sliderInput("fullobjc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("fullobjc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("fullobjc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("fullobjc1oup.ui"),  
            downloadButton("fullobjc1oup.pdf", "Download PDF"),  
            downloadButton("fullobjc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("fullobjc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("fullobjc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("fullobjc2inp1", "Cell information to plot (X-axis):", 
                  choices = fullobjconf[grp == TRUE]$UI, 
                  selected = fullobjdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("fullobjc2inp2", "Cell information to group / colour by:", 
                  choices = fullobjconf[grp == TRUE]$UI, 
                  selected = fullobjdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("fullobjc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("fullobjc2flp", "Flip X/Y", value = FALSE), 
      actionButton("fullobjc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.fullobjc2togL % 2 == 1", 
        selectInput("fullobjc2sub1", "Cell information to subset:", 
                    choices = fullobjconf[grp == TRUE]$UI, 
                    selected = fullobjdef$grp1), 
        uiOutput("fullobjc2sub1.ui"), 
        actionButton("fullobjc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("fullobjc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("fullobjc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.fullobjc2tog % 2 == 1", 
        radioButtons("fullobjc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("fullobjc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("fullobjc2oup.ui"),  
           downloadButton("fullobjc2oup.pdf", "Download PDF"),  
           downloadButton("fullobjc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("fullobjc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("fullobjc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("fullobjd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(fullobjdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("fullobjd1grp", "Group by:", 
                    choices = fullobjconf[grp == TRUE]$UI, 
                    selected = fullobjconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("fullobjd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("fullobjd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("fullobjd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("fullobjd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("fullobjd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.fullobjd1togL % 2 == 1", 
          selectInput("fullobjd1sub1", "Cell information to subset:", 
                      choices = fullobjconf[grp == TRUE]$UI, 
                      selected = fullobjdef$grp1), 
          uiOutput("fullobjd1sub1.ui"), 
          actionButton("fullobjd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("fullobjd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("fullobjd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.fullobjd1tog % 2 == 1", 
          radioButtons("fullobjd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("fullobjd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("fullobjd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("fullobjd1oupTxt")), 
             uiOutput("fullobjd1oup.ui"), 
             downloadButton("fullobjd1oup.pdf", "Download PDF"), 
             downloadButton("fullobjd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("fullobjd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("fullobjd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("T-cells subset",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("tcellsobja1drX", "X-axis:", choices = tcellsobjconf[dimred == TRUE]$UI, 
                           selected = tcellsobjdef$dimred[1]), 
            selectInput("tcellsobja1drY", "Y-axis:", choices = tcellsobjconf[dimred == TRUE]$UI, 
                        selected = tcellsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("tcellsobja1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.tcellsobja1togL % 2 == 1", 
          selectInput("tcellsobja1sub1", "Cell information to subset:", 
                      choices = tcellsobjconf[grp == TRUE]$UI, 
                      selected = tcellsobjdef$grp1), 
          uiOutput("tcellsobja1sub1.ui"), 
          actionButton("tcellsobja1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("tcellsobja1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("tcellsobja1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.tcellsobja1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("tcellsobja1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("tcellsobja1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("tcellsobja1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("tcellsobja1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("tcellsobja1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("tcellsobja1inp1", "Cell information:", 
                           choices = tcellsobjconf$UI, 
                           selected = tcellsobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("tcellsobja1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.tcellsobja1tog1 % 2 == 1", 
              radioButtons("tcellsobja1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("tcellsobja1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("tcellsobja1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("tcellsobja1oup1.ui"))), 
        downloadButton("tcellsobja1oup1.pdf", "Download PDF"), 
        downloadButton("tcellsobja1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("tcellsobja1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("tcellsobja1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("tcellsobja1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.tcellsobja1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("tcellsobja1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("tcellsobja1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("tcellsobja1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("tcellsobja1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.tcellsobja1tog2 % 2 == 1", 
              radioButtons("tcellsobja1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("tcellsobja1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("tcellsobja1oup2.ui"))), 
        downloadButton("tcellsobja1oup2.pdf", "Download PDF"), 
        downloadButton("tcellsobja1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("tcellsobja1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("tcellsobja1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("tcellsobja2drX", "X-axis:", choices = tcellsobjconf[dimred == TRUE]$UI, 
                           selected = tcellsobjdef$dimred[1]), 
            selectInput("tcellsobja2drY", "Y-axis:", choices = tcellsobjconf[dimred == TRUE]$UI, 
                        selected = tcellsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("tcellsobja2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.tcellsobja2togL % 2 == 1", 
          selectInput("tcellsobja2sub1", "Cell information to subset:", 
                      choices = tcellsobjconf[grp == TRUE]$UI, 
                      selected = tcellsobjdef$grp1), 
          uiOutput("tcellsobja2sub1.ui"), 
          actionButton("tcellsobja2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("tcellsobja2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("tcellsobja2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.tcellsobja2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("tcellsobja2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("tcellsobja2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("tcellsobja2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("tcellsobja2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("tcellsobja2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("tcellsobja2inp1", "Cell information:", 
                           choices = tcellsobjconf$UI, 
                           selected = tcellsobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("tcellsobja2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.tcellsobja2tog1 % 2 == 1", 
              radioButtons("tcellsobja2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("tcellsobja2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("tcellsobja2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("tcellsobja2oup1.ui"))), 
        downloadButton("tcellsobja2oup1.pdf", "Download PDF"), 
        downloadButton("tcellsobja2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("tcellsobja2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("tcellsobja2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("tcellsobja2inp2", "Cell information:", 
                           choices = tcellsobjconf$UI, 
                           selected = tcellsobjdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("tcellsobja2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.tcellsobja2tog2 % 2 == 1", 
              radioButtons("tcellsobja2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("tcellsobja2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("tcellsobja2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("tcellsobja2oup2.ui"))), 
        downloadButton("tcellsobja2oup2.pdf", "Download PDF"), 
        downloadButton("tcellsobja2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("tcellsobja2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("tcellsobja2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("tcellsobja3drX", "X-axis:", choices = tcellsobjconf[dimred == TRUE]$UI, 
                           selected = tcellsobjdef$dimred[1]), 
            selectInput("tcellsobja3drY", "Y-axis:", choices = tcellsobjconf[dimred == TRUE]$UI, 
                        selected = tcellsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("tcellsobja3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.tcellsobja3togL % 2 == 1", 
          selectInput("tcellsobja3sub1", "Cell information to subset:", 
                      choices = tcellsobjconf[grp == TRUE]$UI, 
                      selected = tcellsobjdef$grp1), 
          uiOutput("tcellsobja3sub1.ui"), 
          actionButton("tcellsobja3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("tcellsobja3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("tcellsobja3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.tcellsobja3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("tcellsobja3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("tcellsobja3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("tcellsobja3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("tcellsobja3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("tcellsobja3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("tcellsobja3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("tcellsobja3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.tcellsobja3tog1 % 2 == 1", 
              radioButtons("tcellsobja3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("tcellsobja3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("tcellsobja3oup1.ui"))), 
        downloadButton("tcellsobja3oup1.pdf", "Download PDF"), 
        downloadButton("tcellsobja3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("tcellsobja3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("tcellsobja3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("tcellsobja3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("tcellsobja3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.tcellsobja3tog2 % 2 == 1", 
              radioButtons("tcellsobja3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("tcellsobja3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("tcellsobja3oup2.ui"))), 
        downloadButton("tcellsobja3oup2.pdf", "Download PDF"), 
        downloadButton("tcellsobja3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("tcellsobja3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("tcellsobja3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("tcellsobjb2drX", "X-axis:", choices = tcellsobjconf[dimred == TRUE]$UI, 
                           selected = tcellsobjdef$dimred[1]), 
           selectInput("tcellsobjb2drY", "Y-axis:", choices = tcellsobjconf[dimred == TRUE]$UI, 
                       selected = tcellsobjdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("tcellsobjb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.tcellsobjb2togL % 2 == 1", 
         selectInput("tcellsobjb2sub1", "Cell information to subset:", 
                     choices = tcellsobjconf[grp == TRUE]$UI, 
                    selected = tcellsobjdef$grp1), 
         uiOutput("tcellsobjb2sub1.ui"), 
         actionButton("tcellsobjb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("tcellsobjb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("tcellsobjb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.tcellsobjb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("tcellsobjb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("tcellsobjb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("tcellsobjb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("tcellsobjb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("tcellsobjb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("tcellsobjb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("tcellsobjb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("tcellsobjb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.tcellsobjb2tog1 % 2 == 1", 
         radioButtons("tcellsobjb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("tcellsobjb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("tcellsobjb2oup1.ui"), 
       downloadButton("tcellsobjb2oup1.pdf", "Download PDF"), 
       downloadButton("tcellsobjb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("tcellsobjb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("tcellsobjb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("tcellsobjb2oup2.ui"), 
       downloadButton("tcellsobjb2oup2.pdf", "Download PDF"), 
       downloadButton("tcellsobjb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("tcellsobjb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("tcellsobjc1inp1", "Cell information (X-axis):", 
                   choices = tcellsobjconf[grp == TRUE]$UI, 
                   selected = tcellsobjdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("tcellsobjc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("tcellsobjc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("tcellsobjc1pts", "Show data points", value = FALSE), 
       actionButton("tcellsobjc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.tcellsobjc1togL % 2 == 1", 
         selectInput("tcellsobjc1sub1", "Cell information to subset:", 
                     choices = tcellsobjconf[grp == TRUE]$UI, 
                     selected = tcellsobjdef$grp1), 
         uiOutput("tcellsobjc1sub1.ui"), 
         actionButton("tcellsobjc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("tcellsobjc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("tcellsobjc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.tcellsobjc1tog % 2 == 1", 
         sliderInput("tcellsobjc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("tcellsobjc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("tcellsobjc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("tcellsobjc1oup.ui"),  
            downloadButton("tcellsobjc1oup.pdf", "Download PDF"),  
            downloadButton("tcellsobjc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("tcellsobjc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("tcellsobjc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("tcellsobjc2inp1", "Cell information to plot (X-axis):", 
                  choices = tcellsobjconf[grp == TRUE]$UI, 
                  selected = tcellsobjdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("tcellsobjc2inp2", "Cell information to group / colour by:", 
                  choices = tcellsobjconf[grp == TRUE]$UI, 
                  selected = tcellsobjdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("tcellsobjc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("tcellsobjc2flp", "Flip X/Y", value = FALSE), 
      actionButton("tcellsobjc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.tcellsobjc2togL % 2 == 1", 
        selectInput("tcellsobjc2sub1", "Cell information to subset:", 
                    choices = tcellsobjconf[grp == TRUE]$UI, 
                    selected = tcellsobjdef$grp1), 
        uiOutput("tcellsobjc2sub1.ui"), 
        actionButton("tcellsobjc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("tcellsobjc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("tcellsobjc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.tcellsobjc2tog % 2 == 1", 
        radioButtons("tcellsobjc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("tcellsobjc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("tcellsobjc2oup.ui"),  
           downloadButton("tcellsobjc2oup.pdf", "Download PDF"),  
           downloadButton("tcellsobjc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("tcellsobjc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("tcellsobjc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("tcellsobjd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(tcellsobjdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("tcellsobjd1grp", "Group by:", 
                    choices = tcellsobjconf[grp == TRUE]$UI, 
                    selected = tcellsobjconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("tcellsobjd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("tcellsobjd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("tcellsobjd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("tcellsobjd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("tcellsobjd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.tcellsobjd1togL % 2 == 1", 
          selectInput("tcellsobjd1sub1", "Cell information to subset:", 
                      choices = tcellsobjconf[grp == TRUE]$UI, 
                      selected = tcellsobjdef$grp1), 
          uiOutput("tcellsobjd1sub1.ui"), 
          actionButton("tcellsobjd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("tcellsobjd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("tcellsobjd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.tcellsobjd1tog % 2 == 1", 
          radioButtons("tcellsobjd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("tcellsobjd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("tcellsobjd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("tcellsobjd1oupTxt")), 
             uiOutput("tcellsobjd1oup.ui"), 
             downloadButton("tcellsobjd1oup.pdf", "Download PDF"), 
             downloadButton("tcellsobjd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("tcellsobjd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("tcellsobjd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("Macrophages subset",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("macrophagesobja1drX", "X-axis:", choices = macrophagesobjconf[dimred == TRUE]$UI, 
                           selected = macrophagesobjdef$dimred[1]), 
            selectInput("macrophagesobja1drY", "Y-axis:", choices = macrophagesobjconf[dimred == TRUE]$UI, 
                        selected = macrophagesobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("macrophagesobja1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.macrophagesobja1togL % 2 == 1", 
          selectInput("macrophagesobja1sub1", "Cell information to subset:", 
                      choices = macrophagesobjconf[grp == TRUE]$UI, 
                      selected = macrophagesobjdef$grp1), 
          uiOutput("macrophagesobja1sub1.ui"), 
          actionButton("macrophagesobja1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("macrophagesobja1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("macrophagesobja1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.macrophagesobja1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("macrophagesobja1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("macrophagesobja1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("macrophagesobja1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("macrophagesobja1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("macrophagesobja1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("macrophagesobja1inp1", "Cell information:", 
                           choices = macrophagesobjconf$UI, 
                           selected = macrophagesobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("macrophagesobja1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.macrophagesobja1tog1 % 2 == 1", 
              radioButtons("macrophagesobja1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("macrophagesobja1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("macrophagesobja1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("macrophagesobja1oup1.ui"))), 
        downloadButton("macrophagesobja1oup1.pdf", "Download PDF"), 
        downloadButton("macrophagesobja1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("macrophagesobja1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.macrophagesobja1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("macrophagesobja1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("macrophagesobja1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("macrophagesobja1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("macrophagesobja1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.macrophagesobja1tog2 % 2 == 1", 
              radioButtons("macrophagesobja1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("macrophagesobja1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("macrophagesobja1oup2.ui"))), 
        downloadButton("macrophagesobja1oup2.pdf", "Download PDF"), 
        downloadButton("macrophagesobja1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("macrophagesobja2drX", "X-axis:", choices = macrophagesobjconf[dimred == TRUE]$UI, 
                           selected = macrophagesobjdef$dimred[1]), 
            selectInput("macrophagesobja2drY", "Y-axis:", choices = macrophagesobjconf[dimred == TRUE]$UI, 
                        selected = macrophagesobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("macrophagesobja2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.macrophagesobja2togL % 2 == 1", 
          selectInput("macrophagesobja2sub1", "Cell information to subset:", 
                      choices = macrophagesobjconf[grp == TRUE]$UI, 
                      selected = macrophagesobjdef$grp1), 
          uiOutput("macrophagesobja2sub1.ui"), 
          actionButton("macrophagesobja2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("macrophagesobja2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("macrophagesobja2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.macrophagesobja2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("macrophagesobja2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("macrophagesobja2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("macrophagesobja2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("macrophagesobja2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("macrophagesobja2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("macrophagesobja2inp1", "Cell information:", 
                           choices = macrophagesobjconf$UI, 
                           selected = macrophagesobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("macrophagesobja2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.macrophagesobja2tog1 % 2 == 1", 
              radioButtons("macrophagesobja2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("macrophagesobja2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("macrophagesobja2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("macrophagesobja2oup1.ui"))), 
        downloadButton("macrophagesobja2oup1.pdf", "Download PDF"), 
        downloadButton("macrophagesobja2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("macrophagesobja2inp2", "Cell information:", 
                           choices = macrophagesobjconf$UI, 
                           selected = macrophagesobjdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("macrophagesobja2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.macrophagesobja2tog2 % 2 == 1", 
              radioButtons("macrophagesobja2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("macrophagesobja2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("macrophagesobja2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("macrophagesobja2oup2.ui"))), 
        downloadButton("macrophagesobja2oup2.pdf", "Download PDF"), 
        downloadButton("macrophagesobja2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("macrophagesobja3drX", "X-axis:", choices = macrophagesobjconf[dimred == TRUE]$UI, 
                           selected = macrophagesobjdef$dimred[1]), 
            selectInput("macrophagesobja3drY", "Y-axis:", choices = macrophagesobjconf[dimred == TRUE]$UI, 
                        selected = macrophagesobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("macrophagesobja3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.macrophagesobja3togL % 2 == 1", 
          selectInput("macrophagesobja3sub1", "Cell information to subset:", 
                      choices = macrophagesobjconf[grp == TRUE]$UI, 
                      selected = macrophagesobjdef$grp1), 
          uiOutput("macrophagesobja3sub1.ui"), 
          actionButton("macrophagesobja3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("macrophagesobja3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("macrophagesobja3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.macrophagesobja3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("macrophagesobja3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("macrophagesobja3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("macrophagesobja3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("macrophagesobja3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("macrophagesobja3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("macrophagesobja3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("macrophagesobja3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.macrophagesobja3tog1 % 2 == 1", 
              radioButtons("macrophagesobja3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("macrophagesobja3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("macrophagesobja3oup1.ui"))), 
        downloadButton("macrophagesobja3oup1.pdf", "Download PDF"), 
        downloadButton("macrophagesobja3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("macrophagesobja3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("macrophagesobja3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.macrophagesobja3tog2 % 2 == 1", 
              radioButtons("macrophagesobja3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("macrophagesobja3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("macrophagesobja3oup2.ui"))), 
        downloadButton("macrophagesobja3oup2.pdf", "Download PDF"), 
        downloadButton("macrophagesobja3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("macrophagesobja3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("macrophagesobjb2drX", "X-axis:", choices = macrophagesobjconf[dimred == TRUE]$UI, 
                           selected = macrophagesobjdef$dimred[1]), 
           selectInput("macrophagesobjb2drY", "Y-axis:", choices = macrophagesobjconf[dimred == TRUE]$UI, 
                       selected = macrophagesobjdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("macrophagesobjb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.macrophagesobjb2togL % 2 == 1", 
         selectInput("macrophagesobjb2sub1", "Cell information to subset:", 
                     choices = macrophagesobjconf[grp == TRUE]$UI, 
                    selected = macrophagesobjdef$grp1), 
         uiOutput("macrophagesobjb2sub1.ui"), 
         actionButton("macrophagesobjb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("macrophagesobjb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("macrophagesobjb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.macrophagesobjb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("macrophagesobjb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("macrophagesobjb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("macrophagesobjb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("macrophagesobjb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("macrophagesobjb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("macrophagesobjb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("macrophagesobjb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("macrophagesobjb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.macrophagesobjb2tog1 % 2 == 1", 
         radioButtons("macrophagesobjb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("macrophagesobjb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("macrophagesobjb2oup1.ui"), 
       downloadButton("macrophagesobjb2oup1.pdf", "Download PDF"), 
       downloadButton("macrophagesobjb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("macrophagesobjb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("macrophagesobjb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("macrophagesobjb2oup2.ui"), 
       downloadButton("macrophagesobjb2oup2.pdf", "Download PDF"), 
       downloadButton("macrophagesobjb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("macrophagesobjb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("macrophagesobjc1inp1", "Cell information (X-axis):", 
                   choices = macrophagesobjconf[grp == TRUE]$UI, 
                   selected = macrophagesobjdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("macrophagesobjc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("macrophagesobjc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("macrophagesobjc1pts", "Show data points", value = FALSE), 
       actionButton("macrophagesobjc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.macrophagesobjc1togL % 2 == 1", 
         selectInput("macrophagesobjc1sub1", "Cell information to subset:", 
                     choices = macrophagesobjconf[grp == TRUE]$UI, 
                     selected = macrophagesobjdef$grp1), 
         uiOutput("macrophagesobjc1sub1.ui"), 
         actionButton("macrophagesobjc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("macrophagesobjc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("macrophagesobjc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.macrophagesobjc1tog % 2 == 1", 
         sliderInput("macrophagesobjc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("macrophagesobjc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("macrophagesobjc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("macrophagesobjc1oup.ui"),  
            downloadButton("macrophagesobjc1oup.pdf", "Download PDF"),  
            downloadButton("macrophagesobjc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("macrophagesobjc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("macrophagesobjc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("macrophagesobjc2inp1", "Cell information to plot (X-axis):", 
                  choices = macrophagesobjconf[grp == TRUE]$UI, 
                  selected = macrophagesobjdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("macrophagesobjc2inp2", "Cell information to group / colour by:", 
                  choices = macrophagesobjconf[grp == TRUE]$UI, 
                  selected = macrophagesobjdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("macrophagesobjc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("macrophagesobjc2flp", "Flip X/Y", value = FALSE), 
      actionButton("macrophagesobjc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.macrophagesobjc2togL % 2 == 1", 
        selectInput("macrophagesobjc2sub1", "Cell information to subset:", 
                    choices = macrophagesobjconf[grp == TRUE]$UI, 
                    selected = macrophagesobjdef$grp1), 
        uiOutput("macrophagesobjc2sub1.ui"), 
        actionButton("macrophagesobjc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("macrophagesobjc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("macrophagesobjc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.macrophagesobjc2tog % 2 == 1", 
        radioButtons("macrophagesobjc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("macrophagesobjc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("macrophagesobjc2oup.ui"),  
           downloadButton("macrophagesobjc2oup.pdf", "Download PDF"),  
           downloadButton("macrophagesobjc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("macrophagesobjc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("macrophagesobjc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("macrophagesobjd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(macrophagesobjdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("macrophagesobjd1grp", "Group by:", 
                    choices = macrophagesobjconf[grp == TRUE]$UI, 
                    selected = macrophagesobjconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("macrophagesobjd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("macrophagesobjd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("macrophagesobjd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("macrophagesobjd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("macrophagesobjd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.macrophagesobjd1togL % 2 == 1", 
          selectInput("macrophagesobjd1sub1", "Cell information to subset:", 
                      choices = macrophagesobjconf[grp == TRUE]$UI, 
                      selected = macrophagesobjdef$grp1), 
          uiOutput("macrophagesobjd1sub1.ui"), 
          actionButton("macrophagesobjd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("macrophagesobjd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("macrophagesobjd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.macrophagesobjd1tog % 2 == 1", 
          radioButtons("macrophagesobjd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("macrophagesobjd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("macrophagesobjd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("macrophagesobjd1oupTxt")), 
             uiOutput("macrophagesobjd1oup.ui"), 
             downloadButton("macrophagesobjd1oup.pdf", "Download PDF"), 
             downloadButton("macrophagesobjd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("macrophagesobjd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("macrophagesobjd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("Dendritic cells subset",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("dendriticobja1drX", "X-axis:", choices = dendriticobjconf[dimred == TRUE]$UI, 
                           selected = dendriticobjdef$dimred[1]), 
            selectInput("dendriticobja1drY", "Y-axis:", choices = dendriticobjconf[dimred == TRUE]$UI, 
                        selected = dendriticobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("dendriticobja1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.dendriticobja1togL % 2 == 1", 
          selectInput("dendriticobja1sub1", "Cell information to subset:", 
                      choices = dendriticobjconf[grp == TRUE]$UI, 
                      selected = dendriticobjdef$grp1), 
          uiOutput("dendriticobja1sub1.ui"), 
          actionButton("dendriticobja1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("dendriticobja1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("dendriticobja1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.dendriticobja1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("dendriticobja1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("dendriticobja1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("dendriticobja1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("dendriticobja1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("dendriticobja1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("dendriticobja1inp1", "Cell information:", 
                           choices = dendriticobjconf$UI, 
                           selected = dendriticobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("dendriticobja1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.dendriticobja1tog1 % 2 == 1", 
              radioButtons("dendriticobja1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("dendriticobja1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("dendriticobja1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("dendriticobja1oup1.ui"))), 
        downloadButton("dendriticobja1oup1.pdf", "Download PDF"), 
        downloadButton("dendriticobja1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("dendriticobja1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("dendriticobja1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("dendriticobja1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.dendriticobja1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("dendriticobja1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("dendriticobja1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("dendriticobja1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("dendriticobja1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.dendriticobja1tog2 % 2 == 1", 
              radioButtons("dendriticobja1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("dendriticobja1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("dendriticobja1oup2.ui"))), 
        downloadButton("dendriticobja1oup2.pdf", "Download PDF"), 
        downloadButton("dendriticobja1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("dendriticobja1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("dendriticobja1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("dendriticobja2drX", "X-axis:", choices = dendriticobjconf[dimred == TRUE]$UI, 
                           selected = dendriticobjdef$dimred[1]), 
            selectInput("dendriticobja2drY", "Y-axis:", choices = dendriticobjconf[dimred == TRUE]$UI, 
                        selected = dendriticobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("dendriticobja2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.dendriticobja2togL % 2 == 1", 
          selectInput("dendriticobja2sub1", "Cell information to subset:", 
                      choices = dendriticobjconf[grp == TRUE]$UI, 
                      selected = dendriticobjdef$grp1), 
          uiOutput("dendriticobja2sub1.ui"), 
          actionButton("dendriticobja2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("dendriticobja2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("dendriticobja2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.dendriticobja2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("dendriticobja2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("dendriticobja2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("dendriticobja2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("dendriticobja2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("dendriticobja2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("dendriticobja2inp1", "Cell information:", 
                           choices = dendriticobjconf$UI, 
                           selected = dendriticobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("dendriticobja2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.dendriticobja2tog1 % 2 == 1", 
              radioButtons("dendriticobja2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("dendriticobja2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("dendriticobja2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("dendriticobja2oup1.ui"))), 
        downloadButton("dendriticobja2oup1.pdf", "Download PDF"), 
        downloadButton("dendriticobja2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("dendriticobja2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("dendriticobja2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("dendriticobja2inp2", "Cell information:", 
                           choices = dendriticobjconf$UI, 
                           selected = dendriticobjdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("dendriticobja2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.dendriticobja2tog2 % 2 == 1", 
              radioButtons("dendriticobja2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("dendriticobja2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("dendriticobja2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("dendriticobja2oup2.ui"))), 
        downloadButton("dendriticobja2oup2.pdf", "Download PDF"), 
        downloadButton("dendriticobja2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("dendriticobja2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("dendriticobja2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("dendriticobja3drX", "X-axis:", choices = dendriticobjconf[dimred == TRUE]$UI, 
                           selected = dendriticobjdef$dimred[1]), 
            selectInput("dendriticobja3drY", "Y-axis:", choices = dendriticobjconf[dimred == TRUE]$UI, 
                        selected = dendriticobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("dendriticobja3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.dendriticobja3togL % 2 == 1", 
          selectInput("dendriticobja3sub1", "Cell information to subset:", 
                      choices = dendriticobjconf[grp == TRUE]$UI, 
                      selected = dendriticobjdef$grp1), 
          uiOutput("dendriticobja3sub1.ui"), 
          actionButton("dendriticobja3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("dendriticobja3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("dendriticobja3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.dendriticobja3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("dendriticobja3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("dendriticobja3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("dendriticobja3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("dendriticobja3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("dendriticobja3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("dendriticobja3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("dendriticobja3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.dendriticobja3tog1 % 2 == 1", 
              radioButtons("dendriticobja3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("dendriticobja3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("dendriticobja3oup1.ui"))), 
        downloadButton("dendriticobja3oup1.pdf", "Download PDF"), 
        downloadButton("dendriticobja3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("dendriticobja3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("dendriticobja3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("dendriticobja3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("dendriticobja3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.dendriticobja3tog2 % 2 == 1", 
              radioButtons("dendriticobja3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("dendriticobja3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("dendriticobja3oup2.ui"))), 
        downloadButton("dendriticobja3oup2.pdf", "Download PDF"), 
        downloadButton("dendriticobja3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("dendriticobja3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("dendriticobja3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("dendriticobjb2drX", "X-axis:", choices = dendriticobjconf[dimred == TRUE]$UI, 
                           selected = dendriticobjdef$dimred[1]), 
           selectInput("dendriticobjb2drY", "Y-axis:", choices = dendriticobjconf[dimred == TRUE]$UI, 
                       selected = dendriticobjdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("dendriticobjb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.dendriticobjb2togL % 2 == 1", 
         selectInput("dendriticobjb2sub1", "Cell information to subset:", 
                     choices = dendriticobjconf[grp == TRUE]$UI, 
                    selected = dendriticobjdef$grp1), 
         uiOutput("dendriticobjb2sub1.ui"), 
         actionButton("dendriticobjb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("dendriticobjb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("dendriticobjb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.dendriticobjb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("dendriticobjb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("dendriticobjb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("dendriticobjb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("dendriticobjb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("dendriticobjb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("dendriticobjb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("dendriticobjb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("dendriticobjb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.dendriticobjb2tog1 % 2 == 1", 
         radioButtons("dendriticobjb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("dendriticobjb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("dendriticobjb2oup1.ui"), 
       downloadButton("dendriticobjb2oup1.pdf", "Download PDF"), 
       downloadButton("dendriticobjb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("dendriticobjb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("dendriticobjb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("dendriticobjb2oup2.ui"), 
       downloadButton("dendriticobjb2oup2.pdf", "Download PDF"), 
       downloadButton("dendriticobjb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("dendriticobjb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("dendriticobjc1inp1", "Cell information (X-axis):", 
                   choices = dendriticobjconf[grp == TRUE]$UI, 
                   selected = dendriticobjdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("dendriticobjc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("dendriticobjc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("dendriticobjc1pts", "Show data points", value = FALSE), 
       actionButton("dendriticobjc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.dendriticobjc1togL % 2 == 1", 
         selectInput("dendriticobjc1sub1", "Cell information to subset:", 
                     choices = dendriticobjconf[grp == TRUE]$UI, 
                     selected = dendriticobjdef$grp1), 
         uiOutput("dendriticobjc1sub1.ui"), 
         actionButton("dendriticobjc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("dendriticobjc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("dendriticobjc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.dendriticobjc1tog % 2 == 1", 
         sliderInput("dendriticobjc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("dendriticobjc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("dendriticobjc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("dendriticobjc1oup.ui"),  
            downloadButton("dendriticobjc1oup.pdf", "Download PDF"),  
            downloadButton("dendriticobjc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("dendriticobjc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("dendriticobjc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("dendriticobjc2inp1", "Cell information to plot (X-axis):", 
                  choices = dendriticobjconf[grp == TRUE]$UI, 
                  selected = dendriticobjdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("dendriticobjc2inp2", "Cell information to group / colour by:", 
                  choices = dendriticobjconf[grp == TRUE]$UI, 
                  selected = dendriticobjdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("dendriticobjc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("dendriticobjc2flp", "Flip X/Y", value = FALSE), 
      actionButton("dendriticobjc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.dendriticobjc2togL % 2 == 1", 
        selectInput("dendriticobjc2sub1", "Cell information to subset:", 
                    choices = dendriticobjconf[grp == TRUE]$UI, 
                    selected = dendriticobjdef$grp1), 
        uiOutput("dendriticobjc2sub1.ui"), 
        actionButton("dendriticobjc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("dendriticobjc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("dendriticobjc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.dendriticobjc2tog % 2 == 1", 
        radioButtons("dendriticobjc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("dendriticobjc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("dendriticobjc2oup.ui"),  
           downloadButton("dendriticobjc2oup.pdf", "Download PDF"),  
           downloadButton("dendriticobjc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("dendriticobjc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("dendriticobjc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("dendriticobjd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(dendriticobjdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("dendriticobjd1grp", "Group by:", 
                    choices = dendriticobjconf[grp == TRUE]$UI, 
                    selected = dendriticobjconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("dendriticobjd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("dendriticobjd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("dendriticobjd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("dendriticobjd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("dendriticobjd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.dendriticobjd1togL % 2 == 1", 
          selectInput("dendriticobjd1sub1", "Cell information to subset:", 
                      choices = dendriticobjconf[grp == TRUE]$UI, 
                      selected = dendriticobjdef$grp1), 
          uiOutput("dendriticobjd1sub1.ui"), 
          actionButton("dendriticobjd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("dendriticobjd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("dendriticobjd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.dendriticobjd1tog % 2 == 1", 
          radioButtons("dendriticobjd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("dendriticobjd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("dendriticobjd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("dendriticobjd1oupTxt")), 
             uiOutput("dendriticobjd1oup.ui"), 
             downloadButton("dendriticobjd1oup.pdf", "Download PDF"), 
             downloadButton("dendriticobjd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("dendriticobjd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("dendriticobjd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("Neutrophils subset",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("neutrophilsobja1drX", "X-axis:", choices = neutrophilsobjconf[dimred == TRUE]$UI, 
                           selected = neutrophilsobjdef$dimred[1]), 
            selectInput("neutrophilsobja1drY", "Y-axis:", choices = neutrophilsobjconf[dimred == TRUE]$UI, 
                        selected = neutrophilsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("neutrophilsobja1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.neutrophilsobja1togL % 2 == 1", 
          selectInput("neutrophilsobja1sub1", "Cell information to subset:", 
                      choices = neutrophilsobjconf[grp == TRUE]$UI, 
                      selected = neutrophilsobjdef$grp1), 
          uiOutput("neutrophilsobja1sub1.ui"), 
          actionButton("neutrophilsobja1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("neutrophilsobja1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("neutrophilsobja1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.neutrophilsobja1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("neutrophilsobja1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("neutrophilsobja1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("neutrophilsobja1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("neutrophilsobja1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("neutrophilsobja1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("neutrophilsobja1inp1", "Cell information:", 
                           choices = neutrophilsobjconf$UI, 
                           selected = neutrophilsobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("neutrophilsobja1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.neutrophilsobja1tog1 % 2 == 1", 
              radioButtons("neutrophilsobja1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("neutrophilsobja1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("neutrophilsobja1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("neutrophilsobja1oup1.ui"))), 
        downloadButton("neutrophilsobja1oup1.pdf", "Download PDF"), 
        downloadButton("neutrophilsobja1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("neutrophilsobja1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.neutrophilsobja1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("neutrophilsobja1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("neutrophilsobja1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("neutrophilsobja1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("neutrophilsobja1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.neutrophilsobja1tog2 % 2 == 1", 
              radioButtons("neutrophilsobja1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("neutrophilsobja1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("neutrophilsobja1oup2.ui"))), 
        downloadButton("neutrophilsobja1oup2.pdf", "Download PDF"), 
        downloadButton("neutrophilsobja1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("neutrophilsobja2drX", "X-axis:", choices = neutrophilsobjconf[dimred == TRUE]$UI, 
                           selected = neutrophilsobjdef$dimred[1]), 
            selectInput("neutrophilsobja2drY", "Y-axis:", choices = neutrophilsobjconf[dimred == TRUE]$UI, 
                        selected = neutrophilsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("neutrophilsobja2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.neutrophilsobja2togL % 2 == 1", 
          selectInput("neutrophilsobja2sub1", "Cell information to subset:", 
                      choices = neutrophilsobjconf[grp == TRUE]$UI, 
                      selected = neutrophilsobjdef$grp1), 
          uiOutput("neutrophilsobja2sub1.ui"), 
          actionButton("neutrophilsobja2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("neutrophilsobja2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("neutrophilsobja2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.neutrophilsobja2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("neutrophilsobja2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("neutrophilsobja2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("neutrophilsobja2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("neutrophilsobja2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("neutrophilsobja2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("neutrophilsobja2inp1", "Cell information:", 
                           choices = neutrophilsobjconf$UI, 
                           selected = neutrophilsobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("neutrophilsobja2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.neutrophilsobja2tog1 % 2 == 1", 
              radioButtons("neutrophilsobja2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("neutrophilsobja2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("neutrophilsobja2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("neutrophilsobja2oup1.ui"))), 
        downloadButton("neutrophilsobja2oup1.pdf", "Download PDF"), 
        downloadButton("neutrophilsobja2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("neutrophilsobja2inp2", "Cell information:", 
                           choices = neutrophilsobjconf$UI, 
                           selected = neutrophilsobjdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("neutrophilsobja2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.neutrophilsobja2tog2 % 2 == 1", 
              radioButtons("neutrophilsobja2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("neutrophilsobja2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("neutrophilsobja2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("neutrophilsobja2oup2.ui"))), 
        downloadButton("neutrophilsobja2oup2.pdf", "Download PDF"), 
        downloadButton("neutrophilsobja2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("neutrophilsobja3drX", "X-axis:", choices = neutrophilsobjconf[dimred == TRUE]$UI, 
                           selected = neutrophilsobjdef$dimred[1]), 
            selectInput("neutrophilsobja3drY", "Y-axis:", choices = neutrophilsobjconf[dimred == TRUE]$UI, 
                        selected = neutrophilsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("neutrophilsobja3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.neutrophilsobja3togL % 2 == 1", 
          selectInput("neutrophilsobja3sub1", "Cell information to subset:", 
                      choices = neutrophilsobjconf[grp == TRUE]$UI, 
                      selected = neutrophilsobjdef$grp1), 
          uiOutput("neutrophilsobja3sub1.ui"), 
          actionButton("neutrophilsobja3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("neutrophilsobja3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("neutrophilsobja3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.neutrophilsobja3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("neutrophilsobja3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("neutrophilsobja3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("neutrophilsobja3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("neutrophilsobja3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("neutrophilsobja3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("neutrophilsobja3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("neutrophilsobja3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.neutrophilsobja3tog1 % 2 == 1", 
              radioButtons("neutrophilsobja3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("neutrophilsobja3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("neutrophilsobja3oup1.ui"))), 
        downloadButton("neutrophilsobja3oup1.pdf", "Download PDF"), 
        downloadButton("neutrophilsobja3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("neutrophilsobja3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("neutrophilsobja3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.neutrophilsobja3tog2 % 2 == 1", 
              radioButtons("neutrophilsobja3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("neutrophilsobja3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("neutrophilsobja3oup2.ui"))), 
        downloadButton("neutrophilsobja3oup2.pdf", "Download PDF"), 
        downloadButton("neutrophilsobja3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("neutrophilsobja3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("neutrophilsobjb2drX", "X-axis:", choices = neutrophilsobjconf[dimred == TRUE]$UI, 
                           selected = neutrophilsobjdef$dimred[1]), 
           selectInput("neutrophilsobjb2drY", "Y-axis:", choices = neutrophilsobjconf[dimred == TRUE]$UI, 
                       selected = neutrophilsobjdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("neutrophilsobjb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.neutrophilsobjb2togL % 2 == 1", 
         selectInput("neutrophilsobjb2sub1", "Cell information to subset:", 
                     choices = neutrophilsobjconf[grp == TRUE]$UI, 
                    selected = neutrophilsobjdef$grp1), 
         uiOutput("neutrophilsobjb2sub1.ui"), 
         actionButton("neutrophilsobjb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("neutrophilsobjb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("neutrophilsobjb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.neutrophilsobjb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("neutrophilsobjb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("neutrophilsobjb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("neutrophilsobjb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("neutrophilsobjb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("neutrophilsobjb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("neutrophilsobjb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("neutrophilsobjb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("neutrophilsobjb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.neutrophilsobjb2tog1 % 2 == 1", 
         radioButtons("neutrophilsobjb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("neutrophilsobjb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("neutrophilsobjb2oup1.ui"), 
       downloadButton("neutrophilsobjb2oup1.pdf", "Download PDF"), 
       downloadButton("neutrophilsobjb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("neutrophilsobjb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("neutrophilsobjb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("neutrophilsobjb2oup2.ui"), 
       downloadButton("neutrophilsobjb2oup2.pdf", "Download PDF"), 
       downloadButton("neutrophilsobjb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("neutrophilsobjb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("neutrophilsobjc1inp1", "Cell information (X-axis):", 
                   choices = neutrophilsobjconf[grp == TRUE]$UI, 
                   selected = neutrophilsobjdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("neutrophilsobjc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("neutrophilsobjc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("neutrophilsobjc1pts", "Show data points", value = FALSE), 
       actionButton("neutrophilsobjc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.neutrophilsobjc1togL % 2 == 1", 
         selectInput("neutrophilsobjc1sub1", "Cell information to subset:", 
                     choices = neutrophilsobjconf[grp == TRUE]$UI, 
                     selected = neutrophilsobjdef$grp1), 
         uiOutput("neutrophilsobjc1sub1.ui"), 
         actionButton("neutrophilsobjc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("neutrophilsobjc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("neutrophilsobjc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.neutrophilsobjc1tog % 2 == 1", 
         sliderInput("neutrophilsobjc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("neutrophilsobjc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("neutrophilsobjc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("neutrophilsobjc1oup.ui"),  
            downloadButton("neutrophilsobjc1oup.pdf", "Download PDF"),  
            downloadButton("neutrophilsobjc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("neutrophilsobjc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("neutrophilsobjc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("neutrophilsobjc2inp1", "Cell information to plot (X-axis):", 
                  choices = neutrophilsobjconf[grp == TRUE]$UI, 
                  selected = neutrophilsobjdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("neutrophilsobjc2inp2", "Cell information to group / colour by:", 
                  choices = neutrophilsobjconf[grp == TRUE]$UI, 
                  selected = neutrophilsobjdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("neutrophilsobjc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("neutrophilsobjc2flp", "Flip X/Y", value = FALSE), 
      actionButton("neutrophilsobjc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.neutrophilsobjc2togL % 2 == 1", 
        selectInput("neutrophilsobjc2sub1", "Cell information to subset:", 
                    choices = neutrophilsobjconf[grp == TRUE]$UI, 
                    selected = neutrophilsobjdef$grp1), 
        uiOutput("neutrophilsobjc2sub1.ui"), 
        actionButton("neutrophilsobjc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("neutrophilsobjc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("neutrophilsobjc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.neutrophilsobjc2tog % 2 == 1", 
        radioButtons("neutrophilsobjc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("neutrophilsobjc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("neutrophilsobjc2oup.ui"),  
           downloadButton("neutrophilsobjc2oup.pdf", "Download PDF"),  
           downloadButton("neutrophilsobjc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("neutrophilsobjc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("neutrophilsobjc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("neutrophilsobjd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(neutrophilsobjdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("neutrophilsobjd1grp", "Group by:", 
                    choices = neutrophilsobjconf[grp == TRUE]$UI, 
                    selected = neutrophilsobjconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("neutrophilsobjd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("neutrophilsobjd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("neutrophilsobjd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("neutrophilsobjd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("neutrophilsobjd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.neutrophilsobjd1togL % 2 == 1", 
          selectInput("neutrophilsobjd1sub1", "Cell information to subset:", 
                      choices = neutrophilsobjconf[grp == TRUE]$UI, 
                      selected = neutrophilsobjdef$grp1), 
          uiOutput("neutrophilsobjd1sub1.ui"), 
          actionButton("neutrophilsobjd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("neutrophilsobjd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("neutrophilsobjd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.neutrophilsobjd1tog % 2 == 1", 
          radioButtons("neutrophilsobjd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("neutrophilsobjd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("neutrophilsobjd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("neutrophilsobjd1oupTxt")), 
             uiOutput("neutrophilsobjd1oup.ui"), 
             downloadButton("neutrophilsobjd1oup.pdf", "Download PDF"), 
             downloadButton("neutrophilsobjd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("neutrophilsobjd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("neutrophilsobjd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("Mast cells subset",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("mastcellsobja1drX", "X-axis:", choices = mastcellsobjconf[dimred == TRUE]$UI, 
                           selected = mastcellsobjdef$dimred[1]), 
            selectInput("mastcellsobja1drY", "Y-axis:", choices = mastcellsobjconf[dimred == TRUE]$UI, 
                        selected = mastcellsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("mastcellsobja1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.mastcellsobja1togL % 2 == 1", 
          selectInput("mastcellsobja1sub1", "Cell information to subset:", 
                      choices = mastcellsobjconf[grp == TRUE]$UI, 
                      selected = mastcellsobjdef$grp1), 
          uiOutput("mastcellsobja1sub1.ui"), 
          actionButton("mastcellsobja1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("mastcellsobja1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("mastcellsobja1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.mastcellsobja1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("mastcellsobja1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("mastcellsobja1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("mastcellsobja1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("mastcellsobja1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("mastcellsobja1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("mastcellsobja1inp1", "Cell information:", 
                           choices = mastcellsobjconf$UI, 
                           selected = mastcellsobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("mastcellsobja1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.mastcellsobja1tog1 % 2 == 1", 
              radioButtons("mastcellsobja1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("mastcellsobja1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("mastcellsobja1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("mastcellsobja1oup1.ui"))), 
        downloadButton("mastcellsobja1oup1.pdf", "Download PDF"), 
        downloadButton("mastcellsobja1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("mastcellsobja1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.mastcellsobja1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("mastcellsobja1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("mastcellsobja1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("mastcellsobja1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("mastcellsobja1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.mastcellsobja1tog2 % 2 == 1", 
              radioButtons("mastcellsobja1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("mastcellsobja1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("mastcellsobja1oup2.ui"))), 
        downloadButton("mastcellsobja1oup2.pdf", "Download PDF"), 
        downloadButton("mastcellsobja1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("mastcellsobja2drX", "X-axis:", choices = mastcellsobjconf[dimred == TRUE]$UI, 
                           selected = mastcellsobjdef$dimred[1]), 
            selectInput("mastcellsobja2drY", "Y-axis:", choices = mastcellsobjconf[dimred == TRUE]$UI, 
                        selected = mastcellsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("mastcellsobja2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.mastcellsobja2togL % 2 == 1", 
          selectInput("mastcellsobja2sub1", "Cell information to subset:", 
                      choices = mastcellsobjconf[grp == TRUE]$UI, 
                      selected = mastcellsobjdef$grp1), 
          uiOutput("mastcellsobja2sub1.ui"), 
          actionButton("mastcellsobja2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("mastcellsobja2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("mastcellsobja2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.mastcellsobja2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("mastcellsobja2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("mastcellsobja2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("mastcellsobja2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("mastcellsobja2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("mastcellsobja2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("mastcellsobja2inp1", "Cell information:", 
                           choices = mastcellsobjconf$UI, 
                           selected = mastcellsobjdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("mastcellsobja2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.mastcellsobja2tog1 % 2 == 1", 
              radioButtons("mastcellsobja2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("mastcellsobja2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("mastcellsobja2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("mastcellsobja2oup1.ui"))), 
        downloadButton("mastcellsobja2oup1.pdf", "Download PDF"), 
        downloadButton("mastcellsobja2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("mastcellsobja2inp2", "Cell information:", 
                           choices = mastcellsobjconf$UI, 
                           selected = mastcellsobjdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("mastcellsobja2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.mastcellsobja2tog2 % 2 == 1", 
              radioButtons("mastcellsobja2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("mastcellsobja2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("mastcellsobja2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("mastcellsobja2oup2.ui"))), 
        downloadButton("mastcellsobja2oup2.pdf", "Download PDF"), 
        downloadButton("mastcellsobja2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("mastcellsobja3drX", "X-axis:", choices = mastcellsobjconf[dimred == TRUE]$UI, 
                           selected = mastcellsobjdef$dimred[1]), 
            selectInput("mastcellsobja3drY", "Y-axis:", choices = mastcellsobjconf[dimred == TRUE]$UI, 
                        selected = mastcellsobjdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("mastcellsobja3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.mastcellsobja3togL % 2 == 1", 
          selectInput("mastcellsobja3sub1", "Cell information to subset:", 
                      choices = mastcellsobjconf[grp == TRUE]$UI, 
                      selected = mastcellsobjdef$grp1), 
          uiOutput("mastcellsobja3sub1.ui"), 
          actionButton("mastcellsobja3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("mastcellsobja3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("mastcellsobja3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.mastcellsobja3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("mastcellsobja3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("mastcellsobja3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("mastcellsobja3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("mastcellsobja3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("mastcellsobja3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("mastcellsobja3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("mastcellsobja3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.mastcellsobja3tog1 % 2 == 1", 
              radioButtons("mastcellsobja3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("mastcellsobja3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("mastcellsobja3oup1.ui"))), 
        downloadButton("mastcellsobja3oup1.pdf", "Download PDF"), 
        downloadButton("mastcellsobja3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("mastcellsobja3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("mastcellsobja3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.mastcellsobja3tog2 % 2 == 1", 
              radioButtons("mastcellsobja3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("mastcellsobja3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("mastcellsobja3oup2.ui"))), 
        downloadButton("mastcellsobja3oup2.pdf", "Download PDF"), 
        downloadButton("mastcellsobja3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("mastcellsobja3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("mastcellsobjb2drX", "X-axis:", choices = mastcellsobjconf[dimred == TRUE]$UI, 
                           selected = mastcellsobjdef$dimred[1]), 
           selectInput("mastcellsobjb2drY", "Y-axis:", choices = mastcellsobjconf[dimred == TRUE]$UI, 
                       selected = mastcellsobjdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("mastcellsobjb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.mastcellsobjb2togL % 2 == 1", 
         selectInput("mastcellsobjb2sub1", "Cell information to subset:", 
                     choices = mastcellsobjconf[grp == TRUE]$UI, 
                    selected = mastcellsobjdef$grp1), 
         uiOutput("mastcellsobjb2sub1.ui"), 
         actionButton("mastcellsobjb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("mastcellsobjb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("mastcellsobjb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.mastcellsobjb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("mastcellsobjb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("mastcellsobjb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("mastcellsobjb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("mastcellsobjb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("mastcellsobjb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("mastcellsobjb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("mastcellsobjb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("mastcellsobjb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.mastcellsobjb2tog1 % 2 == 1", 
         radioButtons("mastcellsobjb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("mastcellsobjb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("mastcellsobjb2oup1.ui"), 
       downloadButton("mastcellsobjb2oup1.pdf", "Download PDF"), 
       downloadButton("mastcellsobjb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("mastcellsobjb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("mastcellsobjb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("mastcellsobjb2oup2.ui"), 
       downloadButton("mastcellsobjb2oup2.pdf", "Download PDF"), 
       downloadButton("mastcellsobjb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("mastcellsobjb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("mastcellsobjc1inp1", "Cell information (X-axis):", 
                   choices = mastcellsobjconf[grp == TRUE]$UI, 
                   selected = mastcellsobjdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("mastcellsobjc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("mastcellsobjc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("mastcellsobjc1pts", "Show data points", value = FALSE), 
       actionButton("mastcellsobjc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.mastcellsobjc1togL % 2 == 1", 
         selectInput("mastcellsobjc1sub1", "Cell information to subset:", 
                     choices = mastcellsobjconf[grp == TRUE]$UI, 
                     selected = mastcellsobjdef$grp1), 
         uiOutput("mastcellsobjc1sub1.ui"), 
         actionButton("mastcellsobjc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("mastcellsobjc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("mastcellsobjc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.mastcellsobjc1tog % 2 == 1", 
         sliderInput("mastcellsobjc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("mastcellsobjc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("mastcellsobjc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("mastcellsobjc1oup.ui"),  
            downloadButton("mastcellsobjc1oup.pdf", "Download PDF"),  
            downloadButton("mastcellsobjc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("mastcellsobjc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("mastcellsobjc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("mastcellsobjc2inp1", "Cell information to plot (X-axis):", 
                  choices = mastcellsobjconf[grp == TRUE]$UI, 
                  selected = mastcellsobjdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("mastcellsobjc2inp2", "Cell information to group / colour by:", 
                  choices = mastcellsobjconf[grp == TRUE]$UI, 
                  selected = mastcellsobjdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("mastcellsobjc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("mastcellsobjc2flp", "Flip X/Y", value = FALSE), 
      actionButton("mastcellsobjc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.mastcellsobjc2togL % 2 == 1", 
        selectInput("mastcellsobjc2sub1", "Cell information to subset:", 
                    choices = mastcellsobjconf[grp == TRUE]$UI, 
                    selected = mastcellsobjdef$grp1), 
        uiOutput("mastcellsobjc2sub1.ui"), 
        actionButton("mastcellsobjc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("mastcellsobjc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("mastcellsobjc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.mastcellsobjc2tog % 2 == 1", 
        radioButtons("mastcellsobjc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("mastcellsobjc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("mastcellsobjc2oup.ui"),  
           downloadButton("mastcellsobjc2oup.pdf", "Download PDF"),  
           downloadButton("mastcellsobjc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("mastcellsobjc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("mastcellsobjc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("mastcellsobjd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(mastcellsobjdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("mastcellsobjd1grp", "Group by:", 
                    choices = mastcellsobjconf[grp == TRUE]$UI, 
                    selected = mastcellsobjconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("mastcellsobjd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("mastcellsobjd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("mastcellsobjd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("mastcellsobjd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("mastcellsobjd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.mastcellsobjd1togL % 2 == 1", 
          selectInput("mastcellsobjd1sub1", "Cell information to subset:", 
                      choices = mastcellsobjconf[grp == TRUE]$UI, 
                      selected = mastcellsobjdef$grp1), 
          uiOutput("mastcellsobjd1sub1.ui"), 
          actionButton("mastcellsobjd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("mastcellsobjd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("mastcellsobjd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.mastcellsobjd1tog % 2 == 1", 
          radioButtons("mastcellsobjd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("mastcellsobjd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("mastcellsobjd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("mastcellsobjd1oupTxt")), 
             uiOutput("mastcellsobjd1oup.ui"), 
             downloadButton("mastcellsobjd1oup.pdf", "Download PDF"), 
             downloadButton("mastcellsobjd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("mastcellsobjd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("mastcellsobjd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

   
br(), 
p(strong("Reference: "),strong("An Interactive Web Resource for Exploring Equine BAL Cell scRNA-seq Data, "),"(2025) Sulyaeva J., Fegraeus K., Riihimäki M., Nordlund J., Raine A. ","(Funding: FORMAS) ",style = "font-size: 125%;"), 
p(em("This webpage was made using "), a("ShinyCell", 
  href = "https://github.com/SGDDNB/ShinyCell",target="_blank")), 
br(),br(),br(),br(),br() 
))) 
 
 
 
 