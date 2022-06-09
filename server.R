#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tools)


# server

shinyServer(function(input, output) {
  
  ###########
  
  # input csv or use example data
  
  ###########
  
  values <- reactiveValues(df_data = NULL)
  display_df <- NULL
  
  observeEvent(input$file1, {
    
    values$df_data <- NULL
    
    output$message1 <- renderPrint({
      ext <- file_ext(input$file1$datapath)
      validate(need(ext == "csv", "Incorrect file type"))
      values$df_data  <- read.csv(input$file1$datapath)
      
      dimen <- dim(values$df_data)
      
      str1 <- cat(dimen[1], 'rows,', dimen[2], ' dimensions(variables).')
      
      str2 <- NULL
      if(!all(sapply(values$df_data, is.numeric))){
        # values$df_data <- NULL
        str2 <- 'Data contains non-numeric inputs. Re-upload to proceed.'
      }
      
      invisible(cat(str1, '\n', str2, sep=''))
    })
    
  })
  
  observeEvent(input$use_example, {
    
    
    if (input$use_example){
      
      observeEvent(input$selected_example, {
        if(input$selected_example==1){
          values$df_data <- precip
          
        }
        else if (input$selected_example==2){
          values$df_data <- two_dim_example
          
        }else{
          values$df_data <- three_dim_example
          
        }
        
        output$message1 <- renderPrint({
          
          dimen <- dim(values$df_data)
          
          str1 <- cat(dimen[1], 'rows,', dimen[2], 'dimensions(variables).')
          invisible(str1)
        })
        
      })
    }else{
      values$df_data <- NULL
      output$message1 <- renderPrint({

        invisible(NULL)
        
      })
    }
    
  })
  
  output$table1 <- renderTable({
    head(values$df_data)
  })
  
  #########
  
  # data plot
  
  #########
  
  gmm_result <- reactive({

    wrapped_gmm(values$df_data, input$num_k)
  })
  
  
  upper.panel <- function(x, y){
    
    points(x, y, pch=19)
    
    d <- dim(values$df_data)[2]
    n <- dim(values$df_data)[1]
    
    x.point <- seq(min(x), max(x), length=100)
    y.point <- seq(min(y), max(y), length=100)
    
    idx1 <- idx2 <- 0
    for (i in 1:d){
      if(sum(x == values$df_data[, i]) == n){
        idx1 <- i
      }
      if(sum(y == values$df_data[, i]) == n){
        idx2 <- i
      }
    }
    for (k in 1:input$num_k){
      small_mu <- gmm_result()$para$mu[c(idx1, idx2), k]
      small_sigma <- matrix(0, 2, 2)
      small_sigma[1, 1] <- gmm_result()$para$sigma[[k]][idx1, idx1]
      small_sigma[1, 2] <- small_sigma[2, 1] <- gmm_result()$para$sigma[[k]][idx1, idx2]
      small_sigma[2, 2] <- gmm_result()$para$sigma[[k]][idx2, idx2]

      dens <- dmvnorm(as.matrix(expand.grid(x.point, y.point)), small_mu, small_sigma)
      dens <- matrix(dens, 100, 100, byrow = T)
      contour(x.point,y.point, dens, col=k+1, add=T, drawlabels = F, nlevels=5)
    }
  }
  
  upper.panel2 <- function(x, y){
    points(x, y, pch=19)
  }
  
  output$plot2 <- renderPlot({
    
    if (!is.null(values$df_data) & all(sapply(values$df_data, is.numeric))){
      
      if(is.na(input$num_k) | input$num_k%%1 != 0 | input$num_k<=0){
        
        if(dim(values$df_data)[2] == 1){
          
          xx <- as.matrix(values$df_data)
          hist(xx, probability = T, breaks=20, main = 'Histogram of data')
        }
        else{
          pairs(values$df_data, lower.panel=NULL, upper.panel = upper.panel2, main = 'Paired scatter plot of data')
        }
        
      }
      else{
        
        if(dim(values$df_data)[2] == 1){
          
          mix_pdf <- function(x, para, k){
            pdf <- 0
            for (i in 1:k){
              pdf <- pdf + para$para$weight[i] * dnorm(x, para$para$mu[, i], para$para$sigma[[i]]^0.5)
            }
            return(pdf)
          }
          xx <- as.matrix(values$df_data)
          if (gmm_result()$success == F){
            xx <- as.matrix(values$df_data)
            hist(xx, probability = T, breaks=20, main = 'Histogram of data')
          }
          else{
            hist(xx, probability = T, breaks=20, main = 'Histogram of data')
            curve(mix_pdf(x, para= gmm_result(), k=input$num_k), from = min(xx), to = max(xx), col='red', add=T)
          }
        }
        
        else{
          if(gmm_result()$success == F){
            pairs(values$df_data, lower.panel=NULL, upper.panel = upper.panel2, main = 'Paired scatter plot of data')
          }
          else{
            pairs(values$df_data, lower.panel=NULL, upper.panel = upper.panel, main = 'Paired scatter plot of data')
          }
        }
        
        
        
      }
    }
    
    
  })
  
  output$plot3 <- renderPlot({ # BIC
    
    if (!is.null(values$df_data) & all(sapply(values$df_data, is.numeric))){
      
      lb <- input$bic_lb
      ub <- input$bic_ub
      
      if(is.na(ub) | ub <= 0 | ub%%1 != 0 | is.na(lb) | lb <= 0 | lb%%1 != 0 | lb >= ub){
        
        plot(0,type='n',axes=T, xlab = 'k component', ylab='bic score', main = 'BIC score')
        
        
      }else{
        
        print('bic plot start')
        
        n <- dim(values$df_data)[1]
        
        bic_matrix <- matrix(0, 5, ub-lb+1)
        
        unable_to_plot <- F
        
        for (t in 1:5){
          bic <- c()
          set.seed(1234* 0.5 *t)
          for (k in c(lb:ub)){
            
            result <- wrapped_gmm(values$df_data, k)
            if (result$success == F){
              unable_to_plot <- T
              break
            }
            lld <- wrapped_gmm(values$df_data, k)$lld
            bic <- c(bic, -2*lld + k*log(n))
            
          }
          if(unable_to_plot == T){
            break
          }
          bic_matrix[t ,] <- bic
        }
        if (unable_to_plot == F){
          df_bic <- data.frame(k=c(lb:ub), aic=colSums(bic_matrix)/5)
          plot(df_bic, type='b', xlab = 'k component', ylab='bic score', main = 'BIC score')
          print('bic plot complete')
        }else{
          plot(NULL, xlim=c(lb,ub), ylim=c(-1,1), xlab = 'k component', ylab='bic score', main='Error: unable to cluster!')
        }
        
      }
    }else{
      return(NULL)
    }
  })  
  
  
  output$plot4 <- renderPlot({ # AIC
    
    if (!is.null(values$df_data) & all(sapply(values$df_data, is.numeric))){
      
      lb <- input$aic_lb
      ub <- input$aic_ub
      
      if(is.na(ub) | ub <= 0 | ub%%1 != 0 | is.na(lb) | lb <= 0 | lb%%1 != 0 | lb >= ub){
        
        plot(0,type='n',axes=T, xlab = 'k component', ylab='aic score', main = 'AIC score')
        
        
      }else{
        
        print('aic plot start')
        
        n <- dim(values$df_data)[1]
        
        aic_matrix <- matrix(0, 5, ub-lb+1)
        
        unable_to_plot <- F
        
        for (t in 1:5){
          aic <- c()
          set.seed(1234* 0.5 *t)
          for (k in c(lb:ub)){
            result <- wrapped_gmm(values$df_data, k)
            if (result$success == F){
              unable_to_plot <- T
              break
            }

            lld <- result$lld
            aic <- c(aic, -2*lld + k*2)
            
          }
          if(unable_to_plot == T){
            break
          }
          aic_matrix[t ,] <- aic
        }
        if (unable_to_plot == F){
          df_aic <- data.frame(k=c(lb:ub), aic=colSums(aic_matrix)/5)
          plot(df_aic, type='b', xlab = 'k component', ylab='aic score', main = 'AIC score')
          print('aic plot complete')
        }else{
          plot(NULL, xlim=c(lb,ub), ylim=c(-1,1), xlab = 'k component', ylab='bic score', main='Error: unable to cluster!')
        }
      }
    }else{
      return(NULL)
    }
  })
  ##########
  
  # text messages
  
  ##########
  
  output$message2 <- renderPrint({
    
    if (is.null(values$df_data)){
      cat('Please choose a table to import.')
    }else if (!all(sapply(values$df_data, is.numeric))){
      cat('Data is not appropriate for the model.', '\n','Please re-upload.', sep = '')
    }
  })
  
  output$message3 <- renderPrint({
    
    if (is.null(values$df_data)){
      cat('Please choose a table to import.')
    }else if (!all(sapply(values$df_data, is.numeric))){
      cat('Data is not appropriate for the model.', '\n','Please re-upload.', sep='')
    }
  })
  
  output$message4 <- renderPrint({
    
    if (is.null(values$df_data)){
      cat('Please choose a table to import.')
    }else if (!all(sapply(values$df_data, is.numeric))){
      cat('Data is not appropriate for the model.', '\n','Please re-upload.', sep='')
    }
  })
  
  output$message5 <- renderPrint({ # BIC error message
    
    if (!is.null(values$df_data) & all(sapply(values$df_data, is.numeric))){
      if(is.na(input$bic_ub) | is.na(input$bic_lb)){
        cat('Please input upper bound and lower bound.')
      }
      else if (input$bic_ub <= input$bic_lb){
        cat('Upper bound must be greater than lower bound.')
      }
      else if(input$bic_ub%%1 != 0 | input$bic_lb%%1 != 0){
        cat('Upper and lower bound must be integers.')
      }
      else if (input$bic_ub <= 0 | input$bic_lb <= 0){
        cat('Upper and lower bound must be greater than 0.')
      }
      
    }
  })
  
  output$message6 <- renderPrint({ # AIC error message
    
    if (!is.null(values$df_data) & all(sapply(values$df_data, is.numeric))){
      if(is.na(input$aic_ub) | is.na(input$aic_lb)){
        cat('Please input upper bound and lower bound.')
      }
      else if (input$aic_ub <= input$aic_lb){
        cat('Upper bound must be greater than lower bound.')
      }
      else if(input$aic_ub%%1 != 0 | input$aic_lb%%1 != 0){
        cat('Upper and lower bound must be integers.')
      }
      else if (input$aic_ub <= 0 | input$aic_lb <= 0){
        cat('Upper and lower bound must be greater than 0.')
      }
      
    }
  })
  
  output$para1 <- renderPrint({
    
    if(!is.null(values$df_data) & all(sapply(values$df_data, is.numeric))){
      if(is.na(input$num_k)){
        cat('Please input k.')
      }else if (input$num_k <=0){
        cat('k must be greater than 0.')
      }else if(input$num_k %% 1 != 0){
        cat('k must be an integer.')
      }else{
        gmm_result()$para
      }
      
    }
  })
  
  
  #########
  
  output$select_example <- renderUI({
      req(input$use_example)
      selectInput(inputId = "selected_example",
                  label = 'Choose one',
                  choices =list('precip (1-dim)'=1, 'X1, X2 (2-dim)'=2, 'x, y, z (3-dim)'=3))
  })
  
  output$showfile <- renderUI({
    includeHTML('readme_final.html')
  })

})
