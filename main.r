# app.R
library(shiny)
library(imputeTS)
library(ggplot2)
library(Metrics)
library(reshape2)

ui <- fluidPage(
  titlePanel("MIMeC – Modelo de Interpolação Meteorológica Computacional"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Carregar CSV (uma coluna, vírgula decimal)", accept = ".csv"),
      checkboxGroupInput("methods", "Escolha os métodos:",
                         choices = c("Linear" = "linear",
                                     "Spline" = "spline",
                                     "Kalman" = "kalman",
                                     "STL (sazonal)" = "stl"),
                         selected = c("linear","kalman")),
      actionButton("run", "Rodar Imputação"),
      downloadButton("download", "Exportar resultados (CSV)"),
      br(),
      h4("Parâmetros automáticos utilizados:"),
      verbatimTextOutput("params")
    ),
    
    mainPanel(
      plotOutput("plot", height = "400px"),
      tableOutput("metrics"),
      br(),
      tableOutput("preview")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$run, {
    req(input$file)
    
    # 1) Ler CSV com uma coluna de temperatura, vírgula como decimal
    data <- read.csv(input$file$datapath, header = FALSE, dec = ",")
    x <- as.numeric(data[,1])
    
    if (all(is.na(x))) {
      showNotification("Arquivo inválido: todos os valores são NA.", type = "error")
      return(NULL)
    }
    
    n <- length(x)
    idx_known <- which(!is.na(x))
    
    # Seed automático = tamanho da série
    seed <- n
    set.seed(seed)
    
    # Porcentagem de validação adaptativa
    valpct <- if (n < 500) 5 else if (n <= 10000) 10 else 20
    n_val <- max(1, ceiling(valpct/100 * length(idx_known)))
    val_idx <- sample(idx_known, size = n_val)
    
    # Frequência sazonal automática (horária)
    dias <- floor(n/24)
    freq <- 24 * dias
    
    # Mostrar parâmetros na interface
    output$params <- renderText({
      paste0("Frequência sazonal (STL): ", freq, "\n",
             "Porcentagem para validação: ", valpct, "% (", n_val, " pontos mascarados)\n",
             "Seed (reprodutibilidade): ", seed)
    })
    
    x_masked <- x
    x_masked[val_idx] <- NA
    
    imputations <- list()
    metrics <- data.frame(Metodo = character(), MAE = numeric(), RMSE = numeric(), stringsAsFactors = FALSE)
    
    # Linear
    if ("linear" %in% input$methods) {
      x_lin <- na_interpolation(x_masked, option = "linear")
      imputations$Linear <- x_lin
      metrics <- rbind(metrics, data.frame(
        Metodo = "Linear",
        MAE = mae(x[val_idx], x_lin[val_idx]),
        RMSE = rmse(x[val_idx], x_lin[val_idx])
      ))
    }
    
    # Spline
    if ("spline" %in% input$methods) {
      x_spline <- na_interpolation(x_masked, option = "spline")
      imputations$Spline <- x_spline
      metrics <- rbind(metrics, data.frame(
        Metodo = "Spline",
        MAE = mae(x[val_idx], x_spline[val_idx]),
        RMSE = rmse(x[val_idx], x_spline[val_idx])
      ))
    }
    
    # Kalman
    if ("kalman" %in% input$methods) {
      x_kalman <- na_kalman(x_masked, model = "auto.arima")
      imputations$Kalman <- x_kalman
      metrics <- rbind(metrics, data.frame(
        Metodo = "Kalman",
        MAE = mae(x[val_idx], x_kalman[val_idx]),
        RMSE = rmse(x[val_idx], x_kalman[val_idx])
      ))
    }
    
    # STL
    if ("stl" %in% input$methods) {
      ts_x <- ts(x_masked, frequency = freq)
      x_stl <- na_seadec(ts_x, algorithm = "interpolation")
      x_stl <- as.numeric(x_stl)
      imputations$`STL (sazonal)` <- x_stl
      metrics <- rbind(metrics, data.frame(
        Metodo = "STL (sazonal)",
        MAE = mae(x[val_idx], x_stl[val_idx]),
        RMSE = rmse(x[val_idx], x_stl[val_idx])
      ))
    }
    
    # Plot
    output$plot <- renderPlot({
      df <- data.frame(Index = 1:n, Original = x, Mascarada = x_masked)
      for (m in names(imputations)) df[[m]] <- imputations[[m]]
      df_long <- melt(df, id.vars = "Index")
      
      ggplot(df_long, aes(x = Index, y = value, color = variable)) +
        geom_line(alpha = 0.9) +
        geom_point(data = subset(df_long, variable != "Original" & Index %in% val_idx),
                   size = 1.2, alpha = 0.7) +
        labs(title = "Comparação dos Métodos de Imputação (com validação)",
             y = "Temperatura", x = "Índice") +
        theme_minimal()
    })
    
    # Métricas
    output$metrics <- renderTable({
      metrics[order(metrics$MAE), ]
    }, digits = 4)
    
    # Prévia
    output$preview <- renderTable({
      df_out <- data.frame(Index = 1:n, Original = x, Mascarada = x_masked)
      for (m in names(imputations)) df_out[[m]] <- imputations[[m]]
      head(df_out, 12)
    }, digits = 3)
    
    # Exportar CSV
    output$download <- downloadHandler(
      filename = function() paste0("imputacoes_mimec_", Sys.Date(), ".csv"),
      content = function(file) {
        df_out <- data.frame(Index = 1:n, Original = x, Mascarada = x_masked)
        for (m in names(imputations)) df_out[[m]] <- imputations[[m]]
        write.csv(df_out, file, row.names = FALSE)
      }
    )
  })
}

shinyApp(ui = ui, server = server)
