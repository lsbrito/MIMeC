# app.R
library(shiny)
library(imputeTS)
library(ggplot2)
library(Metrics)
library(reshape2)

ui <- fluidPage(
  titlePanel("Preenchimento de Dados de Temperatura - INMET"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Carregar CSV (uma coluna, vírgula decimal)", accept = ".csv"),
      checkboxGroupInput("methods", "Escolha os métodos:",
                         choices = c("Linear" = "linear",
                                     "Spline" = "spline",
                                     "Kalman" = "kalman",
                                     "STL (sazonal)" = "stl"),
                         selected = c("linear","kalman")),
      numericInput("freq", "Frequência sazonal (ex.: diária=365, horária=24*365)", value = 365, min = 1),
      sliderInput("valpct", "Porcentagem para validação (mascarar):", min = 5, max = 30, value = 10, step = 1),
      numericInput("seed", "Seed (reprodutibilidade)", value = 123, min = 1),
      actionButton("run", "Rodar Imputação"),
      downloadButton("download", "Exportar resultados (CSV)")
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
    # Não definir sep, pois há apenas uma coluna
    data <- read.csv(input$file$datapath, header = FALSE, dec = ",")
    x <- as.numeric(data[,1])
    
    # Verificações simples
    if (all(is.na(x))) {
      showNotification("Arquivo inválido: todos os valores são NA. Verifique o formato (uma coluna, vírgula decimal).", type = "error")
      return(NULL)
    }
    
    n <- length(x)
    idx_known <- which(!is.na(x))
    if (length(idx_known) < 10) {
      showNotification("Poucos pontos conhecidos. Adicione mais dados ou reduza a porcentagem de validação.", type = "warning")
    }
    
    # 2) Criar máscara de validação: retirar aleatoriamente uma porcentagem dos pontos conhecidos
    set.seed(input$seed)
    n_val <- max(1, ceiling(input$valpct/100 * length(idx_known)))
    val_idx <- sample(idx_known, size = n_val)
    
    x_masked <- x
    x_masked[val_idx] <- NA
    
    imputations <- list()
    metrics <- data.frame(Metodo = character(), MAE = numeric(), RMSE = numeric(), stringsAsFactors = FALSE)
    
    # 3) Aplicar métodos selecionados
    if ("linear" %in% input$methods) {
      x_lin <- tryCatch(na_interpolation(x_masked, option = "linear"), error = function(e) rep(NA_real_, n))
      imputations$Linear <- x_lin
      # Métricas calculadas somente nos pontos mascarados (val_idx)
      metrics <- rbind(metrics, data.frame(
        Metodo = "Linear",
        MAE = mae(x[val_idx], x_lin[val_idx]),
        RMSE = rmse(x[val_idx], x_lin[val_idx])
      ))
    }
    
    if ("spline" %in% input$methods) {
      x_spline <- tryCatch(na_interpolation(x_masked, option = "spline"), error = function(e) rep(NA_real_, n))
      imputations$Spline <- x_spline
      metrics <- rbind(metrics, data.frame(
        Metodo = "Spline",
        MAE = mae(x[val_idx], x_spline[val_idx]),
        RMSE = rmse(x[val_idx], x_spline[val_idx])
      ))
    }
    
    if ("kalman" %in% input$methods) {
      x_kalman <- tryCatch(na_kalman(x_masked, model = "auto.arima"), error = function(e) rep(NA_real_, n))
      imputations$Kalman <- x_kalman
      metrics <- rbind(metrics, data.frame(
        Metodo = "Kalman",
        MAE = mae(x[val_idx], x_kalman[val_idx]),
        RMSE = rmse(x[val_idx], x_kalman[val_idx])
      ))
    }
    
    if ("stl" %in% input$methods) {
      # STL precisa de um ts com frequência definida
      ts_x <- ts(x_masked, frequency = input$freq)
      x_stl <- tryCatch(na_seadec(ts_x, algorithm = "interpolation"), error = function(e) rep(NA_real_, n))
      x_stl <- as.numeric(x_stl)
      imputations$`STL (sazonal)` <- x_stl
      metrics <- rbind(metrics, data.frame(
        Metodo = "STL (sazonal)",
        MAE = mae(x[val_idx], x_stl[val_idx]),
        RMSE = rmse(x[val_idx], x_stl[val_idx])
      ))
    }
    
    # 4) Plot: original, mascarada e imputações
    output$plot <- renderPlot({
      df <- data.frame(Index = 1:n, Original = x, Mascarada = x_masked)
      for (m in names(imputations)) df[[m]] <- imputations[[m]]
      df_long <- melt(df, id.vars = "Index")
      
      ggplot(df_long, aes(x = Index, y = value, color = variable)) +
        geom_line(alpha = 0.9) +
        # Destacar pontos mascarados e estimados
        geom_point(data = subset(df_long, variable != "Original" & Index %in% val_idx),
                   size = 1.2, alpha = 0.7) +
        labs(title = "Comparação dos Métodos de Imputação (com validação)",
             y = "Temperatura", x = "Índice") +
        theme_minimal()
    })
    
    # 5) Tabela de métricas
    output$metrics <- renderTable({
      metrics[order(metrics$MAE), ]
    }, digits = 4)
    
    # 6) Prévia dos dados imputados
    output$preview <- renderTable({
      df_out <- data.frame(Index = 1:n, Original = x, Mascarada = x_masked)
      for (m in names(imputations)) df_out[[m]] <- imputations[[m]]
      head(df_out, 12)
    }, digits = 3)
    
    # 7) Exportar CSV com imputações
    output$download <- downloadHandler(
      filename = function() paste0("imputacoes_inmet_", Sys.Date(), ".csv"),
      content = function(file) {
        df_out <- data.frame(Index = 1:n, Original = x, Mascarada = x_masked)
        for (m in names(imputations)) df_out[[m]] <- imputations[[m]]
        write.csv(df_out, file, row.names = FALSE)
      }
    )
  })
}

shinyApp(ui = ui, server = server)
