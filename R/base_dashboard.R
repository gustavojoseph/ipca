################################################################################
############################## MODELO PARA O IPCA ##############################
################################################################################

# PACOTES -----------------------------------------------------------------

library(tidyverse)
library(rbcb)
library(ipeadatar)
library(sidrar)
library(gtrendsR)
library(tsibble)
library(forecast)
library(feasts)
library(fabletools)
library(timetk)
library(caret)
library(glmnet)
library(plotly)
library(scales)
library(ggtext)
library(DT)
library(rmarkdown)
library(flexdashboard)

# Pacotes secundários

library(bit)
library(class)
library(codetools)
library(cpp11)
library(datawizard)
library(HDeconometrics)
library(jsonlite)
library(KernSmooth)
library(lattice)
library(MASS)
library(Matrix)
library(mgcv)
library(nlme)
library(nnet)
library(readxl)
library(rpart)
library(sparsevctrs)
library(survival)

#remotes::install_github("gabrielrvsc/HDeconometrics")
library(HDeconometrics)

# FUNÇÕES -----------------------------------------------------------------

#1ª Função

report_ndiffs2 <- function(
    x,
    test = c("kpss", "adf", "pp"),
    term = c("level", "trend"),
    alpha = 0.05,
    na_rm = TRUE
) {

  # Garantir que x seja lista de séries
  if (!is.list(x)) x <- list(x)

  # Remover NAs se necessário
  if (na_rm) x <- purrr::map(x, stats::na.omit)

  # Nomes para variáveis se não existirem
  if (is.null(names(x))) names(x) <- paste0("var", seq_along(x))

  # Definir combinações de testes e tipos
  test_grid <- tidyr::expand_grid(
    test = test,
    type = term
  )

  # Aplicar testes em cada série
  resultados <- purrr::map_dfr(
    .x = x,
    .f = function(serie) {
      purrr::pmap_dfc(
        test_grid,
        function(test, type) {
          ndiff <- forecast::ndiffs(
            serie,
            alpha = alpha,
            test = test,
            type = type
          )
          tibble::tibble(!!paste(test, type, sep = "_") := ndiff)
        }
      )
    },
    .id = "variable"
  )

  # Calcular o ndiffs mais frequente por série
  resultados <- resultados |>
    dplyr::rowwise() |>
    dplyr::mutate(
      ndiffs = as.numeric(
        names(
          sort(
            table(dplyr::c_across(-variable)),
            decreasing = TRUE
          )[1]
        )
      )
    ) |>
    dplyr::ungroup()

  return(resultados)
}

#2ª Função

#' Format x-axislabels (dates) as year/month in two lines
#'
#' @param x Date vector, usually takes as input the output of `breaks` in `ggplot2::scale_x_date`
#'
#' @return
#' @export
#'
#' @examples
ym_label <- function(x) {

  x <- lubridate::as_date(x)

  dplyr::if_else(
    is.na(dplyr::lag(x)) | tsibble::yearmonth(dplyr::lag(x)) != tsibble::yearmonth(x),
    paste(lubridate::month(x, label = TRUE), "\n", lubridate::year(x)),
    paste(lubridate::month(x, label = TRUE))
  )

}

#3ª Função

#' CSR and Bagging: estimates model with cross-validation and reports accuracy by forecast horizon
#'
#' @param model Model to be estimated, possible values are `csr` or `bagging`, see HDeconometrics functions
#' @param data A data frame
#' @param y_target Column name of the variable of interest (used to report accuracy)
#' @param date_col Date class column name
#' @param init_window Number of initial observations to be used in the first cross-validation subsample
#' @param step A positive integer for incremental step (see `tsibble::stretch_tsibble`)
#' @param horizon Forecast horizon
#' @param ... Additional arguments to `HDeconometrics::csr` or `HDeconometrics::bagging`
#'
#' @return tibble with the RMSE per forecast horizon.
get_cv_rmse_hdecon <- function (
    model,
    data,
    y_target,
    date_col,
    init_window = 150,
    step        = 1,
    horizon     = 12,
    ...
) {

  cv_train_index <- data |>
    dplyr::slice(1:(dplyr::n() - horizon)) |>
    nrow() |>
    seq_len() |>
    # function not exported, use with caution!
    tsibble:::stretcher2(.step = step, .init = init_window)

  n_fcst <- length(cv_train_index)

  point_fcst <- list()

  for (i in seq_len(n_fcst)) {

    cat(paste0("\nIteration: ", i, "/", n_fcst))

    curr_index <- cv_train_index[[i]]
    data_train <- data[curr_index, ]

    yy_in <- dplyr::pull(data_train, dplyr::all_of(y_target))
    xx_in <- dplyr::select(data_train, !dplyr::any_of(c(date_col, y_target)))

    xx_out <- as.matrix(data[-curr_index, names(xx_in)][1:horizon, ])

    if (model == "csr") {

      fit_csr <- HDeconometrics::csr(
        x = xx_in,
        y = yy_in,
        ...
      )

      fcsts <- predict(object = fit_csr, newdata = xx_out)

    } else if (model == "bagging") {

      fit_bagging <- HDeconometrics::bagging(
        x = as.matrix(xx_in),
        y = yy_in,
        ...
      )

      fcsts <- predict(object = fit_bagging, newdata = xx_out)

    } else stop("model must be 'csr' or 'bagging'.")

    point_fcst[[i]] <- dplyr::tibble(
      {{ date_col }} := seq.Date(
        from       = max(dplyr::pull(data_train, {{ date_col }})) + months(1),
        by         = "month",
        length.out = length(fcsts)
      ),
      fcst = fcsts
    )

  }

  fc <- point_fcst |>
    dplyr::bind_rows(.id = ".id") |>
    dplyr::mutate(model = dplyr::if_else(model == "csr", "csr", "bagging"))

  rmse_tbl <- dplyr::left_join(
    x  = fc,
    y  = dplyr::select(data, dplyr::all_of(c(date_col, y_target))),
    by = {{ date_col }}
  ) |>
    dplyr::group_by(.id) |>
    dplyr::mutate(h = dplyr::row_number()) |>
    dplyr::group_by(h, model) |>
    dplyr::summarise(
      rmse = sqrt(mean((!!rlang::sym(y_target) - fcst)^2, na.rm = TRUE)),
      .groups = "drop"
    )

  return(
    list("rmse" = rmse_tbl, "forecasts" = fc)
  )

}

# COLETA DE DADOS ---------------------------------------------------------

# Metadados/códigos de coleta (CSV)
metadata <- readr::read_csv2(
  file      = "./R/metadados.csv",
  col_names = TRUE,
  col_types = "c"
)

# Data inicial da coleta de dados
init_date <- lubridate::ymd("2002-12-01")

# 1. Dados do Banco Central ----

# Códigos de coleta
codes_bcb <- metadata |>
  dplyr::filter(fonte == "BCB") |>
  dplyr::reframe(
    purrr::set_names(x = as.numeric(codigo), nm = acronimo)
  ) |>
  dplyr::pull()


# Coleta de dados do SGS/BCB
raw_bcb <- purrr::map(codes_bcb, function(code) {
  tryCatch({
    rbcb::get_series(
      code       = code,
      start_date = init_date,
      end_date   = lubridate::today()
    )
  }, error = function(e) {
    message("Erro ao baixar a série ", code, ": ", e$message)
    return(NULL)  # Retorna NULL em caso de erro
  })
}) |>
  purrr::discard(is.null)

# 2. Dados do IPEA ----

# Códigos de coleta
codes_ipeadata <- metadata |>
  dplyr::filter(fonte == "IPEADATA") |>
  dplyr::reframe(
    purrr::set_names(x = codigo, nm = acronimo)
  ) |>
  dplyr::pull()

# Coleta de dados do IPEADATA
raw_ipeadata <- ipeadatar::ipeadata(code = codes_ipeadata)

# 3. Dados do IBGE ----

# Códigos de coleta
codes_ibge <- metadata |>
  dplyr::filter(fonte == "IBGE") |>
  dplyr::reframe(
    purrr::set_names(x = codigo, nm = acronimo)
  ) |>
  dplyr::pull()

# Coleta de dados do IBGE
raw_ibge <- purrr::map(
  .x = codes_ibge,
  .f = ~sidrar::get_sidra(api = .x)
)

# 4. Dados do Google Trends ----

# Códigos de coleta
codes_google <- metadata |>
  dplyr::filter(fonte == "Google Trends") |>
  dplyr::reframe(
    purrr::set_names(x = codigo, nm = acronimo)
  ) |>
  dplyr::pull()


# Coleta de dados do Google Trends
raw_google <- gtrendsR::gtrends(
  keyword      = codes_google,
  geo          = "BR",
  time         = "all",
  onlyInterest = TRUE
)


# 5. Dados do Noletim Focus/BCB ----

# Códigos de coleta
codes_focus <- metadata |>
  dplyr::filter(fonte == "Focus/BCB") |>
  dplyr::reframe(
    purrr::set_names(x = codigo, nm = acronimo)
  ) |>
  dplyr::pull()


# Coleta de dados do Focus/BCB
raw_focus <- rbcb::get_market_expectations(
  type       = "monthly",
  indic      = codes_focus,
  start_date = init_date,
  end_date   = lubridate::today()
)

# Salvar dados brutos para reprodução
readr::write_rds(x = mget(ls(pattern = "raw_")), file = "raw_data.rds")

# TRATAMENTO DE DADOS -----------------------------------------------------

# Dados do SGS/BCB
df_bcb <- raw_bcb |> purrr::reduce(.f = dplyr::full_join, by = "date")

# Dados do IPEADATA
df_ipeadata <- raw_ipeadata |>
  dplyr::select("date", "variable" = "code", "value") |>
  dplyr::left_join(
    y  = dplyr::select(metadata, "codigo", "acronimo"),
    by = c("variable" = "codigo")
  ) |>
  tidyr::pivot_wider(
    id_cols     = "date",
    names_from  = "acronimo",
    values_from = "value"
  ) |>
  # Obter média mensal (para séries com freq. diária)
  dplyr::group_by(date = tsibble::yearmonth(.data$date)) |>
  dplyr::summarise(
    dplyr::across(where(is.numeric), ~mean(.x, na.rm = TRUE))
  ) |>
  dplyr::mutate(date = lubridate::as_date(.data$date))


# Dados do IBGE
df_ibge <- raw_ibge |>
  purrr::map_dfr(
    .f  = ~dplyr::select(.x, "date" = "Mês (Código)", "value" = "Valor"),
    .id = "variable"
  ) |>
  tidyr::pivot_wider(
    id_cols     = "date",
    names_from  = "variable",
    values_from = "value"
  ) |>
  dplyr::mutate(date = lubridate::ym(.data$date)) |>
  dplyr::filter(date >= init_date) |>
  dplyr::relocate("date", "ipca")


# Dados do Google
df_google <- raw_google |>
  purrr::pluck("interest_over_time") |>
  dplyr::mutate(
    date    = lubridate::as_date(.data$date),
    keyword = stringr::str_replace(keyword, " ", "_")
  ) |>
  tidyr::pivot_wider(
    id_cols     = "date",
    names_from  = "keyword",
    values_from = "hits"
  )


# Dados do Focus/BCB
df_focus <- raw_focus |>
  dplyr::arrange(Data) |>
  # Calcula horizonte em meses da expectativa
  dplyr::mutate(
    monthyear = tsibble::yearmonth(Data),
    horizon   = tsibble::yearmonth(DataReferencia, format = "%m/%Y") - monthyear
  ) |>
  # Agrupar por mês de estatísticas registradas no Focus
  dplyr::group_by(monthyear) |>
  # Filtra estatísticas tipo 0 (últimos 30 dias) no horizonte de 1 ano e
  # de datas próxima ao dia 15
  dplyr::filter(
    baseCalculo == 0,
    horizon > 0 & horizon < 13,
    lubridate::day(Data) < 16
  ) |>
  dplyr::filter(
    abs(lubridate::day(Data) - 15) == min(abs(lubridate::day(Data) - 15))
  ) |>
  dplyr::ungroup() |>
  dplyr::distinct(
    "date" = lubridate::floor_date(x = Data, unit = "month"),
    horizon,
    .keep_all = TRUE
  ) |>
  tidyr::pivot_wider(
    id_cols      = date,
    names_from   = horizon,
    values_from  = Mediana,
    names_prefix = "expectativa_ipca_h_"
  )


# Cruzamento de dados
df_ipca <- purrr::reduce(
  .x = list(df_ibge, df_bcb, df_ipeadata, df_google, df_focus),
  .f = dplyr::left_join,
  by = "date"
)


# 1. Teste de Estacionariedade ----

# Aplicar testes de estacionariedade
vars_ndiffs <- df_ipca |>
  dplyr::select(-"date") |>
  report_ndiffs2()


# Diferenciar séries que são não estacionárias
df_ipca_diff <- df_ipca |>
  dplyr::mutate(dplyr::across(
    .cols = vars_ndiffs$variable[vars_ndiffs$ndiffs > 0 & vars_ndiffs$variable != "ipca"],
    .fns  = ~tsibble::difference(
      x           = .x,
      differences = vars_ndiffs$ndiffs[vars_ndiffs$variable == dplyr::cur_column()]
    )
  )
  )

# 2. Valores Ausentes ----

# Expandir base criando defasagens e dummies, preenchendo valores NA (de baixo),
# filtrar amostra e selecionar variáveis de interesse
df_ipca_lagged <- df_ipca_diff |>
  timetk::tk_augment_lags(.value = !dplyr::any_of("date"), .lags = 1:4) |> #defasagens
  tidyr::fill(!dplyr::any_of(c("date", "ipca")), .direction = "down") |> #Na's de baixo
  tidyr::drop_na() |>
  dplyr::select(
    !dplyr::any_of(                           # filtrando variáveis de interesse
      stringr::str_subset(
        string  = names(df_ipca_diff),
        pattern = "date|ipca|expectativa",
        negate  = TRUE
      )
    )
  )

# Criando um vetor time-series
yy_ts <- stats::ts(
  data = df_ipca_lagged$ipca,
  start = c(
    lubridate::year(min(df_ipca_lagged$date)),
    lubridate::month(min(df_ipca_lagged$date))
  ),
  frequency = 12
)

# Criando um vetor de dummies sazonais
seasonal_dummies <- forecast::seasonaldummy(yy_ts)

# Incluindo as dummies na base de dados
df_ipca_lagged <- df_ipca_lagged |>
  dplyr::bind_cols(seasonal_dummies)

# VALIDAÇÃO CRUZADA -------------------------------------------------------

# Número de observações iniciais
init_obs <- 100

# Horizonte de previsão
horizon <- 12
t1 <- 1

# 1. Modelo CSR ----

# Aplica função criada para reportar RMSE por horizonte preditivo e pontos de previsão
acc_csr <- get_cv_rmse_hdecon(
  model          = "csr",
  data           = df_ipca_lagged,
  y_target       = "ipca",
  date_col       = "date",
  init_window    = init_obs,
  step           = 1,
  horizon        = horizon,
  K              = 20,
  k              = 15,
  fixed.controls = colnames(seasonal_dummies)
)

acc_csr_t1 <- get_cv_rmse_hdecon(
  model          = "csr",
  data           = df_ipca_lagged,
  y_target       = "ipca",
  date_col       = "date",
  init_window    = init_obs,
  step           = 1,
  horizon        = t1,
  K              = 20,
  k              = 15,
  fixed.controls = colnames(seasonal_dummies)
)

# Tratar dados para obter previsões por amostra de validação cruzada do modelo
df_csr <- acc_csr$forecasts |>
  tidyr::pivot_wider(
    id_cols     = c(".id", "date"),
    names_from  = "model",
    values_from = "fcst"
  )

df_csr_t1 <- acc_csr_t1$forecasts |>
  tidyr::pivot_wider(
    id_cols     = c(".id", "date"),
    names_from  = "model",
    values_from = "fcst"
  )


# 2. RW ----


# Calcular acurácia do modelo RW
acc_rw <- df_csr |>
  dplyr::left_join(
    y  = dplyr::select(df_ipca_lagged, "date", "ipca", "ipca_lag1"),
    by = "date"
  ) |>
  dplyr::group_by(.id) |>
  dplyr::mutate(h = dplyr::row_number()) |>
  dplyr::group_by(h) |>
  dplyr::summarise(
    rw       = sqrt(mean((ipca - ipca_lag1)^2, na.rm = TRUE)),
    .groups  = "drop"
  ) |>
  tidyr::pivot_longer(cols = -"h", names_to = "model", values_to = "rmse")



# 3. Acurácia do Focus/BCB ----

# Calcular acurácia do Focus por horizonte preditivo
acc_focus <- df_ipca_lagged |>
  dplyr::select(
    "date",
    "ipca",
    dplyr::matches("expectativa_ipca_h_\\d{1,2}$")
  ) |>
  tidyr::pivot_longer(
    cols      = -c("date", "ipca"),
    names_to  = "h",
    values_to = "focus"
  ) |>
  dplyr::mutate(
    h     = as.integer(stringr::str_remove(h, "expectativa_ipca_h_")),
    model = "focus"
  ) |>
  dplyr::left_join(
    y = df_csr |>
      dplyr::group_by(.id) |>
      dplyr::mutate(h = dplyr::row_number()) |>
      dplyr::select(-csr),
    by = c("date", "h")
  ) |>
  tidyr::drop_na() |>
  dplyr::group_by(h, model) |>
  dplyr::summarise(
    rmse = sqrt(mean((ipca - focus)^2, na.rm = TRUE)),
    .groups  = "drop"
  )


# AVALIAÇÃO DA ACURÁCIA DOS MODELOS ---------------------------------------

# Gerar gráfico de linha do RMSE por horizonte preditivo de cada modelo
plotly::ggplotly(

  dplyr::bind_rows(
    acc_csr$rmse,
    acc_rw,
    acc_focus
  ) |>
    ggplot2::ggplot() +

    ggplot2::aes(x = h, y = rmse, color = model) +

    scale_x_continuous(n.breaks = 10) +

    ggplot2::geom_line() +

    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.caption = element_text(hjust = 0, face = "bold"),
      plot.caption.position = "plot",
      legend.position = "top") +

    labs(
      title    = "Avaliando a Acurácia dos Modelos",
      y        = "RMSE",
      x        = NULL,
      color    = "Modelos",
      caption  = NULL)
)

# PREVISÃO FORA DA AMOSTRA ------------------------------------------------

# Criar objetos com vetor da variável dependente e matriz das independentes
yy <- dplyr::pull(df_ipca_lagged, "ipca")
xx <- df_ipca_lagged |>
  dplyr::select(!dplyr::any_of(c("date", "ipca"))) |>
  as.matrix()


# Estimar modelo CSR
fit_csr <- HDeconometrics::csr(
  x              = xx,
  y              = yy,
  K              = 20,
  k              = 15,
  fixed.controls = colnames(seasonal_dummies)
)

# 1. Avaliação dos resíduos: CSR ----
resids_csr <- fit_csr |>
  purrr::pluck("residuals") |>
  dplyr::as_tibble(.name_repair = ~"resid") |>
  dplyr::mutate(
    date = seq.Date(
      from = min(df_ipca_lagged$date),
      by   = "month",
      to   = max(df_ipca_lagged$date)
    ) |> tsibble::yearmonth()
  )  |>
  tsibble::as_tsibble(index = "date")

# Autocorrelação
resids_csr |> # correlograma ACF
  feasts::ACF(resid) |>
  fabletools::autoplot()

# Normalidade
resids_csr |> # histograma
  ggplot2::ggplot() +
  ggplot2::aes(x = resid) +
  ggplot2::geom_histogram()


# 2. Cenários para variáveis independentes ----

# Modelos univariados - auto ARIMA
xreg_arima <- purrr::map_df(
  .x = dplyr::select(
    df_ipca_lagged,
    !dplyr::any_of(c("date", "ipca", colnames(seasonal_dummies)))
  ),
  .f = ~{
    forecast::auto.arima(.x) |>
      forecast::forecast(h = horizon) |>
      magrittr::extract2("mean") |>
      as.numeric()
  }
)

# Dummies sazonais de fora da amostra
seasonal_dummies_oos <- forecast::seasonaldummy(x = yy_ts, h = horizon)


# Variáveis exógenas para gerar previsões fora da amostra
xx_oos <- as.matrix(dplyr::bind_cols(xreg_arima, seasonal_dummies_oos))

# 3. Previsões fora da amostra ----

# Previsão CSR
fcst_csr <- predict(object = fit_csr, newdata = xx_oos)

# Tabela
previsao <- dplyr::tibble(csr = fcst_csr) |>
  dplyr::mutate(
    date = seq.Date(
      from       = max(df_ipca_lagged$date) + months(1),
      by         = "month",
      length.out = length(fcst_csr)
    ),
    .before = 1
  )

# VISUALIZAÇÃO ------------------------------------------------------------

# Cores para gráfico
colors <- c(
  "#325d88", # blue
  "#289386", # green
  "#fc9849", # red
  "#e8bb4f", # yellow
  "black"
)

# Juntar dados observados com pontos de previsão
df_fanchart <- df_ipca_lagged |>
  dplyr::select("date", "ipca") |>
  dplyr::full_join(y = previsao, by = "date") |>
  dplyr::mutate(
    csr = dplyr::if_else(date == max(df_ipca_lagged$date), ipca, csr)
  )


# 1. Gráfico Previsão ----
grafico_previsao <- df_fanchart |>

  ggplot2::ggplot() +

  ggplot2::aes(x = date) +

  ggplot2::geom_line(
    mapping = ggplot2::aes(y = ipca),
    size    = 2,
    color   = colors[5],
    data    = dplyr::slice_tail(df_ipca_lagged, n = 36)
  ) +

  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(ymin = -Inf, ymax = Inf),
    fill    = colors[1],
    alpha   = 0.35,
    data    = dplyr::filter(df_fanchart, date >= max(df_ipca_lagged$date))
  ) +

  ggplot2::geom_line(
    mapping = ggplot2::aes(y = csr),
    size    = 2,
    color   = colors[1],
    data    = dplyr::filter(df_fanchart, date >= max(df_ipca_lagged$date))
  ) +

  ggplot2::geom_vline(
    xintercept = max(df_ipca_lagged$date),
    linetype   = "dashed"
  ) +

  ggplot2::scale_y_continuous(
    n.breaks = 6,
    labels = scales::number_format(
      suffix       = "%",
      accuracy     = 0.1,
      decimal.mark = ","
    )
  ) +

  ggplot2::scale_x_date(
    breaks = scales::breaks_width("4 months"),
    labels = ym_label
  ) +

  ggplot2::theme_light() +
  ggplot2::labs(
    title    = "Previsão para os próximos 12 meses",
    subtitle = "Indice de Preços ao Consumidor Amplo - IPCA",
    y        = "% a.m.",
    x        = NULL,
    caption  = "**Modelo:** Complete Subset Regression - CSR"
  ) +

  ggplot2::theme(
    plot.title       = ggtext::element_markdown(size = 20, colour = colors[1]),
    plot.subtitle    = ggtext::element_markdown(size = 14),
    axis.text        = ggtext::element_markdown(size = 12, face = "bold"),
    axis.title       = ggtext::element_markdown(size = 12, face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    plot.caption     = ggtext::element_textbox_simple(
      size   = 10,
      colour = "grey20",
      margin = ggplot2::margin(10, 5.5, 10, 5.5)
    )
  )


# 2. Gráfico de Acurácia ----

# Juntar dados de acurácia de todos os modelos
acc_rmse <- dplyr::bind_rows(
  acc_csr$rmse,
  acc_rw,
  acc_focus
) |>
  dplyr::mutate(
    model = dplyr::recode(
      model,
      "csr"      = "Complete Subset Regression",
      "rw"       = "Random Walk",
      "focus"    = "Focus"
    )
  ) |>
  dplyr::arrange(h, rmse, model) |>
  dplyr::select("Horizonte" = "h", "Modelo" = "model", "RMSE" = "rmse")


# Gráfico do RMSE por horizonte de previsão
grafico_rmse <- acc_rmse |>
  ggplot2::ggplot(
    ggplot2::aes(x = Horizonte, y = RMSE, colour = Modelo)
  ) +
  ggplot2::geom_line(size = 1.5) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_colour_manual(values = colors[1:3]) +
  ggplot2::scale_y_continuous(
    breaks = scales::breaks_extended(n = 8),
    labels = scales::number_format(
      accuracy     = 0.001,
      decimal.mark = ",",
      big.mark     = "."
    )
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::number_format(accuracy = 1),
    breaks = 1:horizon
  ) +
  ggplot2::theme_light() +
  ggplot2::labs(
    title    = "Performance por horizonte preditivo",
    subtitle = "Modelos de previsão do IPCA",
    x        = "Horizonte (meses)",
    color    = NULL,
    caption  = NULL
  ) +
  ggplot2::theme(
    plot.title       = ggtext::element_markdown(size = 20, colour = colors[1]),
    plot.subtitle    = ggtext::element_markdown(size = 14),
    axis.text        = ggtext::element_markdown(size = 12, face = "bold"),
    axis.title       = ggtext::element_markdown(size = 12, face = "bold"),
    legend.position  = "bottom",
    legend.text      = ggplot2::element_text(size = 12, face = "bold"),
    legend.key.width = ggplot2::unit(1, "cm"),
    panel.grid.minor = ggplot2::element_blank(),
    plot.caption     = ggtext::element_textbox_simple(
      size   = 12,
      colour = "grey20",
      margin = ggplot2::margin(10, 5.5, 10, 5.5)
    )
  )

# 3. Gráfico últimos 12 meses ----

# Tabelas
tabela_focus <- df_ipca_lagged |>
  dplyr::select(
    "date",
    "ipca",
    dplyr::matches("expectativa_ipca_h_\\d{1,2}$")
  ) |>
  tidyr::pivot_longer(
    cols      = -c("date", "ipca"),
    names_to  = "h",
    values_to = "focus"
  ) |>
  dplyr::mutate(
    h     = as.integer(stringr::str_remove(h, "expectativa_ipca_h_")),
    model = "focus"
  ) |>
  dplyr::filter( h == 1) |>
  dplyr:: mutate(date  = date,
                 IPCA  = ipca,
                 Focus = round(lag(focus), 2),
                 .keep = "none") |>
  dplyr::arrange(desc(date))


tabela_modelo <- df_csr_t1 |>
  dplyr::mutate(date   = date,
                Modelo = round(lag(csr), 2),
                .keep  = "none") |>
  dplyr::arrange(desc(date))

# Gráfico
grafico_12m <- tabela_focus |>
  dplyr::left_join(y = tabela_modelo,
                   by = "date") |>
  tidyr::drop_na() |>
  dplyr::mutate(Data   = date,
                IPCA   = IPCA,
                Focus  = Focus,
                Modelo = Modelo,
                .keep  = "none"
  ) |>
  dplyr::select(4,1,3,2) |>
  tidyr::pivot_longer(
    cols      = 2:4,
    names_to  = "Modelo",
    values_to = "Valor"
  ) |>
  head(36) |>

  ggplot2::ggplot() +

  ggplot2::aes(x = Data, y = Valor, colour = Modelo) +

  ggplot2::geom_line(size = 1.5) +

  ggplot2::geom_point(size = 3) +

  ggplot2::scale_colour_manual(values = c(
    "Modelo" = colors[1],
    "IPCA"   = colors[5],
    "Focus"  = colors[2])
  ) +

  ggplot2::scale_y_continuous(
    n.breaks = 6,
    labels   = scales::number_format(
      suffix       = "%",
      accuracy     = 0.1,
      decimal.mark = ","
    )) +

  ggplot2::scale_x_date(
    breaks = scales::breaks_width("1 months"),
    labels = ym_label
  ) +

  ggplot2::theme_light() +

  ggplot2::labs(
    title    = "Análise dos últimos 12 meses",
    subtitle = "Comparação do IPCA, Previsão do Modelo e Expectativa de Mercado",
    x        = NULL,
    y        = "% a.m.",
    color    = NULL,
    caption  = NULL
  ) +

  ggplot2::theme(
    plot.title       = ggtext::element_markdown(size = 20, colour = colors[1]),
    plot.subtitle    = ggtext::element_markdown(size = 14),
    axis.text        = ggtext::element_markdown(size = 12, face = "bold"),
    axis.title       = ggtext::element_markdown(size = 12, face = "bold"),
    legend.position  = "bottom",
    legend.text      = ggplot2::element_text(size = 12, face = "bold"),
    legend.key.width = ggplot2::unit(1, "cm"),
    panel.grid.minor = ggplot2::element_blank(),
    plot.caption     = ggtext::element_textbox_simple(
      size   = 12,
      colour = "grey20",
      margin = ggplot2::margin(10, 5.5, 10, 5.5)
    )
  )

# 4. Tabelas ----

# Tabela de pontos de previsão
tabela_previsao <- df_fanchart |>
  dplyr::filter(date > max(df_ipca_lagged$date)) |>
  dplyr::mutate(date = lubridate::as_date(date) |> format("%b/%Y")) |>
  dplyr::select("Mes" = "date", "Previsao" = "csr") |>
  DT::datatable(
    options = list(dom = "tip", pageLength = 10, scrollX = TRUE, scrollY = TRUE),
    rownames = FALSE
  ) |>
  DT::formatRound(columns = 2, digits = 2, dec.mark = ",", mark = ".") |>
  DT::formatStyle(columns = 2, fontWeight = "bold")


# Tabela com valores do RMSE vs. horizonte/modelos
tabela_rmse <- acc_rmse |>
  dplyr::mutate(Horizonte = as.character(Horizonte)) |>
  DT::datatable(
    options = list(dom = "tip", pageLength = 10, scrollX = TRUE, scrollY = TRUE),
    rownames = FALSE
  ) |>
  DT::formatRound(columns = 3, digits = 2, dec.mark = ",", mark = ".")

# Tabelas Modelo x Focus
tabela_versus <- tabela_focus |>
  dplyr::left_join(y = tabela_modelo,
                  by = "date") |>
  tidyr::drop_na() |>
  dplyr::mutate(Data   = lubridate::as_date(date) |> format("%b/%Y"),
                IPCA   = IPCA,
                Focus  = Focus,
                Modelo = Modelo,
                .keep  = "none"
  ) |>
  dplyr::select(4,1,3,2) |>
  DT::datatable(
    options = list(dom = "tip", pageLength = 10, scrollX = TRUE, scrollY = TRUE),
    rownames = FALSE
  ) |>
  DT::formatRound(columns = 2, digits = 2, dec.mark = ",", mark = ".") |>
  DT::formatStyle(columns = 2, fontWeight = "bold")


# Dashboard ---------------------------------------------------------------

# Verificar se pasta "docs" existe no projeto
if(!dir.exists("docs")){dir.create("docs")}

# Renderizar dashboard
rmarkdown::render(
  input       = "./docs/dash_ipca.Rmd",
  output_file = "index.html"
)



