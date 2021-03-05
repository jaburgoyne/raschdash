.digits <- function(x) {
  floor(pmax(0, 2 - log10(diff(range(x, na.rm = TRUE)))))
}

.formatS <- function(x, digits = NULL) {
  digits <- if (is_null(digits)) .digits(x) else digits
  formatC(x, digits = digits, format = "f")
}

.formatP <- function(p) {
  dplyr::case_when(
    p < 0.01 ~ " < 0.01",
    p > 0.99 ~ " > 0.99",
    TRUE ~
    stringr::str_c(
      " = ",
      formatC(p, digits = 2, format = "f")
    )
  )
}

.join_calibration <- function(df, rdfit, par, mean, sd) {
  stan_par <-
    if (vec_in(par, c("group", "person"))) {
      stringr::str_c(par, "_ability")
    } else if (vec_in(par, c("testlet", "item"))) {
      stringr::str_c(par, "_difficulty")
    }
  if (has_name(df, par) && !has_name(df, stan_par)) {
    ## There must be a way to avoid defining an extra symbol
    ## variable, but I cannot figure it out.
    parsym <- sym(par)
    df %>%
      dplyr::left_join(
        dplyr::select(
          summary(rdfit, {{ parsym }}),
          par, stan_par
        ),
        by = par
      ) %>%
      dplyr::mutate(
        dplyr::across(
          stan_par,
          function(x) mean + sd * x
        )
      )
  } else {
    df
  }
}

.plot_bubble <- function(df, mark_limits) {
  df %>%
    plotly::plot_ly(
      x = ~item_difficulty,
      y = ~person_ability
    ) %>%
    plotly::add_markers(
      text =
        ~ str_c(
          person, " (θ = ",
          .formatS(person_ability), ")\n",
          item, " (δ = ",
          .formatS(item_difficulty), ")\n",
          "\n",
          "Observed Score = ",
          observed_score, " / ",
          max_score, "\n",
          "Expected Score = ",
          .formatS(expected_score), " / ",
          max_score, "\n",
          "P(Score ≤ Observed)",
          .formatP(p_score), "\n",
          "\n",
          "Information Content = ",
          .formatS(information_content, 1), "\n",
          "Entropy = ",
          .formatS(entropy, 1), "\n",
          "P(Information ≤ Observed)",
          .formatP(p_information)
        ),
      hoverinfo = "text",
      size = ~ pmax(0.25, information_content),
      sizes = ~ c(0.25, max(information_content)),
      marker = list(sizemode = "area", sizeref = 0.01),
      color = ~observed_mark,
      stroke = ~expected_mark, span = I(2)
    ) %>%
    plotly::colorbar(which = 1, limits = mark_limits, title = "Mark") %>%
    plotly::colorbar(which = 2, limits = mark_limits, title = "Expected") %>%
    ## Who knows why the expected colour scale is listed as trace 4?
    plotly::style(traces = 3, marker.showscale = FALSE) %>%
    plotly::layout(
      autosize = TRUE,
      showlegend = FALSE,
      xaxis =
        list(
          zeroline = FALSE,
          showline = FALSE,
          showgrid = FALSE,
          title = "",
          ticks = ""
        ),
      yaxis =
        list(
          zeroline = FALSE,
          showline = FALSE,
          showgrid = FALSE,
          title = "",
          ticks = ""
        )
    )
}

#' @export
plot.rdfit <- function(x, ..., mark = identity, mean = 0, sd = 1) {
  df <-
    summary(x, ..., mark = mark, mean = mean, sd = sd) %>%
    dplyr::mutate(
      dplyr::across(
        where(~ vec_is(.x, new_eloo())),
        field, "value"
      )
    ) %>%
    .join_calibration(x, "group", mean, sd) %>%
    .join_calibration(x, "person", mean, sd) %>%
    .join_calibration(x, "testlet", mean, sd) %>%
    .join_calibration(x, "item", mean, sd)
  mark_limits <- c(mark(0), mark(1))
  plotfun <-
    if ((has_name(df, "group") || has_name(df, "person"))
    && (has_name(df, "testlet") || has_name(df, "item"))) {
      .plot_bubble
    } else if (has_name(df, "group") || has_name(df, "person")) {
      .plot_horizontal
    } else if (has_name(df, "testlet") || has_name(df, "item")) {
      .plot_vertical
    }
  plotlist <-
    df %>%
    tidyr::nest(dat = !any_of("cohort")) %>%
    deframe() %>%
    map(plotfun, mark_limits)
  if (length(plotlist) == 1) plotlist[[1]] else plotlist
}
