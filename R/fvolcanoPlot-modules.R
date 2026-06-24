#' Reusable volcano plot server module
#'
#' Renders a prepared volcano table and returns brushed feature keys. Data
#' preparation and downstream interpretation stay with the caller.
#'
#' @param id Module id.
#' @param data Reactive data frame containing the x-axis, y-metric, and key
#'   columns.
#' @param ... Reserved for future module options.
#' @param x Column name for the x-axis value.
#' @param y Named character vector of y-axis metrics.
#' @param key Column name for the selected feature key.
#' @param label_priority Ordered label columns to use before falling back to
#'   `key`.
#' @param hover Columns to include in Plotly hover text.
#' @param event_source Plotly event source.
#' @param active Reactive or logical value controlling whether the plot renders.
#' @param .selected_event Internal reactive override for tests.
#' @param label_limit Maximum number of labeled points.
#' @param width,height,webgl Passed to [FacileViz::fscatterplot()].
#' @param debug Reserved for future diagnostics.
#' @return A list with reactive `selected`, `selected_keys`, `plot_data`, and
#'   `label_choices` members.
#'
#' @export
fvolcanoPlotServer <- function(
  id,
  data,
  ...,
  x = "logFC",
  y = c("p-value" = "pval", "FDR" = "FDR"),
  key = "feature_id",
  label_priority = c("symbol", "feature_id"),
  hover = c("symbol", "logFC", "FDR"),
  event_source = NULL,
  active = shiny::reactive(TRUE),
  .selected_event = NULL,
  label_limit = 20L,
  width = NULL,
  height = NULL,
  webgl = TRUE,
  debug = FALSE
) {
  if (!shiny::is.reactive(data)) {
    stop("`data` must be a reactive expression", call. = FALSE)
  }
  if (!shiny::is.reactive(active)) {
    active <- shiny::reactive(isTRUE(active))
  }
  y <- .fvolcano_y_choices(y)
  checkmate::assert_string(x)
  checkmate::assert_string(key)
  checkmate::assert_character(label_priority, min.len = 1L)
  checkmate::assert_count(label_limit, positive = TRUE)

  shiny::moduleServer(id, function(input, output, session) {
    if (is.null(event_source)) {
      event_source <- session$ns("selection")
    }
    selected_event <- .selected_event
    if (is.null(selected_event)) {
      selected_event <- shiny::reactive({
        plotly::event_data("plotly_selected", source = event_source)
      })
    } else if (!shiny::is.reactive(selected_event)) {
      stop("`.selected_event` must be a reactive expression", call. = FALSE)
    }

    selected_keys <- shiny::reactiveVal(character())

    plot_data <- shiny::reactive({
      dat. <- shiny::req(data())
      .fvolcano_require_columns(dat., c(x, unname(y), key))
      dat.
    })

    y_metric <- shiny::reactive({
      .fvolcano_y_metric(input$y_metric, y)
    })

    label_choices <- shiny::reactive({
      .fvolcano_label_choices(plot_data(), label_priority, key)
    })

    selected <- shiny::reactive({
      selected_keys()
    })

    shiny::observeEvent(
      plot_data(),
      {
        dat. <- shiny::req(plot_data())
        keep <- input$labels
        if (is.null(keep)) {
          keep <- character()
        }
        keep <- intersect(keep, as.character(dat.[[key]]))
        keep <- utils::head(keep, label_limit)
        current <- selected_keys()
        current <- intersect(current, as.character(dat.[[key]]))
        if (!setequal(current, selected_keys())) {
          selected_keys(current)
        }
        shiny::updateSelectizeInput(
          session,
          "labels",
          choices = label_choices(),
          selected = keep,
          server = TRUE
        )
      },
      ignoreInit = FALSE
    )

    shiny::observeEvent(
      input$labels,
      {
        keep <- input$labels
        if (length(keep) <= label_limit) {
          return()
        }
        keep <- keep[seq_len(label_limit)]
        shiny::updateSelectizeInput(
          session,
          "labels",
          selected = keep,
          server = TRUE
        )
        shinyWidgets::sendSweetAlert(
          session,
          title = "Too many volcano labels",
          text = sprintf(
            "Only the first %d selected targets will be labeled.",
            label_limit
          ),
          type = "warning"
        )
      },
      ignoreNULL = FALSE
    )

    shiny::observeEvent(
      active(),
      {
        if (!isTRUE(active()) && length(selected_keys())) {
          selected_keys(character())
        }
      },
      ignoreInit = FALSE
    )

    output$plot <- plotly::renderPlotly({
      shiny::req(isTRUE(active()))
      dat. <- shiny::req(plot_data())
      metric <- shiny::isolate(y_metric())
      dat.$yaxis <- .fvolcano_y_values(dat., metric)
      label.selected <- shiny::isolate(input$labels)
      label.dat <- .fvolcano_label_data(
        dat., dat.$yaxis, label.selected, x, key, label_priority, label_limit
      )
      point.style <- .fvolcano_point_style(
        dat.,
        dat.$yaxis,
        shiny::isolate(input$x_cutoff),
        shiny::isolate(input$y_cutoff),
        label.selected,
        x,
        key,
        label_limit
      )
      axis.ranges <- .fvolcano_axis_ranges(dat., x)
      hover. <- intersect(hover, colnames(dat.))

      fplot <- FacileViz::fscatterplot(
        dat.,
        c(x, "yaxis"),
        xlabel = x,
        ylabel = .fvolcano_y_label(metric),
        width = width,
        height = height,
        hover = hover.,
        webgl = webgl,
        event_source = event_source,
        key = key
      )
      plt <- plot(fplot)
      plt <- .fvolcano_style_points(plt, point.style)
      plotly::layout(
        plt,
        annotations = .fvolcano_label_annotations(label.dat, dat., x),
        xaxis = list(range = axis.ranges$x),
        yaxis = list(range = axis.ranges$y)
      )
    })

    shiny::observeEvent(
      {
        list(
          input$x_cutoff,
          input$y_cutoff,
          input$y_metric,
          input$labels,
          active(),
          plot_data()
        )
      },
      {
        shiny::req(isTRUE(active()))
        dat. <- shiny::req(plot_data())
        metric <- y_metric()
        y.values <- .fvolcano_y_values(dat., metric)
        dat.$yaxis <- y.values
        point.style <- .fvolcano_point_style(
          dat.,
          y.values,
          input$x_cutoff,
          input$y_cutoff,
          input$labels,
          x,
          key,
          label_limit
        )
        axis.ranges <- .fvolcano_axis_ranges(dat., x)
        label.dat <- .fvolcano_label_data(
          dat., y.values, input$labels, x, key, label_priority, label_limit
        )
        proxy <- plotly::plotlyProxy("plot", session)
        plotly::plotlyProxyInvoke(
          proxy,
          "restyle",
          list(
            y = list(y.values),
            `marker.color` = list(point.style$color),
            `marker.line.color` = list(point.style$line.color),
            `marker.line.width` = list(point.style$line.width)
          ),
          0
        )
        plotly::plotlyProxyInvoke(
          proxy,
          "relayout",
          list(
            "yaxis.title.text" = .fvolcano_y_label(metric),
            "xaxis.range" = axis.ranges$x,
            "yaxis.range" = axis.ranges$y,
            annotations = .fvolcano_label_annotations(label.dat, dat., x)
          )
        )
      },
      ignoreInit = FALSE
    )

    shiny::observeEvent(
      selected_event(),
      {
        shiny::req(isTRUE(active()))
        next_keys <- .fvolcano_selected_keys(selected_event())
        if (!setequal(next_keys, selected_keys())) {
          selected_keys(next_keys)
        }
      },
      ignoreNULL = FALSE,
      ignoreInit = TRUE
    )

    vals <- list(
      selected = selected,
      selected_keys = selected,
      plot_data = plot_data,
      label_choices = label_choices,
      .ns = session$ns
    )
    class(vals) <- c("FVolcanoPlot", "ReactiveFacileViz")
    vals
  })
}

#' Reusable volcano plot UI module
#'
#' Provides y-axis metric, label, cutoff, and Plotly output controls for
#' [fvolcanoPlotServer()].
#'
#' @param id Module id.
#' @param ... Reserved for future module options.
#' @param y Named character vector of y-axis metrics.
#' @param label_limit Maximum number of labeled points.
#' @param debug Reserved for future diagnostics.
#' @return A Shiny tag list.
#'
#' @export
fvolcanoPlotUI <- function(
  id,
  ...,
  y = c("p-value" = "pval", "FDR" = "FDR"),
  label_limit = 20L,
  debug = FALSE
) {
  ns <- shiny::NS(id)
  y <- .fvolcano_y_choices(y)
  label.selectize.options <- .assay_feature_selectize_options(list(
    maxItems = label_limit,
    placeholder = sprintf("Paste or select up to %d targets", label_limit)
  ))
  controls.style <- paste(
    "border: 1px solid #b8c7ce;",
    "border-radius: 4px;",
    "padding: 10px 12px 0;",
    "margin: 10px 0;",
    "background-color: #f9fafb;"
  )
  controls.subtitle.style <- paste(
    "font-weight: normal;",
    "margin-left: 0.5em;",
    "color: #666;"
  )
  feedback.style <- paste(
    "display:none;",
    "margin-top: 0.35em;",
    "font-size: 0.9em;"
  )

  shiny::tagList(
    shiny::tags$div(
      style = controls.style,
      shiny::tags$div(
        style = "font-weight: bold; margin-bottom: 8px; color: #444;",
        shiny::tags$span("Volcano plot controls"),
        shiny::tags$span(
          style = controls.subtitle.style,
          "for labels and Blue / Red coloring"
        )
      ),
      shiny::fluidRow(
        shiny::column(
          width = 4,
          shinyWidgets::radioGroupButtons(
            ns("y_metric"),
            "Y-axis",
            choices = y,
            selected = unname(y)[1L],
            size = "xs",
            justified = TRUE
          )
        ),
        shiny::column(
          width = 8,
          shiny::selectizeInput(
            ns("labels"),
            "Labels",
            choices = NULL,
            multiple = TRUE,
            options = label.selectize.options
          ),
          shiny::tags$div(
            id = ns("labels_feedback"),
            class = "text-muted",
            role = "status",
            style = feedback.style
          )
        )
      ),
      shiny::fluidRow(
        shiny::column(
          width = 6,
          shiny::numericInput(
            ns("x_cutoff"),
            "X cutoff",
            value = NULL,
            min = 0,
            step = 0.1
          )
        ),
        shiny::column(
          width = 6,
          shiny::numericInput(
            ns("y_cutoff"),
            "Y cutoff",
            value = NULL,
            min = 0,
            step = 0.1
          )
        )
      )
    ),
    shinycssloaders::withSpinner(plotly::plotlyOutput(ns("plot")))
  )
}

.fvolcano_y_choices <- function(y) {
  checkmate::assert_character(y, min.len = 1L, any.missing = FALSE)
  if (is.null(names(y))) {
    names(y) <- y
  }
  y
}

.fvolcano_y_metric <- function(metric, y) {
  y <- .fvolcano_y_choices(y)
  if (length(metric) != 1L || !metric %in% unname(y)) {
    metric <- unname(y)[1L]
  }
  metric
}

.fvolcano_y_label <- function(metric) {
  sprintf("-log10(%s)", metric)
}

.fvolcano_y_values <- function(dat., metric) {
  -log10(dat.[[metric]])
}

.fvolcano_feature_names <- function(dat., label_priority, key) {
  labels <- rep(NA_character_, nrow(dat.))
  for (col in label_priority) {
    if (!col %in% colnames(dat.)) {
      next
    }
    vals <- as.character(dat.[[col]])
    missing <- is.na(labels) | !nzchar(labels)
    labels[missing] <- vals[missing]
  }
  missing <- is.na(labels) | !nzchar(labels)
  labels[missing] <- as.character(dat.[[key]][missing])
  labels
}

.fvolcano_label_choices <- function(dat., label_priority, key) {
  labels <- .fvolcano_feature_names(dat., label_priority, key)
  duplicate.labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicate.labels] <- paste0(
    labels[duplicate.labels],
    " (",
    dat.[[key]][duplicate.labels],
    ")"
  )
  stats::setNames(as.character(dat.[[key]]), labels)
}

.fvolcano_label_data <- function(
  dat., y.values, selected, x, key, label_priority, label_limit
) {
  if (is.null(selected)) {
    selected <- character()
  }
  selected <- utils::head(selected, label_limit)
  idx <- match(selected, as.character(dat.[[key]]), nomatch = 0L)
  idx <- idx[idx > 0L]
  out <- dat.[idx, , drop = FALSE]
  out$yaxis <- y.values[idx]
  out <- out[is.finite(out[[x]]) & is.finite(out$yaxis), , drop = FALSE]
  out$volcano_label <- .fvolcano_feature_names(out, label_priority, key)
  out
}

.fvolcano_axis_ranges <- function(dat., x) {
  padded_range <- function(vals, pad) {
    rng <- suppressWarnings(range(vals, finite = TRUE))
    if (!all(is.finite(rng))) {
      rng <- c(-0.5, 0.5)
    } else if (diff(rng) == 0) {
      rng <- rng + c(-0.5, 0.5)
    }
    span <- diff(rng)
    rng + c(-span * pad[1L], span * pad[2L])
  }
  list(
    x = padded_range(dat.[[x]], c(0.08, 0.18)),
    y = padded_range(dat.$yaxis, c(0.02, 0.14))
  )
}

.fvolcano_label_positions <- function(label.dat, plot.dat, x) {
  if (nrow(label.dat) == 0L) {
    return(label.dat)
  }

  plot.dat <- plot.dat[
    is.finite(plot.dat[[x]]) & is.finite(plot.dat$yaxis),
    ,
    drop = FALSE
  ]
  if (nrow(plot.dat) == 0L) {
    label.dat$label_x <- label.dat[[x]]
    label.dat$label_y <- label.dat$yaxis
    return(label.dat)
  }

  axis.ranges <- .fvolcano_axis_ranges(plot.dat, x)
  xlim <- axis.ranges$x
  ylim <- axis.ranges$y
  grid::pushViewport(grid::viewport(
    width = grid::unit(3.5, "in"),
    height = grid::unit(4, "in"),
    xscale = xlim,
    yscale = ylim
  ))
  on.exit(grid::popViewport(), add = TRUE)

  label.width <- grid::convertWidth(
    grid::stringWidth(label.dat$volcano_label) + grid::unit(14, "pt"),
    "native",
    TRUE
  )
  label.height <- grid::convertHeight(
    grid::stringHeight(label.dat$volcano_label) + grid::unit(11, "pt"),
    "native",
    TRUE
  )
  point.size <- grid::convertWidth(grid::unit(3, "pt"), "native", TRUE)
  point.padding.x <- grid::convertWidth(grid::unit(5, "pt"), "native", TRUE)
  point.padding.y <- grid::convertHeight(grid::unit(5, "pt"), "native", TRUE)

  label.points <- as.matrix(label.dat[, c(x, "yaxis")])
  repel <- utils::getFromNamespace("repel_boxes2", "ggrepel")(
    data_points = label.points,
    point_size = rep(point.size, nrow(label.dat)),
    point_padding_x = point.padding.x,
    point_padding_y = point.padding.y,
    boxes = cbind(
      x1 = label.dat[[x]] - label.width / 2,
      y1 = label.dat$yaxis - label.height / 2,
      x2 = label.dat[[x]] + label.width / 2,
      y2 = label.dat$yaxis + label.height / 2
    ),
    xlim = xlim,
    ylim = ylim,
    hjust = rep(0.5, nrow(label.dat)),
    vjust = rep(0.5, nrow(label.dat)),
    force_push = 1e-5,
    force_pull = 5e-4,
    max_time = 1.5,
    max_iter = 50000L,
    max_overlaps = Inf,
    direction = "both",
    verbose = 0L
  )

  label.dat$label_x <- ifelse(is.finite(repel$x), repel$x, label.dat[[x]])
  label.dat$label_y <- ifelse(is.finite(repel$y), repel$y, label.dat$yaxis)
  label.dat
}

.fvolcano_label_annotations <- function(label.dat, plot.dat, x) {
  if (nrow(label.dat) == 0L) {
    return(list())
  }

  label.dat <- .fvolcano_label_positions(label.dat, plot.dat, x)
  lapply(seq_len(nrow(label.dat)), function(idx) {
    list(
      x = label.dat[[x]][idx],
      y = label.dat$yaxis[idx],
      xref = "x",
      yref = "y",
      text = htmltools::htmlEscape(label.dat$volcano_label[idx]),
      showarrow = TRUE,
      ax = label.dat$label_x[idx],
      ay = label.dat$label_y[idx],
      axref = "x",
      ayref = "y",
      xanchor = "center",
      yanchor = "middle",
      arrowhead = 0,
      arrowwidth = 1,
      arrowcolor = "#555555",
      bgcolor = "rgba(238, 238, 238, 0.86)",
      bordercolor = "#9a9a9a",
      borderwidth = 1,
      borderpad = 1,
      font = list(size = 10, color = "#222222"),
      opacity = 0.98,
      captureevents = FALSE
    )
  })
}

.fvolcano_point_palette <- list(
  background = "rgba(105, 105, 105, 0.18)",
  left = "rgba(0, 82, 204, 0.62)",
  right = "rgba(220, 35, 35, 0.62)",
  left_selected = "rgb(0, 82, 204)",
  right_selected = "rgb(220, 35, 35)",
  selected = "rgb(0, 0, 0)",
  outline = "rgba(35, 35, 35, 0.7)",
  no_outline = "rgba(0, 0, 0, 0)"
)

.fvolcano_point_style <- function(
  dat., y.values, x.cutoff, y.cutoff, selected = NULL, x, key, label_limit
) {
  colors <- rep(.fvolcano_point_palette$background, nrow(dat.))
  line.colors <- rep(.fvolcano_point_palette$no_outline, nrow(dat.))
  line.widths <- rep(0, nrow(dat.))
  left.idx <- right.idx <- rep(FALSE, nrow(dat.))
  x.cutoff <- suppressWarnings(as.numeric(x.cutoff))
  y.cutoff <- suppressWarnings(as.numeric(y.cutoff))
  has.cutoffs <- length(x.cutoff) == 1L && length(y.cutoff) == 1L &&
    is.finite(x.cutoff) && is.finite(y.cutoff)
  if (has.cutoffs) {
    x.cutoff <- abs(x.cutoff)
    left.idx <- dat.[[x]] <= -x.cutoff & y.values >= y.cutoff
    right.idx <- dat.[[x]] >= x.cutoff & y.values >= y.cutoff
    colors[left.idx] <- .fvolcano_point_palette$left
    colors[right.idx] <- .fvolcano_point_palette$right
  }
  selected.idx <- rep(FALSE, nrow(dat.))
  if (length(selected)) {
    selected <- utils::head(selected, label_limit)
    selected.idx <- as.character(dat.[[key]]) %in% selected
    colors[selected.idx & left.idx] <- .fvolcano_point_palette$left_selected
    colors[selected.idx & right.idx] <- .fvolcano_point_palette$right_selected
    colors[selected.idx & !left.idx & !right.idx] <-
      .fvolcano_point_palette$selected
  }
  outline.idx <- left.idx | right.idx | selected.idx
  line.colors[outline.idx] <- .fvolcano_point_palette$outline
  line.widths[outline.idx] <- 0.5
  list(color = colors, line.color = line.colors, line.width = line.widths)
}

.fvolcano_style_points <- function(plt, point.style) {
  if (length(plt$x$data) == 0L) {
    return(plt)
  }
  marker <- plt$x$data[[1L]]$marker
  if (is.null(marker)) {
    marker <- list()
  }
  marker$color <- point.style$color
  marker$size <- 4
  marker$opacity <- 1
  marker$line <- list(
    color = point.style$line.color,
    width = point.style$line.width
  )
  plt$x$data[[1L]]$marker <- marker
  plt
}

.fvolcano_selected_keys <- function(selected) {
  if (is.null(selected)) {
    return(character())
  }
  selected <- selected$key
  selected <- selected[!is.na(selected) & nzchar(selected)]
  as.character(selected)
}

.fvolcano_require_columns <- function(dat., cols) {
  missing <- setdiff(cols, colnames(dat.))
  if (length(missing)) {
    stop(
      sprintf(
        "Volcano data is missing required columns: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(dat.)
}
