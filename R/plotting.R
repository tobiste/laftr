#' Adaptive kernel density estimates
#'
#' Modification of [ggplot2::stat_density()] for kernel density estimates using a combination of the
#' Botev (2010) bandwidth selector and the Abramson (1982) adaptive kernel
#' bandwidth modifier.
#'
#' @inheritParams ggplot2::stat_density
#' @inheritParams ggplot2::geom_density
#' @param from,to the left and right-most points of the grid at which the density is to be estimated
#' @param bw the bandwidth of the KDE. If `NULL`, `bw` will be calculated automatically using the algorithm by Botev et al. (2010).
#' @param adaptive logical flag controlling if the adaptive KDE modifier of Abramson (1982) is used
#' @import ggplot2
#' @import IsoplotR
#' @name aKDE
#' @importFrom provenance botev
#' @source Algorithm for adaptive kernel is modified from [IsoplotR]. The
#' algorithm for the optimal kernel bandwidth is from [provenance::botev()].
#' @examples
#' data("sample")
#' example <- age_ICP(sample, zeta = c(0.1188, 0.0119))
#' # IsoplotR::kde(example$ages$t)
#' ggplot2::ggplot(data = example$ages, mapping = ggplot2::aes(x = t)) +
#'   stat_aKDE(adaptive = TRUE) +
#'   stat_aKDE(adaptive = FALSE, color = "red", fill = NA)
#'
#' ggplot2::ggplot(
#'   data = example$ages,
#'   mapping = ggplot2::aes(x = t, weight = t / st)
#' ) +
#'   geom_aKDE(
#'     ggplot2::aes(y = ggplot2::after_stat(scaled)),
#'     kernel = "epanechnikov", fill = "steelblue", alpha = .75
#'   ) +
#'   ggplot2::geom_histogram(
#'     ggplot2::aes(y = ggplot2::after_stat(ncount)),
#'     color = "grey", fill = "grey", alpha = .5
#'   ) +
#'   ggplot2::geom_rug(alpha = 0.5)
NULL

#' @rdname aKDE
#' @export
stat_aKDE <- function(mapping = NULL, data = NULL, stat = "DensityAdaptive", geom = "area",
                      position = "stack",
                      ...,
                      from = NA, to = NA, bw = NA, adjust = 1, kernel = "gaussian", n = 512,
                      adaptive = TRUE,
                      na.rm = FALSE, bounds = c(-Inf, Inf), show.legend = NA, orientation = NA,
                      inherit.aes = TRUE) {
  layer(
    stat = StatDensityAdaptive, data = data, mapping = mapping, geom = geom, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      from = from, to = to, bw = bw, adaptive = adaptive,
      adjust = adjust, kernel = kernel, n = n, na.rm = na.rm, bounds = bounds, orientation = orientation, ...
    )
  )
}

#' @rdname aKDE
#' @export
geom_aKDE <- function(mapping = NULL, data = NULL, stat = "DensityAdaptive", position = "identity",
                      ..., na.rm = FALSE, orientation = NA, show.legend = NA, inherit.aes = TRUE,
                      outline.type = "upper") {
  outline.type <- rlang::arg_match0(outline.type, c(
    "both", "upper",
    "lower", "full"
  ))
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomDensity,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm, orientation = orientation,
      outline.type = outline.type, ...
    )
  )
}



compute_density2 <- function(x, w = NULL, from, to, bw = NA, adaptive = TRUE, adjust = 1, kernel = "gaussian",
                             n = 512, bounds = c(-Inf, Inf)) {
  # nx <- length(x)
  dens <- getkde(x, from = from, to = to, bw = bw, adaptive = adaptive, n = n, adjust = adjust, kernel = kernel, weights = w)

  if (any(is.finite(bounds))) {
    dens <- ggplot2:::reflect_density(
      dens = dens, bounds = bounds,
      from = from, to = to
    )
  }

  dx <- dens$x
  dy <- dens$y

  nx <- length(stats::na.omit(x))
  maxdy <- max(dy, na.rm = TRUE)

  ggplot2:::data_frame0(
    x = dx, y = dy, density = dy, scaled = dy / maxdy,
    ndensity = dy / maxdy,
    count = dy * nx, n = nx, .size = length(dx)
  )
}

StatDensityAdaptive <- ggplot2::ggproto(
  "_class" = "StatDensityAdaptive", "_inherit" = Stat,
  compute_group = function(data, scales, from = NA, to = NA, bw = NA, adaptive = TRUE, adjust = 1, kernel = "gaussian",
                           n = 512, na.rm = FALSE, bounds = c(-Inf, Inf),
                           flipped_aes = FALSE) {
    data <- ggplot2::flip_data(data, flipped_aes)

    density <- compute_density2(data$x,
      w = data$weight,
      from = from, to = to, bw = bw, adaptive = adaptive, adjust = adjust, kernel = kernel,
      n = n, bounds = bounds
    )
    density$flipped_aes <- flipped_aes
    ggplot2::flip_data(density, flipped_aes)
  },
  required_aes = c("x"),
  dropped_aes = c("weight")
)

#' Adaptive kernel density estimates
#'
#' Modification of [ggplot2::stat_density()] for kernel density estimates using a combination of the
#' Botev (2010) bandwidth selector and the Abramson (1982) adaptive kernel
#' bandwidth modifier.
#' @param x numeric vector. Ages
#' @param ... optional arguments to be passed on to R's density function.
#' @param from numeric. minimum age of the time axis. If NULL, this is set automatically
#' @param to numeric. maximum age of the time axis. If NULL, this is set automatically
#' @param bw numeric. the bandwidth of the KDE. If NULL, bw will be calculated
#' automatically using the algorithm by Botev et al. (2010).
#' @param adaptive logical flag controlling if the adaptive KDE modifier of Abramson (1982) is used
#' @param log	logical. transform the ages to a log scale if TRUE
#' @param  n integer. horizontal resolution (i.e., the number of segments) of the density estimate.
#' @param weights numeric vector of non-negative observation weights, hence of
#' same length as `x`. The default `NULL` is equivalent to
#' `weights = rep(1/nx, nx)` where `nx` is the length of (the finite entries of)
#' `x[]`. If `na.rm = TRUE` and there are `NA`'s in `x`, they and the
#' corresponding weights are removed before computations.
#' In that case, when the original weights have summed to one, they are
#' re-scaled to keep doing so. Will be ignored is `adaptive` is `TRUE`.
getkde <- function(x, from = NA, to = NA, bw = NA, adaptive = TRUE, log = FALSE,
                   n = 512, weights = NULL, ...) {
  out <- list()
  class(out) <- "KDE"
  out$name <- deparse(substitute(x))
  out$log <- log

  # weights <- weights[which(!is.na(x))]
  # x <- stats::na.omit(x)

  if (log) {
    d <- log(x)
  } else {
    d <- x
  }

  if (is.na(bw)) {
    bw <- provenance::botev(d)
  }
  if (is.na(from) | is.na(to)) {
    mM <- IsoplotR:::getmM(x, from, to, log)
    to <- mM$M + bw
    from <- mM$m
    if (mM$m > bw) {
      from <- from - bw
    }
  }
  if (log) {
    from <- log(from)
    to <- log(to)
  }



  out$x <- seq(from = from, to = to, length.out = n)
  if (adaptive) {
    out$y <- IsoplotR:::Abramson(d,
      from = from, to = to, bw = bw, n = n, weights = NULL,
      ...
    )
  } else {
    if (!is.null(weights)) {
      weights <- weights / sum(weights)
    }
    out$y <- stats::density(d, bw,
      from = from, to = to,
      n = n, weights = weights, ...
    )$y
  }
  if (log) {
    out$x <- exp(out$x)
  }
  out$y <- out$y / (sum(out$y) * (to - from) / n)
  out$x <- c(out$x[1], out$x, out$x[n])
  out$y <- c(0, out$y, 0)
  out$bw <- bw
  out$ages <- x
  out
}


# Abramson_weighted <- function(dat, from, to, bw, n = 512, weights = NULL, ...) {
#   nn <- length(dat)
#   if (is.null(weights)) {
#     weights <- rep.int(1, nn)
#   } else {
#     weights / nn
#   }
#   d <- stats::na.omit(dat)
#   nn <- length(d)
#
#   w <- weights[which(!is.na(dat))]
#
#   pdens <- IsoplotR:::pilotdensity(d, bw)
#   G <- IsoplotR:::getG(pdens)
#   lambda <- 0
#   dens <- rep(0, n)
#   for (i in 1:nn) {
#     lambda <- sqrt(G / pdens[i])
#     suppressWarnings(
#       dens <- dens + stats::density(d[i],
#         bw = bw * lambda,
#         from = from, to = to, n = n, weights = w[i], ...
#       )$y
#     )
#   }
#   dens
# }
