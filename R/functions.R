# LA-ICP-MS FT ages and stats
#
## Constants
pkg.env <- new.env()
pkg.env$e <- exp(1)
pkg.env$lambda <- c(
  1.55125e-10,
  0.000000083
) # U238 decay constant of Jaffey et al. (1971)
pkg.env$g <- 0.5 # geometry factor
pkg.env$lambda.fission <- c(8.45e-17, 0.1e-17) # fission decay constant (Holden and Hoffmann 2000)

## Standards
pkg.env$standards <- data.frame(
  mineral = c("zircon", "zircon", "apatite", "apatite"),
  standard = c(NA, "FCT", "Dur", NA),
  t = c(NA, 28200000, 31440000, NA),
  st = c(NA, NA, NA, NA),
  ref = c(NA, NA, NA, NA)
)
pkg.env$standards$zeta.factor <-
  pkg.env$e^(pkg.env$lambda[1] * pkg.env$standards$t) - 1


#' Area of the laser spot and count area
#' @param x radius
get_area <- function(x) {
  pi * x^2
}


#' @title spontaneous fission-track density
#' @param Ns spontaneous fission track count
#' @param area Area
get_rhos <- function(Ns, area) {
  Ns / area
}

#' Helper function for U times A
#' @param r radius
#' @param U Uranium concentration
get_UxA <- function(r, U) {
  get_area(r) * U
}

#' Single-grain FT age
get_age <-
  function(U,
           rhos,
           zeta,
           lambda = pkg.env$lambda[1],
           g = pkg.env$g) {
    (1 / lambda * log(1 + g * lambda * zeta[1] * rhos / U)) / 1000000
  }

get_age_error <- function(t, Ns, U, sU, zeta) {
  if (is.null(zeta[2])) {
    zeta[2] <- 0
  }
  t * (sqrt((1 / Ns) + (zeta[2] / zeta[1])^2 + (sU / U)^2))
}

get_age_error_perc <- function(t, st) {
  st / t * 100
}

#' Zeta error
get_zeta_error <- function(x, zeta) {
  Ns <- U <- sU <- NULL
  sqrt(1 / x$Ns + (x$st / x$t)^2 + (x$sU / x$U)^2) * zeta
}


#' Single-grain FT age FT ages
#'
#' blabla bla
#'
#' @param Ns  count of spontaneous fission-tracks
#' @param rhos Spontaneous fission track density
#' @param U U-concentration, given as \eqn{^{238}}{238}U/\eqn{^{29}}{29}Ca based on cps (for zircon), \eqn{^{238}}{238}U/\eqn{^{43}}{43}Ca based on cps (apatite), or \eqn{^{238}}{238}U in ppm
#' @param sU standard deviation of \code{U}
#' @param zeta calibration factor based on FCT2 fission-track age standard
#' @note You can either use U ratio cps, or U ppm. Whatever you choose for zeta calculation you must use for age calculation too
get_ages <- function(x, rhos, zeta, ...) {
  t <- get_age(
    U = x$U,
    rhos = rhos,
    zeta = zeta
  )
  st <-
    get_age_error(
      t = t,
      Ns = x$Ns,
      U = x$U,
      sU = x$sU,
      zeta = zeta
    )
  st.perc <- get_age_error_perc(t, st)
  t.min <- t - st
  t.max <- t + st

  data.frame(t, st, st.perc, t.max, t.min, rhos)
}

#' Chi-squared statistic
#' @inheritParams chi_stats
get_chisq <- function(t, st, na.rm = TRUE) {
  z.i <- log(t)
  sigma.i <- st / t

  a <- (z.i / sigma.i)^2
  b <- z.i / sigma.i^2
  c <- 1 / sigma.i^2

  sum(a, na.rm = na.rm) - sum(b, na.rm = na.rm)^2 / sum(c, na.rm = na.rm)
}

#' Probability of chi-squared
#'
#' @param chisq Chi-squared statistic
#' @param n Number of grains
#' @importFrom stats pchisq
get_prob_chisq <- function(chisq, n) {
  stats::pchisq(chisq, df = n - 1, lower.tail = FALSE)
}

#' Chisq stats
#' @param t ages
#' @param st analytical uncertainties
#' @param na.rm logical. wether NA values should be removed (TRUE) or not (FALSE)
chi_stats <- function(t, st, na.rm = TRUE) {
  if (na.rm) {
    x <- data.frame(t, st)
    x <- subset(x, !is.na(t))
    t <- x$t
    st <- x$st
  }
  chisq <- get_chisq(t, st, na.rm = pkg.env$e)
  P.chisq <- get_prob_chisq(chisq, n = length(t))
  return(c(
    "chisq" = chisq,
    "probability" = P.chisq
  ))
}


get_age_pooled0 <-
  function(x,
           r,
           zeta,
           lambda = pkg.env$lambda[1],
           g = pkg.env$g,
           ...) {
    UxA <- get_UxA(r, x$U)
    sums <- sum(x$Ns, na.rm = TRUE) / sum(UxA, na.rm = TRUE)
    age <- 1 / lambda *
      log(1 + (g * lambda * zeta * sums))
    age / 1000000
  }

#' Pooled age
#'
#' @param x dataset
#' @param r Radius
#' @param zeta two-element vector with the zeta-factor and its standard error.
pooled_age <- function(x, r, zeta, ...) {
  t.pool <- get_age_pooled0(x, r = r, zeta = zeta[1])
  st.pool <- t.pool * sqrt(1 / sum(x$Ns, na.rm = TRUE) + (zeta[2] / zeta[1])^
    2 + (sum(x$sU, na.rm = TRUE) / sum(x$U, na.rm = TRUE))^2)

  names(t.pool) <- NULL
  names(st.pool) <- NULL

  c(age = t.pool, std = st.pool)
}

#' Absolute FT ages
#'
#' Calculation of the ages from LA-ICP-MS using zeta calibration
#'
#' @param x \code{data.frame} containing counts of spontaneous tracks,
#' U concentrations, and the uncertainties of U concentrations
#' @param spotsize the laser ablation spot size (diameter in micrometer)
#' @param zeta two-element vector with the zeta-factor and its standard error.
#' @param ... additional arguments
#' @return \code{list} with the single grain ages and their analytical uncertainties,
#' the pooled age, and the chisq statistics
#' @details The fission track age using the zeta calibration approach
#' (Hasebe et al. 2009) is given by
#' \deqn{t = \frac{1}{\lambda} \ln \left(1+\frac{1}{2} \lambda \zeta_{icp} \frac{N)s}{A_s [ ^{238}U/^xX]} \right)}
#' where \eqn{\lambda}{lambda} is the 238U decay constant,
#' \eqn{\zeta_{icp}}{zeta_icp} is the zeta calibration factor,
#' \eqn{N_s}{Ns} is the count of spontaneous fission tracks,
#' \eqn{A_s}{As} is the size of the ablated area, and
#' \eqn{[ ^{238}U/^xX]}{238U} either stands for the
#' \eqn{^{238}}{238}U-concentration (in ppm) or for the U/Ca (for apatite) or
#' U/Si (for zircon) ratio measurement.
#' @references Hasebe, N., Carter, A., Hurford, A.J., Arai, S., 2009. The effect of chemical etching on LA-ICP-MS analysis in determining uranium concentration for fission-track chronometry. Geol. Soc. Lond. Spec. Publ. 324 (1), 37--46. \doi{10.1144/SP324.3}
#'
#' Vermeesch, Pieter, 2017. Statistics for LA-ICP-MS based fission track dating. Chem. Geol. 456, 19--27. \doi{10.1016/j.chemgeo.2017.03.002}
#' @export
#' @examples
#' data("sample")
#' age_ICP(sample, spotsize = 40, zeta = c(0.1188, 0.0119))
age_ICP <- function(x, spotsize = 40, zeta, ...) {
  r <- spotsize / 2 / 10000 # radius in cm
  rhos <- get_rhos(Ns = x$Ns, area = get_area(r))
  res <- get_ages(x, rhos = rhos, zeta = zeta)
  df <- cbind(x, res)
  t.pooled <- pooled_age(df, r = r, zeta = zeta)
  t.stats <- chi_stats(res$t, res$st)
  zeta.err <- get_zeta_error(df, zeta[1])

  list(
    data = x,
    ages = data.frame(res, zeta.error = zeta.err),
    age.pooled = t.pooled,
    stats = t.stats
  )
}
#' Zeta calibration coefficient for LA-ICP-MS fission track dating
#'
#' Calibrate the ICP-session zeta factor using the number of spontaneous tracks,
#' the area over which these were counted and one ore more U/Ca, U/Si or
#' U-concentration measurements (in ppm) and their analytical uncertainties.
#'
#' @param x \code{data.frame} containing counts of spontaneous tracks, U concentrations, and the uncertainties of U concentrations
#' @param spotsize the laser ablation spot size (diameter in micrometer)
#' @param tst (optional) two-column vector giving the true age and its standard error. overrides mineral and standard
#' @param mineral Mineral of the standard, one of 'zircon' or 'apatite'
#' @param standard Name of the standard, one of 'Dur' or 'FCT'
#' @param lambda (optional) decay constant
#' @param g (optional) geometry factor (default: 0.5)
#' @param ... additional arguments
#' @return list containing: two-element vector with the session zeta and its standard error
#' @references Hasebe, N., Carter, A., Hurford, A.J., Arai, S., 2009. The effect of chemical etching on LA-ICP-MS analysis in determining uranium concentration for fission-track chronometry. Geol. Soc. Lond. Spec. Publ. 324 (1), 37--46. \doi{10.1144/SP324.3}
#' @export
#' @examples
#' data("standard")
#' zeta_ICP(standard, spotsize = 40, mineral = "apatite", standard = "Dur")
zeta_ICP <-
  function(x,
           spotsize = 40,
           mineral = c("zircon", "apatite"),
           standard = c("Dur", "FCT"),
           tst,
           lambda = pkg.env$lambda[1],
           g = pkg.env$g,
           ...) {
    min <- match.arg(mineral)
    name <- match.arg(standard)
    zeta.factor <- NULL
    stopifnot(is.numeric(spotsize))

    if (missing(tst)) {
      std <-
        subset(pkg.env$standards,
          mineral == min & standard == name,
          select = c(t, zeta.factor)
        )
    } else {
      stopifnot(is.numeric(tst))
      std <- data.frame(t = tst[1], zeta.factor = tst[2])
    }

    r <- spotsize / 2 / 10000 # radius in cm
    # area <- get_area(r)
    rhos <- get_rhos(Ns = x$Ns, area = get_area(r))

    zeta <- (std$zeta.factor * x$U) / (g * lambda * rhos)


    # final zeta calc
    UxA <- get_UxA(r, x$U)
    sumNS <- sum(x$Ns, na.rm = TRUE)
    sumUxA <- sum(UxA, na.rm = TRUE)
    sumU <- sum(x$U, na.rm = TRUE)
    sumsU <- sum(x$sU, na.rm = TRUE)

    st <- (std$t / 3) / 100

    zeta.icp <- std$zeta.factor * sumUxA / (g * lambda * sumNS)
    szeta.icp <-
      zeta.icp * sqrt((1 / sumNS) + (st / std$t)^2 + (sumsU / sumU)^2)

    res <- data.frame(rhos = rhos, zeta = zeta)

    list(
      data = x,
      results = res,
      zeta = c(zeta = zeta.icp, std = szeta.icp)
    )
  }


#' Export to IsoplotR
#'
#' transform into a \code{"fissiontracks"} object
#'
#' @param x a data.frame of spontaneous and fission track counts and the U-concentration
#' @param spotsize the laser ablation spot size
#' @param zeta the zeta calibration constant extracted from the input data
#' @param ages logical. \code{FALSE} if \code{x} represents unprocessed data
#' (the default), or \code{TRUE} if processed ages
#' @return An object of class \code{fissiontracks}
#' @export
#' @examples
#' data("sample")
#' as.fissiontracks(sample, spotsize = 40, zeta = session.zeta$zeta)
as.fissiontracks <- function(x, spotsize = 40, zeta, ages = FALSE) {
  if (!ages) {
    dfU <- cbind(x$U, NA, NA)
    lU <- split(dfU, 1:NROW(dfU))

    dfsU <- cbind(x$sU, NA, NA)
    lsU <- split(dfsU, 1:NROW(dfsU))

    x2 <- list(
      format = 2,
      zeta = zeta,
      spotSize = spotsize,
      Ns = x$Ns,
      U = lU,
      sU = lsU
    )
    class(x2) <- "fissiontracks"
  } else {
    x2 <- cbind("t" = x$t, "s[t]" = x$st)
  }
  return(x2)
}
