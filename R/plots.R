matchit.covplot <- function(object, type = "qq", interactive = TRUE, which.xs = NULL, ...) {

  #Create covariate matrix; include exact and mahvars
  if (!is.null(which.xs)) {
    if (!is.character(which.xs)) {
      stop("'which.xs' should be a vector of variables for which balance is to be displayed.", call. = FALSE)
    }
    if (!all(which.xs %in% names(object$X))) {
      missing.vars <- setdiff(which.xs, names(object$X))
      stop(paste0("The requested ", ngettext(length(missing.vars), "variable ", "variables "),
                  word_list(missing.vars, is.are = TRUE, quotes = 2),
                  " not present in the supplied matchit object."), call. = FALSE)
    }

    which.xs.f <- terms(reformulate(which.xs))

    X <- get.covs.matrix(which.xs.f, data = object$X)
  }
  else {
    #Create covariate matrix; include exact and mahvars
    X <- get.covs.matrix(object$formula, data = object$X)

    if (!is.null(object$exact)) {
      Xexact <- get.covs.matrix(object$exact, data = object$X)
      X <- cbind(X, Xexact[,setdiff(colnames(Xexact), colnames(X)), drop = FALSE])
    }

    if (!is.null(object$mahvars)) {
      Xmahvars <- get.covs.matrix(object$mahvars, data = object$X)
      X <- cbind(X, Xmahvars[,setdiff(colnames(Xmahvars), colnames(X)), drop = FALSE])
    }
  }

  t <- object$treat

  sw <- if (is.null(object$s.weights)) rep(1, length(t)) else object$s.weights
  w <- object$weights * sw
  if (is.null(w)) w <- rep(1, length(t))

  split(w, t) <- lapply(split(w, t), function(x) x/sum(x))
  split(sw, t) <- lapply(split(sw, t), function(x) x/sum(x))

  varnames <- colnames(X)

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  oma <- c(2.25, 0, 3.75, 1.5)
  if (type == "qq") {
  opar <- par(mfrow = c(3, 3), mar = rep.int(1/2, 4), oma = oma)
  }
  else if (type == "ecdf") {
    opar <- par(mfrow = c(3, 3), mar = c(1.5,.5,1.5,.5), oma = oma)
  }

  for (i in seq_along(varnames)){
    x <- X[,i]

    plot(x, type= "n" , axes=FALSE)
    if (((i-1)%%3)==0) {

      if (type == "qq") {
        htext <- "eQQ Plots"
        mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
        mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
        mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
        mtext("Control Units", 1, 0, TRUE, 2/3, cex=1, font = 1)
        mtext("Treated Units", 4, 0, TRUE, 0.5, cex=1, font = 1)
      }
      else if (type == "ecdf") {
        htext <- "eCDF Plots"
        mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
        mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
        mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
      }

    }

    par(usr = c(0, 1, 0, 1))
    l.wid <- strwidth(varnames, "user")
    cex.labels <- max(0.75, min(1.45, 0.85/max(l.wid)))
    text(0.5, 0.5, varnames[i], cex = cex.labels)

    if (type == "qq") {
      qqplot_match(x = x, t = t, w = w, sw = sw, ...)
    }
    else if (type == "ecdf") {
      ecdfplot_match(x = x, t = t, w = w, sw = sw, ...)
    }

    devAskNewPage(ask = interactive)
  }
  devAskNewPage(ask = FALSE)
}

matchit.covplot.subclass <- function(object, type = "qq", which.subclass = NULL,
                                    interactive = TRUE, which.xs = NULL, ...) {

  #Create covariate matrix; include exact and mahvars
  if (!is.null(which.xs)) {
    if (!is.character(which.xs)) {
      stop("'which.xs' should be a vector of variables for which balance is to be displayed.", call. = FALSE)
    }
    if (!all(which.xs %in% names(object$X))) {
      missing.vars <- setdiff(which.xs, names(object$X))
      stop(paste0("The requested ", ngettext(length(missing.vars), "variable ", "variables "),
                  word_list(missing.vars, is.are = TRUE, quotes = 2),
                  " not present in the supplied matchit object."), call. = FALSE)
    }

    which.xs.f <- terms(reformulate(which.xs))

    X <- get.covs.matrix(which.xs.f, data = object$X)
  }
  else {
    #Create covariate matrix; include exact and mahvars
    X <- get.covs.matrix(object$formula, data = object$X)

    if (!is.null(object$exact)) {
      Xexact <- get.covs.matrix(object$exact, data = object$X)
      X <- cbind(X, Xexact[,setdiff(colnames(Xexact), colnames(X)), drop = FALSE])
    }

    if (!is.null(object$mahvars)) {
      Xmahvars <- get.covs.matrix(object$mahvars, data = object$X)
      X <- cbind(X, Xmahvars[,setdiff(colnames(Xmahvars), colnames(X)), drop = FALSE])
    }
  }

  discrete.cutoff <- 5

  t <- object$treat

  if (!is.atomic(which.subclass)) {
    stop("The argument to 'subclass' must be NULL or the indices of the subclasses for which to display covariate distributions.", call. = FALSE)
  }
  if (!all(which.subclass %in% object$subclass[!is.na(object$subclass)])) {
    stop("The argument supplied to 'subclass' is not the index of any subclass in the matchit object.", call. = FALSE)
  }

  varnames <- colnames(X)

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  oma <- c(2.25, 0, 3.75, 1.5)

  for (s in which.subclass) {
    if (type == "qq") {
      opar <- par(mfrow = c(3, 3), mar = rep.int(1/2, 4), oma = oma)
    }
    else if (type == "ecdf") {
      opar <- par(mfrow = c(3, 3), mar = c(1.5,.5,1.5,.5), oma = oma)
    }

    sw <- if (is.null(object$s.weights)) rep(1, length(t)) else object$s.weights
    w <- sw*(!is.na(object$subclass) & object$subclass == s)

    split(w, t) <- lapply(split(w, t), function(x) x/sum(x))
    split(sw, t) <- lapply(split(sw, t), function(x) x/sum(x))

    for (i in seq_along(varnames)){

      x <- X[,i]

      plot(x, type = "n" , axes = FALSE)

      if (((i-1)%%3)==0) {

        if (type == "qq") {
          htext <- paste0("eQQ Plots (Subclass ", s,")")
          mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
          mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
          mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
          mtext("Control Units", 1, 0, TRUE, 2/3, cex=1, font = 1)
          mtext("Treated Units", 4, 0, TRUE, 0.5, cex=1, font = 1)
        }
        else if (type == "ecdf") {
          htext <- paste0("eCDF Plots (Subclass ", s,")")
          mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
          mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
          mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
        }

      }

      #Empty plot with variable name
      par(usr = c(0, 1, 0, 1))
      l.wid <- strwidth(varnames, "user")
      cex.labels <- max(0.75, min(1.45, 0.85/max(l.wid)))
      text(0.5, 0.5, varnames[i], cex = cex.labels)

      if (type == "qq") {
        qqplot_match(x = x, t = t, w = w, sw = sw, ...)
      }
      else if (type == "ecdf") {
        ecdfplot_match(x = x, t = t, w = w, sw = sw, ...)
      }

      devAskNewPage(ask = interactive)
    }
  }
  devAskNewPage(ask = FALSE)
}

qqplot_match <- function(x, t, w, sw, discrete.cutoff = 5, ...) {

  ord <- order(x)
  x_ord <- x[ord]
  t_ord <- t[ord]

  u <- unique(x_ord)

  #Need to interpolate larger group to be same size as smaller group

  #Unmatched sample
  sw_ord <- sw[ord]
  sw1 <- sw_ord[t_ord == 1]
  sw0 <- sw_ord[t_ord != 1]

  x1 <- x_ord[t_ord == 1][sw1 > 0]
  x0 <- x_ord[t_ord != 1][sw0 > 0]

  swn1 <- sum(sw[t==1] > 0)
  swn0 <- sum(sw[t==0] > 0)

  if (swn1 < swn0) {
    if (length(u) <= discrete.cutoff) {
      x0probs <- vapply(u, function(u_) weighted.mean(x0 == u_, sw0[sw0 > 0]), numeric(1L))
      x0cumprobs <- c(0, cumsum(x0probs)[-length(u)], 1)
      x0 <- u[findInterval(cumsum(sw1[sw1 > 0]), x0cumprobs, rightmost.closed = TRUE)]
    }
    else {
      x0 <- approx(cumsum(sw0[sw0 > 0]), y = x0, xout = cumsum(sw1[sw1 > 0]), rule = 2,
                   method = "constant", ties = "ordered")$y
    }
  }
  else {
    if (length(u) <= discrete.cutoff) {
      x1probs <- vapply(u, function(u_) weighted.mean(x1 == u_, sw1[sw1 > 0]), numeric(1L))
      x1cumprobs <- c(0, cumsum(x1probs)[-length(u)], 1)
      x1 <- u[findInterval(cumsum(sw0[sw0 > 0]), x1cumprobs, rightmost.closed = TRUE)]
    }
    else {
      x1 <- approx(cumsum(sw1[sw1 > 0]), y = x1, xout = cumsum(sw0[sw0 > 0]), rule = 2,
                   method = "constant", ties = "ordered")$y
    }
  }

  if (length(u) <= discrete.cutoff) {
    md <- min(diff(u))
    x0 <- jitter(x0, amount = .1*md)
    x1 <- jitter(x1, amount = .1*md)
  }

  rr <- range(c(x0, x1))
  plot(x0, x1, xlab = "", ylab = "", xlim = rr, ylim = rr, axes = FALSE, ...)
  abline(a = 0, b = 1)
  abline(a = (rr[2]-rr[1])*0.1, b = 1, lty = 2)
  abline(a = -(rr[2]-rr[1])*0.1, b = 1, lty = 2)
  axis(2)
  box()

  #Matched sample
  w_ord <- w[ord]
  w1 <- w_ord[t_ord == 1]
  w0 <- w_ord[t_ord != 1]

  x1 <- x_ord[t_ord == 1][w1 > 0]
  x0 <- x_ord[t_ord != 1][w0 > 0]

  wn1 <- sum(w[t==1] > 0)
  wn0 <- sum(w[t==0] > 0)

  if (wn1 < wn0) {
    if (length(u) <= discrete.cutoff) {
      x0probs <- vapply(u, function(u_) weighted.mean(x0 == u_, w0[w0 > 0]), numeric(1L))
      x0cumprobs <- c(0, cumsum(x0probs)[-length(u)], 1)
      x0 <- u[findInterval(cumsum(w1[w1 > 0]), x0cumprobs, rightmost.closed = TRUE)]
    }
    else {
      x0 <- approx(cumsum(w0[w0 > 0]), y = x0, xout = cumsum(w1[w1 > 0]), rule = 2,
                   method = "constant", ties = "ordered")$y
    }
  }
  else {
    if (length(u) <= discrete.cutoff) {
      x1probs <- vapply(u, function(u_) weighted.mean(x1 == u_, w1[w1 > 0]), numeric(1L))
      x1cumprobs <- c(0, cumsum(x1probs)[-length(u)], 1)
      x1 <- u[findInterval(cumsum(w0[w0 > 0]), x1cumprobs, rightmost.closed = TRUE)]
    }
    else {
      x1 <- approx(cumsum(w1[w1 > 0]), y = x1, xout = cumsum(w0[w0 > 0]), rule = 2,
                   method = "constant", ties = "ordered")$y
    }
  }

  if (length(u) <= discrete.cutoff) {
    md <- min(diff(u))
    x0 <- jitter(x0, amount = .1*md)
    x1 <- jitter(x1, amount = .1*md)
  }

  plot(x0, x1, xlab = "", ylab = "", xlim = rr, ylim = rr, axes = FALSE, ...)
  abline(a = 0, b = 1)
  abline(a = (rr[2]-rr[1])*0.1, b = 1, lty = 2)
  abline(a = -(rr[2]-rr[1])*0.1, b = 1, lty = 2)
  box()
}

ecdfplot_match <- function(x, t, w, sw, ...) {
  ord <- order(x)
  x.min <- x[ord][1]
  x.max <- x[ord][length(x)]
  x.range <- x.max - x.min

  #Unmatched samples
  plot(x = x, y = w, type= "n" , xlim = c(x.min - .02 * x.range, x.max + .02 * x.range),
       ylim = c(0, 1), axes = TRUE, ...)

  for (tr in 0:1) {
    in.tr <- t[ord] == tr
    ordt <- ord[in.tr]
    cswt <- c(0, cumsum(sw[ordt]), 1)
    xt <- c(x.min - .02 * x.range, x[ordt], x.max + .02 * x.range)

    lines(x = xt, y = cswt, type = "s", col = if (tr == 0) "grey60" else "black")

  }

  abline(h = 0:1)
  box()

  #Matched sample
  plot(x = x, y = w, type= "n" , xlim = c(x.min - .02 * x.range, x.max + .02 * x.range),
       ylim = c(0, 1), axes = FALSE, ...)
  for (tr in 0:1) {
    in.tr <- t[ord] == tr
    ordt <- ord[in.tr]
    cwt <- c(0, cumsum(w[ordt]), 1)
    xt <- c(x.min - .02 * x.range, x[ordt], x.max + .02 * x.range)

    lines(x = xt, y = cwt, type = "s", col = if (tr == 0) "grey60" else "black")

  }

  abline(h = 0:1)
  axis(1)
  box()
}

