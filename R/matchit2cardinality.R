#' Cardinality Matching
#' @name method_cardinality
#' @aliases method_cardinality
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "cardinality"` performs cardinality
#' matching and other forms of matching that use mixed integer programming.
#' Rather than forming pairs, cardinality matching selects the largest subset
#' of units that satisfies user-supplied balance constraints on mean
#' differences. One of several available optimization programs can be used to
#' solve the mixed integer program. The default is the GLPK library as
#' implemented in the *Rglpk* package, but performance can be dramatically
#' improved using the HiGHS and the *highs* package, which are free, or Gurobi and the *gurobi* package, for which there is a
#' free academic license.
#'
#' This page details the allowable arguments with `method =
#' "cardinality"`. See [matchit()] for an explanation of what each argument
#' means in a general context and how it can be specified.
#'
#' Below is how `matchit()` is used for cardinality matching:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "cardinality",
#'         estimand = "ATT",
#'         exact = NULL,
#'         mahvars = NULL,
#'         s.weights = NULL,
#'         ratio = 1,
#'         verbose = FALSE,
#'         tols = .05,
#'         std.tols = TRUE,
#'         solver = "glpk",
#'         ...) }
#'
#' @param formula a two-sided [formula] object containing the treatment and
#' covariates to be balanced.
#' @param data a data frame containing the variables named in `formula`.
#' If not found in `data`, the variables will be sought in the
#' environment.
#' @param method set here to `"cardinality"`.
#' @param estimand a string containing the desired estimand. Allowable options
#' include `"ATT"`, `"ATC"`, and `"ATE"`. See Details.
#' @param exact for which variables exact matching should take place. Separate
#' optimization will occur within each subgroup of the exact matching
#' variables.
#' @param mahvars which variables should be used for pairing after subset selection. Can only be set when `ratio` is a whole number. See Details.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into the optimization. The balance constraints refer to the product of the
#' sampling weights and the matching weights, and the sum of the product of the
#' sampling and matching weights will be maximized.
#' @param ratio the desired ratio of control to treated units. Can be set to
#' `NA` to maximize sample size without concern for this ratio. See
#' Details.
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console.
#' @param \dots additional arguments that control the matching specification:
#' \describe{
#' \item{`tols`}{ `numeric`; a vector of imbalance
#' tolerances for mean differences, one for each covariate in `formula`.
#' If only one value is supplied, it is applied to all. See `std.tols`
#' below. Default is `.05` for standardized mean differences of at most
#' .05 for all covariates between the treatment groups in the matched sample.
#' }
#' \item{`std.tols`}{ `logical`; whether each entry in `tols`
#' corresponds to a raw or standardized mean difference. If only one value is
#' supplied, it is applied to all. Default is `TRUE` for standardized mean
#' differences. The standardization factor is the pooled standard deviation
#' when `estimand = "ATE"`, the standard deviation of the treated group
#' when `estimand = "ATT"`, and the standard deviation of the control
#' group when `estimand = "ATC"` (the same as used in
#' [summary.matchit()]).}
#' \item{`solver`}{ the name of solver to use to
#' solve the optimization problem. Available options include `"highs"`, `"glpk"`,
#' `"symphony"`, and `"gurobi"` for HiGHS (implemented in the *highs* package), GLPK (implemented in the
#' *Rglpk* package), SYMPHONY (implemented in the *Rsymphony*
#' package), and Gurobi (implemented in the *gurobi* package),
#' respectively. The differences between them are in speed and solving ability.
#' GLPK (the default) and HiGHS are the easiest to install, but Gurobi is recommended as
#' it consistently outperforms other solvers and can find solutions even when
#' others can't, and in less time. Gurobi is proprietary but can be used with a
#' free trial or academic license. SYMPHONY may not produce reproducible
#' results, even with a seed set.  }
#' \item{`time`}{ the maximum amount of
#' time before the optimization routine aborts, in seconds. Default is 120 (2
#' minutes). For large problems, this should be set much higher.  }
#' }
#'
#' The arguments `distance` (and related arguments), `replace`, `m.order`, and `caliper` (and related arguments) are ignored with a warning.
#'
#' @section Outputs:
#'
#' Most outputs described in [matchit()] are returned with
#' `method = "cardinality"`. Unless `mahvars` is specified, the `match.matrix` and `subclass`
#' components are omitted because no pairing or subclassification is done. When
#' `include.obj = TRUE` in the call to `matchit()`, the output of the
#' optimization function will be included in the output. When `exact` is
#' specified, this will be a list of such objects, one for each stratum of the
#' exact variables.
#'
#' @details
#' ## Cardinality and Profile Matching
#'
#' Two types of matching are
#' available with `method = "cardinality"`: cardinality matching and
#' profile matching.
#'
#' **Cardinality matching** finds the largest matched set that satisfies the
#' balance constraints between treatment groups, with the additional constraint
#' that the ratio of the number of matched control to matched treated units is
#' equal to `ratio` (1 by default), mimicking k:1 matching. When not all
#' treated units are included in the matched set, the estimand no longer
#' corresponds to the ATT, so cardinality matching should be avoided if
#' retaining the ATT is desired. To request cardinality matching,
#' `estimand` should be set to `"ATT"` or `"ATC"` and
#' `ratio` should be set to a positive integer. 1:1 cardinality matching
#' is the default method when no arguments are specified.
#'
#' **Profile matching** finds the largest matched set that satisfies balance
#' constraints between each treatment group and a specified target sample. When
#' `estimand = "ATT"`, it will find the largest subset of the control
#' units that satisfies the balance constraints with respect to the treated
#' group, which is left intact. When `estimand = "ATE"`, it will find the
#' largest subsets of the treated group and of the control group that are
#' balanced to the overall sample. To request profile matching for the ATT,
#' `estimand` should be set to `"ATT"` and `ratio` to `NA`.
#' To request profile matching for the ATE, `estimand` should be set to
#' `"ATE"` and `ratio` can be set either to `NA` to maximize the
#' size of each sample independently or to a positive integer to ensure that
#' the ratio of matched control units to matched treated treats is fixed,
#' mimicking k:1 matching. Unlike cardinality matching, profile matching
#' retains the requested estimand if a solution is found.
#'
#' Neither method involves creating pairs in the matched set, but it is
#' possible to perform an additional round of pairing within the matched sample
#' after cardinality matching or profile matching for the ATE with a fixed whole number
#' sample size ratio by supplying the desired pairing variables to `mahvars`. Doing so will trigger [optimal matching][method_optimal] using `optmatch::pairmatch()` on the Mahalanobis distance computed using the variables supplied to `mahvars`. The balance or composition of the matched sample will not change, but additional
#' precision and robustness can be gained by forming the pairs.
#'
#' The weights are scaled so that the sum of the weights in each group is equal
#' to the number of matched units in the smaller group when cardinality
#' matching or profile matching for the ATE, and scaled so that the sum of the
#' weights in the control group is equal to the number of treated units when
#' profile matching for the ATT. When the sample sizes of the matched groups
#' is the same (i.e., when `ratio = 1`), no scaling is done. Robust
#' standard errors should be used in effect estimation after cardinality or
#' profile matching (and cluster-robust standard errors if additional pairing
#' is done in the matched sample). See `vignette("estimating-effects")`
#' for more information.
#'
#' ## Specifying Balance Constraints
#'
#' The balance constraints are on
#' the (standardized) mean differences between the matched treatment groups for
#' each covariate. Balance constraints should be set by supplying arguments to
#' `tols` and `std.tols`. For example, setting `tols = .1` and
#' `std.tols = TRUE` requests that all the mean differences in the matched
#' sample should be within .1 standard deviations for each covariate. Different
#' tolerances can be set for different variables; it might be beneficial to
#' constrain the mean differences for highly prognostic covariates more tightly
#' than for other variables. For example, one could specify `tols = c(.001, .05), std.tols = c(TRUE, FALSE)`
#' to request that the standardized
#' mean difference for the first covariate is less than .001 and the raw mean
#' difference for the second covariate is less than .05. The values should be
#' specified in the order they appear in `formula`, except when
#' interactions are present. One can run the following code:
#'
#' \preformatted{MatchIt:::get_assign(model.matrix(~X1*X2 + X3, data = data))[-1]}
#'
#' which will output a vector of numbers and the variable to which each number
#' corresponds; the first entry in `tols` corresponds to the variable
#' labeled 1, the second to the variable labeled 2, etc.
#'
#' ## Dealing with Errors and Warnings
#'
#' When the optimization cannot be
#' solved at all, or at least within the time frame specified in the argument
#' to `time`, an error or warning will appear. Unfortunately, it is hard
#' to know exactly the cause of the failure and what measures should be taken
#' to rectify it.
#'
#' A warning that says `"The optimizer failed to find an optimal solution
#' in the time alotted. The returned solution may not be optimal."` usually
#' means that an optimal solution may be possible to find with more time, in
#' which case `time` should be increased or a faster solver should be
#' used. Even with this warning, a potentially usable solution will be
#' returned, so don't automatically take it to mean the optimization failed.
#' Sometimes, when there are multiple solutions with the same resulting sample
#' size, the optimizers will stall at one of them, not thinking it has found
#' the optimum. The result should be checked to see if it can be used as the
#' solution.
#'
#' An error that says `"The optimization problem may be infeasible."`
#' usually means that there is a issue with the optimization problem, i.e.,
#' that there is no possible way to satisfy the constraints. To rectify this,
#' one can try relaxing the constraints by increasing the value of `tols`
#' or use another solver. Sometimes Gurobi can solve problems that the other
#' solvers cannot.
#'
#' @seealso [matchit()] for a detailed explanation of the inputs and outputs of
#' a call to `matchit()`.
#'
#' *\CRANpkg{designmatch}*, which performs cardinality and profile matching with many more options and
#' more flexibility. The implementations of cardinality matching differ between
#' *MatchIt* and *designmatch*, so their results might differ.
#'
#' *\CRANpkg{optweight}*, which offers similar functionality but in the context of weighting rather
#' than matching.
#'
#' @references In a manuscript, you should reference the solver used in the
#' optimization. For example, a sentence might read:
#'
#' *Cardinality matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R with the optimization performed by HiGHs (Huangfu & Hall, 2018).*
#'
#' See `vignette("matching-methods")` for more literature on cardinality
#' matching.
#'
#' @examplesIf requireNamespace("highs", quietly = TRUE)
#' data("lalonde")
#'
#' #Choose your solver; "gurobi" is best, "highs" is free and
#' #easy to install
#' solver <- "highs"
#'
#' # 1:1 cardinality matching
#' m.out1 <- matchit(treat ~ age + educ + re74,
#'                   data = lalonde, method = "cardinality",
#'                   estimand = "ATT", ratio = 1,
#'                   tols = .2, solver = solver)
#' m.out1
#' summary(m.out1)
#'
#' # Profile matching for the ATT
#' m.out2 <- matchit(treat ~ age + educ + re74,
#'                   data = lalonde, method = "cardinality",
#'                   estimand = "ATT", ratio = NA,
#'                   tols = .2, solver = solver)
#' m.out2
#' summary(m.out2, un = FALSE)
#'
#' # Profile matching for the ATE
#' m.out3 <- matchit(treat ~ age + educ + re74,
#'                   data = lalonde, method = "cardinality",
#'                   estimand = "ATE", ratio = NA,
#'                   tols = .2, solver = solver)
#' m.out3
#' summary(m.out3, un = FALSE)
#' @examplesIf (requireNamespace("highs", quietly = TRUE) && requireNamespace("optmatch", quietly = TRUE))
#' # Pairing after 1:1 cardinality matching:
#' m.out1b <- matchit(treat ~ age + educ + re74,
#'                    data = lalonde, method = "cardinality",
#'                    estimand = "ATT", ratio = 1,
#'                    tols = .15, solver = solver,
#'                    mahvars = ~ age + educ + re74)
#'
#' # Note that balance doesn't change but pair distances
#' # are lower for the paired-upon variables
#' summary(m.out1b, un = FALSE)
#' summary(m.out1, un = FALSE)
#'
#' # In these examples, a high tol was used and
#' # few covariate matched on in order to not take too long;
#' # with real data, tols should be much lower and more
#' # covariates included if possible.
NULL

matchit2cardinality <-  function(treat, data, discarded, formula,
                                 ratio = 1, focal = NULL, s.weights = NULL,
                                 replace = FALSE, mahvars = NULL, exact = NULL,
                                 estimand = "ATT", verbose = FALSE,
                                 tols = .05, std.tols = TRUE,
                                 solver = "glpk", time = 1*60, ...){

  if (verbose) {
    cat("Cardinality matching... \n")
  }

  tvals <- unique(treat)
  nt <- length(tvals)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
  if (!is.null(focal)) {
    if (!focal %in% tvals) .err("`focal` must be a value of the treatment")
  }
  else if (estimand == "ATC") {
    focal <- min(tvals)
  }
  else {
    focal <- max(tvals)
  }

  lab <- names(treat)

  weights <- setNames(rep(0, length(treat)), lab)

  X <- get.covs.matrix(formula, data = data)

  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data), nam = lab, sep = ", ", include_vars = TRUE))

    cc <- do.call("intersect", lapply(tvals, function(t) as.integer(ex)[treat == t]))
    if (length(cc) == 0) .err("no matches were found")
  }
  else {
    ex <- gl(1, length(treat))
    cc <- 1
  }

  #Process mahvars
  if (!is.null(mahvars)) {
    if (!is.finite(ratio) || !chk::vld_whole_number(ratio)) {
      .err("`mahvars` can only be used with `method = \"cardinality\"` when `ratio` is a whole number")
    }
    rlang::check_installed("optmatch")
    mahcovs <- transform_covariates(mahvars, data = data, method = "mahalanobis",
                                    s.weights = s.weights, treat = treat,
                                    discarded = discarded)
    pair <- setNames(rep(NA_character_, length(treat)), lab)

    #Set max problem size to Inf and return to original value after match
    omps <- getOption("optmatch_max_problem_size")
    on.exit(options(optmatch_max_problem_size = omps))
    options(optmatch_max_problem_size = Inf)
  }
  else {
    pair <- NULL
  }

  #Process tols
  assign <- get_assign(X)

  chk::chk_numeric(tols)
  if (length(tols) == 1) {
    tols <- rep(tols, ncol(X))
  }
  else if (length(tols) == max(assign)) {
    tols <- tols[assign]
  }
  else if (length(tols) != ncol(X)) {
    .err("`tols` must have length 1 or the number of covariates. See `?method_cardinality` for details")
  }

  chk::chk_logical(std.tols)
  if (length(std.tols) == 1) {
    std.tols <- rep(std.tols, ncol(X))
  }
  else if (length(std.tols) == max(assign)) {
    std.tols <- std.tols[assign]
  }
  else if (length(std.tols) != ncol(X)) {
    .err("`std.tols` must have length 1 or the number of covariates. See `?method_cardinality` for details")
  }

  #Apply std.tols
  if (any(std.tols)) {
    sds <- {
      if (estimand == "ATE") {
        pooled_sd(X[, std.tols, drop = FALSE], t = treat,
                  w = s.weights, contribution = "equal")
      }
      else {
        sqrt(apply(X[treat == focal, std.tols, drop = FALSE], 2,
                   wvar, w = s.weights[treat == focal]))
      }
    }

    zero.sds <- sds < 1e-10

    X[,std.tols][,!zero.sds] <- scale(X[, std.tols, drop = FALSE][,!zero.sds, drop = FALSE],
                                      center = FALSE, scale = sds[!zero.sds])
  }

  opt.out <- setNames(vector("list", nlevels(ex)), levels(ex))

  for (e in levels(ex)[cc]) {
    if (nlevels(ex) > 1 && verbose) {
      cat(sprintf("Matching subgroup %s/%s: %s...\n",
                  match(e, levels(ex)[cc]), length(cc), e))
    }

    in.exact <- which(!discarded & ex == e)

    treat_in.exact <- treat[in.exact]
    out <- cardinality_matchit(treat = treat_in.exact,
                               X = X[in.exact,, drop = FALSE],
                               estimand = estimand, tols = tols,
                               s.weights = s.weights[in.exact],
                               ratio = ratio,
                               focal = focal, tvals = tvals,
                               solver = solver, time = time,
                               verbose = verbose)

    weights[in.exact] <- out[["weights"]]
    opt.out[[e]] <- out[["opt.out"]]

    if (!is.null(mahvars)) {
      mo <- eucdist_internal(mahcovs[in.exact[out[["weights"]] > 0],, drop = FALSE],
                             treat_in.exact[out[["weights"]] > 0])

      pm <- optmatch::pairmatch(mo,
                                controls = ratio,
                                data = data.frame(treat_in.exact))

      pair[names(pm)[!is.na(pm)]] <- paste(as.character(pm[!is.na(pm)]), e, sep = "|")
    }
  }

  if (!is.null(pair)) {
    psclass <- factor(pair)
    levels(psclass) <- seq_len(nlevels(psclass))
    names(psclass) <- names(treat)

    mm <- nummm2charmm(subclass2mmC(psclass, treat, focal = switch(estimand, "ATC" = 0, 1)),
                       treat)
  }
  else {
    mm <- psclass <- NULL
  }

  if (length(opt.out) == 1L) out <- out[[1]]

  res <- list(match.matrix = mm,
              subclass = psclass,
              weights = weights,
              obj = opt.out)

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"

  res
}

## Function to actually do the matching
cardinality_matchit <- function(treat, X, estimand = "ATT", tols = .05, s.weights = NULL,
                                ratio = 1, focal = NULL, tvals = NULL, solver = "highs",
                                time = 2*60, verbose = FALSE) {

  n <- length(treat)
  if (is.null(tvals)) tvals <- if (is.factor(treat)) levels(treat) else sort(unique(treat))
  nt <- length(tvals)

  #Check inputs
  if (is.null(s.weights)) s.weights <- rep(1, n)
  else for (i in tvals) s.weights[treat == i] <- s.weights[treat == i]/mean(s.weights[treat == i])

  if (is.null(focal)) focal <- tvals[length(tvals)]

  chk::chk_number(time)
  chk::chk_gt(time, 0)

  chk::chk_string(solver)
  solver <- match_arg(solver, c("highs", "glpk", "symphony", "gurobi"))
  rlang::check_installed(switch(solver, glpk = "Rglpk", symphony = "Rsymphony", gurobi = "gurobi",
                                highs = "highs"))

  #Select match type
  if (estimand == "ATE") match_type <- "profile_ate"
  else if (!is.finite(ratio)) match_type <- "profile_att"
  else match_type <- "cardinality"

  #Set objective and constraints
  if (match_type == "profile_ate") {
    #Find largest sample that matches full sample

    #Objective function: total sample size
    O <- c(
      s.weights, #weight for each unit
      rep(0, nt)       #slack coefs for each sample size (n1, n0)
    )

    #Constraint matrix
    target.means <- apply(X, 2, wm, w = s.weights)

    C <- matrix(0, nrow = nt * (1 + 2*ncol(X)), ncol = length(O))
    Crhs <- rep(0, nrow(C))
    Cdir <- rep("==", nrow(C))

    for (i in seq_len(nt)) {
      #Num in group i = ni
      C[i, seq_len(n)] <- s.weights * (treat == tvals[i])
      C[i, n + i] <- -1

      #Cov means must be less than target.means+tols/2
      r1 <- nt + (i - 1)*2*ncol(X) + 1:ncol(X)
      C[r1, seq_len(n)] <- t((treat==tvals[i])*s.weights*X)
      C[r1, n + i] <- -target.means-tols/2
      Cdir[r1] <- "<"

      #Cov means must be greater than target.means-tols/2
      r2 <- r1 + ncol(X)
      C[r2, seq_len(n)] <- t((treat==tvals[i])*s.weights*X)
      C[r2, n + i] <- -target.means+tols/2
      Cdir[r2] <- ">"
    }

    #If ratio != 0, constrain n0 to be ratio*n1
    if (nt == 2L && is.finite(ratio)) {
      C_ratio <- c(rep(0, n), rep(-1, nt))
      C_ratio[n + which(tvals == focal)] <- ratio
      C <- rbind(C, C_ratio)
      Crhs <- c(Crhs, 0)
      Cdir <- c(Cdir, "==")
    }

    #Coef types
    types <- c(rep("B", n), #Matching weights
               rep("C", nt)) #Slack coefs for matched group size

    lower.bound <- c(rep(0, n),
                     rep(1, nt))
    upper.bound <- c(rep(1, n),
                     rep(Inf, nt))
  }
  else if (match_type == "profile_att") {
    #Find largest control group that matches treated group

    nonf <- which(treat != focal)
    n0 <- length(nonf)
    tvals_ <- setdiff(tvals, focal)

    #Objective function: size of matched control group
    O <- c(
      rep(1, n0), #weights for each non-focal unit
      rep(0, nt - 1)  #slack coef for size of non-focal groups
    )

    #Constraint matrix
    target.means <- apply(X[treat==focal,,drop=FALSE], 2, wm, w = s.weights[treat==focal])
    #One row per constraint, one column per coef

    C <- matrix(0, nrow = (nt - 1) * (1 + 2*ncol(X)), ncol = length(O))
    Crhs <- rep(0, nrow(C))
    Cdir <- rep("==", nrow(C))

    for (i in seq_len(nt - 1)) {
      #Num in group i = ni
      C[i, seq_len(n0)] <- s.weights[nonf] * (treat[nonf] == tvals_[i])
      C[i, n0 + i] <- -1

      #Cov means must be less than target.means+tols
      r1 <- nt - 1 + (i - 1)*2*ncol(X) + 1:ncol(X)
      C[r1, seq_len(n0)] <- t((treat[nonf]==tvals_[i])*s.weights[nonf]*X[nonf,,drop = FALSE])
      C[r1, n0 + i] <- -target.means-tols
      Cdir[r1] <- "<"

      #Cov means must be greater than target.means-tols
      r2 <- r1 + ncol(X)
      C[r2, seq_len(n0)] <- t((treat[nonf]==tvals_[i])*s.weights[nonf]*X[nonf,,drop = FALSE])
      C[r2, n0 + i] <- -target.means+tols
      Cdir[r2] <- ">"
    }

    #Coef types
    types <- c(rep("B", n0), #Matching weights
               rep("C", nt - 1))  #Slack for num control matched

    lower.bound <- c(rep(0, n0),
                     rep(0, nt - 1))
    upper.bound <- c(rep(1, n0),
                     rep(Inf, nt - 1))
  }
  else if (match_type == "cardinality") {
    #True cardinality matching: find largest balanced sample
    if (nt > 2) ratio <- 1

    #Objective function: total sample size
    O <- c(
      s.weights, #weight for each unit
      0          #coef for treated sample size (n1)
    )

    #Constraint matrix
    t_combs <- combn(tvals, 2, simplify = FALSE)

    C <- matrix(0, nrow = nt + 2*ncol(X)*length(t_combs), ncol = length(O))
    Crhs <- rep(0, nrow(C))
    Cdir <- rep("==", nrow(C))

    for (i in seq_len(nt)) {
      #Num in group i = ni
      C[i, seq_len(n)] <- s.weights * (treat == tvals[i])
      C[i, n + 1] <- if (tvals[i] == focal) -1 else -ratio
    }

    for (j in seq_along(t_combs)) {
      t_comb <- t_combs[[j]]
      if (t_comb[2] == focal) t_comb <- rev(t_comb)

      r1 <- nt + (j - 1)*2*ncol(X) + 1:ncol(X)
      C[r1, seq_len(n)] <- t(((treat==t_comb[1]) - (treat==t_comb[2])/ratio)*s.weights*X)
      C[r1, n + 1] <- -tols
      Cdir[r1] <- "<"

      r2 <- r1 + ncol(X)
      C[r2, seq_len(n)] <- t(((treat==t_comb[1]) - (treat==t_comb[2])/ratio)*s.weights*X)
      C[r2, n + 1] <- tols
      Cdir[r2] <- ">"
    }

    #Coef types
    types <- c(rep("B", n), #Matching weights
               rep("C", 1)) #Slack coef for treated group size (n1)

    lower.bound <- c(rep(0, n),
                     rep(0, 1))
    upper.bound <- c(rep(1, n),
                     rep(min(tabulateC(treat)), 1))
  }

  weights <- NULL

  opt.out <- dispatch_optimizer(solver = solver, obj = O, mat = C, dir = Cdir,
                                rhs = Crhs, types = types, max = TRUE, lb = lower.bound,
                                ub = upper.bound, time = time, verbose = verbose)

  cardinality_error_report(opt.out, solver)

  sol <- switch(solver,
                "glpk" = opt.out$solution,
                "symphony" = opt.out$solution,
                "gurobi" = opt.out$x,
                "highs" = opt.out$primal_solution)

  if (match_type %in% c("profile_ate", "cardinality")) {
    weights <- round(sol[seq_len(n)])
  }
  else if (match_type %in% c("profile_att")) {
    weights <- rep(1, n)
    weights[treat != focal] <- round(sol[seq_len(n0)])
  }

  #Make sure sum of weights in both groups is the same (important for exact matching)
  if (match_type == "profile_att" && (is.na(ratio) || ratio != 1)) {
    for (t in setdiff(tvals, focal)) {
      weights[treat == t] <- weights[treat == t]*sum(weights[treat == focal])/sum(weights[treat == t])
    }
  }
  else {
    smallest.group <- tvals[which.min(vapply(tvals, function(t) sum(treat == t), numeric(1L)))]
    for (t in setdiff(tvals, smallest.group)) {
      weights[treat == t] <- weights[treat == t]*sum(weights[treat == smallest.group])/sum(weights[treat == t])
    }
  }

  list(weights = weights, opt.out = opt.out)
}

cardinality_error_report <- function(out, solver) {
  if (solver == "glpk") {
    if (out$status == 1) {
      if (all(out$solution == 0)) {
        .err("the optimization problem may be infeasible. Try increasing the value of `tols`.\nSee `?method_cardinality` for additional details")
      }
      .wrn("the optimizer failed to find an optimal solution in the time alotted. The returned solution may not be optimal.\nSee `?method_cardinality` for additional details")
    }
  }
  else if (solver == "symphony") {
    if (names(out$status) %in% c("TM_TIME_LIMIT_EXCEEDED") && !all(out$solution == 0) && all(out$solution <= 1)) {
      .wrn("the optimizer failed to find an optimal solution in the time alotted. The returned solution may not be optimal")
    }
    else if (names(out$status) != "TM_OPTIMAL_SOLUTION_FOUND") {
      .err("the optimizer failed to find an optimal solution in the time alotted. The optimization problem may be infeasible. Try increasing the value of 'tols'.\nSee `?method_cardinality` for additional details")
    }
  }
  else if (solver == "gurobi") {
    if (out$status %in% c("TIME_LIMIT", "SUBOPTIMAL") && !all(out$x == 0)) {
      .wrn("the optimizer failed to find an optimal solution in the time alotted. The returned solution may not be optimal.\nSee `?method_cardinality` for additional details")
    }
    else if (out$status %in% c("INFEASIBLE", "INF_OR_UNBD", "NUMERIC") || all(out$x == 0)) {
      .err("The optimization problem may be infeasible. Try increasing the value of `tols`.\nSee `?method_cardinality` for additional details")
    }
  }
  else if (solver == "highs") {
    if (out$status_message %in% c("Infeasible", "Primal infeasible or unbounded")) {
      # if (out$status_message %in% c("Infeasible", "Primal infeasible or unbounded") ||
      #     all(abs(out$primal_solution) < 1e-8)) {
      .err("the optimization problem may be infeasible. Try increasing the value of `tols`.\nSee `?method_cardinality` for additional details")
    }
    if (out$status_message %in% c("Time limit reached", "Iteration limit reached")) {
      .err("the optimizer failed to find an optimal solution in the time alotted. Try increasing the value of `time`.\nSee `?method_cardinality` for additional details")
    }
  }
}

dispatch_optimizer <- function(solver = "glpk", obj, mat, dir, rhs, types, max = TRUE, lb = NULL, ub = NULL, time = NULL, verbose = FALSE) {
  if (solver == "glpk") {
    dir[dir == "="] <- "=="
    opt.out <- Rglpk::Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs, max = max,
                                     types = types,
                                     # bounds = list(lower = lb, upper = ub), #Spurious warning when using bounds
                                     control = list(tm_limit = time*1000, verbose = verbose))
  }
  else if (solver == "symphony") {
    dir[dir == "<"] <- "<="
    dir[dir == ">"] <- ">="
    dir[dir == "="] <- "=="
    opt.out <- Rsymphony::Rsymphony_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs, max = TRUE,
                                             types = types, verbosity = verbose - 2,
                                             # bounds = list(lower = lb, upper = ub), #Spurious warning when using bounds
                                             time_limit = time)
  }
  else if (solver == "gurobi") {
    dir[dir == "<="] <- "<"
    dir[dir == ">="] <- ">"
    dir[dir == "=="] <- "="
    opt.out <- gurobi::gurobi(list(A = mat, obj = obj, sense = dir, rhs = rhs, vtype = types,
                                   modelsense = "max", lb = lb, ub = ub),
                              params = list(OutputFlag = as.integer(verbose), TimeLimit = time))
  }
  else if (solver == "highs") {
    rhs_h <- lhs_h <- rhs

    rhs_h[dir == ">"] <- Inf
    lhs_h[dir == "<"] <- -Inf

    types[types == "B"] <- "I"

    opt.out <- highs::highs_solve(L = obj, lower = lb, upper = ub,
                                  A = mat, lhs = lhs_h, rhs = rhs_h,
                                  types = types,
                                  maximum = max,
                                  control = list(time_limit = time,
                                                 log_to_console = verbose))
  }

  opt.out
}
