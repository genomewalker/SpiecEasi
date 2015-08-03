# Added multicore support for glasso method
#' @export
huge.mc <- function (x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, 
                     method = "mb", scr = NULL, scr.num = NULL, cov.output = FALSE, 
                     sym = "or", verbose = TRUE, ncores = 1, mc.progress = FALSE) 
{
  gcinfo(FALSE)
  est = list()
  est$method = method
  if (method == "glasso") {
    fit = huge.glasso.mc(x, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, 
                      lambda = lambda, scr = scr, cov.output = cov.output, 
                      verbose = verbose, ncores = ncores, mc.progress = mc.progress)
    est$path = fit$path
    est$lambda = fit$lambda
    est$icov = fit$icov
    est$df = fit$df
    est$sparsity = fit$sparsity
    est$loglik = fit$loglik
    if (cov.output) 
      est$cov = fit$cov
    est$cov.input = fit$cov.input
    est$cov.output = fit$cov.output
    est$scr = fit$scr
    rm(fit)
    gc()
  }
  est$data = x
  rm(x, scr, lambda, lambda.min.ratio, nlambda, cov.output, 
     verbose)
  gc()
  class(est) = "huge"
  return(est)
}

#' @export
huge.glasso.mc<-function (x, lambda = NULL, lambda.min.ratio = NULL, nlambda = NULL, 
                          scr = NULL, cov.output = FALSE, verbose = TRUE, ncores = 1, mc.progress = FALSE) 
{
  gcinfo(FALSE)
  n = nrow(x)
  d = ncol(x)
  fit = list()
  fit$cov.input = isSymmetric(x)
  if (fit$cov.input) {
    if (verbose) 
      cat("The input is identified as the covriance matrix.\n")
    S = x
  }
  if (!fit$cov.input) {
    x = scale(x)
    S = cor(x)
  }
  #rm(x)
  gc()
  if (is.null(scr)) {
    scr = FALSE
  }
  fit$scr = scr
  if (!is.null(lambda)) {
    nlambda = length(lambda)
  }
  if (is.null(lambda)) {
    if (is.null(nlambda)) {
      nlambda = 10
    }
    if (is.null(lambda.min.ratio)) {
      lambda.min.ratio = 0.1
    }
    lambda.max = max(max(S - diag(d)), -min(S - diag(d)))
    lambda.min = lambda.min.ratio * lambda.max
    lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }
  fit$lambda = lambda
  fit$loglik = rep(-d, nlambda)
  fit$sparsity = rep(0, nlambda)
  fit$df = rep(0, nlambda)
  fit$path = list()
  fit$icov = list()
  fit$cov.output = cov.output
  if (cov.output) {
    fit$cov = list()
  }
  out.glasso = NULL
  fit1<-list()
  fit.mc<-mclapply.pb(nlambda:1, function(i) {
    fit1$sparsity<-0
    fit1$loglik<--d
    fit1$df<-0
    z = which(rowSums(abs(S) > lambda[i]) > 1)
    q = length(z)
    if (q > 0) {
      #       if (verbose) {
      #         if (scr) {
      #           cat(paste(c("Conducting the graphical lasso (glasso) wtih lossy screening....in progress:", 
      #                       floor(100 * (1 - i/nlambda)), "%"), collapse = ""), 
      #               "\r")
      #         }
      #         if (!scr) {
      #           cat(paste(c("Conducting the graphical lasso (glasso) with lossless screening....in progress:", 
      #                       floor(100 * (1 - i/nlambda)), "%"), collapse = ""), 
      #               "\r")
      #         }
      #         flush.console()
      #       }
      if (scr) {
        if (!is.null(out.glasso)) {
          out.glasso = .C("hugeglassoscr",
                          S = as.double(S[z,z]),
                          W = as.double(tmp.cov[z, z]), 
                          T = as.double(tmp.icov[z,z]), 
                          dd = as.integer(q), 
                          lambda = as.double(lambda[i]), 
                          df = as.integer(0), 
                          PACKAGE = "huge")
        }
        if (is.null(out.glasso)) {
          out.glasso = .C("hugeglassoscr",
                          S = as.double(S[z,z]), 
                          W = as.double(S[z, z]),
                          T = as.double(diag(q)), 
                          dd = as.integer(q),
                          lambda = as.double(lambda[i]), 
                          df = as.integer(0),
                          PACKAGE = "huge")
        }
      }else {
        if (!is.null(out.glasso)) {
          out.glasso = .C("hugeglasso", 
                          S = as.double(S[z,z]),
                          dd = as.integer(q), 
                          lambda = as.double(lambda[i]), 
                          df = as.integer(0), 
                          PACKAGE = "huge")
        }
        if (is.null(out.glasso)) {
          out.glasso = .C("hugeglasso", 
                          S = as.double(S[z,z]),
                          W = as.double(S[z, z]), 
                          T = as.double(diag(q)), 
                          dd = as.integer(q), 
                          lambda = as.double(lambda[i]), 
                          df = as.integer(0), 
                          PACKAGE = "huge")
        }
      }
      out.glasso$T = matrix(out.glasso$T, ncol = q)
      out.glasso$W = matrix(out.glasso$W, ncol = q)
    }
    if (q == 0){ 
      out.glasso = NULL
    }
    tmp.icov = matrix(0, d, d)
    diag(tmp.icov) = 1/(diag(S) + lambda[i])
    tmp.cov = matrix(0, d, d)
    diag(tmp.cov) = diag(S) + lambda[i]
    fit1$path <- Matrix(0, d, d)
    if (!is.null(out.glasso)) {
      tmp.icov[z, z] = out.glasso$T
      tmp.cov[z, z] = out.glasso$W
      fit1$path[z, z] <- abs(sign(out.glasso$T))
      diag(fit1$path) <- 0
      fit1$sparsity <- as.double(out.glasso$df)/d/(d - 1)
      fit1$df <- out.glasso$df/2
      fit1$loglik <- (log(det(out.glasso$T)) - sum(diag(out.glasso$T %*% S[z, z])) - (d - q))
    }
    fit1$icov <- Matrix(tmp.icov)
    if (cov.output) 
      fit1$cov <- Matrix(tmp.cov)
    return(fit1)
  }, mc.cores = ncores, mc.progress = mc.progress)
  
  fit.l<-do.call(Map, c(c, fit.mc))
  fit$loglik<-rev(fit.l$loglik)
  fit$sparsity<-rev(fit.l$sparsity)
  fit$df<-rev(fit.l$df)
  fit$path<-rev(fit.l$path)
  fit$icov<-rev(fit.l$icov)
  fit$cov<-rev(fit.l$cov)
  
  rm(S, out.glasso, fit.l, fit.mc)
  gc()
  if (verbose) {
    cat("Conducting the graphical lasso (glasso)....done.                                          \r")
    cat("\n")
    flush.console()
  }
  return(fit)
}


##------------------------------------------------------------------------------
##' Wrapper around mclapply to track progress
##' 
##' Based on http://stackoverflow.com/questions/10984556
##' 
##' @param X         a vector (atomic or list) or an expressions vector. Other
##'                  objects (including classed objects) will be coerced by
##'                  ‘as.list’
##' @param FUN       the function to be applied to
##' @param ...       optional arguments to ‘FUN’
##' @param mc.preschedule see mclapply
##' @param mc.set.seed see mclapply
##' @param mc.silent see mclapply
##' @param mc.cores see mclapply
##' @param mc.cleanup see mclapply
##' @param mc.allow.recursive see mclapply
##' @param mc.progress track progress?
##' @param mc.style    style of progress bar (see txtProgressBar)
##'
##' @examples
##' x <- mclapply2(1:1000, function(i, y) Sys.sleep(0.01))
##' x <- mclapply2(1:3, function(i, y) Sys.sleep(1), mc.cores=1)
##------------------------------------------------------------------------------
#' @export
mclapply.pb <- function(X, FUN, ..., 
                      mc.preschedule = TRUE, mc.set.seed = TRUE,
                      mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                      mc.cleanup = TRUE, mc.allow.recursive = TRUE,
                      mc.progress=TRUE, mc.style=3) 
{
  if (!is.vector(X) || is.object(X)) X <- as.list(X)
  
  if (mc.progress) {
    f <- fifo(tempfile(), open="w+b", blocking=T)
    p <- parallel:::mcfork()
    pb <- txtProgressBar(0, length(X), style=mc.style)
    setTxtProgressBar(pb, 0) 
    progress <- 0
    if (inherits(p, "masterProcess")) {
      while (progress < length(X)) {
        readBin(f, "double")
        progress <- progress + 1
        setTxtProgressBar(pb, progress) 
      }
      cat("\n")
      parallel:::mcexit()
    }
  }
  tryCatch({
    result <- mclapply(X, function(...) {
      res <- FUN(...)
      if (mc.progress) writeBin(1, f)
      res
    }, 
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
    mc.silent = mc.silent, mc.cores = mc.cores,
    mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
    )
    
  }, finally = {
    if (mc.progress) close(f)
  })
  result
}