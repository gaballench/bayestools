#' Function for estimation of lognormal parameters through quadratic optimization
#' 
#' @param p1,p2,p3 Quantile probabilities, 0.05, 0.50 and 0.95 by default
#'
#' @param q1,q2,q3 Observed quantiles, might come from a HPD from a previous study (along with a median), or from biostratigraphic constrains.
#'
#' @param params Internal params to be optimized by the optim funciton on the internal .quadraticFun function
#'
#' @param output One of two possible values (unconstrained internally): "complete" and "parameters". For the latter the complete output of the optim function is returned with information on convergence and squared errors (that might be useless for simple cases) or just the parameters
#'
#' @return Either a list with the complete output of convergence, squared errors and parameter values, or just a vector of parameter values. Depends on the value of \code{output}
#'
#' @examples
#' # Find the best parameters for a lognormal density that fit the observed quantiles q1 = 1, q2 = 5, and q3 = 10, providing full output for the calculations in the form of a list
#' LNparams(q1 = 1, q2 = 5, q3 = 10)
#' # As above but return only the parameter values instead
#' LNparams(q1 = 1, q2 = 5, q3 = 10, output = "parameters")
#'
#' @note Notes on distributions
#' pbeta(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#' pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
#' pcauchy(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
#' pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#' pexp(q, rate = 1, lower.tail = TRUE, log.p = FALSE)
#' pf(q, df1, df2, ncp, lower.tail = TRUE, log.p = FALSE)
#' pgamma(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
#' pgeom(q, prob, lower.tail = TRUE, log.p = FALSE)
#' phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
#' plnorm(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
#' pnbinom(q, size, prob, mu, lower.tail = TRUE, log.p = FALSE)
#' pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
#' ppois(q, lambda, lower.tail = TRUE, log.p = FALSE)
#' pt(q, df, ncp, lower.tail = TRUE, log.p = FALSE)
#' punif(q, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)
#' pweibull(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
#' #for the dirichlet distribution there is the Compositional::diri.est
LNparams <- function(p1 = 0.05, p2 = 0.50, p3 = 0.95, q1, q2, q3, params, output = "complete") {
    .quadraticFun <- function(pars){
        (plnorm(q1, pars[1], pars[2]) - p1)^2 + (plnorm(q2, pars[1], pars[2]) - p2)^2 + (plnorm(q3, pars[1], pars[2]) - p3)^2
    }
    res <- optim(c(1, 1), .quadraticFun)
    # if output is set to "complete", return the whole res object, otherwise, return just the parameters 
    if (output == "parameters") {
        return(res$par)
    }
    return(res)
}

# prototype for using ellipsis args, it worked!
encap <- function(FUN, ...) {
    sexpr <- substitute(FUN)
    dots <- substitute(...())
    expr <- as.call(c(as.list(sexpr), dots))
    print(expr)
    print(length(dots))
    cat(eval(expr), "\n")
    cat(eval(substitute(expr - 1)), "\n")
}

# params MUST be a list of CHARACTER mode in the same sense as in aggregate()
findParams <- function(p1 = 0.05, p2 = 0.50, p3 = 0.95, q1, q2, q3, output = "complete", densit, params) {
    densit <- substitute(densit)
    params <- sapply(params, as.name)
    densiCall1 <- as.call(c(as.list(densit), as.list(q1), params))
    densiCall2 <- as.call(c(as.list(densit), as.list(q2), params))
    densiCall3 <- as.call(c(as.list(densit), as.list(q3), params))
    quadratEq <- substitute((densiCall1 - p1)^2 +
                            (densiCall2 - p2)^2 +
                            (densiCall3 - p3)^2)
    quadraticFun <- function(params) {
        eval(quadratEq)
    }    
    initVals <- rep(1, times = length(params))
    res <- optim(initVals, quadraticFun) # ERROR SOMEWHERE AROUD HERE...
    # if output is set to "complete", return the whole res object, otherwise, return just the parameters 
    if (output == "parameters") {
        return(res$par)
    }
    return(res)
}
