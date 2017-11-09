#' Function for estimation of lognormal parameters through quadratic optimization
#' 
#' @param p1,p2,p3 Quantile probabilities, 0.05, 0.50 and 0.95 by default
#'
#' @param q1,q2,q3 Observed quantiles, might come from a HPD from a previous study (along with a median), or from biostratigraphic constrains.
#'
#' @param params Internal params to be optimized by the optim funciton on the internal .quadraticFun function. THIS ARGUMENT DOES NOT SHOW ANYWHERE IN THE FUNCTION BODY!
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
