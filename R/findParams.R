#' Function for estimation of probability density function parameters through quadratic optimization
#' 
#' @param p A numeric vector of quantile probabilities, 0.05, 0.50 and 0.95 by default
#'
#' @param q A numeric vector of observed quantiles, might come from a HPD from a previous study (along with a median), or from other sources of prior information. See \link{Notes}.
#'
#' @param pdfunction A character vector (of length one) with the name of the PDF function of interest. Technically this argument supports any PDF function of the form pDIST (e.g., pnorm, ppois, pexp).
#' 
#' @param params A character vector with the name of the parameter(s) to optimize in the probability density function. These should match the parameter names of the respective PDF function, e.g., \code{"lambda"} in the function \code{ppois}
#'
#' @param output One of two possible values: \code{"complete"} and \code{"parameters"}. For the latter the complete output of the \code{optim} function is returned with information on convergence and squared errors (that might be useless for simple cases) or just the parameters.
#'
#' @return Either a list with the complete output of convergence, squared errors and parameter values, or just a vector of parameter values. Depends on the value of \code{output}.
#' Warnings may be triggered by the function \code{optim} since the optimization is a heuristic process, whenever a given iteration results in an invalid value for a given combination of parameters, the \code{optim} function tries another combination of values but inform the user about the problem through a warning. In general these can be safely disregarded.
#'
#' @author Main code by Gustavo A. Ballen with important contributions in expression call structure and vectorized design by Klaus Schliep (\email{Klaus.Schliep@umb.edu}).
#' 
#' @examples
#' # Find the best parameters for a standard normal density that fit the observed quantiles -1.644854, 0, and 1.644854, providing full output for the calculations in the form of a list
#' findParams(q = c(-1.644854, 0, 1.644854), p = c(0.05,  0.50, 0.95), output = "complete", pdfunction = "pnorm", params = c("mean", "sd"))
#'
#' # Given that we have prior on the age of a fossil to be 1 - 10 Ma and that we want to model it with a lognormal distribution, fin the parameters of the PDF that best reflect the uncertainty in question (i.e., the parameters  for which the observed quantiles are 1, 5.5, and 10, assuming that we want the midpoint to reflect the mean of the PDF.
#' findParams(q = c(1, 5.5, 10), p = c(0.05,  0.50, 0.95), output = "complete", pdfunction = "plnorm", params = c("meanlog", "sdlog"))
#' 
#' @note Notes on distributions

findParams <- function(q, p, output = "complete", pdfunction, params) {
    # define the initial values for optimization
    initVals = rep(mean(q), times = length(params))
    # construct the call as a list
    l <- length(params)
    cl <- vector("list", 2  + length(params))
    cl[[1]] <- as.name(pdfunction)
    cl[[2]] <- q
    names(cl) <- c(NA, "q", params)
    mode(cl) <- "call"
    # the quadratic function will minimize the error around estimates
    quadraticFun <- function(x) {
        cl[3:(l+2)] <- x
        res <- eval(cl)
        sum((res - p)^2)
    }    
    # pick the method for optimization based on the number of parameters
    if (l == 1) {
        res <- optim(initVals, quadraticFun, method = "BFGS")
    } else {
        res <- optim(initVals, quadraticFun)
    }
    # return depending on the desired output    
    if (output == "parameters") {
        return(res$par)
    } else {
        return(res)
    }
}
