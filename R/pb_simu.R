#' Simulate biomarkers from Normal or Log-Normal Distribution
#'
#'
#'
#'
#'
#' @export
#'
pbSimuBiom <- function(n, mu = 0, sd = 1, type = c("lognormal", "normal"),
                       trunc.lb = -Inf, trunc.ub = Inf) {

    type <- match.arg(type);

    f.lognormal <- function() {
        stopifnot(trunc.ub > 0);
        trunc.ub <- log(trunc.ub);

        if (trunc.lb > 0) {
            trunc.lb <- log(trunc.lb);
        }

        rst <- truncnorm::rtruncnorm(n, mean = mu, sd = sd, a = trunc.lb, b = trunc.ub);

        exp(rst);
    }

    rst <- switch(type,
                  normal    = truncnorm::rtruncnorm(n, mean = mu, sd = sd,
                                                    a = trunc.lb, b = trunc.ub),
                  lognormal = f.lognormal());

    rst;

}


#' Simulate responses given biomarkers
#'
#'
#'
#'
#'
#' @export
#'
pbSimuResp <- function(x, betas, type = c("fourpar", "logit")) {

    type <- match.arg(type);
    py1  <- switch(type,
                   fourpar = betas[1] + (betas[2] - betas[1]) / (1 + (betas[3]/x)^betas[4]),
                   logit   = betas[1] + betas[2]*x
                   )
    py1  <- py1 / (1 + exp(py1));
    rst  <- rbinom(length(x), 1, py1);

    list(y   = rst,
         py1 = py1);
}


