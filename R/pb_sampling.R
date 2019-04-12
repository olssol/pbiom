#' Cut Biomarker using candidate cut points
#'
#'
#'
#' @export
#'
pbCutBiom <- function(x, cuts = NULL, probs = c(0.25, 0.5, 0.75)) {
    if (is.null(cuts)) {
        cuts <- quantile(x, probs = probs);
    }

    cuts.inf <- c(-Inf, cuts, Inf);
    cut.x    <- cut(x, breaks = cuts.inf, labels = 1:(length(cuts.inf)-1));
    cou.x    <- table(cut.x)

    rst <- list(x.cut    = as.numeric(cut.x),
                x.count  = cou.x,
                cuts     = cuts);
}



#' Posterior distribution of biomarker as intervals
#'
#'
#' @export
#'
pbSmpBiom <- function(x.count, prior.p = NULL, iter = 4000) {
    n.cat <- length(x.count);
    if (is.null(prior.p)) {
        prior.p <- rep(0.5, n.cat);
    }

    post.p <- prior.p + x.count;
    rst    <- rdirichlet(iter, post.p);
}


#' Posterior joint sample cut points and step response rates
#'
#' @param x.cut     biomarker levels
#' @param y         response status
#' @param cand.cuts candidate cut points
#'
#'
#' @export
pbSmpResp <- function(x.cut, y, cand.cuts, type = c("simplebin"), ...) {
    type <- match.arg(type);
    rst  <- switch(type,
                   simplebin = prvSmpSimplebin(x.cut, y, cand.cuts, ...));

    rst
}

#' Private Function: Sampling from simple binomial of response rates
#'
#'
#'
prvSmpSimplebin <- function(x.cut, y, cand.cuts, iter = 4000, pri.a = 0.5, pri.b = 0.5) {
    rst <- NULL;
    for (ct in cand.cuts) {
        cur.inx <- which(ct <= x.cut);
        cur.n   <- length(cur.inx);
        if (0 == cur.n) {
            cur.a <- 0;
            cur.b <- 0;
        } else {
            cur.a <- length(which(0 == y[cur.inx]));
            cur.b <- cur.n - cur.a;
        }

        cur.rst <- rbeta(n = iter, cur.a + pri.a, cur.b + pri.b);
        rst     <- cbind(rst, cur.rst);
    }
    colnames(rst) <- cand.cuts;

    rst
}
