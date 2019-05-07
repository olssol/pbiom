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

    cuts.inf <- unique(c(-Inf, cuts, Inf));
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


#' Posterior response rates for each cut interval
#'
#' @param x.cut     observed biomarker levels
#' @param y         response status
#' @param cand.cuts candidate cut points
#'
#'
#' @export
pbSmpResp <- function(x.cut, y, cand.cuts = NULL, type = c("simplebin"), iter = 4000, ...) {
    type <- match.arg(type);

    if (is.null(cand.cuts))
        cand.cuts <- 1:max(x.cut);

    rst  <- switch(type,
                   simplebin = prvSmpSimplebin(x.cut, y, cand.cuts,
                                               iter = iter, ...));

    rst
}


#' Posterior cumulative response rate
#'
#' @param post.q biomarker interval probabilities
#' @param post.p response rates
#'
#' @return biomarker probabilities and response rates conditioning on biomarker
#'     values bigger than cut points
#'
#' @export
pbCumuPQ <- function(post.q, post.p) {

    nc <- ncol(post.p);
    stopifnot(nc == ncol(post.q));

    post.q <- post.q[,nc:1];
    post.p <- post.p[,nc:1];
    pq     <- post.q * post.p;

    cumu.q  <- apply(post.q, 1, cumsum);
    cumu.pq <- apply(pq,     1, cumsum);
    cumu.p  <- cumu.pq / cumu.q;

    list(cumu.q = t(cumu.q)[,nc:1],
         cumu.p = t(cumu.p)[,nc:1]);
}

#' Private Function: Sampling from simple binomial of response rates
#'
#'
#'
prvSmpSimplebin <- function(x.cut, y, cand.cuts, iter = 4000,
                            prior.q = c(a = 0.5, b = 0.5), ...) {
    rst <- NULL;
    for (ct in cand.cuts) {
        cur.inx <- which(ct == x.cut);
        cur.n   <- length(cur.inx);
        if (0 == cur.n) {
            cur.a <- 0;
            cur.b <- 0;
        } else {
            cur.a <- sum(y[cur.inx]);
            cur.b <- cur.n - cur.a;
        }

        cur.rst <- rbeta(n = iter, cur.a + prior.q["a"], cur.b + prior.q["b"]);
        rst     <- cbind(rst, cur.rst);
    }
    colnames(rst) <- cand.cuts;

    rst
}
