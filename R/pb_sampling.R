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
pbSmpBiom <- function(x.count, prior.q = NULL, iter = 4000) {
    n.cat <- length(x.count);
    if (is.null(prior.q)) {
        prior.q <- rep(0.5, n.cat);
    }

    post.q <- rdirichlet(iter, prior.q + x.count);

    ## get cumulative probabilities
    nq        = ncol(post.q);
    post.cumu = apply(post.q, 1, function(x) get.cumu(x, nq));

    list(post.q = post.q,
         cumu.q = t(post.cumu));
}


#' Posterior response rates for each cut interval
#'
#' @param x.cut     observed biomarker levels
#' @param y         response status
#' @param cand.cuts candidate cut points
#'
#'
#' @export
pbSmpResp <- function(x.cut, y, cand.cuts = NULL, type = c("cumubin", "simplebin"),
                      iter = 4000, ...) {
    type <- match.arg(type);

    if (is.null(cand.cuts))
        cand.cuts <- 1:max(x.cut);

    rst  <- switch(type,
                   simplebin = prvSmpSimpleBin(x.cut, y, cand.cuts, cumu = FALSE,
                                               iter = iter, ...),
                   cumubin   = prvSmpSimpleBin(x.cut, y, cand.cuts, cumu = TRUE,
                                             iter = iter, ...)
                   );

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

    cumu.q = post.q$cumu.q;
    post.q = post.q$post.q;

    if (!attr(post.p, "cumu")) {
        nc = ncol(post.p);
        stopifnot(nc == ncol(post.q));

        pq      <- post.q * post.p;
        cumu.pq <- apply(pq, 1, function(x) get.cumu(x, nc));
        cumu.p  <- t(cumu.pq) / cumu.q;
    } else {
        cumu.p = post.p;
    }

    list(cumu.q = cumu.q,
         cumu.p = cumu.p);
}


#' Private Function: Sampling from simple binomial of response rates
#'
#'
#'
prvSmpSimpleBin <- function(x.cut, y, cand.cuts, iter = 4000, cumu = FALSE,
                            prior.p = c(a = 0.0025, b = 0.0025), ...) {
    rst <- NULL;
    for (ct in cand.cuts) {
        if (cumu) {
            cur.inx <- which(x.cut >= ct);
        } else {
            cur.inx <- which(ct == x.cut);
        }

        cur.n <- length(cur.inx);
        if (0 == cur.n) {
            cur.a <- 0;
            cur.b <- 0;
        } else {
            cur.a <- sum(y[cur.inx]);
            cur.b <- cur.n - cur.a;
        }

        cur.rst <- rbeta(n = iter, cur.a + prior.p["a"], cur.b + prior.p["b"]);
        rst     <- cbind(rst, cur.rst);
    }

    attr(rst, "cumu") = cumu;
    colnames(rst)     = cand.cuts;
    rst
}


## get cumulative sums
get.cumu  <-  function(q, nq) {
    q = q[nq:1];
    q = cumsum(q);
    q[nq:1]
}
