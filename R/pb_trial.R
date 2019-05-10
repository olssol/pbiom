#' Simulate a single trial
#'
#' @param par.biom list of parameters for simulating biomarkers
#' @param par.resp list of parameters for simulating responses
#' @param iter number of iterations for interim analysis
#' @param nlarge large n for evaluating the truth
#'
#' @export
#'
pbSimuSingleTrial <- function(par.biom, par.resp, n2, theta0,
                              prior.p = NULL, prior.q = c(a = 0.5, b = 0.5),
                              cut.x = NULL, cut.quants = c(0.25, 0.5, 0.75), cut.known = TRUE,
                              cand.cuts = 1:(length(cut.quants)+1),
                              resp.mdl = "simplebin",
                              uti.f, alpha = 0.05, uti.cut = 0.1,
                              B1 = 1, C1 = 1, C2 = 0, C3 = 0,
                              iter = 4000, nlarge = 50000,
                              repeach = 1) {
    ## find cut points
    if (is.null(cut.quants)) {
        stop("Biomarker level cut quantiles (cut.quants) needs to be provided.");
    }

    ## truth
    true.par.biom    <- par.biom;
    true.par.biom$n  <- nlarge;
    true.x           <- do.call(pbSimuBiom, true.par.biom);
    true.y           <- do.call(pbSimuResp, c(list(x = true.x), par.resp))$y;

    if (is.null(cut.x) & cut.known) {
        cut.x <- quantile(true.x, probs = cut.quants);
    }

    ## simulate first stage;
    s1x    <- do.call(pbSimuBiom, par.biom);
    s1y    <- do.call(pbSimuResp, c(list(x = s1x), par.resp))$y;
    n1     <- par.biom$n;
    nresp1 <- sum(s1y);

    ## interim analysis
    s1x.cut <- pbCutBiom(s1x, cuts = cut.x, probs = cut.quants);
    post.q  <- pbSmpBiom(s1x.cut$x.count, prior.p = prior.p, iter = iter);
    post.p  <- pbSmpResp(s1x.cut$x.cut, s1y, cand.cuts = cand.cuts, type = resp.mdl);
    cumu.pq <- pbCumuPQ(post.q, post.p);

    ## truth pq
    ## use the current cut points for the truth
    true.x.cut   <- pbCutBiom(true.x, cuts = s1x.cut$cuts, probs = NULL);
    true.post.q  <- pbSmpBiom(true.x.cut$x.count, prior.p = prior.p, iter = iter);
    true.post.p  <- pbSmpResp(true.x.cut$x.cut,   true.y, cand.cuts = cand.cuts, type = resp.mdl);
    true.cumu.pq <- pbCumuPQ(post.q, post.p);

    browser();

    ## predict outcomes
    rst <- NULL;
    for (i in cand.cuts) {
        cur.estt <- c(0, cut.quants)[i];

        for (j in n2) {
            ## truth
            true.cur.pred <- pbCfPred(true.cumu.pq$cumu.q[,i,drop = F],
                                      true.cumu.pq$cumu.p[,i,drop = F],
                                      j, theta0, alpha, repeach);

            true.cur.sum <- apply(true.cur.pred, 2, mean);

            ## predict based on posterior
            cur.pred <- pbCfPred(cumu.pq$cumu.q[,i,drop = F],
                                 cumu.pq$cumu.p[,i,drop = F],
                                 j, theta0, alpha, repeach);

            cur.sum <- apply(cur.pred, 2, mean);

            ## get result
            cur.rst <- c(n1, nresp1, i, j,
                         cur.sum[2], cur.sum[1], cur.sum[3],cur.sum[2]/j,
                         true.cur.sum[2], true.cur.sum[1], true.cur.sum[3], true.cur.sum[2]/j);

            ## utility
            for (k in uti.f) {
                cur.uti      <- pbCfUti(cbind(cur.pred, n1, nresp1),
                                        utif = k, theta0 = theta0, estt = cur.estt,
                                        B1 = B1, C1 = C1, C2 = C2, C3 = C3);
                true.cur.uti <- pbCfUti(cbind(true.cur.pred, n1, nresp1),
                                        utif = k, theta0 = theta0, estt = cur.estt,
                                        B1 = B1, C1 = C1, C2 = C2, C3 = C3);
                rst <- rbind(rst,
                             c(cur.rst, k,
                               mean(cur.uti), mean(cur.uti > uti.cut),
                               mean(true.cur.uti), mean(true.cur.uti > uti.cut)));
            }
        }
    }

    ## return
    colnames(rst) <- c("N1", "N1Resp", "Cut", "N2",
                       "N2Resp", "NScr", "PrRej", "ETheta",
                       "TrN2Resp", "TrNScr", "TrPrRej", "TrETheta",
                       "UtiF", "MeanUti", "PrLgU", "TrMeanUti", "TrPrLgU");
    rownames(rst) <- NULL;

    rst
}

