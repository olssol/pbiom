#' Simulate a single trial
#'
#' @param par.biom list of parameters for simulating biomarkers
#' @param par.resp list of parameters for simulating responses
#' @param iter number of iterations for interim analysis
#' @param nlarge large n for evaluating the truth
#'
#' @export
#'
pbSimuSingleTrial <- function(par.biom, par.resp, n2,
                              prior.p = NULL, prior.q = c(a = 0.5, b = 0.5),
                              theta0,
                              cuts = NULL, probs = c(0.25, 0.5, 0.75),
                              cand.cuts = 1:(length(probs)+1),
                              resp.mdl = "simplebin",
                              uti.f, alpha = 0.05, uti.cut = 0.1,
                              B1 = 1, C1 = 1, C2 = 0, C3 = 0,
                              iter = 4000, nlarge = 10000, repeach = 1) {
    ## find cut points
    if (is.null(cuts)) {
        par.biom.truth   <- par.biom;
        par.biom.truth$n <- nlarge;
        truth.x          <- do.call(pbSimuBiom, par.biom);
        truth.cuts       <- quantile(truth.x, probs = probs);
    } else {
        truth.cuts <- cuts;
    }

    ## simulate first stage;
    s1x  <- do.call(pbSimuBiom, par.biom);
    s1y  <- do.call(pbSimuResp, c(list(x = s1x), par.resp))$y;

    ## interim analysis
    s1x.cut <- pbCutBiom(s1x, cuts = truth.cuts);
    post.q  <- pbSmpBiom(s1x.cut$x.count, prior.p = prior.p, iter = iter);
    post.p  <- pbSmpResp(s1x.cut$x.cut, s1y, cand.cuts = cand.cuts, type = resp.mdl);
    cumu.pq <- pbCumuPQ(post.q, post.p);

    ## predict outcomes
    rst <- NULL;
    for (i in cand.cuts) {
        cur.estt <- 1 - c(0, probs)[i];
        for (j in n2) {
            cur.pred <- pbCfPred(cumu.pq$cumu.q[,i,drop = F],
                                 cumu.pq$cumu.p[,i,drop = F],
                                 j, theta0, alpha, repeach);
            for (k in uti.f) {
                cur.uti <- pbCfUti(cur.pred, j, theta0 = theta0, estt = cur.estt,
                                   B1 = B1, C1 = C1, C2 = C2, C3 = C3);
                cur.rst <- c(i, j, k, mean(cur.uti), mean(cur.uti > uti.cut));
                rst     <- rbind(rst, cur.rst);
            }
        }
    }
    colnames(rst) <- c("Cut", "N2", "Uti", "Mean", "PLgU");
    rst
}

