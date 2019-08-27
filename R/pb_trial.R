#' Simulate a single trial
#'
#' @param par.biom list of parameters for simulating biomarkers
#' @param par.resp list of parameters for simulating responses
#' @param iter number of iterations for interim analysis
#' @param nlarge large n for evaluating the truth
#' @param true.cumu.pq truth of cumu.pq
#'
#' @return A dataframe with the following columns
#' \itemize{
#'    \item{N1All}{Total number patients in stage 1}
#'    \item{N1RespAll}{Total number responders in stage 1}
#'    \item{N1}{Number of stage 1 patients satisfying the biomarker cut criteria}
#'    \item{N1Resp}{Number of stage 1 responders satisfying the biomarker cut criteria}
#' }
#'
#' @export
#'
pbSimuSingleTrial <- function(par.biom, par.resp, n2, theta0,
                              prior.p = NULL, prior.q = c(a = 0.5, b = 0.5),
                              cut.x = NULL, cut.quants = c(0.25, 0.5, 0.75), cut.known = TRUE,
                              cand.cuts = 1:(length(cut.quants)+1), true.cumu.pq = NULL,
                              resp.mdl = "simplebin",
                              uti.f, alpha = 0.05, uti.cut = 0.1,
                              B1 = 1, C1 = 1, C2 = 0, C3 = 0,
                              iter = 4000, nlarge = 50000,
                              repeach = 1) {

    ## utility columns
    uti.cols <- list(c("nscreenout", "n2", "nresp2", "rej2",  "etheta2"),
                     c("nscreenout", "n2", "nresp2", "rej12", "etheta12"));

    ## find cut points
    if (is.null(cut.quants)) {
        stop("Biomarker level cut quantiles (cut.quants) needs to be provided.");
    }

    ## simulate first stage;
    s1x        <- do.call(pbSimuBiom, par.biom);
    s1y        <- do.call(pbSimuResp, c(list(x = s1x), par.resp))$y;
    n1.all     <- par.biom$n;
    nresp1.all <- sum(s1y);

    ## cut points
    if (is.null(cut.x)) {
        if (cut.known) {
            true.par.biom    <- par.biom;
            true.par.biom$n  <- nlarge;
            true.x           <- do.call(pbSimuBiom, true.par.biom);
            cut.x            <- quantile(true.x, probs = cut.quants);
        } else {
            cut.x <- quantile(s1x, probs = cut.quants);
        }
    }

    ## interim analysis
    s1x.cut <- pbCutBiom(s1x, cuts = cut.x, probs = cut.quants);

    if (is.null(true.cumu.pq)) {
        post.q  <- pbSmpBiom(s1x.cut$x.count, prior.p = prior.p, iter = iter);
        post.p  <- pbSmpResp(s1x.cut$x.cut, s1y, cand.cuts = cand.cuts, type = resp.mdl);
        cumu.pq <- pbCumuPQ(post.q, post.p);
    } else {
        cumu.pq <- true.cumu.pq;
    }

    ## predict outcomes
    rst <- NULL;
    for (i in cand.cuts) {
        ## cut proportion
        cur.estt <- c(0, cut.quants)[i];

        ## stage 1 pts satisfying cuts;
        s1y.cut <- s1y[which(s1x.cut$x.cut >= i)];

        ## predict stage 2
        for (j in n2) {
            ## predict based on posterior
            cur.pred <- pbCfPred(cumu.pq$cumu.q[,i,drop = F],
                                 cumu.pq$cumu.p[,i,drop = F],
                                 j, theta0, alpha, repeach, s1y.cut);
            cur.sum <- apply(cur.pred, 2, mean);

            ## get result
            cur.rst <- c(n1.all, nresp1.all, i, cur.sum);

            ## utility
            for (k in uti.f) {
                cur.uti <- k;
                for (l in 1:2) {
                    tmp <- pbCfUti(cbind(cur.pred[,uti.cols[[l]], drop = FALSE],
                                         n1.all, nresp1.all),
                                   utif = k,
                                   theta0 = theta0, estt = cur.estt,
                                   B1 = B1, C1 = C1, C2 = C2, C3 = C3);

                    cur.uti <- c(cur.uti, mean(tmp), mean(tmp > uti.cut));
                }

                ## use stage 1&2
                rst <- rbind(rst, c(cur.rst,cur.uti));
            }
        }
    }

    ## return
    colnames(rst) <- c("N1All", "N1RespAll", "Cut", "NScr",
                       "N2", "N2Resp", "ETheta2",  "PrRej2",
                       "N1", "N1Resp", "ETheta12", "PrRej12",
                       "UtiF",
                       "MeanUti2", "PrLgU2", "MeanUti12", "PrLgU12");
    rownames(rst) <- NULL;

    list(cumu.pq  = cumu.pq,
         simu.rst = rst);
}

