#' Logit function
#'
#'
#'
#' @export
#'
logit  <- function(x) {
    log(x / (1-x));
}

#' Compute utility
#'
#'
#'
#' @description
#' theta0      <- ["theta0"];
#' resp.all    <- ["resp.all"];
#' resp2       <- ["resp2"];
#' r           <- ["Stage1prop"];
#' est.t       <- ["est.t"];
#' etheta      <- ["etheta"];
#' rej         <- ["rej"];
#' n.ScreenOut <- ["n.ScreenOut"];
#' n1          <- ["n1"]; 	n2   <- rst.summary["n2"];
#' p.extra     <- n.ScreenOut/(n1+n2)
#' p.resp.all  <- resp.all/(n1+n2)
#' p.resp2     <- resp2/(n1+n2)
#'
#' @export
#'
pbUtility  <- function(uti, ...) {
    par.lst <- list(...);
    par.lst <- tkListAddValue(par.lst,
                              list(B1 = 1, C1 = 1, C2 = 0, C3 = 0),
                              overwrite = FALSE);

    bc <- do.call(tkCallFun, c(list(vec = c("prvUti.", uti)), par.lst));

    c(utility = bc[1] - bc[2],
      benefit = bc[1],
      cost    = bc[2]);
}

#' Compute utility using a matrix of inputs
#'
#'
#'
#'
#' @export
#'
pbUtilityMat <- function(uti, mat, ...) {
    fu <- function(x) {
        rst <- NULL;
        for (i in uti) {
            cur.par <- c(list(uti = uti[i]), x, list(...));
            crst    <- do.call(pbUtility, cur.par);
            rst     <- c(rst, crst[1]);
        }
        rst
    }

    rst.uti           <- t(apply(mat, 1, fu));
    colnames(rst.uti) <- paste("utility", uti, sep ="");
    cbind(mat, rst.uti);
}

prvUti.1 <- function(n2, nscreenout, nresp, theta0, rej, B1, C1, ...) {
    etheta  <- nresp / n2;
    benefit <- B1 * (etheta - theta0) * rej;
    cost    <- C1 * nscreenout / n2;

    c(benefit, cost);
}

prvUti.2 <- function(n2, nscreenout, nresp, theta0, rej, estt, B1, C1, ...) {
    etheta  <- nresp / n2;
    benefit <- B1 * (etheta - theta0) * (1-estt) * rej;
    cost    <- C1 * nscreenout / n2;

    c(benefit, cost);
}

prvUti.5 <- function(rej, prespall, pextra, B1, C1, ...) {
    benefit <- B1 * rej + prespall;
    cost    <- (1 - prespall) + C1 * pextra;

    c(benefit, cost);
}

prvUti.10 <- function(rej, ...) {
    benefit <- rej;
    cost    <- 0;

    c(benefit, cost);
}


