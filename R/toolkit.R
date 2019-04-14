
#' Export results into a template file
#'
#' @param numbers    vector of results
#' @param template.f template file name
#' @param out.f      output file name
#' @param sub.str    pattern of string to be replaced
#'
#'
#' @export
#'
tkExpRst <- function(numbers, template.f,  out.f="rst.txt", sub.str="AA") {
    if (!file.exists(template.f)) {
        return;
    }
    ##read
    tpla <- readChar(template.f, file.info(template.f)$size);

    ##substitute
    for (i in 1:length(numbers)) {
        tpla <- sub(sub.str, numbers[i], tpla);
    }

    ##write out
    write(tpla, file=out.f);
}


#' Import objects in a list into a designated environment
#'
#' @param alist list of objects
#' @param dest.env designated environment
#'
#' @export
#'
tkMakeLocal <- function(alist, dest.env='.GlobalEnv') {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest.env);
    }
}


#' Call function by its name organized as a vector
#'
#' @param vec function names as a vector
#'
#' @export
#'
tkCallFun <- function(vec, ...) {
    eval(parse(text=paste("rst <- ",
                          paste(vec, collapse = ""),
                          "(...)",
                          sep="")
               )
         );
    rst
}


#' Set default values in a list
#'
#'
#' @param tar.lst target list
#' @param value.lst list with values
#' @param overwrite Whether overwrite tar.lst with value.lst
#'
#' @export
#'
tkListAddValue <- function(tar.lst, value.lst, overwrite = FALSE) {

    v.names <- names(value.lst);
    for (i in 1:length(value.lst)) {
        cur.name <- v.names[i];

        if (overwrite | is.null(tar.lst[[cur.name]])) {
            tar.lst[[cur.name]] <- value.lst[[i]];
        }
    }

    tar.lst
}


#' Random draw from Dirichlet
#'
#' @param n     number of samples
#' @param alpha dirichlet parameters
#'
#' @export
#'
rdirichlet <- function (n, alpha) {
    l  <- length(alpha);
    x  <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE);
    sm <- x %*% rep(1, l);

    return(x/as.vector(sm));
}


#' Call STAN models
#'
#'
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#'
#' @export
#'
tkSTAN <- function(lst.data, stan.mdl = NULL,
                   chains = 4, iter = 2000, warmup = 1000,
                   control = list(adapt_delta=0.95), ...) {

    if (is.null(stan.mdl))
        return(NULL);

    stan.rst <- rstan::sampling(stanmodels[[stan.mdl]],
                                data    = lst.data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                control = control,
                                ...);

    stan.rst;
}
