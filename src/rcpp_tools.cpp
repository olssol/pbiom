
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]




// [[Rcpp::init]]
void my_package_init(DllInfo *dll) {
  // initialization code here
  R_useDynamicSymbols(dll, TRUE);
}


//' Final Hypothesis testing
//'
//' @param y observed binary outcome data
//' @param theta0 reference theta for H0
//'
//' @export
// [[Rcpp::export]]
double pbCfYtest(NumericVector y, double theta0) {
  int    n = y.size();
  double sumy = 0, phat, zscore, pval;

  for (int i = 0; i < n; i++) {
    sumy += y[i];
  }

  phat   = sumy/n;
  zscore = (phat - theta0)/sqrt(phat*(1-phat)/n);
  pval   = 1 - R::pnorm(zscore, 0.0, 1.0, 1, 0);
	return pval;
}

//' Predict for a single set of value for rep times
//'
//' @param rep number of replications
//' 
//'
//' @export
// [[Rcpp::export]]
void pbCfPredSingle(NumericVector vq, NumericVector vp, int n2,
                    double theta0, double alpha, int rep, NumericMatrix rst) {
  int           nc = vq.size(), rej;
  double        nscreen, sumvq = 0, pval;
  IntegerVector x(nc);
  NumericVector y(n2);
  double        cy, nresp;
  int           i, j, k, inx;

  // total probability of eligible biomarkers
  for (i = 0; i < nc; i++) {
    sumvq += vq[i];
  }

  for (i = 0; i < nc; i++) {
    vq[i] /= sumvq;
  }

  for (i = 0; i < rep; i++) {
    // negative binomial to get screened out numbers
    if (sumvq > 0.9999999999) {
      nscreen = 0;
    } else {
      nscreen = R::rnbinom(n2, sumvq);
    }

    // sample biomarkers
    R::rmultinom(n2, vq.begin(), nc, x.begin());

    // sample responses
    nresp = 0;
    inx   = 0;
    for (j = 0; j < nc; j++) {
      for (k = 0; k < x[j]; k++) {
        cy       = R::rbinom(1, vp[j]);
        y[inx]   = cy;
        nresp    += cy;
        inx++;
      }
    }

    //test pvalues
    pval = pbCfYtest(y, theta0);
    rej  = pval < alpha;

    //save
    rst(i, 0) = nscreen;
    rst(i, 1) = nresp;
    rst(i, 2) = rej;
    rst(i, 3) = n2;
    //rst(i, 5) = nresp / n2;
    //rst(i, 6) = theta0;
  }
}

//' Predict for posterior samples from the interim analysis
//'
//' 
//' 
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pbCfPred(NumericMatrix vq, NumericMatrix vp, IntegerVector n2,
                       double theta0, double alpha, int repeach) {
  int           nr = vq.nrow(), ncr = 4, nn = n2.size();
  NumericMatrix rst(repeach * nr * nn, ncr), crst(repeach, ncr);
  int           i, j, k, l;

  for (l = 0; l < nn; l++) {
    for (i = 0; i < nr; i++) {
      // predict single 
      pbCfPredSingle(vq(i,_), vp(i,_), n2(l), theta0, alpha, repeach, crst);

      for (j = 0; j < repeach; j++) {
        for (k = 0; k < ncr; k++) {
          rst(l*(repeach*nr) + i*repeach + j, k) = crst(j, k);
        }
      }
    }
  }

  colnames(rst) = CharacterVector::create("nscreenout", "nresp", "rej", "n2");
  return rst;
}

//' Get utilities
//'
//'
//' @export
// [[Rcpp::export]]
NumericVector pbCfUti(NumericMatrix prst, int utif, double theta0, double estt,
                      double B1, double C1, double C2, double C3) {
  int           nr = prst.nrow();;
  NumericVector rst(nr);
  double        etheta, benefit, cost, extra, rej;
  int           i;

  for (i = 0; i < nr; i++) {
    etheta = prst(i,1) / prst(i,3);
    extra  = prst(i,0) / prst(i,3);
    rej    = prst(i,2);
    switch (utif)
      {
      case 1:
        benefit = B1 * (etheta - theta0) * rej;
        cost    = C1 * extra;
        break;
      case 2:
        benefit = B1 * (etheta - theta0) * (1 - estt) * rej;
        cost    = C1 * extra;
        break;
      default:
        benefit = rej;
        cost    = 0;
      }
    rst(i) = benefit - cost;
  }

  return rst;
}

//' Test function
//'
//' @param u input value
//' @return u squared 
//' @export
// [[Rcpp::export]]
double ptemp(double u) {
  double rst;
  rst = pow(u, 2);
  return(rst);
}
