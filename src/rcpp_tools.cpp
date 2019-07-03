
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
                    double theta0, double alpha, int rep, IntegerVector s1y,
                    NumericMatrix rst) {

  int           nc = vq.size(), n1 = s1y.size(), rej2, rej12;
  double        nscreen = 0, sumvq = 0;
  IntegerVector x(nc);
  NumericVector y2(n2), y12(n1 + n2);
  double        nresp2 = 0, nresp1 = 0, cy;
  int           i, j, k, inx;

  // total probability of eligible biomarkers
  for (i = 0; i < nc; i++) {
    sumvq += vq[i];
  }

  for (i = 0; i < nc; i++) {
    vq[i] /= sumvq;
  }

  for (i = 0; i < rep; i++) {
    //stage 1 data
    for (j = 0; j < n1; j++) {
      y12[j] =  s1y[j];
      nresp1 += s1y[j];
    }

    if (n2 > 0) {
      // negative binomial to get screened out numbers
      if (sumvq > 0.9999999999) {
        nscreen = 0;
      } else {
        nscreen = R::rnbinom(n2, sumvq);
      }

      // sample biomarkers
      R::rmultinom(n2, vq.begin(), nc, x.begin());

      // sample responses
      inx = 0;
      for (j = 0; j < nc; j++) {
        for (k = 0; k < x[j]; k++) {
          cy       = R::rbinom(1, vp[j]);
          y2[inx]  = cy;
          nresp2   += cy;
          inx++;
        }
      }

      //test pvalues
      rej2  = pbCfYtest(y2, theta0) < alpha;

      //add stage 2 to stage 1
      for (j = 0; j < n2; j++) {
        y12[n1 + j] = y2[j];
      }
    }

    //test pvalues using stage 1 and 2 data
    rej12  = pbCfYtest(y12, theta0)  < alpha;

    if (0 == n2) {
      rej2 = rej12;
    }

    //save
    rst(i, 0) = nscreen;
    rst(i, 1) = n2;
    rst(i, 2) = nresp2;
    rst(i, 3) = nresp2 / n2;
    rst(i, 4) = rej2;
    rst(i, 5) = n1;
    rst(i, 6) = nresp1;
    rst(i, 7) = (nresp1 + nresp2)/(n1 + n2);
    rst(i, 8) = rej12;
  }
}

//' Predict for posterior samples from the interim analysis
//'
//' 
//' @param s1y sub-vector from stage 1 for patients satisfying biomarker cut point
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pbCfPred(NumericMatrix vq, NumericMatrix vp, IntegerVector n2,
                       double theta0, double alpha, int repeach, IntegerVector s1y) {

  int           nr = vq.nrow(), ncr = 9, nn = n2.size();
  NumericMatrix rst(repeach * nr * nn, ncr), crst(repeach, ncr);
  int           i, j, k, l;

  for (l = 0; l < nn; l++) {
    for (i = 0; i < nr; i++) {
      // predict single 
      pbCfPredSingle(vq(i,_), vp(i,_), n2(l), theta0, alpha, repeach, s1y, crst);

      for (j = 0; j < repeach; j++) {
        for (k = 0; k < ncr; k++) {
          rst(l*(repeach*nr) + i*repeach + j, k) = crst(j, k);
        }
      }
    }
  }

  colnames(rst) = CharacterVector::create("nscreenout",
                                          "n2", "nresp2", "etheta2", "rej2",
                                          "n1", "nresp1", "etheta12", "rej12");
  return rst;
}

//' Get utilities
//'
//' @param prst predicted results plus stage 1 total size and number of responders
//'
//' @export
// [[Rcpp::export]]
NumericVector pbCfUti(NumericMatrix prst, int utif, double theta0, double estt,
                      double B1, double C1, double C2, double C3) {
  int           nr = prst.nrow();;
  NumericVector rst(nr);
  double        n1, n2, nresp1, nresp2, nsout, rej;
  double        etheta, respall;
  double        benefit, cost;
  int           i;

  for (i = 0; i < nr; i++) {
    nsout   = prst(i,0);
    n2      = prst(i,1);
    nresp2  = prst(i,2);
    rej     = prst(i,3);
    etheta  = prst(i,4);
    n1      = prst(i,5);
    nresp1  = prst(i,6);

    respall = nresp1 + nresp2;

    switch (utif)
      {
      case 1:
        benefit = B1 * (n1 + n2) * (etheta - theta0) * rej;
        cost    = C2 * nsout + C1 * n1 + C3 * n2;
        break;
      case 2:
        benefit = B1 * (n1 + n2) * (etheta - theta0) * (1 - estt) * rej;
        cost    = C2 * nsout + C1 * n1 + C3 * n2;
        break;
      case 5:
        benefit = B1 * (n1 + n2) * rej + respall;
        cost    = (n1 + n2 - respall) + C2 * nsout + C1 * n1 + C3 * n2;
        break;
      case 7:
        benefit = B1 * (n1 + n2) * (etheta - theta0) * rej + respall;
        cost    = (n1 + n2 - respall) + C2 * nsout + C1 * n1 + C3 * n2;
        break;
      case 12:
        benefit = B1 * (n1 + n2) * rej;
        cost    = C2 * nsout + C1 * n1 + C3 * n2;
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
  Rcout << "value " << u << "\n";

  rst = pow(u, 2);
  return(rst);
}
