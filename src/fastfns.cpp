#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>
#include "fastfns.h"
using namespace Rcpp;
using namespace arma;
// using std::vector;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double fastMean(const arma::vec & X){
    double out = arma::mean(X);
    return(out);
}

// [[Rcpp::export]]
double fastVar(const arma::vec & X){
    double out = arma::var(X);
    return(out);
}

// [[Rcpp::export]]
arma::vec fastcolMeans(const arma::mat & X){
    arma::vec out = arma::mean(X,0); //Armadillo mean(X,dim)
    return(out);
}

// [[Rcpp::export]]
arma::vec fastrowMeans(const arma::mat & X){
    //vec out = mean(X,1); //Armadillo mean(X,dim)
    //return(out);
    return mean(X,1);
}

// [[Rcpp::export]]
arma::vec fastrowMads(const arma::mat & X){
    int n = X.n_rows;
    vec res(n);
    // this next is the scale factor imagining normal data
    // https://en.wikipedia.org/wiki/Median_absolute_deviation#:~:text=In%20statistics%2C%20the%20median%20absolute,MAD%20calculated%20from%20a%20sample.
    // using it for unit testing. Shouldn't matter.
    const double aconstant = 1.4826;
    for(int i=0;i<n;++i){
        res(i) = median( abs( X.row(i) - median( X.row(i) ) ) );
    }
    return(res*aconstant) ;
    // arma::mat meddev = abs(X.each_col() - rowmeds);
    //return median(meddev,1);
}

// [[Rcpp::export]]
arma::vec fastrowMads2(const arma::mat & X){
    // this next is the scale factor imagining normal data
    // https://en.wikipedia.org/wiki/Median_absolute_deviation#:~:text=In%20statistics%2C%20the%20median%20absolute,MAD%20calculated%20from%20a%20sample.
    // using it for unit testing. Shouldn't matter.
    //const double aconstant = 1.4826;
    vec meddev = 1.4826 * median( abs(X.each_col() - median(X,1)), 1);
    return meddev;
}

// [[Rcpp::export]]
arma::vec fastrowMaxs(const arma::mat & X){
    return max(X,1);
}

// [[Rcpp::export]]
SEXP fastrowMaxs2(SEXP x){
    // copied from https://github.com/RfastOfficial/Rfast/blob/master/src/col_row_utilities.cpp
    int ncol=Rf_ncols(x),nrow=Rf_nrows(x);
    SEXP F;
    F=PROTECT(Rf_allocVector(REALSXP,nrow));
    double *xx=REAL(x),*end=xx+ncol*nrow,*f=REAL(F),*x3,*ff;
    const double *endf=f+LENGTH(F);
    for(ff=f;ff!=endf;++ff,++xx)
        *ff=*xx;
    for(;xx!=end;)
        for(ff=f,x3=xx,xx+=nrow;x3!=xx;++ff,++x3){
            *ff=std::max(*ff,*x3);
        }
    UNPROTECT(1);
    return F;
}


//[[Rcpp::export]]
arma::mat fastcova(const arma::mat & X){
    int n;
    rowvec m;
    mat s;
    n = X.n_rows;
    m = arma::mean(X,0) * sqrt(n);
    // arma::vec m = sqrt(n) * fastcolMeans(x);
    s =  ( X.t() * X - m.t() * m ) / (n - 1);
    // arma::mat s = (arma::cross(x,x) - arma::dot(m,m))/(n - 1);
    return s;
}


//[[Rcpp::export]]
Rcpp::NumericVector zscore_vec(const arma::vec & X){
    // mahalanobis distance between units in a 1d vector is a z-score
    // https://en.wikipedia.org/wiki/Mahalanobis_distance#Definition
    // https://codereview.stackexchange.com/questions/9554/implementation-of-mahalanobis-distances-using-eigen
    //https://learn.stat.ubc.ca/~andy.leung/files/seminars-talks/2012/03/RcppDemo.pdf
    //https://github.com/RfastOfficial/Rfast/blob/master/src/maha.cpp
    // This is just a z-score in the end.
    arma::vec res1 =  square( X - mean(X))/var(X) ;
    Rcpp::NumericVector res2 = Rcpp::wrap(res1);
    res2.attr("dim") = R_NilValue;
    // why doesn't .as_col() or  vectorise() work?
    //std::vector<double> res2 = conv_to<vector<double>>::from(res1)
    return res2;
    // return res;
}

//[[Rcpp::export]]
arma::mat vecdist_arma(const arma::vec & x){
    int n=x.n_elem, i = 0, j = 0;
    double d;
    arma::mat out(n,n);
    for (i = 0; i < n - 1; i++){
        for (j = i + 1; j < n ; j ++){
            // d = sqrt(pow( (x(i)-x(j)),2.0));
           d = std::abs( x(i)-x(j) );
           out(j,i)=d;
           out(i,j)=d;
        }
    }
    return(out);
}


//[[Rcpp::export]]
arma::mat vecdist3(const arma::vec & A){
    /// This is just to compute the pairwise distances of a vector
      //Z is not used here. But in the API for distance functions.
      // Calc euclidean distance of x -- return a square matrix n x n with all distances
      // http://nonconditional.com/2014/04/on-the-trick-for-computing-the-squared-euclidian-distances-between-two-s  ets-of-vectors/
      // see also https://stackoverflow.com/questions/35273292/eigen-calculate-the-distance-matrix-between-2-sets-o  f-vectors
      // and https://github.com/RoyiAvital/Projects/tree/master/CalcDistanceMatrix
    // int n = A.nrow();
    // arma::mat A = arma::mat(x.begin(), n, 1, false);
    // int n = A.n_elem;

    // arma::colvec An =  sum(square(A),1);
    arma::colvec An =  square(A);

    arma::mat C = -2 * (A * A.t());
    C.each_col() += An;
    C.each_row() += An.t();

    return sqrt(C);
}


// from mn.h in Rfast void minus_c(double f[],double &,double *,int,int &);
  // and from mn.cpp
  //vecdist
void minus_c(double f[],double &x,double *y,int offset,int &len){
      double *ff=f;
      for(int i=0;i<len;++i,ff+=offset,++y){
          *ff=std::abs(x-*y);
      }
}

//[[Rcpp::export]]
Rcpp::NumericMatrix vecdist2(const Rcpp::NumericVector & x){
    int nrow,ncol;
    nrow=ncol=x.length();
    SEXP F=PROTECT(Rf_allocMatrix(REALSXP,ncol,nrow));
    double *xx=REAL(x),*end=xx+nrow,*f=REAL(F),*y=xx;
    for(;xx!=end;++xx,f+=nrow)
        minus_c(f,*xx,y,1,ncol);
    UNPROTECT(1);
      return F;
  }


// [[Rcpp::export]]
arma::vec avg_rank_arma(const arma::vec & x)
{
    R_xlen_t sz = x.n_elem;
    arma::vec w = linspace(0, sz-1, sz);
    std::sort(w.begin(), w.end(),
        [&x](std::size_t i, std::size_t j) { return x[i] < x[j]; });

    arma::vec r(sz);
    for ( R_xlen_t n, i = 0; i < sz; i += n) {
        n = 1;
        while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
        for ( R_xlen_t k = 0; k < n; k++) {
            r[w[i + k]] = i + (n + 1) / 2.;
        }
    }
    return r;
}

// see http://arma.sourceforge.net/docs.html#conv_to
typedef std::vector<double> stdvec;


template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
   // https://stackoverflow.com/questions/14253069/convert-rcpparmadillo-vector-to-rcpp-vector
    return Rcpp::NumericVector(x.begin(), x.end());
}

//[[Rcpp::export]]
Rcpp::List fast_dists_and_trans(const arma::vec & x, const arma::vec & Z){
    //https://stackoverflow.com/questions/24997510/how-to-speed-up-this-rcpp-function
    // Z is not used
    int n=x.n_elem;
    arma::mat dx(n , n), dxRank0(n , n);
    arma::vec mndx(n), mndrx(n), maddx(n), maddrx(n), maxdx(n), maxdrx(n), zx(n), rankx(n);

    dx = vecdist_arma(x);
    //  Rcpp::NumericMatrix dxTmp = vecdist2(Rcpp::NumericVector(x.begin(),x.end()));
    // dx = as<arma::mat>(dxTmp);
    // arma::mat dx(dxTmp.begin(), n, n, false);

    rankx = avg_rank_arma(x);
    dxRank0 = vecdist_arma(rankx);
    // Rcpp::NumericMatrix drxTmp = vecdist2(Rcpp::NumericVector(rankx.begin(),rankx.end()));
    // arma::mat dxRank0(drxTmp.begin(), n, n, false);

    mndx = mean(dx,1);
    mndrx = mean(dxRank0,1);
    maddx = fastrowMads2(dx);
    maddrx = fastrowMads2(dxRank0);
    maxdx = fastrowMaxs(dx);
    maxdrx = fastrowMaxs(dxRank0);
    zx = zscore_vec(x);
    // Order here matters:
//  outcome_names <- c(theresponse, "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY")
    List res = List::create(mndx, mndrx, maddx, maddrx, maxdx, maxdrx, zx, rankx);
    // List L = List::create(Named("name1") = v1 , _["name2"] = v2);
     return(res);
}

