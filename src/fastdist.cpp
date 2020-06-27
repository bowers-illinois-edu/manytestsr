#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <R.h>
#include <Rinternals.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppEigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

// next from https://yiransblog.com/2019/05/08/article1/

// [[Rcpp::export]]
SEXP armaMatMult(arma::mat & A, arma::mat & B){
    arma::mat C = A * B;
    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd & A, Eigen::MatrixXd & B){
    Eigen::MatrixXd C = A * B;
    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(Eigen::Map<Eigen::MatrixXd> & A, Eigen::Map<Eigen::MatrixXd> & B){
    Eigen::MatrixXd C = A * B;
    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
double fastMean(const arma::vec & X){
	double out = arma::mean(X);
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

//[[Rcpp::export]]
arma::vec row_means(NumericMatrix & x){
     // https://github.com/RfastOfficial/Rfast/blob/master/src/col_row_utilities.cpp
     mat X = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
     return mean(X, 1);
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


/////Working on implementing the whole distance creation step in C
// Mostly by copying code from https://github.com/RfastOfficial/Rfast
// with some small adaptations

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
NumericMatrix vecdist2(NumericVector & x){
    // double x=<as std::vector<double> >(x_);
    int nrow,ncol;
    nrow=ncol=LENGTH(x);
    SEXP F=PROTECT(Rf_allocMatrix(REALSXP,ncol,nrow));
    double *xx=REAL(x),*end=xx+nrow,*f=REAL(F),*y=xx;
    for(;xx!=end;++xx,f+=nrow)
        minus_c(f,*xx,y,1,ncol);
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
    //return(s);
   return(s);
}

// fastmahal{}
// https://codereview.stackexchange.com/questions/9554/implementation-of-mahalanobis-distances-using-eigen
//
//https://learn.stat.ubc.ca/~andy.leung/files/seminars-talks/2012/03/RcppDemo.pdf
//https://github.com/RfastOfficial/Rfast/blob/master/src/maha.cpp

//
// dists_and_trans <- function(x, Z) {
// # I bet we can speed this up by just doing it in cpp
// ## if(!is.numeric(x){x <- as.numeric(x)}
//     dx <- vecdist(x)
//     dxRank0 <- vecdist(Rank(x)) # distance among the ranks
//     mnx <- fastMean(x)
//     sigx <- cova(matrix(x, ncol = 1))
//     res <- list(
//             mndist = rowmeans(dx),
//             mndistRank0 = rowmeans(dxRank0),
//             maddist = rowMads(dx),
//             maddistRank0 = rowMads(dxRank0),
//             maxdist = rowMaxs(dx, value = TRUE),
//             maxdistRank0 = rowMaxs(dxRank0, value = TRUE),
//             mhdist = mahala(xmat, mu = mnx, sigma = sigx)
//     )
//     return(res)
// }
//
//


// see also https://stackoverflow.com/questions/35273292/eigen-calculate-the-distance-matrix-between-2-sets-of-vectors
// and https://github.com/RoyiAvital/Projects/tree/master/CalcDistanceMatrix


// [[Rcpp::export]]
SEXP eigenDist(Eigen::Map<Eigen::VectorXd> & X){
    /// This is just to compute the pairwise distances of a vector
    //Z is not used here. But in the API for distance functions.
    // Calc euclidean distance of x -- return a square matrix n x n with all distances
    // http://nonconditional.com/2014/04/on-the-trick-for-computing-the-squared-euclidian-distances-between-two-sets-of-vectors/

    const int N = X.rows();

    Eigen::MatrixXd XX(N,N);
    Eigen::MatrixXd YY(N,N);
    Eigen::MatrixXd XY(N,N);
    Eigen::MatrixXd D(N,N);

    XX = X.array().square().rowwise().sum();
    YY = XX.transpose();
    XY = (2*X)*X.transpose();
    D = XX * Eigen::MatrixXd::Ones(1,N);
    D = D + Eigen::MatrixXd::Ones(N,1) * YY;
    D = D - XY;
    Eigen::MatrixXd sqrtD = D.array().sqrt();

    return Rcpp::wrap(sqrtD);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix calcPWD1 (const Rcpp::NumericVector & x){
    // from https://stackoverflow.com/questions/36829700/rcpp-my-distance-matrix-program-is-slower-than-the-function-in-package
    unsigned int outrows = x.length(), i = 0, j = 0;
    double d;
    Rcpp::NumericMatrix out(outrows,outrows);

    for (i = 0; i < outrows - 1; i++){
        for (j = i + 1; j < outrows ; j ++){
           // d = sqrt(pow(x(i)-x(j), 2.0));
            d = sqrt(pow(x(i)-x(j),2.0));
            out(j,i)=d;
            out(i,j)=d;
        }
    }

    return out;
}

