#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// next from https://yiransblog.com/2019/05/08/article1/

// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
    arma::mat C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
double fastMean(const arma::vec X){
	double out = arma::mean(X);
	return(out);
}

// [[Rcpp::export]]
arma::vec fastcolMeans(const arma::mat & X){
	arma::vec out = arma::mean(X,0); //Armadillo mean(X,dim)
	return(out);
}

