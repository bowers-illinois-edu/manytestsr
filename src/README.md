
Ideas for implementing the distance function in cpp: <https://stackoverflow.com/questions/3088650/how-do-i-create-a-list-of-vectors-in-rcpp>

<https://stackoverflow.com/questions/44508366/efficient-distance-calculations-in-armadillo>

```c
sp_mat arma_distmat_LT2(const arma::mat& x) { // input expected X_{n x p} n << p                               
  int nr, nc;
  Col<double> col0, col1;

  nr = x.n_rows;
  nc = x.n_cols;

  sp_mat out(nc, nc);
  out = trimatl(x.t() * x, k=-1);
  return out;
}

```

<http://nonconditional.com/2014/04/on-the-trick-for-computing-the-squared-euclidian-distances-between-two-sets-of-vectors/>


