#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

//' Compute image denoising.
//'
//' @param img the original image matrix
//' @param N the dimension of texture replacement search region
//' @param k the dimension of the texture region of the pixel to be processed
//' @param sigma the smoothing parameters
//' @return the denoised image
// [[Rcpp::export]]
arma::mat nlm(arma::mat img, int N, int K, double sigma){
  int m, n, pad_len;
  m = img.n_rows;
  n = img.n_cols;
  pad_len = N+K;
  arma::mat ssd;

  arma::mat arraypad(2 * pad_len + m, 2 * pad_len + n, fill::zeros);
  arraypad(span(pad_len, pad_len - 1 + m), span(pad_len, pad_len - 1 + n)) = img;
  arma::mat yy(m, n, fill::zeros);
  arma::mat B(m, n, fill::zeros);

  for(int nx = -N; nx < N; nx++)
    for(int ny = -N; ny < N; ny++){
      arma::mat ssd(m, n, fill::zeros);
      for(int kx = -K; kx < K; kx++)
        for(int ky = -K; ky < K; ky++){
          arma::mat s = arraypad(span(pad_len + ny + ky, m + pad_len + ny + ky - 1), span(pad_len + nx + kx, n + pad_len + nx + kx - 1))
          - arraypad(span(pad_len + ky, m + pad_len + ky - 1), span(pad_len + kx, n + pad_len + kx - 1));
          ssd = ssd + s % s;
          arma::mat ex = exp(- ssd / (2 * sigma * sigma));
          B = B + ex;
          yy = yy + ex % arraypad(span(pad_len + ny, m + pad_len + ny - 1), span(pad_len + nx, n + pad_len + nx - 1));
        }
    }
  return yy / B;
}
