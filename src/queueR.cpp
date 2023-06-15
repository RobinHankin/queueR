#include "queueR.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pochhammer(const NumericVector a_real, const NumericVector a_imag,
                               const NumericVector q_real, const NumericVector q_imag,
			       const NumericVector n) {
  const size_t r = a_real.size();
  NumericVector result_real(r);
  NumericVector result_imag(r);
  const std::complex<double> one (1,0);

  for (int i = 0; i < r; i++) {
          std::complex<double> out = one;
          std::complex<double> outold (0,0);
    const std::complex<double> a (a_real[i], q_imag[i]);
    const std::complex<double> q (q_real[i], q_imag[i]);
          std::complex<double> s (a_real[i], a_imag[i]);  // s == subt

	  unsigned int f=0;

    while(
	  ((std::real(outold) != std::real(out)) || (std::imag(outold) != std::imag(out))) && (f < n[i])
	  ){
      outold = out;
      out *= (one-s);
      s *= q;
      f++;
    }
    result_real[i] = std::real(out);
    result_imag[i] = std::imag(out);
  } // i loop closes
  return Rcpp::cbind(result_real, result_imag);
}
