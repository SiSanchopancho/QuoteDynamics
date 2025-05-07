// file: src/build_info.cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppEigen.h>
#include <string>

// [[Rcpp::export]]
Rcpp::List build_info() {

  // ----- Eigen -----------------------------------------------------------
  int eig_world = EIGEN_WORLD_VERSION;   // 3
  int eig_major = EIGEN_MAJOR_VERSION;   // 4
  int eig_minor = EIGEN_MINOR_VERSION;   // 0

  // ----- Compiler --------------------------------------------------------
  std::string cc, cc_ver;

#if defined(__clang__)
  cc = "Clang/LLVM";
  cc_ver = std::to_string(__clang_major__) + "." +
    std::to_string(__clang_minor__) + "." +
    std::to_string(__clang_patchlevel__);
#elif defined(__GNUC__)
  cc = "GCC";
  cc_ver = __VERSION__;                 // z. B. "13.2.1 20240210"
#elif defined(_MSC_VER)
  cc = "MSVC";
  cc_ver = std::to_string(_MSC_VER);    // z. B. 1936
#else
  cc = "Unknown";
  cc_ver = "n/a";
#endif

  // ----- Rückgabe an R ---------------------------------------------------
  return Rcpp::List::create(
    Rcpp::Named("Eigen") =
    std::to_string(eig_world) + "." +
    std::to_string(eig_major) + "." +
    std::to_string(eig_minor),
    Rcpp::Named("Compiler") = cc,
    Rcpp::Named("Version") = cc_ver
  );
}
