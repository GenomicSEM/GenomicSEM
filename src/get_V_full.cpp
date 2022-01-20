#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.get_V_full)]]
arma::mat get_V_full(int k,
                        arma::mat V_LD,
                        double varSNPSE2,
                        arma::mat SE_SNP,
                        arma::mat I_LD,
                        arma::vec varSNP,
                        std::string GC,
                        arma::mat coords,
                        int i) {
    arma::mat V_SNP(k, k, arma::fill::zeros);
    //loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
    if (GC.compare(std::string("conserv")) == 0) {
        for (unsigned int p = 0; p < size(coords)[0]; p += 1) {
            int x = coords(p, 0) - 1;
            int y = coords(p, 1) - 1;
            if (x != y) {
                V_SNP(x, y) = SE_SNP(i - 1, y) * SE_SNP(i - 1, x) * I_LD(x, y) * I_LD(x, x) * I_LD(y, y) * (varSNP[i - 1] * varSNP[i - 1]);
            }
            if (x == y) {
                V_SNP(x, x) = (SE_SNP(i - 1, x) * I_LD(x, x) * varSNP[i - 1]) * (SE_SNP(i - 1, x) * I_LD(x, x) * varSNP[i - 1]);
            }
        }
    } else if (GC.compare(std::string("standard")) == 0) {
        for (unsigned int p = 0; p < size(coords)[0]; p += 1) {
            int x = coords(p, 0) - 1;
            int y = coords(p, 1) - 1;
            if (x != y) {
                V_SNP(x, y) = SE_SNP(i - 1, y) * SE_SNP(i - 1, x) * I_LD(x, y) * sqrt(I_LD(x, x)) * sqrt(I_LD(y, y)) * (varSNP[i - 1] * varSNP[i - 1]);
            }
            if (x == y) {
                V_SNP(x, x) = (SE_SNP(i - 1, x) * sqrt(I_LD(x, x)) * varSNP[i - 1]) * (SE_SNP(i - 1, x) * sqrt(I_LD(x, x)) * varSNP[i - 1]);
            }
        }
    } else if (GC.compare(std::string("none")) == 0) {
        for (unsigned int p = 0; p < size(coords)[0]; p += 1) {
            int x = coords(p, 0) - 1;
            int y = coords(p, 1) - 1;
            if (x != y) {
                V_SNP(x, y) = SE_SNP(i - 1, y) * SE_SNP(i - 1, x) * I_LD(x, y) * (varSNP[i - 1] * varSNP[i - 1]);
            }
            if (x == y) {
                V_SNP(x, x) = (SE_SNP(i - 1, x) * varSNP[i - 1]) * (SE_SNP(i - 1, x) * varSNP[i - 1]);
            }
        }
    }
    arma::vec v(((k+1)*(k+2))/2, arma::fill::ones);
    arma::mat V_full = diagmat(v);

    //input the ld-score regression region of sampling covariance from ld-score regression SEs
    V_full.submat(k+1, k+1, arma::size(V_full)[0]-1, arma::size(V_full)[0]-1) = V_LD;

    //add in SE of SNP variance as first observation in sampling covariance matrix
    V_full(0,0) = varSNPSE2;

    //add in SNP region of sampling covariance matrix
    V_full.submat(1, 1, k, k) = V_SNP;
    return V_full;
}