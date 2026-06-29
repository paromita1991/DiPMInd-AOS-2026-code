// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cstring>
using namespace Rcpp;

inline double thres(double x) { return (x == 0.0) ? 1.0 : x; }

// Build sorted unique values of a
static inline std::vector<double> unique_sorted(const NumericVector& a) {
  std::vector<double> u(a.begin(), a.end());
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  return u;
}

// Map each entry of a to its tie-group index (0..m-1), where groups are unique_sorted(a).
// We use lower_bound here (exact group), and for queries with "≤" we use upper_bound on the unique array.
static inline void map_to_groups(const NumericVector& a,
                                 const std::vector<double>& uniq,
                                 std::vector<int>& idx) {
  const int n = a.size();
  idx.resize(n);
  for (int t = 0; t < n; ++t) {
    idx[t] = (int)(std::lower_bound(uniq.begin(), uniq.end(), a[t]) - uniq.begin());
  }
}

// [[Rcpp::export]]
std::vector<double> teststatv2Rcpp_sorted2D(const NumericMatrix& Dx, const NumericMatrix& Dy) {
  const int n = Dx.nrow();
  const double invn = 1.0 / n;
  
  double jt_stat = 0.0, prod_stat = 0.0, jt_AD = 0.0, prod_AD = 0.0,
    jt_ADsum = 0.0, prod_ADsum = 0.0, jt_W = 0.0, prod_W = 0.0;
  
  // Reused buffers across i
  std::vector<int> grp_u(n), grp_v(n);
  std::vector<int> ru_idx(n), rv_idx(n);  // ≤-group indices for j,k (via upper_bound)
  std::vector<double> Fx(n), Fy(n);
  
  for (int i = 0; i < n; ++i) {
    NumericVector u = Dx(i, _);
    NumericVector v = Dy(i, _);
    
    // Unique sorted values (tie-groups)
    std::vector<double> su = unique_sorted(u);
    std::vector<double> sv = unique_sorted(v);
    const int nu = (int)su.size();
    const int nv = (int)sv.size();
    
    // Map each sample t to its group (exact equality group)
    map_to_groups(u, su, grp_u);
    map_to_groups(v, sv, grp_v);
    
    // Build contingency H (nu x nv) and convert to 2D prefix sum P in-place.
    // Use 32-bit ints (n fits); flatten as row-major [ru*nv + rv].
    std::vector<int> P((size_t)nu * nv, 0);
    for (int t = 0; t < n; ++t) {
      P[(size_t)grp_u[t] * nv + grp_v[t]] += 1;
    }
    // 2D prefix sum: P[r,c] = sum_{r'<=r,c'<=c} H[r',c']
    for (int r = 0; r < nu; ++r) {
      int rowsum = 0;
      for (int c = 0; c < nv; ++c) {
        rowsum += P[(size_t)r * nv + c];
        if (r > 0) P[(size_t)r * nv + c] = P[(size_t)(r - 1) * nv + c] + rowsum;
        else       P[(size_t)r * nv + c] = rowsum;
      }
    }
    
    // Precompute ≤-indices for every j,k using upper_bound on unique arrays
    for (int j = 0; j < n; ++j) {
      // index of last group with su[ru] <= u[j]
      int ru = (int)(std::upper_bound(su.begin(), su.end(), u[j]) - su.begin()) - 1;
      ru_idx[j] = ru; // 0..nu-1
    }
    for (int k = 0; k < n; ++k) {
      int rv = (int)(std::upper_bound(sv.begin(), sv.end(), v[k]) - sv.begin()) - 1;
      rv_idx[k] = rv; // 0..nv-1
    }
    
    // Marginals via the 2D prefix sums:
    // Fx(j) = P[ru(j), nv-1]/n ; Fy(k) = P[nu-1, rv(k)]/n
    for (int j = 0; j < n; ++j) {
      Fx[j] = P[(size_t)ru_idx[j] * nv + (nv - 1)] * invn;
    }
    for (int k = 0; k < n; ++k) {
      Fy[k] = P[(size_t)(nu - 1) * nv + rv_idx[k]] * invn;
    }
    
    // Accumulate statistics over all (j,k)
    double sum_prod = 0.0, sum_AD = 0.0, sum_ADsum = 0.0, sum_W = 0.0;
    double sum_diag_prod = 0.0, sum_diag_AD = 0.0, sum_diag_ADsum = 0.0, sum_diag_W = 0.0;
    
    for (int j = 0; j < n; ++j) {
      const double Fxj = Fx[j];
      const double axj = Fxj * (1.0 - Fxj);
      const int    ru  = ru_idx[j];
      
      for (int k = 0; k < n; ++k) {
        const double Fyk = Fy[k];
        const double ayk = Fyk * (1.0 - Fyk);
        const int    rv  = rv_idx[k];
        
        // Fxy via 2D prefix sum
        const double Fxy = P[(size_t)ru * nv + rv] * invn;
        
        const double diff  = Fxy - Fxj * Fyk;
        const double stat  = diff * diff;
        const double denAD = axj * ayk;
        const double denAS = axj + ayk;
        const double denW  = Fxj * Fyk;
        
        sum_prod  += stat;
        sum_AD    += stat / thres(denAD);
        sum_ADsum += stat / thres(denAS);
        sum_W     += stat / thres(denW);
        
        if (j == k) {
          sum_diag_prod  += stat;
          sum_diag_AD    += stat / thres(denAD);
          sum_diag_ADsum += stat / thres(denAS);
          sum_diag_W     += stat / thres(denW);
        }
      }
    }
    
    // per-i normalization (same as your code)
    const double n2 = (double)n * n;
    prod_stat  += sum_prod / n2;
    jt_stat    += sum_diag_prod / n;
    prod_AD    += sum_AD / n2;
    jt_AD      += sum_diag_AD / n;
    prod_ADsum += sum_ADsum / n2;
    jt_ADsum   += sum_diag_ADsum / n;
    prod_W     += sum_W / n2;
    jt_W       += sum_diag_W / n;
  }
  
  return {jt_stat, jt_AD, jt_W, jt_ADsum, prod_stat, prod_AD, prod_W, prod_ADsum};
}
