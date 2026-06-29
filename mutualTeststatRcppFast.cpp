// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <numeric>
using namespace Rcpp;

inline double thres(double x) { return (x == 0.0) ? 1.0 : x; }

// unique sorted copy
static inline std::vector<double> unique_sorted(const NumericVector& a) {
  std::vector<double> u(a.begin(), a.end());
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  return u;
}

// map values to group indices in uniq (0..m-1)
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
std::vector<double> mutualTeststatRcpp_sorted3D(const NumericMatrix& Dx,
                                                const NumericMatrix& Dy,
                                                const NumericMatrix& Dz) {
  const int n = Dx.nrow();
  const double invn = 1.0 / n;
  
  // ---- Global accumulators (NO per-i normalization) ----
  // Diagonal (j==k==l)
  double jt_stat_sum = 0.0, jt_AD_sum = 0.0, jt_ADsum_sum = 0.0, jt_W_sum = 0.0;
  // Product (all j,k,l)
  double prod_stat_sum = 0.0, prod_AD_sum = 0.0, prod_ADsum_sum = 0.0, prod_W_sum = 0.0;
  
  // reusable buffers
  std::vector<int> gu(n), gv(n), gw(n);
  
  for (int i = 0; i < n; ++i) {
    NumericVector u = Dx(i, _);
    NumericVector v = Dy(i, _);
    NumericVector w = Dz(i, _);
    
    // tie groups for u,v,w (closed balls)
    std::vector<double> su = unique_sorted(u);
    std::vector<double> sv = unique_sorted(v);
    std::vector<double> sw = unique_sorted(w);
    const int nu = (int)su.size();
    const int nv = (int)sv.size();
    const int nw = (int)sw.size();
    
    // map samples to groups
    map_to_groups(u, su, gu);
    map_to_groups(v, sv, gv);
    map_to_groups(w, sw, gw);
    
    // group counts (exact group) and ≤-prefix (for marginals)
    std::vector<int> cu(nu, 0), cv(nv, 0), cw(nw, 0);
    for (int t = 0; t < n; ++t) { ++cu[gu[t]]; ++cv[gv[t]]; ++cw[gw[t]]; }
    
    std::vector<int> pcu = cu, pcv = cv, pcw = cw; // make prefixes
    std::partial_sum(pcu.begin(), pcu.end(), pcu.begin());
    std::partial_sum(pcv.begin(), pcv.end(), pcv.begin());
    std::partial_sum(pcw.begin(), pcw.end(), pcw.begin());
    
    // Fx,Fy by group + Ax,Ay
    std::vector<double> Fx_u(nu), Fy_v(nv), Ax_u(nu), Ay_v(nv);
    for (int r = 0; r < nu; ++r) { Fx_u[r] = pcu[r] * invn; Ax_u[r] = Fx_u[r] * (1.0 - Fx_u[r]); }
    for (int c = 0; c < nv; ++c) { Fy_v[c] = pcv[c] * invn; Ay_v[c] = Fy_v[c] * (1.0 - Fy_v[c]); }
    
    // bucket indices by w-group
    std::vector<std::vector<int>> bucket_w(nw);
    for (int t = 0; t < n; ++t) bucket_w[gw[t]].push_back(t);
    
    // cumulative 2D table S(r,c) and its 2D prefix P2D for XY up to current d
    std::vector<int> S((size_t)nu * nv, 0);
    std::vector<int> planeD((size_t)nu * nv, 0); // counts at EXACT current d for diagonal
    std::vector<int> P2D((size_t)nu * nv, 0);
    
    for (int d = 0; d < nw; ++d) {
      // add this w-group into S; track exact plane for diagonal weights
      std::fill(planeD.begin(), planeD.end(), 0);
      for (int idx : bucket_w[d]) {
        const int r = gu[idx], c = gv[idx];
        S[(size_t)r * nv + c] += 1;
        planeD[(size_t)r * nv + c] += 1;
      }
      
      // 2D prefix of S → P2D
      for (int r = 0; r < nu; ++r) {
        int acc = 0;
        int* prow = &P2D[(size_t)r * nv];
        const int* srow = &S[(size_t)r * nv];
        if (r == 0) {
          for (int c = 0; c < nv; ++c) { acc += srow[c]; prow[c] = acc; }
        } else {
          const int* prow_up = &P2D[(size_t)(r - 1) * nv];
          for (int c = 0; c < nv; ++c) { acc += srow[c]; prow[c] = prow_up[c] + acc; }
        }
      }
      
      const double Fz_d  = pcw[d] * invn;
      const double Az_d  = Fz_d * (1.0 - Fz_d);
      const int    cw_d  = cw[d];               // exact group size at d
      
      for (int r = 0; r < nu; ++r) {
        const double Fx  = Fx_u[r];
        const double Ax  = Ax_u[r];
        const double wr  = (double)cu[r];       // #j in group r
        const int*   prow = &P2D[(size_t)r * nv];
        
        for (int c = 0; c < nv; ++c) {
          const double Fy  = Fy_v[c];
          const double Ay  = Ay_v[c];
          const double wc  = (double)cv[c];     // #k in group c
          
          const double Fxyz = prow[c] * invn;   // ≤ in all three
          const double diff = Fxyz - Fx * Fy * Fz_d;
          const double stat = diff * diff;
          
          // keep EXACTLY your previous weights (thres)
          const double denAD  = Ax * Ay * Az_d;
          const double denAS  = Ax + Ay + Az_d;
          const double denW   = Fx * Fy * Fz_d;
          
          const double w_all  = wr * wc * (double)cw_d;                 // all (j,k,l)
          const double w_diag = (double)planeD[(size_t)r * nv + c];     // j==k==l
          
          // product sums (all triplets)
          prod_stat_sum  += w_all * stat;
          prod_AD_sum    += w_all * (stat / thres(denAD));
          prod_ADsum_sum += w_all * (stat / thres(denAS));
          prod_W_sum     += w_all * (stat / thres(denW));
          
          // diagonal sums (j==k==l)
          jt_stat_sum    += w_diag * stat;
          jt_AD_sum      += w_diag * (stat / thres(denAD));
          jt_ADsum_sum   += w_diag * (stat / thres(denAS));
          jt_W_sum       += w_diag * (stat / thres(denW));
        }
      }
    }
  }
  
  // ---- Normalize ONCE at the end (kbcov-style counts) ----
  const double n2 = (double)n * (double)n;     // centers × radii (diagonal)
  const double n4 = n2 * n2;                   // centers × all triplets
  
  double jt_stat   = jt_stat_sum   / n2;
  double jt_AD     = jt_AD_sum     / n2;
  double jt_ADsum  = jt_ADsum_sum  / n2;
  double jt_W      = jt_W_sum      / n2;
  
  double prod_stat = prod_stat_sum / n4;
  double prod_AD   = prod_AD_sum   / n4;
  double prod_ADsum= prod_ADsum_sum/ n4;
  double prod_W    = prod_W_sum    / n4;
  
  return {jt_stat, jt_AD, jt_W, jt_ADsum, prod_stat, prod_AD, prod_W, prod_ADsum};
}
