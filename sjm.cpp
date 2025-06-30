// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Soft‐thresholding: returns sign(x) * max(0, |x| − delta)
static arma::rowvec soft_threshold(const arma::rowvec &x, double delta) {
  arma::rowvec absx = arma::abs(x);
  arma::rowvec sub = absx - delta;
  sub.for_each([](arma::mat::elem_type &val) { if (val < 0) val = 0; });
  return arma::sign(x) % sub;
}

// Calculate new feature weights: normalize soft‐thresholded objective
static arma::rowvec calc_new_feature_weights(const arma::rowvec &objective, double delta) {
  arma::rowvec soft = soft_threshold(objective, delta);
  double norm_soft = arma::norm(soft, 2);
  if (norm_soft == 0.0) {
    int p = objective.n_elem;
    return arma::rowvec(p).fill(1.0 / std::sqrt((double)p));
  }
  return soft / norm_soft;
}

// Binary search to find delta so that ∥soft_threshold(objective, δ)∥₁/∥soft_threshold(objective, δ)∥₂ ≤ norm_constraint
static double binary_search(const arma::rowvec &objective, double norm_constraint, int max_iter = 15) {
  double l2n = arma::norm(objective, 2);
  if (l2n == 0.0 || arma::accu(arma::abs(objective) / l2n) <= norm_constraint) {
    return 0.0;
  }
  double lam1 = 0.0;
  double lam2 = arma::abs(objective).max() - 1e-5;
  double mid = 0.0;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    mid = 0.5 * (lam1 + lam2);
    arma::rowvec su = soft_threshold(objective, mid);
    double norm_su = arma::norm(su, 2);
    if (norm_su == 0.0) {
      lam2 = mid;
      continue;
    }
    double l1frac = arma::accu(arma::abs(su) / norm_su);
    if (l1frac < norm_constraint) {
      lam2 = mid;
    } else {
      lam1 = mid;
    }
    if (std::abs(lam2 - lam1) < 1e-4) {
      break;
    }
  }
  return 0.5 * (lam1 + lam2);
}

// Compute Between‐Cluster Sum of Squares (BCSS) per feature
// [[Rcpp::export]]
static arma::rowvec get_BCSS(const arma::mat &Y, const arma::uvec &states, int n_states) {
  int p = Y.n_cols;
  arma::rowvec global_mean = arma::mean(Y, 0);
  arma::mat diff = Y.each_row() - global_mean;
  arma::rowvec TSS = arma::sum(arma::square(diff), 0);
  
  arma::rowvec WCSS(p, arma::fill::zeros);
  for (int i = 0; i < n_states; ++i) {
    arma::uvec idx = arma::find(states == (unsigned)i);
    if (idx.n_elem > 0) {
      arma::mat Yc    = Y.rows(idx);
      arma::rowvec cm = arma::mean(Yc, 0);
      arma::mat d     = Yc.each_row() - cm;
      WCSS           += arma::sum(arma::square(d), 0);
    }
  }
  return TSS - WCSS;
}

// Wrapper to get feature weights given Y and a state sequence
static arma::rowvec get_weights(const arma::mat &Y, const arma::uvec &states, double max_features, int n_states) {
  arma::rowvec BCSS   = get_BCSS(Y, states, n_states);
  double      delta   = binary_search(BCSS, max_features);
  return calc_new_feature_weights(BCSS, delta);
}

// K‐means++ style initialization of states (0‐based labels)
static arma::uvec init_states(const arma::mat &Y, int n_states) {
  int n_obs = Y.n_rows;
  int p     = Y.n_cols;
  arma::mat centers(n_states, p, arma::fill::zeros);
  
  // 1) Choose first center at random
  int first_idx = (int)std::floor(R::runif(0.0, (double)n_obs));
  centers.row(0) = Y.row(first_idx);
  
  int              n_local_trials = 2 + (int)std::log((double)n_states);
  arma::vec        closest_dist_sq(n_obs);
  for (int i = 0; i < n_obs; ++i) {
    arma::rowvec diff = Y.row(i) - centers.row(0);
    closest_dist_sq[i] = arma::dot(diff, diff);
  }
  double current_pot = arma::accu(closest_dist_sq);
  
  for (int c = 1; c < n_states; ++c) {
    arma::vec rand_vals(n_local_trials);
    for (int t = 0; t < n_local_trials; ++t) {
      rand_vals[t] = R::runif(0.0, current_pot);
    }
    arma::vec cumsum = arma::cumsum(closest_dist_sq);
    
    std::vector<int> candidate_ids(n_local_trials);
    for (int t = 0; t < n_local_trials; ++t) {
      double rv = rand_vals[t];
      int    idx = (int)(std::lower_bound(cumsum.begin(), cumsum.end(), rv) - cumsum.begin());
      if (idx >= n_obs) idx = n_obs - 1;
      candidate_ids[t] = idx;
    }
    
    arma::mat best_dist_sq(n_local_trials, n_obs, arma::fill::zeros);
    for (int t = 0; t < n_local_trials; ++t) {
      arma::rowvec cand = Y.row(candidate_ids[t]);
      for (int i = 0; i < n_obs; ++i) {
        arma::rowvec diff = Y.row(i) - cand;
        best_dist_sq(t, i) = arma::dot(diff, diff);
      }
    }
    
    arma::vec       best_pot_vec(n_local_trials);
    arma::mat       new_closest(n_local_trials, n_obs, arma::fill::zeros);
    for (int t = 0; t < n_local_trials; ++t) {
      arma::rowvec new_dist = arma::min(closest_dist_sq.t(), best_dist_sq.row(t));
      new_closest.row(t)    = new_dist;
      best_pot_vec[t]       = arma::accu(new_dist);
    }
    arma::uword best_trial     = best_pot_vec.index_min();
    int         best_candidate = candidate_ids[best_trial];
    
    centers.row(c) = Y.row(best_candidate);
    current_pot    = best_pot_vec[best_trial];
    closest_dist_sq = new_closest.row(best_trial).t();
  }
  
  arma::uvec states(n_obs);
  for (int i = 0; i < n_obs; ++i) {
    arma::rowvec dists(n_states);
    for (int j = 0; j < n_states; ++j) {
      arma::rowvec diff = Y.row(i) - centers.row(j);
      dists[j]          = arma::dot(diff, diff);
    }
    states[i] = dists.index_min();
  }
  return states;
}

// The core "jump" algorithm: piecewise‐constant fit with n_states
// Returns pair<state‐vector, best_loss>
static std::pair<arma::uvec, double> jump_model(const arma::mat &Y,
                                                int               n_states,
                                                double            jump_penalty,
                                                int               max_iter,
                                                int               n_init,
                                                double            tol,
                                                bool              verbose,
                                                bool              has_initial,
                                                const arma::uvec &initial_states) {
  int n_obs = Y.n_rows;
  int p     = Y.n_cols;
  
  arma::mat Gamma(n_states, n_states, arma::fill::ones);
  Gamma -= arma::eye<arma::mat>(n_states, n_states);
  Gamma *= jump_penalty;
  
  arma::uvec best_s;
  double    best_loss = R_PosInf;
  
  arma::uvec s(n_obs);
  if (has_initial) {
    arma::uvec unique_init = arma::unique(initial_states);
    if ((int)unique_init.n_elem == n_states) {
      s = initial_states;
    } else {
      s = init_states(Y, n_states);
    }
  } else {
    s = init_states(Y, n_states);
  }
  
  for (int init = 0; init < n_init; ++init) {
    if (init > 0) {
      s = init_states(Y, n_states);
    }
    arma::mat mu(n_states, p, arma::fill::zeros);
    double   loss_old = 1e10;
    
    for (int it = 0; it < max_iter; ++it) {
      // 1) Update cluster means μ_i
      for (int i = 0; i < n_states; ++i) {
        arma::uvec idx = arma::find(s == (unsigned)i);
        if (idx.n_elem > 0) {
          arma::mat Yc    = Y.rows(idx);
          arma::rowvec cm = arma::mean(Yc, 0);
          mu.row(i)       = cm;
        }
      }
      
      // 2) Compute loss_by_state[t,i] = || Y[t] − μ_i ||²
      arma::mat loss_by_state(n_obs, n_states, arma::fill::zeros);
      for (int i = 0; i < n_states; ++i) {
        arma::rowvec center = mu.row(i);
        for (int t = 0; t < n_obs; ++t) {
          arma::rowvec diff = Y.row(t) - center;
          loss_by_state(t, i) = arma::dot(diff, diff);
        }
      }
      
      // 3) Dynamic programming backward pass to fill V
      arma::mat V = loss_by_state;
      for (int t = n_obs - 1; t > 0; --t) {
        for (int i = 0; i < n_states; ++i) {
          double min_val = R_PosInf;
          for (int j = 0; j < n_states; ++j) {
            double val = V(t, j) + Gamma(i, j);
            if (val < min_val) {
              min_val = val;
            }
          }
          V(t - 1, i) = loss_by_state(t - 1, i) + min_val;
        }
      }
      
      // 4) Reconstruct state sequence via argmin
      arma::uvec s_old = s;
      
      {
        arma::rowvec first_row = V.row(0);
        arma::uword  argmin_i  = first_row.index_min();
        s[0] = argmin_i;
      }
      for (int t = 1; t < n_obs; ++t) {
        int    prev   = s[t - 1];
        double best_v = R_PosInf;
        int    best_j = 0;
        for (int j = 0; j < n_states; ++j) {
          double val = Gamma(prev, j) + V(t, j);
          if (val < best_v) {
            best_v = val;
            best_j = j;
          }
        }
        s[t] = best_j;
      }
      
      {
        arma::uvec uniq_s = arma::unique(s);
        if ((int)uniq_s.n_elem == 1) {
          break;
        }
      }
      
      double loss = V.row(0).min();
      if (verbose) {
        Rcpp::Rcout << "    inner iter " << it << "  loss = " << loss << "\n";
      }
      
      if (tol > 0) {
        double epsilon = loss_old - loss;
        if (epsilon < tol) {
          break;
        }
      } else {
        if (arma::all(s == s_old)) {
          break;
        }
      }
      loss_old = loss;
    }
    
    if (loss_old < best_loss) {
      best_loss = loss_old;
      best_s    = s;
    }
  }
  
  return std::make_pair(best_s, best_loss);
}

// [[Rcpp::export]]
List sparse_jump(NumericMatrix Y_in,
                 int           n_states,
                 double        max_features,
                 double        jump_penalty  = 1e-5,
                 int           max_iter      = 10,
                 double        tol           = 1e-4,
                 int           n_init        = 10,
                 bool          verbose       = false) {
  arma::mat Y = as<arma::mat>(Y_in);
  int       p = Y.n_cols;
  
  // Clip max_features between [1, sqrt(p)]
  double upper = std::sqrt((double)p);
  if (max_features < 1.0) max_features = 1.0;
  if (max_features > upper) max_features = upper;
  
  // Initialize feature weights uniformly: 1/sqrt(p)
  arma::rowvec feat_w(p);
  feat_w.fill(1.0 / std::sqrt((double)p));
  
  arma::uvec states;      // will store 0-based labels
  bool      have_states = false;
  double    last_loss   = NA_REAL;
  
  for (int it = 0; it < max_iter; ++it) {
    // 1) Weight each column j by sqrt(feat_w[j])
    arma::rowvec sqw = arma::sqrt(feat_w);
    arma::mat   Yw  = Y;
    for (int j = 0; j < p; ++j) {
      Yw.col(j) *= sqw[j];
    }
    
    // 2) Run the jump_model on weighted data, capturing returned loss
    std::pair<arma::uvec,double> jl;
    if (!have_states) {
      jl           = jump_model(Yw, n_states, jump_penalty, max_iter, n_init, tol, verbose, false, states);
      have_states  = true;
    } else {
      jl           = jump_model(Yw, n_states, jump_penalty, max_iter, n_init, tol, verbose, true, states);
    }
    states    = jl.first;
    last_loss = jl.second;
    
    {
      arma::uvec uniq_s = arma::unique(states);
      if ((int)uniq_s.n_elem == 1) {
        break;
      }
    }
    
    // 3) Compute new feature weights on un‐weighted Y
    arma::rowvec new_w = get_weights(Y, states, max_features, n_states);
    
    double diff = arma::accu(arma::abs(new_w - feat_w)) / arma::accu(arma::abs(feat_w));
    if (diff < tol) {
      feat_w = new_w;
      break;
    }
    if (verbose) {
      Rcpp::Rcout << "  outer iter " << it << "  weight‐diff = " << diff << "\n";
    }
    feat_w = new_w;
  }
  
  // Compute weighted BCSS on Y using final feat_w
  arma::mat Y_weighted = Y;
  for (int j = 0; j < p; ++j) {
    Y_weighted.col(j) *= feat_w[j];
  }
  //for (int j = 0; j < p; ++j) {
    // ogni elemento di colonna j viene elevato a feat_w[j]
    //Y_weighted.col(j) = arma::pow(Y_weighted.col(j), feat_w[j]);
  //}
  arma::rowvec bcss_vec = get_BCSS(Y_weighted, states, n_states);
  double weighted_BCSS = arma::accu(bcss_vec);
  
  
  // Prepare R return
  int n_obs = Y.n_rows;
  IntegerVector R_states(n_obs);
  for (int i = 0; i < n_obs; ++i) {
    R_states[i] = (int)states[i] + 1; // convert to 1‐based
  }
  NumericVector R_feat_w(p);
  for (int j = 0; j < p; ++j) {
    R_feat_w[j] = feat_w[j];
  }
  
  return List::create(
    Named("states")  = R_states,
    Named("feat_w")  = R_feat_w,
    Named("loss") = weighted_BCSS
  );
}
