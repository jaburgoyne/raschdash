functions {
  // Log likelihood under a Rasch partial-credit model
  real pcm_lpmf(int y, vector beta_raw) {
    int K = num_elements(beta_raw);
    vector[K+1] beta = append_row(0, cumulative_sum(beta_raw));
    return categorical_logit_lpmf(y + 1 | beta);
  }

  // Log likelihood under a Rasch partial-credit model, dropping constant
  // additive terms (when rstan updates to latest version)
  // real pcm_lupmf(int y, vector beta_raw) {
  //   int K = num_elements(beta_raw);
  //   vector[K+1] beta = append_row(0, cumulative_sum(beta_raw));
  //   return categorical_logit_lupmf(y + 1 | beta);
  // }

  // Random responses under a Rasch partial-credit model
  int pcm_rng(vector beta_raw) {
    int K = num_elements(beta_raw);
    vector[K+1] beta = append_row(0, cumulative_sum(beta_raw));
    return categorical_logit_rng(beta) - 1;
  }
}

data {
  int<lower=1> L;               // number of testlets
  int<lower=1> I;               // number of items
  // int<lower=0> J;               // number of facets
  int<lower=1> K;               // max score on rating scale
  int<lower=1> M;               // number of groups
  int<lower=1> N;               // number of persons
  int<lower=1> O;               // number of observations
  int<lower=1,upper=M> mm[N];   // group for person n
  int<lower=-M,upper=N> nn[O];  // person for observation m
                                //   negative value => group observation
  int<lower=1,upper=L> ll[I];   // testlet for item i
  int<lower=1,upper=I> ii[O];   // item for observation o
  int<lower=1> kk[I];           // maximum score for item i
                                //   Rasch if 1
                                //   rating scale/partial credit if K (and > 1)
                                //   binomial otherwise
                                //   TODO: Add more flexible typing.
  int<lower=0> y[O];            // score for observation o
}

parameters {
  // The parameters are mostly in LISREL notation, with some abuse
  // of κ, δ, ε, and θ for items and testlets.
  real kappa;                      // intercept
  vector<offset=mu,[L] epsilon_raw;           // relative difficulty for testlet l
  vector[I] upsilon_raw;           // relative difficulty for item i
  real d_tau_raw;                  // spacing of last two thresholds.
  vector[max(0, K-2)] d_tau_err;   // relative spacing for threshold t.
                                   //   Vanishes if no rating scales
                                   //   or less than two thresholds.
  vector[M] xi_raw;                // ability for group j
  vector[N] zeta_raw;              // relative ability for person n
  real<lower=0> theta_epsilon_raw; // scale of testlet difficulties
  real<lower=0> theta_upsilon;     // scale of relative item difficulties
  real<lower=0> psi_raw;           // scale of group abilities
  real<lower=0> phi;               // scale of relative person abilities
}

transformed parameters {
  vector[L] epsilon;
  vector[I] delta;
  vector[K] tau;
  vector[M] xi;
  vector[N] eta;
  real theta_epsilon = L == 1 ? 0 : theta_epsilon_raw;
  real psi = M == 1 ? 0 : psi_raw;
  epsilon = nu + theta_epsilon * epsilon_raw;
  delta = epsilon[ll] + theta_upsilon * upsilon_raw;
  if (K == 1) {
    tau = [0]';
  } else if (K > 1) {
    // Our prior is that thresholds should be ascending,
    // but disordered thresholds are possible. Use the
    // traditional log(2) separation requirement
    // (Linacre, 2002) as a prior mean.
    vector[K-1] d_tau =
      log(2) + d_tau_raw + append_row(d_tau_err, 0);
    vector[K] tau_raw = cumulative_sum(append_row(0, d_tau));
    // Centre between the last two thresholds.
    tau = tau_raw - tau_raw[K] + 0.5 * d_tau_raw;
  }
  xi = psi * xi_raw;
  eta = xi[mm] + phi * zeta_raw;
}

model {
  // Gelman's standard normal prior everywhere is sensible here
  // and runs fast; testing for failure with normal(0, 100) has
  // essentially the same result. The Nuffic distribution
  // is closer to t_10, which allows for some outliers, but it runs
  // slower and fits equivalently to the standard normal.
  nu                ~ std_normal();
  epsilon_raw       ~ std_normal();
  upsilon_raw       ~ std_normal();
  d_tau_raw         ~ std_normal();
  d_tau_err         ~ std_normal();
  xi_raw            ~ std_normal();
  zeta_raw          ~ std_normal();
  theta_epsilon_raw ~ std_normal();
  theta_upsilon     ~ std_normal();
  psi_raw           ~ std_normal();
  phi               ~ std_normal();
  for (o in 1:O) {
    int i = ii[o];
    int k = kk[i];
    int n = nn[o];
    real beta = (n < 0 ? xi[-n] : eta[n]) - delta[i];
    // When K == 1, using the rsm function is correct but
    // inefficient. Checking for k == 1 first routes these items
    // to the built-in bernoulli_logit().
    if (k == 1)      y[o] ~ bernoulli_logit(beta);
    else if (k == K) y[o] ~ rsm(beta - tau);
    else             y[o] ~ binomial_logit(k, beta);
  }
}

generated quantities {
  real<lower=0> sigma = sqrt(square(psi) + square(phi));
  real<lower=0> theta =
    sqrt(square(theta_epsilon) + square(theta_upsilon));
  vector[L] testlet_difficulty = epsilon / sigma;
  vector[I] item_difficulty = delta / sigma;
  vector[K] thresholds = tau / sigma;
  vector[M] group_ability = xi / sigma;
  vector[N] person_ability = eta / sigma;
  real prior_testlet_difficulty = normal_rng(nu, theta_epsilon / sigma);
  real prior_item_difficulty =
    normal_rng(prior_testlet_difficulty, theta_upsilon / sigma);
  real prior_group_ability = normal_rng(0, psi / sigma);
  real prior_person_ability =
    normal_rng(prior_group_ability, phi / sigma);
  int<lower=0> y_rep[O];
  vector[O] log_lik;
  vector[O] log_lik_rep;
  for (o in 1:O) {
    int i = ii[o];
    int k = kk[i];
    int n = nn[o];
    real beta = (n < 0 ? xi[-n] : eta[n]) - delta[i];
    if (k == 1) {
      y_rep[o] = bernoulli_logit_rng(beta);
      log_lik[o] = bernoulli_logit_lpmf(y[o] | beta);
      log_lik_rep[o] = bernoulli_logit_lpmf(y_rep[o] | beta);
    } else if (k == K) {
      y_rep[o] = rsm_rng(beta - tau);
      log_lik[o] = rsm_lpmf(y[o] | beta - tau);
      log_lik_rep[o] = rsm_lpmf(y_rep[o] | beta - tau);
    } else {
      y_rep[o] = binomial_rng(k, inv_logit(beta));
      log_lik[o] = binomial_logit_lpmf(y[o] | k, beta);
      log_lik_rep[o] =
        binomial_logit_lpmf(y_rep[o] | k, beta);
    }
  }
}
