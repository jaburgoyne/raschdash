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
  vector<lower=0>[O] w;         // weight of observation o
}

parameters {
  // The parameters are mostly in LISREL notation, with some abuse
  // of ν, δ, ε, and θ for items and testlets.
  real nu;                                // intercept
  vector[L > 1 ? L : 0] epsilon_raw;      // relative difficulty for testlet l
  vector[I] upsilon_raw;                  // relative difficulty for item i
  vector[K > 1 ? K : 0] tau_raw;          // relative difficulty of step t
  vector[M > 1 ? M : 0] xi_raw;           // ability for group j
  vector[N] zeta_raw;                     // relative ability for person n
  vector<lower=0>[L > 1 ? 1 : 0] theta_epsilon;
  // scale of testlet difficulties
  real<lower=0> theta_upsilon_raw;        // scale of relative item difficulties
  vector<lower=0>[K > 1 ? 1 : 0] theta_tau;
  // scale of step offsets
  vector<lower=0>[M > 1 ? 1 : 0] psi;     // scale of group abilities
  real<lower=0> phi;                      // scale of relative person abilities
}

transformed parameters {
  real<lower=0> theta_upsilon =
    K > 1
    ? sqrt(square(theta_upsilon_raw) + square(theta_tau[1]))
    : theta_upsilon_raw;
  vector[L > 1 ? L : 0] epsilon;
  vector[I] delta;
  vector[K > 1 ? K : 0] tau;
  vector[M > 1 ? M : 0] xi;
  vector[N] eta;
  if (L > 1) epsilon = nu + theta_epsilon[1] * epsilon_raw;
  delta = L > 1 ? epsilon[ll] : rep_vector(nu, I);
  // We want *thresholds* to be distributed the same as binary difficulties.
  for (i in 1:I) {
    delta += (kk[i] == K ? theta_upsilon_raw : theta_upsilon) * upsilon_raw;
  }
  if (K > 1) tau = theta_tau[1] * tau_raw;
  if (M > 1) {
    xi = psi[1] * xi_raw;
    eta = xi[mm] + phi * zeta_raw;
  } else {
    eta = phi * zeta_raw;
  }
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
  tau_raw           ~ std_normal();
  xi_raw            ~ std_normal();
  zeta_raw          ~ std_normal();
  theta_epsilon     ~ std_normal();
  theta_upsilon_raw ~ std_normal();
  theta_tau         ~ std_normal();
  psi               ~ std_normal();
  phi               ~ std_normal();
  for (o in 1:O) {
    int i = ii[o];
    int k = kk[i];
    int n = nn[o];
    real beta = (n < 0 ? xi[-n] : eta[n]) - delta[i];
    // When K == 1, using the rsm function is correct but
    // inefficient. Checking for k == 1 first routes these items
    // to the built-in bernoulli_logit().
    if (k == 1)      target += w[o] * bernoulli_logit_lpmf(y[o] | beta);
    else if (k == K) target += w[o] * pcm_lpmf(y[o] | beta - tau);
    else             target += w[o] * binomial_logit_lpmf(y[o] | k, beta);
  }
}

generated quantities {
  real<lower=0> sigma =
    M > 1 ? sqrt(square(psi[1]) + square(phi)) : phi;
  vector[L > 1 ? L : 0] testlet_difficulty = epsilon / sigma;
  vector[I] item_difficulty = delta / sigma;
  vector[K > 1 ? K : 0] thresholds = tau / sigma;
  vector[M > 1 ? M : 0] group_ability = xi / sigma;
  vector[N] person_ability = eta / sigma;
  int<lower=0> y_rep[O];
  vector[O] log_lik;
  vector[O] log_lik_rep;
  vector[O] log_lik_prior_person;
  vector[O] log_lik_prior_item;
  for (o in 1:O) {
    int i = ii[o];
    int k = kk[i];
    int n = nn[o];
    real prior_epsilon = L > 1 ? normal_rng(nu, theta_epsilon[1]) : nu;
    real prior_delta;
    real prior_xi = M > 1 ? normal_rng(0, psi[1]) : 0;
    real prior_eta = normal_rng(prior_xi, phi);
    real beta = (n < 0 ? xi[-n] : eta[n]) - delta[i];
    real beta_prior_person = (n < 0 ? prior_xi : prior_eta) - delta[i];
    real beta_prior_item = (n < 0 ? xi[-n] : eta[n]) - prior_delta;
    int y_prior_person;
    int y_prior_item;
    if (k == 1) {
      prior_delta = normal_rng(prior_epsilon, theta_upsilon);
      y_rep[o] = bernoulli_logit_rng(beta);
      y_prior_person = bernoulli_logit_rng(beta_prior_person);
      y_prior_item = bernoulli_logit_rng(beta_prior_item);
      log_lik[o] = bernoulli_logit_lpmf(y[o] | beta);
      log_lik_rep[o] = bernoulli_logit_lpmf(y_rep[o] | beta);
      log_lik_prior_person[o] =
        bernoulli_logit_lpmf(y_prior_person | beta_prior_person);
      log_lik_prior_item[o] =
        bernoulli_logit_lpmf(y_prior_item | beta_prior_item);
    } else if (k == K) {
      vector[K] prior_tau = to_vector(normal_rng(0, theta_tau));
      prior_delta = normal_rng(prior_epsilon, theta_upsilon_raw);
      y_rep[o] = pcm_rng(beta - tau);
      y_prior_person = pcm_rng(beta_prior_person - tau);
      y_prior_item = pcm_rng(beta_prior_item - prior_tau);
      log_lik[o] = pcm_lpmf(y[o] | beta - tau);
      log_lik_rep[o] = pcm_lpmf(y_rep[o] | beta - tau);
      log_lik_prior_person[o] =
        pcm_lpmf(y_prior_person | beta_prior_person - tau);
      log_lik_prior_item[o] =
        pcm_lpmf(y_prior_item | beta_prior_item - prior_tau);
    } else {
      prior_delta = normal_rng(prior_epsilon, theta_upsilon);
      y_rep[o] = binomial_rng(k, inv_logit(beta));
      y_prior_person = binomial_rng(k, inv_logit(beta_prior_person));
      y_prior_item = binomial_rng(k, inv_logit(beta_prior_item));
      log_lik[o] = binomial_logit_lpmf(y[o] | k, beta);
      log_lik_rep[o] =
        binomial_logit_lpmf(y_rep[o] | k, beta);
      log_lik_prior_person[o] =
        binomial_logit_lpmf(y_prior_person | k, beta_prior_person);
      log_lik_prior_item[o] =
        binomial_logit_lpmf(y_prior_item | k, beta_prior_item);
    }
  }
}
