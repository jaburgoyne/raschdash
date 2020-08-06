/*
* Rasch rating-scale model for four-level criteria, with missing data possible.
*/

functions {
        // Rasch rating-scale model
        real rsm_lpmf(int y, vector beta_raw) {
                int K = num_elements(beta_raw);
                vector[K+1] beta = append_row(0, cumulative_sum(beta_raw));
                return categorical_logit_lpmf(y + 1 | beta);
        }
        int rsm_rng(vector beta_raw) {
                int K = num_elements(beta_raw);
                vector[K+1] beta = append_row(0, cumulative_sum(beta_raw));
                return categorical_logit_rng(beta) - 1;
        }
}

data {
        int<lower=1> J;               // number of groups
        int<lower=1> N;               // number of students
        int<lower=1> L;               // number of testlets
        int<lower=1> I;               // number of items
        int<lower=1,upper=3> K_raw;   // number of rating-scale thresholds
        int<lower=1> M;               // number of observations
        int<lower=1,upper=J> jj[N];   // group for student n
        int<lower=-J,upper=N> nn[M];  // student for observation m
                                      //   negative value => group observation
        int<lower=1,upper=L> ll[I];   // testlet for item i
        int<lower=1,upper=I> ii[M];   // item for observation m
        int<lower=1> kk[I];           // max score of item i
                                      //   TODO: Add error checking
                                      //   1   = Rasch
                                      //   2-3 = rating scale
                                      //   4+  = binomial
        int<lower=0> y[M];        // score for observation m
        real alpha;                   // centre of ability scale
        real gamma;                   // standard deviation of ability scale
}

transformed data {
        // Eliminate tau if there are no rating scales
        int K = K_raw < 2 ? 0 : K_raw;
}

// The parameters are mostly in LISREL notation, with some abuse of
// ν, δ, ε, and θ for items and testlets.
parameters {
        vector[J] xi_raw;               // ability for group j
        vector[N] zeta_raw;             // relative ability for student n
        real nu;                        // mean difficulty
        vector[L] epsilon_raw;          // relative difficulty for testlet l
        vector[I] upsilon_raw;          // relative difficulty for item i
        vector[K] tau_raw;              // (rating thresholds) - delta
        real<lower=0> psi;              // scale of group abilities
        real<lower=0> phi;              // scale of relative student abilities
        real<lower=0> theta_epsilon;    // scale of testlet difficulties
        real<lower=0> theta_upsilon;    // scale of relative item difficulties
}

transformed parameters {
        vector[J] xi;
        vector[N] eta;
        vector[L] epsilon;
        vector[I] delta;
        vector[K] tau;
        xi = psi * xi_raw;
        eta = xi[jj] + phi * zeta_raw;
        epsilon = nu + theta_epsilon * epsilon_raw;
        delta = epsilon[ll] + theta_upsilon * upsilon_raw;
        tau = K < 2 ? tau_raw : tau_raw - 0.5 * (tau_raw[K-1] + tau_raw[K]);
        // Check for illegal item types.
        // TODO: Make more flexible item types, including Poisson
        //       and partial credit.
        for (i in 1:I) {
                if (kk[i] > 1 && kk[i] < 4 && kk[i] != K) {
                        reject(
                                "All rating scales must have ",
                                K,
                                " thresholds."
                                );
                }
        }
}

model {
        // Gelman's standard normal prior everywhere is sensible here
        // and runs fast; testing for failure with normal(0, 100) has
        // essentially the same result. The Nuffic distribution
        // is closer to t_10, which allows for some outliers, but it runs
        // slower and fits equivalently to the standard normal.
        xi_raw        ~ std_normal();
        zeta_raw      ~ std_normal();
        nu            ~ std_normal();
        epsilon_raw   ~ std_normal();
        upsilon_raw   ~ std_normal();
        tau_raw       ~ std_normal();
        psi           ~ std_normal();
        phi           ~ std_normal();
        theta_epsilon ~ std_normal();
        theta_upsilon ~ std_normal();
        {
                for (m in 1:M) {
                        int n = nn[m];
                        int i = ii[m];
                        int j = n < 0 ? -n : jj[n];
                        int k = kk[i];
                        real beta = (n < 0 ? xi[j] : eta[n]) - delta[i];
                        if (k < 2)      y[m] ~ bernoulli_logit(beta);
                        else if (k < 4) y[m] ~ rsm(beta - tau);
                        else            y[m] ~ binomial_logit(k, beta);
                }
        }
}

generated quantities {
        real<lower=0> gamma_logit = gamma / sqrt(dot_self([psi, phi]));
        vector[J] group_ability = alpha + gamma_logit * xi;
        real prior_group_ability = gamma_logit * normal_rng(0, psi);
        vector[N] person_ability = alpha + gamma_logit * eta;
        real prior_person_ability =
                prior_group_ability + gamma_logit * normal_rng(0, phi);
        vector[L] testlet_difficulty = alpha + gamma_logit * epsilon;
        real prior_testlet_difficulty =
                alpha + gamma_logit * normal_rng(nu, theta_epsilon);
        vector[I] item_difficulty = alpha + gamma_logit * delta;
        real prior_item_difficulty =
                prior_testlet_difficulty
                + gamma_logit * normal_rng(0, theta_upsilon);
        vector[K] thresholds = gamma_logit * tau;
        // Expected ratings and entropies can be computed from
        // replications, but we need two prior versions, one with
        // persons freed and one with items freed.
        int<lower=0> y_rep[M];
        vector[M] log_lik;
        vector[M] log_lik_rep;
        for (m in 1:M) {
                int n = nn[m];
                int i = ii[m];
                int j = n < 0 ? -n : jj[n];
                int k = kk[i];
                real beta = (n < 0 ? xi[j] : eta[n]) - delta[i];
                if (k < 2) {
                        y_rep[m] = bernoulli_logit_rng(beta);
                        log_lik[m] = bernoulli_logit_lpmf(y[m] | beta);
                        log_lik_rep[m] = bernoulli_logit_lpmf(y_rep[m] | beta);
                }
                else if (k < 4) {
                        y_rep[m] = rsm_rng(beta - tau);
                        log_lik[m] = rsm_lpmf(y[m] | beta - tau);
                        log_lik_rep[m] = rsm_lpmf(y_rep[m] | beta - tau);
                } else {
                        y_rep[m] = binomial_rng(k, inv_logit(beta));
                        log_lik[m] = binomial_logit_lpmf(y[m] | k, beta);
                        log_lik_rep[m] =
                                binomial_logit_lpmf(y_rep[m] | k, beta);
                }
        }
}
