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
        int<lower=1> L;               // number of testlets
        int<lower=1> I;               // number of items
        // J is reserved for facets in a future implementation
        int<lower=1> M;               // number of groups
        int<lower=1> N;               // number of students
        int<lower=1> O;               // number of observations
        int<lower=1,upper=M> mm[N];   // group for student n
        int<lower=-M,upper=N> nn[O];  // student for observation m
                                      //   negative value => group observation
        int<lower=1,upper=L> ll[I];   // testlet for item i
        int<lower=1,upper=I> ii[O];   // item for observation o
        int<lower=1> kk[I];           // maximum score for item i
                                      //   1   = Rasch
                                      //   2-4 = rating scale
                                      //   5+  = binomial
                                      //   TODO: Add more flixible typing.
        int<lower=0> y[O];            // score for observation o
}

transformed data {
        int K_MAX = 4;                  // maximum allowable rating-scale score
        int<lower=0,upper=K_MAX> K = 0; // number of rating-scale thresholds
        for (i in 1:I) {
                if (kk[i] > 1 && kk[i] <= K_MAX) {
                        if (K == 0) {
                                K = kk[i];
                        } else if (K == kk[i]) {
                                continue;
                        } else {
                                reject(
                                        "All rating-scale items must have ",
                                        "the same number of categories."
                                )
                        }
                }
        }
}

parameters {
        // The parameters are mostly in LISREL notation, with some abuse
        // of ν, δ, ε, and θ for items and testlets.
        real nu;                        // mean difficulty
        vector[L] epsilon_raw;          // relative difficulty for testlet l
        vector[I] upsilon_raw;          // relative difficulty for item i
        vector[K > 1] d_tau_raw;        // spacing of last two thresholds.
                                        //   Vanishes if no rating scales.
        vector[max(0, K-2)] d_tau_err;  // relative spacing for threshold t.
                                        //   Vanishes if no rating scales
                                        //   or a single threshold.
        vector[M] xi_raw;               // ability for group j
        vector[N] zeta_raw;             // relative ability for student n
        real<lower=0> theta_epsilon;    // scale of testlet difficulties
        real<lower=0> theta_upsilon;    // scale of relative item difficulties
        real<lower=0> psi;              // scale of group abilities
        real<lower=0> phi;              // scale of relative student abilities
}

transformed parameters {
        vector[L] epsilon;
        vector[I] delta;
        vector[K] tau;
        vector[M] xi;
        vector[N] eta;
        epsilon = nu + theta_epsilon * epsilon_raw;
        delta = epsilon[ll] + theta_upsilon * upsilon_raw;
        if (K > 1) {
                // Our prior is that thresholds should be ascending,
                // disordered thresholds are possible. Use the
                // traditional log(2) separation requirement
                // (Linacre, 2002) as a prior mean.
                vector[K-1] d_tau =
                        log(2) + d_tau_raw[1] + append_row(d_tau_err, 0);
                vector[K] tau_raw = cumulative_sum(append_row(0, d_tau));
                // Centre between the last two thresholds.
                tau =  tau_raw - tau_raw[K] + 0.5 * d_tau[K-1];
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
        nu             ~ std_normal();
        epsilon_raw    ~ std_normal();
        upsilon_raw    ~ std_normal();
        d_tau_raw      ~ std_normal();
        d_tau_err      ~ std_normal();
        xi_raw         ~ std_normal();
        zeta_raw       ~ std_normal();
        theta_epsilon  ~ std_normal();
        theta_upsilon  ~ std_normal();
        psi            ~ std_normal();
        phi            ~ std_normal();
        for (o in 1:O) {
                int i = ii[o];
                int k = kk[i];
                int n = nn[o];
                real beta = (n < 0 ? xi[-n] : eta[n]) - delta[i];
                if (k < 2)           y[o] ~ bernoulli_logit(beta);
                else if (k <= K_MAX) y[o] ~ rsm(beta - tau);
                else                 y[o] ~ binomial_logit(k, beta);
        }
}

generated quantities {
        real<lower=0> gamma = inv_sqrt(square(psi) + square(phi));
        vector[L] testlet_difficulty = gamma * epsilon;
        real prior_testlet_difficulty = gamma * normal_rng(nu, theta_epsilon);
        vector[I] item_difficulty = gamma * delta;
        real prior_item_difficulty =
                prior_testlet_difficulty + gamma * normal_rng(0, theta_upsilon);
        vector[K] thresholds = gamma * tau;
        vector[M] group_ability = gamma * xi;
        real prior_group_ability = gamma * normal_rng(0, psi);
        vector[N] person_ability = gamma * eta;
        real prior_person_ability =
                prior_group_ability + gamma * normal_rng(0, phi);
        // Expected ratings and entropies can be computed from
        // replications, but we need two prior versions, one with
        // persons freed and one with items freed.
        int<lower=0> y_rep[O];
        vector[O] log_lik;
        vector[O] log_lik_rep;
        for (o in 1:O) {
                int i = ii[o];
                int k = kk[i];
                int n = nn[o];
                real beta = (n < 0 ? xi[-n] : eta[n]) - delta[i];
                if (k < 2) {
                        y_rep[o] = bernoulli_logit_rng(beta);
                        log_lik[o] = bernoulli_logit_lpmf(y[o] | beta);
                        log_lik_rep[o] = bernoulli_logit_lpmf(y_rep[o] | beta);
                }
                else if (k < K_MAX) {
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
