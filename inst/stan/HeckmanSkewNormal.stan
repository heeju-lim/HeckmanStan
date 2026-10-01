functions {
  // Inverse Mills ratio m = phi(x) / Phi(x) and V = m * (x + m), 0 < V < 1.
  // For x < -30 an asymptotic expansion is used, because the direct formula
  // loses precision there (x + m suffers cancellation).
  vector mills(real x) {
    vector[2] out;
    if (x < -30) {
      real u = inv_square(x);
      out[1] = -x * (1 + u - 2 * square(u) + 10 * u * square(u));
      out[2] = 1 - u + 6 * square(u);
    } else {
      out[1] = exp(std_normal_lpdf(x) - std_normal_lcdf(x));
      out[2] = out[1] * (x + out[1]);
    }
    return out;
  }

  // For l(z) = log phi(z) + log Phi((a - r z) / s) returns [ l'(z), -l''(z) ].
  // l is concave and -l'' lies in [1, 1 / s^2].
  vector lgrad(real z, real a, real r, real s) {
    real k = r / s;
    vector[2] m = mills((a - r * z) / s);
    return [ -z - k * m[1], 1 + square(k) * m[2] ]';
  }

  // Distance over which a log-concave integrand drops by D (log units), given
  // the outward slope g >= 0 and a lower bound c on the curvature.
  // Written in the cancellation-free form of (-g + sqrt(g^2 + 2 c D)) / c.
  real tail_width(real g, real c, real D) {
    return 2 * D / (g + sqrt(square(g) + 2 * c * D));
  }

  // log Phi2(a, b; rho), robust version used for |rho| > 0.75 or tiny probabilities.
  //
  // Evaluated as log int_{-inf}^{b} phi(z) Phi((a - rho z) / s) dz,
  // s = sqrt(1 - rho^2), entirely on the log scale.  The integrand is
  // log-concave, so Gauss-Legendre panels are placed at its mode and at the
  // "shoulder" z = a / rho (where Phi(.) switches from ~1 to ~0), and the
  // outer panels end where the integrand has dropped by exp(-40).
  real log_Phi2_quad(real a, real b, real rho, data matrix gl) {
    int n = rows(gl);
    real D = 40;  // panel coverage: integrand within exp(-D) of its maximum
    real K = 6;   // shoulder half-width in units of s / |rho|
    real r = fmin(fmax(rho, -1 + 1e-12), 1 - 1e-12);
    real s = sqrt((1 - r) * (1 + r));

    // 1) mode of the integrand on (-inf, b]  (Newton on a concave function)
    real zm = b;
    vector[2] g = lgrad(b, a, r, s);
    if (g[1] < 0) {
      for (it in 1:100) {
        real step = g[1] / g[2];
        zm += step;
        g = lgrad(zm, a, r, s);
        if (abs(step) < 1e-10) break;
      }
      zm = fmin(zm, b);
      g = lgrad(zm, a, r, s);
    }

    // 2) breakpoints: mode, shoulder, and edge of the shoulder on the "on" side
    real sup_lo = zm - tail_width(fmax(g[1], 0), 1, D);
    real zs = a / r;
    real zon = zs - K * s / r;
    int use_zs = (zs > sup_lo) && (zs < b);
    int use_zon = (zon > sup_lo) && (zon < b);
    real pL = zm;
    real pR = zm;
    if (use_zs) { pL = fmin(pL, zs); pR = fmax(pR, zs); }
    if (use_zon) { pL = fmin(pL, zon); pR = fmax(pR, zon); }

    // 3) outer panels: curvature grows to the right when r > 0, to the left when r < 0
    vector[2] gL = lgrad(pL, a, r, s);
    vector[2] gR = lgrad(pR, a, r, s);
    real cL = r < 0 ? gL[2] : 1;
    real cR = r > 0 ? gR[2] : 1;
    real lo = pL - tail_width(fmax(gL[1], 0), cL, D);
    real hi = fmin(b, pR + tail_width(fmax(-gR[1], 0), cR, D));

    vector[5] e = sort_asc([ lo, zm, use_zs ? zs : zm, use_zon ? zon : zm, hi ]');

    // 4) Gauss-Legendre on each panel, summed on the log scale
    vector[4] lp = rep_vector(negative_infinity(), 4);
    vector[n] term;
    for (j in 1:4) {
      real w = e[j + 1] - e[j];
      if (w > 0) {
        for (i in 1:n) {
          real z = e[j] + 0.5 * w * (gl[i, 1] + 1);
          term[i] = gl[i, 2] - 0.5 * square(z) + std_normal_lcdf((a - r * z) / s);
        }
        lp[j] = log(0.5 * w) + log_sum_exp(term);
      }
    }
    return log_sum_exp(lp) - 0.5 * log(2 * pi());
  }

  // log Phi2(a, b; rho) via the Drezner-Wesolowsky representation (as in the
  // previous version), but summed on the log scale:
  //   Phi2 = Phi(a) Phi(b) + 1/(2 pi) int_0^{asin rho} exp(-(a^2 + b^2 - 2ab sin t) / (2 cos^2 t)) dt
  // Fast and accurate for |rho| <= 0.75 unless Phi2 is tiny.
  real log_Phi2_drezner(real a, real b, real rho, data matrix gl) {
    int n = rows(gl);
    real sg = rho > 0 ? 1 : -1;
    real al = asin(abs(rho));
    vector[n] term;
    for (i in 1:n) {
      real t = 0.5 * al * (gl[i, 1] + 1);
      term[i] = gl[i, 2]
                - (square(a) + square(b) - 2 * sg * a * b * sin(t)) / (2 * square(cos(t)));
    }
    real lI = log(0.5 * al) - log(2 * pi()) + log_sum_exp(term);
    real lb = std_normal_lcdf(a) + std_normal_lcdf(b);
    if (rho > 0) return log_sum_exp(lb, lI);
    if (lI < lb) return log_diff_exp(lb, lI);
    return negative_infinity();
  }

  // log Phi2(a, b; rho) = log P(Z1 <= a, Z2 <= b), corr(Z1, Z2) = rho.
  // Fast path (Drezner, 10 nodes) when it is accurate to < 1e-8 relative;
  // otherwise the robust log-scale quadrature (16 nodes).
  real log_Phi2(real a0, real b0, real rho, data matrix gl10, data matrix gl16) {
    if (is_nan(a0) || is_nan(b0) || is_nan(rho)) return not_a_number();
    // Phi(1e6) == 1 and log Phi(-1e6) = -5e11 in double precision, so clamping
    // changes nothing but keeps +-inf / absurd values from producing NaN.
    real a = fmin(fmax(a0, -1e6), 1e6);
    real b = fmin(fmax(b0, -1e6), 1e6);
    real v;
    if (abs(rho) < 1e-12)
      return std_normal_lcdf(a) + std_normal_lcdf(b);
    if (abs(rho) <= 0.75) {
      v = log_Phi2_drezner(a, b, rho, gl10);
      if (v > log(1e-6)) return fmin(v, 0);
    }
    v = log_Phi2_quad(a, b, rho, gl16);
    return fmin(v, 0);
  }

  // Pointwise log-likelihood (length N)
  vector heckman_sn_loglik(vector y, vector Xb, vector Zg, array[] int D,
                           real omega, real auxmt, real auxkt,
                           real Omega11, real rhoaaux,
                           data matrix gl10, data matrix gl16) {
    int N = num_elements(D);
    vector[N] ll;
    real sOm = sqrt(Omega11);
    int ny = 1;
    for (n in 1:N) {
      if (D[n] == 1) {
        real res = y[ny] - Xb[ny];
        real mut1 = Zg[n] + auxmt * res;
        real mut2 = auxkt * res;
        // log(aa - pp) = log Phi2(mut1 / sqrt(Omega11), mut2; -rhoaaux)
        ll[n] = normal_lpdf(res | 0, omega) + log(2)
                + log_Phi2(mut1 / sOm, mut2, -rhoaaux, gl10, gl16);
        ny += 1;
      } else {
        ll[n] = std_normal_lcdf(-Zg[n]);
      }
    }
    return ll;
  }
}

data {
  // dimensions
  int<lower=1> N;
  int<lower=1, upper=N> N_y;
  int<lower=1> p;
  int<lower=1> q;
  // covariates
  matrix[N_y, p] X;   // rows = selected observations (D == 1), in order
  matrix[N, q] Z;
  // responses
  array[N] int<lower=0, upper=1> D;
  vector[N_y] y;
}

transformed data {
  // Gauss-Legendre nodes (column 1) and log-weights (column 2) on [-1, 1]
  matrix[10, 2] gl10;
  matrix[16, 2] gl16;
  gl10[ : , 1] = [ -0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                   -0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
                   0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
                   0.9739065285171717 ]';
  gl10[ : , 2] = log([ 0.0666713443086881, 0.1494513491505806, 0.2190863625159820,
                       0.2692667193099963, 0.2955242247147529, 0.2955242247147529,
                       0.2692667193099963, 0.2190863625159820, 0.1494513491505806,
                       0.0666713443086881 ]');
  gl16[ : , 1] = [ -0.9894009349916499, -0.9445750230732326, -0.8656312023878318,
                   -0.755404408355003, -0.6178762444026438, -0.45801677765722737,
                   -0.2816035507792589, -0.09501250983763744, 0.09501250983763744,
                   0.2816035507792589, 0.45801677765722737, 0.6178762444026438,
                   0.755404408355003, 0.8656312023878318, 0.9445750230732326,
                   0.9894009349916499 ]';
  gl16[ : , 2] = log([ 0.027152459411754176, 0.062253523938647456, 0.0951585116824926,
                       0.12462897125553407, 0.1495959888165767, 0.16915651939500265,
                       0.18260341504492364, 0.18945061045506864, 0.18945061045506864,
                       0.18260341504492364, 0.16915651939500265, 0.1495959888165767,
                       0.12462897125553407, 0.0951585116824926, 0.062253523938647456,
                       0.027152459411754176 ]');
  if (sum(D) != N_y)
    reject("sum(D) = ", sum(D), " but N_y = ", N_y,
           ". y and X must contain exactly the selected observations, in order.");
}

parameters {
  vector[p] beta;
  vector[q] gamma;
  real<lower=-1, upper=1> rho;
  real<lower=0> sigma;
  real lambda;
}

transformed parameters {
  real<lower=0> sigma2 = square(sigma);
  real<lower=0> omega = sqrt(sigma2 + square(lambda));  // scale of the skew-normal error
  real lambdat = -lambda * rho / omega;
  real auxmt = sigma * rho / square(omega);
  real auxkt = lambda / (sigma * omega);                 // = -lambdat / (sigma * rho), no 0/0
  real<lower=0> Omega11 = (1 - square(rho)) + square(lambdat);
  real rhoaaux = -lambdat / sqrt(Omega11);
}

model {
  // priors (normal(0, 10) = the previous multi_normal(0, 100 * I))
  beta ~ normal(0, 10);
  gamma ~ normal(0, 10);
  sigma ~ normal(0, 5);    // half-normal because of <lower=0>; adjust to the scale of y
  lambda ~ normal(0, 5);   // adjust to the scale of y
  // rho: uniform(-1, 1)

  target += sum(heckman_sn_loglik(y, X * beta, Z * gamma, D,
                                  omega, auxmt, auxkt, Omega11, rhoaaux,
                                  gl10, gl16));
}

generated quantities {
  vector[N] log_lik = heckman_sn_loglik(y, X * beta, Z * gamma, D,
                                        omega, auxmt, auxkt, Omega11, rhoaaux,
                                        gl10, gl16);
  real log_lik_total = sum(log_lik);
  // E[eps] is not zero: the intercept in beta is shifted by err_mean
  real err_mean = lambda * sqrt(2 / pi());
  real err_sd = sqrt(sigma2 + square(lambda) * (1 - 2 / pi()));
}
