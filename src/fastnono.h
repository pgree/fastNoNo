#ifndef __FASTNONO_H__
#define __FASTNONO_H__

#include <RcppEigen.h>
#include <stdlib.h>
#include <cmath>
#include <chrono>
// #include <iostream>

extern "C" {
#include "legeexps.h"
}

struct fit_out {
  Eigen::VectorXd means;
  Eigen::VectorXd sds;
  Eigen::MatrixXd cov;
  double time;
};

fit_out mixed_2group(int nnt, int nn, int n, int k1, int k2, int k,
                     Eigen::MatrixXd a, Eigen::VectorXd y,
                     Eigen::VectorXd ss, double sigy, double sig1);

void lege_nodes_whts(int nn, double t0, double t1,
                     Eigen::VectorXd &ts, Eigen::VectorXd &whts);

void rescale_a(double t, Eigen::MatrixXd &a, Eigen::MatrixXd &asca,
               int n, int k, int k1, int k2);

void eval_inner(int nn, int n, int k1, int k2, int k, double d1, double d2,
                Eigen::MatrixXd b, double t, double resid,
                Eigen::VectorXd ynew, Eigen::VectorXd &dsumsi, double &dsumi,
                double &ss1i, double &ss2i, Eigen::MatrixXd &dsums_covi,
                Eigen::VectorXd &dsum_xsi, Eigen::MatrixXd &xxti, double &fmi);

void get_beta_alpha_prefact(double phi, Eigen::VectorXd ys, Eigen::VectorXd ys2,
                            Eigen::VectorXd s, Eigen::VectorXd s2, int n, int k,
                            double resid, double &alpha, double &beta,
                            double &prefact, double &exp_fact);

void get_phi(int nn, int n, int k, double t, Eigen::VectorXd ysmall,
             Eigen::VectorXd ys2, Eigen::VectorXd s, Eigen::VectorXd s2, double resid,
             double d1, double d2, double &phi0i, double &phi1i,
             double &fmax);

void eval_jac_det(double rho, double phi, double t, double &f);

void eval_logdens_rho(int n, double rho, double a1, double c1, double &f);

void get_int_bds(int n, Eigen::VectorXd phis, Eigen::VectorXd fs, int &i0, int &i1);

void get_mjs(int k, Eigen::VectorXd s2, Eigen::VectorXd ys, double phi, Eigen::VectorXd &vmoms);

void get_xs_to_ws_matrix(int k, int k1, int k2, Eigen::MatrixXd v, double t,
                         Eigen::MatrixXd &a);

void get_xs_from_ws(Eigen::VectorXd ws, int k, int k1, int k2, Eigen::MatrixXd vt,
                    double t, Eigen::VectorXd &xs);


////////////////////////////////////////////////////////////////////

fit_out mixed_2group(int nnt, int nn, int n, int k1, int k2, int k,
                     Eigen::MatrixXd a, Eigen::VectorXd y,
                     Eigen::VectorXd ss, double sigy, double sig1) {
  Eigen::VectorXd dsums = Eigen::VectorXd::Zero(k+2);
  Eigen::VectorXd stds = Eigen::VectorXd::Zero(k+2);
  Eigen::MatrixXd dsums_cov = Eigen::MatrixXd::Zero(k, k);
  double d1, d2, t, fi, wt, fm, ss1, ss2, dsum;
  fit_out fit;

  // time this function run
  std::chrono::high_resolution_clock::time_point tt1, tt2;
  tt1 = std::chrono::high_resolution_clock::now();

  // adjust for fixed priors on coefficients k1+1 to k1+k2
  a.rightCols(k2) = a.rightCols(k2).array().rowwise() * ss.array().square().transpose().array();

  // hyperpriors -- variances, not standard deviations
  d1 = pow(sigy, 2);
  d2 = pow(sig1, 2);

  // before computing integral find least squares solution to
  // ax = y as well as projection of y onto the columns of a
  // an improved (and uglier) a^t * a
  Eigen::MatrixXd ata(Eigen::MatrixXd(k, k).setZero().
		      selfadjointView<Eigen::Lower>().rankUpdate(a.adjoint()));

  // eigendecomposition of a^t * a
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ata);

  Eigen::VectorXd s0(k), s(k);
  s = es.eigenvalues();
  s0 = s.cwiseAbs().cwiseSqrt();
  // eigenvalues are in increasing order, check how many are bigger than 0
  int nnull = 0;
  double tol = 1e-12;
  for (int i=0; i<k; i++) {
    if (std::abs(s[i] / s[k-1]) < tol) {
      nnull = i+1;
    }
  }
  int k0 = k - nnull;
  Eigen::VectorXd s_inv(k0);
  s_inv = s.tail(k0).cwiseInverse();
  Eigen::MatrixXd v0(k, k0);
  v0 = es.eigenvectors().rightCols(k0);

  s_inv = s.tail(k0).cwiseInverse();

  //std::cout << "s0: " << s0 << std::endl;
  Eigen::MatrixXd b(k, k);
  b = s0.asDiagonal() * es.eigenvectors().transpose();

  // v * s_inv * vt * at * y
  Eigen::VectorXd x1(k);
  x1 = v0 * s_inv.asDiagonal() * v0.transpose() * a.transpose() * y;

  Eigen::VectorXd ynew(k);
  // component-wise multiplication of s0 and vt * x1
  ynew = s0.cwiseProduct(es.eigenvectors().transpose() * x1);

  double resid;
  resid = (a * x1 - y).norm();
  resid = pow(resid, 2);

  double t0=0, t1 = M_PI / 2.0;
  Eigen::VectorXd ts(nnt), whts_ts(nnt);
  lege_nodes_whts(nnt, t0, t1, ts, whts_ts);

  // initialize sums (integrals) to be computed
  // dsums - expectations, dsums_cov - covariance
  // ss1 - E[sig1**2], ss2 - E[sig2**2]
  // fm - scaling constant
  dsum = 0;
  ss1 = 0;
  ss2 = 0;
  fm = -1.0e250;
  double dsumi, ss1i, ss2i;
  Eigen::VectorXd dsumsi = Eigen::VectorXd::Zero(k+2);
  Eigen::VectorXd dsum_xsi = Eigen::VectorXd::Zero(k+2);
  Eigen::MatrixXd dsums_covi = Eigen::MatrixXd::Zero(k, k);
  Eigen::MatrixXd xxti = Eigen::MatrixXd::Zero(k, k);

  // theta integral
  for (int i=0; i<nnt; i++) {
    t = ts[nnt - 1 - i];
    wt = whts_ts[nnt - 1 - i];

    // compute phi integral
    eval_inner(nn, n, k1, k2, k, d1, d2, b, t, resid, ynew,
               dsumsi, dsumi, ss1i, ss2i, dsums_covi, dsum_xsi, xxti, fi);

    // due to underflow issues, integrate over theta by computing
    // a sum of the form \sum_i exp(fi)*gi such that at the end
    // we have an expression exp(fm)*dsum
    if (fi > fm) {
      dsum = dsum * exp(fm-fi) + dsumi*wt;

      dsums.head(k) *= exp(fm - fi);
      dsums.head(k) += dsum_xsi.head(k)*wt;
      dsums(k) = dsums(k) * exp(fm-fi) + dsumsi(k)*wt;
      dsums(k+1) = dsums(k+1) * exp(fm-fi) + dsumsi(k+1)*wt;

      ss1 = ss1 * exp(fm-fi) + ss1i*wt;
      ss2 = ss2 * exp(fm-fi) + ss2i*wt;

      dsums_cov *= exp(fm - fi);
      dsums_cov += (dsums_covi + xxti) * wt;

      fm = fi;
     } else {

      dsum = dsum + exp(fi-fm) * dsumi*wt;

      dsums.head(k) += exp(fi - fm)*dsum_xsi.head(k)*wt;
      dsums(k) = dsums(k) + exp(fi-fm)*dsumsi(k)*wt;
      dsums(k+1) = dsums(k+1) + exp(fi-fm)*dsumsi(k+1)*wt;

      ss1 = ss1 + exp(fi-fm)*ss1i*wt;
      ss2 = ss2 + exp(fi-fm)*ss2i*wt;

      dsums_cov += exp(fi - fm) * (dsums_covi + xxti) * wt;
    }
    // if we've gotten to the point that the function is small, break
    if (fi/log(10.0) < fm/log(10.0) - 20.0) break;
  }
  //std::cout << "stds: " << stds << std::endl;

  // scale first and second moments by normalizing constant
  dsums /= dsum;
  dsums_cov /= dsum;

  // dsums_cov now contains E[xx^t], adjust to get covariance
  dsums_cov -= (dsums.head(k) * dsums.head(k).transpose());

  // get stds
  //std::cout << "dsums_cov.diagonal: " << dsums_cov.diagonal() << std::endl;
  stds.head(k) = dsums_cov.diagonal().cwiseSqrt();

  // get variances of sig1, sig2
  stds(k) = sqrt(ss1/dsum - pow(dsums(k), 2));
  stds(k+1) = sqrt(ss2/dsum - pow(dsums(k+1), 2));

  // readjust for scale parameter priors
  for (int j=0; j<k2; j++) {
    dsums(k1 + j) *= pow(ss[j], 2);
    stds(k1 + j) *= pow(ss[j], 2);
  }

  // end timing
  tt2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = tt2 - tt1;

  // fill struct with output
  fit.means = dsums;
  fit.sds = stds;
  fit.cov = dsums_cov;
  fit.time = duration.count()/1e3;

  return fit;
}


void eval_inner(int nn, int n, int k1, int k2, int k, double d1,
                double d2, Eigen::MatrixXd b, double t, double resid,
                Eigen::VectorXd ynew, Eigen::VectorXd &dsumsi, double &dsumi,
                double &ss1i, double &ss2i, Eigen::MatrixXd &dsums_covi,
                Eigen::VectorXd &dsum_xsi, Eigen::MatrixXd &xxti, double &fmi) {
  double phi, alpha, beta, prefact, exp_fact, sig22, sig12, rho,
    a1, f, djac, fi, wt, gi, gi1, gi2, coef;
  Eigen::VectorXd phis(nn), whts_phis(nn);
  Eigen::MatrixXd vt(k, k), v(k, k), c(k, k), ct(k, k), asca(k, k);


  // compute the diagonal form of a^t*a for each theta which
  // will be used to convert the integrand to a diagonal
  // gaussian for each (theta, phi, rho)
  rescale_a(t, b, asca, k, k, k1, k2);
  //std::cout << "asca: " << asca(k-1, k-1) << std::endl;
  Eigen::BDCSVD<Eigen::MatrixXd> svd(asca, Eigen::ComputeThinU | Eigen::ComputeThinV);
  //std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;

  // take svd of asca, that is asca  = u * s * v^t
  Eigen::VectorXd s = svd.singularValues();
  Eigen::MatrixXd u = svd.matrixU();
  v = svd.matrixV();

  // compute the residual of least squares solution
  Eigen::VectorXd ysmall = ynew.transpose() * u;
  // double res2 = ysmall.norm() - ysmall.norm();

  // square entries of two vectors
  Eigen::VectorXd ys2 = ysmall.array().square();
  Eigen::VectorXd s2 = s.array().square();

  // initialize sums taken in inner loop
  // wwti is E[w*w^t] where w is in the coordinate
  // system that depends on theta
  dsumi = 0;
  ss1i = 0;
  ss2i = 0;
  dsumsi *= 0;
  Eigen::MatrixXd wwti = Eigen::MatrixXd::Zero(k, k);
  Eigen::VectorXd dsum_vars = Eigen::VectorXd::Zero(k);
  fmi = -1e250;

  // for each theta, get upper integration bound for phi
  double phi0i, phi1i, fmax_theta;
  get_phi(nn, n, k, t, ysmall, ys2, s, s2, resid, d1, d2, phi0i, phi1i, fmax_theta);
  lege_nodes_whts(nn, phi0i, phi1i, phis, whts_phis);

  // phi integral
  Eigen::VectorXd vmoms(k);
  for (int j=0; j<nn; j++) {
    phi = phis[j];
    wt = whts_phis[j];
    get_beta_alpha_prefact(phi, ysmall, ys2, s, s2, n, k,
                           resid, alpha, beta, prefact, exp_fact);
    get_mjs(k, s2, ysmall, phi, vmoms);

    // compute rho as a function of theta and phi
    sig22 = 1 / pow(tan(t), 2);
    sig12 = -pow(cos(phi), 2) * (sig22 + 1)/(pow(cos(phi), 2) - 1);
    rho = sqrt(sig12 + sig22 + 1);

    a1 = pow(cos(phi), 2) / d1 + pow(sin(phi), 2)*pow(cos(t), 2)/d2;
    eval_logdens_rho(n, rho, a1, exp_fact, f);

    // compute determinant of jacobian
    eval_jac_det(rho, phi, t, djac);

    // for fixed theta, compute integral over phi by a sum of
    // the form \sum_i exp(fi)*gi so that
    // we have an expression exp(fmi)*dsum, the reason we do
    // this is that exp(fi) is often smaller than 10^{-250}
    fi = prefact + f - log(abs(djac));
    gi = 1.0;
    gi1 = rho;
    gi2 = pow(rho, 2);
    if (fi > fmi) {
      dsumi = dsumi * exp(fmi-fi) + gi*wt;

      dsumsi.head(k) *= exp(fmi - fi);
      dsumsi.head(k) += gi * wt * vmoms;
      dsumsi(k) = dsumsi(k)*exp(fmi-fi) + gi1*cos(phi)*wt;
      dsumsi(k+1) = dsumsi(k+1)*exp(fmi-fi) + gi1*sin(phi)*cos(t)*wt;

      ss1i = ss1i*exp(fmi-fi) + gi2*pow(cos(phi), 2)*wt;
      ss2i = ss2i*exp(fmi-fi) + gi2*pow(sin(phi), 2)*pow(cos(t), 2)*wt;

      for (int ijk=0; ijk<k; ijk++) {
        coef = 1/(s2(ijk)/pow(cos(phi), 2) + 1/pow(sin(phi), 2));
        dsum_vars(ijk) = dsum_vars(ijk) * exp(fmi-fi) + coef*gi2*wt;
      }

      wwti.array() *= exp(fmi - fi);
      wwti = wwti + gi*(vmoms * vmoms.transpose())*wt;

      fmi = fi;

    } else {

      dsumi = dsumi + exp(fi-fmi) * gi*wt;

      dsumsi.head(k) += exp(fi - fmi) * gi * wt * vmoms;
      dsumsi(k) = dsumsi(k) + exp(fi-fmi)*gi1*cos(phi)*wt;
      dsumsi(k+1) = dsumsi(k+1) + exp(fi-fmi)*gi1*sin(phi)*cos(t)*wt;

      ss1i = ss1i + exp(fi-fmi)*gi2*pow(cos(phi), 2)*wt;
      ss2i = ss2i + exp(fi-fmi)*gi2*pow(sin(phi), 2)*pow(cos(t), 2)*wt;

      for (int ijk=0; ijk<k; ijk++) {
        coef = 1/(s2(ijk)/pow(cos(phi), 2) + 1/pow(sin(phi), 2));
        dsum_vars(ijk) = dsum_vars(ijk) + exp(fi-fmi)*coef*gi2*wt;
      }

      wwti = wwti + exp(fi - fmi)*gi*(vmoms * vmoms.transpose())*wt;
    }

  }
  // for covariance we need to convert back to original coordinate
  // system. we have a different change of variables for each value
  // of theta, so do this for each theta
  get_xs_to_ws_matrix(k, k1, k2, v, t, c);
  dsums_covi = c * dsum_vars.asDiagonal() * c.transpose();

  // recover matrix E[x*x^t]
  xxti = c * wwti * c.transpose();
  get_xs_from_ws(dsumsi, k, k1, k2, v.transpose(), t, dsum_xsi);
}


void get_xs_from_ws(Eigen::VectorXd ws,
                    int k, int k1, int k2,
                    Eigen::MatrixXd vt, double t,
                    Eigen::VectorXd &xs) {

  // convert expectations back from the diagonal coordinate (ws)
  // system to the original coordinate system (xs)

  for (int i=0; i<k; i++) {
    xs(i) = 0;
    for (int j=0; j<k; j++) {
      xs(i) = xs(i) + vt(j, i) * ws(j);
    }
  }

  for (int i=0; i<k1; i++) {
    xs(i) = xs(i) * cos(t);
  }

  for (int i=0; i<k2; i++) {
    xs(i + k1) = xs(i + k1) * sin(t);
  }
}


void get_xs_to_ws_matrix(int k, int k1, int k2,
                         Eigen::MatrixXd v, double t, Eigen::MatrixXd &a) {
  Eigen::VectorXd s1(k);

  // construct matrix that goes from diagonal coordinate
  // system (depending on t) to original coordinate system
  for (int i=0; i<k1; i++) {
    s1(i) = cos(t);
  }

  for (int i=0; i<k2; i++) {
    s1(i+k1) = sin(t);
  }

  a = s1.asDiagonal() * v;

}

void get_mjs(int k, Eigen::VectorXd s2, Eigen::VectorXd ys, double phi,
             Eigen::VectorXd &vmoms) {

  // compute conditional expectations of gaussian
  double sp = sin(phi);
  double sp2 = sp*sp;
  double cp = pow(cos(phi), 2);

  for (int i=0; i<k; i++) {
    vmoms(i) = ys(i) * sqrt(s2(i)) * sp2 / (s2(i)*sp2 + cp);
  }

}


void get_phi(int nn, int n, int k, double t, Eigen::VectorXd ysmall,
             Eigen::VectorXd ys2, Eigen::VectorXd s, Eigen::VectorXd s2, double resid,
             double d1, double d2, double &phi0i, double &phi1i, double &fmax) {
  double phi, alpha, beta, prefact, exp_fact, a1, sig22, sig12, rho, djac, f;
  Eigen::VectorXd fs(nn), phis(nn), whts_phis(nn);

  /*
     find the upper integration bound of the integral with
     respect to phi. do this by evaluating the integral on
     a sparse grid and checking when the value decreases below
     10^{-18} of its maximum. also return the maximum value of
     the integrand
  */

  // lay down nodes
  double phi0 = 0.0;
  double phi1 = M_PI / 2.0;
  lege_nodes_whts(nn, phi0, phi1, phis, whts_phis);

  // tabulate function
  for (int j=0; j<nn; j++) {
    phi = phis[j];
    get_beta_alpha_prefact(phi, ysmall, ys2, s, s2, n, k,
                           resid, alpha, beta, prefact, exp_fact);
    a1 = pow(cos(phi), 2) / d1;
    a1 = a1 + pow(sin(phi), 2)*pow(cos(t), 2) / d2;
    //std::cout << "a1: " << a1 << std::endl;

    // compute rho as a function of theta and phi
    sig22 = 1/pow(tan(t), 2);
    sig12 = -pow(cos(phi), 2) * (sig22+1)/(pow(cos(phi), 2) - 1);
    rho = sqrt(sig12 + sig22 + 1);

    eval_jac_det(rho, phi, t, djac);
    eval_logdens_rho(n, rho, a1, exp_fact, f);

    fs(j) = prefact + f - log(djac);
  }

  // get bounds of integration
  int i0, i1;
  get_int_bds(nn, phis, fs, i0, i1);
  phi0i = phis[i0];
  phi1i = phis[i1];
  if (i0 == 0) {
    phi0i = phi0;
  }
  if (i1 == nn) {
    phi1i = phi1;
  }

  // check that bounds are reasonable
  // fmax = fs.maxCoeff();
  // double r_dd = fmax - fs[i0];
  // double l_dd = fmax - fs[i1];
  //std::cout << "r_dd: " << r_dd << std::endl;
  //std::cout << "l_dd: " << l_dd << std::endl;
}


void get_int_bds(int nn, Eigen::VectorXd phis, Eigen::VectorXd fs, int &i0, int &i1) {
  // this subroutine takes as input an array of the form
  // (1, 1, ..., 1, 0, ..., 0, 1, 1, ..., 1) and returns
  // the indices of the last 1 before the 0s and the first
  // 1 after the 0s

  // tolerance
  double tol = 16;

  // get maximum of function
  Eigen::Index ind_max;
  double fmax = fs.maxCoeff(&ind_max);
  double fmax_phi = phis[ind_max];

  // find points at which function is within 1e-16 of its maximum
  Eigen::VectorXi tmpvec1 = (fs.array() > (fmax - tol)).cast<int>();

  // get maximum bound
  //Eigen::Map<Eigen::VectorXd> phis_eig(&phis[0], nn);
  Eigen::VectorXi tmpvec = ((tmpvec1.array() == 1) && (phis.array() >= fmax_phi)) \
    .cast<int>();
  for (int i=0; i<nn; i++) {
    if (tmpvec[i] > 0) {
      i1 = i;
    }
  }
  i1 = std::min(i1 + 1, nn);

  // get minimum bound
  tmpvec = ((tmpvec1.array() == 1) && (phis.array() <= fmax_phi)) \
    .cast<int>();
  for (int i=0; i<nn; i++) {
    if (tmpvec[i] > 0) {
      i0 = i;
      break;
    }
  }
  i0 = std::max(i0 - 1, 0);
}



void eval_logdens_rho(int n, double rho, double a1, double c1, double &f) {
  f = -n*log(rho) - pow(rho, 2)/2.0*a1 - c1/(2*pow(rho, 2));
}


void eval_jac_det(double rho, double phi, double t, double &f) {
  // compute the determinant of the change of variables
  // jacobian for the mixed effects model

  double sig1, sig2, x1, y1, alph, djac;

  sig1 = rho*cos(phi);
  sig2 = rho*sin(phi)*cos(t);
  x1 = sig1;
  y1 = sig2;

  alph = pow(x1, 2) + pow(y1, 2) + 1;
  djac = 1 / (1 + pow(y1, 2));
  djac *= (1/sqrt(alph)-pow(x1, 2)/pow(alph, 1.5))/sqrt(1- pow(x1, 2)/alph);
  f = abs(djac);
}



void get_beta_alpha_prefact(double phi, Eigen::VectorXd ys, Eigen::VectorXd ys2,
                            Eigen::VectorXd s, Eigen::VectorXd s2, int n, int k,
                            double resid, double &alpha, double &beta,
                            double &prefact, double &exp_fact) {

  // for convenience, compute quantities that appear in
  // several locations in the integrand. prefact is the
  // part of the integrand that doesn't depend on rho

  double sp = pow(sin(phi), 2);
  double cp = pow(cos(phi), 2);

  // alpha
  alpha = 0;
  for (int i=0; i<k; i++) {
    alpha = alpha - log(cp+sp*s2(i));
  }
  alpha = alpha/2.0 + k/2.0 * log(2 * M_PI);

  // beta
  beta = 0;
  for (int i=0; i<k; i++) {
    beta = beta + ys2(i)/(s2(i)*sp+cp);
  }

  // prefact
  prefact = -(n-k)/2.0*log(cp)+alpha;
  exp_fact= resid/cp+beta;
}



void rescale_a(double t, Eigen::MatrixXd &a, Eigen::MatrixXd &asca,
               int n, int k, int k1, int k2) {

  // multiply a by a diagonal matrix
  double tc = cos(t);
  double ts = sin(t);

  for (int i=0; i<n; i++) {
    for (int j=0; j<k1; j++) {
      asca(i, j) = a(i, j) * tc;
    }
    for (int j=0; j<k2; j++) {
      asca(i, k1 + j) = a(i, k1 + j) * ts;
    }
  }
}


void lege_nodes_whts(int nn, double t0, double t1,
                     Eigen::VectorXd &ts, Eigen::VectorXd &whts) {
  // this function wraps legeexps.c, which was constructed from
  // Vladimir Rokhlin's legeexps.f via a call to f2c

  // convert to appropriate data types of calling sequence
  long int nn0[1] = {nn};
  long int itype[1] = {1};
  double ts0[nn], whts0[nn], dummy[1];
  legeexps_(itype, nn0, ts0, dummy, dummy, whts0);

  // copy nodes and weights into Eigen::VectorXd and scale to fit interval
  for (int i=0; i<nn; i++) {
    ts(i) = t0 + (t1 - t0) * (ts0[i] + 1)/2.0;
    whts(i) = whts0[i] * (t1 - t0) / 2.0;
  }
}


# endif
