#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
extern "C" {
#include "legeexps.c"
}
#include <chrono>
using namespace std::chrono;
//#include <Rcpp.h>
//using namespace Rcpp;

using namespace std;
using namespace Eigen;


void lege_nodes_whts(int nn, double t0, double t1, 
		     VectorXd &ts, VectorXd &whts);

void mixed_read_params(string filename, int &n, int &k1, int &k2, MatrixXd &a, 
		       VectorXd &y, VectorXd &ss, double &sigy, 
		       double &sig1);

void rescale_a(double t, MatrixXd &a, MatrixXd &asca, int n, int k, 
	       int k1, int k2);

void eval_inner(int nn, int n, int k1, int k2, int k, double d1, 
		double d2, MatrixXd b, double t, double resid,
		VectorXd ynew, VectorXd &dsumsi, double &dsumi,
		double &ss1i, double &ss2i, MatrixXd &dsums_covi,
		VectorXd &dsum_xsi, MatrixXd &xxti, double &fmi);

void get_beta_alpha_prefact(double phi, VectorXd ys, VectorXd ys2,
			    VectorXd s, VectorXd s2, int n, int k,
			    double resid, double &alpha, double &beta,
			    double &prefact, double &exp_fact);

void get_phi(int nn, int n, int k, double t, VectorXd ysmall, 
	     VectorXd ys2, VectorXd s, VectorXd s2, double resid,
	     double d1, double d2, double &phi0i, double &phi1i,
	     double &fmax);

void eval_jac_det(double rho, double phi, double t, double &f);

void eval_logdens_rho(int n, double rho, double a1, double c1, double &f);

void get_int_bds(int n, VectorXd phis, VectorXd fs, int &i0, int &i1);

void get_mjs(int k, VectorXd s2, VectorXd ys, double phi, VectorXd &vmoms);

void mixed_2group(int nnt, int nn, int n, int k1, int k2, int k, 
		  MatrixXd a, VectorXd y, VectorXd ss, 
		  double sigy, double sig1, 
		  VectorXd &dsums, double &dsum, VectorXd &stds, 
		  MatrixXd &cov);

void get_xs_to_ws_matrix(int k, int k1, int k2, MatrixXd v, double t,
			 MatrixXd &a);

void get_xs_from_ws(VectorXd ws, int k, int k1, int k2, MatrixXd vt, 
		    double t, VectorXd &xs);


////////////////////////////////////////////////////////////////////


void test_mixed_effects() {
  string filename = "params.dat";
  //cout << "filename: " << filename << endl;
  int n, k1, k2, k3, k, nn, nnt, nn2, nnt2;
  MatrixXd a;
  double sigy, sig1, dsum;
  VectorXd y, ss;

  // dimensions of problem
  mixed_read_params(filename, n, k1, k2, a, y, ss, sigy, sig1);
  k = k1 + k2;
  //cout << "n: " << n << endl;
  //cout << "k1: " << k1 << endl;
  //cout << "k: " << k << endl;

  // allocate memory accordingly
  VectorXd dsums(k+2), dsums2(k+2), stds(k+2), stds2(k+2), dds(k+2);
  MatrixXd cov(k, k);

  // run algorithm
  high_resolution_clock::time_point t1, t2;
  t1 = high_resolution_clock::now();
  nn = 80;
  nnt = 40;
  mixed_2group(nnt, nn, n, k1, k2, k, a, y, ss, sigy, sig1,
	       dsums, dsum, stds, cov);
  t2 = high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = t2 - t1;
  cout << "total time: " << duration.count()/1e6 << endl;
  //cout << "dsums: " << dsums << endl;
  //cout << "stds: " << stds << endl;

  // double number of nodes
  nn2 = 2 * nn;
  nnt2 = 2 * nnt;
  mixed_2group(nnt2, nn2, n, k1, k2, k, a, y, ss, sigy, sig1,
	       dsums2, dsum, stds2, cov);

  // difference
  dds = dsums - dsums2;
  cout << "max posterior mean merror: " << dds.cwiseAbs().maxCoeff() << endl;
  dds = stds - stds2;
  cout << "max posterior mean merror: " << dds.cwiseAbs().maxCoeff() << endl;

}


void mixed_2group(int nnt, int nn, int n, int k1, int k2, int k, 
		  MatrixXd a, VectorXd y, VectorXd ss, 
		  double sigy, double sig1, 
		  VectorXd &dsums, double &dsum, VectorXd &stds, 
		  MatrixXd &dsums_cov) {
  double d1, d2, t, wht_t, fi, wt, fm, ss1, ss2, tmp;

  // adjust for fixed priors on coefficients k1+1 to k1+k2
  MatrixXd a2 = a;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k2; j++) {
      a2(i, k1 + j) = a2(i, k1 + j) * pow(ss[j], 2);
    }
  }

  // hyperpriors -- variances, not standard deviations
  d1 = pow(sigy, 2);
  d2 = pow(sig1, 2);

  // before computing integral find least squares solution to
  // ax = y as well as projection of y onto the columns of a
  MatrixXd ata(k, k);
  ata = a2.transpose() * a2;

  // eigendecomposition of a^t * a
  SelfAdjointEigenSolver<MatrixXd> es(ata);
  VectorXd s_inv(k), s0(k), s(k);
  s = es.eigenvalues();
  for (int i=0; i<k; i++) {
    s_inv[i] = 1 / es.eigenvalues()[i];
    s0[i] = sqrt(abs(s[i]));
  }
  //cout << "s0: " << s0 << endl;

  MatrixXd b(k, k);
  b = s0.asDiagonal() * es.eigenvectors().transpose();

  // v * s_inv * vt * a2t * y
  VectorXd x1(k);
  x1 = es.eigenvectors() * s_inv.asDiagonal() * es.eigenvectors().transpose()\
    * a2.transpose() * y;

  VectorXd ynew(k), tmp_vec(k);
  tmp_vec = es.eigenvectors().transpose() * x1;
  for (int i=0; i<k; i++) {
    ynew[i] = s0[i] * tmp_vec[i];
  }

  double resid;
  resid = (a2 * x1 - y).norm();
  resid = pow(resid, 2);
  //cout << "resid: " << resid << endl;

  double t0=0, t1 = M_PI / 2.0;
  VectorXd ts(nnt), whts_ts(nnt);
  lege_nodes_whts(nnt, t0, t1, ts, whts_ts);

  // initialize sums (integrals) to be computed
  // dsums - expectations
  // dsums_cov - covariance
  // ss1 - E[sig1**2]
  // ss2 - E[sig2**2]
  // fm - scaling constant
  dsums *= 0;
  dsum = 0;
  dsums_cov *= 0;
  ss1 = 0;
  ss2 = 0;
  fm = -1.0e250;
  double dsumi, ss1i, ss2i, fmi;
  VectorXd dsumsi(k+2), dsum_xsi(k+2); 
  MatrixXd dsums_covi(k, k), xxti(k, k);

  // theta integral
  for (int i=0; i<nnt; i++) {
    t = ts[nnt - 1 - i];
    wt = whts_ts[nnt - 1 - i];

    // compute phi integral
    eval_inner(nn, n, k1, k2, k, d1, d2, b, t, resid, ynew, dsumsi, 
	       dsumi, ss1i, ss2i, dsums_covi, dsum_xsi, xxti, fi);

    // due to underflow issues, integrate over theta by computing 
    // a sum of the form \sum_i exp(fi)*gi such that at the end
    // we have an expression exp(fm)*dsum 
    if (fi > fm) {
      dsum = dsum * exp(fm-fi) + dsumi*wt;

      for (int ijk=0; ijk<k; ijk++) {
	dsums(ijk) = dsums(ijk) * exp(fm-fi) + dsum_xsi(ijk)*wt;
      }

      dsums(k) = dsums(k) * exp(fm-fi) + dsumsi(k)*wt;
      dsums(k+1) = dsums(k+1) * exp(fm-fi) + dsumsi(k+1)*wt;

      ss1 = ss1 * exp(fm-fi) + ss1i*wt;
      ss2 = ss2 * exp(fm-fi) + ss2i*wt;

      for (int i1=0; i1<k; i1++) {
	for (int i2=0; i2<k; i2++) {
          tmp = dsums_covi(i1,i2) + xxti(i1,i2);
          dsums_cov(i1,i2) = dsums_cov(i1,i2)*exp(fm-fi) + tmp*wt;
	}
      }

      fm = fi;
     } else {

      dsum = dsum + exp(fi-fm) * dsumi*wt;

      for (int ijk=0; ijk<k; ijk++) {
	dsums(ijk) = dsums(ijk) + exp(fi-fm)*dsum_xsi(ijk)*wt;
      }

      dsums(k) = dsums(k) + exp(fi-fm)*dsumsi(k)*wt;
      dsums(k+1) = dsums(k+1) + exp(fi-fm)*dsumsi(k+1)*wt;

      ss1 = ss1 + exp(fi-fm)*ss1i*wt;
      ss2 = ss2 + exp(fi-fm)*ss2i*wt;

      for (int i1=0; i1<k; i1++) {
	for (int i2=0; i2<k; i2++) {
          tmp = dsums_covi(i1,i2)+xxti(i1,i2);
          dsums_cov(i1,i2)=dsums_cov(i1,i2) + exp(fi-fm)*tmp*wt;
	}
      }
    }
    // if we've gotten to the point that the function is small, break
    if (fi/log(10.0) < fm/log(10.0) - 20.0) break;
  }


  // scale first and second moments by normalizing constant
  dsums /= dsum;
  dsums_cov /= dsum;

  // dsums_cov now contains E[xx^t], adjust to get covariance 
  dsums_cov -= (dsums.head(k) * dsums.head(k).transpose());

  // get stds
  stds.head(k) = dsums_cov.diagonal().cwiseSqrt();

  // get variances of sig1, sig2
  stds(k) = sqrt(ss1/dsum - pow(dsums(k), 2));
  stds(k+1) = sqrt(ss2/dsum - pow(dsums(k+1), 2));

  // readjust for scale parameter priors
  for (int j=0; j<k2; j++) {
    dsums(k1 + j) *= pow(ss[j], 2);
    stds(k1 + j) *= pow(ss[j], 2);
  }


}


void eval_inner(int nn, int n, int k1, int k2, int k, double d1, 
		double d2, MatrixXd b, double t, double resid,
		VectorXd ynew, VectorXd &dsumsi, double &dsumi,
		double &ss1i, double &ss2i, MatrixXd &dsums_covi,
		VectorXd &dsum_xsi, MatrixXd &xxti, double &fmi) {
  double phi, alpha, beta, prefact, exp_fact, sig22, sig12, rho,
    a1, f, djac, fi, wt, gi, gi1, gi2, coef;
  VectorXd phis(nn), whts_phis(nn);
  MatrixXd vt(k, k), v(k, k), c(k, k), ct(k, k), tmp_mat(k, k), 
    asca(k, k);


  // compute the diagonal form of a^t*a for each theta which
  // will be used to convert the integrand to a diagonal 
  // gaussian for each (theta, phi, rho) 
  rescale_a(t, b, asca, k, k, k1, k2);
  //cout << "asca: " << asca(k-1, k-1) << endl;

  JacobiSVD<MatrixXd> svd(asca, ComputeThinU | ComputeThinV);
  //cout << "Its singular values are:" << endl << svd.singularValues() << endl;

  // take svd of asca, that is asca  = u * s * v^t
  VectorXd s = svd.singularValues();
  MatrixXd u = svd.matrixU();
  v = svd.matrixV();

  // compute the residual of least squares solution 
  VectorXd ysmall = ynew.transpose() * u;
  double res2 = ysmall.norm() - ysmall.norm();

  // square entries of two vectors
  VectorXd ys2 = ysmall.array().square();   
  VectorXd s2 = s.array().square();   

  // initialize sums taken in inner loop
  // wwti is E[w*w^t] where w is in the coordinate 
  // system that depends on theta
  dsumi = 0;
  ss1i = 0;
  ss2i = 0;
  dsumsi *= 0;
  MatrixXd wwti = MatrixXd::Zero(k, k);
  VectorXd dsum_vars = VectorXd::Zero(k);
  fmi = -1e250;

  // for each theta, get upper integration bound for phi
  double phi0i, phi1i, fmax_theta;
  get_phi(nn, n, k, t, ysmall, ys2, s, s2, resid, d1, d2, phi0i, 
	  phi1i, fmax_theta);
  lege_nodes_whts(nn, phi0i, phi1i, phis, whts_phis);
  //cout << "phi0i: " << phi0i << endl;
  //cout << "phi1i: " << phi1i << endl;


  // phi integral
  VectorXd vmoms(k);
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

      for (int ijk=0; ijk<k; ijk++) {
	dsumsi(ijk) = dsumsi(ijk) * exp(fmi - fi) + gi*wt*vmoms(ijk);
      }

      dsumsi(k) = dsumsi(k)*exp(fmi-fi) + gi1*cos(phi)*wt;
      dsumsi(k+1) = dsumsi(k+1)*exp(fmi-fi) + gi1*sin(phi)*cos(t)*wt;

      ss1i = ss1i*exp(fmi-fi) + gi2*pow(cos(phi), 2)*wt;
      ss2i = ss2i*exp(fmi-fi) + gi2*pow(sin(phi), 2)*pow(cos(t), 2)*wt;

      for (int ijk=0; ijk<k; ijk++) {
	coef = 1/(s2(ijk)/pow(cos(phi), 2) + 1/pow(sin(phi), 2));
	dsum_vars(ijk) = dsum_vars(ijk) * exp(fmi-fi) + coef*gi2*wt;
      }

      for (int i1=0; i1<k; i1++) {
	for (int i2=0; i2<k; i2++) {
          wwti(i1, i2)=wwti(i1,i2)*exp(fmi-fi)+gi*vmoms(i1)*vmoms(i2)*wt;
	}
      }

      fmi = fi;

    } else { 
      
      dsumi = dsumi + exp(fi-fmi) * gi*wt;

      for (int ijk=0; ijk<k; ijk++) {
	dsumsi(ijk) = dsumsi(ijk) + exp(fi - fmi) * gi*wt*vmoms(ijk);
      }

      dsumsi(k) = dsumsi(k) + exp(fi-fmi)*gi1*cos(phi)*wt;
      dsumsi(k+1) = dsumsi(k+1) + exp(fi-fmi)*gi1*sin(phi)*cos(t)*wt;

      ss1i = ss1i + exp(fi-fmi)*gi2*pow(cos(phi), 2)*wt;
      ss2i = ss2i + exp(fi-fmi)*gi2*pow(sin(phi), 2)*pow(cos(t), 2)*wt;

      for (int ijk=0; ijk<k; ijk++) {
	coef = 1/(s2(ijk)/pow(cos(phi), 2) + 1/pow(sin(phi), 2));
        dsum_vars(ijk) = dsum_vars(ijk) + exp(fi-fmi)*coef*gi2*wt;
      }

      for (int i1=0; i1<k; i1++) {
	for (int i2=0; i2<k; i2++) {
          wwti(i1, i2)=wwti(i1,i2)+exp(fi-fmi)*gi*vmoms(i1)*vmoms(i2)*wt;
	}
      }
    }

  }
  //cout << "dsum_vars: " << dsum_vars << endl;

  // for covariance we need to convert back to original coordinate
  // system. we have a different change of variables for each value
  // of theta, so do this for each theta
  get_xs_to_ws_matrix(k, k1, k2, v, t, c);
  dsums_covi = c * dsum_vars.asDiagonal() * c.transpose();

  // recover matrix E[x*x^t]
  xxti = c * wwti * c.transpose();
  get_xs_from_ws(dsumsi, k, k1, k2, v.transpose(), t, dsum_xsi);

}


void get_xs_from_ws(VectorXd ws, int k, int k1, int k2, MatrixXd vt, 
		    double t, VectorXd &xs) {

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


void get_xs_to_ws_matrix(int k, int k1, int k2, MatrixXd v, double t,
			 MatrixXd &a) {
  VectorXd s1(k);

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

void get_mjs(int k, VectorXd s2, VectorXd ys, double phi, VectorXd &vmoms) {

  // compute conditional expectations of gaussian
  double sp = sin(phi);
  double sp2 = sp*sp;
  double cp = pow(cos(phi), 2);

  for (int i=0; i<k; i++) {
    vmoms(i) = ys(i) * sqrt(s2(i)) * sp2 / (s2(i)*sp2 + cp);
  }

}


void get_phi(int nn, int n, int k, double t, VectorXd ysmall, 
	     VectorXd ys2, VectorXd s, VectorXd s2, double resid,
	     double d1, double d2, double &phi0i, double &phi1i,
	     double &fmax) {
  double phi, alpha, beta, prefact, exp_fact, a1, sig22, sig12, rho,
    djac, f, fi;
  VectorXd fs(nn), tmp_mat(nn);
  VectorXd phis(nn), whts_phis(nn);

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
    //cout << "a1: " << a1 << endl;

    // compute rho as a function of theta and phi
    sig22 = 1/pow(tan(t), 2);
    sig12 = -pow(cos(phi), 2) * (sig22+1)/(pow(cos(phi), 2) - 1);
    rho = sqrt(sig12 + sig22 + 1);

    eval_jac_det(rho, phi, t, djac);
    //cout << "djac: " << djac << endl;
    eval_logdens_rho(n, rho, a1, exp_fact, f);
    //cout << "f: " << f << endl;

    fs(j) = prefact + f - log(djac);
    //cout << "fs(j): " << fs(j) << endl;
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
  fmax = fs.maxCoeff();  
  double r_dd = fmax - fs[i0];
  double l_dd = fmax - fs[i1];
  //cout << "r_dd: " << r_dd << endl;
  //cout << "l_dd: " << l_dd << endl;

}


void get_int_bds(int nn, VectorXd phis, VectorXd fs, int &i0, int &i1) {
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
  VectorXi tmpvec1 = (fs.array() > (fmax - tol)).cast<int>();

  // get maximum bound
  Map<VectorXd> phis_eig(&phis[0], nn); 
  VectorXi tmpvec = ((tmpvec1.array() == 1) && (phis_eig.array() > fmax_phi))\
    .cast<int>();
  for (int i=0; i<nn; i++) {
    if (tmpvec[i] > 0) {
      i1 = i;
    }
  }
  i1 = min(i1 + 1, nn);

  // get minimum bound
  tmpvec = ((tmpvec1.array() == 1) && (phis_eig.array() < fmax_phi))\
    .cast<int>();
  for (int i=0; i<nn; i++) {
    if (tmpvec[i] > 0) {
      i0 = i;
      break;
    }
  }
  i0 = max(i0 - 1, 0);

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



void get_beta_alpha_prefact(double phi, VectorXd ys, VectorXd ys2,
			    VectorXd s, VectorXd s2, int n, int k,
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



void rescale_a(double t, MatrixXd &a, MatrixXd &asca, int n, int k, 
	       int k1, int k2) {

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



void lege_nodes_whts(int nn, double t0, double t1, VectorXd &ts, 
		     VectorXd &whts) {
  // this function wraps legeexps.c, which was constructed from
  // Vladimir Rokhlin's legeexps.f via a call to f2c

  // convert to appropriate data types of calling sequence 
  integer nn0[1] = {nn};
  integer itype[1] = {1};
  double ts0[nn], whts0[nn], dummy[1], tmp;
  legeexps_(itype, nn0, ts0, dummy, dummy, whts0);
  
  // copy nodes and weights into VectorXd and scale to fit 
  // interval 
  for (int i=0; i<nn; i++) {
    ts(i) = t0 + (t1 - t0) * (ts0[i] + 1)/2.0;
    whts(i) = whts0[i] * (t1 - t0) / 2.0;
  }

}


void mixed_read_params(string filename, int &n, int &k1, int &k2, MatrixXd &a, VectorXd &y, VectorXd &ss, double &sigy, double &sig1) {
  int k;

  ifstream input_file(filename);

  // dimensions of matrix
  input_file >> n >> k1 >> k2;
  //cout << "n: " << n << endl;
  //cout << "k1: " << k1 << endl;
  //cout << "k2: " << k2 << endl;

  // read data matrix, a
  k = k1 + k2;
  a.resize(n, k);
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < n; j++) {
      input_file >> a(j, i);
    }
  }
  //cout << "a: " << a(0,0) << endl;

  // read y
  y.resize(n);
  for (int i = 0; i < n; i++) {
    input_file >> y(i);
  }
  //cout << "y: " << y[99] << endl;

  // read ss
  ss.resize(k2);
  for (int i = 0; i < k2; i++) {
    input_file >> ss(i);
  }

  // read scale hyperpriors
  input_file >> sigy >> sig1;

  input_file.close();

}


int main() {

  test_mixed_effects();

  return 0;
}



