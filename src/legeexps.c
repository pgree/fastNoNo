/* legeexps.f -- translated by f2c (version 20190311).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublecomplex c_b125 = {1.,0.};
static integer c__2 = 2;





/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*        this is the end of the debugging code and the beginning */
/*        of the legendre expansion routines */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */



/*        This file contains a set of subroutines for the handling */
/*        of Legendre expansions. It contains 19 subroutines that are */
/*        user-callable. Following is a brief description of these */
/*        subroutines. */

/*   legeexps - constructs Legendre nodes, and  corresponding Gaussian */
/*        weights. Also constructs the matrix v converting the */
/*         coefficients of a legendre expansion into its values at */
/*         the n Gaussian nodes, and its inverse u, converting the */
/*         values of a function at n Gaussian nodes into the */
/*         coefficients of the corresponding Legendre series. */

/*   legepol - evaluates a single Legendre polynomial (together */
/*         with its derivative) at the user-provided point */

/*   legepols - evaluates a bunch of Legendre polynomials */
/*         at the user-provided point */
/*   legepls2 - an accelerated version of legepols, evaluating a */
/*         Legendre polynomials at the user-provided point; maximum */
/*         order of the polynomials to be evaluated is 290 */

/*   legepolders - evaluates a bunch of Legendre polynomials */
/*         at the user-provided point, and the derivatives of the */
/*         said polynomials */
/*   legeinmt - for the user-specified n, constructs the matrices of */
/*        spectral indefinite integration differentiation on the n */
/*        Gaussian nodes on the interval [-1,1]. */

/*   legeinte - computes the indefinite integral of the legendre */
/*        expansion polin getting the expansion polout */

/*   legediff -  differentiates the legendre expansion polin getting */
/*        the expansion polout */

/*   legefder - computes the value and the derivative of a Legendre */
/*        expansion at point X in interval [-1,1]; this subroutine */
/*        is not designed to be very efficient, but it does not */
/*        use any exdternally supplied arrays */

/*   legefde2 - the same as legefder, except it is desigmed to be */
/*        fairly efficient; it uses externally supplied arrays */
/*        that are precomputed */

/*   legeexev - computes the value of a Legendre expansion with */
/*        at point X in interval [-1,1]; same as legefder, but does */
/*        not compute the derivative of the expansion */

/*   legeexe2 - the same as legeexev, except it is desigmed to be */
/*        fairly efficient; it uses externally supplied arrays */
/*        that are precomputed */

/*   lematrin - constructs the matrix interpolating functions from */
/*        the n-point Gaussian grid on the interval [-1,1] to an */
/*        arbitrary m-point grid (the nodes of the latter are */
/*        user-provided) */

/*   levecin - constructs the coefficients of the standard */
/*        interpolation formula connecting the values of a */
/*        function at n Gaussian nodes on the interval [a,b] with */
/*        its value at the point x \in R^1 */

/*   legeodev - evaluates at the point x a Legendre expansion */
/*        having only odd-numbered elements; this is a fairly */
/*        efficient code, using external arrays that are */
/*        precomputed */

/*   legeevev - evaluates at the point x a Legendre expansion */
/*        having only even-numbered elements; this is a fairly */
/*        efficient code, using external arrays that are */
/*        precomputed */

/*   legepeven - evaluates even-numbered Legendre polynomials */
/*        of the argument x; this is a fairly efficient code, */
/*        using external arrays that are precomputed */

/*   legepodd - evaluates odd-numbered Legendre polynomials */
/*        of the argument x; this is a fairly efficient code, */
/*        using external arrays that are precomputed */

/*   legefdeq - computes the value and the derivative of a */
/*        Legendre Q-expansion with coefficients coefs */
/*     at point X in interval (-1,1); please note that this is */
/*     the evil twin of the subroutine legefder, evaluating the */
/*     proper (P-function) Legendre expansion; this subroutine */
/*        is not designed to be very efficient, but it does not */
/*        use any exdternally supplied arrays */

/*   legeq - calculates the values and derivatives of a bunch */
/*        of Legendre Q-functions at the user-specified point */
/*        x on the interval (-1,1) */

/*   legeqs - calculates the value and the derivative of a single */
/*        Legendre Q-function at the user-specified point */
/*        x on the interval (-1,1) */

/*   legecfde - computes the value and the derivative of a Legendre */
/*        expansion with complex coefficients at point X in interval */
/*        [-1,1]; this subroutine is not designed to be very efficient, */
/*        but it does not use any exdternally supplied arrays. This is */
/*        a complex version of the subroutine legefder. */

/*   legecfd2 - the same as legecfde, except it is designed to be */
/*        fairly efficient; it uses externally supplied arrays */
/*        that are precomputed. This is a complex version of the */
/*        subroutine legefde2. */

/*   legecva2 - the same as legecfd2, except it is does not evaluate */
/*        the derivative of the function */
/*   legerts - an improved code for the construction of Gaussian */
/*        quadratures. Its asymptotic CPU time requirements are of */
/*        the order $O(n)$; it has been tested for n \leq 100 000. */


/* Subroutine */ int legeexps_(integer *itype, integer *n, doublereal *x, 
	doublereal *u, doublereal *v, doublereal *whts)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int legerts2_(integer *, integer *, doublereal *, 
	    doublereal *), legepols_(doublereal *, integer *, doublereal *);
    static doublereal d__;
    static integer i__, j, itype_rts__;


/*         this subroutine constructs the gaussiaqn nodes */
/*         on the interval [-1,1], and the weights for the */
/*         corresponding order n quadrature. it also constructs */
/*         the matrix v converting the coefficients */
/*         of a legendre expansion into its values at the n */
/*         gaussian nodes, and its inverse u, converting the */
/*         values of a function at n gaussian nodes into the */
/*         coefficients of the corresponding legendre series. */
/*         no attempt has been made to make this code efficient, */
/*         but its speed is normally sufficient, and it is */
/*         mercifully short. */

/*                 input parameters: */

/*  itype - the type of the calculation to be performed */
/*          itype=0 means that only the gaussian nodes are */
/*                  to be constructed. */
/*          itype=1 means that only the nodes and the weights */
/*                  are to be constructed */
/*          itype=2 means that the nodes, the weights, and */
/*                  the matrices u, v are to be constructed */
/*  n - the number of gaussian nodes and weights to be generated */

/*                 output parameters: */

/*  x - the order n gaussian nodes - computed independently */
/*          of the value of itype. */
/*  u - the n*n matrix converting the  values at of a polynomial of order */
/*         n-1 at n legendre nodes into the coefficients of its */
/*         legendre expansion - computed only in itype=2 */
/*  v - the n*n matrix converting the coefficients */
/*         of an n-term legendre expansion into its values at */
/*         n legendre nodes (note that v is the inverse of u) */
/*          - computed only in itype=2 */
/*  whts - the corresponding quadrature weights - computed only */
/*         if itype .ge. 1 */

/*       . . . construct the nodes and the weights of the n-point gaussian */
/*             quadrature */

/* ccc        ifwhts=0 */
/* ccc        if(itype. gt. 0) ifwhts=1 */
/* ccc        call legewhts(n,x,whts,ifwhts) */

    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    u_dim1 = *n;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --x;
    --whts;

    /* Function Body */
    itype_rts__ = 0;
    if (*itype > 0) {
	itype_rts__ = 1;
    }

    legerts2_(&itype_rts__, n, &x[1], &whts[1]);

/*       construct the matrix of values of the legendre polynomials */
/*       at these nodes */

    if (*itype != 2) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n - 1;
	legepols_(&x[i__], &i__2, &u[i__ * u_dim1 + 1]);
/* L1400: */
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    v[i__ + j * v_dim1] = u[j + i__ * u_dim1];
/* L1600: */
	}
/* L1800: */
    }

/*       now, v converts coefficients of a legendre expansion */
/*       into its values at the gaussian nodes. construct its */
/*       inverse u, converting the values of a function at */
/*       gaussian nodes into the coefficients of a legendre */
/*       expansion of that function */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__ = 1.;
	d__ = d__ * ((i__ << 1) - 1) / 2;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    u[i__ + j * u_dim1] = v[j + i__ * v_dim1] * whts[j] * d__;
/* L2600: */
	}
/* L2800: */
    }
    return 0;
} /* legeexps_ */






/* Subroutine */ int legerts2_(integer *itype, integer *n, doublereal *ts, 
	doublereal *whts)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal d__;
    extern /* Subroutine */ int legetayl2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *);
    static doublereal h__;
    static integer i__, k;
    static doublereal d2, h2;
    static integer n2;
    static doublereal x0, x1;
    static integer ii, kk;
    static doublereal pi, der, eps, pol, der3, pol3, done, xold, d_den__;
    static integer ifodd;
    static doublereal theta, derold, polold;
    static integer ifstop;
    extern /* Subroutine */ int legepol_(doublereal *, integer *, doublereal *
	    , doublereal *);


/*        This subroutine constructs the Gaussian quadrature */
/*        or order n. Its claim to fame is the fact that the */
/*        cost of the calculation is proportional to n; in */
/*        practice, with n=10 000 the calculation is more or */
/*        less instantaneous. PLEASE NOTE THAT THIS IS A */
/*        MILDLY OPTIMIZED - and much less readable - VERSION */
/*        OF THE SUBROUTINE LEGERTS (SEE) */


/*                 Input parameters: */

/*  itype - the type of calculation desired: */
/*     itype=1 will cause both the roots and the weights to be returned */
/*     itype=0 will cause only the roots to be returned */
/*  n - the number of nodes to be returned */

/*                 Output parameters: */

/*  ts - the n Gaussian nodes on the interval [-1,1] */
/*  whts - the n Gaussian weights on the interval [-1,1] */


/*        . . . determine the number of Taylor coefficients */
/*              to be used */

    /* Parameter adjustments */
    --whts;
    --ts;

    /* Function Body */
    k = 30;
    eps = 1e-8;

    d__ = 1.;
    d2 = d__ + 1e-24;
    if (d2 != d__) {
	k = 54;
	eps = 1e-13;
    }

/*       . . . construct the array of initial approximations */
/*             to the roots of the n-th legendre polynomial */

    i__ = *n / 2;
    ifodd = *n - (i__ << 1);

    done = 1.;
    pi = atan(done) * 4;
    h__ = pi / (*n << 1);
    ii = 0;

    d_den__ = pi / ((*n << 2) + 2);
    i__1 = *n;
    for (i__ = *n / 2 + 1; i__ <= i__1; ++i__) {

	++ii;
	theta = ((i__ << 2) - done) * d_den__;
	ts[ii] = -cos(theta);
/* L1100: */
    }

/*       starting from the center, find roots one after another */

    pol = 1.;
    der = 0.;

    x0 = 0.;
    legepol_(&x0, n, &pol, &der);
    x1 = ts[1];

    n2 = (*n + 1) / 2;

    pol3 = pol;
    der3 = der;

    i__1 = n2;
    for (kk = 1; kk <= i__1; ++kk) {

	if (kk == 1 && ifodd == 1) {
	    ts[kk] = x0;
	    whts[kk] = der;
	    x0 = x1;
	    x1 = ts[kk + 1];
	    pol3 = pol;
	    der3 = der;
	    goto L2000;
	}

/*        conduct newton */

	ifstop = 0;
	for (i__ = 1; i__ <= 10; ++i__) {

	    if (i__ != 1) {
		h2 = x1 - xold;

		if (abs(h2) < 1e-36) {
		    goto L1600;
		}

		if (i__ == 2) {
		    i__2 = k / 2;
		    legetayl2_(&polold, &derold, &xold, &h2, n, &i__2, &pol, &
			    der);
		}
		if (i__ != 2) {
		    i__2 = k / 5;
		    legetayl2_(&polold, &derold, &xold, &h2, n, &i__2, &pol, &
			    der);
		}

		polold = pol;
		derold = der;
		xold = x1;
	    }

	    if (i__ == 1) {
		h__ = x1 - x0;
		legetayl2_(&pol3, &der3, &x0, &h__, n, &k, &pol, &der);

		polold = pol;
		derold = der;
		xold = x1;
	    }

	    x1 -= pol / der;

	    if (abs(pol) < eps) {
		++ifstop;
	    }
	    if (ifstop == 3) {
		goto L1600;
	    }

/* L1400: */
	}
L1600:

	ts[kk] = x1;
	if (*itype > 0) {
	    whts[kk] = der;
	}

	x0 = x1;
	x1 = ts[kk + 1];
	pol3 = pol;
	der3 = der;
L2000:
	;
    }

/*        put the obtained roots in the proper order */

    for (i__ = (*n + 1) / 2; i__ >= 1; --i__) {

	ts[i__ + *n / 2] = ts[i__];
/* L2200: */
    }

    i__1 = *n / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	ts[i__] = -ts[*n - i__ + 1];
/* L2400: */
    }

    if (*itype <= 0) {
	return 0;
    }

/*        put the obtained roots in the proper order */

    for (i__ = (*n + 1) / 2; i__ >= 1; --i__) {

	whts[i__ + *n / 2] = whts[i__];
/* L2600: */
    }

    i__1 = *n / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	whts[i__] = whts[*n - i__ + 1];
/* L2800: */
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

/* Computing 2nd power */
	d__1 = ts[i__];
/* Computing 2nd power */
	d__2 = whts[i__];
	whts[i__] = 2 / (1 - d__1 * d__1) / (d__2 * d__2);
/* L3600: */
    }

    return 0;
} /* legerts2_ */






/* Subroutine */ int legetayl2_(doublereal *pol, doublereal *der, doublereal *
	x, doublereal *h__, integer *n, integer *k, doublereal *sum, 
	doublereal *sumder)
{
    /* Initialized data */

    static integer ifcalled = 0;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer i__;
    static doublereal d7, q0, q1, q2, dn, qi, ddd, two, qip1, dddd, half, 
	    done, hinv, prods[100], rnums[100], prodinv[100], squares[100];


/*       initialize */

    if (ifcalled == 1) {
	goto L2000;
    }

    done = 1.;
    two = 2.;
    half = 1 / two;
    for (i__ = 1; i__ <= 90; ++i__) {

	prods[i__ - 1] = i__ * done * (i__ + 1);
/* Computing 2nd power */
	i__1 = i__;
	squares[i__ - 1] = (doublereal) (i__1 * i__1);
	prodinv[i__ - 1] = 1 / prods[i__ - 1];
	rnums[i__ - 1] = (doublereal) i__;
/* L1200: */
    }

    ifcalled = 1;
L2000:

/*        . . . evaluate the derivatives of P_n scaled by h^n/n!, */
/*              and sum the taylor series for P_n and its */
/*              derivative */

    dn = *n * (*n + done);
/* Computing 2nd power */
    d__1 = *x;
    d7 = done / (done - d__1 * d__1);
    q0 = *pol;
    q1 = *der * *h__;
    q2 = (two * *x * *der - dn * *pol) * d7;
/* Computing 2nd power */
    d__1 = *h__;
    q2 = q2 * (d__1 * d__1) * half;

    hinv = done / *h__;
    *sum = q0 + q1 + q2;
    *sumder = (q1 + q2 + q2) * hinv;

/* ccc        if(k .le. 2) return */

    qi = q1;
    qip1 = q2;

/* Computing 2nd power */
    d__1 = *h__;
    ddd = d__1 * d__1 * d7;
    dddd = two * *x * hinv;

    i__1 = *k - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	d__ = dddd * squares[i__] * qip1 - (dn - prods[i__ - 1]) * qi;
	d__ = d__ * prodinv[i__] * ddd;

	*sum += d__;
	*sumder += d__ * rnums[i__ + 1] * hinv;

	qi = qip1;
	qip1 = d__;
/* L2200: */
    }

    return 0;
} /* legetayl2_ */






/* Subroutine */ int legewhts_old__(integer *n, doublereal *ts, doublereal *
	whts, integer *ifwhts)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal a, b, h__;
    static integer i__, k;
    static doublereal t, fm, fp, pi, xk, der, eps, pol, done, zero, delta;
    static integer ifout;
    static doublereal deltold;
    extern /* Subroutine */ int legepol_(doublereal *, integer *, doublereal *
	    , doublereal *), prodend_(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *);


/*        this subroutine constructs the nodes and the */
/*        weights of the n-point gaussian quadrature on */
/*        the interval [-1,1] */

/*                input parameters: */

/*  n - the number of nodes in the quadrature */

/*                output parameters: */

/*  ts - the nodes of the n-point gaussian quadrature */
/*  w - the weights of the n-point gaussian quadrature */

/*       . . . construct the array of initial approximations */
/*             to the roots of the n-th legendre polynomial */

    /* Parameter adjustments */
    --whts;
    --ts;

    /* Function Body */
    eps = 1e-14;
    zero = 0.;
    done = 1.;
    pi = atan(done) * 4;
    h__ = pi / (*n << 1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = ((i__ << 1) - 1) * h__;
	ts[*n - i__ + 1] = cos(t);
/* L1200: */
    }

/*         use newton to find all roots of the legendre polynomial */

    ts[*n / 2 + 1] = 0.;
    i__1 = *n / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	xk = ts[i__];
	ifout = 0;
	deltold = 1.;
	for (k = 1; k <= 10; ++k) {
	    legepol_(&xk, n, &pol, &der);
	    delta = -pol / der;
/* cccc         call prin2('delta=*',delta,1) */
	    xk += delta;
	    if (abs(delta) < eps) {
		++ifout;
	    }

/* ccc        call prin2('delta=*',delta,1) */
	    if (ifout == 3) {
		goto L1600;
	    }
/* L1400: */
	}
L1600:
	ts[i__] = xk;
	ts[*n - i__ + 1] = -xk;
/* L2000: */
    }

/*       now, use the explicit integral formulae */
/*       to obtain the weights */

    if (*ifwhts == 0) {
	return 0;
    }
    a = -1.;
    b = 1.;
    i__1 = *n / 2 + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	prodend_(&a, &ts[1], n, &i__, &fm);
	prodend_(&b, &ts[1], n, &i__, &fp);
	whts[i__] = fp - fm;
	whts[*n - i__ + 1] = whts[i__];
/* L2200: */
    }
    return 0;
} /* legewhts_old__ */






/* Subroutine */ int legewhts_(integer *n, doublereal *ts, doublereal *whts, 
	integer *ifwhts)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal h__;
    static integer i__, k;
    static doublereal t, pi, xk, der, eps, pol, sum;
    extern /* Subroutine */ int legepol_sum__(doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal done, zero, delta;
    static integer ifout;
    static doublereal deltold;


/*        this subroutine constructs the nodes and the */
/*        weights of the n-point gaussian quadrature on */
/*        the interval [-1,1] */

/*                input parameters: */

/*  n - the number of nodes in the quadrature */

/*                output parameters: */

/*  ts - the nodes of the n-point gaussian quadrature */
/*  w - the weights of the n-point gaussian quadrature */

/*       . . . construct the array of initial approximations */
/*             to the roots of the n-th legendre polynomial */

    /* Parameter adjustments */
    --whts;
    --ts;

    /* Function Body */
    eps = 1e-14;
    zero = 0.;
    done = 1.;
    pi = atan(done) * 4;
    h__ = pi / (*n << 1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = ((i__ << 1) - 1) * h__;
	ts[*n - i__ + 1] = cos(t);
/* L1200: */
    }

/*         use newton to find all roots of the legendre polynomial */

    ts[*n / 2 + 1] = 0.;
    i__1 = *n / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	xk = ts[i__];
	ifout = 0;
	deltold = 1.;
	for (k = 1; k <= 10; ++k) {
	    legepol_sum__(&xk, n, &pol, &der, &sum);
	    delta = -pol / der;
	    xk += delta;
	    if (abs(delta) < eps) {
		++ifout;
	    }

	    if (ifout == 3) {
		goto L1600;
	    }
/* L1400: */
	}
L1600:
	ts[i__] = xk;
	ts[*n - i__ + 1] = -xk;
/* L2000: */
    }

/*        construct the weights via the orthogonality relation */

    if (*ifwhts == 0) {
	return 0;
    }

    i__1 = (*n + 1) / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	legepol_sum__(&ts[i__], n, &pol, &der, &sum);
	whts[i__] = 1 / sum;
	whts[*n - i__ + 1] = whts[i__];
/* L2400: */
    }

    return 0;
} /* legewhts_ */






/* Subroutine */ int legepol_sum__(doublereal *x, integer *n, doublereal *pol,
	 doublereal *der, doublereal *sum)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k;
    static doublereal pk, pkm1, pkp1, done;


    done = 1.;
    *sum = 0.;

    pkm1 = 1.;
    pk = *x;
/* Computing 2nd power */
    d__1 = pkm1;
    *sum += d__1 * d__1 / 2;
/* Computing 2nd power */
    d__1 = pk;
    *sum += d__1 * d__1 * (done / 2 + 1);

    pk = 1.;
    pkp1 = *x;

/*        if n=0 or n=1 - exit */

    if (*n >= 2) {
	goto L1200;
    }
    *sum = 0.;

    *pol = 1.;
    *der = 0.;
/* Computing 2nd power */
    d__1 = *pol;
    *sum += d__1 * d__1 / 2;
    if (*n == 0) {
	return 0;
    }

    *pol = *x;
    *der = 1.;
/* Computing 2nd power */
    d__1 = *pol;
    *sum += d__1 * d__1 * (done / 2 + 1);
    return 0;
L1200:

/*       n is greater than 1. conduct recursion */

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	pkm1 = pk;
	pk = pkp1;
	pkp1 = (((k << 1) + 1) * *x * pk - k * pkm1) / (k + 1);
/* Computing 2nd power */
	d__1 = pkp1;
	*sum += d__1 * d__1 * (k + 1 + done / 2);
/* L2000: */
    }

/*        calculate the derivative */

    *pol = pkp1;
/* Computing 2nd power */
    d__1 = *x;
    *der = *n * (*x * pkp1 - pk) / (d__1 * d__1 - 1);
    return 0;
} /* legepol_sum__ */






/* Subroutine */ int legepol_(doublereal *x, integer *n, doublereal *pol, 
	doublereal *der)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k;
    static doublereal pk, pkm1, pkp1;


    pkm1 = 1.;
    pk = *x;

    pk = 1.;
    pkp1 = *x;

/*        if n=0 or n=1 - exit */

    if (*n >= 2) {
	goto L1200;
    }
    *pol = 1.;
    *der = 0.;
    if (*n == 0) {
	return 0;
    }

    *pol = *x;
    *der = 1.;
    return 0;
L1200:

/*       n is greater than 1. conduct recursion */

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	pkm1 = pk;
	pk = pkp1;
	pkp1 = (((k << 1) + 1) * *x * pk - k * pkm1) / (k + 1);
/* L2000: */
    }

/*        calculate the derivative */

    *pol = pkp1;
/* Computing 2nd power */
    d__1 = *x;
    *der = *n * (*x * pkp1 - pk) / (d__1 * d__1 - 1);
    return 0;
} /* legepol_ */






/* Subroutine */ int prodend_(doublereal *x, doublereal *xs, integer *n, 
	integer *i__, doublereal *f)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer j;
    static doublereal dd, d10000;
    static integer large;
    static doublereal dlarge, dsmall;


/*      evaluate the product */

    /* Parameter adjustments */
    --xs;

    /* Function Body */
    *f = 1.;
    dlarge = 1e20;
    dsmall = *f / dlarge;

    large = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dd = abs(*f);
	if (dd > dsmall) {
	    goto L1200;
	}
	*f *= 10000;
	--large;
L1200:

	if (dd < dlarge) {
	    goto L1400;
	}
	*f /= 10000;
	++large;
L1400:
	if (j == *i__) {
	    goto L2000;
	}
	*f = *f * (*x - xs[j]) / (xs[*i__] - xs[j]);
L2000:
	;
    }
    d10000 = 1e4;
    *f *= pow_di(&d10000, &large);
/* Computing 2nd power */
    d__1 = *f;
    *f = d__1 * d__1 * (*x - xs[*i__]);
    return 0;
} /* prodend_ */






/* Subroutine */ int legepols_(doublereal *x, integer *n, doublereal *pols)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    static doublereal pk, pkm1, pkp1;


    /* Parameter adjustments */
    --pols;

    /* Function Body */
    pkm1 = 1.;
    pk = *x;

    pk = 1.;
    pkp1 = *x;


/*        if n=0 or n=1 - exit */

    if (*n >= 2) {
	goto L1200;
    }
    pols[1] = 1.;
    if (*n == 0) {
	return 0;
    }

    pols[2] = *x;
    return 0;
L1200:

    pols[1] = 1.;
    pols[2] = *x;

/*       n is greater than 2. conduct recursion */

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	pkm1 = pk;
	pk = pkp1;
	pkp1 = (((k << 1) + 1) * *x * pk - k * pkm1) / (k + 1);
	pols[k + 2] = pkp1;
/* L2000: */
    }

    return 0;
} /* legepols_ */






/* Subroutine */ int legepolders_(doublereal *x, doublereal *vals, doublereal 
	*ders, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal pj, pjm1, pjm2, derj, done, derjm1, derjm2;


/*     This subroutine computes the values and the derivatives */
/*     of n+1 first Legendre polynomials at the point x */
/*     in interval [-1,1]. */

/*                input parameters: */

/*     X = evaluation point */
/*     N  = order of expansion */
/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*     VALs = computed values of Legendre polynomials */
/*     ders = computed values of the derivatives */


    /* Parameter adjustments */
    --ders;
    --vals;

    /* Function Body */
    done = 1.;
    pjm2 = 1.;
    pjm1 = *x;
    derjm2 = 0.;
    derjm1 = 1.;

    vals[1] = 1.;
    ders[1] = 0.;

    vals[2] = *x;
    ders[2] = 1.;

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {

	pj = (((j << 1) - 1) * *x * pjm1 - (j - 1) * pjm2) / j;
	derj = ((j << 1) - 1) * (pjm1 + *x * derjm1) - (j - 1) * derjm2;

	derj /= j;
	vals[j + 1] = pj;
	ders[j + 1] = derj;

	pjm2 = pjm1;
	pjm1 = pj;
	derjm2 = derjm1;
	derjm1 = derj;
/* L600: */
    }

    return 0;
} /* legepolders_ */






/* Subroutine */ int legepls2_(doublereal *x, integer *n, doublereal *pols)
{
    /* Initialized data */

    static integer ifcalled = 0;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;
    static doublereal pk, pkm1, pkp1, done;
    static integer ninit;
    static doublereal pjcoefs1[2000], pjcoefs2[300];

    /* Parameter adjustments */
    --pols;

    /* Function Body */

/*        if need be - initialize the arrays pjcoefs1, pjcoefs2 */

    if (ifcalled == 1) {
	goto L1100;
    }

    done = 1.;
    ninit = 290;
    i__1 = ninit;
    for (j = 2; j <= i__1; ++j) {

	pjcoefs1[j - 1] = ((j << 1) - done) / j;
	pjcoefs2[j - 1] = -(j - done) / j;

/* L1050: */
    }

    ifcalled = 1;
L1100:
    pkm1 = 1.;
    pk = *x;

    pk = 1.;
    pkp1 = *x;

/*        if n=0 or n=1 - exit */

    if (*n >= 2) {
	goto L1200;
    }
    pols[1] = 1.;
    if (*n == 0) {
	return 0;
    }

    pols[2] = *x;
    return 0;
L1200:

    pols[1] = 1.;
    pols[2] = *x;

/*       n is greater than 2. conduct recursion */

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	pkm1 = pk;
	pk = pkp1;
	pkp1 = *x * pk * pjcoefs1[k] + pkm1 * pjcoefs2[k];
	pols[k + 2] = pkp1;
/* L2000: */
    }

    return 0;
} /* legepls2_ */






/* Subroutine */ int legeinmt_(integer *n, doublereal *ainte, doublereal *
	adiff, doublereal *x, doublereal *whts, doublereal *endinter, integer 
	*itype, doublereal *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iu, iv, iw, lu, lv, lw, ltot, ipolin, lpolin, ipolout, 
	    lpolout;
    extern /* Subroutine */ int legeinm0_(integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    ;



/*        for the user-specified n, this subroutine constructs */
/*        the matrices of spectral indefinite integration and/or */
/*        spectral differentiation on the n Gaussian nodes */
/*        on the interval [-1,1]. Actually, this is omnly a */
/*        memory management routine. All the actual work is done */
/*        by the subroutine legeinm0 (see) */

/*                           input parameters: */

/*  n - the number of Gaussian nodes on the interval [-1,1] */
/*  itype - the type of the calculation to be performed */
/*          EXPLANATION: */
/*       itype=1 means that only the matrix ainte will */
/*               be constructed */
/*       itype=2 means that only the matrix adiff will */
/*               be constructed */
/*       itype=3 means that both matrices ainte and adiff */
/*               will be constructed */

/*                           output paramaters: */

/*  ainte - the matrix of spectral indefinite integration on */
/*          the Gaussian nodes */
/*  adiff - the matrix of spectral differentiation on */
/*          the Gaussian nodes */
/*  x - the n Gaussian nodes on the intervl [-1,1] */
/*  whts - the n Gaussian weights on the interval [-1,1] */
/*  endinter - the interpolation coefficients converting the */
/*          values of a function at n Gaussian nodes into its */
/*          value at 1 (the right end of the interval) */

/*                           work arrays: */

/*  w - must be 3* n**2 + 2*n +50 *8 locations long */

/*        . . . allocate memory for the construction of the integrating */
/*              matrix */

    /* Parameter adjustments */
    --w;
    --endinter;
    --whts;
    --x;
    --adiff;
    --ainte;

    /* Function Body */
    ipolin = 1;
    lpolin = *n + 5;

    ipolout = ipolin + lpolin;
    lpolout = *n + 5;

    iu = ipolout + lpolout;
/* Computing 2nd power */
    i__1 = *n;
    lu = i__1 * i__1 + 1;

    iv = iu + lu;
/* Computing 2nd power */
    i__1 = *n;
    lv = i__1 * i__1 + 1;

    iw = iv + lv;
/* Computing 2nd power */
    i__1 = *n;
    lw = i__1 * i__1 + 1;

    ltot = iw + lw;

/*        construct the integrating matrix */

    legeinm0_(n, &ainte[1], &adiff[1], &w[ipolin], &w[ipolout], &x[1], &whts[
	    1], &w[iu], &w[iv], &w[iw], itype, &endinter[1]);

    return 0;
} /* legeinmt_ */






/* Subroutine */ int legeinm0_(integer *n, doublereal *ainte, doublereal *
	adiff, doublereal *polin, doublereal *polout, doublereal *x, 
	doublereal *whts, doublereal *u, doublereal *v, doublereal *w, 
	integer *itype, doublereal *endinter)
{
    /* System generated locals */
    integer ainte_dim1, ainte_offset, u_dim1, u_offset, v_dim1, v_offset, 
	    w_dim1, w_offset, adiff_dim1, adiff_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int legediff_(doublereal *, integer *, doublereal 
	    *), legeinte_(doublereal *, integer *, doublereal *), legeexps_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal d__;
    static integer i__, j, itype2;
    extern /* Subroutine */ int matmul_(doublereal *, doublereal *, 
	    doublereal *, integer *);


/*        for the user-specified n, this subroutine constructs */
/*        the matrices of spectral indefinite integration and/or */
/*        spectral differentiation on the n Gaussian nodes */
/*        on the interval [-1,1] */

/*                           input parameters: */

/*  n - the number of Gaussian nodes on the interval [-1,1] */
/*  itype - the type of the calculation to be performed */
/*          EXPLANATION: */
/*       itype=1 means that only the matrix ainte will */
/*               be constructed */
/*       itype=2 means that only the matrix adiff will */
/*               be constructed */
/*       itype=3 means that both matrices ainte and adiff */
/*               will be constructed */

/*                           output paramaters: */

/*  ainte - the matrix of spectral indefinite integration on */
/*          the Gaussian nodes */
/*  adiff - the matrix of spectral differentiation on */
/*          the Gaussian nodes */
/*  x - the n Gaussian nodes on the intervl [-1,1] */
/*  whts - the n Gaussian weights on the interval [-1,1] */

/*                           work arrays: */

/*  polin, polout - must be n+3 real *8 locations each */

/*  u, v, w - must be n**2+1 real *8 locations each */

/*        . . . construct the matrices of the forward and inverse */
/*              Legendre transforms */

    /* Parameter adjustments */
    w_dim1 = *n;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    u_dim1 = *n;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --whts;
    --x;
    --polout;
    --polin;
    adiff_dim1 = *n;
    adiff_offset = 1 + adiff_dim1;
    adiff -= adiff_offset;
    ainte_dim1 = *n;
    ainte_offset = 1 + ainte_dim1;
    ainte -= ainte_offset;
    --endinter;

    /* Function Body */
    itype2 = 2;
    legeexps_(&itype2, n, &x[1], &u[u_offset], &v[v_offset], &whts[1]);

/* ccc         call prin2('after legeexps, u=*',u,n*n) */

/*        if the user so requested, */
/*        construct the matrix converting the coefficients of */
/*        the Legendre series of a function into the coefficients */
/*        of the indefinite integral of that function */

    if (*itype == 2) {
	goto L2000;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n + 2;
	for (j = 1; j <= i__2; ++j) {
	    polin[j] = 0.;
/* L1200: */
	}

	polin[i__] = 1.;

	legeinte_(&polin[1], n, &polout[1]);

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    ainte[j + i__ * ainte_dim1] = polout[j];
/* L1400: */
	}

/* L1600: */
    }

/* ccc         call prin2('ainte initially is*',ainte,n*n) */

/*        multiply the three, obtaining the integrating matrix */

    matmul_(&ainte[ainte_offset], &u[u_offset], &w[w_offset], n);
    matmul_(&v[v_offset], &w[w_offset], &ainte[ainte_offset], n);

L2000:

/*        if the user so requested, */
/*        construct the matrix converting the coefficients of */
/*        the Legendre series of a function into the coefficients */
/*        of the derivative of that function */

    if (*itype == 1) {
	goto L3000;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n + 2;
	for (j = 1; j <= i__2; ++j) {
	    polin[j] = 0.;
/* L2200: */
	}

	polin[i__] = 1.;

	legediff_(&polin[1], n, &polout[1]);

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    adiff[j + i__ * adiff_dim1] = polout[j];
/* ccc        ainte(i,j)=polout(j) */
/* L2400: */
	}

/* L2600: */
    }

/* ccc         call prin2('adiff initially is*',adiff,n*n) */

/*        multiply the three, obtaining the integrating matrix */

    matmul_(&adiff[adiff_offset], &u[u_offset], &w[w_offset], n);
    matmul_(&v[v_offset], &w[w_offset], &adiff[adiff_offset], n);

L3000:

/*        construct the vector of interpolation coefficients */
/*        converting the values of a polynomial at the Gaussian */
/*        nodes into its value at the right end of the interval */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	d__ = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    d__ += u[j + i__ * u_dim1];
/* L3200: */
	}
	endinter[i__] = d__;
/* L3400: */
    }

    return 0;
} /* legeinm0_ */






/* Subroutine */ int legeinte_(doublereal *polin, integer *n, doublereal *
	polout)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal dd, sss;


/*       this subroutine computes the indefinite integral of the */
/*       legendre expansion polin getting the expansion polout */


/*                       input parameters: */

/*  polin - the legendre expansion to be integrated */
/*  n - the order of the expansion polin */
/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */
/*         also nothe that the order of the integrated expansion is */
/*         n+1 (who could think!) */

/*                       output parameters: */

/*  polout - the legendre expansion of the integral of the function */
/*         represented by the expansion polin */

    /* Parameter adjustments */
    --polout;
    --polin;

    /* Function Body */
    i__1 = *n + 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	polout[i__] = 0.;
/* L1200: */
    }

    i__1 = *n + 1;
    for (k = 2; k <= i__1; ++k) {
	j = k - 1;

/* ccc        polout(k+1)=polin(k)/(2*j+1)+polout(k+1) */
	polout[k + 1] = polin[k] / ((j << 1) + 1);
	polout[k - 1] = -polin[k] / ((j << 1) + 1) + polout[k - 1];

/* L2000: */
    }

    polout[2] = polin[1] + polout[2];

    dd = 0.;
    sss = -1.;
    i__1 = *n + 1;
    for (k = 2; k <= i__1; ++k) {

	dd += polout[k] * sss;
	sss = -sss;
/* L2200: */
    }

/* cc        call prin2('dd=*',dd,1) */
    polout[1] = -dd;

    return 0;
} /* legeinte_ */






/* Subroutine */ int legediff_(doublereal *polin, integer *n, doublereal *
	polout)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;
    static doublereal pk, pkm1, pkm2;


/*       this subroutine differentiates the legendre */
/*       expansion polin getting the expansion polout */


/*                       input parameters: */

/*  polin - the legendre expansion to be differentiated */
/*  n - the order of the expansion polin */
/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */
/*         also nothe that the order of the integrated expansion is */
/*         n+1 (who could think!) */

/*                       output parameters: */

/*  polout - the legendre expansion of the derivative of the function */
/*         represented by the expansion polin */

    /* Parameter adjustments */
    --polout;
    --polin;

    /* Function Body */
    i__1 = *n + 1;
    for (k = 1; k <= i__1; ++k) {
	polout[k] = 0.;
/* L1200: */
    }

    pk = polin[*n + 1];
    pkm1 = polin[*n];
    pkm2 = 0.;
    for (k = *n + 1; k >= 2; --k) {

	j = k - 1;

	polout[k - 1] = pk * ((j << 1) - 1);
	if (k >= 3) {
	    pkm2 = polin[k - 2] + pk;
	}

	pk = pkm1;
	pkm1 = pkm2;

/* L2000: */
    }
    return 0;
} /* legediff_ */






/* Subroutine */ int legefder_(doublereal *x, doublereal *val, doublereal *
	der, doublereal *pexp, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal pj, pjm1, pjm2, derj, done, derjm1, derjm2;


/*     This subroutine computes the value and the derivative */
/*     of a gaussian expansion with coefficients PEXP */
/*     at point X in interval [-1,1]. */

/*                input parameters: */

/*     X = evaluation point */
/*     PEXP = expansion coefficients */
/*     N  = order of expansion */
/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*     VAL = computed value */
/*     der = computed value of the derivative */


    /* Parameter adjustments */
    --pexp;

    /* Function Body */
    done = 1.;
    pjm2 = 1.;
    pjm1 = *x;
    derjm2 = 0.;
    derjm1 = 1.;

    *val = pexp[1] * pjm2 + pexp[2] * pjm1;
    *der = pexp[2];

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {

	pj = (((j << 1) - 1) * *x * pjm1 - (j - 1) * pjm2) / j;
	*val += pexp[j + 1] * pj;

	derj = ((j << 1) - 1) * (pjm1 + *x * derjm1) - (j - 1) * derjm2;

	derj /= j;
	*der += pexp[j + 1] * derj;

	pjm2 = pjm1;
	pjm1 = pj;
	derjm2 = derjm1;
	derjm1 = derj;
/* L600: */
    }

    return 0;
} /* legefder_ */






/* Subroutine */ int legefde2_(doublereal *x, doublereal *val, doublereal *
	der, doublereal *pexp, integer *n, doublereal *pjcoefs1, doublereal *
	pjcoefs2, integer *ninit)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ifcalled, j;
    static doublereal pj, pjm1, pjm2, derj, done, derjm1, derjm2;


/*     This subroutine computes the value and the derivative */
/*     of a gaussian expansion with coefficients PEXP */
/*     at point X in interval [-1,1]. */

/*                input parameters: */

/*  X - evaluation point */
/*  PEXP - expansion coefficients */
/*  N  - order of expansion */
/*  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call */
/*      on a previous call to this subroutine. Please note that this */
/*      is only an input parameter if the parameter ninit (see below) */
/*      has been set to 0; otherwise, these are output parameters */
/*  ninit - tells the subroutine whether and to what maximum order the */
/*       arrays coepnm1,coepnp1,coexpnp1 should be initialized. */
/*     EXPLANATION: The subroutine will initialize the first ninit */
/*       elements of each of the arrays pjcoefs1, pjcoefs2. On the first */
/*       call to this subroutine, ninit should be set to the maximum */
/*       order n for which this subroutine might have to be called; */
/*       on subsequent calls, ninit should be set to 0. PLEASE NOTE */
/*       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE */
/*       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE */
/*       SUBROUTINE LEGEEXE2. If these arrays have been initialized */
/*       by one of these two subroutines, they do not need to be */
/*       initialized by the other one. */

/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*  VAL - computed value */
/*  der - computed value of the derivative */


    /* Parameter adjustments */
    --pjcoefs2;
    --pjcoefs1;
    --pexp;

    /* Function Body */
    if (*ninit == 0) {
	goto L1400;
    }

    done = 1.;
    i__1 = *ninit;
    for (j = 2; j <= i__1; ++j) {

	pjcoefs1[j] = ((j << 1) - done) / j;
	pjcoefs2[j] = -(j - done) / j;

/* L1200: */
    }

    ifcalled = 1;
L1400:

    pjm2 = 1.;
    pjm1 = *x;
    derjm2 = 0.;
    derjm1 = 1.;

    *val = pexp[1] * pjm2 + pexp[2] * pjm1;
    *der = pexp[2];

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {

/* ccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j */

	pj = pjcoefs1[j] * *x * pjm1 + pjcoefs2[j] * pjm2;
	*val += pexp[j + 1] * pj;

/* ccc        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2 */
	derj = pjcoefs1[j] * (pjm1 + *x * derjm1) + pjcoefs2[j] * derjm2;
/* cc         call prin2('derj=*',derj,1) */
/* ccc        derj=derj/j */
	*der += pexp[j + 1] * derj;

	pjm2 = pjm1;
	pjm1 = pj;
	derjm2 = derjm1;
	derjm1 = derj;
/* L1600: */
    }

    return 0;
} /* legefde2_ */






/* Subroutine */ int legeexev_(doublereal *x, doublereal *val, doublereal *
	pexp, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal pj, der, pjm1, pjm2, done;


/*     This subroutine computes the value o a Legendre */
/*     expansion with coefficients PEXP at point X in interval [-1,1] */

/*                input parameters: */

/*     X = evaluation point */
/*     PEXP = expansion coefficients */
/*     N  = order of expansion */
/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*     VAL = computed value */

    /* Parameter adjustments */
    --pexp;

    /* Function Body */
    done = 1.;
    pjm2 = 1.;
    pjm1 = *x;

    *val = pexp[1] * pjm2 + pexp[2] * pjm1;
    der = pexp[2];

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {

	pj = (((j << 1) - 1) * *x * pjm1 - (j - 1) * pjm2) / j;
	*val += pexp[j + 1] * pj;

	pjm2 = pjm1;
	pjm1 = pj;
/* L600: */
    }

    return 0;
} /* legeexev_ */






/* Subroutine */ int legeexe2_(doublereal *x, doublereal *val, doublereal *
	pexp, integer *n, doublereal *pjcoefs1, doublereal *pjcoefs2, integer 
	*ninit)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ifcalled, j;
    static doublereal pj, der, pjm1, pjm2, done;


/*     This subroutine computes the value o a Legendre */
/*     expansion with coefficients PEXP at point X in interval [-1,1] */

/*                input parameters: */

/*     X = evaluation point */
/*     PEXP = expansion coefficients */
/*     N  = order of expansion */
/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*     VAL = computed value */

    /* Parameter adjustments */
    --pjcoefs2;
    --pjcoefs1;
    --pexp;

    /* Function Body */
    done = 1.;
    if (*ninit == 0) {
	goto L1400;
    }

    done = 1.;
    i__1 = *ninit;
    for (j = 2; j <= i__1; ++j) {

	pjcoefs1[j] = ((j << 1) - done) / j;
	pjcoefs2[j] = -(j - done) / j;

/* L1200: */
    }

    ifcalled = 1;
L1400:

    pjm2 = 1.;
    pjm1 = *x;

    *val = pexp[1] * pjm2 + pexp[2] * pjm1;
    der = pexp[2];

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {

	pj = pjcoefs1[j] * *x * pjm1 + pjcoefs2[j] * pjm2;
/* ccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j */
	*val += pexp[j + 1] * pj;

	pjm2 = pjm1;
	pjm1 = pj;
/* L600: */
    }

    return 0;
} /* legeexe2_ */






/* Subroutine */ int lematrin_(integer *n, integer *m, doublereal *xs, 
	doublereal *amatrint, doublereal *ts, doublereal *w)
{
    /* System generated locals */
    integer amatrint_dim1, amatrint_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, iu, iv, lu, icoefs, lcoefs, ifinit;
    extern /* Subroutine */ int levecin_(integer *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *);



/*        This subroutine constructs the matrix interpolating */
/*        functions from the n-point Gaussian grid on the interval [-1,1] */
/*        to an arbitrary m-point grid (the nodes of the latter are */
/*        user-provided) */

/*                 Input parameters: */

/*  n - the number of interpolation nodes */
/*  m - the number of nodes to which the functions will be interpolated */
/*  xs - the points at which the function is to be interpolated */

/*                  Output parameters: */

/*  amatrint - the m \times n matrix conerting the values of a function */
/*        at the n Legendre nodes into its values at m user-specified */
/*        (arbitrary) nodes */
/*  ts - the n Gaussian nodes on the interval [-1,1] */

/*                  Work arrays: */

/*  w - must be at least 2*n**2+n + 100 real *8 locations long */

    /* Parameter adjustments */
    amatrint_dim1 = *m;
    amatrint_offset = 1 + amatrint_dim1;
    amatrint -= amatrint_offset;
    --xs;
    --ts;
    --w;

    /* Function Body */
    icoefs = 1;
    lcoefs = *n + 2;

    iu = icoefs + lcoefs;
/* Computing 2nd power */
    i__1 = *n;
    lu = i__1 * i__1 + 10;

    iv = iu + lu;

    ifinit = 1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {

	levecin_(n, &xs[i__], &ts[1], &w[iu], &w[iv], &w[icoefs], &ifinit);

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    amatrint[i__ + j * amatrint_dim1] = w[j];
/* L1400: */
	}

	ifinit = 0;
/* L2000: */
    }

    return 0;
} /* lematrin_ */






/* Subroutine */ int levecin_(integer *n, doublereal *x, doublereal *ts, 
	doublereal *u, doublereal *v, doublereal *coefs, integer *ifinit)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int lematvec_(doublereal *, doublereal *, 
	    doublereal *, integer *), legepols_(doublereal *, integer *, 
	    doublereal *), legeexps_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer itype;


/*        This subroutine constructs the coefficients of the */
/*        standard interpolation formula connecting the values of a */
/*        function at n Gaussian nodes on the interval [a,b] with */
/*        its value at the point x \in R^1 */

/*                 Input parameters: */

/*  n - the number of interpolation nodes */
/*  x - the points at which the function is to be interpolated */
/*  ts - the n Gaussian nodes on the interval [-1,1]; please note that */
/*        it is an input parameter only if the parameter ifinit (see */
/*        below) has been set to 1; otherwise, it is an output parameter */
/*  u - the n*n matrix converting the  values at of a polynomial of order */
/*         n-1 at n legendre nodes into the coefficients of its */
/*        legendre expansion; please note that */
/*        it is an input parameter only if the parameter ifinit (see */
/*        below) has been set to 1; otherwise, it is an output parameter */
/*  ifinit - an integer parameter telling the subroutine whether it should */
/*        initialize the Legendre expander; */
/*     ifinit=1 will cause the subroutine to perform the initialization */
/*     ifinit=0 will cause the subroutine to  skip the initialization */

/*                  Output parameters: */

/*  coefs - the interpolation coefficients */

/*                 Work arrays: */

/*  v - must be at least n*n real *8 locations long */

/*       . . . construct the n Gausian nodes on the interval [-1,1]; */
/*             also the corresponding Gaussian expansion-evaluation */
/*             matrices */

    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    u_dim1 = *n;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --ts;
    --coefs;

    /* Function Body */
    itype = 2;
    if (*ifinit != 0) {
	legeexps_(&itype, n, &ts[1], &u[u_offset], &v[v_offset], &coefs[1]);
    }

/*       evaluate the n Legendre polynomials at the point where the */
/*       functions will have to be interpolated */

    i__1 = *n + 1;
    legepols_(x, &i__1, &v[v_offset]);

/*       apply the interpolation matrix to the ector of values */
/*       of polynomials from the right */

    lematvec_(&u[u_offset], &v[v_offset], &coefs[1], n);
    return 0;
} /* levecin_ */






/* Subroutine */ int lematvec_(doublereal *a, doublereal *x, doublereal *y, 
	integer *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static doublereal d__;
    static integer i__, j;


    /* Parameter adjustments */
    --y;
    --x;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__ = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    d__ += a[j + i__ * a_dim1] * x[j];
/* L1200: */
	}
	y[i__] = d__;
/* L1400: */
    }
    return 0;
} /* lematvec_ */






/* Subroutine */ int matmul_0_(int n__, doublereal *a, doublereal *b, 
	doublereal *c__, integer *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;


    /* Parameter adjustments */
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *n;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    switch(n__) {
	case 1: goto L_matmua;
	}

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    d__ = 0.;
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		d__ += a[i__ + k * a_dim1] * b[k + j * b_dim1];
/* L1600: */
	    }
	    c__[i__ + j * c_dim1] = d__;
/* L1800: */
	}
/* L2000: */
    }
    return 0;





L_matmua:
/* cc          call prin2('in matmua, a=*',a,n**2) */
/* cc          call prin2('in matmua, b=*',b,n**2) */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    d__ = 0.;
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		d__ += a[i__ + k * a_dim1] * b[j + k * b_dim1];
/* L2600: */
	    }
	    c__[i__ + j * c_dim1] = d__;
/* L2800: */
	}
/* L3000: */
    }
/* cc          call prin2('exiting, c=*',c,n**2) */
    return 0;
} /* matmul_ */

/* Subroutine */ int matmul_(doublereal *a, doublereal *b, doublereal *c__, 
	integer *n)
{
    return matmul_0_(0, a, b, c__, n);
    }

/* Subroutine */ int matmua_(doublereal *a, doublereal *b, doublereal *c__, 
	integer *n)
{
    return matmul_0_(1, a, b, c__, n);
    }






/* Subroutine */ int legeodev_(doublereal *x, integer *nn, doublereal *coefs, 
	doublereal *val, integer *ninit, doublereal *coepnm1, doublereal *
	coepnp1, doublereal *coexpnp1)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer i__, n;
    static doublereal pi, x22;
    static integer nnn;
    static doublereal pip1, pip2, done;



/*       This subroutine evaluates at the point x a Legendre expansion */
/*       having only odd-numbered elements */

/*                  Input parameters: */

/*  x - point on the interval [-1,1] at which the Legendre expansion */
/*       is to be evaluated */
/*  nn - order of the expansion to be evaluated */
/*  coefs - odd-numbered coefficients of the Legendre expansion */
/*       to be evaluated at the point x (nn/2+2 of them things) */
/*  ninit - tells the subroutine whether and to what maximum order the */
/*       arrays coepnm1,coepnp1,coexpnp1 should be initialized. */
/*     EXPLANATION: The subroutine will initialize the first ninit/2+2 */
/*                  (or so) elements of each of the arrays  coepnm1, */
/*       coepnp1, coexpnp1. On the first call to this subroutine, ninit */
/*       should be set to the maximum order nn for which this subroutine */
/*       might have to be called; on subsequent calls, ninit should be */
/*       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE */
/*       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE */
/*       SUBROUTINE LEGEPODD. IF these arrays have been initialized */
/*       by one of these two subroutines, they do not need to be */
/*       initialized by the other one. */
/*  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long */
/*       each. Please note that these are input arrays only if ninit */
/*       (see above) has been set to 0; otherwise, these are output arrays. */

/*                  Output parameters: */

/*  val - the value at the point x of the Legendre expansion with */
/*       coefficients coefs (see above) */
/*    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x), */
/*       pols(2) = P_2(x), pols(3) =P_4(x),  etc. */
/*  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long */
/*       each. Please note that these are output parameters only if ninit */
/*       (see above) has not been set to 0; otherwise, these are input */
/*       parameters */


    /* Parameter adjustments */
    --coexpnp1;
    --coepnp1;
    --coepnm1;
    --coefs;

    /* Function Body */
    if (*ninit == 0) {
	goto L1400;
    }
    done = 1.;
    n = 0;
    i__ = 0;

    i__1 = *ninit;
    for (nnn = 2; nnn <= i__1; nnn += 2) {

	n += 2;
	++i__;

/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coepnm1[i__] = -(n * 5 + d__1 * d__1 * 7 + d__2 * (d__2 * d__2) * 2);
/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coepnp1[i__] = -(n * 24 + 9 + d__1 * d__1 * 18 + d__2 * (d__2 * d__2) 
		* 4);
/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coexpnp1[i__] = n * 46 + 15 + d__1 * d__1 * 36 + d__2 * (d__2 * d__2) 
		* 8;

	d__ = (n * done + 2) * (n * done + 3) * ((n << 1) * done + 1);
	coepnm1[i__] /= d__;
	coepnp1[i__] /= d__;
	coexpnp1[i__] /= d__;

/* L1200: */
    }

L1400:

/* Computing 2nd power */
    d__1 = *x;
    x22 = d__1 * d__1;

    pi = *x;
    pip1 = *x * (x22 * 2.5 - 1.5);

    *val = coefs[1] * pi + coefs[2] * pip1;
    i__1 = *nn / 2 - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	pip2 = coepnm1[i__] * pi + (coepnp1[i__] + coexpnp1[i__] * x22) * 
		pip1;

	*val += coefs[i__ + 2] * pip2;

	pi = pip1;
	pip1 = pip2;
/* L2000: */
    }

    return 0;
} /* legeodev_ */






/* Subroutine */ int legeevev_(doublereal *x, integer *nn, doublereal *coefs, 
	doublereal *val, integer *ninit, doublereal *coepnm1, doublereal *
	coepnp1, doublereal *coexpnp1)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer i__, n;
    static doublereal pi, x22;
    static integer nnn;
    static doublereal pip1, pip2, done;



/*       This subroutine evaluates at the point x a Legendre expansion */
/*       having only even-numbered elements */

/*                  Input parameters: */

/*  x - point on the interval [-1,1] at which the Legendre expansion */
/*       is to be evaluated */
/*  nn - order of the expansion to be evaluated */
/*  coefs - even-numbered coefficients of the Legendre expansion */
/*       to be evaluated at the point x (nn/2+2 of them things) */
/*  ninit - tells the subroutine whether and to what maximum order the */
/*       arrays coepnm1,coepnp1,coexpnp1 should be initialized. */
/*     EXPLANATION: The subroutine will initialize the first ninit/2+2 */
/*                  (or so) elements of each of the arrays  coepnm1, */
/*       coepnp1, coexpnp1. On the first call to this subroutine, ninit */
/*       should be set to the maximum order nn for which this subroutine */
/*       might have to be called; on subsequent calls, ninit should be */
/*       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE */
/*       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE */
/*       SUBROUTINE LEGEPEVEN. IF these aqrrays have been initialized */
/*       by one of these two subroutines, they do not need to be */
/*       initialized by the other one. */
/*  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long */
/*       each. Please note that these are input arrays only if ninit */
/*       (see above) has been set to 0; otherwise, these are output arrays. */

/*                  Output parameters: */

/*  val - the value at the point x of the Legendre expansion with */
/*       coefficients coefs (see above) */
/*    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x), */
/*       pols(2) = P_2(x), pols(3) =P_4(x),  etc. */
/*  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long */
/*       each. Please note that these are output parameters only if ninit */
/*       (see above) has not been set to 0; otherwise, these are input */
/*       parameters */

    /* Parameter adjustments */
    --coexpnp1;
    --coepnp1;
    --coepnm1;
    --coefs;

    /* Function Body */
    if (*ninit == 0) {
	goto L1400;
    }

    done = 1.;
    n = -1;
    i__ = 0;
    i__1 = *ninit;
    for (nnn = 1; nnn <= i__1; nnn += 2) {

	n += 2;
	++i__;

/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coepnm1[i__] = -(n * 5 + d__1 * d__1 * 7 + d__2 * (d__2 * d__2) * 2);
/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coepnp1[i__] = -(n * 24 + 9 + d__1 * d__1 * 18 + d__2 * (d__2 * d__2) 
		* 4);
/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coexpnp1[i__] = n * 46 + 15 + d__1 * d__1 * 36 + d__2 * (d__2 * d__2) 
		* 8;

	d__ = (n * done + 2) * (n * done + 3) * ((n << 1) * done + 1);
	coepnm1[i__] /= d__;
	coepnp1[i__] /= d__;
	coexpnp1[i__] /= d__;

/* L1200: */
    }

L1400:

/* Computing 2nd power */
    d__1 = *x;
    x22 = d__1 * d__1;

    pi = 1.;
    pip1 = x22 * 1.5 - .5;

    *val = coefs[1] + coefs[2] * pip1;

/*       n is greater than 2. conduct recursion */

    i__1 = *nn / 2 - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	pip2 = coepnm1[i__] * pi + (coepnp1[i__] + coexpnp1[i__] * x22) * 
		pip1;
	*val += coefs[i__ + 2] * pip2;

	pi = pip1;
	pip1 = pip2;

/* L2000: */
    }

    return 0;
} /* legeevev_ */






/* Subroutine */ int legepeven_(doublereal *x, integer *nn, doublereal *pols, 
	integer *ninit, doublereal *coepnm1, doublereal *coepnp1, doublereal *
	coexpnp1)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer i__, n;
    static doublereal x22;
    static integer nnn;
    static doublereal done;


/*       This subroutine evaluates even-numbered Legendre polynomials */
/*       of the argument x, up to order nn+1 */

/*                  Input parameters: */

/*  x - the argument for which the Legendre polynomials are */
/*       to be evaluated */
/*  nn - the maximum order for which the Legendre polynomials are */
/*       to be evaluated */
/*  ninit - tells the subroutine whether and to what maximum order the */
/*       arrays coepnm1,coepnp1,coexpnp1 should be initialized. */
/*     EXPLANATION: The subroutine ill initialize the first ninit/2+2 */
/*                  (or so) elements of each of the arrays  coepnm1, */
/*       coepnp1, coexpnp1. On the first call to this subroutine, ninit */
/*       should be set to the maximum order nn for which this subroutine */
/*       might have to be called; on subsequent calls, ninit should be */
/*       set to 0. */
/*  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long */
/*       each. Please note that these are input arrays only if ninit */
/*       (see above) has been set to 0; otherwise, these are output arrays. */
/*       PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE */
/*       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE */
/*       SUBROUTINE LEGEEVEV. IF these aqrrays have been initialized */
/*       by one of these two subroutines, they do not need to be */
/*       initialized by the other one. */

/*                  Output parameters: */

/*  pols - even-numbered Legendre polynomials of the input parameter x */
/*         (nn/2+2 of them things) */
/*    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x), */
/*       pols(2) = P_2(x), pols(3) =P_4(x),  etc. */
/*  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long */
/*       each. Please note that these are output parameters only if ninit */
/*       (see above) has not been set to 0; otherwise, these are input */
/*       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE */
/*       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE */
/*       SUBROUTINE LEGEEVEV. If these arrays have been initialized */
/*       by one of these two subroutines, they do not need to be */
/*       initialized by the other one. */


    /* Parameter adjustments */
    --coexpnp1;
    --coepnp1;
    --coepnm1;
    --pols;

    /* Function Body */
    if (*ninit == 0) {
	goto L1400;
    }

    done = 1.;
    n = -1;
    i__ = 0;
    i__1 = *ninit;
    for (nnn = 1; nnn <= i__1; nnn += 2) {

	n += 2;
	++i__;

/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coepnm1[i__] = -(n * 5 + d__1 * d__1 * 7 + d__2 * (d__2 * d__2) * 2);
/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coepnp1[i__] = -(n * 24 + 9 + d__1 * d__1 * 18 + d__2 * (d__2 * d__2) 
		* 4);
/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coexpnp1[i__] = n * 46 + 15 + d__1 * d__1 * 36 + d__2 * (d__2 * d__2) 
		* 8;

	d__ = (n * done + 2) * (n * done + 3) * ((n << 1) * done + 1);
	coepnm1[i__] /= d__;
	coepnp1[i__] /= d__;
	coexpnp1[i__] /= d__;

/* L1200: */
    }

L1400:

/* Computing 2nd power */
    d__1 = *x;
    x22 = d__1 * d__1;

    pols[1] = 1.;
    pols[2] = x22 * 1.5 - .5;

/*       n is greater than 2. conduct recursion */

    i__1 = *nn / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	pols[i__ + 2] = coepnm1[i__] * pols[i__] + (coepnp1[i__] + coexpnp1[
		i__] * x22) * pols[i__ + 1];

/* L2000: */
    }

    return 0;
} /* legepeven_ */






/* Subroutine */ int legepodd_(doublereal *x, integer *nn, doublereal *pols, 
	integer *ninit, doublereal *coepnm1, doublereal *coepnp1, doublereal *
	coexpnp1)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer i__, n;
    static doublereal x22;
    static integer nnn;
    static doublereal done;


/*       This subroutine evaluates odd-numbered Legendre polynomials */
/*       of the argument x, up to order nn+1 */

/*                  Input parameters: */

/*  x - the argument for which the Legendre polynomials are */
/*       to be evaluated */
/*  nn - the maximum order for which the Legendre polynomials are */
/*       to be evaluated */
/*  ninit - tells the subroutine whether and to what maximum order the */
/*       arrays coepnm1,coepnp1,coexpnp1 should be initialized. */
/*     EXPLANATION: The subroutine will initialize the first ninit/2+2 */
/*                  (or so) elements of each of the arrays  coepnm1, */
/*       coepnp1, coexpnp1. On the first call to this subroutine, ninit */
/*       should be set to the maximum order nn for which this subroutine */
/*       might have to be called; on subsequent calls, ninit should be */
/*       set to 0. */
/*  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long */
/*       each. Please note that these are input arrays only if ninit */
/*       (see above) has been set to 0; otherwise, these are output arrays. */

/*                  Output parameters: */

/*  pols - the odd-numbered Legendre polynomials of the input parameter x */
/*         (nn/2+2 of them things) */
/*    EXPLANATION: On exit from the subroutine, pols(1) = P_1(x), */
/*       pols(2) = P_3(x), pols(3) = P_5 (x), etc. */
/*  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long */
/*       each. Please note that these are output parameters only if ninit */
/*       (see above) has not been set to 0; otherwise, these are input */
/*       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS */
/*       SUBROUTINE ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES */
/*       SUSED BY THE UBROUTINE LEGEODEV. IF these arrays have been */
/*       initialized by one of these two subroutines, they do not need */
/*       to be initialized by the other one. */

    /* Parameter adjustments */
    --coexpnp1;
    --coepnp1;
    --coepnm1;
    --pols;

    /* Function Body */
    if (*ninit == 0) {
	goto L1400;
    }
    done = 1.;
    n = 0;
    i__ = 0;

    i__1 = *ninit;
    for (nnn = 2; nnn <= i__1; nnn += 2) {

	n += 2;
	++i__;

/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coepnm1[i__] = -(n * 5 + d__1 * d__1 * 7 + d__2 * (d__2 * d__2) * 2);
/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coepnp1[i__] = -(n * 24 + 9 + d__1 * d__1 * 18 + d__2 * (d__2 * d__2) 
		* 4);
/* Computing 2nd power */
	d__1 = n * done;
/* Computing 3rd power */
	d__2 = n * done;
	coexpnp1[i__] = n * 46 + 15 + d__1 * d__1 * 36 + d__2 * (d__2 * d__2) 
		* 8;

	d__ = (n * done + 2) * (n * done + 3) * ((n << 1) * done + 1);
	coepnm1[i__] /= d__;
	coepnp1[i__] /= d__;
	coexpnp1[i__] /= d__;

/* L1200: */
    }

L1400:

/* Computing 2nd power */
    d__1 = *x;
    x22 = d__1 * d__1;

    pols[1] = *x;
    pols[2] = *x * (x22 * 2.5 - 1.5);

    i__1 = *nn / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	pols[i__ + 2] = coepnm1[i__] * pols[i__] + (coepnp1[i__] + coexpnp1[
		i__] * x22) * pols[i__ + 1];

/* L2000: */
    }

    return 0;
} /* legepodd_ */






/* Subroutine */ int legefdeq_(doublereal *x, doublereal *val, doublereal *
	der, doublereal *coefs, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer k;
    static doublereal pk, pkm1, pkp1, derk, derkm1, derkp1;


/*     This subroutine computes the value and the derivative */
/*     of a Legendre Q-expansion with coefficients coefs */
/*     at point X in interval (-1,1); please note that this is */
/*     the evil twin of the subroutine legefder, evaluating the */
/*     proper (P-function) Legendre expansion */

/*                input parameters: */

/*  X = evaluation point */
/*  coefs = expansion coefficients */
/*  N  = order of expansion */

/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*     VAL = computed value */
/*     der = computed value of the derivative */


    /* Parameter adjustments */
    --coefs;

    /* Function Body */
    *val = 0.;
    *der = 0.;

    d__ = log((*x + 1) / (1 - *x)) / 2;
    pkm1 = d__;
    pk = d__ * *x - 1;

    pk = d__;
    pkp1 = d__ * *x - 1;
    derk = (1 / (*x + 1) + 1 / (1 - *x)) / 2;
    derkp1 = d__ + derk * *x;

    *val = coefs[1] * pk + coefs[2] * pkp1;
    *der = coefs[1] * derk + coefs[2] * derkp1;

/*        if n=0 or n=1 - exit */

    if (*n >= 2) {
	goto L1200;
    }

    if (*n == 0) {
	return 0;
    }

    return 0;
L1200:

/*       n is greater than 2. conduct recursion */

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	pkm1 = pk;
	pk = pkp1;

	pkp1 = (((k << 1) + 1) * *x * pk - k * pkm1) / (k + 1);

	derkm1 = derk;
	derk = derkp1;

	derkp1 = (((k << 1) + 1) * pk + ((k << 1) + 1) * *x * derk - k * 
		derkm1) / (k + 1);

	*val += coefs[k + 2] * pkp1;
	*der += coefs[k + 2] * derkp1;

/* L2000: */
    }

    return 0;
} /* legefdeq_ */






/* Subroutine */ int legeqs_(doublereal *x, integer *n, doublereal *pols, 
	doublereal *ders)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer k;
    static doublereal pk, pkm1, pkp1, derk, derkm1, derkp1;


/*       This subroutine calculates the values and derivatives of */
/*       a bunch of Legendre Q-functions at the user-specified point */
/*       x on the interval (-1,1) */

/*                     Input parameters: */

/*  x - the point on the interval [-1,1] where the Q-functions and */
/*       their derivatives are to be evaluated */
/*  n - the highest order for which the functions are to be evaluated */

/*                     Output parameters: */

/*  pols - the values of the Q-functions (the evil twins of the */
/*       Legeendre polynomials) at the point x (n+1 of them things) */
/*  ders - the derivatives of the Q-functions (the evil twins of the */
/*       Legeendre polynomials) at the point x (n+1 of them things) */


    /* Parameter adjustments */
    --ders;
    --pols;

    /* Function Body */
    d__ = log((*x + 1) / (1 - *x)) / 2;
    pkm1 = d__;
    pk = d__ * *x - 1;

    pk = d__;
    pkp1 = d__ * *x - 1;
    derk = (1 / (*x + 1) + 1 / (1 - *x)) / 2;
    derkp1 = d__ + derk * *x;

/*        if n=0 or n=1 - exit */

    if (*n >= 2) {
	goto L1200;
    }
    pols[1] = pk;
    ders[1] = derk;
    if (*n == 0) {
	return 0;
    }

    pols[2] = pkp1;
    ders[2] = derkp1;
    return 0;
L1200:

    pols[1] = pk;
    pols[2] = pkp1;

/*       n is greater than 2. conduct recursion */

    ders[1] = derk;
    ders[2] = derkp1;

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	pkm1 = pk;
	pk = pkp1;

	pkp1 = (((k << 1) + 1) * *x * pk - k * pkm1) / (k + 1);
	pols[k + 2] = pkp1;

	derkm1 = derk;
	derk = derkp1;

	derkp1 = (((k << 1) + 1) * pk + ((k << 1) + 1) * *x * derk - k * 
		derkm1) / (k + 1);
	ders[k + 2] = derkp1;
/* L2000: */
    }

    return 0;
} /* legeqs_ */






/* Subroutine */ int legeq_(doublereal *x, integer *n, doublereal *pol, 
	doublereal *der)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer k;
    static doublereal pk, pkm1, pkp1;


/*       This subroutine calculates the value and derivative of */
/*       a Legendre Q-function at the user-specified point */
/*       x on the interval (-1,1) */


/*                     Input parameters: */

/*  x - the point on the interval [-1,1] where the Q-functions and */
/*       their derivatives are to be evaluated */
/*  n - the order for which the function is to be evaluated */

/*                     Output parameters: */

/*  pol - the value of the n-th Q-function (the evil twin of the */
/*       Legeendre polynomial) at the point x */
/*  ders - the derivatives of the Q-function at the point x */


    d__ = log((*x + 1) / (1 - *x)) / 2;
    pk = d__;
    pkp1 = d__ * *x - 1;

/*        if n=0 or n=1 - exit */

    if (*n >= 2) {
	goto L1200;
    }
    *pol = d__;
    *der = (1 / (*x + 1) + 1 / (1 - *x)) / 2;
    if (*n == 0) {
	return 0;
    }

    *pol = pkp1;
    *der = d__ + *der * *x;
    return 0;
L1200:

/*       n is greater than 1. conduct recursion */

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	pkm1 = pk;
	pk = pkp1;
	pkp1 = (((k << 1) + 1) * *x * pk - k * pkm1) / (k + 1);
/* L2000: */
    }

/*        calculate the derivative */

    *pol = pkp1;
/* Computing 2nd power */
    d__1 = *x;
    *der = *n * (*x * pkp1 - pk) / (d__1 * d__1 - 1);
    return 0;
} /* legeq_ */






/* Subroutine */ int clegeq_(doublecomplex *x, integer *n, doublecomplex *pol,
	 doublecomplex *der)
{
    /* Initialized data */

    static doublecomplex ima = {0.,1.};

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double atan(doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_log(
	    doublecomplex *, doublecomplex *), pow_zi(doublecomplex *, 
	    doublecomplex *, integer *);

    /* Local variables */
    static doublecomplex d__;
    static integer k;
    static doublereal d3, pi;
    static doublecomplex pk, pkm1, pkp1;
    static doublereal done;


/*       This subroutine calculates the value and derivative of */
/*       a Legendre Q-function at the user-specified point */
/*       x on the interval (-1,1) */


/*                     Input parameters: */

/*  x - the point on the interval [-1,1] where the Q-functions and */
/*       their derivatives are to be evaluated */
/*  n - the order for which the function is to be evaluated */

/*                     Output parameters: */

/*  pol - the value of the n-th Q-function (the evil twin of the */
/*       Legeendre polynomial) at the point x */
/*  ders - the derivatives of the Q-function at the point x */



    done = 1.;
    pi = atan(done) * 4;

    z__3.r = x->r + 1, z__3.i = x->i;
    z__4.r = 1 - x->r, z__4.i = -x->i;
    z_div(&z__2, &z__3, &z__4);
    z_log(&z__1, &z__2);
    d__.r = z__1.r, d__.i = z__1.i;
    d__1 = 2.;
    z__1.r = d__.r / d__1, z__1.i = d__.i / d__1;
    d__.r = z__1.r, d__.i = z__1.i;

    z__2.r = -ima.r, z__2.i = -ima.i;
    z__1.r = z__2.r * x->r - z__2.i * x->i, z__1.i = z__2.r * x->i + z__2.i * 
	    x->r;
    d3 = z__1.r;
    if (d3 > 0.) {
	d__1 = pi / 2;
	z__2.r = d__1 * ima.r, z__2.i = d__1 * ima.i;
	z__1.r = d__.r - z__2.r, z__1.i = d__.i - z__2.i;
	d__.r = z__1.r, d__.i = z__1.i;
    }
    if (d3 < 0.) {
	d__1 = pi / 2;
	z__2.r = d__1 * ima.r, z__2.i = d__1 * ima.i;
	z__1.r = d__.r + z__2.r, z__1.i = d__.i + z__2.i;
	d__.r = z__1.r, d__.i = z__1.i;
    }

    pk.r = d__.r, pk.i = d__.i;
    z__2.r = d__.r * x->r - d__.i * x->i, z__2.i = d__.r * x->i + d__.i * 
	    x->r;
    z__1.r = z__2.r - 1, z__1.i = z__2.i;
    pkp1.r = z__1.r, pkp1.i = z__1.i;

/*        if n=0 or n=1 - exit */

    if (*n >= 2) {
	goto L1200;
    }
    pol->r = d__.r, pol->i = d__.i;
    z__4.r = x->r + 1, z__4.i = x->i;
    z_div(&z__3, &c_b125, &z__4);
    z__6.r = 1 - x->r, z__6.i = -x->i;
    z_div(&z__5, &c_b125, &z__6);
    z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
    d__1 = 2.;
    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
    der->r = z__1.r, der->i = z__1.i;

    if (*n == 0) {
	return 0;
    }

    pol->r = pkp1.r, pol->i = pkp1.i;
    z__2.r = der->r * x->r - der->i * x->i, z__2.i = der->r * x->i + der->i * 
	    x->r;
    z__1.r = d__.r + z__2.r, z__1.i = d__.i + z__2.i;
    der->r = z__1.r, der->i = z__1.i;
    return 0;
L1200:

/*       n is greater than 1. conduct recursion */

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	pkm1.r = pk.r, pkm1.i = pk.i;
	pk.r = pkp1.r, pk.i = pkp1.i;
	i__2 = (k << 1) + 1;
	d__1 = (doublereal) i__2;
	z__4.r = d__1 * x->r, z__4.i = d__1 * x->i;
	z__3.r = z__4.r * pk.r - z__4.i * pk.i, z__3.i = z__4.r * pk.i + 
		z__4.i * pk.r;
	d__2 = (doublereal) k;
	z__5.r = d__2 * pkm1.r, z__5.i = d__2 * pkm1.i;
	z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	i__3 = k + 1;
	d__3 = (doublereal) i__3;
	z__1.r = z__2.r / d__3, z__1.i = z__2.i / d__3;
	pkp1.r = z__1.r, pkp1.i = z__1.i;
/* L2000: */
    }

/*        calculate the derivative */

    pol->r = pkp1.r, pol->i = pkp1.i;
    z__4.r = x->r * pkp1.r - x->i * pkp1.i, z__4.i = x->r * pkp1.i + x->i * 
	    pkp1.r;
    z__3.r = z__4.r - pk.r, z__3.i = z__4.i - pk.i;
    d__1 = (doublereal) (*n);
    z__2.r = d__1 * z__3.r, z__2.i = d__1 * z__3.i;
    pow_zi(&z__6, x, &c__2);
    z__5.r = z__6.r - 1, z__5.i = z__6.i;
    z_div(&z__1, &z__2, &z__5);
    der->r = z__1.r, der->i = z__1.i;
    return 0;
} /* clegeq_ */






/* Subroutine */ int legecfde_(doublereal *x, doublecomplex *val, 
	doublecomplex *der, doublecomplex *pexp, integer *n)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer j;
    static doublereal pj, pjm1, pjm2, derj, done, derjm1, derjm2;


/*     This subroutine computes the value and the derivative */
/*     of a gaussian expansion with complex coefficients PEXP */
/*     at point X in interval [-1,1]. */

/*                input parameters: */

/*     X = evaluation point */
/*     PEXP = expansion coefficients */
/*     N  = order of expansion */
/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*     VAL = computed value */
/*     der = computed value of the derivative */


    /* Parameter adjustments */
    --pexp;

    /* Function Body */
    done = 1.;
    pjm2 = 1.;
    pjm1 = *x;
    derjm2 = 0.;
    derjm1 = 1.;

    z__2.r = pjm2 * pexp[1].r, z__2.i = pjm2 * pexp[1].i;
    z__3.r = pjm1 * pexp[2].r, z__3.i = pjm1 * pexp[2].i;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    val->r = z__1.r, val->i = z__1.i;
    der->r = pexp[2].r, der->i = pexp[2].i;

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {

	pj = (((j << 1) - 1) * *x * pjm1 - (j - 1) * pjm2) / j;
	i__2 = j + 1;
	z__2.r = pj * pexp[i__2].r, z__2.i = pj * pexp[i__2].i;
	z__1.r = val->r + z__2.r, z__1.i = val->i + z__2.i;
	val->r = z__1.r, val->i = z__1.i;

	derj = ((j << 1) - 1) * (pjm1 + *x * derjm1) - (j - 1) * derjm2;

	derj /= j;
	i__2 = j + 1;
	z__2.r = derj * pexp[i__2].r, z__2.i = derj * pexp[i__2].i;
	z__1.r = der->r + z__2.r, z__1.i = der->i + z__2.i;
	der->r = z__1.r, der->i = z__1.i;

	pjm2 = pjm1;
	pjm1 = pj;
	derjm2 = derjm1;
	derjm1 = derj;
/* L600: */
    }

    return 0;
} /* legecfde_ */






/* Subroutine */ int legecfd2_(doublereal *x, doublecomplex *val, 
	doublecomplex *der, doublecomplex *pexp, integer *n, doublereal *
	pjcoefs1, doublereal *pjcoefs2, integer *ninit)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer ifcalled, j;
    static doublereal pj, pjm1, pjm2, derj, done, derjm1, derjm2;


/*     This subroutine computes the value and the derivative */
/*     of a Legendre expansion with complex coefficients PEXP */
/*     at point X in interval [-1,1]. */

/*                input parameters: */

/*  X - evaluation point */
/*  PEXP - expansion coefficients */
/*  N  - order of expansion */
/*  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call */
/*      on a previous call to this subroutine. Please note that this */
/*      is only an input parameter if the parameter ninit (see below) */
/*      has been set to 0; otherwise, these are output parameters */
/*  ninit - tells the subroutine whether and to what maximum order the */
/*       arrays coepnm1,coepnp1,coexpnp1 should be initialized. */
/*     EXPLANATION: The subroutine will initialize the first ninit */
/*       elements of each of the arrays pjcoefs1, pjcoefs2. On the first */
/*       call to this subroutine, ninit should be set to the maximum */
/*       order n for which this subroutine might have to be called; */
/*       on subsequent calls, ninit should be set to 0. PLEASE NOTE */
/*       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE */
/*       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE */
/*       SUBROUTINE LEGEEXE2. If these arrays have been initialized */
/*       by one of these two subroutines, they do not need to be */
/*       initialized by the other one. */

/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*  VAL - computed value */
/*  der - computed value of the derivative */

    /* Parameter adjustments */
    --pjcoefs2;
    --pjcoefs1;
    --pexp;

    /* Function Body */
    if (*ninit == 0) {
	goto L1400;
    }

    done = 1.;
    i__1 = *ninit;
    for (j = 2; j <= i__1; ++j) {

	pjcoefs1[j] = ((j << 1) - done) / j;
	pjcoefs2[j] = -(j - done) / j;

/* L1200: */
    }

    ifcalled = 1;
L1400:

    pjm2 = 1.;
    pjm1 = *x;
    derjm2 = 0.;
    derjm1 = 1.;

    z__2.r = pjm2 * pexp[1].r, z__2.i = pjm2 * pexp[1].i;
    z__3.r = pjm1 * pexp[2].r, z__3.i = pjm1 * pexp[2].i;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    val->r = z__1.r, val->i = z__1.i;
    der->r = pexp[2].r, der->i = pexp[2].i;

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {

/* ccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j */

	pj = pjcoefs1[j] * *x * pjm1 + pjcoefs2[j] * pjm2;
	i__2 = j + 1;
	z__2.r = pj * pexp[i__2].r, z__2.i = pj * pexp[i__2].i;
	z__1.r = val->r + z__2.r, z__1.i = val->i + z__2.i;
	val->r = z__1.r, val->i = z__1.i;

/* ccc        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2 */
	derj = pjcoefs1[j] * (pjm1 + *x * derjm1) + pjcoefs2[j] * derjm2;
/* cc         call prin2('derj=*',derj,1) */
/* ccc        derj=derj/j */
	i__2 = j + 1;
	z__2.r = derj * pexp[i__2].r, z__2.i = derj * pexp[i__2].i;
	z__1.r = der->r + z__2.r, z__1.i = der->i + z__2.i;
	der->r = z__1.r, der->i = z__1.i;

	pjm2 = pjm1;
	pjm1 = pj;
	derjm2 = derjm1;
	derjm1 = derj;
/* L1600: */
    }

    return 0;
} /* legecfd2_ */






/* Subroutine */ int legecva2_(doublereal *x, doublecomplex *val, 
	doublecomplex *pexp, integer *n, doublereal *pjcoefs1, doublereal *
	pjcoefs2, integer *ninit)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer ifcalled, j;
    static doublereal pj, pjm1, pjm2, done;


/*     This subroutine computes the value of a Legendre expansion */
/*     with complex coefficients PEXP at point X in interval [-1,1]. */

/*                input parameters: */

/*  X - evaluation point */
/*  PEXP - expansion coefficients */
/*  N  - order of expansion */
/*  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call */
/*      on a previous call to this subroutine. Please note that this */
/*      is only an input parameter if the parameter ninit (see below) */
/*      has been set to 0; otherwise, these are output parameters */
/*  ninit - tells the subroutine whether and to what maximum order the */
/*       arrays coepnm1,coepnp1,coexpnp1 should be initialized. */
/*     EXPLANATION: The subroutine will initialize the first ninit */
/*       elements of each of the arrays pjcoefs1, pjcoefs2. On the first */
/*       call to this subroutine, ninit should be set to the maximum */
/*       order n for which this subroutine might have to be called; */
/*       on subsequent calls, ninit should be set to 0. PLEASE NOTE */
/*       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE */
/*       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE */
/*       SUBROUTINE LEGEEXE2. If these arrays have been initialized */
/*       by one of these two subroutines, they do not need to be */
/*       initialized by the other one. */

/*   IMPORTANT NOTE: n is {\bf the order of the expansion, which is */
/*         one less than the number of terms in the expansion!!} */

/*                output parameters: */

/*  VAL - computed value */

    /* Parameter adjustments */
    --pjcoefs2;
    --pjcoefs1;
    --pexp;

    /* Function Body */
    if (*ninit == 0) {
	goto L1400;
    }

    done = 1.;
    i__1 = *ninit;
    for (j = 2; j <= i__1; ++j) {

	pjcoefs1[j] = ((j << 1) - done) / j;
	pjcoefs2[j] = -(j - done) / j;

/* L1200: */
    }

    ifcalled = 1;
L1400:

    pjm2 = 1.;
    pjm1 = *x;

    z__2.r = pjm2 * pexp[1].r, z__2.i = pjm2 * pexp[1].i;
    z__3.r = pjm1 * pexp[2].r, z__3.i = pjm1 * pexp[2].i;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    val->r = z__1.r, val->i = z__1.i;

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {

	pj = pjcoefs1[j] * *x * pjm1 + pjcoefs2[j] * pjm2;
	i__2 = j + 1;
	z__2.r = pj * pexp[i__2].r, z__2.i = pj * pexp[i__2].i;
	z__1.r = val->r + z__2.r, z__1.i = val->i + z__2.i;
	val->r = z__1.r, val->i = z__1.i;

	pjm2 = pjm1;
	pjm1 = pj;
/* L1600: */
    }

    return 0;
} /* legecva2_ */






/* Subroutine */ int legerts_(integer *itype, integer *n, doublereal *ts, 
	doublereal *whts)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal);

    /* Local variables */
    extern /* Subroutine */ int legetayl_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *);
    static doublereal d__, h__;
    static integer i__, k;
    static doublereal t, d2;
    static integer n2;
    static doublereal x0, x1;
    static integer ii, kk;
    static doublereal pi, der, pol, der3, pol3, done;
    static integer ifodd, ifstop;
    extern /* Subroutine */ int legepol_(doublereal *, integer *, doublereal *
	    , doublereal *);


/*        This subroutine constructs the Gaussian quadrature */
/*        or order n. Its claim to fame is the fact that the */
/*        cost of the calculation is proportional to n; in */
/*        practice, with n=10 000 the calculation is more or */
/*        less instantaneous */

/*                 Input parameters: */

/*  itype - the type of calculation desired: */
/*     itype=1 will cause both the roots and the weights to be returned */
/*     itype=0 will cause only the roots to be returned */
/*  n - the number of nodes to be returned */

/*                 Output parameters: */

/*  ts - the n Gaussian nodes on the interval [-1,1] */
/*  whts - the n Gaussian weights on the interval [-1,1] */


/*        . . . determine the number of Taylor coefficients */
/*              to be used */

    /* Parameter adjustments */
    --whts;
    --ts;

    /* Function Body */
    k = 30;
    d__ = 1.;
    d2 = d__ + 1e-24;
    if (d2 != d__) {
	k = 54;
    }

/*       . . . construct the array of initial approximations */
/*             to the roots of the n-th legendre polynomial */

    i__ = *n / 2;
    ifodd = *n - (i__ << 1);

    done = 1.;
    pi = atan(done) * 4;
    h__ = pi / (*n << 1);
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	if (i__ < *n / 2 + 1) {
	    goto L1100;
	}
	++ii;
	t = ((i__ << 1) - 1) * h__;
	ts[ii] = -cos(t);
L1100:
	;
    }

/*       starting from the center, find roots one after another */

    pol = 1.;
    der = 0.;

    x0 = 0.;
    legepol_(&x0, n, &pol, &der);
    x1 = ts[1];

    n2 = (*n + 1) / 2;

    pol3 = pol;
    der3 = der;
    i__1 = n2;
    for (kk = 1; kk <= i__1; ++kk) {

	if (ifodd == 1 && kk == 1) {
	    ts[kk] = x0;
	    whts[kk] = der;
	    x0 = x1;
	    x1 = ts[kk + 1];
	    pol3 = pol;
	    der3 = der;
	    goto L2000;
	}

/*        conduct newton */

	ifstop = 0;
	for (i__ = 1; i__ <= 10; ++i__) {

	    h__ = x1 - x0;
	    legetayl_(&pol3, &der3, &x0, &h__, n, &k, &pol, &der);

	    x1 -= pol / der;

	    if (abs(pol) < 1e-12) {
		++ifstop;
	    }
	    if (ifstop == 3) {
		goto L1600;
	    }

/* L1400: */
	}
L1600:

	ts[kk] = x1;
	if (*itype > 0) {
	    whts[kk] = der;
	}

	x0 = x1;
	x1 = ts[kk + 1];
	pol3 = pol;
	der3 = der;
L2000:
	;
    }

/*        put the obtained roots in the proper order */

    for (i__ = (*n + 1) / 2; i__ >= 1; --i__) {

	ts[i__ + *n / 2] = ts[i__];
/* L2200: */
    }

    i__1 = *n / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	ts[i__] = -ts[*n - i__ + 1];
/* L2400: */
    }

    if (*itype <= 0) {
	return 0;
    }

/*        put the obtained roots in the proper order */

    for (i__ = (*n + 1) / 2; i__ >= 1; --i__) {

	whts[i__ + *n / 2] = whts[i__];
/* L2600: */
    }

    i__1 = *n / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

	whts[i__] = whts[*n - i__ + 1];
/* L2800: */
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

/* Computing 2nd power */
	d__1 = ts[i__];
/* Computing 2nd power */
	d__2 = whts[i__];
	whts[i__] = 2 / (1 - d__1 * d__1) / (d__2 * d__2);
/* L3600: */
    }

    return 0;
} /* legerts_ */






/* Subroutine */ int legetayl_(doublereal *pol, doublereal *der, doublereal *
	x, doublereal *h__, integer *n, integer *k, doublereal *sum, 
	doublereal *sumder)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer i__;
    static doublereal q0, q1, q2, qi, qip1, qip2, done;


/*        This subroutine evaluates the Taylor series for the */
/*        Legendre polynomial and its derivative at the point */
/*        x+h, starting with the value of the polynomial at the */
/*        point x, and the value of the derivative of that */
/*        polynomial. It uses the obvious three-term recursion */
/*        for the derivatives of Legendre polynomials. */

/*                 Input parameters: */

/*  pol - the value of the polynomial at the pount x */
/*  der - the derivative of the polynomial at the pount x */
/*  x - the point where the values pol, der are specified */
/*  h - the polynomial and its derivative will be evaluated at */
/*        the point x+h */
/*  n - the order of the Legendre polynomial */
/*  k - the order of the Taylor series to be used */

/*                 Output parameters: */

/*  sum - the value of P_n at the point x+h */
/*  sumder - the value of the derivative of P_n at the point x+h */

/*        . . . evaluate the derivatives of P_n scaled by h^n/n!, */
/*              and sum the taylor series for P_n and its */
/*              derivative */

    done = 1.;
    q0 = *pol;
    q1 = *der * *h__;
/* Computing 2nd power */
    d__1 = *x;
    q2 = (*x * 2 * *der - *n * (*n + done) * *pol) / (1 - d__1 * d__1);
/* Computing 2nd power */
    d__1 = *h__;
    q2 = q2 * (d__1 * d__1) / 2;

    *sum = q0 + q1 + q2;
    *sumder = q1 / *h__ + q2 * 2 / *h__;

    if (*k <= 2) {
	return 0;
    }

    qi = q1;
    qip1 = q2;

    i__1 = *k - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {

/* Computing 2nd power */
	i__2 = i__ + 1;
	d__ = *x * 2 * (i__2 * i__2) / *h__ * qip1 - (*n * (*n + done) - i__ *
		 (i__ + 1)) * qi;
/* Computing 2nd power */
	d__1 = *h__;
/* Computing 2nd power */
	d__2 = *x;
	d__ = d__ / (i__ + 1) / (i__ + 2) * (d__1 * d__1) / (1 - d__2 * d__2);
	qip2 = d__;

	*sum += qip2;
	*sumder += d__ * (i__ + 2) / *h__;

	qi = qip1;
	qip1 = qip2;
/* L1200: */
    }

    return 0;
} /* legetayl_ */

