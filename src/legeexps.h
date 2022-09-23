/* Copyright (C) 2019-2020 The Simons Foundation, Inc., Zydrunas Gimbutas,
 and Vladimir Rokhlin - All Rights Reserved. */

#ifndef __FASTNONO_LEGEEXPS_H__
#define __FASTNONO_LEGEEXPS_H__

/* legeexps.f -- translated by f2c (version 20190311). */

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


/* Subroutine */ int legeexps_(long int *itype, long int *n, double *x,
	double *u, double *v, double *whts)
{
    /* System generated locals */
    long int u_dim1, u_offset, v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int legerts2_(long int *, long int *, double *,
	    double *), legepols_(double *, long int *, double *);
    static double d__;
    static long int i__, j, itype_rts__;


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






/* Subroutine */ int legerts2_(long int *itype, long int *n, double *ts,
	double *whts)
{
    /* System generated locals */
    long int i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
    double atan(double), cos(double);

    /* Local variables */
    static double d__;
    extern /* Subroutine */ int legetayl2_(double *, double *,
	    double *, double *, long int *, long int *, double *,
	    double *);
    static double h__;
    static long int i__, k;
    static double d2, h2;
    static long int n2;
    static double x0, x1;
    static long int ii, kk;
    static double pi, der, eps, pol, der3, pol3, done, xold, d_den__;
    static long int ifodd;
    static double theta, derold, polold;
    static long int ifstop;
    extern /* Subroutine */ int legepol_(double *, long int *, double *
	    , double *);


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






/* Subroutine */ int legetayl2_(double *pol, double *der, double *
	x, double *h__, long int *n, long int *k, double *sum,
	double *sumder)
{
    /* Initialized data */

    static long int ifcalled = 0;

    /* System generated locals */
    long int i__1;
    double d__1;

    /* Local variables */
    static double d__;
    static long int i__;
    static double d7, q0, q1, q2, dn, qi, ddd, two, qip1, dddd, half,
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
	squares[i__ - 1] = (double) (i__1 * i__1);
	prodinv[i__ - 1] = 1 / prods[i__ - 1];
	rnums[i__ - 1] = (double) i__;
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








/* Subroutine */ int legepol_(double *x, long int *n, double *pol,
	double *der)
{
    /* System generated locals */
    long int i__1;
    double d__1;

    /* Local variables */
    static long int k;
    static double pk, pkm1, pkp1;


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








/* Subroutine */ int legepols_(double *x, long int *n, double *pols)
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    static long int k;
    static double pk, pkm1, pkp1;


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











/* Subroutine */ int matmul_0_(int n__, double *a, double *b,
	double *c__, long int *n)
{
    /* System generated locals */
    long int a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
	    i__3;

    /* Local variables */
    static double d__;
    static long int i__, j, k;


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

/* Subroutine */ int matmul_(double *a, double *b, double *c__,
	long int *n)
{
    return matmul_0_(0, a, b, c__, n);
    }

/* Subroutine */ int matmua_(double *a, double *b, double *c__,
	long int *n)
{
    return matmul_0_(1, a, b, c__, n);
    }


# endif
