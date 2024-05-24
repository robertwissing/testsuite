#ifdef GASOLINE
#ifndef NOCOOLING


/*
 * Integrate "stiff" equations similar to chemical equilibrium
 * equations.
 *
 * The source code below is taken from Mott, D.R. & Oran, E.S., 2001,
 * "CHEMEQ2: A Solver for the Stiff Ordinary Differential Equations of
 * Chemical Kinetics", Naval Research Laboratory,
 * NRL/MR/6400-01-8553.  The report documentation page states under
 * distribution/availability statement states:
 * "Approved for public release; distribution is unlimited."
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "stiff_cooling_cosmo.h"

inline double max(double x, double y) 
{
    if(x > y) return x;
    else return y;
    }

inline double min(double a, double b)
{
    if(a < b) return a;
    else return b;
    }

/* implement fortran sign function: return a with the sign of b */
inline double sign(double a, double b) 
{
    double aabs = fabs(a);
    if(b >= 0.0) return aabs;
    else return -aabs;
    }

/*
 * Set integration parameters and allocate scratch arrays.
 */
STIFF *StiffInit( double eps, int nv, void *Data,
		  void (*derivs)(void *Data, double, const double *, double *,
				 double *)
		  ) 
{
  STIFF *s;
  int i;

  s = (STIFF *) malloc(sizeof(STIFF));
  assert(s!=NULL);

  s->nv = nv;
  s->epsmin = eps;
  s->sqreps = 5.0*sqrt(eps);
  s->epscl = 1.0/eps;
  s->epsmax = 10.0;
  s->dtmin = 1e-15;
  s->itermax = 3; /*Increased from 1 to 3 to speed integration by calculating more correctors.  Adjustable parameter*/
  s->ymin = malloc(nv*sizeof(*(s->ymin)));
  for(i = 0; i < nv; i++)
      s->ymin[i] = 1e-300;
  s->y0 = malloc(nv*sizeof(*(s->y0)));
  s->y1 = malloc(nv*sizeof(*(s->y1)));
  s->q = malloc(nv*sizeof(*(s->q)));
  s->d = malloc(nv*sizeof(*(s->d)));
  s->rtau = malloc(nv*sizeof(*(s->rtau)));
  s->ys = malloc(nv*sizeof(*(s->ys)));
  s->qs = malloc(nv*sizeof(*(s->qs)));
  s->rtaus = malloc(nv*sizeof(*(s->rtaus)));
  s->scrarray = malloc(nv*sizeof(*(s->scrarray)));

  s->Data = Data;
  s->derivs = derivs;

  s->epsmax = 10.0;

  return s;
}

/*
 * Specify minimum values of quantities.
 */
void StiffSetYMin(STIFF *s, const double *ymin) 
{
    int i;
    
    for(i = 0; i < s->nv; i++)
	s->ymin[i] = ymin[i];
    }

void StiffFinalize( STIFF *s ) 
{
    free(s->ymin);
    free(s->y0);
    free(s->y1);
    free(s->q);
    free(s->d);
    free(s->rtau);
    free(s->ys);
    free(s->qs);
    free(s->rtaus);
    free(s->scrarray);
    free(s);
}

void StiffStep(STIFF *s,
	       double y[],	/* dependent variables */
	       double tstart, 	/* start time */
	       double dtg) 	/* time step */
{
/*
C
cd* * *** ** ********** * * * *** **** * ** ** * *
cd
cd chemeq2(dtg, gsub, n, y)
cd
cd original chemeq development:
cd originators: t.r. young nrl 1982
cd vax version: t.r. young nrl code 4040 may 1983
cd workstation: g. patnaik berkeley research jun 1995
cd
cd chemeq2 development: d.r. mott nrl code 6404 may 1999

   Conversion to C: T. Quinn, UW, Dec, 2011
cd
cd Description: Subroutine chemeq2 solves a class of "stiff" ODEs
cd associated with reactive flow problems that cannot be readily
cd solved by the standard classical methods. In contrast to the
cd original chemeq subroutine, this version uses the same
cd quasi-steady-state update for every species regardless of the
cd timescale for that species. An adaptive stepsize is chosen to
cd give accurate results for the fastest changing quantity, and a
cd stability check on the timestep is also available when the
cd corrector is iterated.
cd
cd NOTE: The accuracy-based timestep calculation can be augmented
cd with a stability-based check when at least three corrector
cd iterations are performed. To include this check, "uncomment"
cd the lines that start with "D", or use the compiler flag -d-lines"
cd if available to compile the code including these lines.  If the
cd lines are manually uncommented, the continuation characters
cd must be placed in the correct column. For most problems, the
cd stability check is not needed, and eliminating the calculations
cd and logic associated with the check enhances performance.
cd 
cd The routine assumes that all of the integrated quantites and the
cd time step are positive.
cd
cd argument list definition (name, type, description, input vs. output):
cd dtg  real    the interval of integration or the i
cd              range of the independent variable.
cd              0.0 <= t <= dtg. (global timestep)
cd gsub real    the name of the derivitive function i
cd              evaluator subroutine.
cd n    integer the number of equations to be i
cd              integrated. an error exisis if n is
cd              greater than nd set by the parameter
cd              statement.
cd y(n) real    the initial values at call time i/o
cd              and the final values at return time.
cd
cd Language and limitations: This subroutine is written in standard
cd FORTRAN 77. For high accuracy, this routine should be compiled
cd using whatever "double precision" flag is appropriate for the
cd platform being used (such as "f77 -r8 . . . .)
cd
cd subroutines referenced:
cd
cd gsub;    whose actual name and definition are supplied by the user
cd 	    is called to obtain the derivitive functions.
cd
cd call gsub(y, q, d, t)
cd argument list to gsub;
cd y(n) real current values of the dependent  i
cd           variable.
cd q(n) real calculated formation rates.      o
cd d(n) real calculated loss rates.           o
cd t    real current value of the independent i 
cd           variable. 
cd
*/

    double tn;			/* time within step */
    int i;
    /*
     * Local copies of Stiff context
     */
    int n = s->nv;
    double *y0 = s->y0;
    double *ymin = s->ymin;
    double *q = s->q;
    double *d = s->d;
    double *rtau = s->rtau;
    double *ys = s->ys;
    double *qs = s->qs;
    double *rtaus = s->rtaus;
    double *scrarray = s->scrarray;
    double *y1 = s->y1;
    double epsmin = s->epsmin;
    double sqreps = s->sqreps;
    double epscl = s->epscl;
    double epsmax = s->epsmax;
    double dtmin = s->dtmin;
    int itermax = s->itermax;
    int gcount = 0;		/* count calls to derivs */
    int rcount = 0;		/* count restart steps */
    double scrtch;
    double ascr;
    double scr1;
    double scr2;
    double dt;			/* timestep used by the integrator */
    double ts;			/* t at start of the chemical timestep */
    double alpha;		/* solution parameter used in update */
    int iter;			/* counter for corrector iterations */
    double eps;			/* maximum correction term */
    double rtaub;
    double qt;			/* alpha weighted average of q */
    double pb;
    const double tfd = 1.000008; /* fudge for completion of timestep */
    double rteps;		/* estimate of sqrt(eps) */
    double dto;			/* old timestep; to rescale rtaus */
    double temp;                
    
    tn = 0.0;
    for(i = 0; i < n; i++) {
	q[i] = 0.0;
	d[i] = 0.0;
	y0[i] = y[i];
	y[i] = max(y[i], ymin[i]);
	}
    
    s->derivs(s->Data, tn + tstart, y, q, d);
    gcount++;
    
    /*
    C
    c estimate the initial stepsize.
    C
    c strongly increasing functions(q >>> d assumed here) use a step-
    c size estimate proportional to the step needed for the function to
    c reach equilibrium where as functions decreasing or in equilibrium
    c use a stepsize estimate directly proportional to the character-
    c istic stepsize of the function. convergence of the integration
    c scheme is likely since the smallest estimate is chosen for the
    c initial stepsize.
    */
    scrtch  = 1.0e-25;
    for(i = 0; i < n; i++) {
	ascr = fabs(q[i]);
	scr2 = sign(1./y[i],.1*epsmin*ascr - d[i]);
	scr1 = scr2 * d[i];
	temp = -fabs(ascr-d[i])*scr2;
	if (y[i] == ymin[i]) temp = 0; /*If the species is already at the minimum, disregard destruction when calculating step size */
	scrtch = max(scr1,max(temp,scrtch));
	}
    dt = min(sqreps/scrtch,dtg);
    while(1) {
	/*
	  c the starting values are stored.
	*/
	ts = tn;
	for(i = 0; i < n; i++) {
	    rtau[i] = dt*d[i]/y[i];
	    ys[i] = y[i];
	    qs[i] = q[i];
	    rtaus[i] = rtau[i];
	    }

	/*
	 * find the predictor terms.
	 */
     restart:
	for(i = 0; i < n; i++) {
	    /*
	     * prediction
	     */
	    double rtaui = rtau[i];
	    /*
	    c note that one of two approximations for alpha is chosen:
	    c 1) Pade b for all rtaui (see supporting memo report)
	    c or
	    c 2) Pade a for rtaui<=rswitch,
	    c linear approximation for rtaui > rswitch
	    c (again, see supporting NRL memo report (Mott et al., 2000))
	    c
	    c Option 1): Pade b
	    */
	    alpha = (180.+rtaui*(60.+rtaui*(11.+rtaui)))
		/(360.+ rtaui*(60. + rtaui*(12. + rtaui)));
	    /*
	    c Option 2): Pade a or linear
	    c
	    c if(rtaui.le.rswitch) then
	    c      alpha = (840.+rtaui*(140.+rtaui*(20.+rtaui)))
	    c    &         / (1680. + 40. * rtaui*rtaui)
	    c else
	    c    alpha = 1.-1./rtaui
	    c end if
	    */
	    scrarray[i] = (q[i]-d[i])/(1.0 + alpha*rtaui);
	    }

	iter = 1;
	while(iter <= itermax) {
	    for(i = 0; i < n; i++) {
		/*
		C ym2(i) = ym1(i)
		C ym1(i) = y(i)
		*/
		y[i] = max(ys[i] + dt*scrarray[i], ymin[i]);
		}
	    /*	    if(iter == 1) {  Removed from original algorithm so that previous, rather than first, corrector is compared to.  Results in faster integration. */
		/*
		c the first corrector step advances the time (tentatively) and
		c saves the initial predictor value as y1 for the timestep
		check later.
		*/
		tn = ts + dt;
		for(i = 0; i < n; i++)
		    y1[i] = y[i];
		/*		} Close for "if(iter == 1)" above */
	    /*
	      evaluate the derivitives for the corrector.
	    */
	    s->derivs(s->Data, tn + tstart, y, q, d);
	    gcount++;
	    eps = 1.0e-10;
	    for(i = 0; i < n; i++) {
		rtaub = .5*(rtaus[i]+dt*d[i]/y[i]);
		/*
		c Same options for calculating alpha as in predictor:
		c
		c Option 1): Pade b
		*/
		alpha = (180.+rtaub*(60.+rtaub*(11.+rtaub)))
		    / (360. + rtaub*(60. + rtaub*(12. + rtaub)));
		/*
		c Option 2): Pade a or linear
		c
		c if(rtaub.le.rswitch)
		c then
		c alpha = (840.+rtaub*(140.+rtaub*(20.+rtaub)))
		c & / (1680. + 40.*rtaub*rtaub)
		c else
		c alpha = 1.- 1./rtaub
		c end if
		*/
		qt = qs[i]*(1. - alpha) + q[i]*alpha;
		pb = rtaub/dt;
		scrarray[i] = (qt - ys[i]*pb) / (1.0 + alpha*rtaub);
		}
	    iter++;
	    }
	/*
	c calculate new f, check for convergence, and limit decreasing
	c functions. the order of the operations in this loop is important.
	*/
	for(i = 0; i < n; i++) {
	    scr2 = max(ys[i] + dt*scrarray[i], 0.0);
	    scr1 = fabs(scr2 - y1[i]);
	    y[i] = max(scr2, ymin[i]);
	    /*
	    C ym2(i) = ymi(i)
	    C yml(i) = y(i)
	    */
	    if(.25*(ys[i] + y[i]) > ymin[i]) {
		scr1 = scr1/y[i];
		eps = max(.5*(scr1+
			      min(fabs(q[i]-d[i])/(q[i]+d[i]+1.0e-30),scr1)),eps);
		}
	    }
	eps = eps*epscl;
	/* 
	   print out dianostics if stepsize becomes too small.
	*/
	if(dt <= dtmin + 1.0e-16*tn) {
	    fprintf(stderr, "stiffchem: step size too small\n");
	    assert(0);
	    }
	/*
	c check for convergence.
	c
	c The following section is used for the stability check
	C       stab = 0.01
	C if(itermax.ge.3) then
	C       do i=1,n
	C           stab = max(stab, abs(y(i)-yml(i))/
	C       &       (abs(ymi(i)-ym2(i))+1.e-20*y(i)))
	C end do
	C endif
	*/
	if(eps <= epsmax) {
	    /*
	      & .and.stab.le.1.
	    c
	    c Valid step. Return if dtg has been reached.
	    */
	    if(dtg <= tn*tfd) return;
	    }
	else {
	    /*
	      Invalid step; reset tn to ts
	    */
	    tn = ts;
	    }
	/*
	  perform stepsize modifications.
	  estimate sqrt(eps) by newton iteration.
	*/
	rteps = 0.5*(eps + 1.0);
	rteps = 0.5*(rteps + eps/rteps);
	rteps = 0.5*(rteps + eps/rteps);

	dto = dt;
	dt = min(dt*(1.0/rteps+.005), tfd*(dtg - tn));
	if(dt <= dtmin + 1.0e-16*tn) {
	    fprintf(stderr, "stiffchem: step size too small\n");
	    assert(0);
	    }
	/* & ,dto/(stab+.001) */
	/*
	  begin new step if previous step converged.
	*/
	if(eps > epsmax) {
	    /*    & .or. stab. gt. 1 */
	    rcount++;
	    /*
	    c After an unsuccessful step the initial timescales don't
	    c change, but dt does, requiring rtaus to be scaled by the
	    c ratio of the new and old timesteps.
	    */
	    dto = dt/dto;
	    for(i = 0; i < n; i++) {
		rtaus[i] = rtaus[i]*dto;
		}
	    /*
	     * Unsuccessful steps return to line 101 so that the initial
	     * source terms do not get recalculated.
	    */
	    goto restart;
	    }
	/*
	  Successful step; get the source terms for the next step
	  and continue back at line 100
	*/
	s->derivs(s->Data, tn + tstart, y, q, d);
	gcount++;
	}
    }

/* test case */
void derivs(double x, double y[], double dydx[])
{
  dydx[0] = -0.013*y[0]-1000.0*y[0]*y[2];
  dydx[1] = -2500.0*y[1]*y[2];
  dydx[2] = -0.013*y[0]-1000.0*y[0]*y[2]-2500.0*y[1]*y[2];
}

/*
 * The following is an implementation of the Brent root finding
 * algorithm.  This is based on code from the GSL which has the
 * following copyright.  The code was modified to avoid the overhead
 * that is somewhat unnecessary for such a simple routine.
 */
/* roots/brent.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

const int itmaxRoot = 100;
const double epsRoot = 3.0e-8;

double RootFind(double (*func)(void *Data, double), void *Data, double x1,
		double x2, double tol)
{
    int i;
    double a = x1;
    double fa = (*func)(Data, a);
    double b = x2;
    double fb = (*func)(Data, b);
    double c = x2;
    double fc = fb;
    double d = x2 - x1;
    double e = x2 - x1;
    
    if ((fa < 0.0 && fb < 0.0) || (fa > 0.0 && fb > 0.0))
    {
	fprintf(stderr, "RootFind: endpoints do not straddle y=0");
	assert(0);
	}
    
    for(i = 0; i < itmaxRoot; i++) {
	double m;
	double dTol;

	int ac_equal = 0;

	if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {
	    ac_equal = 1;
	    c = a;
	    fc = fa;
	    d = b - a;
	    e = b - a;
	    }
  
	if (fabs (fc) < fabs (fb)) {
	    ac_equal = 1;
	    a = b;
	    b = c;
	    c = a;
	    fa = fb;
	    fb = fc;
	    fc = fa;
	    }
  
	dTol = 0.5 * epsRoot * fabs (b) + 0.5*tol;
	m = 0.5 * (c - b);
  
	if (fb == 0) {
	    return b;	/* SUCCESS */
	    }
	if (fabs (m) <= dTol) {
	    return b; /* SUCCESS */
	    }
  
	if (fabs (e) < dTol || fabs (fa) <= fabs (fb)) {
	    d = m;            /* use bisection */
	    e = m;
	    }
	else {
	    double p, q, r;   /* use inverse cubic interpolation */
	    double s = fb / fa;
      
	    if (ac_equal) {
		p = 2 * m * s;
		q = 1 - s;
		}
	    else {
		q = fa / fc;
		r = fb / fc;
		p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
		q = (q - 1) * (r - 1) * (s - 1);
		}
	    if (p > 0) {
		q = -q;
		}
	    else {
		p = -p;
		}
      
	    if (2 * p < min(3 * m * q - fabs (dTol * q), fabs (e * q))) {
		e = d;
		d = p / q;
		}
	    else {
		/* interpolation failed, fall back to bisection */
		d = m;
		e = m;
		}
	    }
	a = b;
	fa = fb;
  
	if (fabs (d) > dTol) {
	    b += d;
	    }
	else {
	    b += (m > 0 ? +dTol : -dTol);
	    }
	fb = (*func)(Data, b);
	}
    
    fprintf(stderr, "brent: number of interations exceeded");
    assert(0);
    return 0.0;
    }

#endif
#endif

