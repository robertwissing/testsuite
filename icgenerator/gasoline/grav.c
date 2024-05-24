
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mdl.h"
#include "pkd.h"
#include "grav.h"
#include "meval.h"
#include "qeval.h"

/*
 * Architectures for which the native sqrt() is faster.
 */

#define NATIVE_SQRT (defined(_MIPS_ISA) && (_MIPS_ISA == _MIPS_ISA_MIPS4) \
	     || defined(__i486__) || defined(CCC) || defined(__ia64__) \
	     || defined(__crayx1) || defined(OPTERON))

#if !(NATIVE_SQRT)
void v_sqrt1(int,double *,double *);
#endif

int pkdBucketInteract(PKD pkd,int iBucket,int iOrder)
{
	PARTICLE *p;
	KDN *pkdn;
	ILP *ilp;
	ILCS *ilcs;
	ILCN *ilcn;
	int n,i,j;
	double fPot,ax,ay,az;
	double x,y,z,dx,dy,dz,d2,h,twoh,a,b,c,d;
	double dir2,qirx,qiry,qirz,qir,tr,qir3;
	double gam[6];
	double idt2; /* reciprocal square of symmetric timestep */
	int nFlop;
	int nActive = 0;
#ifdef COMPLETE_LOCAL
	int nMultiFlop[5] = MEVAL_FLOP;
#else
	int nMultiFlop[5] = QEVAL_FLOP;
#endif
#if (NATIVE_SQRT)
	double dir;
#endif
	int nPart, nCellSoft, nCellNewt;
	double *d2a; 
	double *sqrttmp; 

#ifdef GR_DRAG
	const double prefac = -5.163e-21; /* c = 1.006e4 in pkdgrav units */
	FLOAT vx,vy,vz;
	double ct2,e;
#endif

	/*
	 ** Now process the two interaction lists for each active particle.
	 */
	pkdn = &pkd->kdNodes[iBucket];
	p = &pkd->pStore[pkdn->pLower];
	n = pkdn->pUpper - pkdn->pLower + 1;
	ilp = pkd->ilp;
	ilcs = pkd->ilcs;
	ilcn = pkd->ilcn;
	d2a = pkd->d2a;
	sqrttmp = pkd->sqrttmp;
	nPart = pkd->nPart;
	nCellSoft = pkd->nCellSoft;
	nCellNewt = pkd->nCellNewt;
	for (i=0;i<n;++i) {
		if (!TYPEQueryACTIVE(&(p[i]))) continue;
		++nActive;
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		fPot = 0.0;
		x = p[i].r[0];
		y = p[i].r[1];
		z = p[i].r[2];
		h = p[i].fSoft;
		/*
		 ** Scoring for Part (+,*)
		 ** 	Without sqrt = (10,8)
		 **     1/sqrt est.  = (6,11)
		 **     SPLINEM      = (0,3)  for soft = (8,30)
		 **     Total        = (16,22)           (24,49)
		 **     			 = 38	  for soft = 73
		 */
#if !(NATIVE_SQRT)
		for (j=0;j<nPart;++j) {
		        dx = x - ilp[j].x;
			dy = y - ilp[j].y;
			dz = z - ilp[j].z;
			d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (nPart>0) v_sqrt1(nPart,d2a,sqrttmp);
#endif
		for (j=0;j<nPart;++j) {
			dx = x - ilp[j].x;
			dy = y - ilp[j].y;
			dz = z - ilp[j].z;
			twoh = h + ilp[j].h;
#if (NATIVE_SQRT)
			d2 = dx*dx + dy*dy + dz*dz;
			SPLINE(d2,twoh,a,b);
#else
			SPLINEM(sqrttmp[j],d2a[j],twoh,a,b);
#endif			
			idt2 = (p[i].fMass + ilp[j].m)*b;
			if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;
			a *= ilp[j].m;
			b *= ilp[j].m;
			fPot -= a;
			ax -= dx*b;
			ay -= dy*b;
			az -= dz*b;
			}
		/*
		 ** Scoring for CellSoft (+,*)
		 ** 	Without sqrt = (27,29)
		 **     1/sqrt est.  = (6,11)
		 **     SPLINEQ      = (0,9)  for soft = (13,62)
		 **     Total        = (33,49)           (46,102)
		 **     			 = 82	  for soft = 148
		 */
#if !(NATIVE_SQRT)
		for (j=0;j<nCellSoft;++j) {
			dx = x - ilcs[j].x;
			dy = y - ilcs[j].y;
			dz = z - ilcs[j].z;
			d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (nCellSoft>0) v_sqrt1(nCellSoft,d2a,sqrttmp);
#endif
		for (j=0;j<nCellSoft;++j) {
			dx = x - ilcs[j].x;
			dy = y - ilcs[j].y;
			dz = z - ilcs[j].z;
			twoh = h + ilcs[j].h;
#if (NATIVE_SQRT)
			d2 = dx*dx + dy*dy + dz*dz;
			dir = 1.0/sqrt(d2);
			SPLINEQ(dir,d2,twoh,a,b,c,d);
#else
			SPLINEQ(sqrttmp[j],d2a[j],twoh,a,b,c,d);
#endif
			qirx = ilcs[j].xx*dx + ilcs[j].xy*dy + ilcs[j].xz*dz;
			qiry = ilcs[j].xy*dx + ilcs[j].yy*dy + ilcs[j].yz*dz;
			qirz = ilcs[j].xz*dx + ilcs[j].yz*dy + ilcs[j].zz*dz;
			qir = 0.5*(qirx*dx + qiry*dy + qirz*dz);
			tr = 0.5*(ilcs[j].xx + ilcs[j].yy + ilcs[j].zz);
			qir3 = b*ilcs[j].m + d*qir - c*tr;
			fPot -= a*ilcs[j].m + c*qir - b*tr;
			ax -= qir3*dx - c*qirx;
			ay -= qir3*dy - c*qiry;
			az -= qir3*dz - c*qirz;
			idt2 = (p[i].fMass + ilcs[j].m)*b;
			if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;
			}
		/*
		 ** Try a cache check to improve responsiveness?
		 */
		mdlCacheCheck(pkd->mdl);
		/*
		 ** Scoring for CellNewt (+,*)
		 ** 	Without sqrt = (5,13)
		 **     1/sqrt est.  = (6,11)
		 **     Subtotal	 = (11,24)  = 35
		 **     Qeval        (Hex)		= 277 (see qeval.h)
		 **     Total        = (85,227) = 312 Flops/Newt-Interact
		 */
#if !(NATIVE_SQRT)
		for (j=0;j<nCellNewt;++j) {
			dx = x - ilcn[j].x;
			dy = y - ilcn[j].y;
			dz = z - ilcn[j].z;
			d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (nCellNewt>0) v_sqrt1(nCellNewt,d2a,sqrttmp);
#endif
		for (j=0;j<nCellNewt;++j) {
			dx = x - ilcn[j].x;
			dy = y - ilcn[j].y;
			dz = z - ilcn[j].z;
#if (NATIVE_SQRT)
			d2 = dx*dx + dy*dy + dz*dz;
			gam[0] = 1.0/sqrt(d2);
#else
			gam[0] = sqrttmp[j];
#endif
			dir2 = gam[0]*gam[0];
			gam[1] = gam[0]*dir2;
			gam[2] = 3*gam[1]*dir2;
			gam[3] = 5*gam[2]*dir2;
			gam[4] = 7*gam[3]*dir2;
			gam[5] = 9*gam[4]*dir2;
#ifdef COMPLETE_LOCAL
			MEVAL(iOrder,ilcn[j],gam,dx,dy,dz,ax,ay,az,fPot);
#else
			QEVAL(iOrder,ilcn[j],gam,dx,dy,dz,ax,ay,az,fPot);
#endif
			idt2 = (p[i].fMass + ilcn[j].m)*gam[1];
			if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;
			}
		p[i].fPot += fPot;
		p[i].a[0] += ax;
		p[i].a[1] += ay;
		p[i].a[2] += az;
		/*
		 ** Try a cache check to improve responsiveness?
		 */
		mdlCacheCheck(pkd->mdl);
		}
	/*
	 ** Do the intra-bucket interactions.
	 ** Scoring (+,*):
	 ** 	without sqrt = (14,17)
	 **     sqrt est.    = (6,11)
	 **     SPLINE       = (0,3)  for soft = (8,30)
	 **     Total        = (20,31)           (28,58)
	 **                  = 51     for soft = 86
	 ** Multiplied by (n*(n-1)/2)!
	 */
	for (i=0;i<n-1;++i) {
		for (j=i+1;j<n;++j) {
			if (!TYPEQueryACTIVE(&(p[i])) 
				&& !TYPEQueryACTIVE(&(p[j]))) continue;
			dx = p[j].r[0] - p[i].r[0];
			dy = p[j].r[1] - p[i].r[1];
			dz = p[j].r[2] - p[i].r[2];
			d2 = dx*dx + dy*dy + dz*dz;
#ifdef COLLISIONS
			/* determine closest approach */
			if (d2<p[i].mindist2) 
				p[i].mindist2=d2;
#endif
			twoh = p[i].fSoft + p[j].fSoft;
			SPLINE(d2,twoh,a,b);
#ifdef GR_DRAG
			vx = p[j].v[0] - p[i].v[0];
			vy = p[j].v[1] - p[i].v[1];
			vz = p[j].v[2] - p[i].v[2];
			ct2 = (dx*vx + dy*vy + dz*vz)*(dx*vx + dy*vy + dz*vz)/(d2*(vx*vx + vy*vy + vz*vz)); /* cos^2(theta) */
			e = prefac*(p[i].fMass + p[j].fMass)*p[i].fMass*p[j].fMass*(12 - 11*ct2)/(d2*d2);
#endif
			idt2 = (p[i].fMass + p[j].fMass)*b;
			if (TYPEQueryACTIVE(&(p[j]))) {
				p[j].fPot -= a*p[i].fMass;
				p[j].a[0] -= dx*b*p[i].fMass;
				p[j].a[1] -= dy*b*p[i].fMass;
				p[j].a[2] -= dz*b*p[i].fMass;
#ifdef GR_DRAG
				p[j].a[0] += vx*e;
				p[j].a[1] += vy*e;
				p[j].a[2] += vz*e;
#endif
				if (idt2 > p[j].dtGrav) p[j].dtGrav = idt2;
				}
			if (TYPEQueryACTIVE(&(p[i]))) {
				p[i].fPot -= a*p[j].fMass;
				p[i].a[0] += dx*b*p[j].fMass;
				p[i].a[1] += dy*b*p[j].fMass;
				p[i].a[2] += dz*b*p[j].fMass;
#ifdef GR_DRAG
				p[i].a[0] -= vx*e;
				p[i].a[1] -= vy*e;
				p[i].a[2] -= vz*e;
#endif
				if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;
				}
			}
		}
	/*
	 ** Compute the nFlop estimate.
	 */
	nFlop = nActive*((nPart + n)*38 + nCellSoft*82 +
					 nCellNewt*(35 + nMultiFlop[iOrder]));
	return(nFlop);
	}


