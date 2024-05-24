
#ifdef STARFORM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "feedback.h"
#include "supernova.h"

void snCalcWindFeedback(SN sn, SFEvent sfEvent,
                        double dTimeYr, /* current time in years */
                        double dDelta, /* length of timestep (years) */
                        FBEffects *fbEffects);

void snCalcUVFeedback(SN sn, SFEvent sfEvent,
		      double dTimeYr, /* current time in years */
		      double dDelta, /* length of timestep (years) */
		      FBEffects *fbEffects);

/*
 * Feedback module for GASOLINE
 */

double mod(double a, int b) {
    /* printf("MOD: a=%g, b=%d, floor(a/b)=%g\n",a,b,floor(a/b));*/
    return (a-b*floor(a/b));
    }

void fbInitialize(FB *pfb)
{
    FB fb;
    
    fb = (FB) malloc(sizeof(struct fbContext));
    assert(fb != NULL);
    
    fb->dGmUnit = 0.0;
    fb->dSecUnit = 0.0;
    fb->dErgPerGmUnit = 0.0;
    *pfb = fb;
    }

void pkdFeedback(PKD pkd, FB fb, SN sn, double dTime, double dDelta,
		 FBEffects *fbTotals)
    {
    int i;
    PARTICLE *p;
    int n = pkdLocal(pkd);
    SFEvent sfEvent;
    FBEffects fbEffects;
    double dTotMassLoss;
    double dTotMetals;
    double dTotMOxygen;
    double dTotMIron;
    double dSNIaMassStore;
    double dNSNII, dProb, dStarAge, dMinAge;
    double dRandomNum;
    double dTimeYr = dTime*fb->dSecUnit/SEC_YR;
    double dDeltaYr = dDelta*fb->dSecUnit/SEC_YR;
    int j;
        
    for(i = 0; i < FB_NFEEDBACKS; i++) {
        fbTotals[i].dMassLoss = 0.0;
        fbTotals[i].dEnergy = 0.0;
        fbTotals[i].dMetals = 0.0;
        fbTotals[i].dMIron = 0.0;
        fbTotals[i].dMOxygen = 0.0;
	}
    
    for(i = 0; i < n; ++i) {
        p = &pkd->pStore[i];
#ifdef FBPARTICLE
        if (TYPETest(p, TYPE_FEEDBACK)) printf("FBP status: %d: %g %g %g %g\n",p->iOrder,dTimeYr,p->fTimeCoolIsOffUntil*fb->dSecUnit/SEC_YR,p->u,p->uDotFB);
#endif
        if (pkdIsStar(pkd, p) && p->fTimeForm >= 0.0) {
            dTotMassLoss = 0.0;
            dTotMetals = 0.0;
            dTotMOxygen = 0.0;
            dTotMIron = 0.0;
            p->uDotFB = 0.0;
            p->fNSN = 0.0;
#ifdef FBPARTICLE
            {
            double tstart = 4e6, tend = 30e6; /* yr */
            double mdotonmstar = 25/100./(tend-tstart); /* gm / gm / yr */
            double edotonmstar = 1e51/(100*2e33)/(tend-tstart); /* erg / gm / yr */
            double mFB = fb->dFBMassRatio*p->fMassForm; /* mass of fb particles (code units) */
            double tFB0,tFB1,nFac,dDeltaFB;
            int nFB0,nFB1;
            
            tFB1 =  dTimeYr-p->fTimeForm*fb->dSecUnit/SEC_YR; /* yr */
            tFB0 =  tFB1 - dDeltaYr;

            if (tFB1 > tstart && tFB0 < tend) {
                nFac = mdotonmstar*p->fMassForm/mFB; /* yr^-1 */
                nFB0 = floor(nFac*(tFB0 > tstart ? tFB0-tstart : 0)+0.9);
                nFB1 = floor(nFac*(tFB1 < tend ? tFB1-tstart : tend-tstart)+0.9);
                printf("FBP: %d %g %g %g %d %d\n",p->iOrder,dTimeYr,tFB0,tFB1,nFB0,nFB1);
                while (nFB1 > nFB0) {
                    PARTICLE pNew;
                    pNew = *p;
                    TYPEReset(&pNew, TYPE_STAR|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE|TYPE_ACTIVE);
                    TYPESet(&pNew, TYPE_GAS|TYPE_FEEDBACK);
                    
                    /* For mass conservation would need to do this ... */
                    /*p->fMass -= mFB;
                      assert(p->fMass > 0);*/
                    dDeltaFB = 1/(nFac*fb->dSecUnit/SEC_YR); /* code time */
                    if (dDeltaFB < dDelta) dDeltaFB = dDelta;

                    pNew.fMass = mFB;
                    pNew.u = (edotonmstar/mdotonmstar)/fb->dErgPerGmUnit;
                    pNew.uPred = pNew.u;
                    pNew.uDotFB = (pNew.u*0.999)/dDeltaFB;
                    pNew.fTimeCoolIsOffUntil = dTime + 0.9999*dDeltaFB;
                    pNew.c = sqrt((5./3.*2/3.)*pNew.u); /* Estimate final C */
                    {
                    double ph,dt_diff, dt;
                    struct GasPressureContext gpc;

                    ph = sqrt(p->fBall2*0.25); /* gas h at star */
                    dt = fb->dtFacCourant*ph/(2*pNew.c);	/* dt Courant */
#ifdef THERMALCOND
                        { 
                    gpc.dThermalCondCoeffCode = fb->dThermalCondCoeffCode*fb->a;
                    gpc.dThermalCondSatCoeff = fb->dThermalCondSatCoeff/fb->a;
                    gpc.dThermalCond2CoeffCode = fb->dThermalCond2CoeffCode*fb->a;
                    gpc.dThermalCond2SatCoeff = fb->dThermalCond2SatCoeff/fb->a;
                        }
                    pNew.fDensity = p->fDensity;
                    pkdSetThermalCond(pkd, &gpc, &pNew);
                    if (p->fThermalCond > 0) {
                        dt_diff = fb->dtFacDiffusion*ph*ph*pNew.fDensity/pNew.fThermalCond; /* gas density at star */
                        if (dt_diff < dt) dt = dt_diff;   
                        }
#endif
#if (1)
                    pNew.dt = dt; 
                    pNew.iRung = pkdOneParticleDtToRung( 0,dDelta,pNew.dt );
#else
                    pNew.dt = 1.5625e-9; /* horrible hack! -- should be gas min dt maybe? */
                    pNew.iRung = 6; /* horrible hack! */
#endif
                    }
                    pNew.u *= 0.001;
                    pNew.uPred = pNew.u;

                    pNew.PoverRho2 = (2/3.)*pNew.uPred/p->fDensity;
                    pNew.uDot = 0;
                    pNew.uDotPdV = 0;
                    pNew.uDotAV = 0;
                    pNew.uDotDiff = 0;
                    pkdNewParticle(pkd, pNew);
                    printf("FBP: Particle made %g %g %g %g %d\n",mFB,pNew.u,dDeltaFB*fb->dSecUnit/SEC_YR,pNew.dt,pNew.iRung);
                    nFB1--;
                    }
                }
            }
        continue; /* skip normal feedback */
#endif                

	    sfEvent.dMass = p->fMassForm*fb->dGmUnit/MSOLG;
	    sfEvent.dTimeForm = p->fTimeForm*fb->dSecUnit/SEC_YR;
	    dStarAge = dTime - sfEvent.dTimeForm;
	    sfEvent.dMetals = p->fMetals;
	    sfEvent.dMFracOxygen = p->fMFracOxygen;
	    sfEvent.dMFracIron = p->fMFracIron;

	    /*
	     * Call all the effects in order and accumulate them.
	     */
	    dSNIaMassStore=0.0;  /* Stores mass loss of Ia so as
				    not to double count it in 
				    wind feedback */
                                    
	    for(j = 0; j < FB_NFEEDBACKS; j++) {
		dNSNII = 0;
		switch (j) {
		case FB_SNII:
		    snCalcSNIIFeedback(sn, sfEvent, dTimeYr,
				       dDeltaYr, &fbEffects);
		    if( sn->dESN > 0.0)
			dNSNII = fbEffects.dEnergy * MSOLG*fbEffects.dMassLoss/
			    sn->dESN;
		    /* Blow winds before SN (power of winds ~ power of SN) */
		    dMinAge = dSTLtimeMStar(&sn->ppdva, sn->dMSNIImax, 
					    sfEvent.dMetals); 
		    if (dNSNII > 0 && sn->iNSNIIQuantum > 0 && 
			dStarAge > dMinAge) {
			/* Make sure only a iNSNIIQuantum number of
			 * SNII go off at a time */
			dProb = mod(dNSNII, sn->iNSNIIQuantum)/
			    sn->iNSNIIQuantum;
			dRandomNum = (rand()/((double) RAND_MAX));
			/*	  printf("Random Number = %g\n",dRandomNum);*/
			if(dRandomNum < dProb) /* SN occurred */ {
			    /* Adds missing part to make up quantum */
			    p->fNSN = dNSNII + (1.-dProb)*sn->iNSNIIQuantum;
/*			    printf("NSN: +factor=%g   dNSNII=%g  result=%g fNSN=%g\n",(1.-dProb)*sn->iNSNIIQuantum,dNSNII,dNSNII + (1.-dProb)*sn->iNSNIIQuantum,p->fNSN);*/
			    } 
			else {
			    p->fNSN = dNSNII - dProb*sn->iNSNIIQuantum;
/*			    printf("NSN: -factor=%g   dNSNII=%g  result=%g fNSN=%g\n",dProb*sn->iNSNIIQuantum,dNSNII,dNSNII - dProb*sn->iNSNIIQuantum,p->fNSN);*/
			    }
			if(p->fNSN < sn->iNSNIIQuantum) p->fNSN = 0;
			fbEffects.dEnergy = p->fNSN*sn->dESN/(MSOLG*fbEffects.dMassLoss);   
			} 
		    else if(dStarAge < dMinAge && sn->iNSNIIQuantum) 
			p->fNSN = dNSNII;
		    else p->fNSN += dNSNII;
		    break;
		case FB_SNIA:
		    snCalcSNIaFeedback(sn, sfEvent, dTimeYr,
				       dDeltaYr, &fbEffects);
		    dSNIaMassStore=fbEffects.dMassLoss;
		    break;
		case FB_WIND:
		    snCalcWindFeedback(sn, sfEvent, dTimeYr,
				       dDeltaYr, &fbEffects);
		    if(dSNIaMassStore < fbEffects.dMassLoss)
			fbEffects.dMassLoss -= dSNIaMassStore;
		    break;
		case FB_UV:
		    snCalcUVFeedback(sn, sfEvent, dTimeYr, dDeltaYr,
				     &fbEffects);
		    break;
		default:
		    assert(0);
		    }

		fbEffects.dMassLoss *= MSOLG/fb->dGmUnit;
		fbEffects.dEnergy /= fb->dErgPerGmUnit;
#ifdef FBPARTICLEMIMIC
		double tstart = 4e6, tend = 30e6; 
		double mdotonmstar = 25/100./(tend-tstart); 
		double edotonmstar = 1e51/(100*2e33)/(tend-tstart); 
		double tFB0,tFB1,dDeltaFB;
		
		tFB1 =  dTimeYr-p->fTimeForm*fb->dSecUnit/SEC_YR; 
		tFB0 =  tFB1 - dDeltaYr;

		if (tFB1 > tstart && tFB0 < tend) {
			fbEffects.dMassLoss = dDeltaYr*p->fMassForm*mdotonmstar/FB_NFEEDBACKS;
			fbEffects.dEnergy = (edotonmstar/mdotonmstar)/fb->dErgPerGmUnit;
		}
		else {
			fbEffects.dMassLoss = 0;
			fbEffects.dEnergy = 0;
		}
#endif

		dTotMassLoss += fbEffects.dMassLoss;
		p->uDotFB += fbEffects.dEnergy*fbEffects.dMassLoss;
		dTotMetals += fbEffects.dMetals*fbEffects.dMassLoss;
		dTotMOxygen += fbEffects.dMOxygen*fbEffects.dMassLoss;
		dTotMIron += fbEffects.dMIron*fbEffects.dMassLoss;

		fbTotals[j].dMassLoss += fbEffects.dMassLoss;
		fbTotals[j].dEnergy += fbEffects.dEnergy*fbEffects.dMassLoss;
		fbTotals[j].dMetals += fbEffects.dMetals*fbEffects.dMassLoss;
		fbTotals[j].dMIron += fbEffects.dMIron*fbEffects.dMassLoss;
		fbTotals[j].dMOxygen += fbEffects.dMOxygen*fbEffects.dMassLoss;
		}

	    /*
	     * Modify star particle
	     */
	    /*        fprintf(stderr,"Mass dTotMassLoss %d %g %g  %g %g\n",p->iOrder,p->fMass,dTotMassLoss,dTime/(fb->dSecUnit/SEC_YR),p->fTimeForm);*/
	    assert(p->fMass > dTotMassLoss);

	    p->fMass -= dTotMassLoss;
	    p->fMSN = dTotMassLoss;
	    /* The SNMetals and uDotFB (SN rate) used to be specific
	       quantities, but we are making them totals as
	       they leave the stars so that they are easier
	       to divvy up among the gas particles in 
	       distSNEnergy in smoothfcn.c.  These quantities
	       will be converted back to specific quantities when
	       they are parts of gas particles. */
	    p->fSNMetals = dTotMetals;
	    p->fMIronOut = dTotMIron;
	    p->fMOxygenOut = dTotMOxygen;
	    p->uDotFB /= dDelta; /* convert to rate */

#ifdef  RADIATIVEBOX 
/* Calculates LW radiation from a stellar particle of a given age and mass (assumes Kroupa IMF), CC */
	    double dAge1, dAge2, dLW1, dLW2;
	    if (dStarAge != 0) {
	      dAge1 = dStarAge;
	    }
	    else {
	      /*Avoids log of zero error*/
	      dAge1 = dDeltaYr;
	    }
	    dAge2 = dStarAge + dDeltaYr;
	    dLW1 = CoolLymanWerner(pkd->Cool, dAge1);
	    dLW2 = CoolLymanWerner(pkd->Cool, dAge2);
  
	    p->CoolParticle.dLymanWerner = (dLW1 + dLW2)/2;
#endif
	    }
	else if(pkdIsGas(pkd, p)){
	    assert(p->u >= 0.0);
	    assert(p->uPred >= 0.0);
#ifdef FBPARTICLE
        if (dTime > p->fTimeCoolIsOffUntil) p->uDotFB = 0;
#else
	    p->uDotFB = 0.0;	/* reset SN heating rate of gas to zero */
	    p->uDotESF = 0.0;	/* reset SN heating rate of gas to zero */
#endif
	    }
	}

    }


void snCalcWindFeedback(SN sn, SFEvent sfEvent,
                        double dTimeYr, /* current time in years */
                        double dDelta, /* length of timestep (years) */
                        FBEffects *fbEffects)
{
    double dMStarMin, dMStarMax;
    double dStarLtimeMin, dStarLtimeMax;
    double dMCumMin, dMCumMax,dMTot;
    double dMmin, dMmax;
    double dMassFracReturned;
    double dMDying;

    /* First determine if dying stars are between 1-8 Msolar

    * stellar lifetimes corresponding to beginning and end of 
    * current timestep with respect to starbirth time in yrs */
    dMmin=1.0;
    dMmax=8.0;
    dStarLtimeMin = dTimeYr - sfEvent.dTimeForm; 
    dStarLtimeMax = dStarLtimeMin + dDelta;

    dMStarMin = dSTMStarLtime(&sn->ppdva, dStarLtimeMax, sfEvent.dMetals); 
    dMStarMax = dSTMStarLtime(&sn->ppdva, dStarLtimeMin, sfEvent.dMetals); 
    assert(dMStarMax >= dMStarMin);

    if (((dMStarMin < dMmax) && (dMStarMax > dMmin)) && dMStarMax > dMStarMin) {

	/* Mass Fraction returned to ISM taken from Weidemann, 1987, A&A 188 74 
	   then fit to function: MFreturned = 0.86 - exp(-Mass/1.1) */
	dMassFracReturned=0.86-exp(-((dMStarMax+dMStarMin)/2.)/1.1);

	dMCumMin = dMSCumMass(&sn->MSparam, dMStarMin);
	dMCumMax = dMSCumMass(&sn->MSparam, dMStarMax);
	dMTot = dMSCumMass(&sn->MSparam,0.0);
	/* Find out mass fraction of dying stars, then multiply by the original
	   mass of the star particle */
	if (dMTot == 0.0){
	    dMDying = 0.0;
	    } else { 
		dMDying = (dMCumMin - dMCumMax)/dMTot;
		}
	dMDying *= sfEvent.dMass;

	/* Figure out feedback effects */
	fbEffects->dMassLoss = dMDying * dMassFracReturned;
	fbEffects->dEnergy = 0.0;    
	/* Use star's metallicity for gas returned */
	fbEffects->dMetals = sfEvent.dMetals; 
	fbEffects->dMIron = sfEvent.dMFracIron; 
	fbEffects->dMOxygen = sfEvent.dMFracOxygen; 

	} else {
	    fbEffects->dMassLoss = 0.0;
	    fbEffects->dEnergy = 0.0;    
	    fbEffects->dMetals = 0.0; 
	    fbEffects->dMIron = 0.0;
	    fbEffects->dMOxygen = 0.0;
	    }
    }
void snCalcUVFeedback(SN sn, SFEvent sfEvent,
		      double dTimeYr, /* current time in years */
		      double dDelta, /* length of timestep (years) */
		      FBEffects *fbEffects)
{
    fbEffects->dMassLoss = 0.0;
    fbEffects->dEnergy = 0.0;
    fbEffects->dMetals = 0.0;
    fbEffects->dMIron = 0.0;
    fbEffects->dMOxygen = 0.0;
    }

#endif
