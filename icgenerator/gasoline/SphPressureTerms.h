/*
   Macros for function:
       void SphPressureTermsSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
   in smoothfcn.c

   This code is used 3 times.  Constrained to use Macros (not functions) because
   Macro definitions change between uses as follows:
   1) with p and q active
     hash define PACTIVE(xxx) xxx
     hash define QACTIVE(xxx) xxx
   2) with just p active
     hash define PACTIVE(xxx) xxx
     hash define QACTIVE(xxx)
   3) with just q active
     hash define PACTIVE(xxx)
     hash define QACTIVE(xxx) xxx

   All Macros and Variables not defined here are defined in smoothfcn.c
 */
#ifndef MHDCLEANONLY
#ifdef DRHODT
    #define DRHODTACTIVE(xxx) xxx
    #ifdef GDFORCE
        #define RHO_DIVV(a,b) (b)
    #else
        #define RHO_DIVV(a,b) (a)
    #endif
#else
    #define DRHODTACTIVE(xxx)
#endif

#ifdef DIFFUSION

    #ifdef FEEDBACKDIFFLIMIT
        #define DIFFUSIONLimitTest() (diffSum == 0 || smf->dTime < p->fTimeCoolIsOffUntil || smf->dTime < q->fTimeCoolIsOffUntil)
    #else
        #define DIFFUSIONLimitTest() (diffSum == 0)
    #endif


    #ifdef DIFFUSIONHARMONIC
        #define DIFFUSIONBase() double diffSum = (p->diff+q->diff); \
                                double diffBase = (DIFFUSIONLimitTest() ? 0 : 4*p->diff*q->diff/diffSum);
    #else
        #define DIFFUSIONBase() double diffSum = (p->diff+q->diff); \
                                double diffBase = (DIFFUSIONLimitTest() ? 0 : diffSum);
    #endif
    #ifdef MASSDIFF
        #define MASSDIFFFAC(pother_) ((pother_)->fMass)
        #define DIFFUSIONMetalsBase() double diffMetalsBase = 4*smf->dMetalDiffusionCoeff*diffBase   \
             /((p->fDensity+q->fDensity)*(p->fMass+q->fMass));
        #define DIFFUSIONMass() \
            { double diff = diffMetalsBase*(p->fMass - q->fMass); \
              PACTIVE( p->fMassDot += diff*p->fMass*rq ); \
              QACTIVE( q->fMassDot -= diff*q->fMass*rp ); \
            }
        #define DIFFUSIONVelocity() \
            { double diff0 = diffMetalsBase*(p->v[0] - q->v[0]); \
              double diff1 = diffMetalsBase*(p->v[1] - q->v[1]); \
              double diff2 = diffMetalsBase*(p->v[2] - q->v[2]); \
              PACTIVE( ACCEL(p,0) += diff0*rq*MASSDIFFFAC(q) ); \
              QACTIVE( ACCEL(q,0) -= diff0*rp*MASSDIFFFAC(p) ); \
              PACTIVE( ACCEL(p,1) += diff1*rq*MASSDIFFFAC(q) ); \
              QACTIVE( ACCEL(q,1) -= diff1*rp*MASSDIFFFAC(p) ); \
              PACTIVE( ACCEL(p,2) += diff2*rq*MASSDIFFFAC(q) ); \
              QACTIVE( ACCEL(q,2) -= diff2*rp*MASSDIFFFAC(p) ); \
            }
    #else
        #define MASSDIFFFAC(pother_) 1
        #define DIFFUSIONMetalsBase() double diffMetalsBase = 2*smf->dMetalDiffusionCoeff*diffBase \
             /(p->fDensity+q->fDensity);
        #define DIFFUSIONMass()
        #define DIFFUSIONVelocity()
    #endif
#else
    #define DIFFUSIONBase()
    #define DIFFUSIONMetalsBase()
    #define DIFFUSIONMass()
    #define DIFFUSIONVelocity()
#endif


#ifdef DIFFUSION
    #if defined(UNONCOOL) && !defined(TWOPHASE)
        #define DIFFUSIONThermaluHot() \
                { double diffuNc = diffTh*(p->uHotPred-q->uHotPred); \
                PACTIVE( p->uHotDotDiff += diffuNc*rq );        \
                QACTIVE( q->uHotDotDiff -= diffuNc*rp );        \
                }
    #else
        #define DIFFUSIONThermaluHot()
    #endif
    #ifdef DIFFUSIONPRICE
        #define DIFFUSIONThermal(dt_) \
            { double irhobar = 2/(p->fDensity+q->fDensity);     \
             double vsig = sqrt(fabs(qPoverRho2*q->fDensity*q->fDensity - pPoverRho2*p->fDensity*p->fDensity)*irhobar); \
             double diffTh = smf->dThermalDiffusionCoeff*0.5*(ph+sqrt(0.25*BALL2(q)))*irhobar*vsig*drterm2; \
             double diffu = diffTh*(p->uPred-q->uPred);             \
             PACTIVE( p->uDotDiff += diffu*rq );                     \
             QACTIVE( q->uDotDiff-= diffu*rp );                     \
             DIFFUSIONThermaluHot(); }
        #define DIFFUSIONShockCondBase() double dShockCond = 0;
        #define DIFFUSIONShockCond() double dShockCond = 0;
    #else
        #ifndef NODIFFUSIONTHERMAL
            /* Default -- thermal diffusion */
            #ifdef THERMALCOND
                #if (1)
                    /* Harmonic average coeff */
                    #define DIFFUSIONThermalCondBase(dt_) double dThermalCondSum = p->fThermalCond + q->fThermalCond; \
                        double dThermalCond = ( dThermalCondSum <= 0 ? 0 : 4*p->fThermalCond*q->fThermalCond/(dThermalCondSum*p->fDensity*q->fDensity) ); \
                        if (dThermalCond > 0 && (dt_diff = smf->dtFacDiffusion*ph*ph/(dThermalCond*p->fDensity)) < dt_) dt_ = dt_diff;
                #else
                    /* Arithmetic average coeff */
                    #define DIFFUSIONThermalCondBase(dt_) \
                          double dThermalCond = (p->fThermalCond + q->fThermalCond)/(p->fDensity*q->fDensity); \
                          if (dThermalCond > 0 && (dt_diff = smf->dtFacDiffusion*ph*ph/(dThermalCond*p->fDensity)) < dt_) dt_ = dt_diff;

                #endif
            #else
                #define DIFFUSIONThermalCondBase(dt_) double dThermalCond=0;
            #endif
            #ifdef DIFFUSIONSHOCKCOND
                #define DIFFUSIONShockCondBase() double dShockCond = 0;
                #define DIFFUSIONShockCond() double havrg=0.5*(ph+sqrt(0.25*BALL2(q))); \
                        absmu = abs(havrg*dvdotdr*smf->a/(nnList[i].fDist2+0.01*havrg*havrg)); \
                        double constg = 2.0; \
                        dShockCond = SWITCHCOMBINE(p,q)*havrg*constg*(8*absmu)/(p->fDensity+q->fDensity);
            #else
                #define DIFFUSIONShockCondBase() double dShockCond = 0;
                #define DIFFUSIONShockCond() double dShockCond = 0;
            #endif
            #ifdef DTTEST
                #define DIFFUSIONThermal(dt_) \
                    { double diffTh = (2*smf->dThermalDiffusionCoeff*diffBase*drterm2/(p->fDensity+q->fDensity)); \
                      double dt_diff, diffu;                                                  \
                      DIFFUSIONThermalCondBase(dt_);                                    \
                      p->dt_Sph_cond = dt_; \
                      q->dt_Sph_cond = dt_; \
                      if (diffTh > 0 && (dt_diff= smf->dtFacDiffusion*ph*ph/(diffTh*p->fDensity)) < dt_) dt_ = dt_diff; \
                      p->dt_Sph_diff = dt_; \
                      q->dt_Sph_diff = dt_; \
                      diffu = (diffTh+dThermalCond+dShockCond)*(p->uPred-q->uPred);              \
                      PACTIVE( p->uDotDiff += diffu*rq*MASSDIFFFAC(q) );                \
                      QACTIVE( q->uDotDiff -= diffu*rp*MASSDIFFFAC(p) );                \
                      DIFFUSIONThermaluHot(); }
            #else
                #define DIFFUSIONThermal(dt_) \
                    { double diffTh = (2*smf->dThermalDiffusionCoeff*diffBase*drterm2/(p->fDensity+q->fDensity)); \
                      double dt_diff, diffu;				\
                      DIFFUSIONThermalCondBase(dt_);                                    \
                      if ((diffTh+dShockCond+dThermalCond) > 0 && (dt_diff= smf->dtFacDiffusion*ph*ph/((diffTh+dShockCond+dThermalCond)*p->fDensity)) < dt_) dt_ = dt_diff; \
                      diffu = (diffTh+dThermalCond+dShockCond)*(p->uPred-q->uPred);              \
                      PACTIVE( p->uDotDiff += diffu*rq*MASSDIFFFAC(q) );                \
                      QACTIVE( q->uDotDiff -= diffu*rp*MASSDIFFFAC(p) );                \
                      DIFFUSIONThermaluHot(); }
            #endif
        #else
            #define DIFFUSIONThermal(dt_)
        #endif
    #endif

    #define DIFFUSIONMetals() \
        { double diff = diffMetalsBase*(p->fMetals - q->fMetals); \
          PACTIVE( p->fMetalsDot += diff*rq*MASSDIFFFAC(q) ); \
          QACTIVE( q->fMetalsDot -= diff*rp*MASSDIFFFAC(p) ); }
    #ifdef STARFORM
        #define DIFFUSIONMetalsOxygen() \
            { double diff = diffMetalsBase*(p->fMFracOxygen - q->fMFracOxygen); \
              PACTIVE( p->fMFracOxygenDot += diff*rq*MASSDIFFFAC(q) ); \
              QACTIVE( q->fMFracOxygenDot -= diff*rp*MASSDIFFFAC(p) ); }
        #define DIFFUSIONMetalsIron() \
            { double diff = diffMetalsBase*(p->fMFracIron - q->fMFracIron); \
              PACTIVE( p->fMFracIronDot += diff*rq*MASSDIFFFAC(q) ); \
              QACTIVE( q->fMFracIronDot -= diff*rp*MASSDIFFFAC(p) ); }
    #else
        #define DIFFUSIONMetalsOxygen()
        #define DIFFUSIONMetalsIron()
    #endif /* STARFORM */
#else /* No diffusion */
    #define DIFFUSIONShockCondBase()
    #define DIFFUSIONShockCond()
    #define DIFFUSIONThermal(dt_)
    #define DIFFUSIONMetals()
    #define DIFFUSIONMetalsOxygen()
    #define DIFFUSIONMetalsIron()
#endif

#if defined(VARALPHA) || defined(CULLENDEHNEN)
    #define ALPHA (smf->alpha*0.5*(p->alpha+q->alpha))
#ifdef PRICEVISCBETA
    #define BETA  (smf->beta)
#else
    #define BETA  (smf->beta*0.5*(p->alpha+q->alpha))
#endif
#else
    #define ALPHA (smf->alpha)
    #define BETA  (smf->beta)
#endif
#endif

#define SETDTNEW_PQ(dt_)  { if (dt_ < p->dtNew) p->dtNew=dt_; \
                            if (dt_ < q->dtNew) q->dtNew=dt_; \
                            if (4*q->dt < p->dtNew) p->dtNew = 4*q->dt; \
                            if (4*p->dt < q->dtNew) q->dtNew = 4*p->dt; }
#ifndef CONDUCTIONONLY
#ifndef MHDCLEANONLY
#ifdef VSIGVISC
    #define ARTIFICIALVISCOSITY(visc_,dt_) { absmu = -dvdotdr*smf->a            \
                /sqrt(nnList[i].fDist2); /* mu multiply by a to be consistent with physical c */ \
            if (absmu>p->mumax) p->mumax=absmu; /* mu terms for gas time step */ \
    		if (absmu>q->mumax) q->mumax=absmu; \
    		visc_ = (ALPHA*(vSigp + vSigq) + BETA*1.5*absmu); \
    		dt_ = smf->dtFacCourant*ph/(0.625*(vSigp + vSigq)+0.375*visc_);     \
    		visc_ = SWITCHCOMBINE(p,q)*visc_ \
    		    *absmu/(pDensity + q->fDensity); }
#else
    #define ARTIFICIALVISCOSITY(visc_,dt_) { double hav=0.5*(ph+sqrt(0.25*BALL2(q)));  /* h mean */ \
    		absmu = -hav*dvdotdr*smf->a  \
    		    /(nnList[i].fDist2+0.01*hav*hav); /* mu multiply by a to be consistent with physical c */ \
    		if (absmu>p->mumax) p->mumax=absmu; /* mu terms for gas time step */ \
    		if (absmu>q->mumax) q->mumax=absmu; \
    		visc_ = (ALPHA*(vSigp + vSigq) + BETA*2*absmu);	\
    		dt_ = smf->dtFacCourant*hav/(0.625*(vSigp + vSigq)+0.375*visc_); \
    		visc_ = SWITCHCOMBINE(p,q)*visc_ \
    		    *absmu/(pDensity + q->fDensity); }
#endif

    /* Force Calculation between particles p and q */
DRHODTACTIVE( PACTIVE( p->fDivv_PdV -= rq/p->fDivv_Corrector/RHO_DIVV(pDensity,q->fDensity)*dvdotdr; ));
DRHODTACTIVE( QACTIVE( q->fDivv_PdV -= rp/p->fDivv_Corrector/RHO_DIVV(q->fDensity,pDensity)*dvdotdr; ));
DRHODTACTIVE( PACTIVE( p->fDivv_PdVcorr -= rq/RHO_DIVV(pDensity,q->fDensity)*dvdotdr; ));
DRHODTACTIVE( QACTIVE( q->fDivv_PdVcorr -= rp/RHO_DIVV(q->fDensity,pDensity)*dvdotdr; ));
PACTIVE( p->uDotPdV += rq*PRES_PDV(pPoverRho2,qPoverRho2)*dvdotdrI; );
QACTIVE( q->uDotPdV += rp*PRES_PDV(qPoverRho2,pPoverRho2)*dvdotdrI; );
//if (p->iOrder == 865177) printf("sphp %d: %g %g %g %g\n",p->iOrder,p->u,pPoverRho2,dvdotdr,rq*PRES_PDV(pPoverRho2,qPoverRho2)*dvdotdr );
//if (q->iOrder == 865177) printf("sphq %q: %g %g %g %g\n",q->iOrder,q->u,qPoverRho2,dvdotdr,rp*PRES_PDV(qPoverRho2,pPoverRho2)*dvdotdr );
double dvec[3] = {dx,dy,dz};
#ifdef USEVEC
PACTIVE( multvs(Accp,dvec,PRES_ACC(pPoverRho2f,qPoverRho2f)); );
QACTIVE( multvs(Accq,dvec,PRES_ACC(qPoverRho2f,pPoverRho2f)); );
#else
PACTIVE( Accp = (PRES_ACC(pPoverRho2f,qPoverRho2f)); );
QACTIVE( Accq = (PRES_ACC(pPoverRho2f,qPoverRho2f)); );
#endif
//DIFFUSIONShockCondBase();
if (dvdotdr>=0.0) {
    dt = smf->dtFacCourant*ph/(2*(vSigp > vSigq ? vSigp : vSigq));
    #ifdef DTTEST
        p->dt_Sph_dvdotdr = dt;
        q->dt_Sph_dvdotdr = dt;
        p->dt_Sph_av = 0;
        q->dt_Sph_av = 0;
    #endif
    }
else {
    ARTIFICIALVISCOSITY(visc,dt); /* Calculate Artificial viscosity terms and associated dt */
    #ifdef DTTEST
        p->dt_Sph_dvdotdr = 0;
        q->dt_Sph_dvdotdr = 0;
        p->dt_Sph_av = dt;
        q->dt_Sph_av = dt;
    #endif
        PACTIVE( p->uDotAV += rq*(0.5*visc)*dvdotdrI; );
        QACTIVE( q->uDotAV += rp*(0.5*visc)*dvdotdrI; );
        PACTIVE( Accp += visc );
        QACTIVE( Accq += visc );
#ifdef SETPRANDTLVISC
        PACTIVE( p->AccAV[0] -= visc*rq*aFac*dxterm );
        QACTIVE( q->AccAV[0] += visc*rp*aFac*dxterm );
        PACTIVE( p->AccAV[1] -= visc*rq*aFac*dyterm );
        QACTIVE( q->AccAV[1] += visc*rp*aFac*dyterm );
        PACTIVE( p->AccAV[2] -= visc*rq*aFac*dzterm );
        QACTIVE( q->AccAV[2] += visc*rp*aFac*dzterm );
#endif
	}
#endif
#ifdef MHD
#ifndef MHDCLEANONLY
#ifndef GLASSMHDBETA
  /* Induction Eq */
  PACTIVE( p->BDotRho[0] += (p->BPred[0]*dvdotdrI-dvx*Bdotdrp)*rqoverRhoq;);
  QACTIVE( q->BDotRho[0] += (q->BPred[0]*dvdotdrI-dvx*Bdotdrq)*rpoverRhop;);
  PACTIVE( p->BDotRho[1] += (p->BPred[1]*dvdotdrI-dvy*Bdotdrp)*rqoverRhoq;);
  QACTIVE( q->BDotRho[1] += (q->BPred[1]*dvdotdrI-dvy*Bdotdrq)*rpoverRhop;);
  PACTIVE( p->BDotRho[2] += (p->BPred[2]*dvdotdrI-dvz*Bdotdrp)*rqoverRhoq;);
  QACTIVE( q->BDotRho[2] += (q->BPred[2]*dvdotdrI-dvz*Bdotdrq)*rpoverRhop;);
//if (p->iOrder == 0) printf("dvdotdr %f Bdotdrp %f iorder %d \n", dvdotdr,Bdotdrp,q->iOrder);
#endif
  /* Add magnetic pressure term*/
  PACTIVE( Accp += 0.5*(PRES_ACC(pB2*igDensity2p,qB2*igDensity2q)););
  QACTIVE( Accq += 0.5*(PRES_ACC(pB2*igDensity2p,qB2*igDensity2q)););
  /* Add magnetic tensor terms*/
  PACTIVE( Accpx = -(PRES_ACC(p->BPred[0]*Bdotdrp*igDensity2p,q->BPred[0]*Bdotdrq*igDensity2q))*rq*aFac;);
  QACTIVE( Accqx = -(PRES_ACC(p->BPred[0]*Bdotdrp*igDensity2p,q->BPred[0]*Bdotdrq*igDensity2q))*rp*aFac;);
  PACTIVE( Accpy = -(PRES_ACC(p->BPred[1]*Bdotdrp*igDensity2p,q->BPred[1]*Bdotdrq*igDensity2q))*rq*aFac;);
  QACTIVE( Accqy = -(PRES_ACC(p->BPred[1]*Bdotdrp*igDensity2p,q->BPred[1]*Bdotdrq*igDensity2q))*rp*aFac;);
  PACTIVE( Accpz = -(PRES_ACC(p->BPred[2]*Bdotdrp*igDensity2p,q->BPred[2]*Bdotdrq*igDensity2q))*rq*aFac;);
  QACTIVE( Accqz = -(PRES_ACC(p->BPred[2]*Bdotdrp*igDensity2p,q->BPred[2]*Bdotdrq*igDensity2q))*rp*aFac;);
  /* Add tensile correction */
  PACTIVE( Accpx += Bfacp*p->BPred[0]*(PRES_ACC(Bdotdrp*igDensity2p,Bdotdrq*igDensity2q))*rq*aFac;);
  QACTIVE( Accqx += Bfacq*q->BPred[0]*(PRES_ACC(Bdotdrp*igDensity2p,Bdotdrq*igDensity2q))*rp*aFac;);
  PACTIVE( Accpy += Bfacp*p->BPred[1]*(PRES_ACC(Bdotdrp*igDensity2p,Bdotdrq*igDensity2q))*rq*aFac;);
  QACTIVE( Accqy += Bfacq*q->BPred[1]*(PRES_ACC(Bdotdrp*igDensity2p,Bdotdrq*igDensity2q))*rp*aFac;);
  PACTIVE( Accpz += Bfacp*p->BPred[2]*(PRES_ACC(Bdotdrp*igDensity2p,Bdotdrq*igDensity2q))*rq*aFac;);
  QACTIVE( Accqz += Bfacq*q->BPred[2]*(PRES_ACC(Bdotdrp*igDensity2p,Bdotdrq*igDensity2q))*rp*aFac;);
  #ifdef MHDRESISTIVITY
  /*Artificial resistivity terms*/
  PACTIVE( p->uDotBdiss -= 0.25*dB2*(PRES_ACC(alphaBp*vSigBp*igDensity2p,alphaBq*vSigBq*igDensity2q))*rq*drterm;);
  PACTIVE( p->BDotDiss[0] += 0.5*dBx*p->fDensity*(PRES_ACC(alphaBp*igDensity2p*vSigBp,alphaBq*igDensity2q*vSigBq))*rq*drterm;);
  PACTIVE( p->BDotDiss[1] += 0.5*dBy*p->fDensity*(PRES_ACC(alphaBp*igDensity2p*vSigBp,alphaBq*igDensity2q*vSigBq))*rq*drterm;);
  PACTIVE( p->BDotDiss[2] += 0.5*dBz*p->fDensity*(PRES_ACC(alphaBp*igDensity2p*vSigBp,alphaBq*igDensity2q*vSigBq))*rq*drterm;);
  QACTIVE( q->uDotBdiss -= 0.25*dB2*(PRES_ACC(alphaBp*vSigBp*igDensity2p,alphaBq*vSigBq*igDensity2q))*rp*drterm;);
  QACTIVE( q->BDotDiss[0] -= 0.5*dBx*q->fDensity*(PRES_ACC(alphaBp*igDensity2p*vSigBp,alphaBq*igDensity2q*vSigBq))*rp*drterm;);
  QACTIVE( q->BDotDiss[1] -= 0.5*dBy*q->fDensity*(PRES_ACC(alphaBp*igDensity2p*vSigBp,alphaBq*igDensity2q*vSigBq))*rp*drterm;);
  QACTIVE( q->BDotDiss[2] -= 0.5*dBz*q->fDensity*(PRES_ACC(alphaBp*igDensity2p*vSigBp,alphaBq*igDensity2q*vSigBq))*rp*drterm;);
  #endif
#endif
  PACTIVE( p->divB -= dBdotdr*rqoverRhoq;);
  QACTIVE( q->divB -= dBdotdr*rpoverRhop;);
  PACTIVE(p->curlB[0] -= (dBz*dyterm - dBy*dzterm)*rqoverRhoq;);
  PACTIVE(p->curlB[1] -= (dBx*dzterm - dBz*dxterm)*rqoverRhoq;);
  PACTIVE(p->curlB[2] -= (dBy*dxterm - dBx*dyterm)*rqoverRhoq;);
  QACTIVE(q->curlB[0] -= (dBz*dyterm - dBy*dzterm)*rpoverRhop;);
  QACTIVE(q->curlB[1] -= (dBx*dzterm - dBz*dxterm)*rpoverRhop;);
  QACTIVE(q->curlB[2] -= (dBy*dxterm - dBx*dyterm)*rpoverRhop;);

//SETDTNEW_PQ(dt);
  #ifdef MHDCLEANING
  /* Cleaning equations*/
  PACTIVE( p->BDotBClean[0] -= p->fDensity*(PRES_ACC(p->BCleanPred*chp*igDensity2p,q->BCleanPred*chq*igDensity2q))*rq*dxterm;);
  QACTIVE( q->BDotBClean[0] += q->fDensity*(PRES_ACC(p->BCleanPred*chp*igDensity2p,q->BCleanPred*chq*igDensity2q))*rp*dxterm;);
  PACTIVE( p->BDotBClean[1] -= p->fDensity*(PRES_ACC(p->BCleanPred*chp*igDensity2p,q->BCleanPred*chq*igDensity2q))*rq*dyterm;);
  QACTIVE( q->BDotBClean[1] += q->fDensity*(PRES_ACC(p->BCleanPred*chp*igDensity2p,q->BCleanPred*chq*igDensity2q))*rp*dyterm;);
  PACTIVE( p->BDotBClean[2] -= p->fDensity*(PRES_ACC(p->BCleanPred*chp*igDensity2p,q->BCleanPred*chq*igDensity2q))*rq*dzterm;);
  QACTIVE( q->BDotBClean[2] += q->fDensity*(PRES_ACC(p->BCleanPred*chp*igDensity2p,q->BCleanPred*chq*igDensity2q))*rp*dzterm;);
  /*Apply to x,y,z using times dx,dy,dz*/
  PACTIVE( p->BCleanDot += (chp*dBdotdr+p->BCleanPred*dvdotdrI*0.5)*rqoverRhoq;);
  QACTIVE( q->BCleanDot += (chq*dBdotdr+q->BCleanPred*dvdotdrI*0.5)*rpoverRhop;);
dt=smf->dtFacCourant*ph/(2*chp);
SETDTNEW_PQ(dt);
#endif //MHDCLEANING
#endif
#endif //CONDUCTIONONLY
#ifndef MHDCLEANONLY
#ifdef THERMALCONDANI
double TCIso = 0.0;
double TCAni = 1.0;
double TCondp = 0.0;
double TCondq = 0.0;
TCondp += TCIso*gradTdotdrp*igDensity2p;
TCondq += TCIso*gradTdotdrq*igDensity2q;
//if(p->iOrder==1974) printf ("Tisocond %f",TCondp);
#ifdef MHD
TCondp += TCAni*BgradTdotBdotdrp*igDensity2p;
TCondq += TCAni*BgradTdotBdotdrq*igDensity2q;
#endif
//if(p->iOrder==1974) printf ("TAnicond %f",TCondp);
double dThermalAniCond = TCondp+TCondq;
PACTIVE( p->uDotAniCond += (TCondp+TCondq)*rq*MASSDIFFFAC(q););
QACTIVE( q->uDotAniCond -= (TCondp+TCondq)*rp*MASSDIFFFAC(q););
//dt=(TCondp+TCondq) > 0.0 ? 25000*ph*ph/(TCondp+TCondq) : 100000.0;
//SETDTNEW_PQ(dt);
#endif
#ifndef CONDUCTIONONLY
PACTIVE( Accp *= rq*aFac; );/* aFac - convert to comoving acceleration */
QACTIVE( Accq *= rp*aFac; );
PACTIVE( ACCEL(p,0) -= (Accp * dxterm) + Accpx; );
PACTIVE( ACCEL(p,1) -= (Accp * dyterm) + Accpy; );
PACTIVE( ACCEL(p,2) -= (Accp * dzterm) + Accpz; );
QACTIVE( ACCEL(q,0) += (Accq * dxterm) + Accqx; );
QACTIVE( ACCEL(q,1) += (Accq * dyterm) + Accqy; );
QACTIVE( ACCEL(q,2) += (Accq * dzterm) + Accqz; );

DIFFUSIONBase();
DIFFUSIONShockCondBase();
if (dvdotdr < 0.0) {DIFFUSIONShockCond();}
DIFFUSIONThermal(dt);
DIFFUSIONMetalsBase();
DIFFUSIONMetals();
DIFFUSIONMetalsOxygen();
DIFFUSIONMetalsIron();
DIFFUSIONMass();
DIFFUSIONVelocity();
SETDTNEW_PQ(dt);
#else
PACTIVE( ACCEL(p,0) -= 0.0; );
PACTIVE( ACCEL(p,1) -= 0.0; );
PACTIVE( ACCEL(p,2) -= 0.0; );
QACTIVE( ACCEL(q,0) += 0.0; );
QACTIVE( ACCEL(q,1) += 0.0; );
QACTIVE( ACCEL(q,2) += 0.0; );
#endif //ONLYMHDCLEAN
#endif //CONDUCTIONONLY
