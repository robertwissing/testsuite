#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "pkd.h"
#include "outtype.h"
#include "floattype.h"
#ifdef CRAY_XT3
#include "../xdr/types.h"
#include "../xdr/xdr.h"
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

#ifdef COLLISIONS
#include "collision.h"
#endif /* COLLISIONS */

FLOAT VecType(PKD pkd, PARTICLE *p,int iDim,int iType)
    {
#ifdef GASOLINE
    FLOAT vTemp;
#ifdef TWOPHASE
    double frac = p->fMassHot/p->fMass;
#endif
#ifdef COOLING_MOLECULARH
    /*Define the correlation length used for shielding in H2 calculation from the gas shear CC*/
    double correL = 1.0;
#ifdef NEWSHEAR
    /* Calculated same way as for diffusion*/
    if (p->diff != 0) correL = 0.25*p->fBall2*p->c/p->diff;
    /*Minimum correlation length is the smoothing*/
    if (correL > sqrt(0.25*p->fBall2) || p->diff == 0) correL = sqrt(0.25*p->fBall2);
#else /*NEWSHEAR*/
    /* Shear from curl */
    double shear = sqrt(p->curlv[0]*p->curlv[0] + p->curlv[1]*p->curlv[1] + p->curlv[2]*p->curlv[2]);
    if (shear != 0) correL = p->c/shear;
#endif /*NEWSHEAR*/
#endif  /*COOLING_MOLECULARH */
#endif
    //printf("%f",(p->E0x*p->E0x+p->E0y*p->E0y+p->E0z*p->E0z));
	switch (iType) {
	case OUT_IORDER_ARRAY:
	    return((FLOAT) p->iOrder);
#ifdef GASOLINE
    case OUT_GASDENSITY_ARRAY:
	    return(p->fDensity);
	case OUT_DENSITY_ARRAY:
	    return(p->curlv[0]);
#else
	case OUT_DENSITY_ARRAY:
	    return(p->fDensity);
#endif
	case OUT_DENSITYRFC_ARRAY:
	    return(p->fDensity);
	case OUT_COLOR_ARRAY:
#ifdef COLORCODE
	    return(p->fColor);
#endif
	case OUT_POT_ARRAY:
	    return(p->fPot);
	case OUT_AMAG_ARRAY:
	    return(sqrt(p->a[0]*p->a[0] + p->a[1]*p->a[1] + p->a[2]*p->a[2]));
	case OUT_RUNG_ARRAY:
	    return(p->iRung);
	case OUT_SPHDT_ARRAY:
	case OUT_DT_ARRAY:
	    return(p->dt);
	case OUT_SOFT_ARRAY:
	    return(p->fSoft);
#ifdef GASOLINE
#ifdef DENSITYU
	case OUT_DENSITYU_ARRAY:
	    return(p->fDensityU);
#endif
	case OUT_PRES_ARRAY:
	    return(p->fDensity*p->fDensity*p->PoverRho2);
	case OUT_U_ARRAY:
	    return(p->u);
#ifdef MHD
  case OUT_MHDX_ARRAY:
    return(p->B[0]);
  case OUT_MHDY_ARRAY:
    return(p->B[1]);
  case OUT_MHDZ_ARRAY:
    return(p->B[2]);
  case OUT_MHDCURLX_ARRAY:
    return(p->curlB[0]);
  case OUT_MHDCURLY_ARRAY:
    return(p->curlB[1]);
  case OUT_MHDCURLZ_ARRAY:
    return(p->curlB[2]);
  case OUT_MHDDIV_ARRAY:
    return(p->divB);
  case OUT_MHDUDOTBDISS_ARRAY:
    return(p->uDotBdiss);
  case OUT_MHDBDOTRHO_ARRAY:
    return(p->BDotRho[0]*p->BDotRho[0]+p->BDotRho[1]*p->BDotRho[1]+p->BDotRho[2]*p->BDotRho[2]);
  case OUT_MHDBDOTDISS_ARRAY:
    return(p->BDotDiss[0]*p->BDotDiss[0]+p->BDotDiss[1]*p->BDotDiss[1]+p->BDotDiss[2]*p->BDotDiss[2]);
  case OUT_MHDBDOTCLEAN_ARRAY:
    return(p->BDotBClean[0]*p->BDotBClean[0]+p->BDotBClean[1]*p->BDotBClean[1]+p->BDotBClean[2]*p->BDotBClean[2]);
  case OUT_MHDBCLEAN_ARRAY:
    return(p->BClean);
#endif
#ifdef GETVISCPARA
	case OUT_VISCCOEF_ARRAY:
	  return(p->kinvisc);
#ifdef MHD
	case OUT_RESICOEF_ARRAY:
	  return(p->etares);
#endif
#endif
#ifdef CULLENDEHNEN
        case OUT_AX_ARRAY:
          return(p->a[0]);
        case OUT_AY_ARRAY:
          return(p->a[1]);
        case OUT_AZ_ARRAY:
          return(p->a[2]);
	case OUT_GX_ARRAY:
          return(p->grx);
        case OUT_GY_ARRAY:
          return(p->gry);
        case OUT_GZ_ARRAY:
	return(p->grz);
#endif
#ifdef THERMALCONDANI
case OUT_GRADTX_ARRAY:
  return(p->gradTx);
case OUT_GRADTY_ARRAY:
  return(p->gradTy);
case OUT_GRADTZ_ARRAY:
  return(p->gradTz);
case OUT_UDOTANICOND_ARRAY:
  return(p->uDotAniCond);
#endif
#ifdef CHECKQUALITY
case OUT_Q1_ARRAY:
  //return(p->Q1);
  return(p->ICfmax);
case OUT_Q2_ARRAY:
  //return(sqrt((p->Q2x*p->Q2x+p->Q2y*p->Q2y+p->Q2z*p->Q2z)/3));
  return(sqrt(p->ICvel[0]*p->ICvel[0] + p->ICvel[1]*p->ICvel[1] + p->ICvel[2]*p->ICvel[2]));
case OUT_E0_ARRAY:
  return(sqrt((p->E0x*p->E0x+p->E0y*p->E0y+p->E0z*p->E0z)/3));
  // return(p->E0y);
case OUT_Q4_ARRAY:
  return(sqrt((p->Q4xx*p->Q4xx+p->Q4yy*p->Q4yy+p->Q4zz*p->Q4zz)/3));
#endif
	case OUT_TWOPHASE_ARRAY:
#ifdef TWOPHASE
	    return(p->fMassHot);
#else
	    return(0.);
#endif
	case OUT_UNONCOOL_ARRAY:
#ifdef UNONCOOL
	    return(p->uHot);
#else
	    return(0.);
#endif
	case OUT_TEMPINC_ARRAY:
#ifdef UNONCOOL
	    vTemp = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u+p->uHot, p->fDensity, p->fMetals );
#else
#ifndef NOCOOLING
	    vTemp = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u, p->fDensity, p->fMetals );
#else
	    vTemp = pkd->duTFac*p->u;
#endif
#endif
#ifdef TWOPHASE
	    vTemp = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u*(1-frac)+p->uHot*frac, p->fDensity, p->fMetals );
#endif
	    return(vTemp);
	case OUT_TEMP_ARRAY:
#ifndef NOCOOLING
	    vTemp = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u, p->fDensity, p->fMetals );
#else
	    vTemp = pkd->duTFac*p->u;
#endif
	    return(vTemp);
#ifndef NOCOOLING
	case OUT_UDOT_ARRAY:
	    return(p->uDot);
	case OUT_COOL_ARRAY0:
	    return(COOL_ARRAY0(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY1:
	    return(COOL_ARRAY1(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY2:
	    return(COOL_ARRAY2(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY3:
	    return(COOL_ARRAY3(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY4:
	    return(COOL_ARRAY4(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY5:
	    return(COOL_ARRAY5(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY6:
	    return(COOL_ARRAY6(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY7:
	    return(COOL_ARRAY7(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY8:
	    return(COOL_ARRAY8(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY9:
	    return(COOL_ARRAY9(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY10:
	    return(COOL_ARRAY10(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY11:
	    return(COOL_ARRAY11(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY12:
	    return(COOL_ARRAY12(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY13:
	    return(COOL_ARRAY13(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY14:
	    return(COOL_ARRAY14(pkd->Cool, &p->CoolParticle,p->fMetals));
	case OUT_COOL_ARRAY15:
	    return(COOL_ARRAY15(pkd->Cool, &p->CoolParticle,p->fMetals));
#ifdef COOLING_MOLECULARH
	case OUT_CORREL_ARRAY:
        return correL;
#endif
/*Gas shear in terms of mach number, used when calculating column density*/
#ifdef COOLING_METAL_BROKEN
	case OUT_COOL_SHEAR_ARRAY:
	    return(COOL_SHEAR_ARRAY(p->c, p->curlv[0], p->curlv[1], p->curlv[2], p->iOrder));
#endif


#ifdef  RADIATIVEBOX
	case OUT_COOL_LYMANWERNER_ARRAY:
        return(p->CoolParticle.dLymanWerner); /* Lyman Werner Radiation output array*/
#endif /*RADIATIVEBOX*/

#ifdef DENSITYU
#ifdef COOLING_MOLECULARH
	case OUT_COOL_EDOT_ARRAY: /* Cooling array with H2*/
        return( COOL_EDOT( pkd->Cool, &p->CoolParticle, p->u, p->fDensityU, p->fMetals, p->r, correL) );
	case OUT_COOL_COOLING_ARRAY:
        return( COOL_COOLING( pkd->Cool, &p->CoolParticle, p->u, p->fDensityU, p->fMetals, p->r, correL) );
	case OUT_COOL_HEATING_ARRAY:
        return( COOL_HEATING( pkd->Cool, &p->CoolParticle, p->u, p->fDensityU, p->fMetals, p->r, correL) );
#else
	case OUT_COOL_EDOT_ARRAY:
        if(pkdIsGas(pkd, p)) {
            return( COOL_EDOT( pkd->Cool, &p->CoolParticle, p->u, p->fDensityU, p->fMetals, p->r) );
            }
        else {
            return 0;
            }
	case OUT_COOL_COOLING_ARRAY:
        if(pkdIsGas(pkd, p)) {
            return( COOL_COOLING( pkd->Cool, &p->CoolParticle, p->u, p->fDensityU, p->fMetals, p->r) );
            }
        else {
            return 0;
            }
	case OUT_COOL_HEATING_ARRAY:
        if(pkdIsGas(pkd, p)) {
            return( COOL_HEATING( pkd->Cool, &p->CoolParticle, p->u, p->fDensityU, p->fMetals, p->r) );
            }
        else {
            return 0;
            }
#endif
#else
#ifdef COOLING_MOLECULARH
	case OUT_COOL_EDOT_ARRAY:
        return( COOL_EDOT( pkd->Cool, &p->CoolParticle, p->u, p->fDensity, p->fMetals, p->r, correL) );
	case OUT_COOL_COOLING_ARRAY:
	    return( COOL_COOLING( pkd->Cool, &p->CoolParticle, p->u, p->fDensity, p->fMetals,p->r, correL) );
	case OUT_COOL_HEATING_ARRAY:
	    return( COOL_HEATING( pkd->Cool, &p->CoolParticle, p->u, p->fDensity, p->fMetals, p->r, correL) );
#else
	case OUT_COOL_EDOT_ARRAY:
        return( COOL_EDOT( pkd->Cool, &p->CoolParticle, p->u, p->fDensity, p->fMetals, p->r) );
	case OUT_COOL_COOLING_ARRAY:
	    return( COOL_COOLING( pkd->Cool, &p->CoolParticle, p->u, p->fDensity, p->fMetals, p->r) );
	case OUT_COOL_HEATING_ARRAY:
	    return( COOL_HEATING( pkd->Cool, &p->CoolParticle, p->u, p->fDensity, p->fMetals, p->r) );
#endif
#endif
#endif
	case OUT_BALSARASWITCH_ARRAY:
	    return(p->BalsaraSwitch);
#if defined(VARALPHA) || defined(CULLENDEHNEN)
	case OUT_ALPHA_ARRAY:
	    return(p->alpha);
#endif
#ifdef CD_DEBUG
    case OUT_VSIGMAX_ARRAY:
        return(p->vSigMax);
    case OUT_R_CD_ARRAY:
        return(p->R_CD);
    case OUT_SNORM_ARRAY:
        return(p->SNorm);
    case OUT_SFULL_ARRAY:
        return(p->SFull);
    case OUT_DVDSONSFULL_ARRAY:
        return(p->dvdsonSFull);
    case OUT_ALPHALOC_ARRAY:
        return(p->alphaLoc);
    case OUT_ALPHANOISE_ARRAY:
        return(p->alphaNoise);
	case OUT_DIVV_DENS_ARRAY:
	    return(p->divv_dens);
	case OUT_DIVVDOT_ARRAY:
	    return(p->divvDot);
#endif
	case OUT_DIVV_ARRAY:
	    return(p->divv);
#ifdef CD_DEBUG
	case OUT_DVXDX_ARRAY:
        return(p->dvxdx);
	case OUT_DVYDX_ARRAY:
        return(p->dvydx);
	case OUT_DVZDX_ARRAY:
        return(p->dvzdx);
	case OUT_DVXDY_ARRAY:
        return(p->dvxdy);
	case OUT_DVYDY_ARRAY:
        return(p->dvydy);
	case OUT_DVZDY_ARRAY:
        return(p->dvzdy);
	case OUT_DVXDZ_ARRAY:
        return(p->dvxdz);
	case OUT_DVYDZ_ARRAY:
        return(p->dvydz);
	case OUT_DVZDZ_ARRAY:
        return(p->dvzdz);
#endif
#ifdef DODVDS
	case OUT_DVDS_ARRAY:
	    return(p->dvds);
#endif
	case OUT_CSOUND_ARRAY:
	    return(p->c);
	case OUT_MUMAX_ARRAY:
	    return(p->mumax);
	case OUT_DIVONCONH_ARRAY:
	    return(p->divv/(p->c/sqrt(p->fBall2*0.25)));
	case OUT_DIVONCONX_ARRAY:
	    return(p->divv/(p->c/pow(p->fMass/p->fDensity,1./3.)));
    case OUT_UDOTHYDRO_ARRAY:
    case OUT_PDVRFC_ARRAY:
        return(UDOT_HYDRO(p));
	case OUT_METALS_ARRAY:
	    return(p->fMetals);
#ifdef DIFFUSION
	case OUT_METALSDOT_ARRAY:
	    return(p->fMetalsDot);
#endif
	case OUT_UDOTPDV_ARRAY:
	    return(p->uDotPdV);
	case OUT_UDOTAV_ARRAY:
	    return(p->uDotAV);
	case OUT_UDOTDIFF_ARRAY:
	    return(p->uDotDiff);
#ifdef UNONCOOLDEBUG
	case OUT_UHOTDOT_ARRAY:
	    return(p->uHotDot);
	case OUT_UHOTDOTDIFF_ARRAY:
	    return(p->uHotDotDiff);
	case OUT_UHOTDOTPDV_ARRAY:
	    return(p->uHotDotPdV);
	case OUT_UHOTDOTCONV_ARRAY:
	    return(p->uHotDotConv);
	case OUT_UHOTDOTFB_ARRAY:
	    return(p->uHotDotFB);
#endif
#ifdef SHOCKTRACK
	case OUT_SHOCKTRACKER_ARRAY:
	    return(p->ShockTracker);
	case OUT_DIVRHOV_ARRAY:
	    return(p->divrhov);
#endif
#ifdef SFBOUND
	case OUT_SIGMA2_ARRAY:
	    return(p->fSigma2);
#endif
#ifdef STARFORM
	case OUT_IGASORDER_ARRAY:
	    return((FLOAT) p->iGasOrder);
	case OUT_TIMEFORM_ARRAY:
	    return((FLOAT) p->fTimeForm);
	case OUT_MASSFORM_ARRAY:
	    return((FLOAT) p->fMassForm);
	case OUT_COOLTURNONTIME_ARRAY:
	    return((FLOAT) p->fTimeCoolIsOffUntil);
	case OUT_OXYGENMASSFRAC_ARRAY:
	    return((FLOAT) p->fMFracOxygen);
	case OUT_IRONMASSFRAC_ARRAY:
	    return((FLOAT) p->fMFracIron);
#ifdef DIFFUSION
	case OUT_OXYGENMASSFRACDOT_ARRAY:
	    return((FLOAT) p->fMFracOxygenDot);
	case OUT_IRONMASSFRACDOT_ARRAY:
	    return((FLOAT) p->fMFracIronDot);
#endif
	case OUT_UDOTFB_ARRAY:
	    return((FLOAT) p->uDotFB);
#ifdef CHECKSF
	case OUT_TOFF_YR_ARRAY:
	    return((FLOAT) p->tOff );
	case OUT_TCOOL_YR_ARRAY:
	    return((FLOAT) p->tcool );
	case OUT_TDYN_YR_ARRAY:
	    return((FLOAT) p->tdyn );
	case OUT_RATIOSOUNDDYN_ARRAY:
	    return((FLOAT) p->ratiosounddyn );
	case OUT_L_JEANS_ARRAY:
	    return((FLOAT) p->l_jeans );
	case OUT_ISMALL_JEANS_ARRAY:
	    return((FLOAT) p->small_jeans );
#endif
#endif

#ifdef SIMPLESF
	case OUT_TCOOLAGAIN_ARRAY:
	    return(p->fTimeForm);
	case OUT_MSTAR_ARRAY:
	    return(p->fMassStar);
#endif
	case OUT_SPHH_ARRAY:
#endif
	case OUT_H_ARRAY:
	    return(sqrt(p->fBall2*0.25));
#ifdef SURFACEAREA
	case OUT_SURFACEAREA_ARRAY:
	    return(p->fArea);
#endif
#ifdef DRHODT
	case OUT_DIVV_T_ARRAY:
	    return(p->fDivv_t);
	case OUT_DIVV_CORRECTOR_ARRAY:
	    return(p->fDivv_Corrector);
#endif
	case OUT_IACTIVE_ARRAY:
	    return((FLOAT) p->iActive);
#ifdef COLLISIONS
	case OUT_REJECTS_ARRAY:
		/* Rejected particles indicated by their iOrder, otherwise -1 */
	    return(REJECT(p) ? p->iOrder : -1);
#endif /* COLLISIONS */
	case OUT_POS_VECTOR:
	    return(p->r[iDim]);
	case OUT_VEL_VECTOR:
	    return(pkd->dvFac*p->v[iDim]);
	case OUT_MASS_ARRAY:
	    return(p->fMass);
	case OUT_ACCELRFC_VECTOR:
	case OUT_ACCELG_VECTOR:
	case OUT_ACCEL_VECTOR:
	    return(p->a[iDim]);
#ifdef NEED_VPRED
	case OUT_VPRED_VECTOR:
	    return(p->vPred[iDim]);
#endif
#if defined(GASOLINE)
	case OUT_CURLV_VECTOR:
	    return(p->curlv[iDim]);
#ifdef NORMAL
	case OUT_NORMAL_VECTOR:
	    return(p->normal[iDim]);
#endif
#if defined(SHOCKTRACK)
	case OUT_GRADRHO_VECTOR:
	    return(p->gradrho[iDim]);
	case OUT_ACCELPRES_VECTOR:
	    return(p->aPres[iDim]);
#endif
	case OUT_ANGMOM_VECTOR:
	    return(((FLOAT *) (&(SINK_Lx(p))))[iDim]);
#endif
	default:
	    return(0.0);
		}
	}

void VecInType(PKD pkd, PARTICLE *p,int iDim,int iType,FLOAT Data)
    {
#ifdef GASOLINE
    FLOAT vTemp;
#endif

    switch (iType) {
    case OUT_DENSITY_ARRAY:
        p->fDensity = Data; return;
#ifndef NOCOOLING
    case OUT_COOL_ARRAY0:
        COOL_IN_ARRAY0(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY1:
        COOL_IN_ARRAY1(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY2:
        COOL_IN_ARRAY2(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY3:
        COOL_IN_ARRAY3(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return; /*H2*/
    case OUT_COOL_ARRAY4:
        COOL_IN_ARRAY4(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY5:
        COOL_IN_ARRAY5(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY6:
        COOL_IN_ARRAY6(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY7:
        COOL_IN_ARRAY7(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY8:
        COOL_IN_ARRAY8(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY9:
        COOL_IN_ARRAY9(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY10:
        COOL_IN_ARRAY10(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY11:
        COOL_IN_ARRAY11(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY12:
        COOL_IN_ARRAY12(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY13:
        COOL_IN_ARRAY13(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY14:
        COOL_IN_ARRAY14(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
    case OUT_COOL_ARRAY15:
        COOL_IN_ARRAY15(pkd->Cool, &p->CoolParticle,p->fMetals, Data); return;
#endif
#ifdef STARFORM
    case OUT_COOLTURNONTIME_ARRAY:
        p->fTimeCoolIsOffUntil = Data; return;
    case OUT_OXYGENMASSFRAC_ARRAY:
        p->fMFracOxygen = Data; return;
    case OUT_IRONMASSFRAC_ARRAY:
        p->fMFracIron = Data; return;
    case OUT_MASSFORM_ARRAY:
        p->fMassForm = Data; return;
#endif
    default:
        fprintf(stderr,"Bad Option to VecInType %d\n",iType);
        assert(0);
        }
    }

void VecFilename(char *achFile, int iType)
    {
	switch (iType) {
	case OUT_BIG_FILE:
#ifdef COLLISIONS
	    strncat(achFile,"ss",256);
#else
	    strncat(achFile,"tipsy",256);
#endif
	    break;
	case OUT_IORDER_ARRAY:
	    strncat(achFile,"iord",256);
        break;
    case OUT_DENSITY_ARRAY:
		strncat(achFile,"den",256);
        break;
    case OUT_DENSITYRFC_ARRAY:
		strncat(achFile,"denRFC",256);
        break;
	case OUT_COLOR_ARRAY:
#ifdef COLORCODE
	    strncat(achFile,"col",256);
#endif
	case OUT_POT_ARRAY:
        strncat(achFile,"pot",256);
        break;
	case OUT_AMAG_ARRAY:
        strncat(achFile,"amag",256);
        break;
	case OUT_RUNG_ARRAY:
        strncat(achFile,"rung",256);
        break;

	case OUT_MASS_ARRAY:
        strncat(achFile,"mass",256);
        break;
	case OUT_DT_ARRAY:
		strncat(achFile,"dt",256);
        break;
	case OUT_SPHDT_ARRAY:
		strncat(achFile,"SPHdt",256);
        break;
	case OUT_SOFT_ARRAY:
        strncat(achFile,"soft",256);
        break;
#ifdef GASOLINE
#ifdef DENSITYU
    case OUT_DENSITYU_ARRAY:
	    strncat(achFile,"denu",256);
	    break;
#endif
    case OUT_PRES_ARRAY:
	    strncat(achFile,"pres",256);
	    break;
	case OUT_TEMP_ARRAY:
		strncat(achFile,"temperature",256);
        break;
    case OUT_GASDENSITY_ARRAY:
		strncat(achFile,"GasDensity",256);
        break;
#ifdef MHD
    case OUT_MHDX_ARRAY:
      strncat(achFile,"BFieldx",256);
      break;
    case OUT_MHDY_ARRAY:
      strncat(achFile,"BFieldy",256);
      break;
    case OUT_MHDZ_ARRAY:
      strncat(achFile,"BFieldz",256);
      break;
    case OUT_MHDCURLX_ARRAY:
      strncat(achFile,"CurlBx",256);
      break;
    case OUT_MHDCURLY_ARRAY:
      strncat(achFile,"CurlBy",256);
      break;
    case OUT_MHDCURLZ_ARRAY:
      strncat(achFile,"CurlBz",256);
      break;
    case OUT_MHDDIV_ARRAY:
      strncat(achFile,"DivB",256);
      break;
    case OUT_MHDUDOTBDISS_ARRAY:
      strncat(achFile,"uDotBdiss",256);
      break;
    case OUT_MHDBDOTRHO_ARRAY:
      strncat(achFile,"BDotRho",256);
      break;
    case OUT_MHDBDOTDISS_ARRAY:
      strncat(achFile,"BDotDiss",256);
      break;
    case OUT_MHDBDOTCLEAN_ARRAY:
      strncat(achFile,"BDotBClean",256);
      break;
     case OUT_MHDBCLEAN_ARRAY:
      strncat(achFile,"BClean",256);
      break;
#endif
#ifdef THERMALCONDANI
    case OUT_GRADTX_ARRAY:
      strncat(achFile,"gradTx",256);
      break;
    case OUT_GRADTY_ARRAY:
      strncat(achFile,"gradTy",256);
      break;
    case OUT_GRADTZ_ARRAY:
      strncat(achFile,"gradTz",256);
      break;
    case OUT_UDOTANICOND_ARRAY:
    strncat(achFile,"uDotAniCond",256);
    break;
#endif
#ifdef GETVISCPARA
	case OUT_VISCCOEF_ARRAY:
	  strncat(achFile,"kinvisc",256);
	  break;
#ifdef MHD
	case OUT_RESICOEF_ARRAY:
	  strncat(achFile,"etares",256);
	  break;
#endif
#endif
#ifdef CULLENDEHNEN
	case OUT_AX_ARRAY:
	  strncat(achFile,"accx",256);
	  break;
	case OUT_AY_ARRAY:
	  strncat(achFile,"accy",256);
	  break;
	case OUT_AZ_ARRAY:
	  strncat(achFile,"accz",256);
	  break;
	case OUT_GX_ARRAY:
	  strncat(achFile,"grx",256);
	  break;
	case OUT_GY_ARRAY:
	  strncat(achFile,"gry",256);
	  break;
	case OUT_GZ_ARRAY:
	  strncat(achFile,"grz",256);
	  break;
#endif
#ifdef CHECKQUALITY
case OUT_Q1_ARRAY:
  strncat(achFile,"Q1",256);
  break;
case OUT_Q2_ARRAY:
  strncat(achFile,"Q2",256);
  break;
case OUT_E0_ARRAY:
  strncat(achFile,"E0",256);
  break;
case OUT_Q4_ARRAY:
  strncat(achFile,"Q4",256);
  break;
#endif
#ifndef NOCOOLING
	case OUT_U_ARRAY:
		strncat(achFile,"u",256);
		break;
	case OUT_TWOPHASE_ARRAY:
		strncat(achFile,"MassHot",256);
		break;
	case OUT_UNONCOOL_ARRAY:
		strncat(achFile,"uHot",256);
		break;
	case OUT_TEMPINC_ARRAY:
		strncat(achFile,"Tinc",256);
		break;
	case OUT_UDOT_ARRAY:
		strncat(achFile,"uDot",256);
		break;
	case OUT_COOL_ARRAY0:
		strncat(achFile,COOL_ARRAY0_EXT,256);
        break;
	case OUT_COOL_ARRAY1:
		strncat(achFile,COOL_ARRAY1_EXT,256);
        break;
	case OUT_COOL_ARRAY2:
		strncat(achFile,COOL_ARRAY2_EXT,256);
        break;
	case OUT_COOL_ARRAY3:
		strncat(achFile,COOL_ARRAY3_EXT,256);
        break;
	case OUT_COOL_ARRAY4:
		strncat(achFile,COOL_ARRAY4_EXT,256);
        break;
	case OUT_COOL_ARRAY5:
		strncat(achFile,COOL_ARRAY5_EXT,256);
        break;
	case OUT_COOL_ARRAY6:
		strncat(achFile,COOL_ARRAY6_EXT,256);
        break;
	case OUT_COOL_ARRAY7:
		strncat(achFile,COOL_ARRAY7_EXT,256);
        break;
	case OUT_COOL_ARRAY8:
		strncat(achFile,COOL_ARRAY8_EXT,256);
        break;
	case OUT_COOL_ARRAY9:
		strncat(achFile,COOL_ARRAY9_EXT,256);
        break;
	case OUT_COOL_ARRAY10:
		strncat(achFile,COOL_ARRAY10_EXT,256);
        break;
	case OUT_COOL_ARRAY11:
		strncat(achFile,COOL_ARRAY11_EXT,256);
        break;
	case OUT_COOL_ARRAY12:
		strncat(achFile,COOL_ARRAY12_EXT,256);
        break;
	case OUT_COOL_ARRAY13:
		strncat(achFile,COOL_ARRAY13_EXT,256);
        break;
	case OUT_COOL_ARRAY14:
		strncat(achFile,COOL_ARRAY14_EXT,256);
        break;
	case OUT_COOL_ARRAY15:
		strncat(achFile,COOL_ARRAY15_EXT,256);
        break;
#ifdef COOLING_MOLECULARH
	case OUT_CORREL_ARRAY: /*Correlation length determined from Gas Shear in terms of the Mach number -- used when calculating shielding*/
		strncat(achFile,"correL",256);
		break;
#endif
#ifdef  RADIATIVEBOX
	case OUT_COOL_LYMANWERNER_ARRAY:
		strncat(achFile,"lw",256);
		break;
#endif
	case OUT_COOL_EDOT_ARRAY:
		strncat(achFile,"eDot",256);
		break;
	case OUT_COOL_COOLING_ARRAY:
		strncat(achFile,"eCool",256);
		break;
	case OUT_COOL_HEATING_ARRAY:
		strncat(achFile,"eHeat",256);
		break;
#endif
	case OUT_BALSARASWITCH_ARRAY:
		strncat(achFile,"BSw",256);
        break;
	case OUT_ALPHA_ARRAY:
		strncat(achFile,"alpha",256);
        break;
	case OUT_DIVV_ARRAY:
        strncat(achFile,"divv",256);
        break;
	case OUT_DVXDX_ARRAY:
        strncat(achFile,"dvxdx",256);
        break;
	case OUT_DVYDX_ARRAY:
        strncat(achFile,"dvydx",256);
        break;
	case OUT_DVZDX_ARRAY:
        strncat(achFile,"dvzdx",256);
        break;
	case OUT_DVXDY_ARRAY:
        strncat(achFile,"dvxdy",256);
        break;
	case OUT_DVYDY_ARRAY:
        strncat(achFile,"dvydy",256);
        break;
	case OUT_DVZDY_ARRAY:
        strncat(achFile,"dvzdy",256);
        break;
	case OUT_DVXDZ_ARRAY:
        strncat(achFile,"dvxdz",256);
        break;
	case OUT_DVYDZ_ARRAY:
        strncat(achFile,"dvydz",256);
        break;
	case OUT_DVZDZ_ARRAY:
        strncat(achFile,"dvzdz",256);
        break;
	case OUT_DIVV_DENS_ARRAY:
        strncat(achFile,"divv_dens",256);
        break;
	case OUT_DIVVDOT_ARRAY:
        strncat(achFile,"divvdot",256);
        break;
	case OUT_DVDS_ARRAY:
        strncat(achFile,"dvds",256);
        break;
	case OUT_VSIGMAX_ARRAY:
        strncat(achFile,"vsigmax",256);
        break;
	case OUT_R_CD_ARRAY:
        strncat(achFile,"rcd",256);
        break;
	case OUT_SNORM_ARRAY:
        strncat(achFile,"snorm",256);
        break;
	case OUT_SFULL_ARRAY:
        strncat(achFile,"sfull",256);
        break;
	case OUT_DVDSONSFULL_ARRAY:
        strncat(achFile,"dvdsonsfull",256);
        break;
	case OUT_ALPHALOC_ARRAY:
        strncat(achFile,"alphaloc",256);
        break;
	case OUT_ALPHANOISE_ARRAY:
        strncat(achFile,"alphanoise",256);
        break;
	case OUT_SURFACEAREA_ARRAY:
        strncat(achFile,"area",256);
        break;
#ifdef DRHODT
	case OUT_DIVV_T_ARRAY:
        strncat(achFile,"divvt",256);
        break;
	case OUT_DIVV_CORRECTOR_ARRAY:
        strncat(achFile,"divvcorr",256);
        break;
#endif
	case OUT_IACTIVE_ARRAY:
        strncat(achFile,"ia",256);
        break;
	case OUT_CSOUND_ARRAY:
        strncat(achFile,"c",256);
        break;
	case OUT_MUMAX_ARRAY:
        strncat(achFile,"mumax",256);
        break;
	case OUT_DIVONCONH_ARRAY:
        strncat(achFile,"dch",256);
        break;
	case OUT_DIVONCONX_ARRAY:
        strncat(achFile,"dcx",256);
        break;
    case OUT_UDOTHYDRO_ARRAY:
        strncat(achFile,"uDotHydro",256);
        break;
    case OUT_PDVRFC_ARRAY:
        strncat(achFile,"PdVRFC",256);
        break;
	case OUT_UDOTPDV_ARRAY:
        strncat(achFile,"uDotPdV",256);
        break;
	case OUT_UDOTAV_ARRAY:
        strncat(achFile,"uDotAV",256);
        break;
	case OUT_UDOTDIFF_ARRAY:
        strncat(achFile,"uDotDiff",256);
        break;
#ifdef UNONCOOLDEBUG
	case OUT_UHOTDOT_ARRAY:
        strncat(achFile,"uHotDot",256);
        break;
	case OUT_UHOTDOTDIFF_ARRAY:
        strncat(achFile,"uHotDotDiff",256);
        break;
	case OUT_UHOTDOTPDV_ARRAY:
        strncat(achFile,"uHotDotPdV",256);
        break;
	case OUT_UHOTDOTCONV_ARRAY:
        strncat(achFile,"uHotDotConv",256);
        break;
	case OUT_UHOTDOTFB_ARRAY:
        strncat(achFile,"uHotDotFB",256);
        break;
#endif
    case OUT_METALS_ARRAY:
	    strncat(achFile,"Metals",256);
        break;
#ifdef DIFFUSION
	case OUT_METALSDOT_ARRAY:
	    strncat(achFile,"Metalsdot",256);
        break;
#endif
#ifdef SHOCKTRACK
	case OUT_SHOCKTRACKER_ARRAY:
        strncat(achFile,"ST",256);
        break;
	case OUT_DIVRHOV_ARRAY:
        strncat(achFile,"divrhov",256);
        break;
#endif
	case OUT_SIGMA2_ARRAY:
	    strncat(achFile,"sigma2",256);
        break;
#ifdef STARFORM
	case OUT_IGASORDER_ARRAY:
	    strncat(achFile,"igasorder",256);
        break;
	case OUT_TIMEFORM_ARRAY:
	    strncat(achFile,"timeform",256);
        break;
	case OUT_MASSFORM_ARRAY:
	    strncat(achFile,"massform",256);
        break;
#ifndef SUPERBUBBLE
	case OUT_COOLTURNONTIME_ARRAY:
	    strncat(achFile,"coolontime",256);
        break;
#endif
	case OUT_OXYGENMASSFRAC_ARRAY:
	    strncat(achFile,"OxMassFrac",256);
        break;
	case OUT_IRONMASSFRAC_ARRAY:
	    strncat(achFile,"FeMassFrac",256);
        break;
#ifdef DIFFUSION
	case OUT_OXYGENMASSFRACDOT_ARRAY:
	    strncat(achFile,"OxMassFracdot",256);
        break;
	case OUT_IRONMASSFRACDOT_ARRAY:
	    strncat(achFile,"FeMassFracdot",256);
        break;
#endif
	case OUT_UDOTFB_ARRAY:
	    strncat(achFile,"uDotFB",256);
	    break;
#ifdef CHECKSF
	case OUT_TOFF_YR_ARRAY:
	    strncat(achFile,"tcooloff",256);
	    break;
	case OUT_TCOOL_YR_ARRAY:
	    strncat(achFile,"tcoolyr",256);
	    break;
	case OUT_TDYN_YR_ARRAY:
	    strncat(achFile,"tdynyr",256);
	    break;
	case OUT_RATIOSOUNDDYN_ARRAY:
	    strncat(achFile,"rsnddyn",256);
	    break;
	case OUT_L_JEANS_ARRAY:
	    strncat(achFile,"ljeans",256);
	    break;
	case OUT_ISMALL_JEANS_ARRAY:
	    strncat(achFile,"ijeans",256);
	    break;
#endif
#endif

#ifdef SIMPLESF
	case OUT_TCOOLAGAIN_ARRAY:
		strncat(achFile,"tCoolAgain",256);
        break;
	case OUT_MSTAR_ARRAY:
		strncat(achFile,"mStar",256);
        break;
#endif
	case OUT_SPHH_ARRAY:
		strncat(achFile,"SPHH",256);
        break;
#endif
	case OUT_H_ARRAY:
		strncat(achFile,"smoothlength",256);
        break;
	case OUT_POS_VECTOR:
	    strncat(achFile,"pos",256);
        break;
	case OUT_VEL_VECTOR:
	    strncat(achFile,"vel",256);
        break;
	case OUT_ACCEL_VECTOR:
	    strncat(achFile,"acc",256);
        break;
	case OUT_ACCELG_VECTOR:
	    strncat(achFile,"accg",256);
        break;
	case OUT_ACCELRFC_VECTOR:
	    strncat(achFile,"accRFC",256);
        break;
	case OUT_CURLV_VECTOR:
	    strncat(achFile,"curl",256);
        break;
	case OUT_NORMAL_VECTOR:
	    strncat(achFile,"norm",256);
        break;
#ifdef NEED_VPRED
	case OUT_VPRED_VECTOR:
	    strncat(achFile,"vpred",256);
        break;
#endif
#if defined(GASOLINE)
#if defined(SHOCKTRACK)
	case OUT_GRADRHO_VECTOR:
		strncat(achFile,"gradrho",256);
        break;
	case OUT_ACCELPRES_VECTOR:
	    strncat(achFile,"accp",256);
        break;
#endif
	case OUT_ANGMOM_VECTOR:
	    strncat(achFile,"angmom",256);
        break;
#endif
	default:
        assert(1);
		}
	}

void pkdOutNChilada(PKD pkd,char *pszFileName,int nGasStart, int nDarkStart, int nStarStart, int iVecType, int *pnOut, float minValue[3][3], float maxValue[3][3], double duTFac, double dvFac)
    {
    FILE *gasFp = NULL, *darkFp = NULL, *starFp = NULL;
    float fOut, min[3][3], max[3][3];
    int i, iDim, nDim, headerlength;
    int nGas, nDark, nStar, nOut=0;
    char darkFileName[256],gasFileName[256],starFileName[256];
    XDR gasXdrs, darkXdrs, starXdrs;
    for(i=0;i<3;i++){
        for(iDim=0;iDim<3;iDim++){
            min[i][iDim] = 1e30;
            max[i][iDim] = -1e30;
            }
        }

    pkd->duTFac = duTFac;
    pkd->dvFac = dvFac;
    nGas = pkd->nGas; nDark = pkd->nDark; nStar = pkd->nStar;
    switch (iVecType){
	    /* Gas only floats */
    case OUT_COOLTURNONTIME_ARRAY:
    case OUT_DIVV_ARRAY:
    case OUT_DVXDX_ARRAY:
    case OUT_DVYDX_ARRAY:
    case OUT_DVZDX_ARRAY:
    case OUT_DVXDY_ARRAY:
    case OUT_DVYDY_ARRAY:
    case OUT_DVZDY_ARRAY:
    case OUT_DVXDZ_ARRAY:
    case OUT_DVYDZ_ARRAY:
    case OUT_DVZDZ_ARRAY:
    case OUT_DIVV_DENS_ARRAY:
    case OUT_DIVVDOT_ARRAY:
    case OUT_DVDS_ARRAY:
    case OUT_VSIGMAX_ARRAY:
    case OUT_R_CD_ARRAY:
    case OUT_SNORM_ARRAY:
    case OUT_SFULL_ARRAY:
    case OUT_DVDSONSFULL_ARRAY:
    case OUT_ALPHALOC_ARRAY:
    case OUT_ALPHANOISE_ARRAY:
    case OUT_TCOOLAGAIN_ARRAY:
    case OUT_MSTAR_ARRAY:
    case OUT_COOL_ARRAY0:
    case OUT_COOL_ARRAY1:
    case OUT_COOL_ARRAY2:
    case OUT_COOL_ARRAY3:
    case OUT_COOL_ARRAY4:
    case OUT_COOL_ARRAY5:
    case OUT_COOL_ARRAY6:
    case OUT_COOL_ARRAY7:
    case OUT_COOL_ARRAY8:
    case OUT_COOL_ARRAY9:
    case OUT_COOL_ARRAY10:
    case OUT_COOL_ARRAY11:
    case OUT_COOL_ARRAY12:
    case OUT_COOL_ARRAY13:
    case OUT_COOL_ARRAY14:
    case OUT_COOL_ARRAY15:
#ifdef COOLING_MOLECULARH
    case OUT_CORREL_ARRAY:
#endif
#ifdef RADIATIVEBOX
    case OUT_COOL_LYMANWERNER_ARRAY:
#endif
    case OUT_SPHH_ARRAY:
    case OUT_TEMP_ARRAY:
    case OUT_GASDENSITY_ARRAY:
    case OUT_UDOTHYDRO_ARRAY:
    case OUT_UDOTPDV_ARRAY:
    case OUT_UDOTAV_ARRAY:
    case OUT_UDOTDIFF_ARRAY:
#ifdef UNONCOOLDEBUG
    case OUT_UHOTDOT_ARRAY:
    case OUT_UHOTDOTDIFF_ARRAY:
    case OUT_UHOTDOTPDV_ARRAY:
    case OUT_UHOTDOTCONV_ARRAY:
    case OUT_UHOTDOTFB_ARRAY:
#endif
    case OUT_PRES_ARRAY:
    case OUT_BALSARASWITCH_ARRAY:
    case OUT_SPHDT_ARRAY:
    case OUT_CSOUND_ARRAY:
    case OUT_MUMAX_ARRAY:
    case OUT_METALSDOT_ARRAY:
    case OUT_OXYGENMASSFRACDOT_ARRAY:
    case OUT_IRONMASSFRACDOT_ARRAY:
    case OUT_COOL_EDOT_ARRAY:
    case OUT_COOL_COOLING_ARRAY:
    case OUT_COOL_HEATING_ARRAY:
    case OUT_CURLV_VECTOR:
        nDark=nStar=0;
        break;
    case OUT_IGASORDER_ARRAY:
    case OUT_TIMEFORM_ARRAY:
    case OUT_MASSFORM_ARRAY:
    case OUT_ANGMOM_VECTOR:
        nGas=nDark=0;
        break;
    case OUT_METALS_ARRAY:
    case OUT_OXYGENMASSFRAC_ARRAY:
    case OUT_IRONMASSFRAC_ARRAY:
    case OUT_UDOTFB_ARRAY:
        nDark=0;
        break;
        }

    /*
     * N-Chilada has a 28 byte header (see FieldHeader in
     * structures/tree_xdr.h): 4 byte magic + 8 byte time + 8 byte
     * numParticles + 4 byte dimension + 4 byte data type code.  The
     * header is followed by a min and max field 6 numbers for vectors
     * and 2 for scalars.
     */
    /*
     * XXX WARNING:  Seeks below are assuming sizeof(float) = 4, since
     * RFCs for XDR assume 4 byte floats.
     */

    if(iVecType > OUT_1D3DSPLIT) {
        nDim=3;
        headerlength = 28 + 6*4;
        }
    else {
        nDim=1;
        headerlength = 28 + 2*4;
        }

    if(nGas){
        sprintf(gasFileName, "%s/gas/",pszFileName);
        VecFilename(gasFileName,iVecType);
        gasFp = fopen(gasFileName,"r+");
        mdlassert(pkd->mdl,gasFp != NULL);
        xdrstdio_create(&gasXdrs,gasFp,XDR_ENCODE);
        }

    if(nDark){
        sprintf(darkFileName, "%s/dark/",pszFileName);
        VecFilename(darkFileName,iVecType);
        darkFp = fopen(darkFileName,"r+");
        mdlassert(pkd->mdl,darkFp != NULL);
        xdrstdio_create(&darkXdrs,darkFp,XDR_ENCODE);
        }

    if(nStar){
        sprintf(starFileName, "%s/star/",pszFileName);
        VecFilename(starFileName,iVecType);
        starFp = fopen(starFileName,"r+");
        mdlassert(pkd->mdl,starFp != NULL);
        xdrstdio_create(&starXdrs,starFp,XDR_ENCODE);
        }

    switch (iVecType) {
#ifdef STARFORM
	case OUT_IGASORDER_ARRAY:
        if(nStar) {
            pkdGenericSeek(pkd,starFp, nStarStart, headerlength, 4);
            for (i=0;i<pkd->nLocal;++i) {
                if (pkdIsStar(pkd,&pkd->pStore[i])) {
                    xdr_int(&starXdrs,&(pkd->pStore[i].iGasOrder));
                    nOut++;
                    min[2][0] = (pkd->pStore[i].iOrder > min[2][0]) ? min[2][0]: pkd->pStore[i].iOrder;
                    max[2][0] = (pkd->pStore[i].iOrder < max[2][0]) ? max[2][0]: pkd->pStore[i].iOrder;
                    }
                }
            }
	    break;
#endif
    case OUT_IORDER_ARRAY:
        if(nGas) pkdGenericSeek(pkd,gasFp, nGasStart,headerlength, 4);
        if(nDark) pkdGenericSeek(pkd,darkFp, nDarkStart,headerlength, 4);
        if(nStar) pkdGenericSeek(pkd,starFp, nStarStart,headerlength, 4);
        for (i=0;i<pkd->nLocal;++i) {
            if (nGas && pkdIsGas(pkd,&pkd->pStore[i])) {
                xdr_int(&gasXdrs,&(pkd->pStore[i].iOrder));
                nOut++;
                min[0][0] = (pkd->pStore[i].iOrder > min[0][0]) ? min[0][0]: pkd->pStore[i].iOrder;
                max[0][0] = (pkd->pStore[i].iOrder < max[0][0]) ? max[0][0]: pkd->pStore[i].iOrder;
                }
            if (nDark && pkdIsDark(pkd,&pkd->pStore[i])) {
                xdr_int(&darkXdrs,&(pkd->pStore[i].iOrder));
                nOut++;
                min[1][0] = (pkd->pStore[i].iOrder > min[1][0]) ? min[1][0]: pkd->pStore[i].iOrder;
                max[1][0] = (pkd->pStore[i].iOrder < max[1][0]) ? max[1][0]: pkd->pStore[i].iOrder;
                }
            if (nStar && pkdIsStar(pkd,&pkd->pStore[i])) {
                xdr_int(&starXdrs,&(pkd->pStore[i].iOrder));
                nOut++;
                min[2][0] = (pkd->pStore[i].iOrder > min[2][0]) ? min[2][0]: pkd->pStore[i].iOrder;
                max[2][0] = (pkd->pStore[i].iOrder < max[2][0]) ? max[2][0]: pkd->pStore[i].iOrder;
                }
            }
        break;
#ifdef GASOLINE
    case OUT_TEMP_ARRAY:
    case OUT_GASDENSITY_ARRAY:
        if(nGas) {
            pkdGenericSeek(pkd,gasFp, nGasStart,headerlength,sizeof(pkd->pStore[i].iOrder));
            for (i=0;i<pkd->nLocal;++i) {
                if (pkdIsGas(pkd,&pkd->pStore[i])) {
                    fOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    min[0][0] = (fOut > min[0][0]) ? min[0][0]: fOut;
                    max[0][0] = (fOut < max[0][0]) ? max[0][0]: fOut;
                    xdr_float(&gasXdrs,&fOut);
                    nOut++;
                    }
                }
            }
        break;
#endif
    default:
        if(nGas) pkdGenericSeek(pkd,gasFp, ((long)nGasStart)*nDim,
            headerlength,sizeof(fOut));
        if(nDark) pkdGenericSeek(pkd,darkFp, ((long)nDarkStart)*nDim,
            headerlength,sizeof(fOut));
        if(nStar) pkdGenericSeek(pkd,starFp, ((long)nStarStart)*nDim,
            headerlength,sizeof(fOut));
        for (i=0;i<pkd->nLocal;++i) {
            for (iDim = 0; iDim <nDim; iDim++) {
                fOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                if (nGas && pkdIsGas(pkd,&pkd->pStore[i])) {
                    xdr_float(&gasXdrs,&fOut);
                    nOut++;
                    min[0][iDim] = (fOut > min[0][iDim]) ? min[0][iDim]: fOut;
                    max[0][iDim] = (fOut < max[0][iDim]) ? max[0][iDim]: fOut;
                    }
                if (nDark && pkdIsDark(pkd,&pkd->pStore[i])) {
                    xdr_float(&darkXdrs,&fOut);
                    nOut++;
                    min[1][iDim] = (fOut > min[1][iDim]) ? min[1][iDim]: fOut;
                    max[1][iDim] = (fOut < max[1][iDim]) ? max[1][iDim]: fOut;
                    }
                if (nStar && pkdIsStar(pkd,&pkd->pStore[i])) {
                    xdr_float(&starXdrs,&fOut);
                    nOut++;
                    min[2][iDim] = (fOut > min[2][iDim]) ? min[2][iDim]: fOut;
                    max[2][iDim] = (fOut < max[2][iDim]) ? max[2][iDim]: fOut;
                    }
                }
            }
        }
    if(nDark) xdr_destroy(&darkXdrs);
    if(nGas) xdr_destroy(&gasXdrs);
    if(nStar) xdr_destroy(&starXdrs);
    if(nDark) fclose(darkFp);
    if(nGas) fclose(gasFp);
    if(nStar) fclose(starFp);
    for(i=0;i<3;i++){
        for(iDim=0;iDim<nDim;iDim++){
            minValue[i][iDim] = min[i][iDim];
            maxValue[i][iDim] = max[i][iDim];
            }
        }
    *pnOut = nOut;
    }

void xdr_FLOAT(XDR *xdrs, FLOAT *fIn)
    {
#ifdef SINGLE
    xdr_float(xdrs, fIn);
#else
    xdr_double(xdrs, fIn);
#endif
    }

void pkdOutVector(PKD pkd,char *pszFileName,int nStart, int iDim,int iVecType,int iBinaryOutput, int N, int bStandard)
    {
    FILE *fp;
    FLOAT fOut;
    float FloatOut;
    double DoubleOut;
    int IntOut;
    long LongOut;
    int i;
    char vecFileName[256];

    strcpy(vecFileName, pszFileName);
    VecFilename(vecFileName,iVecType);
    if(iBinaryOutput){
        fp = fopen(vecFileName,"r+");
        mdlassert(pkd->mdl,fp != NULL);
        } else {
        fp = fopen(vecFileName,"a");
        mdlassert(pkd->mdl,fp != NULL);
        }
    /*
    ** Write Vector Elements!
    */
    if( bStandard && iBinaryOutput ){
        XDR xdrs;
        xdrstdio_create(&xdrs,fp,XDR_ENCODE);
        switch (iBinaryOutput) {
        case 1:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, ((long)nStart)*3,sizeof(int),sizeof(FloatOut));
                for (i=0;i<pkd->nLocal;++i) {
                    FloatOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    xdr_float(&xdrs,&FloatOut);
                    FloatOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                    xdr_float(&xdrs,&FloatOut);
                    FloatOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                    xdr_float(&xdrs,&FloatOut);
                    }
                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(IntOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_int(&xdrs,&(pkd->pStore[i].iGasOrder));
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp,nStart, sizeof(int), sizeof(IntOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_int(&xdrs,&(pkd->pStore[i].iOrder));
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, ((long)N)*iDim+nStart,sizeof(int),sizeof(FloatOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        FloatOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                        xdr_float(&xdrs,&FloatOut);
                        }
                    }
                }
            break;
        case 2:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, ((long)nStart)*3,sizeof(int),sizeof(DoubleOut));
                for (i=0;i<pkd->nLocal;++i) {
                    DoubleOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    xdr_double(&xdrs,&DoubleOut);
                    DoubleOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                    xdr_double(&xdrs,&DoubleOut);
                    DoubleOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                    xdr_double(&xdrs,&DoubleOut);
                    }

                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        LongOut = pkd->pStore[i].iGasOrder;
                        xdr_longlong_t(&xdrs,(quad_t *) &LongOut);
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        LongOut = pkd->pStore[i].iOrder;
                        xdr_longlong_t(&xdrs,(quad_t *) &LongOut);
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, ((long)N)*iDim+nStart,sizeof(int),sizeof(DoubleOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        DoubleOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                        xdr_double(&xdrs,&DoubleOut);
                        }
                    }
                }
            break;
        case 3:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, ((long)nStart)*3,sizeof(int),sizeof(fOut));
                for (i=0;i<pkd->nLocal;++i) {
                    fOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    xdr_FLOAT(&xdrs,&fOut);
                    fOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                    xdr_FLOAT(&xdrs,&fOut);
                    fOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                    xdr_FLOAT(&xdrs,&fOut);
                    }
                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart, sizeof(int), sizeof(pkd->pStore[i].iGasOrder));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_int(&xdrs,&(pkd->pStore[i].iGasOrder));
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(pkd->pStore[i].iOrder));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_int(&xdrs,&(pkd->pStore[i].iOrder));
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, ((long)N)*iDim+nStart,sizeof(int),sizeof(fOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        fOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                        xdr_FLOAT(&xdrs,&fOut);
                        }
                    }
                }
            break;
            }
        }
    else {
        switch (iBinaryOutput) {
        case 0:
            if (iDim < 0) {
                for (i=0;i<pkd->nLocal;++i) {
                    fOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    fprintf(fp,"%.14g\n",fOut);
                    fOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                    fprintf(fp,"%.14g\n",fOut);
                    fOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                    fprintf(fp,"%.14g\n",fOut);
                    }
                } else {
                for (i=0;i<pkd->nLocal;++i) {
                    fOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                    fprintf(fp,"%.14g\n",fOut);
                    }
                }
            break;
        case 1:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, ((long)nStart)*3,sizeof(int),sizeof(FloatOut));
                for (i=0;i<pkd->nLocal;++i) {
                    FloatOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    fwrite(&FloatOut, sizeof(float), 1, fp );
                    FloatOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                    fwrite(&FloatOut, sizeof(float), 1, fp );
                    FloatOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                    fwrite(&FloatOut, sizeof(float), 1, fp );
                    }
                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(IntOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        fwrite(&(pkd->pStore[i].iGasOrder), sizeof(IntOut), 1, fp );
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp,nStart, sizeof(int), sizeof(IntOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        fwrite(&(pkd->pStore[i].iOrder), sizeof(IntOut), 1, fp );
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, ((long)N)*iDim+nStart,sizeof(int),sizeof(FloatOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        FloatOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                        fwrite(&FloatOut, sizeof(float), 1, fp );
                        }
                    }
                }
            break;
        case 2:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, ((long)nStart)*3,sizeof(int),sizeof(DoubleOut));
                for (i=0;i<pkd->nLocal;++i) {
                    DoubleOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    fwrite(&DoubleOut, sizeof(double), 1, fp );
                    DoubleOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                    fwrite(&DoubleOut, sizeof(double), 1, fp );
                    DoubleOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                    fwrite(&DoubleOut, sizeof(double), 1, fp );
                    }

                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        LongOut = pkd->pStore[i].iGasOrder;
                        fwrite(&LongOut, sizeof(LongOut), 1, fp );
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        LongOut = pkd->pStore[i].iOrder;
                        fwrite(&LongOut, sizeof(LongOut), 1, fp );
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, ((long)N)*iDim+nStart,sizeof(int),sizeof(DoubleOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        DoubleOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                        fwrite(&DoubleOut, sizeof(double), 1, fp );
                        }
                    }
                }
            break;
        case 3:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, ((long)nStart)*3,sizeof(int),sizeof(fOut));
                for (i=0;i<pkd->nLocal;++i) {
                    fOut = VecType(pkd, &pkd->pStore[i],0,iVecType);
                    fwrite(&fOut, sizeof(FLOAT), 1, fp );
                    fOut = VecType(pkd, &pkd->pStore[i],1,iVecType);
                    fwrite(&fOut, sizeof(FLOAT), 1, fp );
                    fOut = VecType(pkd, &pkd->pStore[i],2,iVecType);
                    fwrite(&fOut, sizeof(FLOAT), 1, fp );
                    }
                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart, sizeof(int), sizeof(pkd->pStore[i].iGasOrder));
                    for (i=0;i<pkd->nLocal;++i) {
                        fwrite(&(pkd->pStore[i].iGasOrder), sizeof(pkd->pStore[i].iGasOrder), 1, fp );
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(pkd->pStore[i].iOrder));
                    for (i=0;i<pkd->nLocal;++i) {
                        fwrite(&(pkd->pStore[i].iOrder), sizeof(pkd->pStore[i].iOrder), 1, fp );
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, ((long)N)*iDim+nStart,sizeof(int),sizeof(fOut));
                    for (i=0;i<pkd->nLocal;++i) {
                        fOut = VecType(pkd, &pkd->pStore[i],iDim,iVecType);
                        fwrite(&fOut, sizeof(FLOAT), 1, fp );
                        }
                    }
                }
            break;
            }
        }

	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutVector: could not close file");
		exit(1);
		}
	}

void pkdInVector(PKD pkd,char *pszFileName,int nStart, int nLocal, int iDim,int iVecType,int iBinaryInput, int N, int bStandard)
    {
    FILE *fp;
    FLOAT fIn;
    float FloatIn;
    double DoubleIn;
    int IntIn;
    long LongIn;
    int i;
    char vecFileName[256];
#ifdef CHECKSANITY
    int nProb=0;
#endif

    assert(pkd->nLocal == nLocal);

    strcpy(vecFileName, pszFileName);
    strcat(vecFileName, ".");
    VecFilename(vecFileName,iVecType);
//    printf("id %d opening %s\n",pkd->idSelf,vecFileName);
    fp = fopen(vecFileName,"r");
    mdlassert(pkd->mdl,fp != NULL);
/*
** Read Vector Elements!
*/
    if( bStandard && iBinaryInput ){
        XDR xdrs;
        xdrstdio_create(&xdrs,fp,XDR_DECODE);
        switch (iBinaryInput) {
        case 1:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(FloatIn));
                for (i=0;i<pkd->nLocal;++i) {
                    xdr_float(&xdrs,&FloatIn);
                    VecInType(pkd, &pkd->pStore[i],0,iVecType,FloatIn);
                    xdr_float(&xdrs,&FloatIn);
                    VecInType(pkd, &pkd->pStore[i],1,iVecType,FloatIn);
                    xdr_float(&xdrs,&FloatIn);
                    VecInType(pkd, &pkd->pStore[i],2,iVecType,FloatIn);
                    }
                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(IntIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_int(&xdrs,&(pkd->pStore[i].iGasOrder));
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp,nStart, sizeof(int), sizeof(IntIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_int(&xdrs,&(pkd->pStore[i].iOrder));
/*			    if (i<10) fprintf(stderr,"Proc %d %d read %d,  iOrder %d\n",pkd->idSelf,nStart,i,pkd->pStore[i].iOrder); */
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(FloatIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_float(&xdrs,&FloatIn);
                        CHECKSANEFIX( nProb, FloatIn, ISFINITE(FloatIn) );
                        VecInType(pkd, &pkd->pStore[i],2,iVecType,FloatIn);
                        }
                    }
                }
            break;
        case 2:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(DoubleIn));
                for (i=0;i<pkd->nLocal;++i) {
                    xdr_double(&xdrs,&DoubleIn);
                    VecInType(pkd, &pkd->pStore[i],0,iVecType,DoubleIn);
                    xdr_double(&xdrs,&DoubleIn);
                    VecInType(pkd, &pkd->pStore[i],1,iVecType,DoubleIn);
                    xdr_double(&xdrs,&DoubleIn);
                    VecInType(pkd, &pkd->pStore[i],2,iVecType,DoubleIn);
                    }

                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_longlong_t(&xdrs,(quad_t*) &LongIn);
                        pkd->pStore[i].iGasOrder = LongIn;
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_longlong_t(&xdrs,(quad_t*) &LongIn);
                        pkd->pStore[i].iOrder = LongIn;
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(DoubleIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_double(&xdrs,&DoubleIn);
                        VecInType(pkd, &pkd->pStore[i],iDim,iVecType,DoubleIn);
                        }
                    }
                }
            break;
        case 3:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(fIn));
                for (i=0;i<pkd->nLocal;++i) {
                    xdr_FLOAT(&xdrs,&fIn);
                    VecInType(pkd, &pkd->pStore[i],0,iVecType,fIn);
                    xdr_FLOAT(&xdrs,&fIn);
                    VecInType(pkd, &pkd->pStore[i],1,iVecType,fIn);
                    xdr_FLOAT(&xdrs,&fIn);
                    VecInType(pkd, &pkd->pStore[i],2,iVecType,fIn);
                    }
                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart, sizeof(int), sizeof(pkd->pStore[i].iGasOrder));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_int(&xdrs,&(pkd->pStore[i].iGasOrder));
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(pkd->pStore[i].iOrder));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_int(&xdrs,&(pkd->pStore[i].iOrder));
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(fIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        xdr_FLOAT(&xdrs,&fIn);
                        VecInType(pkd, &pkd->pStore[i],iDim,iVecType,fIn);
                        }
                    }
                }
            break;
            }
        }
    else {
        switch (iBinaryInput) {
        case 0:
            if (iDim < 0) {
                for (i=0;i<pkd->nLocal;++i) {
                    fscanf(fp,"%lg\n",&fIn);
                    VecInType(pkd, &pkd->pStore[i],0,iVecType,fIn);
                    fscanf(fp,"%lg\n",&fIn);
                    VecInType(pkd, &pkd->pStore[i],1,iVecType,fIn);
                    fscanf(fp,"%lg\n",&fIn);
                    VecInType(pkd, &pkd->pStore[i],2,iVecType,fIn);
                    }
                } else {
                for (i=0;i<pkd->nLocal;++i) {
                    fscanf(fp,"%lg\n",&fIn);
                    VecInType(pkd, &pkd->pStore[i],iDim,iVecType,fIn);
                    }
                }
            break;
        case 1:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(FloatIn));
                for (i=0;i<pkd->nLocal;++i) {
                    fread(&FloatIn, sizeof(float), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],0,iVecType, FloatIn);
                    fread(&FloatIn, sizeof(float), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],1, iVecType, FloatIn);
                    fread(&FloatIn, sizeof(float), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],2, iVecType, FloatIn);
                    }
                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(IntIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&(pkd->pStore[i].iGasOrder), sizeof(IntIn), 1, fp );
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp,nStart, sizeof(int), sizeof(IntIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&(pkd->pStore[i].iOrder), sizeof(IntIn), 1, fp );
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(FloatIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&FloatIn, sizeof(float), 1, fp );
                        CHECKSANEFIX( nProb, FloatIn, ISFINITE(FloatIn) );
                        VecInType(pkd, &pkd->pStore[i],iDim,iVecType,FloatIn);
                        }
                    }
                }
            break;
        case 2:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(DoubleIn));
                for (i=0;i<pkd->nLocal;++i) {
                    fread(&DoubleIn, sizeof(double), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],0,iVecType, DoubleIn);
                    fread(&DoubleIn, sizeof(double), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],1,iVecType, DoubleIn);
                    fread(&DoubleIn, sizeof(double), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],2,iVecType, DoubleIn);
                    }

                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&LongIn, sizeof(LongIn), 1, fp );
                        pkd->pStore[i].iGasOrder = LongIn;
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(LongIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&LongIn, sizeof(LongIn), 1, fp );
                        pkd->pStore[i].iOrder = LongIn;
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(DoubleIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&DoubleIn, sizeof(double), 1, fp );
                        VecInType(pkd, &pkd->pStore[i],iDim,iVecType,DoubleIn);
                        }
                    }
                }
            break;
        case 3:
            if (iDim < 0) {
                pkdGenericSeek(pkd,fp, nStart*3,sizeof(int),sizeof(fIn));
                for (i=0;i<pkd->nLocal;++i) {
                    fread(&fIn, sizeof(FLOAT), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],0,iVecType, fIn);
                    fread(&fIn, sizeof(FLOAT), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],1,iVecType, fIn);
                    fread(&fIn, sizeof(FLOAT), 1, fp );
                    VecInType(pkd, &pkd->pStore[i],2,iVecType, fIn);
                    }
                } else {
                switch (iVecType) {
#ifdef STARFORM
                case OUT_IGASORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart, sizeof(int), sizeof(pkd->pStore[i].iGasOrder));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&(pkd->pStore[i].iGasOrder), sizeof(pkd->pStore[i].iGasOrder), 1, fp );
                        }
                    break;
#endif
                case OUT_IORDER_ARRAY:
                    pkdGenericSeek(pkd,fp, nStart,sizeof(int),sizeof(pkd->pStore[i].iOrder));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&(pkd->pStore[i].iOrder), sizeof(pkd->pStore[i].iOrder), 1, fp );
                        }
                    break;
                default:
                    pkdGenericSeek(pkd,fp, N*iDim+nStart,sizeof(int),sizeof(fIn));
                    for (i=0;i<pkd->nLocal;++i) {
                        fread(&fIn, sizeof(FLOAT), 1, fp );
                        VecInType(pkd, &pkd->pStore[i],iDim,iVecType, fIn);
                        }
                    }
                }
            break;
            }
        }

    i = fclose(fp);
    if (i != 0) {
        perror("pkdInVector: could not close file");
        exit(1);
        }

#ifdef CHECKSANITY
    if (nProb > 0) {
        fprintf(stderr,"pkdInVector id %d: %d Inf/Nan Problem of %d Read\n",pkd->idSelf,nProb,pkd->nLocal);
        }
#endif

    }
