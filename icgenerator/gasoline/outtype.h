
#ifndef OUTTYPE_HINCLUDED
#define OUTTYPE_HINCLUDED

#include "pkd.h"

#define NUMOUTPUTS  128

enum DataTypeCode {
	INT8 = 1,
	UINT8,
	INT16,
	UINT16,
	INT32,
	UINT32,
	INT64,
	UINT64,
	FLOAT32,
	FLOAT64
};

enum outtype_arraytype {
	OUT_NULL,
    OUT_BIG_FILE,
	OUT_COLOR_ARRAY,
	OUT_DENSITY_ARRAY,
	OUT_DENSITYU_ARRAY,
	OUT_DENSITYRFC_ARRAY,
	OUT_PRES_ARRAY,
	OUT_POT_ARRAY,
	OUT_AMAG_ARRAY,
	OUT_IMASS_ARRAY,
    OUT_MASS_ARRAY,
	OUT_RUNG_ARRAY,
	OUT_SPHH_ARRAY,
	OUT_U_ARRAY,
	OUT_MHDX_ARRAY,
	OUT_MHDY_ARRAY,
	OUT_MHDZ_ARRAY,
	OUT_MHDCURLX_ARRAY,
	OUT_MHDCURLY_ARRAY,
	OUT_MHDCURLZ_ARRAY,
	OUT_MHDDIV_ARRAY,
	OUT_MHDUDOTBDISS_ARRAY,
	OUT_MHDBDOTRHO_ARRAY,
	OUT_MHDBDOTDISS_ARRAY,
	OUT_MHDBDOTCLEAN_ARRAY,
	OUT_MHDBCLEAN_ARRAY,
	OUT_VISCCOEF_ARRAY,
	OUT_RESICOEF_ARRAY,
	OUT_AX_ARRAY,
	OUT_AY_ARRAY,
	OUT_AZ_ARRAY,
	OUT_GX_ARRAY,
	OUT_GY_ARRAY,
	OUT_GZ_ARRAY,
	OUT_Q1_ARRAY,
	OUT_Q2_ARRAY,
	OUT_E0_ARRAY,
	OUT_Q4_ARRAY,
	OUT_UDOTANICOND_ARRAY,
	OUT_GRADTX_ARRAY,
	OUT_GRADTY_ARRAY,
	OUT_GRADTZ_ARRAY,
    OUT_TWOPHASE_ARRAY,
	OUT_UNONCOOL_ARRAY,
	OUT_TEMPINC_ARRAY,
    OUT_TEMP_ARRAY,
	OUT_COOL_EDOT_ARRAY,
	OUT_COOL_COOLING_ARRAY,
	OUT_COOL_HEATING_ARRAY,
        OUT_GASDENSITY_ARRAY,
	OUT_HSMDIVV_ARRAY,
	OUT_UDOT_ARRAY,
	OUT_COOL_ARRAY0,
	OUT_COOL_ARRAY1,
	OUT_COOL_ARRAY2,
	OUT_COOL_ARRAY3,
	OUT_COOL_ARRAY4,
	OUT_COOL_ARRAY5,
	OUT_COOL_ARRAY6,
	OUT_COOL_ARRAY7,
	OUT_COOL_ARRAY8,
	OUT_COOL_ARRAY9,
	OUT_COOL_ARRAY10,
	OUT_COOL_ARRAY11,
	OUT_COOL_ARRAY12,
	OUT_COOL_ARRAY13,
	OUT_COOL_ARRAY14,
	OUT_COOL_ARRAY15,
	OUT_COOL_LYMANWERNER_ARRAY,
	OUT_CORREL_ARRAY,
	OUT_BALSARASWITCH_ARRAY,
	OUT_ALPHA_ARRAY,
	OUT_DIVV_ARRAY,
	OUT_DVXDX_ARRAY,
	OUT_DVYDX_ARRAY,
	OUT_DVZDX_ARRAY,
	OUT_DVXDY_ARRAY,
	OUT_DVYDY_ARRAY,
	OUT_DVZDY_ARRAY,
	OUT_DVXDZ_ARRAY,
	OUT_DVYDZ_ARRAY,
	OUT_DVZDZ_ARRAY,
	OUT_DIVV_DENS_ARRAY,
	OUT_DVDS_ARRAY,
	OUT_DIVVDOT_ARRAY,
    OUT_ALPHALOC_ARRAY,
    OUT_ALPHANOISE_ARRAY,
    OUT_VSIGMAX_ARRAY,
    OUT_R_CD_ARRAY,
    OUT_SNORM_ARRAY,
    OUT_SFULL_ARRAY,
    OUT_DVDSONSFULL_ARRAY,
	OUT_CSOUND_ARRAY,
	OUT_MUMAX_ARRAY,
	OUT_SHOCKTRACKER_ARRAY,
	OUT_DIVONCONH_ARRAY,
	OUT_DIVONCONX_ARRAY,
	OUT_DIVRHOV_ARRAY,
	OUT_UDOTHYDRO_ARRAY,
	OUT_PDVRFC_ARRAY,
	OUT_UDOTPDV_ARRAY,
	OUT_UDOTAV_ARRAY,
    OUT_UDOTDIFF_ARRAY,
    OUT_UHOTDOT_ARRAY,
    OUT_UHOTDOTDIFF_ARRAY,
    OUT_UHOTDOTPDV_ARRAY,
    OUT_UHOTDOTCONV_ARRAY,
    OUT_UHOTDOTFB_ARRAY,
	OUT_PDVSN_ARRAY,
	OUT_USN_ARRAY,
	OUT_METALS_ARRAY,
	OUT_METALSDOT_ARRAY,
    OUT_SIGMA2_ARRAY,
	OUT_IGASORDER_ARRAY,
    OUT_TIMEFORM_ARRAY,
	OUT_MASSFORM_ARRAY,
	OUT_COOLTURNONTIME_ARRAY,
    OUT_OXYGENMASSFRAC_ARRAY,
    OUT_IRONMASSFRAC_ARRAY,
    OUT_OXYGENMASSFRACDOT_ARRAY,
    OUT_IRONMASSFRACDOT_ARRAY,
	OUT_UDOTFB_ARRAY,
	OUT_TCOOLAGAIN_ARRAY,
	OUT_MSTAR_ARRAY,
	OUT_SOFT_ARRAY,
	OUT_IORDER_ARRAY,
	OUT_H_ARRAY,
        OUT_SPHDT_ARRAY,
	OUT_DT_ARRAY,
	OUT_REJECTS_ARRAY,
	OUT_TOFF_YR_ARRAY,
	OUT_TCOOL_YR_ARRAY,
	OUT_TDYN_YR_ARRAY,
	OUT_RATIOSOUNDDYN_ARRAY,
	OUT_L_JEANS_ARRAY,
	OUT_ISMALL_JEANS_ARRAY,
	OUT_SURFACEAREA_ARRAY,
	OUT_DIVV_T_ARRAY,
	OUT_DIVV_CORRECTOR_ARRAY,
	OUT_IACTIVE_ARRAY,
        OUT_1D3DSPLIT,  /* NOTICE!!
                       * Everything above here is 1D
                       * Everything below here is 3D
                       */
	OUT_POS_VECTOR,
	OUT_VEL_VECTOR,
	OUT_ACCEL_VECTOR,
        OUT_ACCELG_VECTOR,
	OUT_ACCELRFC_VECTOR,
	OUT_CURLV_VECTOR,
	OUT_NORMAL_VECTOR,
	OUT_VPRED_VECTOR,
	OUT_GRADRHO_VECTOR,
	OUT_ACCELPRES_VECTOR,
	OUT_ANGMOM_VECTOR,
	OUT_END_OF_LIST
	};

/*void pkdOutArray(PKD pkd,char *pszFileName,int nStart, int iArrType, int iBinaryOutput);*/
void VecFilename(char *achFile, int iType);
FLOAT VecType(PKD pkd, PARTICLE *p,int iDim,int iType);
void pkdOutVector(PKD pkd,char *pszFileName,int nStart, int iDim,int iVecType,int iBinaryOutput, int N, int bStandard);
void pkdInVector(PKD pkd,char *pszFileName,int nStart, int nLocal, int iDim,int iVecType,int iBinaryInput, int N, int bStandard);
void pkdGenericSeek(PKD pkd,FILE *fp,long nStart,int iHeader, int iElement);
void pkdOutNChilada(PKD pkd,char *pszFileName,int nGasStart, int nDarkStart, int nStarStart, int iVecType, int *pnOut, float minValue[3][3], float maxValue[3][3], double duTFac, double dvFac);
#endif
