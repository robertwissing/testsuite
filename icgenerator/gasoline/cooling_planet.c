
#ifdef GASOLINE
#ifndef NOCOOLING

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* General accuracy target */
#define EPS 1e-5
#define MAXABUNDITERATIONS 20
/* Accuracy for Temperature estimation for E,rho assuming eqm abundances */
#define EPSTEMP 1e-5
#define ETAMINTIMESTEP 1e-4

#include "stiff.h"
#if defined(COOLDEBUG) || defined(STARFORM) || defined(SIMPLESF)
#include "pkd.h"
#else
#include "cooling.h"
#endif
#include "outtype.h"

/* Integrator constants */

/* When to use to just do a first order step and quit */
/* Good value to use I think */
#define ECHANGEFRACSMALL 1e-4  

/* Accuracy target for intergrators */
#define EPSINTEG  1e-5
#define MAXINTEGITS 20000

#define EMUL (1.01)

COOL *CoolInit( )
{
  COOL *cl;
  cl = (COOL *) malloc(sizeof(COOL));
  assert(cl!=NULL);

  cl->nTableRead = 0; /* Internal Tables read from Files */

  cl->DerivsData = malloc(sizeof(clDerivsData));
  assert(cl->DerivsData != NULL);
  ((clDerivsData *) (cl->DerivsData))->IntegratorContext = 
    StiffInit( EPSINTEG, 1, cl->DerivsData, clDerivs);
  
  return cl;
}

void CoolFinalize(COOL *cl ) 
{
  StiffFinalize( ((clDerivsData *) (cl->DerivsData))->IntegratorContext );
  free(cl->DerivsData);
  free(cl);
}

void clInitConstants( COOL *cl, double dGmPerCcUnit, double dComovingGmPerCcUnit, 
		double dErgPerGmUnit, double dSecUnit, double dKpcUnit, COOLPARAM CoolParam) 
{
  assert(cl!=NULL);
  cl->dGmPerCcUnit = dGmPerCcUnit;
  cl->dComovingGmPerCcUnit = dComovingGmPerCcUnit;
  cl->dErgPerGmUnit = dErgPerGmUnit;
  cl->dSecUnit = dSecUnit;
  cl->dErgPerGmPerSecUnit = cl->dErgPerGmUnit / cl->dSecUnit;
  cl->diErgPerGmUnit = 1./dErgPerGmUnit;
  cl->dKpcUnit = dKpcUnit;
  
  cl->dParam1 = CoolParam.dParam1;
  cl->dParam2 = CoolParam.dParam2;
  cl->dParam3 = CoolParam.dParam3;
  cl->dParam4 = CoolParam.dParam4;
  cl->Y_Total = CoolParam.Y_Total;
  cl->Tmin = CoolParam.dCoolingTmin;
  cl->Tmax = CoolParam.dCoolingTmax;

  /* Derivs Data Struct */
  {
    clDerivsData *Data = cl->DerivsData;

    assert(Data != NULL);

    Data->cl = cl;
    Data->dlnE = (log(EMUL)-log(1/EMUL));
  }
}

#define CL_Rgascode         8.2494e7
#define CL_Eerg_gm_degK     CL_Rgascode
#define CL_ev_degK          1.0/1.1604e4
#define CL_Eerg_gm_ev       CL_Eerg_gm_degK/CL_ev_degK
#define CL_Eerg_gm_degK3_2  1.5*CL_Eerg_gm_degK

/* 
 * Though 13.6eV is lost to the Gas as radiation during H recombination, calculating the
 * Energy using u = E per unit mass = 3/2 n/rho k T requires we don't subtract it there.
 * This formulation is useful because then pressure = 2/3 u rho.
 * Instead we subtract the 13.6eV for H collisional ionization which actually
 * removes no energy from the Gas ( similarly for Helium ) 
 * It also means photoionization doesn't add the 13.6eV, only the excess.
 */
double clThermalEnergy( double Y_Total, double T ) {
  return Y_Total*CL_Eerg_gm_degK3_2*T;

}

double clTemperature( double Y_Total, double E ) {
  return E/(Y_Total*CL_Eerg_gm_degK3_2);
}

double clEdotInstant( COOL *cl, double E, double T, double rho, double rFactor,
		      double *dEdotHeat, double *dEdotCool)
{
	double Edot;

	*dEdotHeat = 0.0;
	*dEdotCool = (T < cl->Tmin ? 0 : E*rFactor);

	return *dEdotHeat - *dEdotCool;
	}

/*
 *  We solve the implicit equation:  
 *  Eout = Ein + dt * Cooling ( Yin, Yout, Ein, Eout )
 *
 *  E erg/gm, PdV erg/gm/s, rho gm/cc, dt sec
 * 
 */

void clDerivs(void *Data, double x, const double *y, double *dHeat,
	      double *dCool) {
  clDerivsData *d = Data;

  d->E = y[0];
  d->T = clTemperature( d->Y_Total, d->E );
  clEdotInstant( d->cl, d->E, d->T, d->rho, d->rFactor, dHeat, dCool );
  if(d->PdV > 0.0)
      *dHeat += d->PdV;
  else
      *dCool -= d->PdV;
}

void clIntegrateEnergy(COOL *cl, double *E, 
		       double PdV, double rho, double Y_Total, double radius, double tStep ) {

  double dEdt,dE,Ein = *E,EMin;
  double t=0,dtused,dtnext,tstop = tStep*(1-1e-8),dtEst;
  clDerivsData *d = cl->DerivsData;
  STIFF *sbs = d->IntegratorContext;
  
  if (tStep <= 0) return;

  d->rho = rho;
  d->PdV = PdV;
  d->Y_Total = Y_Total;
  d->rFactor = cl->dParam1*pow(radius,-3./2.);
  
  EMin = clThermalEnergy( d->Y_Total, cl->Tmin );

  {
    int its = 0;
    StiffStep( sbs, E, t, tStep);
#ifdef ASSERTENEG      
      assert(*E > 0.0);
#else
      if (*E < EMin) {
	*E = EMin;
      }
#endif    
    }
  /* 
     Note Stored abundances are not necessary with
     this eqm integrator therefore the following operations
     could be moved to output routines 
  */

  d->E = *E;
  d->T = clTemperature( d->Y_Total, d->E );
  if (d->T < cl->Tmin ) {
	  d->T = cl->Tmin;
	  *E = clThermalEnergy( d->Y_Total, d->T );
	  }
}

/* Module Interface routines */

void CoolAddParams( COOLPARAM *CoolParam, PRM prm ) {
	CoolParam->dParam1 = 0;
	prmAddParam(prm,"CooldParam1",2,&CoolParam->dParam1,
				sizeof(double),"dparam1",
				"<Param1> = 0.0");
	CoolParam->dParam2 = 0;
	prmAddParam(prm,"CooldParam2",2,&CoolParam->dParam2,
				sizeof(double),"dParam2",
				"<Param2> = 0.0");
	CoolParam->dParam3 = 0;
	prmAddParam(prm,"CooldParam3",2,&CoolParam->dParam3,
				sizeof(double),"dParam3",
				"<Param3> = 0.0");
	CoolParam->dParam4 = 0;
	prmAddParam(prm,"CooldParam4",2,&CoolParam->dParam4,
				sizeof(double),"dParam4",
				"<Param4> = 0.0");
	CoolParam->Y_Total = 2;
	prmAddParam(prm,"dY_Total",2,&CoolParam->Y_Total,
				sizeof(double),"Y_Total",
				"<Y_Total> = 2");
	CoolParam->dCoolingTmin = 10;
	prmAddParam(prm,"dCoolingTmin",2,&CoolParam->dCoolingTmin,
				sizeof(double),"ctmin",
				"<Minimum Temperature for Cooling> = 10K");
	CoolParam->dCoolingTmax = 1e9;
	prmAddParam(prm,"dCoolingTmax",2,&CoolParam->dCoolingTmax,
				sizeof(double),"ctmax",
				"<Maximum Temperature for Cooling> = 1e9K");
	}
	
void CoolLogParams( COOLPARAM *CoolParam, LOGGER *lgr) {
    LogParams(lgr, "COOLING", "CooldParam1: %g",CoolParam->dParam1); 
    LogParams(lgr, "COOLING", "CooldParam2: %g",CoolParam->dParam2); 
    LogParams(lgr, "COOLING", "CooldParam3: %g",CoolParam->dParam3); 
    LogParams(lgr, "COOLING", "CooldParam4: %g",CoolParam->dParam4); 
    LogParams(lgr, "COOLING", "Y_Total: %g",CoolParam->Y_Total); 
    LogParams(lgr, "COOLING", "dCoolingTmin: %g",CoolParam->dCoolingTmin); 
    LogParams(lgr, "COOLING", "dCoolingTmax: %g",CoolParam->dCoolingTmax); 
	}

void CoolOutputArray( COOLPARAM *CoolParam, int cnt, int *type, char *suffix ) {
	*type = OUT_NULL;

	}

/* Initialization Routines */

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix )
{
   int localcntTable = 0;

   *nTableColumns = 0;
   }

void CoolTableRead( COOL *Cool, int nData, void *vData)
{
   fprintf(stderr," Attempt to initialize non-exitent table in cooling\n");
   assert(0);

   }

void CoolDefaultParticleData( COOLPARTICLE *cp )
{
	cp->Y_Total = 2;
	}

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double fMetal )
{
	cp->Y_Total = cl->Y_Total;
	*E = clThermalEnergy(cp->Y_Total,dTemp)*cl->diErgPerGmUnit;
	}

void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ) {
	}

void CoolSetTime( COOL *cl, double dTime, double z ) {
	cl->z = z;
	cl->dTime = dTime;
	}

/* Output Conversion Routines */

double CoolEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E ) {
	return clTemperature( cp->Y_Total, E );
	}

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E,
				    double fMetal) {
	return CoolEnergyToTemperature( Cool, cp, E*Cool->dErgPerGmUnit );
	}

/* Integration Routines */

#define CONVERT_CMPERKPC (3.0857e21)

void CoolIntegrateEnergyCode(COOL *cl, COOLPARTICLE *cp, double *ECode, 
		       double PdVCode, double rhoCode, double ZMetal, double *posCode, double tStep ) {
	double radius;
	
	radius= sqrt(posCode[0]*posCode[0]+posCode[1]*posCode[1]+posCode[2]*posCode[2])
		*cl->dKpcUnit*CONVERT_CMPERKPC;

	*ECode = CoolCodeEnergyToErgPerGm( cl, *ECode );
	clIntegrateEnergy(cl,  ECode, CoolCodeWorkToErgPerGmPerSec( cl, PdVCode ), 
					  CodeDensityToComovingGmPerCc(cl, rhoCode), cp->Y_Total, radius, tStep);
	*ECode = CoolErgPerGmToCodeEnergy(cl, *ECode);

	}

/* Star form function -- do not use */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double T, double dDensity ) {
	assert(0);
	return 0;
	}

/* Not implemented */
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			  double rhoCode, double ZMetal, double *posCode ) {
    double T,E,rho,Edot;

    assert(0);
    return CoolErgPerGmPerSecToCodeWork( cl, Edot );
    }

#endif /* NOCOOLING */
#endif /* GASOLINE */

