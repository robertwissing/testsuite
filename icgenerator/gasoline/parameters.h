
#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

#include <sys/param.h> /* for MAXPATHLEN */
#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif
#ifndef MAXLISTLEN
#define MAXLISTLEN 1024
#endif

#include "cosmo.h"
#include "patch.h"

struct parameters {
	/*
    ** Parameters for PKDGRAV.
    */
	int nThreads;
	int bDiag;
	int bOverwrite;
	int bVWarnings;
	int bVStart;
	int bVStep;
	int bVRungStat;
	int bVDetails;
    int bLogTiming;
    int bLogTimingSubStep;
    int bLogTimingStep;
    int bLogTimingSubStepTot;
    int bLogTimingStepTot;
	int bPeriodic;
    int bInflowOutflow;
	int bRestart;
	int bParaRead;
	int bParaWrite;
	int nIOProcessor;
	int bCannonical;
	int bStandard;
	int bKDK;
	int bBinary;
	int bGravStep;
	int bEpsAccStep;
	int bSqrtPhiStep;
	int bAccelStep; /* true if bEpsAccStep or bSqrtPhiStep */
	int bDensityStep;
    int bDeltaAccelStep;
    int bDeltaAccelStepGasTree;
    int bLongRangeStep;
	int nTruncateRung;
	int bNonSymp;
    int iBinaryOutput;
    int bNoReOrder;
    int bPackedVector;
	int bDoDensity;
    int bInitGasDensity;
 	int iReadIOrder;
 	int bDoIOrderOutput;
 	int bDohOutput;
 	int bDoSphhOutput;
	int bDoPressureOutput;
    int bDoHydroOutput;
	int bDoBalsaraOutput;
	int bDoDivvOutput;
	int bDoCurlvOutput;
	int bDoCSoundOutput;
	int bDodtOutput;
	int bDoIonOutput;
#ifdef COOLING_MOLECULARH
    int bDoCorreL; /* Output the correlation length used for calculating H2 shielding*/
#endif
    int bDoCSound;
#ifdef  RADIATIVEBOX
    int bDoStellarLW; /*Turn this on to have the LW radiation outputted*/
#endif
#ifdef PARTICLESPLIT
    double dInitGasMass;
#endif
	int bDoDminOutput;
	int bSymCool;
    int bDoGravity;
    int bDoSelfGravity;
    int bFandG;
    int bHeliocentric;
    double dSunSoft;
    int bLogHalo;
    double dLogHaloVcirc;
    double dLogHaloEps;
    double dLogHaloFlat;
    double dLogHaloSoft;
    int bHernquistSpheroid;
    int bNFWSpheroid;
    double dNFWm200;
    double dNFWr200;
    double dNFWconc;
    double dNFWsoft;
    int bElliptical;
    int bEllipticalDarkNFW;
    int bHomogSpheroid;
    double dHomogSpheroidM;
    double dHomogSpheroidR;
    int bBodyForce;
	double dBodyForceConst;
	int bGalaxyDiskVerticalPotential;
	double dGalaxyDiskVerticalPotentialVc;
	double dGalaxyDiskVerticalPotentialR;
	double dGalaxyDiskVerticalPotentialStarSigma;
	double dGalaxyDiskVerticalPotentialStarH;
	double dGalaxyDiskVerticalPotentialGasSigma;
	double dGalaxyDiskVerticalPotentialGasH;
	int bMiyamotoDisk;
	int bTimeVarying;
	int bRotatingBar;
    ROTBAR  rotbar;

	int bRotFrame;
	double dOmega;
	double dOmegaDot;
	int bPatch;
    int bPatchVerticalGravity;
    PATCH_PARAMS PP;
	int bSimpleGasDrag;
	int bEpstein;
	double dGamma;
	int nBucket;
	int iOutInterval;
    int iOutMinorInterval;
	int iLogInterval;
	int iCheckInterval;
	int iOrder;
	int bEwald;
	int iEwOrder;
	int nReplicas;
	int iStartStep;
	int iStopStep;
	int nSteps;
	int nSmooth;
  int nSmoothMean;
	int nSmoothMax;
	int iMaxRung;
	int nSuperCool;
	int nGrowMass;
    int bGrowGas;
    int bGrowStar;
    int bGrowDark;
	int iWallRunTime;
	int bPhysicalSoft;
	int bSoftMaxMul;
	int bVariableSoft;
	int nSoftNbr;
	int bSoftByType;
    int bVariableSoftStar;
    int bVariableSoftGas;
    int bVariableSoftDark;
	int bDoSoftOutput;
    int bDoSinks;
    int bBHSink;
    int iOrbitOutInterval;
    int bDoSinksAtStart;
    int bSinkThermal;
    int bSinkForm;
    int bSinkFormJeans;
    int bSinkFormDivV;
    int bSinkFormDivAcc;
    int bSinkFormDV;
    int bSinkFormPotMin;
    int bSinkFormSimple;
    int nSinkFormMin;
    int bSinkMerge;
    int bSinkAngMomOutput;
    double dSinkFormDivVCoeff;
    double dSinkFormDivAccCoeff;
    double dSinkFormDensity;
    double dSinkTimeEligible;
    int iSinkRung;
    int iSinkCurrentRung;
    int nJeans;
	double dEta;
    double dEtaDeltaAccel;
	double dExtraStore;
	double dSoft;
	double dSoftMax;
	double dDelta;
	double dEwCut;
	double dEwhCut;
	double dTheta;
	double dTheta2;
	double daSwitchTheta;
	double dAbsPartial;
	double dRelPartial;
	double dAbsTotal;
	double dRelTotal;
	double dPeriod;
	double dxPeriod;
	double dyPeriod;
	double dzPeriod;
    double dxInflow;
    double dxOutflow;
    double drLastInflow;
    double dxLastInflow;
	CSM csm;
	double dRedTo;
	double dCentMass;
    double dSinkRadius;
    double dSinkBoundOrbitRadius;
    double dSinkMustAccreteRadius;
    double dDeltaSink;
    double dSinkMassMin;
    double dJeansConstant;
    double dBHSinkEddEff;
    double dBHSinkFeedbackEff;
    double dBHSinkEddFactor;
    double dBHSinkFeedbackFactor;
    double dBHSinkAlpha;
    double dSinkCurrentDelta;
    int bBHTurnOffCooling;
    int bSmallBHSmooth;
    int bDoBHKick;
    int bBHForm;
    double dBHFormProb;
    double dInitBHMass;
    int bBHMindv;
    int bBHAccreteAll;
	char achDigitMask[MAXPATHLEN];
	char achInFile[MAXPATHLEN];
	char achOutName[MAXPATHLEN];
	char achDataSubPath[MAXPATHLEN];
    char achOutputListGasRed[MAXLISTLEN];
    char achOutputListDarkRed[MAXLISTLEN];
    char achOutputListStarRed[MAXLISTLEN];
    char achOutputListGasInterval[MAXLISTLEN];
    char achOutputListDarkInterval[MAXLISTLEN];
    char achOutputListStarInterval[MAXLISTLEN];
    char achOutputListGasMinorInterval[MAXLISTLEN];
    char achOutputListDarkMinorInterval[MAXLISTLEN];
    char achOutputListStarMinorInterval[MAXLISTLEN];
	double dCoolFac;
	double dCoolDens;
	double dCoolMaxDens;
    double dGrowMinM;
    double dGrowMaxM;
	double dGrowDeltaM;
	double dGrowStartT;
	double dGrowEndT;
	double dFracNoDomainDecomp;
	double dFracNoDomainDimChoice;
	int bSplitWork;
	int    bRungDD;
	double dRungDDWeight;
	double dMHDAlphaB;
	double dMHDPsidecFac;
	double dMHDCleanFac;
  int dICdensprofile;
  int dICdensdir;
  double dICdensouter;
  double dICdensinner;
  double dICdensR;
  double dICdensRsmooth;
  double dICFmax;
  double dICR0;
  double dICR0Rate;
  double dICRLOOP0;
  double dICQ1Avg;
  double dICrhopow;
	double dalp;
	double dalp2;
	double dalp3;
	double dalp4;
	double dalp5;
	double dalp6;
	double dalp7;
	double dalp8;
	double dalp9;
	double dalp10;
	double dalp11;
	double dalp12;
	double dalp13;
	double dalp14;
	double dalp15;
	double dalp16;
	double dalp17;
	double dalp18;
	double dalp19;
	double dalp20;
	double dalp21;
  double dnloc;
  double dnloc2;
  double dnloc3;

	/*
    ** Additional parameters for GASOLINE.
    */
	int bGeometric;
	int bGasAdiabatic;
	int bGasIsothermal;
	int bGasCooling;
	int bGasCoolingNonEqm;
	int iGasModel;
	double dEtaCourant;
	double dEtauDot;
	double dEtaCourantLong;
    double dTMinDt;
    double duMinDt;
	double duDotLimit;
	double dShockTrackerA;
	double dShockTrackerB;
    double dTauAlpha;
    double dAlphaMax;
    double dAlphaMin;
    double dNAlphaNoise;
    double dAFac;
	double dConstAlpha;
	double dConstBeta;
	double dConstGamma;
	double dMeanMolWeight;
	double dGasConst;
    double dTuFac;
	double dMsolUnit;
	double dKpcUnit;
	double ddHonHLimit;
	double dGmPerCcUnit;
	double dComovingGmPerCcUnit;
	double dErgPerGmUnit;
	double dSecUnit;
	int    bViscosityLimiter;
	int    iViscosityLimiter;
	int    bViscosityLimitdt;
	int    bShockTracker;
    int    bVariableAlpha;
	int    bBulkViscosity;
	int    bGasDomainDecomp;
	int    bLowerSoundSpeed;
	int    bFastGas;
	double dFracFastGas;
	double dhMinOverSoft;
	double dResolveJeans;
    double dMetalDiffusionCoeff;
    double dThermalDiffusionCoeff;
	double dFBMassRatio;
    double dThermalCondCoeff;
    double dThermalCondCoeffCode;
    double dThermalCondSatCoeff;
    double dThermalCond2Coeff;
    double dThermalCond2CoeffCode;
    double dThermalCond2SatCoeff;
    double dEvapCoeff;
    double dEvapMinTemp;
    double dEvapCoeffCode;
    double dEtaDiffusion;
    double dDeltaSph;
    int    bConstantDiffusion;
	int    bDoGas;
	int    bSphStep;
    int    bSphSingleStep;
    int    iRungForceCheck;
#if defined(GASOLINE) && !defined(NOCOOLING)
	COOLPARAM CoolParam;
#endif
#ifdef OUTURBDRIVER
    OUTURBPARAM outurbparam;
#endif
	double dDumpFrameStep;
	double dDumpFrameTime;
	int iTreeZipStep;
    int bTreeZipLocal;
    int    iDirector;
	int    bStarForm;
	int    bFeedBack;
	int    bFormOutputs;
    int    bIonize;
    double dIonizeTime;
    double dIonizeMultiple;
    double dIonizeTMin;
    double dIonizeT;
#ifdef SIMPLESF
	double SSF_dEfficiency;
    double SSF_dTMax;
    double SSF_dPhysDenMin;
    double SSF_dComovingDenMin;
    double SSF_dESNPerStarMass;
    double SSF_dInitStarMass;
	double SSF_dtCoolingShutoff;
    int SSF_bdivv;
#endif
#ifdef TWOPHASE
    double dFBInitialMassLoad;
    double dMultiPhaseMinTemp;
    double dMultiPhaseMaxFrac;
    double dMultiPhaseMaxTime;
#endif
#ifdef STARFORM
	STFM   stfm;
	FB     fb;
    SN sn;
    double dDeltaStarForm;
    int bSNTurnOffCooling;
	int bShortCoolShutoff;
	double dExtraCoolShutoff;
	int bSmallSNSmooth;
    int iStarFormRung;
    int iRandomSeed;
	int nSmoothFeedback;
#endif
    double dHotConvTime;
    double dHotConvTimeMul;
    double dHotConvTimeMin;
    double dHotConvVelMin;
    int bESF;
    double dESFEnergy;
    double dESFTime;
	double dKBoltzUnit;
    double dPext;
    double dvturb;
#ifdef GLASS
	/*
    ** Additional parameters for GLASS.
    */
	double dGlassDamper;
	/* Hack for shock tube glasses */
	double dGlassPoverRhoL;
	double dGlassPoverRhoR;
	double dGlassxL;
	double dGlassxR;
	double dGlassVL;
	double dGlassVR;
#endif
#ifdef COLLISIONS
	/*
    ** Additional parameters for collision code...
    */
	int bFindRejects;
	int iCollLogOption;
    int iMinCollRung;
	char achCollLog[MAXPATHLEN];
	double dSmallStep;
	double dxUnifGrav;
	double dyUnifGrav;
	double dzUnifGrav;
	int    iMinBinaryRung;
    double dBallVelFact;
    double dMaxBinaryEcc;
#ifdef SLIDING_PATCH
    int    iRandStep;
    double dLargeMass;
    double dRandBall;
    int iNextRandomization;
#endif
	COLLISION_PARAMS CP;
#endif /* COLLISIONS */
#ifdef SPECIAL_PARTICLES
	int nSpecial;
	int iSpecialId[MAX_NUM_SPECIAL_PARTICLES];
	SPECIAL_PARTICLE_DATA sSpecialData[MAX_NUM_SPECIAL_PARTICLES];
#endif /* SPECIAL_PARTICLES */
#ifdef RUBBLE_ZML
	int bRubbleStep; /* at the moment, this cannot be changed by user */
#endif
	};

#endif
