#include <math.h>
#include <assert.h>
#include "smoothfcn.h"
#include "linalg.h"
#include <stdlib.h>
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))
#define SECONDSPERYEAR   31557600.

#ifdef COLLISIONS
#include "ssdefs.h"
#include "collision.h"
#endif

#ifdef AGGS
#include "aggs.h"
#endif

/* Benz Method (Default) */
#if !defined(PRES_MONAGHAN) && !defined(PRES_HK)
#define PRES_PDV(a,b) (a)
#define PRES_ACC(a,b) (a+b)
#endif

/* Monaghan Method */
#ifdef PRES_MONAGHAN
#define PRES_PDV(a,b) ((a+b)*0.5)
#define PRES_ACC(a,b) (a+b)
#endif

/* HK */
#ifdef PRES_HK
#define PRES_PDV(a,b) sqrt(a*b)
#define PRES_ACC(a,b) (sqrt(a*b)*2)
#endif
/*
 Change the way the Balsara Switch is applied:
*/
/*
#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch+b->BalsaraSwitch))
#define SWITCHCOMBINE(a,b) (a->BalsaraSwitch > b->BalsaraSwitch ? a->BalsaraSwitch : b->BalsaraSwitch)
#define SWITCHCOMBINE(a,b) (a->BalsaraSwitch*b->BalsaraSwitch)
#define SWITCHCOMBINE(a,b) ((a->BalsaraSwitch*b->BalsaraSwitch > 0.5 || \
           (a->BalsaraSwitch > 0.5 && (dx*a->aPres[0]+dy*a->aPres[1]+dz*a->aPres[2]) > 0) || \
           (b->BalsaraSwitch > 0.5 && (dx*b->aPres[0]+dy*b->aPres[1]+dz*b->aPres[2]) < 0)) ? 1 : 0)
#define SWITCHCOMBINEA(a,b) (a->BalsaraSwitch >= 1 || b->BalsaraSwitch >= 1 ? 1 : 0)
#define SWITCHCOMBINEA(a,b) ((a->BalsaraSwitch*b->BalsaraSwitch)*(a->ShockTracker > b->ShockTracker ? a->ShockTracker : b->ShockTracker))
*/

#define ACCEL(p,j) (((PARTICLE *)(p))->a[j])
#define KPCCM 3.085678e21

#ifdef SHOCKTRACK
/* Shock Tracking on: p->ShockTracker and p->aPres are defined */

#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch+b->BalsaraSwitch))
#define SWITCHCOMBINEA(a,b) ((a->BalsaraSwitch*b->BalsaraSwitch)*a->ShockTracker*b->ShockTracker)
#define SWITCHCOMBINEB(a,b) (a->BalsaraSwitch*b->BalsaraSwitch)

#define ACCEL_PRES(p,j) (((PARTICLE *)(p))->aPres[j])
#define ACCEL_COMB_PRES(p,j) ((((PARTICLE *)(p))->a[j])+=(((PARTICLE *)(p))->aPres[j]))

#else

#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch+b->BalsaraSwitch))
/* New idea -- upwind combine
#define SWITCHCOMBINE(a,b) ((a->fDensity < b->fDensity) ? a->BalsaraSwitch : b->BalsaraSwitch)
*/
#define SWITCHCOMBINEA(a,b) SWITCHCOMBINE(a,b)
#define SWITCHCOMBINEB(a,b) SWITCHCOMBINE(a,b)

#define ACCEL_PRES(p,j) (((PARTICLE *)(p))->a[j])
#define ACCEL_COMB_PRES(p,j)

#endif

#ifdef CULLENDEHNEN
#undef SWITCHCOMBINE
#define SWITCHCOMBINE(a,b) (1.0)
#define SWITCHCOMBINEA(a,b) SWITCHCOMBINE(a,b)
#define SWITCHCOMBINEB(a,b) SWITCHCOMBINE(a,b)
#endif

#define multvs(v,u,s)            /*scalar multiply with vector*/	\
  {                                                                     \
    v[0] = u[0] * s; v[1] = u[1] * s; v[2] = u[2] * s;			\
  }

#define multvsadd(v,u,s)            /*scalar multiply with vector*/	\
    {                                                                       \
    v[0] += u[0] * s; v[1] += u[1] * s; v[2] += u[2] * s;		       \
    }

/* KERNEL is the kernel normalized so that only a division by pi is necessary
   DKERNEL is the kernel derivative divided by (r/h) and normalized so that only a division by pi is necessary */

#define BALL2(a) ((a)->fBall2)
#ifdef QUINTIC
/* Quintic spline. Safe Nsmooths: 108-180    Stable against pairing to NS ~ 190  (Dehnen & Aly 2012) */
#define KERNEL(ak,ar2) 	{                                               \
        double ak1 = sqrt(ar2*0.25),ak2;                                \
        ak2 = (1-ak1); ak=(ak2*ak2)*(ak2*ak2)*ak2;                      \
        if (ak1 < (2./3.)) {                                            \
            ak2 = ((2./3.)-ak1); ak -= 6*(ak2*ak2)*(ak2*ak2)*ak2;       \
            if (ak1 < (1./3.)) {                                        \
                ak2 = ((1./3.)-ak1); ak += 15*(ak2*ak2)*(ak2*ak2)*ak2;  \
                }                                                       \
            }                                                           \
        ak *= (2187./40./8.);                                          \
        }
#define DKERNEL(adk,ar2) {                                              \
        double ak1 = sqrt(ar2*0.25),ak2;                                \
        ak2 = (1-ak1); adk=(ak2*ak2)*(ak2*ak2);                         \
        if (ak1 < (2./3.)) {                                            \
            ak2 = ((2./3.)-ak1); adk -= 6*(ak2*ak2)*(ak2*ak2);          \
            if (ak1 < (1./3.)) {                                        \
                if (ak1 <= 0) {                                         \
                    adk = 0; ak1 = 1;                                   \
                    }                                                   \
                else {                                                  \
                    ak2 = ((1./3.)-ak1); adk += 15*(ak2*ak2)*(ak2*ak2); \
                    }                                                   \
                }                                                       \
            }                                                           \
        adk *= (-5*2187./40./8./4.)/ak1;                                  \
        }
#else
#ifdef QUARTICCORE
/* Quartic core kernel.*/
#define KERNEL(ak,ar2) {                                               \
    double ak1 = sqrt(ar2*0.25), ak2;				       \
  ak2 = (1-ak1); ak=(ak2*ak2)*(ak2*ak2)*ak2;                      \
  if (ak1 < (2./3.)) {                                            \
  ak2 = ((2./3.)-ak1); ak -= 6*(ak2*ak2)*(ak2*ak2)*ak2;       \
  if (ak1 < (1./3.)) {                                        \
  ak2 = ((1./3.)-ak1); ak += 15*(ak2*ak2)*(ak2*ak2)*ak2;  \
  }                                                       \
  }                                                           \
  ak *= (2187./40./8.);                                          \
  }
#define DKERNEL(adk,ar2) {                                              \
  double ak1 = sqrt(ar2*0.25),ak2;                                \
  ak2 = (1-ak1); adk=(ak2*ak2)*(ak2*ak2);                         \
  if (ak1 < (2./3.)) {                                            \
  ak2 = ((2./3.)-ak1); adk -= 6*(ak2*ak2)*(ak2*ak2);          \
  if (ak1 < (1./3.)) {                                        \
  if (ak1 <= 0) {                                         \
  adk = 0; ak1 = 1;                                   \
  }                                                   \
  else {                                                  \
    ak2 = ((1./3.)-ak1); adk += 15*(ak2*ak2)*(ak2*ak2); \
  }                                                   \
  }                                                       \
  }                                                           \
  adk *= (-5*2187./40./8./4.)/ak1;                                  \
  }
#else
#ifdef SINC
#define KERNEL(ak,ar2) {       \
  double nloc = smf->nloc;       \
  if (ar2 <= 0.0){ak = 1.0;}	\
  else if (ar2 >= 4.0){ak = 0.0;}       \
  else {       \
  double au = sqrt(ar2*0.25);       \
  double pi2v = M_PI*au;       \
  ak = pow(((sin(pi2v))/(pi2v)),nloc);       \
  }       \
  }
#define DKERNEL(adk,ar2) {                                             \
  double nloc = smf->nloc;       \
  if (ar2 <= 0){adk = 0.0;}       \
  else if (ar2 >= 4.0){adk = 0.0;}                                    \
  else {       \
  double au = sqrt(ar2*0.25);       \
  double pi2v = M_PI*au;       \
  adk = pow(((sin(pi2v))/(pi2v)),nloc);       \
  adk = nloc*0.25*M_PI*(adk/au)*(1.0/tan(pi2v)-1.0/(pi2v));   \
  }       \
  }
#define KERNEL2(ak,ar2) {       \
  double nloc = smf->nloc2;       \
  if (ar2 <= 0.0){ak = 1.0;}    \
  else if (ar2 >= 4.0){ak = 0.0;}       \
  else {       \
  double au = sqrt(ar2*0.25);       \
  double pi2v = M_PI*au;       \
  ak = pow(((sin(pi2v))/(pi2v)),nloc);       \
  }       \
  }
#define DKERNEL2(adk,ar2) {                                             \
  double nloc = smf->nloc2;       \
  if (ar2 <= 0){adk = 0.0;}       \
  else if (ar2 >= 4.0){adk = 0.0;}                                    \
  else {       \
  double au = sqrt(ar2*0.25);       \
  double pi2v = M_PI*au;       \
  adk = pow(((sin(pi2v))/(pi2v)),nloc);       \
  adk = nloc*0.25*M_PI*(adk/au)*(1.0/tan(pi2v)-1.0/(pi2v));   \
  }       \
  }
#define KERNEL3(ak,ar2) {       \
  double nloc = smf->nloc3;       \
  if (ar2 <= 0.0){ak = 1.0;}    \
  else if (ar2 >= 4.0){ak = 0.0;}       \
  else {       \
  double au = sqrt(ar2*0.25);       \
  double pi2v = M_PI*au;       \
  ak = pow(((sin(pi2v))/(pi2v)),nloc);       \
  }       \
  }
#define DKERNEL3(adk,ar2) {                                             \
  double nloc = smf->nloc3;       \
  if (ar2 <= 0){adk = 0.0;}       \
  else if (ar2 >= 4.0){adk = 0.0;}                                    \
  else {       \
  double au = sqrt(ar2*0.25);       \
  double pi2v = M_PI*au;       \
  adk = pow(((sin(pi2v))/(pi2v)),nloc);       \
  adk = nloc*0.25*M_PI*(adk/au)*(1.0/tan(pi2v)-1.0/(pi2v));   \
  }       \
  }
#else
#ifdef COMBOKERNEL
/* Optimized kernel for N24(D=0) V2 fit to 8 poly */
#define KERNEL(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (2.978328957469771/8./M_1_PI)*smf->Wzero; /* ection */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (2.978328957469771/8./M_1_PI)*(1+au*(au*(-15.978583026452117+au*(68.52718671346985+au*(-138.40625996521663+au*(27.309070778331815+au*(470.88213103138025+au*(-1011.1103515402186+au*(974.3708638849415+au*(-466.9086341578398+90.31457628159393*au)))))))))); \
  }                                                           \
  }
#define DKERNEL(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (2.978328957469771/8./M_1_PI/4.)*((-15.978583026452117*2+au*(68.52718671346985*3+au*(-138.40625996521663*4+au*(27.309070778331815*5+au*(470.88213103138025*6+au*(-1011.1103515402186*7+au*(974.3708638849415*8+au*(-466.9086341578398*9+90.31457628159393*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N32(D=0) V3 fit to 7 poly */
#define KERNEL2(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.0335696191263812/8./M_1_PI)*smf->Wzero2; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.0335696191263812/8./M_1_PI)*(1+au*(au*(-12.782425155701672+au*(56.19201508863019+au*(-182.19749782768355+au*(398.3558098227994+au*(-530.5896242699309+au*(413.15969409011626+au*(-173.66137996508962+au*(30.52340821685955)))))))))); \
  }                                                           \
  }
#define DKERNEL2(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.0335696191263812/8./M_1_PI/4.)*((-12.782425155701672*2+au*(56.19201508863019*3+au*(-182.19749782768355*4+au*(398.3558098227994*5+au*(-530.5896242699309*6+au*(413.15969409011626*7+au*(-173.66137996508962*8+au*(30.52340821685955*9))))))))); \
  }                                                           \
  }
/* Optimized kernel for N64(D=0) V2 fit to 8 poly */
#define KERNEL3(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.059891464822678/8./M_1_PI)*smf->Wzero3; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.059891464822678/8./M_1_PI)*(1+au*(au*(-7.146412238800401+au*(8.486638295713002+au*(-4.600790521116535+au*(22.518325576584346+au*(-44.31003177139933+au*(26.628800503958246+au*(6.682529643873523+au*(-13.233388629830701+3.974329141018585*au)))))))))); \
  }                                                           \
  }
#define DKERNEL3(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.059891464822678/8./M_1_PI/4.)*((-7.146412238800401*2+au*(8.486638295713002*3+au*(-4.600790521116535*4+au*(22.518325576584346*5+au*(-44.31003177139933*6+au*(26.628800503958246*7+au*(6.682529643873523*8+au*(-13.233388629830701*9+3.974329141018585*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N128(D=0) V3 fit to 7 poly */
#define KERNEL4(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (2.9693083398478137/8./M_1_PI)*smf->Wzero4; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (2.9693083398478137/8./M_1_PI)*(1+au*(au*(-7.103803076780569+au*(6.248255511888399+au*(12.443184710092742+au*(-28.509020615239596+au*(29.46272319356744+au*(-22.309861056126483+au*(11.288799908178715+au*(-2.5202785755806913)))))))))); \
  }                                                           \
  }
#define DKERNEL4(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (2.9693083398478137/8./M_1_PI/4.)*((-7.103803076780569*2+au*(6.248255511888399*3+au*(12.443184710092742*4+au*(-28.509020615239596*5+au*(29.46272319356744*6+au*(-22.309861056126483*7+au*(11.288799908178715*8+au*(-2.5202785755806913*9))))))))); \
  }                                                           \
  }
/* Optimized kernel for N256(D=0) V2 fit to 8 poly */
#define KERNEL5(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.0302386742621805/8./M_1_PI)*smf->Wzero5; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.0302386742621805/8./M_1_PI)*(1+au*(au*(-7.577958197648637+au*(10.7477849075376+au*(-12.972795637115478+au*(55.63106045361329+au*(-132.84907985143155+au*(162.5874143983151+au*(-110.90708758987748+au*(40.5188941735908-6.178232656983241*au)))))))))); \
  }                                                           \
  }
#define DKERNEL5(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.0302386742621805/8./M_1_PI/4.)*((-7.577958197648637*2+au*(10.7477849075376*3+au*(-12.972795637115478*4+au*(55.63106045361329*5+au*(-132.84907985143155*6+au*(162.5874143983151*7+au*(-110.90708758987748*8+au*(40.5188941735908*9-6.178232656983241*10*au))))))))); \
  }                                                           \
  }
/* Wendland generalize k=0.4 Kernel fit to 8 level polynomial */
#define KERNEL6(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = (2.0436087619911576/8./M_1_PI)*smf->Wzero6; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (2.0436087619911576/8./M_1_PI)*(1-au*au*(12.673489502128177-au*(56.69654094556688-au*(156.9902567239502-au*(269.3345742836768-au*(270.92524943107117-au*(146.36603519184973-32.808154763943335*au))))))); \
  }                                                           \
  }
#define DKERNEL6(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (-2.0436087619911576/8./M_1_PI/4.)*(25.346979004256355-au*(170.08962283670064-au*(627.9610268958008-au*(1346.6728714183837-au*(1625.551496586427-au*(1024.5622463429481-262.4652381115467*au)))))); \
  }                                                           \
  }
/* Wu l=2 k=1 C2 Kernel */
#define KERNEL7(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = 4.*(15./7./8.)*smf->Wzero7; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak;                                                  \
  ak = ak*ak;                                                  \
  ak =  (15./7./8.)*ak*(4.+16.*au+12.*au*au+3.*au*au*au);   \
  }                                                           \
  }
#define DKERNEL7(adk,ar2) {                              \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                                  \
  adk = adk*adk*adk;                                     \
  adk = -(15./7./8./4.)*7.*adk*(8.+9.*au+3.*au*au); \
  }                                                                       \
  }
/* Buh Kernel */
#define KERNEL8(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.65559/8./M_1_PI)*smf->Wzero8; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = au*au;                                                  \
  ak = ak*ak;                                                 \
  ak = (3.65559/8./M_1_PI)*(1.-au*au*(12.-(32.-21.*au)*au)+12.*ak*log(au)); \
  }                                                           \
  }
#define DKERNEL8(adk,ar2) {                              \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                     \
  adk = au*au;                                     \
  adk = (-3.65559/8./4./M_1_PI)*(24.-(96.-72.*au)*au-48.*adk*log(au)); \
}                                   \
}
/* Wendland generalize k=1/2 Kernel */
#define KERNEL9(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = (3.54624/8./M_1_PI)*smf->Wzero9; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = sqrt(1-au*au);                           \
  ak = (3.54624/8./M_1_PI)*(ak+6.5*au*au*ak+1.5*au*au*(4.+au*au)*(log(au)-log(1.+ak))); \
  }                                                           \
  }
#define DKERNEL9(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else if (ar2 >= 4.0) adk = 0.0;                     \
  else {                              \
  double _a2,au = sqrt(ar2*0.25);                     \
  adk = sqrt(1-au*au);                    \
  _a2 = au*au;                                    \
  adk = (3.54624/8./M_1_PI*3./4./adk)*(6.-6.*_a2+adk*((4.+_a2)*log(au/(1.+adk))+_a2*(log(au)-log(1.+adk)))); \
  }                                                           \
  }
/* Wendland C_2 Kernel */
#define KERNEL10(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (21/16.)*smf->Wzero10; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak;                                                 \
  ak = ak*ak;                                                 \
  ak = (21/16.)*ak*(1+4*au);                                  \
  }                                                           \
  }
#define DKERNEL10(adk,ar2) {                              \
  double au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                     \
  adk = (-21/16.*20./4.)*adk*adk*adk;     \
  }
/* Wendland C_4 Kernel */
#define KERNEL11(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = (495/32./8.)*smf->Wzero11; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak*ak;                                                      \
  ak = ak*ak;                                                 \
  ak = (495/32./8.)*ak*(1+6*au+(35/3.)*au*au); \
  }                                                           \
  }
#define DKERNEL11(adk,ar2) {                                  \
  double _a2,au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                         \
  _a2 = adk*adk;                                      \
  adk = (-495/32.*7./3./4.)*_a2*_a2*adk*(1+5*au);        \
  }
/* Wendland C_6 Kernel */
#define KERNEL12(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = (1365/64./8.)*smf->Wzero12; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak*ak*ak;                                                      \
  ak = ak*ak;                                                 \
  ak = (1365/64./8.)*ak*(1+8*au+(25.)*au*au+(32.)*au*au*au);  \
  }                                                           \
  }
#define DKERNEL12(adk,ar2) {                                  \
  double _a2,au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                         \
  _a2 = adk*adk;                                    \
  adk = (-1365/64./8.*11./2.)*_a2*_a2*_a2*adk*(1+(7.)*au+(16.)*au*au);  \
  }
/* Wendland C_8 Kernel */
#define KERNEL13(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = (1785./64./8.)*smf->Wzero13; /* Dehnen & Aly 2012 correction */  \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak*ak*ak*ak;                                                      \
  ak = ak*ak;                                                 \
  ak = (1785./64./8./5.)*ak*(5.+au*(50.+3.*au*(70.+au*(150.+143.*au))));    \
  }                                                           \
  }
#define DKERNEL13(adk,ar2) {                                  \
  double _a4,au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                         \
  _a4 = adk*adk*adk*adk;                                    \
  adk = (-1785./64./8./5.*26./4.)*_a4*_a4*adk*(5.+3.*au*(15.+au*(53.+77.*au))); \
  }
  /* Wu l=3 k=1 C4 Kernel */
#define KERNEL14(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = 6.*(70./33./8.)*smf->Wzero14; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak*ak;                                                  \
  ak = ak*ak;                                                  \
  ak =  (70./33./8.)*ak*(6.+36.*au+82.*au*au+72.*au*au*au+30.*au*au*au*au+5.*au*au*au*au*au); \
  }                                                           \
  }
#define DKERNEL14(adk,ar2) {                              \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                                  \
  adk = adk*adk*adk*adk*adk;                                     \
  adk = -(70./33./8./4.)*11.*adk*(8.+40.*au+48.*au*au+25.*au*au*au+5*au*au*au*au); \
  }                                                                       \
  }
/* Wendland C_0 Kernel */
#define KERNEL15(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = (15./16.)*smf->Wzero8; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak;                                                 \
  ak = (15./16.)*ak;                                  \
  }                                                           \
  }
#define DKERNEL15(adk,ar2) {                              \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                     \
  adk = (-30./16./4.)*adk/au;     \
  }                                                           \
    }
/* Wu l=3 C0 Kernel */
#define KERNEL16(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = 16.*(6./7./8.)*smf->Wzero6; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak;                                                  \
  ak = ak*ak;                                                  \
  ak = (6./7./8.)*ak*(16.+29.*au+20*au*au+5*au*au*au);             \
  }                                                           \
  }
#define DKERNEL16(adk,ar2) {                              \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                     \
  adk = 1-au*au;                                                  \
  adk = adk*adk*adk;                                     \
  adk = -(6./7./8./4.)*35.*adk/au;                      \
  }                                                                       \
  }
/* Wu l=3 k=2 C2 Kernel */
#define KERNEL17(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = 8.*(10./7./8.)*smf->Wzero2; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak;                                                  \
  ak = ak*ak;                                                  \
  ak =  (10./7./8.)*ak*(1.-au)*(8.+40.*au+48.*au*au+25.*au*au*au+5*au*au*au*au); \
  }                                                           \
  }
#define DKERNEL17(adk,ar2) {                              \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                                  \
  adk = adk*adk*adk*adk;                                     \
  adk = -(10./7./8./4.)*9.*adk*(16.+29.*au+20*au*au+5*au*au*au);    \
  }                                                                       \
  }                                               
/* Wu l=2 k=2 C0 Kernel */
#define KERNEL18(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = 8.*(6./5./8.)*smf->Wzero7; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak*ak;                                                  \
  ak =  (6./5./8.)*ak*(8+9*au+3*au*au); \
  }                                                           \
    }
#define DKERNEL18(adk,ar2) {                              \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                     \
  adk = -(6./5./8./4.)*15.*(au*au-1.)*(au*au-1.)/au;                \
  }                                                                       \
  }
/* Optimized kernel for N32G fit to 8 poly */
#define KERNEL19(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.0702399143211334/8./M_1_PI)*smf->Wzero15; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.0702399143211334/8./M_1_PI)*(1+au*(au*(-15.178969768800927+au*(73.51866968356427+au*(-257.63874599152774+au*(645.4937717480914+au*(-1078.589778178365+au*(1171.7510466338265+au*(-797.951908267712+au*(310.1744769215709-52.578562780647644*au)))))))))); \
  }                                                           \
  }
#define DKERNEL19(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.07023991432113348/8./M_1_PI/4.)*((-15.178969768800927*2+au*(73.51866968356427*3+au*(-257.63874599152774*4+au*(645.4937717480914*5+au*(-1078.589778178365*6+au*(1171.7510466338265*7+au*(-797.951908267712*8+au*(310.1744769215709*9-52.578562780647644*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N32GISPH fit to 8 poly */
#define KERNEL20(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.4003319462332073/8./M_1_PI)*smf->Wzero16; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.4003319462332073/8./M_1_PI)*(1+au*(au*(-37.10514490418083+au*(286.59163136982374+au*(-1109.5005318509457+au*(2409.0689549866306+au*(-3009.6818640918887+au*(2067.880326832294+au*(-618.2115694970241+au*(-36.132229318136154+46.09042647341556*au)))))))))); \
  }                                                           \
  }
#define DKERNEL20(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.4003319462332073/8./M_1_PI/4.)*((-37.10514490418083*2+au*(286.59163136982374*3+au*(-1109.5005318509457*4+au*(2409.0689549866306*5+au*(-3009.6818640918887*6+au*(2067.880326832294*7+au*(-618.2115694970241*8+au*(-36.132229318136154*9+46.09042647341556*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N32GISPH fit POSITIVE FOURIER */
#define KERNEL20x(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.443249702546551/8./M_1_PI)*smf->Wzero16; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.443249702546551/8./M_1_PI)*(1+au*(au*(-38.008666264928586+au*(305.101424878325+au*(-1253.7377710752428+au*(2979.081944200849+au*(-4290.0696725349+au*(3775.016595503547+au*(-1957.271894207935+au*(535.1489534631863-56.260913962902464*au)))))))))); \
  }                                                           \
  }
#define DKERNEL20x(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.443249702546551/8./M_1_PI/4.)*((-38.008666264928586*2+au*(305.101424878325*3+au*(-1253.7377710752428*4+au*(2979.081944200849*5+au*(-4290.0696725349*6+au*(3775.016595503547*7+au*(-1957.271894207935*8+au*(535.1489534631863*9-56.260913962902464*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N64G fit to 8 poly */
#define KERNEL21(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.1160348893350283/8./M_1_PI)*smf->Wzero17; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.1160348893350283/8./M_1_PI)*(1+au*(au*(-10.794751909300937+au*(34.923977443679455+au*(-85.66799299517906+au*(148.65227923861258+au*(-124.78847965696686+au*(-11.062047962891707+au*(104.26365151520089+au*(-73.55050556719418+17.023869894040175*au)))))))))); \
  }                                                           \
  }
#define DKERNEL21(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.1160348893350283/8./M_1_PI/4.)*((-10.794751909300937*2+au*(34.923977443679455*3+au*(-85.66799299517906*4+au*(148.65227923861258*5+au*(-124.78847965696686*6+au*(-11.062047962891707*7+au*(104.26365151520089*8+au*(-73.55050556719418*9+17.023869894040175*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N64GISPH fit to 8 poly */
#define KERNEL22(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.4318494020523547/8./M_1_PI)*smf->Wzero18; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.4318494020523547/8./M_1_PI)*(1+au*(au*(-15.217780383519472+au*(66.96785864668652+au*(-205.6618857567194+au*(448.23132682213446+au*(-637.1401597525037+au*(569.1467883955034+au*(-307.2115693839176+au*(91.32532556233177-11.439904149996446*au)))))))))); \
  }                                                           \
  }
#define DKERNEL22(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.4318494020523547/8./M_1_PI/4.)*((-15.217780383519472*2+au*(66.96785864668652*3+au*(-205.6618857567194*4+au*(448.23132682213446*5+au*(-637.1401597525037*6+au*(569.1467883955034*7+au*(-307.2115693839176*8+au*(91.32532556233177*9-11.439904149996446*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N128G fit to 8 poly */
#define KERNEL23(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.279024374031148/8./M_1_PI)*smf->Wzero19; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.279024374031148/8./M_1_PI)*(1+au*(au*(-9.284045192839276+au*(21.16008630720638+au*(-34.54404508703875+au*(40.35435097155496+au*(24.12654086107369+au*(-149.80272853572703+au*(189.8575599789475+au*(-105.16882473827923+22.30110543510182*au)))))))))); \
  }                                                           \
  }
#define DKERNEL23(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.279024374031148/8./M_1_PI/4.)*((-9.284045192839276*2+au*(21.16008630720638*3+au*(-34.54404508703875*4+au*(40.35435097155496*5+au*(24.12654086107369*6+au*(-149.80272853572703*7+au*(189.8575599789475*8+au*(-105.16882473827923*9+22.30110543510182*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N128GISPH fit to 8 poly */
#define KERNEL24(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.2796080864800374/8./M_1_PI)*smf->Wzero20; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.2796080864800374/8./M_1_PI)*(1+au*(au*(-10.385052498240265+au*(23.481577327747694+au*(-17.882479544931122+au*(-44.44973270243745+au*(195.28848747140484+au*(-338.6360198826279+au*(309.5800010858141+au*(-146.15297293589683+28.15619167916885*au)))))))))); \
  }                                                           \
  }
#define DKERNEL24(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.2796080864800374/8./M_1_PI/4.)*((-10.385052498240265*2+au*(23.481577327747694*3+au*(-17.882479544931122*4+au*(-44.44973270243745*5+au*(195.28848747140484*6+au*(-338.6360198826279*7+au*(309.5800010858141*8+au*(-146.15297293589683*9+28.15619167916885*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N256G fit to 8 poly */
#define KERNEL25(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.3158322817560615/8./M_1_PI)*smf->Wzero21; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.3158322817560615/8./M_1_PI)*(1+au*(au*(-6.94709887927258+au*(0.6116847826437315+au*(40.605173332308+au*(-109.39383705074306+au*(199.68569324945807+au*(-269.18551329586074+au*(230.97653909776236+au*(-108.50616849023058+21.153527253935298*au)))))))))); \
  }                                                           \
  }
#define DKERNEL25(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.3158322817560615/8./M_1_PI/4.)*((-6.94709887927258*2+au*(0.6116847826437315*3+au*(40.605173332308*4+au*(-109.39383705074306*5+au*(199.68569324945807*6+au*(-269.18551329586074*7+au*(230.97653909776236*8+au*(-108.50616849023058*9+21.153527253935298*10*au))))))))); \
  }                                                           \
  }
/* Optimized kernel for N256G fit to 8 poly */
#define KERNEL25X(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.317969661358544/8./M_1_PI)*smf->Wzero21; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.317969661358544/8./M_1_PI)*(1+au*(au*(-7.286240304649148+au*(5.818434743290395+au*(7.261024068735384+au*(6.504045269431593+au*(-40.310178334474344+au*(35.64179808739792+au*(-2.192539059615132+au*(-9.825094330475823+3.3887498603591095*au)))))))))); \
  }                                                           \
  }
#define DKERNEL25X(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.317969661358544/8./M_1_PI/4.)*((-7.286240304649148*2+au*(5.818434743290395*3+au*(7.261024068735384*4+au*(6.504045269431593*5+au*(-40.310178334474344*6+au*(35.64179808739792*7+au*(-2.192539059615132*8+au*(-9.825094330475823*9+3.3887498603591095*10*au))))))))); \
  }                                                           \
  }

/* Optimized kernel for N256GISPH fit to 8 poly */
#define KERNEL26(ak,ar2)  {                                               \
  if (ar2 <= 0) ak = (3.181335291127915/8./M_1_PI)*smf->Wzero22; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = (3.181335291127915/8./M_1_PI)*(1+au*(au*(-10.975827966498654+au*(34.022346938397256+au*(-79.62182205737668+au*(145.2028297298656+au*(-156.0395550888878+au*(69.57842722118721+au*(17.65497978738222+au*(-28.52004202463561+7.698663460565884*au)))))))))); \
  }                                                           \
  }
#define DKERNEL26(adk,ar2) {                                  \
  if (ar2 <= 0) adk = 0.0; /* Dehnen & Aly 2012 correction */ \
  else {                                                      \
  double au = sqrt(ar2*0.25);                     \
  adk = (3.181335291127915/8./M_1_PI/4.)*((-10.975827966498654*2+au*(34.022346938397256*3+au*(-79.62182205737668*4+au*(145.2028297298656*5+au*(-156.0395550888878*6+au*(69.57842722118721*7+au*(17.65497978738222*8+au*(-28.52004202463561*9+7.698663460565884*10*au))))))))); \
  }                                                           \
  }
#else
#ifdef WENDLAND
/* Wendland C_2 Kernel */
#ifdef SPH1D
/* Need a new normalization for this Kernel for 1D */
#define KERNEL(ak,ar2) 	{						\
        double au = sqrt(ar2*0.25);					\
	    ak = 1-au;							\
	    ak = ak*ak;							\
	    ak = ak*ak;							\
	    ak = (21/16.)*ak*(1+4*au);					\
	}
#else
#define KERNEL(ak,ar2) 	{						\
    if (ar2 <= 0) ak = (21/16.)*smf->Wzero; /* Dehnen & Aly 2012 correction (1-0.0454684 at Ns=64) */ \
	else {								\
	    double au = sqrt(ar2*0.25);					\
	    ak = 1-au;							\
	    ak = ak*ak;							\
	    ak = ak*ak;							\
	    ak = (21/16.)*ak*(1+4*au);					\
	    }								\
	}
#endif
#define DKERNEL(adk,ar2) {				\
	double au = sqrt(ar2*0.25);			\
	adk = 1-au;					\
	adk = (-21/16.*20./4.)*adk*adk*adk;	\
	}
#else
#ifdef WENDLANDC4
/* Wendland C_4 Kernel */
#define KERNEL(ak,ar2) 	{						\
    if (ar2 <= 0) ak = (495/32./8.)*smf->Wzero; /* Dehnen & Aly 2012 correction */ \
	else {								\
	    double au = sqrt(ar2*0.25);					\
	    ak = 1-au;							\
	    ak = ak*ak*ak;							\
	    ak = ak*ak;							\
	    ak = (495/32./8.)*ak*(1+6*au+(35/3.)*au*au); \
	    }								\
	}
#define DKERNEL(adk,ar2) {                                  \
        double _a2,au = sqrt(ar2*0.25);                     \
        adk = 1-au;                                         \
        _a2 = adk*adk;                                      \
        adk = (-495/32.*7./3./4.)*_a2*_a2*adk*(1+5*au);        \
	}
#else
#ifdef WENDLANDC6
/* Wendland C_6 Kernel */
#define KERNEL(ak,ar2)  {                                               \
    if (ar2 <= 0) ak = (1365/64./8.)*smf->Wzero; /* Dehnen & Aly 2012 correction */ \
  else {                                                          \
  double au = sqrt(ar2*0.25);                                 \
  ak = 1-au;                                                  \
  ak = ak*ak*ak*ak;                                                      \
  ak = ak*ak;                                                 \
  ak = (1365/64./8.)*ak*(1+8*au+(25.)*au*au+(32.)*au*au*au);  \
  }                                                           \
  }
#define DKERNEL(adk,ar2) {                                  \
  double _a2,au = sqrt(ar2*0.25);                     \
  adk = 1-au;                                         \
  _a2 = adk*adk;                                    \
  adk = (-1365/64./8.*11./2.)*_a2*_a2*_a2*adk*(1+(7.)*au+(16.)*au*au);	\
  }
#else
#ifdef PEAKEDKERNEL
/* Standard M_4 Kernel with central peak for dW/du according to Thomas and Couchman 92 (Steinmetz 96) */
#define KERNEL(ak,ar2) { \
		ak = 2.0 - sqrt(ar2); \
		if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
		else ak = 0.25*ak*ak*ak; \
        }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
            if (adk < (2./3.)) { \
               if (adk > 0) adk = -1/adk; \
			   } \
            else {		    \
               adk = -3 + 2.25*adk; \
			   } \
			} \
		else { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
		}
#else
#ifdef M43D
/* M43D Creates a 3D kernel by convolution of 3D tophats the way M4(1D) is made in 1D */
#define KERNEL(ak,ar2) { \
		ak = sqrt(ar2); \
		if (ar2 < 1.0) ak = 6.*0.25/350./3. *(1360+ar2*(-2880 \
			 +ar2*(3528+ak*(-1890+ak*(-240+ak*(270-6*ar2)))))); \
		else ak = 6.*0.25/350./3. *(7040-1152/ak+ak*(-10080+ak*(2880+ak*(4200 \
	                 +ak*(-3528+ak*(630+ak*(240+ak*(-90+2*ar2)))))))); \
                }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) adk = 6.*0.25/350./3. * (-2880*2 \
	                 +ar2*(3528*4+ adk*(-1890*5 + adk*(-240*6+ adk*(270*7-6*9*ar2))))); \
		else adk = 6.*0.25/350./3. *((1152/ar2-10080)/adk+(2880*2+adk*(4200*3 \
	                 +adk*(-3528*4+adk*(630*5+adk*(240*6 +adk*(-90*7+2*9*ar2))))))); \
                }

#else
#ifdef HSHRINK
/* HSHRINK M4 Kernel uses an effective h of (pi/6)^(1/3) times h for nSmooth neighbours */
#define dSHRINKFACTOR        0.805995977
#undef BALL2
#define BALL2(a) ((a)->fBall2*(dSHRINKFACTOR*dSHRINKFACTOR))
#define KERNEL(ak,ar2) { \
		ak = 2.0 - sqrt(ar2); \
		if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
		else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
		else ak = 0.0; \
                }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
			adk = -3 + 2.25*adk; \
			} \
		else if (ar2 < 4.0) { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
		else adk = 0.0; \
                }

#else
/* Standard M_4 Kernel */
#define KERNEL(ak,ar2) { \
		ak = 2.0 - sqrt(ar2); \
		if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
		else ak = 0.25*ak*ak*ak; \
                }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
			adk = -3 + 2.25*adk; \
			} \
		else { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
                }
#endif
#endif
#endif
#endif /* WENDLANDC6 */
#endif /* WENDLANDC4 */
#endif /* WENDLAND */
#endif /* QUINTIC */
#endif /*QUARTIC */
#endif /*SINC*/
#endif

void NullSmooth(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
}

void initDensity(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
#ifdef DENSITYUOLD
	((PARTICLE *)p)->fDensityU = 0.0;
#endif
	}

void combDensity(void *p1,void *p2)
{
	((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
#ifdef DENSITYUOLD
	((PARTICLE *)p1)->fDensityU += ((PARTICLE *)p2)->fDensityU;
#endif
	}

void Density(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    FLOAT ih2,r2,rs,fDensity,fNorm;
#ifdef DENSITYUOLD
	FLOAT fDensityU = 0;
#endif
	int i;

	ih2 = 4.0/BALL2(p);
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
#ifdef FBPARTICLE
        assert(TYPETest(nnList[i].pPart,TYPE_GAS));
#endif
		KERNEL(rs,r2);
		fDensity += rs*nnList[i].pPart->fMass;
#ifdef DENSITYUOLD
		p->fDensityU += rs*nnList[i].pPart->fMass*(nnList[i].pPart->uPred+smf->uMin);
#endif
		}
#ifdef SPH1D
    fNorm = (2/3.)*sqrt(ih2);
#else
    fNorm = M_1_PI*sqrt(ih2)*ih2;
#endif
	p->fDensity = fNorm*fDensity;
#ifdef DENSITYUOLD
	p->fDensityU = fNorm*fDensityU;
#endif
	}

void DensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i;

	ih2 = 4.0/(BALL2(p));
#ifdef SPH1D
	fNorm = (1./3.)*sqrt(ih2);
#else
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
#endif
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rs *= fNorm;
		q = nnList[i].pPart;
		p->fDensity += rs*q->fMass;
		q->fDensity += rs*p->fMass;
#ifdef DENSITYUOLD
		p->fDensityU += rs*q->fMass*(q->uPred+smf->uMin);
		q->fDensityU += rs*p->fMass*(p->uPred+smf->uMin);
#endif
		}
	}

void initDensityTmp(void *p)
{
#ifdef GASOLINE
	((PARTICLE *)p)->curlv[0] = 0.0;
#else
	((PARTICLE *)p)->fDensity = 0.0;
#endif
	}

void combDensityTmp(void *p1,void *p2)
{
#ifdef GASOLINE
	((PARTICLE *)p1)->curlv[0] += ((PARTICLE *)p2)->curlv[0];
#else
	((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
#endif
	}

void DensityTmp(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    FLOAT ih2,r2,rs,fDensity,fNorm;
	int i;

	ih2 = 4.0/BALL2(p);
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		fDensity += rs*nnList[i].pPart->fMass;
		}

#ifdef SPH1D
	fNorm = (2./3.)*sqrt(ih2);
#else
	fNorm = M_1_PI*sqrt(ih2)*ih2;
#endif
#ifdef GASOLINE
	p->curlv[0] = fNorm*fDensity;
#else
	p->fDensity = fNorm*fDensity;
#endif
	}

void DensityTmpSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i;

	ih2 = 4.0/(BALL2(p));
#ifdef SPH1D
	fNorm = (1./3.)*sqrt(ih2);
#else
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
#endif
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rs *= fNorm;
		q = nnList[i].pPart;
#ifdef GASOLINE
		p->curlv[0] += rs*q->fMass;
		q->curlv[0] += rs*p->fMass;
#else
		p->fDensity += rs*q->fMass;
		q->fDensity += rs*p->fMass;
#endif
		}
	}

/* Mark functions:
   These functions calculate density and mark the particles in some way
   at the same time.
   Only SmoothActive particles do this gather operation (as always).

   MarkDensity:
	 Combination: if porig DensZeroed combine porig+pcopy into porig
                  if pcopy DensZeroed set porig = pcopy
     Smooth:
       Active particles get density, nbrofactive
	   All gather neighbours are labelled as DensZeroed, get density
       --> effectively all particles and gather neighbours get density and are labelled DensZeroed
       --> These densities will potentially be lacking scatter neighbours so only correct
           if all particles involved in this operation OR scatter later added
       Gather/Scatter Neighbours of Active Particles get nbrofactive

   MarkSmooth:
     Go through full tree looking for particles than touch a smooth active particle
     and mark them with specified label: eg. TYPE_Scatter

   MarkIIDensity:
     Init:        Densactive particles not dens zeroed get dens zeroed
	 Combination: If pcopy is active make porig active
	              Densactive particles only:
                        if porig DensZeroed combine density porig+pcopy into porig
                        if pcopy DensZeroed set density  porig = pcopy
     Smooth:
       Densactive: get density
                   Densactive gather neighbours get density, DensZeroed (Nbrofactive if reqd)
	   Not Densactive, but Active:
                   get nbrofactive
			       Densactive gather neighbours get density, DensZeroed	(Nbrofactive)
       Not Densactive, Not Active
                   get nbrofactive if gather neighbour active
			       Densactive gather neighbours get density, DensZeroed
*/

void initParticleMarkDensity(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
#ifdef DENSITYUOLD
	((PARTICLE *)p)->fDensityU = 0.0;
#endif
	TYPESet((PARTICLE *) p,TYPE_DensZeroed);
	}

void initMarkDensity(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
#ifdef DENSITYUOLD
	((PARTICLE *)p)->fDensityU = 0.0;
#endif
	}

void combMarkDensity(void *p1,void *p2)
{
        if (TYPETest((PARTICLE *) p1,TYPE_DensZeroed)) {
		((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
#ifdef DENSITYUOLD
		((PARTICLE *)p1)->fDensityU += ((PARTICLE *)p2)->fDensityU;
#endif
        }
	else if (TYPETest((PARTICLE *) p2,TYPE_DensZeroed)) {
		((PARTICLE *)p1)->fDensity = ((PARTICLE *)p2)->fDensity;
#ifdef DENSITYUOLD
		((PARTICLE *)p1)->fDensityU = ((PARTICLE *)p2)->fDensityU;
#endif
		}
    ((PARTICLE *)p1)->iActive |= (((PARTICLE *)p2)->iActive & TYPE_MASK);
	}

void MarkDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	assert(0);
}

void MarkDensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i;
	unsigned int qiActive;

	assert(TYPETest(p,TYPE_GAS));
	ih2 = 4.0/(BALL2(p));
#ifdef SPH1D
	fNorm = (1./3.)*sqrt(ih2);
#else
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
#endif
	if (TYPETest(p,TYPE_ACTIVE)) {
		TYPESet(p,TYPE_NbrOfACTIVE);
		for (i=0;i<nSmooth;++i) {
			r2 = nnList[i].fDist2*ih2;
			KERNEL(rs,r2);
			rs *= fNorm;
			q = nnList[i].pPart;
			assert(TYPETest(q,TYPE_GAS));
			p->fDensity += rs*q->fMass;
#ifdef DENSITYUOLD
			p->fDensityU += rs*q->fMass*(q->uPred+smf->uMin);
#endif
			if (TYPETest(q,TYPE_DensZeroed)) {
				q->fDensity += rs*p->fMass;
#ifdef DENSITYUOLD
				q->fDensityU += rs*p->fMass*(p->uPred+smf->uMin);
#endif
			    }
			else {
				q->fDensity = rs*p->fMass;
#ifdef DENSITYUOLD
				q->fDensityU = rs*p->fMass*(p->uPred+smf->uMin);
#endif
				TYPESet(q, TYPE_DensZeroed);
				}
			TYPESet(q,TYPE_NbrOfACTIVE);
			}
		}
	else {
		qiActive = 0;
		for (i=0;i<nSmooth;++i) {
			r2 = nnList[i].fDist2*ih2;
			KERNEL(rs,r2);
			rs *= fNorm;
			q = nnList[i].pPart;
			assert(TYPETest(q,TYPE_GAS));
			if (TYPETest(p,TYPE_DensZeroed)) {
				p->fDensity += rs*q->fMass;
#ifdef DENSITYUOLD
				p->fDensityU += rs*q->fMass*(q->uPred+smf->uMin);
#endif
			    }
			else {
				p->fDensity = rs*q->fMass;
#ifdef DENSITYUOLD
				p->fDensityU = rs*q->fMass*(q->uPred+smf->uMin);
#endif
				TYPESet(p,TYPE_DensZeroed);
				}
			if (TYPETest(q,TYPE_DensZeroed)) {
				q->fDensity += rs*p->fMass;
#ifdef DENSITYUOLD
				q->fDensityU += rs*p->fMass*(p->uPred+smf->uMin);
#endif
			    }
			else {
				q->fDensity = rs*p->fMass;
#ifdef DENSITYUOLD
				q->fDensityU = rs*p->fMass*(p->uPred+smf->uMin);
#endif
				TYPESet(q, TYPE_DensZeroed);
				}
			qiActive |= q->iActive;
			}
		if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
		}
	}

void initParticleMarkIIDensity(void *p)
{
	if (TYPEFilter((PARTICLE *) p,TYPE_DensACTIVE|TYPE_DensZeroed,
				   TYPE_DensACTIVE)) {
		((PARTICLE *)p)->fDensity = 0.0;
#ifdef DENSITYUOLD
		((PARTICLE *)p)->fDensityU = 0.0;
#endif
		TYPESet((PARTICLE *)p,TYPE_DensZeroed);
/*		if (((PARTICLE *)p)->iOrder == CHECKSOFT) fprintf(stderr,"Init Zero Particle 3031A: %g \n",((PARTICLE *) p)->fDensity);*/
		}
	}
/* copies of remote particles */
void initMarkIIDensity(void *p)
    {
    ((PARTICLE *) p)->fDensity = 0.0;
#ifdef DENSITYUOLD
    ((PARTICLE *)p)->fDensityU = 0.0;
#endif
/*    if (((PARTICLE *)p)->iOrder == CHECKSOFT) fprintf(stderr,"Init Cache Zero Particle 3031A: %g \n",((PARTICLE *) p)->fDensity);*/
    }

void combMarkIIDensity(void *p1,void *p2)
    {
    if (TYPETest((PARTICLE *) p1,TYPE_DensACTIVE)) {
	if (TYPETest((PARTICLE *) p1,TYPE_DensZeroed)) {
	    ((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
#ifdef DENSITYUOLD
	    ((PARTICLE *)p1)->fDensityU += ((PARTICLE *)p2)->fDensityU;
#endif
	    }
	else if (TYPETest((PARTICLE *) p2,TYPE_DensZeroed)) {
	    ((PARTICLE *)p1)->fDensity = ((PARTICLE *)p2)->fDensity;
#ifdef DENSITYUOLD
	    ((PARTICLE *)p1)->fDensityU = ((PARTICLE *)p2)->fDensityU;
#endif
	    TYPESet((PARTICLE *)p1,TYPE_DensZeroed);
	    }
	}
    ((PARTICLE *)p1)->iActive |= (((PARTICLE *)p2)->iActive & TYPE_MASK);
    }

void MarkIIDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    assert(0);
    }

void MarkIIDensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    PARTICLE *q;
    FLOAT fNorm,ih2,r2,rs;
    int i;
    unsigned int qiActive;

    ih2 = 4.0/(BALL2(p));
#ifdef SPH1D
	fNorm = (1./3.)*sqrt(ih2);
#else
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
#endif
    if (TYPETest(p,TYPE_DensACTIVE)) {
	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    qiActive |= q->iActive;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    p->fDensity += rs*q->fMass;
#ifdef DENSITYUOLD
	    p->fDensityU += rs*q->fMass*(q->uPred+smf->uMin);
#endif
/*	    if (p->iOrder == CHECKSOFT) fprintf(stderr,"DensActive Particle %iA: %g %i  %g\n",p->iOrder,p->fDensity,q->iOrder,q->fMass);*/
	    if (TYPETest(q,TYPE_DensACTIVE)) {
		if (TYPETest(q,TYPE_DensZeroed)) {
		    q->fDensity += rs*p->fMass;
#ifdef DENSITYUOLD
		    q->fDensityU += rs*p->fMass*(p->uPred+smf->uMin);
#endif
/*		    if (q->iOrder == CHECKSOFT) fprintf(stderr,"qDensActive Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		    }
		else {
		    q->fDensity = rs*p->fMass;
#ifdef DENSITYUOLD
		    q->fDensityU = rs*p->fMass*(p->uPred+smf->uMin);
#endif
		    TYPESet(q,TYPE_DensZeroed);
/*		    if (q->iOrder == CHECKSOFT) fprintf(stderr,"zero qDensActive Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		    }
		}
	    if (TYPETest(p,TYPE_ACTIVE)) TYPESet(q,TYPE_NbrOfACTIVE);
	    }
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
	}
    else if (TYPETest(p,TYPE_ACTIVE)) {
	TYPESet( p,TYPE_NbrOfACTIVE);
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    TYPESet(q,TYPE_NbrOfACTIVE);
	    if (!TYPETest(q,TYPE_DensACTIVE)) continue;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    if (TYPETest(q,TYPE_DensZeroed)) {
		q->fDensity += rs*p->fMass;
#ifdef DENSITYUOLD
		q->fDensityU += rs*p->fMass*(p->uPred+smf->uMin);
#endif
/*		if (q->iOrder == CHECKSOFT) fprintf(stderr,"qActive Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		}
	    else {
		q->fDensity = rs*p->fMass;
#ifdef DENSITYUOLD
		q->fDensityU = rs*p->fMass*(p->uPred+smf->uMin);
#endif
		TYPESet(q,TYPE_DensZeroed);
/*		if (q->iOrder == CHECKSOFT) fprintf(stderr,"zero qActive Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		}
	    }
	}
    else {
	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    qiActive |= q->iActive;
	    if (!TYPETest(q,TYPE_DensACTIVE)) continue;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    if (TYPETest(q,TYPE_DensZeroed)) {
		q->fDensity += rs*p->fMass;
#ifdef DENSITYUOLD
		q->fDensityU += rs*p->fMass*(p->uPred+smf->uMin);
#endif
/*		if (q->iOrder == CHECKSOFT) fprintf(stderr,"qOther Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		}
	    else {
		q->fDensity = rs*p->fMass;
#ifdef DENSITYUOLD
		q->fDensityU = rs*p->fMass*(p->uPred+smf->uMin);
#endif
		TYPESet(q,TYPE_DensZeroed);
/*		if (q->iOrder == CHECKSOFT) fprintf(stderr,"zero qOther Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		}
	    }
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
	}
    }


void initMark(void *p)
{
        }

void combMark(void *p1,void *p2)
{
    ((PARTICLE *)p1)->iActive |= (((PARTICLE *)p2)->iActive & TYPE_MASK);
	}


void initDeltaAccel(void *p)
{
	}

void combDeltaAccel(void *p1,void *p2)
{
    if (TYPEQueryACTIVE((PARTICLE *) p1) && ((PARTICLE *)p2)->dt < ((PARTICLE *)p1)->dt)
	    ((PARTICLE *)p1)->dt = ((PARTICLE *)p2)->dt;
    }

void DeltaAccel(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	int i;
    FLOAT dax,da2,r2,dt;
	PARTICLE *q;

#ifdef DELTACCELCAP
	FLOAT pSoft2,qSoft2;
	pSoft2 = p->fSoft*p->fSoft;
#endif
	/*	assert(TYPEQueryACTIVE((PARTICLE *) p)); */

	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2;
		if (r2 > 0) {
		  q = nnList[i].pPart;
#ifdef DELTAACCELACTIVE
		  if (!TYPEQueryACTIVE((PARTICLE *) q)) continue;
#endif
		  dax = p->a[0]-q->a[0];
		  da2 = dax*dax;
		  dax = p->a[1]-q->a[1];
		  da2 += dax*dax;
		  dax = p->a[2]-q->a[2];
		  da2 += dax*dax;
		  if (da2 > 0) {
#ifdef DELTACCELCAP
			if (r2 < pSoft2) r2 = pSoft2;
			qSoft2 = q->fSoft*q->fSoft;
			if (r2 < qSoft2) r2 = qSoft2;
#endif

   		    dt = smf->dDeltaAccelFac*sqrt(sqrt(r2/da2));  /* Timestep dt = Eta sqrt(deltar/deltaa) */
		    if (dt < p->dt) p->dt = dt;
		    if (
#ifndef DELTAACCELACTIVE
				TYPEQueryACTIVE((PARTICLE *) q) &&
#endif
				(dt < q->dt)) q->dt = dt;
		    }
		  }
	    }
    }

#define fBindingEnergy(_a)  (((PARTICLE *) (_a))->curlv[0])
#define iOrderSink(_a)      (*((int *) (&((PARTICLE *) (_a))->curlv[1])))
/* Indicator for r,v,a update */
#define bRVAUpdate(_a)         (((PARTICLE *) (_a))->curlv[2])

void initSinkAccreteTest(void *p)
{
#ifdef GASOLINE
    fBindingEnergy(p) = FLT_MAX;
    iOrderSink(p) = -1;
#ifdef SINKING
    if (TYPETest((PARTICLE *) p, TYPE_SINKING)) {
	iOrderSink(p) = ((PARTICLE *) p)->iSinkingOnto;
	}
#endif
#endif
}

void combSinkAccreteTest(void *p1,void *p2)
{
#ifdef GASOLINE
/* Particle p1 belongs to sink iOrderSink(p1) initially but will
   switch to iOrderSink(p2) if more bound to that sink */
#ifdef SINKING
    if (TYPETest((PARTICLE *) p1, TYPE_SINKING)) return;
#endif
    if (fBindingEnergy(p2) < fBindingEnergy(p1)) {
	fBindingEnergy(p1) = fBindingEnergy(p2);
	iOrderSink(p1) = iOrderSink(p2);
	}

    if (TYPETest( ((PARTICLE *) p2), TYPE_NEWSINKING)) TYPESet( ((PARTICLE *) p1), TYPE_NEWSINKING);
#endif
}

#ifdef SINKING
#define TRUEMASS(p__) (p__)->fTrueMass
#else
#define TRUEMASS(p__) (p__)->fMass
#endif

void SinkAccreteTest(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
	int i;
	double dSinkRadius2 = smf->dSinkRadius*smf->dSinkRadius,
	       EBO,Eq,r2,dvx,dv2;
	PARTICLE *q;

#ifdef SINKING
	assert(TYPETest(p, TYPE_SINK));
	assert(!TYPETest(p, TYPE_SINKING));
#endif
	/* G = 1
	 p is sink particle
	 q is gas particle */
        if (smf->dSinkBoundOrbitRadius > 0)
            EBO = -0.5*TRUEMASS(p)/smf->dSinkBoundOrbitRadius;
        else
            EBO = FLT_MAX;

	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2;
	    if (r2 > 0 && r2 <= dSinkRadius2) {
		q = nnList[i].pPart;
#ifdef SINKING
/*		    if (TYPETest(q, TYPE_SINKING)) continue;*/
		if (TYPETest(q, TYPE_SINKING)) {
		    FLOAT r0 = q->rSinking0Mag;
		    FLOAT r2 = r0 + q->vSinkingr0*(smf->dTime-q->fSinkingTime);
/*
  if (p->iOrder == 2000000 && q->iSinkingOnto == 2000000) printf("SINKACCRETETEST %d eats sinking %d %g %g %g\n",p->iOrder,q->iOrder,sqrt((p->r[0]-q->r[0])*(p->r[0]-q->r[0])+(p->r[1]-q->r[1])*(p->r[1]-q->r[1])+(p->r[2]-q->r[2])*(p->r[2]-q->r[2])), r0,r2 );*/
		    continue;
		    }
#endif
		if (TYPETest( q, TYPE_GAS ) && (q->iRung >= smf->iSinkCurrentRung || TYPETest(q, TYPE_NEWSINKING))) {
		    dvx = p->v[0]-q->v[0];
		    dv2 = dvx*dvx;
		    dvx = p->v[1]-q->v[1];
		    dv2 += dvx*dvx;
		    dvx = p->v[2]-q->v[2];
		    dv2 += dvx*dvx;
		    Eq = -TRUEMASS(p)/sqrt(r2) + 0.5*dv2;
		    if (smf->bSinkThermal) Eq+= q->u;
		    /* Being labelled NEWSINKING forces accretion to somebody even if unbound */
		    if (Eq < EBO || r2 < smf->dSinkMustAccreteRadius*smf->dSinkMustAccreteRadius || TYPETest(q, TYPE_NEWSINKING)) {
			if (Eq < fBindingEnergy(q)) {
			    fBindingEnergy(q) = Eq;
			    iOrderSink(q) = p->iOrder; /* Particle q belongs to sink p */
			    TYPESet(q, TYPE_NEWSINKING);
			    }
			}
		    }
		}
	    }
#endif
}

void SinkingAverage(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
#ifdef SINKING
	int i,j;
	PARTICLE *q;
	double dv[3],wt,ih2,rs,r2,norm;

	assert(TYPETest(p,TYPE_NEWSINKING));
	wt = 0;
	for (j=0;j<3;j++) dv[j]=0;

	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    ih2 = 4.0/BALL2(p);
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs*= q->fMass;
	    wt+=rs;
	    for (j=0;j<3;j++) dv[j]+=rs*(q->v[j]-p->v[j]);
	    }

	norm=1./wt;
        /* smoothed velocity */
	for (j=0;j<3;j++) {
	    p->vSinkingTang0Unit[j] = dv[j]*norm;
	    }
#endif
#endif
}


void SinkLogEntry( SINKLOG *pSinkLog, SINKEVENT *pSinkEvent, int iSinkEventType) {

    if (pSinkEvent == NULL) return;

    switch (iSinkEventType) {
    case SINK_EVENT_NULL:
	return;
    case SINK_EVENT_ACCRETE_AT_FORMATION:
	pSinkEvent->iOrdSink = -1; /* Truncated Event -- report Victim iOrder */
	assert(pSinkEvent->iOrdVictim > 0);
	pSinkLog->nAccrete++;
	break;
    case SINK_EVENT_FORM:
	pSinkEvent->iOrdSink = pSinkLog->nForm; /* Actual Sink Formation Event -- Victims listed previously */
	pSinkEvent->iOrdVictim = -1;
	pSinkLog->nForm++;
	break;
    case SINK_EVENT_ACCRETE:
	pSinkEvent->iOrdSink = -2; /* Truncated Event -- report Victim iOrder */
	pSinkLog->nAccrete++;
	break;
    case SINK_EVENT_ACCRETE_UPDATE:
	assert(pSinkEvent->iOrdSink > 0); /* Detailed Event -- report Victim iOrder and Sink Data */
	pSinkLog->nAccrete++;
	break;
    case SINK_EVENT_MERGER: /* Detailed Event -- report deleted and full Sink Data */
	break;
    default:
	assert(0);
	}

    if(pSinkLog->nLog >= pSinkLog->nMaxLog) {
	/* Grow table */
	pSinkLog->nMaxLog *= 1.4;
	pSinkLog->SinkEventTab = realloc(pSinkLog->SinkEventTab,pSinkLog->nMaxLog*sizeof(SINKEVENT));
	assert(pSinkLog->SinkEventTab != NULL);
	}

    pSinkLog->SinkEventTab[pSinkLog->nLog] = *pSinkEvent;
    pSinkLog->nLog++;
    }


void initSinkAccrete(void *p)
    { /* cached copies only */
#ifdef SINKING
    if (TYPETest( ((PARTICLE *) p), TYPE_SINKING )) {
	bRVAUpdate(p) = 0; /* Indicator for r,v,a update */
	}
#endif
    }

void combSinkAccrete(void *p1,void *p2)
{
    if (!(TYPETest( ((PARTICLE *) p1), TYPE_DELETED )) &&
        TYPETest( ((PARTICLE *) p2), TYPE_DELETED ) ) {
	((PARTICLE *) p1)-> fMass = ((PARTICLE *) p2)-> fMass;
	TYPEReset(((PARTICLE *) p1),TYPE_SINKING);
	pkdDeleteParticle( NULL, p1 );
	}
#ifdef SINKING
    else if (TYPETest( ((PARTICLE *) p2), TYPE_SINKING ) && !TYPETest( ((PARTICLE *) p1), TYPE_SINKING ))  {
	PARTICLE *p=p1,*q=p2;
	int j;
	bRVAUpdate(p) = 1;
	assert(!TYPETest(p,TYPE_SINK));
	assert(!TYPETest(q,TYPE_SINK));
	p->fTrueMass = q->fTrueMass;
	TYPESet(p, TYPE_SINKING);
	p->fMetals = -1; /* HACK -- flag in output for sinking state */
	p->iSinkingOnto = q->iSinkingOnto;
	p->dt = q->dt;
	p->iRung = q->iRung;
	p->fSinkingTime = q->fSinkingTime;

	p->vSinkingr0 = q->vSinkingr0;
        p->rSinking0Mag = q->rSinking0Mag;
	p->vSinkingTang0Mag = q->vSinkingTang0Mag;
	assert(q->vSinkingr0 < 0);
	for (j=0;j<3;j++) {
	    p->rSinking0Unit[j] = q->rSinking0Unit[j];
	    p->vSinkingTang0Unit[j] = q->vSinkingTang0Unit[j];
	    }
	}
    if (TYPETest( ((PARTICLE *) p1), TYPE_SINKING )) {
#ifdef SINKINGAVERAGE
	TYPEReset(((PARTICLE *) p1),TYPE_NEWSINKING);
#endif
	if (bRVAUpdate(p2)) {
	    PARTICLE *p=p1,*q=p2;
	    p->r[0] = q->r[0];
	    p->r[1] = q->r[1];
	    p->r[2] = q->r[2];
	    p->v[0] = q->v[0];
	    p->v[1] = q->v[1];
	    p->v[2] = q->v[2];
	    p->vPred[0] = q->vPred[0];
	    p->vPred[1] = q->vPred[1];
	    p->vPred[2] = q->vPred[2];
	    p->a[0] = q->a[0];
	    p->a[1] = q->a[1];
	    p->a[2] = q->a[2];
	    }
	}
#endif
    }

/*#define DBGIORDER 2000000 */
#define DBGIORDER -1

void SinkAccrete(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
    int i,j;
	double ifMass;
	PARTICLE *q;
	int iSinkEventType = SINK_EVENT_NULL;
	SINKEVENT SinkEvent;
#ifdef SINKING
	double drSink[3];
	double r2;

	for (j=0;j<3;j++) drSink[j] = -p->r[j];
#endif
	SinkEvent.time= smf->dTime;

#ifdef STARSINK
	/* Convert to total angular momentum for additions */
	SINK_Lx(p) += TRUEMASS(p)*(p->r[1]*p->v[2] - p->r[2]*p->v[1]);
	SINK_Ly(p) += TRUEMASS(p)*(p->r[2]*p->v[0] - p->r[0]*p->v[2]);
	SINK_Lz(p) += TRUEMASS(p)*(p->r[0]*p->v[1] - p->r[1]*p->v[0]);
#endif
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    if ( iOrderSink(q) == p->iOrder) {
#ifndef SINKING
		FLOAT mp,mq;
		mp = p->fMass; mq = q->fMass;
		ifMass = 1./(mp+mq);
		for (j=0;j<3;j++) {
		    p->r[j] = ifMass*(mp*p->r[j]+mq*q->r[j]);
		    p->v[j] = ifMass*(mp*p->v[j]+mq*q->v[j]);
		    p->a[j] = ifMass*(mp*p->a[j]+mq*q->a[j]);
		    }
#ifdef STARSINK
		SINK_Lz(p) += mq*(q->r[1]*q->v[2] - q->r[2]*q->v[1]); /* add to total L */
        SINK_Ly(p) += mq*(q->r[2]*q->v[0] - q->r[0]*q->v[2]);
		SINK_Lz(p) += mq*(q->r[0]*q->v[1] - q->r[1]*q->v[0]);
#endif
		p->u += ifMass*(mp*p->u+mq*q->u);
		p->fMass += mq;
		assert(q->fMass != 0);
		SinkLogEntry( &smf->pkd->sinkLog, &SinkEvent, iSinkEventType ); /* log past event */
		iSinkEventType = SINK_EVENT_ACCRETE;
		SinkEvent.iOrdSink = p->iOrder;
		SinkEvent.iOrdVictim = q->iOrder;
		q->fMass = 0;
		pkdDeleteParticle(smf->pkd, q);
#endif
#ifdef SINKING
        int bEat = 0;
		assert(TYPETest(q, TYPE_GAS));
		assert(!TYPETest(q, TYPE_SINK));
		r2 = nnList[i].fDist2;
		if (r2 < smf->dSinkMustAccreteRadius*smf->dSinkMustAccreteRadius ||
		    (TYPETest(q, TYPE_SINKING) &&
		     q->rSinking0Mag + q->vSinkingr0*(smf->dTime-q->fSinkingTime) < smf->dSinkMustAccreteRadius)) {
		    if (!TYPETest(q, TYPE_SINKING)) {
			FLOAT mp,mq;
			mp = p->fTrueMass; mq = q->fMass;
			/* Do a standard accrete -- no sinking stage */
			ifMass = 1./(mp+mq);
			for (j=0;j<3;j++) {
			    p->r[j] = ifMass*(mp*p->r[j]+mq*q->r[j]);
			    p->v[j] = ifMass*(mp*p->v[j]+mq*q->v[j]);
			    p->a[j] = ifMass*(mp*p->a[j]+mq*q->a[j]);
			    }
#ifdef STARSINK
			SINK_Lx(p) += mq*(q->r[1]*q->v[2] - q->r[2]*q->v[1]); /* add to total L */
			SINK_Ly(p) += mq*(q->r[2]*q->v[0] - q->r[0]*q->v[2]);
			SINK_Lz(p) += mq*(q->r[0]*q->v[1] - q->r[1]*q->v[0]);
#endif
			p->u += ifMass*(mp*p->u+mq*q->u);
			p->fMass += mq;
			p->fTrueMass += q->fMass;
			if (p->iOrder==DBGIORDER) printf("SINKING %d Direct Delete of %d +%g %g %g\n",p->iOrder,q->iOrder,q->fMass,p->fMass,p->fTrueMass);
			bEat = 1;
			assert(q->fMass != 0);
			SinkLogEntry( &smf->pkd->sinkLog, &SinkEvent, iSinkEventType ); /* log past event */
			iSinkEventType = SINK_EVENT_ACCRETE;
			SinkEvent.iOrdSink = p->iOrder;
			SinkEvent.iOrdVictim = q->iOrder;
			q->fMass = 0;
			pkdDeleteParticle(smf->pkd, q);
			}
		    else if (q->iSinkingOnto == p->iOrder) {
			/* Already sinking -- now do final accrete */
			p->fMass += q->fMass; /* Note: p->fTrueMass NOT modified */
			p->iSinkingOnto--; /* One less sinking onto this sink */
			if (p->iOrder==DBGIORDER) printf("SINKING %d Final Delete of %d +%g %g %g %d\n",p->iOrder,q->iOrder,q->fMass,p->fMass,p->fTrueMass,p->iSinkingOnto);
			bEat = 1;
			assert(q->fMass != 0);
			q->fMass = 0;
			TYPEReset(q,TYPE_SINKING); /* Don't want it to participate later */
			pkdDeleteParticle(smf->pkd, q);
			}
		    }
		else if (!TYPETest(q, TYPE_SINKING)) {
		    FLOAT mp,mq;
		    mp = p->fTrueMass; mq = q->fMass;
		    /* Enter sinking stage */
		    assert(r2 < smf->dSinkRadius*smf->dSinkRadius);
		    /* If this code is active, ignore non-infalling particles */
/*		    vr = 0;
		    for (j=0;j<3;j++) {
			vr +=  (q->r[j]-p->r[j])*(q->v[j]-p->v[j]);
			}
			if (vr >= 0) continue;*/

		    /* All force, velocity, position info associated
		       with particle now belong to sink instead.
		       Note: If we want accurate estimate of L for the sink
		       we should adjust the sink L now */
		    /* Everything short of eating the particle
		       If many additional particles are added sink position could vary wildly
		       -- defer trajectory calculation until after */
		    ifMass = 1./(p->fTrueMass+q->fMass);
		    for (j=0;j<3;j++) {
			p->r[j] = ifMass*(mp*p->r[j]+q->fMass*q->r[j]);
			p->v[j] = ifMass*(mp*p->v[j]+q->fMass*q->v[j]);
			p->a[j] = ifMass*(mp*p->a[j]+q->fMass*q->a[j]);
			}
#ifdef STARSINK
			SINK_Lx(p) += mq*(q->r[1]*q->v[2] - q->r[2]*q->v[1]); /* add to total L */
			SINK_Ly(p) += mq*(q->r[2]*q->v[0] - q->r[0]*q->v[2]);
			SINK_Lz(p) += mq*(q->r[0]*q->v[1] - q->r[1]*q->v[0]);
#endif
		    p->u += ifMass*(mp*p->u+mq*q->u);
		    p->fTrueMass += q->fMass;
		    p->iSinkingOnto++; /* One more sinking onto this sink */
		    SinkLogEntry( &smf->pkd->sinkLog, &SinkEvent, iSinkEventType ); /* log past event */
		    iSinkEventType = SINK_EVENT_ACCRETE;
		    SinkEvent.iOrdSink = p->iOrder;
		    SinkEvent.iOrdVictim = q->iOrder;
		    q->fTrueMass = 0; /* Sinking particles retain mass to act as ghost particles for forces g+SPH */
		    assert(!TYPETest(q,TYPE_SINK));
		    assert(TYPETest(q,TYPE_NEWSINKING));
		    TYPESet(q, TYPE_SINKING);
		    q->fMetals = -1; /* HACK -- flag for sinking state */
		    q->iSinkingOnto = p->iOrder;
		    q->dt = p->dt;
		    q->iRung = p->iRung;
		    q->fSinkingTime = smf->dTime;
		    q->vSinkingr0 = 1e30;
		    if (p->iOrder==DBGIORDER) printf("SINKING %d Add Sinking %d +%g %g %g %d\n",p->iOrder,q->iOrder,q->fMass,p->fMass,p->fTrueMass,p->iSinkingOnto);
		    bEat = 1;
		    }
#endif /*SINKING*/
		}
	    }
#ifdef STARSINK
	/* Store Internal angular momentum only as it is invariant under motions  */
	SINK_Lx(p) -= TRUEMASS(p)*(p->r[1]*p->v[2] - p->r[2]*p->v[1]);
	SINK_Ly(p) -= TRUEMASS(p)*(p->r[2]*p->v[0] - p->r[0]*p->v[2]);
	SINK_Lz(p) -= TRUEMASS(p)*(p->r[0]*p->v[1] - p->r[1]*p->v[0]);
#endif
	if (iSinkEventType != SINK_EVENT_NULL) {
	    iSinkEventType = SINK_EVENT_ACCRETE_UPDATE;  /* Final particle includes full Sink update */
	    SinkEvent.mass = TRUEMASS(p);
	    SinkEvent.r[0] = p->r[0];  SinkEvent.r[1] = p->r[1];  SinkEvent.r[2] = p->r[2];
	    SinkEvent.v[0] = p->v[0];  SinkEvent.v[1] = p->v[1];  SinkEvent.v[2] = p->v[2];
	    SinkEvent.L[0] = SINK_Lx(p); SinkEvent.L[1] = SINK_Ly(p); SinkEvent.L[2] = SINK_Lz(p);
	    SinkLogEntry( &smf->pkd->sinkLog, &SinkEvent, iSinkEventType ); /* log most recent event */
	}

#ifdef SINKING
	if (bEat) {
	    for (j=0;j<3;j++) drSink[j] += p->r[j]; /* If sink ate it moved: old sinking particles must follow */

	    for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		if (TYPETest( q, TYPE_SINKING ) && q->iSinkingOnto == p->iOrder) {
		    if (TYPETest(q,TYPE_NEWSINKING)) {
			double r0, dx[3], dv[3], vr, vrraw, vrmin, vt, norm, dot;
			/* Initialize sinking trajectory */
			r0 = 0;
			for (j=0;j<3;j++) {
			    dx[j] = (q->r[j]-p->r[j]);
			    r0 += dx[j]*dx[j];
			    }
			r0 = sqrt(r0);
			norm = 1/r0;
			for (j=0;j<3;j++) dx[j] *= norm;

			vr = 0;
			for (j=0;j<3;j++) {
			    dv[j] = q->v[j]-p->v[j];  /* dv = raw relative velocity */
			    vr += dv[j]*dx[j];
			    } /* raw infall v */
#ifdef SINKINGAVERAGE
			vrraw = vr;
			    /* See if smoothed velocity is still infalling
			       q->vSinkingTang0Unit temporarily contains mean of (q'-q) */
			vr = 0;
			for (j=0;j<3;j++) {
			    dv[j] = q->vSinkingTang0Unit[j]+q->v[j]-p->v[j]; /* dv now smoothed relative vel */
			    vr += dv[j]*dx[j];
			    }
			if (p->iOrder==DBGIORDER) printf("SINKING %d Add Sinking Trajectory vr %d %g %g %g\n",p->iOrder,q->iOrder,vr,vrraw,-q->c);
/*                    printf("Infall vr avg %g ",vr);*/
#endif
			vt = 0;
			for (j=0;j<3;j++) {
			    dv[j] -= vr*dx[j]; /*subtract off radial part*/
			    vt += dv[j]*dv[j];
			    }
			vt = sqrt(vt); /* Averaged tangential motion */
#ifdef SINKINGAVERAGE
			if (vr >= 0) vr=vrraw;
#endif
			vrmin = -q->c; /* Bondi-Hoyle minimum inflow speed -- perhaps not ideal for rotating case */
			if (vr > vrmin) vr = vrmin;
/*                    printf("min %g : %g   vt %g \n",-vrmin,vr,vt);*/

			q->rSinking0Mag = r0;
			if (r0 > smf->dSinkRadius) printf("WARNING: %i outside sink r=%g %g %g\n",q->iOrder,r0,sqrt(nnList[i].fDist2),smf->dSinkRadius);
			q->vSinkingr0 = vr;
			q->vSinkingTang0Mag = vt;
			norm = 1/(vt+1e-37);
			dot = 0;
			for (j=0;j<3;j++) {
			    q->rSinking0Unit[j] = dx[j];
			    q->vSinkingTang0Unit[j] = dv[j]*norm;
			    dot += q->rSinking0Unit[j]*q->vSinkingTang0Unit[j]; /* sanity check */
			    q->a[j] = p->a[j];
			    q->v[j] = p->v[j]; /* move with sink now */
			    }
			assert(fabs(dot) < 1e-4);
			}
		    else {
			assert(q->vSinkingr0 < 0);
			for (j=0;j<3;j++) {
			    q->r[j] += drSink[j];
			    q->a[j] = p->a[j];
			    q->v[j] = p->v[j];
			    }
			}
		    }
		}
	    }

	/* Make sure sinking particles are still on their trajectories */
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    if (TYPETest(q,TYPE_SINKING) && q->iSinkingOnto == p->iOrder) {
		FLOAT r0 = q->rSinking0Mag;
		FLOAT r2 = r0 + q->vSinkingr0*(smf->dTime-q->fSinkingTime);
		FLOAT thfac, th2, costh2, sinth2, dr2, sqr02;
		FLOAT dr[3];

		assert(q->vSinkingr0 < 0);
		if (r2 < 0.1*r0) r2 = 0.1*r0; /* HACK */
		thfac = q->vSinkingTang0Mag*2/(q->vSinkingr0);
		sqr02 = sqrt(r0/r2);
		th2 = thfac*(1-sqr02);
		costh2 = cos(th2);
		sinth2 = sin(th2);
		for (j=0;j<3;j++) {
		    q->r[j] = p->r[j]+r2*costh2*q->rSinking0Unit[j]+r2*sinth2*q->vSinkingTang0Unit[j];
		    q->vPred[j] = p->v[j]+q->vSinkingTang0Mag*sqr02*
			(-sinth2*q->rSinking0Unit[j]+costh2*q->vSinkingTang0Unit[j])
			+q->vSinkingr0*(costh2*q->rSinking0Unit[j]+sinth2*q->vSinkingTang0Unit[j]);
		    }
		bRVAUpdate(q) = 1; /* Indicator for r,v,a update */
		TYPEReset(q,TYPE_NEWSINKING);
		}

	    }
#endif
#endif
}

void initSinkingForceShare(void *p)
    { /* cached copies only */
#ifdef SINKING
    if (TYPETest( ((PARTICLE *) p), TYPE_SINKING )) {
	bRVAUpdate((PARTICLE *) p) = 0; /* Indicator for r,v,a update */
	}
#endif
    }

void combSinkingForceShare(void *p1,void *p2)
{
#ifdef SINKING
    if (TYPETest( ((PARTICLE *) p1), TYPE_SINKING )) {
	if (bRVAUpdate((PARTICLE *) p2)) {; /* Indicator for r,v,a update */
	    PARTICLE *p=p1,*q=p2;
	    p->r[0] = q->r[0];
	    p->r[1] = q->r[1];
	    p->r[2] = q->r[2];
	    p->v[0] = q->v[0];
	    p->v[1] = q->v[1];
	    p->v[2] = q->v[2];
	    p->vPred[0] = q->vPred[0];
	    p->vPred[1] = q->vPred[1];
	    p->vPred[2] = q->vPred[2];
	    p->a[0] = q->a[0];
	    p->a[1] = q->a[1];
	    p->a[2] = q->a[2];
	    }
	}
#endif
    }
void SinkingForceShare(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#if defined(GASOLINE) && defined(SINKING)
	int i,j,nSinking;
	double totmass,norm,ma[3];
	PARTICLE *q;

	assert(TYPETest(p, TYPE_SINK));
	assert(!TYPETest(p, TYPE_SINKING));

	/* Only do this if sink is active */
	if (!TYPEQueryACTIVE(p)) return;

	totmass = 0;
	nSinking = 0;
	for (j=0;j<3;j++) ma[j] = 0;

	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    if (TYPETest( q, TYPE_SINKING ) && q->iSinkingOnto == p->iOrder) {
		assert(p->iRung == q->iRung);
		assert(!TYPETest(q,TYPE_NEWSINKING));

		totmass += q->fMass;
		nSinking++;
		assert(q->vSinkingr0 < 0);
		if (p->iOrder==DBGIORDER) {
		    FLOAT r0 = q->rSinking0Mag;
		    FLOAT r2 = r0 + q->vSinkingr0*(smf->dTime-q->fSinkingTime);
		    if (p->iOrder == 2000000 && q->iSinkingOnto == 2000000) printf("FORCESHARE %d eats sinking %d %g %g %g\n",p->iOrder,q->iOrder,sqrt((p->r[0]-q->r[0])*(p->r[0]-q->r[0])+(p->r[1]-q->r[1])*(p->r[1]-q->r[1])+(p->r[2]-q->r[2])*(p->r[2]-q->r[2])), r0,r2 );
		    if (q->iOrder == 105639) assert(0);
		    }
		if (p->iOrder==DBGIORDER) printf("%d ForceShare with %d +%g %g  %g\n",p->iOrder,q->iOrder,q->fMass,totmass,sqrt(nnList[i].fDist2)/smf->dSinkRadius);
		ma[0] += q->fMass*q->a[0];
		ma[1] += q->fMass*q->a[1];
		ma[2] += q->fMass*q->a[2];
		}
	    }
	if (p->iOrder==DBGIORDER) printf("%d ForceShCount: %d %d\n",p->iOrder,nSinking,p->iSinkingOnto);
	if (nSinking < p->iSinkingOnto) {
	    p->fBall2 = -1.0; /* redo smooth */
	    return;
	    }
	assert(nSinking == p->iSinkingOnto);

	if (totmass > 0) {
	    totmass += p->fMass;
	    assert(fabs(totmass-p->fTrueMass) < 1e-7*totmass); /* sanity check on true mass
								  Fails if: sinking particle missed or true mass incorrect */
	    norm = 1/totmass;
	    for (j=0;j<3;j++) p->a[j] = (p->a[j]*p->fMass+ma[j])*norm;

	    /* Strictly, only a should need to be set but doing v makes v's exact */
	    for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		if (TYPETest( q, TYPE_SINKING ) && q->iSinkingOnto == p->iOrder) {

		    assert(TYPETest(q, TYPE_INFLOW)); /* Each sinking should appear here exactly once */
		    TYPEReset(q, TYPE_INFLOW ); /* DEBUG -- detect unfound particles */
		    TYPESet(q, TYPE_OUTFLOW );

		    bRVAUpdate( q ) = 1; /* Indicator for r,v,a update */
		    q->v[0] = p->v[0];
		    q->v[1] = p->v[1];
		    q->v[2] = p->v[2];
		    q->a[0] = p->a[0];
		    q->a[1] = p->a[1];
		    q->a[2] = p->a[2];
		    /* set absolute offset and velocity of sinking particles
		       strictly speaking, relative corrections in drift should be enough */
			{
			FLOAT r0 = q->rSinking0Mag;
			FLOAT r2 = r0 + q->vSinkingr0*(smf->dTime-q->fSinkingTime);
			FLOAT thfac, th2, costh2, sinth2, dr2, sqr02;
			FLOAT dr[3];

			if (r2 < 0.1*r0) r2 = 0.1*r0; /* HACK */
			thfac = q->vSinkingTang0Mag*2/(q->vSinkingr0);
			sqr02 = sqrt(r0/r2);
			th2 = thfac*(1-sqr02);
			costh2 = cos(th2);
			sinth2 = sin(th2);
			for (j=0;j<3;j++) {
			    q->r[j] = p->r[j]+r2*costh2*q->rSinking0Unit[j]+r2*sinth2*q->vSinkingTang0Unit[j];
			    q->vPred[j] = p->v[j]+q->vSinkingTang0Mag*sqr02*
				(-sinth2*q->rSinking0Unit[j]+costh2*q->vSinkingTang0Unit[j])
				+q->vSinkingr0*(costh2*q->rSinking0Unit[j]+sinth2*q->vSinkingTang0Unit[j]);
			    }
			}
		    }
		}
	    }
#endif
}

void initTreeParticleBHSinkDensity(void *p)
{
#ifdef GASOLINE

    ((PARTICLE *)p)->curlv[1] = 0.0; /* total mass change */
    ((PARTICLE *)p)->curlv[0] = ((PARTICLE *)p)->fMass; /* initial mass */
    ((PARTICLE *)p)->curlv[2] = 0.0; /* flag if accreted */
#endif
    }

/* Cached Tree Active particles */
void initBHSinkDensity(void *p)
{
#ifdef GASOLINE
    /*
     * Original particle curlv's is trashed (JW)
     */
    ((PARTICLE *)p)->curlv[1] = 0.0; /* total mass change */
    ((PARTICLE *)p)->curlv[0] = ((PARTICLE *)p)->fMass; /* initial mass */
    ((PARTICLE *)p)->curlv[2] = 0.0; /* flag if accreted */
#endif
    }

void combBHSinkDensity(void *p1,void *p2)
{
#ifdef GASOLINE
  if (((PARTICLE *)p2)->curlv[2] > ((PARTICLE *)p1)->curlv[2]) { /* only important if particle has been selected JMB 6/24/09
					 * but necessary - if particle has 2 BH neighbors it can be
					 * overwritten */
    ((PARTICLE *)p1)->curlv[0] = ((PARTICLE *)p2)->curlv[0];
    ((PARTICLE *)p1)->curlv[1] = ((PARTICLE *)p2)->curlv[1]; /* total mass change */
    ((PARTICLE *)p1)->curlv[2] = ((PARTICLE *)p2)->curlv[2];
  }


    /* this way the most recent one wins JMB 5/12/09 */
    /* don't sum up curlv[1], only let one BH eating be recorded.
       curlv[1] will never be greater than curlv[0] this way */
#endif
}

/*
 * Calculate paramters for BH accretion for use in the BHSinkAccrete
 * function.  This is done separately to even competition between
 * neighboring Black Holes.
 * Gas particles are ID'd here to be eaten in BHSinkAccrete.
 */
void BHSinkDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
	PARTICLE *q = NULL;

	FLOAT ih2,r2,rs,fDensity;
	FLOAT v[3],cs,fW,dv2,dv;
	double mdot, mdotEdd, mdotCurr, dmq, dtEff;
	int i,naccreted, ieat;
	FLOAT mdotsum, weat;
	FLOAT weight,wrs;
	double aFac, dCosmoDenFac,dCosmoVel2Fac,aDot;
	FLOAT dvmin,dvx,dvy,dvz,dvcosmo; /* measure mindeltav */

	assert(p->iRung >= smf->iSinkCurrentRung);
	p->curlv[1] = 0.0;
	naccreted = 0;
	aFac = smf->a;
	dCosmoDenFac = aFac*aFac*aFac;
	dCosmoVel2Fac = (smf->bCannonical ? aFac*aFac : 1.0);
	aDot = aFac*smf->H;  /* da/dt */
	/* CosmoVel2Fac converts cannonical velocity^2 to physical velocity^2
	   Cosmo sims tend to use cannonical velocities
	   v_phys = v_cannonical/a + adot*r
	   JMB 9/1/12*/

        mdotsum = 0.;
        weat = -1e37;

	ih2 = 4.0/BALL2(p);
	fDensity = 0.0; cs = 0;
	v[0] = 0; v[1] = 0; v[2] = 0;

	/* for mindv accretion method JMB 8/5/10 */
	if(smf->bBHMindv == 1){
	  dvmin = FLOAT_MAXVAL;
	  for (i=0;i<nSmooth;++i) {
	    dvx = (-p->v[0] + nnList[i].pPart->v[0])/aFac - aDot*nnList[i].dx;
	    dvy = (-p->v[1] + nnList[i].pPart->v[1])/aFac - aDot*nnList[i].dy;
	    dvz = (-p->v[2] + nnList[i].pPart->v[2])/aFac - aDot*nnList[i].dz;
	    dvcosmo = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
	    if (dvcosmo < dvmin) {
	      dvmin=dvcosmo;
	    }
	  }
	}

	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    q = nnList[i].pPart;
	    /*if(q->curlv[2] != 0) {*/
	    q->curlv[0] = q->fMass;
	    q->curlv[1] = 0;
	    q->curlv[2] = 0;
	    /*}*/

	    /* Sink should NOT be part of this list */
	    assert(TYPETest(q,TYPE_GAS));
	    fW = rs*q->fMass;

	    if(smf->bBHMindv == 1) weight = rs*pow(q->c*q->c+(dvmin*dvmin),-1.5)/dCosmoDenFac;
	    else {
	      dvx = (-p->v[0]+q->v[0])/aFac - aDot*nnList[i].dx;
	      dvy = (-p->v[1]+q->v[1])/aFac - aDot*nnList[i].dy;
	      dvz = (-p->v[2]+q->v[2])/aFac - aDot*nnList[i].dz;
	      dvcosmo = sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
	      weight = rs*pow(q->c*q->c+dvcosmo*dvcosmo,-1.5)/dCosmoDenFac; /* weight particles by mdot quantities */
	    /* cosmo factors put in 7/7/09  JMB */
	    }
	      if (weight > weat)  {
		weat = weight;
		ieat = i;
		wrs = rs;
	      }
	      mdotsum += weight*q->fMass;


	    fDensity += fW;
	    v[0] += fW*q->v[0];
	    v[1] += fW*q->v[1];
	    v[2] += fW*q->v[2];
	    cs += fW*q->c;
	    }

	dv2 = 0;
	for (i=0;i<3;i++) {
	    dv = v[i]/fDensity-p->v[i];
	    dv2 += dv*dv;
	    }
	dv2 /= dCosmoVel2Fac;
	/*
	 * Store results in particle.
	 * XXX NB overloading "curlv" field of the BH particle.  I am
	 * assuming it is not used.
	 */
	p->c = cs = cs/fDensity;
	p->fDensity = fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity;
	fDensity /= dCosmoDenFac;

	if(smf->bBHMindv == 1) printf("BHSink %d:  Time: %g Selected Particle Density: %g SPH: %g C_s: %g dv: %g dist: %g\n",p->iOrder,smf->dTime,nnList[ieat].pPart->fDensity,wrs,nnList[ieat].pPart->c,dvmin,sqrt(nnList[ieat].fDist2));
	else printf("BHSink %d:  Time: %g Selected Particle Density: %g SPH: %g C_s: %g dv: %g dist: %g\n",p->iOrder,smf->dTime,nnList[ieat].pPart->fDensity,wrs,nnList[ieat].pPart->c,sqrt(dv2),sqrt(nnList[ieat].fDist2));

	printf("BHSink %d:  Time: %g Mean Particle Density: %g C_s: %g dv: %g\n",p->iOrder,smf->dTime,fDensity,cs,sqrt(dv2));


	mdot = mdotsum*smf->dBHSinkAlphaFactor*p->fMass*p->fMass*M_1_PI*sqrt(ih2)*ih2; /* new mdot! */


	/* Eddington Limit Rate */
	mdotEdd = smf->dBHSinkEddFactor*p->fMass;
	printf("BHSink %d:  Time: %.8f mdot (BH): %g mdot (Edd): %g a: %g\n",p->iOrder,smf->dTime,mdot,mdotEdd, smf->a);

	if (mdot > mdotEdd) mdot = mdotEdd;

	mdotCurr = p->divv = mdot; /* store mdot in divv of sink */

        weight = -1e37;
	weat = -1e37;

	for (;;) {
	    q = NULL;
	    for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2;
		if(TYPETest(nnList[i].pPart,TYPE_DELETED)) continue;
		    /* don't accrete a deleted particle!  JMB 10/1/08 */
		KERNEL(rs,r2);
		fW = rs*nnList[i].pPart->fMass;
		if(smf->bBHAccreteAll != 1 && r2 > 0.25*nnList[i].pPart->fBall2) continue;
		/* don't accrete gas that doesn't have the BH
		 * in its smoothing length  JMB 10/22/08 */
		/* make this an optional parameter JMB 9/21/12 */
		/*if (nnList[i].pPart->iRung < smf->iSinkCurrentRung) continue; // JMB 7/9/09 */

		if(smf->bBHMindv == 1) weight = rs*pow(nnList[i].pPart->c*nnList[i].pPart->c+(dvmin*dvmin),-1.5)/dCosmoDenFac;
		else {
		  dvx = (-p->v[0]+nnList[i].pPart->v[0])/aFac - nnList[i].dx*aDot;
		  dvy = (-p->v[1]+nnList[i].pPart->v[1])/aFac - nnList[i].dy*aDot;
		  dvz = (-p->v[2]+nnList[i].pPart->v[2])/aFac - nnList[i].dz*aDot;
		  dvcosmo = sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
		  weight = rs*pow(nnList[i].pPart->c*nnList[i].pPart->c+dvcosmo*dvcosmo,-1.5)/dCosmoDenFac; /* weight particles by mdot quantities */
	    /* cosmo factors put in 7/7/09  JMB */
		}

		if (weight > weat && nnList[i].pPart->curlv[2] != p->iOrder) {
		  weat = weight;
		  q = nnList[i].pPart;
		}
	    }
	    assert(q != p); /* shouldn't happen because p isn't in
			       tree */


	    if(q != NULL) {
	      dtEff = smf->dSinkCurrentDelta;/* *pow(0.5,iRung-smf->iSinkCurrentRung);*/
	      /* If victim has unclosed kick -- don't actually take the mass
	       If sink has unclosed kick we shouldn't even be here!
	       When victim is active use his timestep if longer
	       Statistically expect to get right effective mdot on average */

	      /* JMB 9/18/12 -- Actually, the method described above
		 does NOT get the right mdot on average, it
		 underestimates mdot.  Instead I have removed all
		 timestep criteria from BH accretion.  Thus, momentum
		 is no longer strictly conserved.  HOWEVER, the
		 amounts of momentum involved are tiny, likely
		 comparable to the errors intrinsic to the code.
		 Thus, we shall not stress about it. */
	      dmq = mdotCurr*dtEff;
	      }
	    else {
	      if(naccreted != 0) printf("BHSink %d:  only %i accreted, %g uneaten mass\n",p->iOrder,naccreted,mdotCurr*dtEff);
	      /* JMB 6/29/09 */
	      if(naccreted == 0) p->curlv[1] = 0.0; /* should be anyway */
	      return;
	    }
	    assert(q != NULL);
	    /* We have our victim */


	      if (dmq < q->fMass) {
		mdotCurr = 0.0;
	      }
	      else {
		mdotCurr -= mdotCurr*(q->fMass/dmq); /* need an additional victim */
		assert(mdotCurr >= 0.0);
		dmq = q->fMass;
	      }
	      q->curlv[2] = p->iOrder; /* flag for pre-used victim particles */
	      /* temporarily store mass lost -- to later check for double dipping */
	      q->curlv[1] += dmq;
	      p->curlv[1] += dmq;
	      naccreted += 1;
	      weat = -1e37; /* reset so other particles can be selected JMB 7/9/09 */

	      printf("BHSink %d:  Time %g %d dmq %g %g %g\n",p->iOrder,smf->dTime,q->iOrder,dmq,q->curlv[1],p->curlv[1]);

	      /*}*/

	    if (mdotCurr == 0.0) break;
	    }
#endif
    }


void initTreeParticleBHSinkAccrete(void *p1)
{
#ifdef STARFORM
    /* Convert energy and metals to non-specific quantities (not per mass)
     * to make it easier to divvy up BH energy.
     */

       if(TYPETest((PARTICLE *)p1, TYPE_GAS))  ((PARTICLE *)p1)->u *= ((PARTICLE *)p1)->fMass;
     if(TYPETest((PARTICLE *)p1, TYPE_GAS))  ((PARTICLE *)p1)->fNSN = 0.0;
  /* fNSN will hold the place of FB energy because uPred is needed for
   *  calculations. it should be zero anyways.  JMB 10/5/09  */

#endif
    }

/* Cached Tree Active particles */
void initBHSinkAccrete(void *p)
{
#ifdef STARFORM
    if (TYPEQueryTREEACTIVE((PARTICLE *) p)) ((PARTICLE *)p)->u = 0.0; /*added 6/10/08*/
  if (TYPEQueryTREEACTIVE((PARTICLE *) p)) ((PARTICLE *)p)->fNSN = 0.0;

	/* Heating due to accretion */
#endif
    }

void combBHSinkAccrete(void *p1,void *p2)
{
#ifdef STARFORM
    PARTICLE *pp1 = p1;
    PARTICLE *pp2 = p2;

    if (!(TYPETest( pp1, TYPE_DELETED )) &&
        TYPETest( pp2, TYPE_DELETED ) ) {
	pp1->fMass = pp2->fMass;
	pkdDeleteParticle( NULL, pp1 );
	}
    else if (TYPEQueryTREEACTIVE(pp1)) {
	/*
	 * See kludgery notice above: record eaten mass in original
	 * particle.
	 */
	FLOAT fEatenMass = pp2->curlv[0] - pp2->fMass;
	pp1->fMass -= fEatenMass;
	if(pp1->fMass < pp1->curlv[0]*1e-3) {
	    /* This could happen if BHs on two
	       different processors are eating
	       gas from a third processor */
	    fprintf(stderr, "ERROR: Overeaten gas particle %d: %g %g\n",
		    pp1->iOrder,
		    pp1->fMass, fEatenMass);
	    if (!(TYPETest( pp1, TYPE_DELETED ))) {
		pkdDeleteParticle( NULL, pp1 );
		}
	    return;
	    }

	if(!(TYPETest(pp1,TYPE_DELETED)))  { /*added 8/21/08 to check -0 mass particles  */
	    assert(pp1->fMass > 0.0);
	}
	pp1->u += pp2->u; /*added 6/10/08 for BH blastwave FB*/
	pp1->fNSN += pp2->fNSN;  /* (this is uPred) JMB 10/5/09  */
        pp1->fTimeCoolIsOffUntil = max( (pp1)->fTimeCoolIsOffUntil,
                (pp2)->fTimeCoolIsOffUntil );
	pp1->fMassForm = max( (pp1)->fMassForm,
			      (pp2)->fMassForm ); /* propagate FB time JMB 2/24/10 */
	}
#endif
}

void BHSinkAccrete(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef STARFORM
	PARTICLE *q = NULL;

	FLOAT ih2,r2,rs;
	double mdot, mdotCurr, dm, dmq, dE, ifMass, dtEff;
	FLOAT fNorm,fNorm_new;
	int i,naccreted;
	FLOAT weat;
	FLOAT weight,fbweight; /* weight is for accretion, fbweight is for FB  */
	double aFac, dCosmoDenFac,aDot;
	FLOAT dvmin, dvx,dvy,dvz,dvcosmo;

        weat = -1e37;
	aFac = smf->a;
	dCosmoDenFac = aFac*aFac*aFac;
	aDot = aFac*smf->H;

	mdot = p->divv;

	mdotCurr = mdot;
	dm = 0;
	naccreted = 0;
       	ih2 = 4.0/BALL2(p);

	for (i=0;i<nSmooth;++i){
	    FLOAT r2min = FLOAT_MAXVAL;
	    q = NULL;
	    if (nnList[i].pPart->curlv[2] == 0) continue;
	    if(TYPETest(nnList[i].pPart,TYPE_DELETED)) continue;
	    if (nnList[i].pPart->curlv[2] == p->iOrder) q = nnList[i].pPart;
	    else continue; /* won't work if selected by another BH */
	    assert(q != NULL);
	    r2min = nnList[i].fDist2;

	    /* Timestep for accretion is larger of sink and victim timestep */
	    dtEff = smf->dSinkCurrentDelta;
	      /* If victim has unclosed kick -- don't actually take the mass
	       If sink has unclosed kick we shouldn't even be here!
	       When victim is active use his timestep if longer
	       Statistically expect to get right effective mdot on average */

	    /* JMB 9/18/12 - see comment in BHSinkDensity above for
	       why the previous comment is not actually what is done. */
	    dmq = mdotCurr*dtEff;

	    /* update mdotCurr */
	    if( fabs(dmq - q->curlv[1])/dmq < 0.0001) mdotCurr = 0.0;
	    /* some floating point issues here JMB 7/9/09 */
	    else {
	      mdotCurr -= q->curlv[1]/dtEff;
	      dmq = q->curlv[1];
	    }
	    assert(mdotCurr >= 0.0);

	    printf("BHSink %d:  Time %g %d dmq %g %g %g\n",p->iOrder,smf->dTime,q->iOrder,dmq,q->curlv[1],p->curlv[1]);
		ifMass = 1./(p->fMass + dmq);
		/* to record angular momentum JMB 11/9/10 */
		printf("BHSink %d:  Gas: %d  dx: %g dy: %g dz: %g \n",p->iOrder,q->iOrder,p->r[0]-q->r[0],p->r[1]-q->r[1],p->r[2]-q->r[2]);
		printf("BHSink %d:  Gas: %d  dvx: %g dvy: %g dvz: %g \n",p->iOrder,q->iOrder,p->v[0]-q->v[0],p->v[1]-q->v[1],p->v[2]-q->v[2]);
		/* Adjust sink properties (conserving momentum etc...) */
		p->r[0] = ifMass*(p->fMass*p->r[0]+dmq*q->r[0]);
		p->r[1] = ifMass*(p->fMass*p->r[1]+dmq*q->r[1]);
		p->r[2] = ifMass*(p->fMass*p->r[2]+dmq*q->r[2]);
		p->v[0] = ifMass*(p->fMass*p->v[0]+dmq*q->v[0]);
		p->v[1] = ifMass*(p->fMass*p->v[1]+dmq*q->v[1]);
		p->v[2] = ifMass*(p->fMass*p->v[2]+dmq*q->v[2]);
		p->a[0] = ifMass*(p->fMass*p->a[0]+dmq*q->a[0]);
		p->a[1] = ifMass*(p->fMass*p->a[1]+dmq*q->a[1]);
		p->a[2] = ifMass*(p->fMass*p->a[2]+dmq*q->a[2]);
		p->fMetals = ifMass*(p->fMass*p->fMetals+dmq*q->fMetals);
        if (dmq > q->fMass) //PRC: this if-else added to avoid negative mass particles
        {
            p->fMass += q->fMass;
            dm += q->fMass;
            q->fMass = 0.;
        }
        else
        {
            q->fMass -= dmq;
            p->fMass += dmq;
            dm += dmq;
            assert(q->fMass >= 0.0);
        }
		naccreted += 1;  /* running tally of how many are accreted JMB 10/23/08 */
		printf("BHSink %d:  Time %g dist2: %d %g gas smooth: %g eatenmass %g \n",p->iOrder,smf->dTime,q->iOrder,r2min,q->fBall2,dmq);
		if (q->fMass <= 1e-3*dmq || q->fMass <= smf->dMinGasMass) { /* = added 8/21/08 */
		    // q->fMass = 0;
		    if(!(TYPETest(q,TYPE_DELETED))) pkdDeleteParticle(smf->pkd, q);
		    /* Particles are getting deleted twice, which is messing
		      up the bookkeeping.  I think this will make them only be
		      deleted once.  JMB 9/23/08*/
		    }

		/*	}*/

	}

	if (mdotCurr == 0.0) goto dofeedback;
	else{
	  /* go looking for more particles, sharing must happen */
	    /****** SHARING CODE HERE ***********/

	  if(smf->bBHMindv == 1){
	    dvmin = FLOAT_MAXVAL;
	    for (i=0;i<nSmooth;++i) {
	      dvx = (-p->v[0] + nnList[i].pPart->v[0])-aDot*nnList[i].dx;
	      dvy = (-p->v[1] + nnList[i].pPart->v[1])-aDot*nnList[i].dy;
	      dvz = (-p->v[2] + nnList[i].pPart->v[2])-aDot*nnList[i].dz;
	      dvcosmo = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
	      if (dvcosmo < dvmin) {
		dvmin=dvcosmo;
	      }
	    }
	  }

	  for (;;) {
	    FLOAT r2min = FLOAT_MAXVAL;
	    q = NULL;
	    for (i=0;i<nSmooth;++i) {
	      if(TYPETest(nnList[i].pPart,TYPE_DELETED)) continue;
	      if(nnList[i].pPart->curlv[2] == p->iOrder) continue;
	      /* don't choose a pre-accreted particle
	         but it can have been accreted by another BH */
	      /*if(nnList[i].pPart->iRung < smf->iSinkCurrentRung) continue;*/
	      /* has to be on the right timestep */
	      r2 = nnList[i].fDist2;
	      if(smf->bBHAccreteAll != 1 && r2 > 0.25*nnList[i].pPart->fBall2) continue;
	      /* has to be nearby! */
	      KERNEL(rs,r2);

	      if(smf->bBHMindv == 1) weight = rs*pow(nnList[i].pPart->c*nnList[i].pPart->c+(dvmin*dvmin),-1.5)/dCosmoDenFac;
	      else {
		dvx = (-p->v[0]+nnList[i].pPart->v[0])-aDot*nnList[i].dx;
		dvy = (-p->v[1]+nnList[i].pPart->v[1])-aDot*nnList[i].dy;
		dvz = (-p->v[2]+nnList[i].pPart->v[2])-aDot*nnList[i].dz;
		dvcosmo = sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
		weight = rs*pow(nnList[i].pPart->c*nnList[i].pPart->c+dvcosmo*dvcosmo,-1.5)/dCosmoDenFac; /* weight particles by mdot quantities */
		/* cosmo factors put in 7/7/09  JMB */
	      }
	      if (weight > weat) {
		r2min = r2; /* note r2min is not really the min r2 anymore */
		weat = weight;
		q = nnList[i].pPart;
	      }
	    }

	    weat = -1e37; /* reset so other particles can be selected JMB 7/9/09 */

	    if(q == NULL) {
	      dtEff = smf->dSinkCurrentDelta;
	      printf("BHSink %d:  WARNING!! Not enough edible particles.  Time: %g Mass not eaten: %g \n",p->iOrder,smf->dTime,mdotCurr*dtEff);
	      if(naccreted == 0) return;
	      else goto dofeedback;
	    }
	    else {
	      assert(q != NULL);
	      dtEff = smf->dSinkCurrentDelta;
	      dmq = mdotCurr*dtEff;


	    if (dmq < q->curlv[0]) mdotCurr = 0.0; /* Original mass in q->curlv[0] */
	    else {
		mdotCurr -= mdotCurr*(q->curlv[0]/dmq); /* need an additional victim */
		dmq = q->curlv[0];
	      }


	    printf("BHSink %d:  Time: %g %d dmq %g %g %g sharing \n",p->iOrder,smf->dTime,q->iOrder,dmq,q->curlv[1],p->curlv[1]);
		ifMass = 1./(p->fMass + dmq);
		/* to record angular momentum JMB 11/9/10 */
		printf("BHSink %d:  Gas: %d  dx: %g dy: %g dz: %g \n",p->iOrder,q->iOrder,p->r[0]-q->r[0],p->r[1]-q->r[1],p->r[2]-q->r[2]);
		printf("BHSink %d:  Gas: %d  dvx: %g dvy: %g dvz: %g \n",p->iOrder,q->iOrder,p->v[0]-q->v[0],p->v[1]-q->v[1],p->v[2]-q->v[2]);
		/* Adjust sink properties (conserving momentum etc...) */
		p->r[0] = ifMass*(p->fMass*p->r[0]+dmq*q->r[0]);
		p->r[1] = ifMass*(p->fMass*p->r[1]+dmq*q->r[1]);
		p->r[2] = ifMass*(p->fMass*p->r[2]+dmq*q->r[2]);
		p->v[0] = ifMass*(p->fMass*p->v[0]+dmq*q->v[0]);
		p->v[1] = ifMass*(p->fMass*p->v[1]+dmq*q->v[1]);
		p->v[2] = ifMass*(p->fMass*p->v[2]+dmq*q->v[2]);
		p->a[0] = ifMass*(p->fMass*p->a[0]+dmq*q->a[0]);
		p->a[1] = ifMass*(p->fMass*p->a[1]+dmq*q->a[1]);
		p->a[2] = ifMass*(p->fMass*p->a[2]+dmq*q->a[2]);
		p->fMetals = ifMass*(p->fMass*p->fMetals+dmq*q->fMetals);
        if (dmq > q->fMass) //PRC: this if-else added to avoid negative mass particles
        {
            p->fMass += q->fMass;
            dm += q->fMass;
            q->fMass = 0.;
        }
        else
        {
            q->fMass -= dmq;
            p->fMass += dmq;
            dm += dmq;
            assert(q->fMass >= 0.0);
        }
		naccreted += 1;  /* running tally of how many are accreted JMB 10/23/08 */
		printf("BHSink %d:  Time %g dist2 %d %g gas smooth: %g eatenmass %g\n",p->iOrder,smf->dTime,q->iOrder,r2min,q->fBall2,dmq);
		if (q->fMass <= 1e-3*dmq || q->fMass <= smf->dMinGasMass) { /* = added 8/21/08 */
		    // q->fMass = 0;
		    if(!(TYPETest(q,TYPE_DELETED))) pkdDeleteParticle(smf->pkd, q);
		    /* Particles are getting deleted twice, which is messing
		      up the bookkeeping.  I think this will make them only be
		      deleted once.  JMB 9/23/08*/
		    }
		/*}*/
	    if (mdotCurr == 0.0) break;

	    }
	  }
	}

	/********************************************************************/
        dofeedback:

	if(smf->dBHSinkFeedbackFactor != 0.0) {  /* allows FB to be turned off JMB 7/20/09 */
	  dE = smf->dBHSinkFeedbackFactor*dm; /* dE based on actual mass eaten */

	  dtEff = smf->dSinkCurrentDelta;
	  printf("BHSink %d:  Delta: %g Time: %.8f dm: %g dE %g\n",p->iOrder,dtEff,smf->dTime,dm,dE);

	  /* Recalculate Normalization */
	  ih2 = 4.0/BALL2(p);
      fNorm_new = 0.0;
	  fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;

	  for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
            fNorm_new += rs;
            rs *= fNorm;
	  }

	  /* spread feedback to nearest 32 particles across the kernel.  JMB 8/29/22 */

	  assert(fNorm_new != 0.0);
          fNorm_new = 1./fNorm_new;

          for(i=0;i<nSmooth;++i) {
            r2 = nnList[i].fDist2*ih2;
            KERNEL(rs,r2);
            fbweight = rs*fNorm_new;
            nnList[i].pPart->u += fbweight*dE;
            nnList[i].pPart->fNSN += fbweight*dE;
            /* now turn off cooling */
            if(  smf->bBHTurnOffCooling) nnList[i].pPart->fTimeCoolIsOffUntil = max(nnList[i].pPart->fTimeCoolIsOffUntil,smf->dTime + nnList[i].pPart->dt);
            printf("BHSink %d: Time %g FB Energy to %i dE %g tCoolOffUntil %g \n",p->iOrder,smf->dTime,nnList[i].pPart->iOrder,fbweight*dE,nnList[i].pPart->fTimeCoolIsOffUntil);
	    if(smf->bBHTurnOffCooling) assert(nnList[i].pPart->fTimeCoolIsOffUntil > smf->dTime);
            nnList[i].pPart->fMassForm = (FLOAT) smf->dTime;
            /* track BHFB time in MassForm  JMB 11/19/08 */
          }
        }

#endif
}

void postBHSinkAccrete(PARTICLE *p1, SMF *smf)
{
#ifdef STARFORM
    /* Convert energy  back to specific quantities (per mass)
       because we are done with our conservative calculations */
    if(TYPETest(p1, TYPE_GAS) && p1->fMass != 0 && !(TYPETest(p1,TYPE_DELETED))) {
    p1->u /= p1->fMass;
    p1->fNSN /= p1->fMass;
    p1->uPred += p1->fNSN; /* now combine fNSN to uPred.  JMB 10/5/09  */
  }

#endif
}



void initBHSinkIdentify(void *p)
{
#ifdef STARFORM
  PARTICLE *pp = p;
  pp->fNSNtot = 0;
  /*  fNSNtot is a placeholder for merging info*/
#endif
}

void combBHSinkIdentify(void *p1, void *p2)
{
#ifdef STARFORM
    PARTICLE *pp1 = p1;
    PARTICLE *pp2 = p2;
    if(pp2->fNSNtot > pp1->fNSNtot) pp1->fNSNtot = pp2->fNSNtot;

#endif
}

void BHSinkIdentify(PARTICLE *p,int nSmooth,NN *nnList, SMF *smf)
{
  /* Identify BHs for merging JMB */

#ifdef STARFORM
	PARTICLE *q = NULL;
	int i;
	FLOAT deltaa, deltar, deltav;
	FLOAT aFac = smf->a;

	for (i=0;i<nSmooth;++i) {
	  q = nnList[i].pPart;

	  if(nnList[i].fDist2 > 4.0*p->fSoft*p->fSoft || q->iOrder == p->iOrder) continue;
	  /* are they close together? (within TWO softenings JMB 10/21/10*/
	  /* don't include yourself!  JMB 12/11/08 */
	  deltaa=sqrt( (p->a[0] - q->a[0])*(p->a[0] - q->a[0]) +  (p->a[1] - q->a[1])*(p->a[1] - q->a[1]) + (p->a[2] - q->a[2])*(p->a[2] - q->a[2]));
	  deltar=sqrt( (p->r[0] - q->r[0])*(p->r[0] - q->r[0]) +  (p->r[1] - q->r[1])*(p->r[1] - q->r[1]) + (p->r[2] - q->r[2])*(p->r[2] - q->r[2]));
	  deltav=sqrt( (p->v[0] - q->v[0])*(p->v[0] - q->v[0]) +  (p->v[1] - q->v[1])*(p->v[1] - q->v[1]) + (p->v[2] - q->v[2])*(p->v[2] - q->v[2]));

	  if ( deltaa*deltar*aFac < 0.5*deltav*deltav ) continue;
	  /* Selects other BH particles that are
	   * within  the criteria
	   * delta_a*delta_r < .5*delta_v^2
	   AND
	   * within the softening  */

	  if(p->iOrder < q->iOrder) {
	    if(p->iOrder < q->fNSNtot || q->fNSNtot == 0) {
	      q->fNSNtot = p->iOrder;
	      printf("BHSink MergeID %d will be eaten by %d \n",q->iOrder,p->iOrder);

	    /* fNSNtot is the place holder for the iord of the
	     eating black hole.  the victim remembers who eats
	     it for the merging step which is next.  JMB 4/30/09 */
	    }
	  }
	}

#endif

}

void initBHSinkMerge(void *p)
{
#ifdef GASOLINE

 /* this init function does nothing.  JMB */

#endif
}

void BHSinkMerge(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
  /* Makes BHs merge together, based on the criteria
     set in BHSinkIdentify above  */

#ifdef STARFORM

	PARTICLE *q = NULL;
	FLOAT ifMass;
	int i;
	/*  for kicks */
	const FLOAT A = 12000.0;
	const FLOAT B = -0.93;
	const FLOAT H = 7300.0;
	const FLOAT K = 60000.0;  /* BH kick parameters in km/s  */
	FLOAT vkick, vm, vpar, vkx, vky, vkz;
	FLOAT vperp, mratio, mfactor, cosine, spin1, spin2, angle1, angle2, xi;
	double xrand, yrand, zrand;
#define SECONDSPERYEAR   31557600.


	for (i=0;i<nSmooth;++i) {
	  q = nnList[i].pPart;

	  /* Give mass to higher iOrd particle, delete lower. */
	  if(q->fNSNtot == p->iOrder) {
	    /* q has already been identified as a victim of p */
	    ifMass = 1./(p->fMass + q->fMass);
	    mratio = p->fMass/q->fMass;
	    if(mratio > 1.0) mratio = q->fMass/p->fMass;
	    /* Adjust sink properties (conserving momentum etc...) */
	    p->r[0] = ifMass*(p->fMass*p->r[0]+q->fMass*q->r[0]);
	    p->r[1] = ifMass*(p->fMass*p->r[1]+q->fMass*q->r[1]);
	    p->r[2] = ifMass*(p->fMass*p->r[2]+q->fMass*q->r[2]);
	    p->v[0] = ifMass*(p->fMass*p->v[0]+q->fMass*q->v[0]);
	    p->v[1] = ifMass*(p->fMass*p->v[1]+q->fMass*q->v[1]);
	    p->v[2] = ifMass*(p->fMass*p->v[2]+q->fMass*q->v[2]);
	    p->a[0] = ifMass*(p->fMass*p->a[0]+q->fMass*q->a[0]);
	    p->a[1] = ifMass*(p->fMass*p->a[1]+q->fMass*q->a[1]);
	    p->a[2] = ifMass*(p->fMass*p->a[2]+q->fMass*q->a[2]);
	    p->fMass += q->fMass;


	      /**** Gravitational Recoil  ****/
	    if(smf->bDoBHKick == 1 && ( fabs(-1.0*p->fTimeForm - smf->dTime) > smf->dDeltaStarForm)) {
	      /* Turn off recoil if BH has just formed. Ejecting them immediately is not helpful. JMB 8/5/09  */
	      mfactor = pow(mratio,2)/pow((1+mratio),5);
	      vm = A * (1.0-mratio) * mfactor * (1.0 + B*mratio/pow((1.0+mratio),2));
	      cosine = (rand()/((double) RAND_MAX )) * 2. * 3.14159; /* angle */
	      spin1 = (rand()/((double) RAND_MAX )); /* spins of BHs 1 and 2 */
	      spin2 = (rand()/((double) RAND_MAX ));
	      angle1 = (rand()/((double) RAND_MAX )) * 2. * 3.14159;  /* orientations of BH spins */
	      angle2 = (rand()/((double) RAND_MAX )) * 2. * 3.14159;
	      xi = 1.5358897; /* 88.0 * 2. * 3.14159 / 360.0  */
	      /* in Campanelli et al xi = 88 degrees */

	      vperp = H * mfactor * (spin2*cos(angle2) - mratio*spin1*cos(angle1));  /* check sin cos! */
	      vpar = K * cos(cosine) * mfactor * (spin2*sin(angle2) - mratio*spin1*sin(angle1));

	      vkick = sqrt ( (vm + vperp*cos(xi))*(vm + vperp*cos(xi)) + (vperp*sin(xi))*(vperp*sin(xi)) + vpar*vpar);
	      vkick = vkick / smf->dKmPerSecUnit;
	      /* comoving coords are important JMB 5/15/09 */
	      if(smf->bCannonical) vkick = vkick*smf->a;

	      /* random direction */

	      for(;;){
		xrand = (rand()/((double) RAND_MAX )) * 2.0 - 1.0;
		yrand = (rand()/((double) RAND_MAX )) * 2.0 - 1.0;
		zrand = (rand()/((double) RAND_MAX )) * 2.0 - 1.0;
		if(sqrt(xrand*xrand+yrand*yrand+zrand*zrand) < 1.0) break;
	      }
	      vkx = vkick * xrand / sqrt(xrand*xrand+yrand*yrand+zrand*zrand);
	      vky = vkick * yrand / sqrt(xrand*xrand+yrand*yrand+zrand*zrand);
	      vkz = vkick * zrand / sqrt(xrand*xrand+yrand*yrand+zrand*zrand);

	      p->v[0] += vkx;
	      p->v[1] += vky;
	      p->v[2] += vkz;
	    }
	    else vkick = 0.0;

	    printf("BHSink Merge %d eating %d  Time %g kick velocity %g mass ratio %g \n",p->iOrder,q->iOrder,smf->dTime,vkick,mratio);

	    if(q->fNSNtot > 0 && !(TYPETest(q,TYPE_DELETED))) {
	      printf("BHSink Merge %d eaten Time %g \n",q->iOrder,smf->dTime);
	      pkdDeleteParticle(smf->pkd,q);
	    }
	  }
	}
#endif
}

void combBHSinkMerge(void *p1,void *p2)
{
#ifdef STARFORM
    PARTICLE *pp1 = p1;
    PARTICLE *pp2 = p2;

    if (!(TYPETest( pp1, TYPE_DELETED )) &&
        TYPETest( pp2, TYPE_DELETED ) ) {
	pkdDeleteParticle( NULL, pp1 );
	}
#endif
}

void initSinkMergeTest(void *p)
{
#ifdef GASOLINE
    fBindingEnergy(p) = FLT_MAX;
    iOrderSink(p) = -1;
#endif
}

void combSinkMergeTest(void *p1, void *p2)
{
#ifdef GASOLINE
/* Particle p1 belongs to sink iOrderSink(p1) initially but will
   switch to iOrderSink(p2) if more bound to that sink
   Mergers only complete if both particles agree on who to merge with
   If any candidate pairs, should always be one "most bound" pair that merge. */
    if (fBindingEnergy(p2) < fBindingEnergy(p1)) {
	fBindingEnergy(p1) = fBindingEnergy(p2);
	iOrderSink(p1) = iOrderSink(p2);
	}
#endif
}

void SinkMergeTest(PARTICLE *p,int nSmooth,NN *nnList, SMF *smf)
{
  /* Identify sinks for merging */

#ifdef GASOLINE
	PARTICLE *q = NULL;
	int i;
	FLOAT dv2, rmerge2, r2, EBind;
	rmerge2 = smf->dSinkRadius*smf->dSinkRadius; /* Centre of sink must be inside other sink */

	for (i=0;i<nSmooth;++i) {
	  q = nnList[i].pPart;
	  r2 = nnList[i].fDist2;
	  if(r2 > rmerge2 || q->iOrder == p->iOrder) continue;
	  dv2=( (p->v[0] - q->v[0])*(p->v[0] - q->v[0]) +  (p->v[1] - q->v[1])*(p->v[1] - q->v[1]) + (p->v[2] - q->v[2])*(p->v[2] - q->v[2]));
	  EBind = -TRUEMASS(p)*TRUEMASS(q)/sqrt(r2)
	      + 0.5*(TRUEMASS(p)*TRUEMASS(q))/(TRUEMASS(p)+TRUEMASS(q))*dv2;

	  /* Reject merger if unbound and outside must accrete radius */
	  if ( EBind > 0 && r2 > smf->dSinkMustAccreteRadius*smf->dSinkMustAccreteRadius) continue;

	  /* Only the most bound pair will merge each step */
	  if (EBind < fBindingEnergy(p) && EBind < fBindingEnergy(q)) {
	      fBindingEnergy(p) = EBind;
	      fBindingEnergy(q) = EBind;
	      iOrderSink(p) = q->iOrder;
	      iOrderSink(q) = p->iOrder;
	      printf("Sink Merge: ID's %d & %d flagged to merge, Ebind=%g\n",p->iOrder,q->iOrder,EBind);
	      }
	}

#endif

}

void initSinkMerge(void *p)
{
}

void combSinkMerge(void *p1, void *p2)
{
    if (!(TYPETest( ((PARTICLE *) p1), TYPE_DELETED )) &&
        TYPETest( ((PARTICLE *) p2), TYPE_DELETED ) ) {
	((PARTICLE *) p1)-> fMass = ((PARTICLE *) p2)-> fMass;
	pkdDeleteParticle( NULL, p1 );
	}
}

void SinkMerge(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
	PARTICLE *q = NULL;
	int i,j;
	FLOAT ifMass;
	SINKEVENT SinkEvent;

	for (i=0;i<nSmooth;++i) {
	  q = nnList[i].pPart;

	  if (iOrderSink(q) == p->iOrder && iOrderSink(p) == q->iOrder
	      && p->iOrder < q->iOrder) {
	      FLOAT mp,mq;
	      mp = p->fMass; mq = q->fMass;
#ifdef STARSINK
	/* Convert to total angular momentum for additions */
	      SINK_Lx(p) += TRUEMASS(p)*(p->r[1]*p->v[2] - p->r[2]*p->v[1]);
	      SINK_Ly(p) += TRUEMASS(p)*(p->r[2]*p->v[0] - p->r[0]*p->v[2]);
	      SINK_Lz(p) += TRUEMASS(p)*(p->r[0]*p->v[1] - p->r[1]*p->v[0]);
        SINK_Lx(p) += mq*(q->r[1]*q->v[2] - q->r[2]*q->v[1]); /* add to total L */
        SINK_Ly(p) += mq*(q->r[2]*q->v[0] - q->r[0]*q->v[2]);
        SINK_Lz(p) += mq*(q->r[0]*q->v[1] - q->r[1]*q->v[0]);
#endif
	      ifMass = 1./(mp+mq);
	      for (j=0;j<3;j++) {
		  p->r[j] = ifMass*(mp*p->r[j]+q->fMass*q->r[j]);
		  p->v[j] = ifMass*(mp*p->v[j]+q->fMass*q->v[j]);
		  p->a[j] = ifMass*(mp*p->a[j]+q->fMass*q->a[j]);
		  }
	      p->u += ifMass*(mp*p->u+mq*q->u);
	      p->fMass += mq;
#ifdef STARSINK
	      /* After mass increased:
		 Store Internal angular momentum only as it is invariant under motions  */
	      SINK_Lx(p) -= TRUEMASS(p)*(p->r[1]*p->v[2] - p->r[2]*p->v[1]);
	      SINK_Ly(p) -= TRUEMASS(p)*(p->r[2]*p->v[0] - p->r[0]*p->v[2]);
	      SINK_Lz(p) -= TRUEMASS(p)*(p->r[0]*p->v[1] - p->r[1]*p->v[0]);
#endif
	      assert(q->fMass != 0);
	      assert(!(TYPETest( ((PARTICLE *) q), TYPE_DELETED )));
	      SinkEvent.time= smf->dTime;
	      SinkEvent.mass = TRUEMASS(p);
	      SinkEvent.iOrdSink = p->iOrder;
	      SinkEvent.iOrdVictim = q->iOrder;
	      SinkEvent.r[0] = p->r[0];  SinkEvent.r[1] = p->r[1];  SinkEvent.r[2] = p->r[2];
	      SinkEvent.v[0] = p->v[0];  SinkEvent.v[1] = p->v[1];  SinkEvent.v[2] = p->v[2];
	      SinkEvent.L[0] = SINK_Lx(p); SinkEvent.L[1] = SINK_Ly(p); SinkEvent.L[2] = SINK_Lz(p);
	      SinkLogEntry( &smf->pkd->sinkLog, &SinkEvent, SINK_EVENT_MERGER );
	      q->fMass = 0;
	      printf("Sinks Merged: ID's %d & %d (deleted) were merged\n",p->iOrder,q->iOrder);
	      pkdDeleteParticle(smf->pkd, q);
	      }
	}
#endif
}

#define fSinkRating(_a)    (((PARTICLE *) (_a))->curlv[0] )

void initSinkFormTest(void *p)
{
#ifdef GASOLINE
    fSinkRating(p) = -FLT_MAX;
    iOrderSink(p) = -1;
#ifdef SINKING
    TYPEReset(((PARTICLE *) p),TYPE_NEWSINKING);
#endif
#endif
	}

void combSinkFormTest(void *p1,void *p2)
{
#ifdef GASOLINE
/* Particle p1 belongs to candidate stored in iOrderSink of p1 initially but
   switch to p2's if that candidate is denser */
#ifdef SINKING
    if (TYPETest((PARTICLE *) p1, TYPE_SINKING)) return;
#endif
    if (fSinkRating(p2) > fSinkRating(p1)) {
	fSinkRating(p1) = fSinkRating(p2);
	iOrderSink(p1) = iOrderSink(p2);
	}
#endif
}


void SinkFormTest(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
	int i;
	double dSinkRadius2 = smf->dSinkRadius*smf->dSinkRadius,r2;
	PARTICLE *q;

	/* Apply Bate tests in next phase
	   For now just decide which sink the particle belongs to:
   	          prefer joining a denser (or deeper potential) sink candidate

		  Must be local extremum to be a candidate or abort
	   Need to prevent double counting particles into two sinks
	*/
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2;
	    if (r2 > 0 && r2 <= dSinkRadius2) {
		q = nnList[i].pPart;
		if (TYPETest( q, TYPE_GAS ) && ((smf->bSinkFormPotMin ? q->fPot < p->fPot : q->fDensity > p->fDensity))) {
		    /* Abort without grabbing any particles -- this isn't an extremum particle */
/*		    printf("Sink aborted %d %g: Denser Neighbour %d %g\n",p->iOrder,p->fDensity,q->iOrder,q->fDensity);*/
		    return;

		    }
		}
	    }

	fSinkRating(p) =  (smf->bSinkFormPotMin ? -p->fPot : p->fDensity); /* rate yourself */

/*	printf("Sink %d %g: looking...\n",p->iOrder,p->fDensity);*/

	/* Identify all nbrs as to be eaten unless those nbrs already belong to a "higher rated" sink
	   Thus sink quality is rated by density or potential depth
	*/
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2;
	    if (r2 > 0 && r2 <= dSinkRadius2) {
		q = nnList[i].pPart;
#ifdef SINKING
		if (TYPETest(q, TYPE_SINKING)) continue;
#endif
		if (TYPETest( q, TYPE_GAS ) && q->iRung >= smf->iSinkCurrentRung) {
		    if (smf->bSinkFormPotMin) {
			if (-p->fPot > fSinkRating(q)) {
			    fSinkRating(q) = -p->fPot;
			    iOrderSink(q) = p->iOrder; /* Particle q belongs to sink p */
			    }

			}
		    else {
			if (p->fDensity > fSinkRating(q)) {
			    fSinkRating(q) = p->fDensity;
			    iOrderSink(q) = p->iOrder; /* Particle q belongs to sink p */
			    }
			}
		    }
		}
	    }
#endif
}

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))

void tred3(double a[][3], double d[], double e[])
{
  int l,k,j,i,n=3;
  double scale,hh,h,g,f;

  for (i=n-1;i>=1;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<=l;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i]=a[i][l];
      else {
        for (k=0;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f > 0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=0;j<=l;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=0;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=0;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=0;k<=j;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[0]=0.0;
  e[0]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
      wanted except for statement d[i]=a[i][i]; */
  for (i=0;i<n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=0;j<=l;j++) {
        g=0.0;
        for (k=0;k<=l;k++)
          g += a[i][k]*a[k][j];
        for (k=0;k<=l;k++)
          a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=0;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }
}

double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) {
      double xxx = absb/absa;
      return absa*sqrt(1.0+xxx*xxx);
//      return absa*sqrt(1.0+SQR(absb/absa));
      }
  else if (absb == 0) return 0.;
  else {
      double xxx = absa/absb;
      return absb*sqrt(1.0+xxx*xxx);
//    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
      }
}

void tqli3(double d[], double e[], double z[][3])
{
  double pythag(double a, double b);
  int m,l,iter,i,k,n=3;
  double s,r,p,g,f,dd,c,b;

  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	  if (iter++ == 30) assert(iter<30); /* "Too many iterations in tqli" */
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=pythag(g,1.0);
        g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;
          for (k=0;k<n;k++) {
            f=z[k][i+1];
            z[k][i+1]=s*z[k][i]+c*f;
            z[k][i]=c*z[k][i]-s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
}




void initSinkForm(void *p)
{
	}

void combSinkForm(void *p1,void *p2)
{
    if (!(TYPETest( ((PARTICLE *) p1), TYPE_DELETED )) &&
        TYPETest( ((PARTICLE *) p2), TYPE_DELETED ) ) {
		((PARTICLE *) p1)-> fMass = ((PARTICLE *) p2)-> fMass;
	    pkdDeleteParticle( NULL, p1 );
	    }
    }

void SinkForm(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
	int i,j,nEaten,nEatNow;
	double mtot,im,Ek,Eth,Eg,dvx,dv2,vsink[3];
	double diva,ih2,r2,rs1;
	PARTICLE *q,*q1,*q2;
	SINKEVENT SinkEvent;
	SinkEvent.time = smf->dTime;

	/* You are accreted */
	if ( iOrderSink(p) != -1 ) {
	    /* What if it was temporarily a candidate -- are it's nbrs confused?
	     Shouldn't be -- they should go to a new better candidate if they are near one */

/*	    printf("Sink aborted %d %g: Accreted by other %d\n",p->iOrder,p->fDensity,iOrderSink(p) );*/
	    return;
	    }
#ifdef SINKING
	assert(!TYPETest(p,TYPE_SINKING));
#endif
	iOrderSink(p) = p->iOrder;

	nEaten=0;
	mtot = 0;
	for (i=0;i<nSmooth;++i) {
	    q1 = nnList[i].pPart;
	    if (iOrderSink(q1) == p->iOrder) {
		nEaten++;
		mtot += q1->fMass;
		}
	    }

	if (nEaten < smf->nSinkFormMin) {
//	    printf("Sink aborted %d %g: np %d Mass %g\n",p->iOrder,p->fDensity,nEaten,mtot);
	    return;
	    }

	if (smf->bSinkFormDivAcc) {
	    diva = 0;
	    ih2 = 1./(smf->dSinkRadius*smf->dSinkRadius);
	    for (i=0;i<nSmooth;++i) {
		q1 = nnList[i].pPart;
		if (iOrderSink(q1) == p->iOrder) {
		    double dax,day,daz,dx,dy,dz;

		    r2 = nnList[i].fDist2*ih2;  /* Calculate div a (un-normalized) */
		    if (r2 < 4) {
			DKERNEL(rs1,r2);
			rs1 *= q1->fMass;
			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dax = (-p->a[0] + q1->a[0]);
			day = (-p->a[1] + q1->a[1]);
			daz = (-p->a[2] + q1->a[2]);
			diva += (dax*dx+day*dy+daz*dz)*rs1;
			}
		    }
		}
	    if (diva > 0) return;
	    }


	if (smf->bSinkFormDV || smf->bSinkFormDivV) {
	    double dvdx[3][3];
	    double d[3],e[3];
	    double vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	    dvdx[0][0]=0; dvdx[0][1]=0; dvdx[0][2]=0;
	    dvdx[1][0]=0; dvdx[1][1]=0; dvdx[1][2]=0;
	    dvdx[2][0]=0; dvdx[2][1]=0; dvdx[2][2]=0;

	    ih2 = 1./(smf->dSinkRadius*smf->dSinkRadius);
	    for (i=0;i<nSmooth;++i) {
		q1 = nnList[i].pPart;
		if (iOrderSink(q1) == p->iOrder) {
		    double dvx,dvy,dvz,dx,dy,dz;

		    r2 = nnList[i].fDist2*ih2;  /* Calculate dv_i/dx_j (un-normalized) */
		    if (r2 < 4) {
			DKERNEL(rs1,r2);
			rs1 *= q1->fMass;
			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = ((-p->vPred[0] + q1->vPred[0])*vFac - dx*smf->H)*rs1; /* NB: dx = px - qx */
			dvy = ((-p->vPred[1] + q1->vPred[1])*vFac - dy*smf->H)*rs1;
			dvz = ((-p->vPred[2] + q1->vPred[2])*vFac - dz*smf->H)*rs1;
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += -1.5*smf->PP.dOrbFreq*smf->PP.dWidth*vFac*rs1;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= -1.5*smf->PP.dOrbFreq*smf->PP.dWidth*vFac*rs1;
#endif
			dvdx[0][0] += dvx*dx;
			dvdx[0][1] += dvx*dy;
			dvdx[0][2] += dvx*dz;
			dvdx[1][0] += dvy*dx;
			dvdx[1][1] += dvy*dy;
			dvdx[1][2] += dvy*dz;
			dvdx[2][0] += dvz*dx;
			dvdx[2][1] += dvz*dy;
			dvdx[2][2] += dvz*dz;
			}
		    }
		}

	    dvdx[0][1] = 0.5*(dvdx[0][1]+dvdx[1][0]); dvdx[1][0]=dvdx[0][1];
	    dvdx[2][1] = 0.5*(dvdx[2][1]+dvdx[1][2]); dvdx[1][2]=dvdx[2][1];
	    dvdx[2][0] = 0.5*(dvdx[2][0]+dvdx[0][2]); dvdx[0][2]=dvdx[2][0];
//	    if (dvdx[0][0] != 0) {
//		printf("Sink dvxd %d: %g %g %g\n",p->iOrder,dvdx[0][0],dvdx[0][1],dvdx[0][2]);
//		printf("Sink dvyd %d: %g %g %g\n",p->iOrder,dvdx[1][0],dvdx[1][1],dvdx[1][2]);
//		printf("Sink dvzd %d: %g %g %g   nEaten %d\n",p->iOrder,dvdx[2][0],dvdx[2][1],dvdx[2][2],nEaten);
//		}

	    if (smf->bSinkFormDivV) {

		if (dvdx[0][0]+dvdx[1][1]+dvdx[2][2] > 0) {
//		    printf("Sink not converging %d %g: divv %g\n",p->iOrder,p->fDensity,dvdx[0][0]+dvdx[1][1]+dvdx[2][2]);
		    return; /* not converging */
		    }
		}
	    if (smf->bSinkFormDV) {
		tred3(dvdx,d,e);
		tqli3(d,e,dvdx);

		if (d[0] > 0 || d[1] > 0 || d[2] > 0) {
//		    printf("Sink not converging in 3ds %d %g: divv %g %g %g\n",p->iOrder,p->fDensity,d[0],d[1],d[2]);
		    return; /* not converging in all directions */
		    }
		}
	    }

	Ek = 0;
	Eth = 0;
	Eg = 0;
	vsink[0] = 0; 	vsink[1] = 0;	vsink[2] = 0;

	for (i=0;i<nSmooth;++i) {
	    q1 = nnList[i].pPart;
/*	    printf("Sink %d %g: trying to use %d %d, rungs %d %d, r %g %g\n",p->iOrder,p->fDensity,q1->iOrder,iOrderSink(q1),q1->iRung,smf->iSinkCurrentRung,nnList[i].fDist2,smf->dSinkRadius);*/
	    if (iOrderSink(q1) == p->iOrder) {
		dvx = q1->v[0]-p->v[0];
		vsink[0] += q1->fMass*dvx;
		dv2 = dvx*dvx;
		dvx = q1->v[1]-p->v[1];
		vsink[1] += q1->fMass*dvx;
		dv2 += dvx*dvx;
		dvx = q1->v[2]-p->v[2];
		vsink[2] += q1->fMass*dvx;
		dv2 += dvx*dvx;
		Ek += 0.5*q1->fMass*dv2;
		Eth += q1->fMass*q1->u;
		for (j=i+1;j<nSmooth;j++) {
		    q2 = nnList[j].pPart;
		    if (iOrderSink(q2) == p->iOrder) {
			dvx = q1->r[0]-q2->r[0];
			dv2 = dvx*dvx;
			dvx = q1->r[1]-q2->r[1];
			dv2 += dvx*dvx;
			dvx = q1->r[2]-q2->r[2];
			dv2 += dvx*dvx;
			Eg -= q1->fMass*q2->fMass/sqrt(dv2);
			}
		    }
		}
	    }

/*	printf("Sink %d Corrected Ek %g %g\n",p->iOrder,Ek,Ek - 0.5*(vsink[0]*vsink[0]+vsink[1]*vsink[1]+vsink[2]*vsink[2])/mtot);*/
	Ek -= 0.5*(vsink[0]*vsink[0]+vsink[1]*vsink[1]+vsink[2]*vsink[2])/mtot;

	/* Apply Bate tests here --
	   1. thermal energy < 1/2 |Grav|,
	   2. thermal + rot E < |Grav|,
	   3. total E < 0 (implies 2.)
	   4. div.acc < 0 (related to rate of change of total E I guess)  (optional see above)
           5. NEW eigenvalues of dvi/dxj all negative (Federrath) (optional see above)
           6. NEW div.v < 0 (Whitworth group)
	*/

	if (Eth < 0.5*fabs(Eg) && Ek + Eth + Eg < 0) {
	    /* Sink approved */
	    if (smf->dSinkTimeEligible > 0) {
		/* Check that formation conditions persisted for minimum time */
		/* This cannot work in its simplest form because the sink candidate
                   for a given density peak changes randomly from step to step as
                   densities etc.. fluctuate.  Need to code in some way to continuously
                   refer to the same pre-sink blob (inherit properties of particles inside 1/2 SinkRadius?) */
		if (p->fTimeForm > 1e36) {
		    p->fTimeForm = smf->dTime;
		    return;
		    }
		if (smf->dTime-p->fTimeForm < smf->dSinkTimeEligible) return;
		}

	    PARTICLE sinkp = *p;
	    sinkp.r[0] = 0;
	    sinkp.r[1] = 0;
	    sinkp.r[2] = 0;
	    sinkp.v[0] = 0;
	    sinkp.v[1] = 0;
	    sinkp.v[2] = 0;
	    sinkp.a[0] = 0;
	    sinkp.a[1] = 0;
	    sinkp.a[2] = 0;
	    sinkp.u = 0;
#ifdef STARSINK
	    SINK_Lx(&sinkp) = 0;
	    SINK_Ly(&sinkp) = 0;
	    SINK_Lz(&sinkp) = 0;
#endif
	    sinkp.fMass = 0;
#ifdef SINKING
	    sinkp.fTrueMass = 0;
	    sinkp.iSinkingOnto = 0; /* Must start at zero -- it records the number of assoc. sinking particles */
#endif
	    nEatNow = 0;
	    for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2;
		if (iOrderSink(q) == p->iOrder) {
#ifdef SINKING
		    if (r2 > smf->dSinkMustAccreteRadius*smf->dSinkMustAccreteRadius) {
			/* Make sinking -- can't force sink choice b/c sink iOrder isn't known yet */
			TYPESet(q, TYPE_NEWSINKING);
			}
		    else
#endif /*SINKING*/
			{
			FLOAT mq,rx,ry,rz,vx,vy,vz;
			mq = q->fMass;
			rx = q->r[0]; ry = q->r[1]; rz = q->r[2];
			vx = q->v[0]; vy = q->v[1]; vz = q->v[2];
			nEatNow++;
			sinkp.r[0] += mq*rx;
			sinkp.r[1] += mq*ry;
			sinkp.r[2] += mq*rz;
			sinkp.v[0] += mq*vx;
			sinkp.v[1] += mq*vy;
			sinkp.v[2] += mq*vz;
			sinkp.a[0] += mq*q->a[0];
			sinkp.a[1] += mq*q->a[1];
			sinkp.a[2] += mq*q->a[2];
			sinkp.u += q->fMass*q->u;
#ifdef STARSINK
			SINK_Lx(&sinkp) += mq*(ry*vz - rz*vy); /* Add to total L */
			SINK_Ly(&sinkp) += mq*(rz*vx - rx*vz);
			SINK_Lz(&sinkp) += mq*(rx*vy - ry*vx);
#endif
			sinkp.fMass += mq;
#ifdef SINKING
			sinkp.fTrueMass += mq;
#endif
			SinkEvent.iOrdVictim = q->iOrder;
			SinkLogEntry( &smf->pkd->sinkLog, &SinkEvent, SINK_EVENT_ACCRETE_AT_FORMATION );
			if (p!=q) {
			    q->fMass = 0;
			    pkdDeleteParticle(smf->pkd, q);
			    }
			}
		    }
		}

	    im = 1/sinkp.fMass;
	    sinkp.r[0] *= im;
	    sinkp.r[1] *= im;
	    sinkp.r[2] *= im;
	    sinkp.v[0] *= im;
	    sinkp.v[1] *= im;
	    sinkp.v[2] *= im;
	    sinkp.a[0] *= im;
	    sinkp.a[1] *= im;
	    sinkp.a[2] *= im;
	    sinkp.u *= im;
#ifdef STARSINK
            /* Store Internal angular momentum only */
	    SINK_Lx(&sinkp) -= sinkp.fMass*(sinkp.r[1]*sinkp.v[2] - sinkp.r[2]*sinkp.v[1]);
	    SINK_Ly(&sinkp) -= sinkp.fMass*(sinkp.r[2]*sinkp.v[0] - sinkp.r[0]*sinkp.v[2]);
	    SINK_Lz(&sinkp) -= sinkp.fMass*(sinkp.r[0]*sinkp.v[1] - sinkp.r[1]*sinkp.v[0]);
#endif
	    TYPEClear(&sinkp);
	    TYPESet(&sinkp,TYPE_SINK|TYPE_STAR|TYPE_ACTIVE);
	    sinkp.fTimeForm = -smf->dTime; /* -ve time is sink indicator */
#ifdef SINKEXTRADATA
	    for(j = 0; j < 3; j++) {
		sinkp.rForm[j] = sinkp.r[j];
		sinkp.vForm[j] = sinkp.v[j];
		}
#endif
	    printf("Sink Formed %d %g: np (%d) %d Mass (%g) %g Ek %g Eth %g Eg %g, %g %g\n",p->iOrder,p->fDensity,nEaten,nEatNow,mtot,sinkp.fMass,Ek,Eth,Eg,4/3.*pow(M_PI,2.5)/50.*sqrt((p->c*p->c*p->c*p->c*p->c*p->c)/(p->fMass*p->fMass*p->fDensity)),pow(Eth/fabs(Eg),1.5) );
#ifndef SINKING
	    assert(fabs(sinkp.fMass/mtot-1) < 1e-4);
#else
/*	    printf("Sink Made %g %g\n",sinkp.fMass,sinkp.fTrueMass);*/
	    assert(fabs(sinkp.fMass/sinkp.fTrueMass-1) < 1e-7);
	    assert(!TYPETest(p, TYPE_SINKING));
#endif
	    p->fMass = 0;
	    pkdDeleteParticle(smf->pkd, p);
	    SinkEvent.mass = TRUEMASS(&sinkp);
	    SinkEvent.r[0] = sinkp.r[0]; SinkEvent.r[1] = sinkp.r[1]; SinkEvent.r[2] = sinkp.r[2];
	    SinkEvent.v[0] = sinkp.v[0]; SinkEvent.v[1] = sinkp.v[1]; SinkEvent.v[2] = sinkp.v[2];
	    SinkEvent.L[0] = SINK_Lx(&sinkp); SinkEvent.L[1] = SINK_Ly(&sinkp); SinkEvent.L[2] = SINK_Lz(&sinkp);
	    SinkLogEntry( &smf->pkd->sinkLog, &SinkEvent, SINK_EVENT_FORM );
	    pkdNewParticle(smf->pkd, sinkp);  /* Sets iOrder to -1 */
	    }
	else {
	    printf("Sink Failed E Tests %d %g: np (%d) Mass (%g) Ek %g Eth %g Eg %g, %g %g\n",p->iOrder,p->fDensity,nEaten,mtot,Ek,Eth,Eg,	    4/3.*pow(M_PI,2.5)/50.*sqrt((p->c*p->c*p->c*p->c*p->c*p->c)/(p->fMass*p->fMass*p->fDensity)),pow(Eth/fabs(Eg),1.5) );
	    }
#endif
}


#ifdef SUPERCOOL
void initMeanVel(void *p)
{
	int j;

	for (j=0;j<3;++j) {
		((PARTICLE *)p)->vMean[j] = 0.0;
		}
	}

void combMeanVel(void *p1,void *p2)
{
	int j;

	for (j=0;j<3;++j) {
		((PARTICLE *)p1)->vMean[j] += ((PARTICLE *)p2)->vMean[j];
		}
	}

void MeanVel(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i,j;

	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (2./3.)*sqrt(ih2);
#else
	fNorm = M_1_PI*sqrt(ih2)*ih2;
#endif
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rs *= fNorm;
		q = nnList[i].pPart;
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
        double dvy = 0;
        if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
            dvy += -1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
        else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
            dvy -= -1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
        p->vMean[0] += rs*q->fMass/q->fDensity*q->v[0];
        p->vMean[1] += rs*q->fMass/q->fDensity*(q->v[1]+dvy);
        p->vMean[2] += rs*q->fMass/q->fDensity*q->v[2];
#else
		for (j=0;j<3;++j) {
			p->vMean[j] += rs*q->fMass/q->fDensity*q->v[j];
			}
#endif
		}
	}

void MeanVelSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i,j;

	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (1./3.)*sqrt(ih2);
#else
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
#endif
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rs *= fNorm;
		q = nnList[i].pPart;
#ifdef SLIDING_PATCH
#error Cannot do mean velocity with SLIDING_PATCH
#endif
		for (j=0;j<3;++j) {
			p->vMean[j] += rs*q->fMass/q->fDensity*q->v[j];
			q->vMean[j] += rs*p->fMass/p->fDensity*p->v[j];
			}
		}
	}
#endif /* SUPER_COOL */


#ifdef GASOLINE
/* Original Particle */
void initSphPressureTermsParticle(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
        ((PARTICLE *)p)->dtNew = FLT_MAX;
        ((PARTICLE *)p)->uDotDiff = 0.0;
		((PARTICLE *)p)->uDotPdV = 0.0;
        ((PARTICLE *)p)->uDotAV = 0.0;
        ((PARTICLE *)p)->uDotBdiss = 0.0;
        ((PARTICLE *)p)->uDotAniCond = 0.0;
#ifdef UNONCOOL
        ((PARTICLE *)p)->uHotDotDiff = 0.0;
#endif
#ifdef MHD
	int i;
	for(i=0;i<3;++i){
      ((PARTICLE *)p)->BDotRho[i] = 0.0;
      ((PARTICLE *)p)->BDotDiss[i] = 0.0;
      ((PARTICLE *)p)->BDotBClean[i] = 0.0;
      ((PARTICLE *)p)->curlB[i] = 0.0;
}
	((PARTICLE *)p)->BCleanDot = 0.0;
	((PARTICLE *)p)->divB = 0.0;
#endif
#ifdef SETPRANDTLVISC
	int j;
	for(j=0;j<3;++j){
	  ((PARTICLE *)p)->AccAV[j] = 0.0;
	}
#endif
#ifdef DRHODT
		((PARTICLE *)p)->fDivv_PdV = 0.0;
		((PARTICLE *)p)->fDivv_PdVcorr = 0.0;
#endif
#ifdef DIFFUSION
		((PARTICLE *)p)->fMetalsDot = 0.0;
#ifdef STARFORM
		((PARTICLE *)p)->fMFracOxygenDot = 0.0;
		((PARTICLE *)p)->fMFracIronDot = 0.0;
#endif /* STARFORM */
#ifdef MASSDIFF
		((PARTICLE *)p)->fMassDot = 0.0;
#endif
#endif /* DIFFUSION */
		}
	}

/* Cached copies of particle */
void initSphPressureTerms(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->dtNew = FLT_MAX;
        ((PARTICLE *)p)->uDotDiff = 0.0;
        ((PARTICLE *)p)->uDotPdV = 0.0;
		((PARTICLE *)p)->uDotAV = 0.0;
    ((PARTICLE *)p)->uDotBdiss = 0.0;
    ((PARTICLE *)p)->uDotAniCond = 0.0;
#ifdef UNONCOOL
        ((PARTICLE *)p)->uHotDotDiff = 0.0;
#endif
#ifdef MHD
	int i;
	for(i=0;i<3;++i){
      ((PARTICLE *)p)->BDotRho[i] = 0.0;
      ((PARTICLE *)p)->BDotDiss[i] = 0.0;
      ((PARTICLE *)p)->BDotBClean[i] = 0.0;
      ((PARTICLE *)p)->curlB[i] = 0.0;
}
	((PARTICLE *)p)->BCleanDot = 0.0;
	((PARTICLE *)p)->divB = 0.0;
#endif
#ifdef SETPRANDTLVISC
	int j;
	for(j=0;j<3;++j){
          ((PARTICLE *)p)->AccAV[j] = 0.0;
	}
#endif
#ifdef DRHODT
		((PARTICLE *)p)->fDivv_PdV = 0.0;
		((PARTICLE *)p)->fDivv_PdVcorr = 0.0;
#endif
		ACCEL(p,0) = 0.0;
		ACCEL(p,1) = 0.0;
		ACCEL(p,2) = 0.0;
#ifdef DIFFUSION
		((PARTICLE *)p)->fMetalsDot = 0.0;
#ifdef STARFORM
		((PARTICLE *)p)->fMFracOxygenDot = 0.0;
		((PARTICLE *)p)->fMFracIronDot = 0.0;
#endif /* STARFORM */
#ifdef MASSDIFF
		((PARTICLE *)p)->fMassDot = 0.0;
#endif
#endif /* DIFFUSION */
		}
	}

void combSphPressureTerms(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
        ((PARTICLE *)p1)->uDotDiff += ((PARTICLE *)p2)->uDotDiff;
        ((PARTICLE *)p1)->uDotPdV += ((PARTICLE *)p2)->uDotPdV;
        ((PARTICLE *)p1)->uDotAV += ((PARTICLE *)p2)->uDotAV;
          ((PARTICLE *)p1)->uDotBdiss += ((PARTICLE *)p2)->uDotBdiss;
          ((PARTICLE *)p1)->uDotAniCond += ((PARTICLE *)p2)->uDotAniCond;
#ifdef UNONCOOL
        ((PARTICLE *)p1)->uHotDotDiff += ((PARTICLE *)p2)->uHotDotDiff;
#endif
#ifdef MHD
        int i;
    for(i=0;i<3;++i){
	  ((PARTICLE *)p1)->BDotRho[i] += ((PARTICLE *)p2)->BDotRho[i];
	  ((PARTICLE *)p1)->BDotDiss[i] += ((PARTICLE *)p2)->BDotDiss[i];
	  ((PARTICLE *)p1)->BDotBClean[i] += ((PARTICLE *)p2)->BDotBClean[i];
	  ((PARTICLE *)p1)->curlB[i] += ((PARTICLE *)p2)->curlB[i];
        }
	((PARTICLE *)p1)->BCleanDot += ((PARTICLE *)p2)->BCleanDot;
	((PARTICLE *)p1)->divB += ((PARTICLE *)p2)->divB;
#endif
#ifdef SETPRANDTLVISC
	int j;
	for(j=0;j<3;++j){
          ((PARTICLE *)p1)->AccAV[j] += ((PARTICLE *)p2)->AccAV[j];
	}
#endif
#ifdef DRHODT
 	        ((PARTICLE *)p1)->fDivv_PdV += ((PARTICLE *)p2)->fDivv_PdV;
 	        ((PARTICLE *)p1)->fDivv_PdVcorr += ((PARTICLE *)p2)->fDivv_PdVcorr;
#endif
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL(p1,0) += ACCEL(p2,0);
		ACCEL(p1,1) += ACCEL(p2,1);
		ACCEL(p1,2) += ACCEL(p2,2);
#ifdef DIFFUSION
		((PARTICLE *)p1)->fMetalsDot += ((PARTICLE *)p2)->fMetalsDot;
#ifdef STARFORM
		((PARTICLE *)p1)->fMFracOxygenDot += ((PARTICLE *)p2)->fMFracOxygenDot;
		((PARTICLE *)p1)->fMFracIronDot += ((PARTICLE *)p2)->fMFracIronDot;
#endif /* STARFORM */
#ifdef MASSDIFF
		((PARTICLE *)p1)->fMassDot += ((PARTICLE *)p2)->fMassDot;
#endif
#endif /* DIFFUSION */
		}
		if (((PARTICLE *)p2)->dtNew < ((PARTICLE *)p1)->dtNew)
			((PARTICLE *)p1)->dtNew = ((PARTICLE *)p2)->dtNew;
	}

/* Gather only version -- never use */
void SphPressureTerms(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
        assert(0);
	}

void SphPressureTermsSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rs,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr,dvdotdrI;
	FLOAT pPoverRho2,pPoverRho2f,pMass;
	FLOAT qPoverRho2,qPoverRho2f;
	FLOAT ph,pc,dt,visc,absmu;
	FLOAT fNorm,fNorm1,aFac,vFac,divvi,divvj;
	FLOAT dxterm,dyterm,dzterm,drterm2=0.0;
	double rpartp,rpartq,rho0q,rho0p;
#ifdef USEVEC
  FLOAT Accp[3],Accq[3];
#else
  FLOAT Accp,Accq,Accpx,Accpy,Accpz,Accqx,Accqy,Accqz;
  Accpx=0.0;
  Accpy=0.0;
  Accpz=0.0;
  Accqx=0.0;
  Accqy=0.0;
  Accqz=0.0;
#endif
  FLOAT vSigp,vSigq,pPeff,qPeff;
  FLOAT pP,qP,igDensity2p,igDensity2q,rqoverRhoq,rpoverRhop;
  double pDensity = p->fDensity;
	int i;
#ifdef STARCLUSTERFORM
    int iDenMax=1;
#endif
	pc = p->c;
	pMass = p->fMass;
  vSigp = pc;
#ifdef PEXT
#define PEXTCORR (-smf->Pext)
#else
#define PEXTCORR 0
#endif
#ifndef GDFORCE
	pPoverRho2 = p->PoverRho2;
  { FLOAT pd2 = p->fDensity*p->fDensity;
    pPoverRho2f = (pPoverRho2*pd2+PEXTCORR)/pd2;
    }
#endif /* ndef GDFORCE */
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (1./3.)/ph;
#else
	fNorm = 0.5*M_1_PI*ih2/ph;
#endif
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);        /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */
	if (p->iOrder==0){printf(" cs %f P %f nsmooth %d ph %f u %f upred %f dens %f", p->PoverRho2*p->fDensity,p->PoverRho2*p->fDensity*p->fDensity,nSmooth,ph,p->u,p->uPred,p->fDensity);}

#ifdef MHD
	/* MHD properties regarding particle p*/
  double pB2 = p->BPred[0]*p->BPred[0]+p->BPred[1]*p->BPred[1]+p->BPred[2]*p->BPred[2];
  vSigp = sqrt(pB2/(pDensity)+pc*pc);
if (p->iOrder==0){printf(" B2 %f bx %f by %f bz %f vSigp %f", pB2,p->BPred[0],p->BPred[1],p->BPred[2],vSigp);}
#ifdef MHDCLEANING
  double chp;
  chp = smf->fch*vSigp;
	if (TYPEQueryACTIVE(p)) {
        p->BCleanDot -= p->BCleanPred*smf->psidecfac*chp/ph;
	}
#endif
#endif
//#ifdef DRHODT
	divvi = 0;
	divvj = 0;
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2*ih2;
	    q = nnList[i].pPart;
	    #ifdef ISPH
	    KERNEL(rs,r2);
            dx = nnList[i].dx;
            dy = nnList[i].dy;
            dz = nnList[i].dz;
            dxterm = -(p->Igrxx*dx+p->Igrxy*dy+p->Igrxz*dz);
            dyterm = -(p->Igrxy*dx+p->Igryy*dy+p->Igryz*dz);
            dzterm = -(p->Igrxz*dx+p->Igryz*dy+p->Igrzz*dz);
	    rs *= (dx*dxterm+dy*dyterm+dz*dzterm)*q->fMass;
	    divvi += rs/p->fDensity;
	    divvj += rs/q->fDensity;
	    #else
	    DKERNEL(rs1,r2);
	    rs1 *= nnList[i].fDist2*q->fMass;
	    divvi += rs1/p->fDensity;
	    divvj += rs1/q->fDensity;
	    #endif
#ifdef STARCLUSTERFORM
        iDenMax &= (q->fDensity <= p->fDensity);
#endif
	    }
#ifdef STARCLUSTERFORM
    if (iDenMax) TYPESet(p,TYPE_DENMAX);
    else TYPEReset(p,TYPE_DENMAX);
#endif

#ifdef GDFORCE
/* The DIVVCORRBAD corrector is better on average but is pathological with very
   uneven particle distributions (very large density gradients) */
#ifdef DIVVCORRBAD
    p->fDivv_Corrector = (divvbad != 0 ? -(3/2.)/(divvbad*fNorm1) : 1); /* fNorm1 normalization includes 0.5 */
#else
    p->fDivv_Corrector = (divvj != 0 ? divvi/divvj : 1); /* GDFORCE CORR */
#endif
#else
#ifdef DIVVCORRBAD
    p->fDivv_Corrector = (divvi != 0 ? -(3/2.)/(divvi*fNorm1) : 1); /* fNorm1 normalization includes 0.5 */
#else
    p->fDivv_Corrector = 1;
#endif
#endif
//#endif

	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    if (!TYPEQueryACTIVE(p) && !TYPEQueryACTIVE(q)) continue;
	    r2 = nnList[i].fDist2*ih2;
      dx = nnList[i].dx;
      dy = nnList[i].dy;
      dz = nnList[i].dz;

      #ifdef ISPH
      KERNEL(rs,r2);
      DKERNEL(rs1,r2); /* rs1 is negative */
      rs *= fNorm;
      //rs *= p->fDivv_Corrector;
      rp = rs * pMass;
      rq = rs * q->fMass;
      dxterm = -(p->Igrxx*dx+p->Igrxy*dy+p->Igrxz*dz);
      dyterm = -(p->Igrxy*dx+p->Igryy*dy+p->Igryz*dz);
      dzterm = -(p->Igrxz*dx+p->Igryz*dy+p->Igrzz*dz);
      #else
      DKERNEL(rs1,r2);
      rs1 *= fNorm1;
      rs1 *= p->fDivv_Corrector;
      rp = rs1 * pMass;
      rq = rs1 * q->fMass;
      dxterm = dx;
      dyterm = dy;
      dzterm = dz;
      #endif
      //if(p->iOrder==0) printf("q->iOrder %d rq %f rp %f dxterm %f dyterm %f dzterm %f pDensity %f p->fDensity %f p->fDivv_Corrector %f \n",q->iOrder,rq,rp,dxterm,dyterm,dzterm,pMass,p->fMass,p->fDivv_Corrector);
	    dvx = p->vPred[0] - q->vPred[0];
	    dvy = p->vPred[1] - q->vPred[1];
	    dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
	    if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
			dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
	    else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
	        dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
	    dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
		+ nnList[i].fDist2*smf->H;
   dvdotdrI = vFac*(dvx*dxterm + dvy*dyterm + dvz*dzterm);
   pP = p->PoverRho2*pDensity*pDensity;
   pPeff = (pP+PEXTCORR) >= 0.0 ? (pP+PEXTCORR) : 0.0;
   qP = q->PoverRho2*q->fDensity*q->fDensity;
   qPeff = (qP+PEXTCORR) >= 0.0 ? (qP+PEXTCORR) : 0.0;
   if(p->iOrder != q->iOrder){drterm2 = (dx*dxterm+dy*dyterm+dz*dzterm)/nnList[i].fDist2; }
#ifdef GDFORCE
        {
       double igDensity2 = 1/(pDensity*q->fDensity);
      igDensity2p = igDensity2*0.6667*p->u*p->fDensity;
      igDensity2q = igDensity2*0.6667*q->u*q->fDensity;;
      rqoverRhoq = rq/q->fDensity;
      rpoverRhop = rp/p->fDensity;
      pPoverRho2f = (pPeff)*igDensity2;
      qPoverRho2f = (qPeff)*igDensity2;
      pPoverRho2 = pP*igDensity2;
      qPoverRho2 = qP*igDensity2;
      }
#else /* now not GDFORCE */
      igDensity2p = 1.0/(p->fDensity*p->fDensity);
      igDensity2q = 1.0/(q->fDensity*q->fDensity);
      rqoverRhoq = rq/p->fDensity;
      rpoverRhop = rp/q->fDensity;
      pPoverRho2 = p->PoverRho2;
      qPoverRho2 = q->PoverRho2;
      pPoverRho2f = (pPeff)*igDensity2p;
      qPoverRho2f = (qPeff)*igDensity2q;
#endif /* GDFORCE */
	vSigq = q->c;
#ifdef THERMALCONDANI
	double gradTdotdrp = p->gradTx*dxterm+p->gradTy*dyterm+p->gradTz*dzterm;
	double gradTdotdrq = q->gradTx*dxterm+q->gradTy*dyterm+q->gradTz*dzterm;
	//if (p->iOrder == 1196) printf("Order %d udot %g gradTpx %g gradTpy %g gradTpz %g gradTpdotdr %g",p->iOrder,p->uDotAniCond,p->gradTx,p->gradTy,p->gradTz,gradTdotdrp);
#endif //THERMALCONDANI
#ifdef MHD
	/*fch and mby tauch should be set in parameter file*/
  double qB2,chq,dBx,dBy,dBz,dB2,Bdotdrp,Bdotdrq,dBdotdr,Bcritp,Bcritq,Bfacp,Bfacq;
	qB2 = q->BPred[0]*q->BPred[0]+q->BPred[1]*q->BPred[1]+q->BPred[2]*q->BPred[2];
	vSigq = sqrt(qB2/q->fDensity+q->c*q->c);
#ifdef MHDCLEANING
chq = smf->fch*vSigq;
#endif
#ifdef GLASSMHDBETA
 double betaglassp = 1.0;
 double betaglassq = 1.0;
 int j;
 if (pB2!=0.0){betaglassp=sqrt((pPeff)/(0.5*pB2)/25.0);}
 if (qB2!=0.0){betaglassq=sqrt((qPeff)/(0.5*qB2)/25.0);}
 if (pPeff!=0.0 && pB2==0.0){p->B[1]=p->BPred[1]=sqrt(pPeff/25.0*0.5);}
 if (qPeff!=0.0 && qB2==0.0){q->B[1]=q->BPred[1]=sqrt(qPeff/25.0*0.5);}
 if(betaglassp<0.01 || betaglassq<0.01){printf("(bfacp: %f bfacq: %f pPeff %f qPeff %f densp %f densq %f iq %d \n",betaglassp,betaglassq,pPeff,qPeff,p->fDensity,q->fDensity,q->iOrder);}
 if(p->iOrder==408){printf("(Before Byp: %f Byq:%f ip %d iq %d \n",p->BPred[1],q->BPred[1],p->iOrder,q->iOrder);}
 for (j=0;j<3;++j) {
   p->BPred[j]=betaglassp*p->BPred[j];
   q->BPred[j]=betaglassq*q->BPred[j];
   p->B[j]=betaglassp*p->BPred[j];
   q->B[j]=betaglassq*q->BPred[j];
   }
 qB2 = q->BPred[0]*q->BPred[0]+q->BPred[1]*q->BPred[1]+q->BPred[2]*q->BPred[2];
 vSigq = sqrt(qB2/q->fDensity+q->c*q->c);
 pB2 = p->BPred[0]*p->BPred[0]+p->BPred[1]*p->BPred[1]+p->BPred[2]*p->BPred[2];
 vSigp = sqrt(pB2/(pDensity)+pc*pc);
#endif
        dBx = p->BPred[0]-q->BPred[0];
        dBy = p->BPred[1]-q->BPred[1];
	dBz = p->BPred[2]-q->BPred[2];
	dB2 = dBx*dBx+dBy*dBy+dBz*dBz;
	Bdotdrp= p->BPred[0]*dxterm+p->BPred[1]*dyterm+p->BPred[2]*dzterm;
	Bdotdrq= q->BPred[0]*dxterm+q->BPred[1]*dyterm+q->BPred[2]*dzterm;
	dBdotdr=dBx*dxterm+dBy*dyterm+dBz*dzterm;
	/* Tensile correction term */
	if (pB2!=0.0){
	  Bcritp=(pPeff)/(0.5*pB2);
	if (Bcritp<1.0) {Bfacp = 1.0;}
        else if(Bcritp<2.0){Bfacp = (2.0-Bcritp);}
        else{Bfacp = 0.0;}}
	else{Bfacp = 0.0;}
	if (qB2!=0.0){
	  Bcritq=(qPeff)/(0.5*qB2);
  if (Bcritq<1.0) {Bfacq = 1.0;}
  else if(Bcritq<2.0){Bfacq = (2.0-Bcritq);}
        else{Bfacq = 0.0;}}
	else{Bfacq = 0.0;}

#ifdef MHDRESISTIVITY
	double dvXdx,dvXdy,dvXdz,vSigBp,vSigBq,alphaBp,alphaBq,drterm=0.0;
/* Signal speed for artificial resistivity */
if(p->iOrder != q->iOrder){
  drterm = (dx*dxterm+dy*dyterm+dz*dzterm)/sqrt(nnList[i].fDist2);
      dvXdx=(-dvz*dy + dvy*dz);
      dvXdy=(-dvx*dz + dvz*dx);
      dvXdz=(-dvy*dx + dvx*dy);
      vSigBp = sqrt(dvXdx*dvXdx+dvXdy*dvXdy+dvXdz*dvXdz)/sqrt(nnList[i].fDist2);
      vSigBq = vSigBp;
      alphaBp=smf->alphaB;
      alphaBq=smf->alphaB;
    }
#endif
#ifdef THERMALCONDANI
double gradTdotBp = (p->BPred[0]*p->gradTx+p->BPred[1]*p->gradTy+p->BPred[2]*p->gradTz);
double gradTdotBq = (q->BPred[0]*q->gradTx+q->BPred[1]*q->gradTy+q->BPred[2]*q->gradTz);
double BgradTdotBdotdrp = (pB2 != 0 ? gradTdotBp*Bdotdrp/pB2 : 0.0);
double BgradTdotBdotdrq = (qB2 != 0 ? gradTdotBq*Bdotdrq/qB2 : 0.0);
//if (p->iOrder == 1196) {printf("Qiorder %d BdotgradTp %g BdotgradTq %g Bgradp %g Bgradq %g",q->iOrder,gradTdotBp,gradTdotBq,BgradTdotBdotdrp,BgradTdotBdotdrq);
// if (BgradTdotBdotdrq != BgradTdotBdotdrq) printf("bx %g,by %g,bz %g ,gtx %g,gty %g,gtz %g",q->BPred[0],q->BPred[1],q->BPred[2],q->gradTx,q->gradTy,q->gradTz);}
#endif //THERMALCONDANI
#endif
	    if (TYPEQueryACTIVE(p)) {
		if (TYPEQueryACTIVE(q)) {
#define PACTIVE(xxx) xxx
#define QACTIVE(xxx) xxx
#include "SphPressureTerms.h"
		    }
		else {
#undef QACTIVE
#define QACTIVE(xxx)
#include "SphPressureTerms.h"
		    }
		}
	    else if (TYPEQueryACTIVE(q)) {
#undef PACTIVE
#define PACTIVE(xxx)
#undef QACTIVE
#define QACTIVE(xxx) xxx
#include "SphPressureTerms.h"
}
	    }
  }

/* NB: ACCEL_PRES used here --
   with shock tracking disabled: #define NOSHOCKTRACK
   it is: a->aPres
   otherwise it is identical to p->a
   The postSphPressure function combines p->a and p->aPres
*/

#if (0)
#define DEBUGFORCE( xxx )  if (p->iOrder == 0 || q->iOrder == 0) { \
        if (p->iOrder == 0) printf("%s  %d-%d %g %g\n",xxx,p->iOrder,q->iOrder,rq,pa[0]+p->a[0]); \
        else printf("%s  %d-%d %g %g\n",xxx,p->iOrder,q->iOrder,rp,q->a[0]); }
#else
#define DEBUGFORCE( xxx )
#endif

void SphPressureTermsSymOld(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pPoverRho2f,pa[3],pMass,pmumax;
	FLOAT qPoverRho2,qPoverRho2f;
	FLOAT ph,pc,pDensity,visc,hav,absmu;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;



#ifdef DEBUGFORCE
	pa[0]=0.0;
#endif

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	pPoverRho2 = p->PoverRho2;
#ifdef PEXT
    {
        FLOAT pd2 = p->fDensity*p->fDensity;
	pPoverRho2f = (pPoverRho2*pd2-smf->Pext)/pd2;
    }
#else
	pPoverRho2f = pPoverRho2;
#endif

	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (1./3.)/ph;
#else
	fNorm = 0.5*M_1_PI*ih2/ph;
#endif
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);    /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	if (TYPEQueryACTIVE(p)) {
		/* p active */
		pmumax = p->mumax;
		pa[0]=0.0;
		pa[1]=0.0;
		pa[2]=0.0;
		for (i=0;i<nSmooth;++i) {
			q = nnList[i].pPart;
			r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rq = rs1 * q->fMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
				+ nnList[i].fDist2*smf->H;

			qPoverRho2 = q->PoverRho2;
#ifdef PEXT
			{
			FLOAT qd2 = q->fDensity*q->fDensity;
			qPoverRho2f = (qPoverRho2*qd2-smf->Pext)/qd2;
			}
#else
			qPoverRho2f = qPoverRho2;
#endif
			if (TYPEQueryACTIVE(q)) {
				/* q active */
				rp = rs1 * pMass;
				if (dvdotdr>0.0) {
					p->uDotPdV += rq*PRES_PDV(pPoverRho2,qPoverRho2)*dvdotdr;
					q->uDotPdV += rp*PRES_PDV(qPoverRho2,pPoverRho2)*dvdotdr;
					rq *= (PRES_ACC(pPoverRho2f,qPoverRho2f));
					rp *= (PRES_ACC(pPoverRho2f,qPoverRho2f));
					rp *= aFac; /* convert to comoving acceleration */
					rq *= aFac;
					pa[0] -= rq * dx;
					pa[1] -= rq * dy;
					pa[2] -= rq * dz;
					ACCEL(q,0) += rp * dx;
					ACCEL(q,1) += rp * dy;
					ACCEL(q,2) += rp * dz;
					DEBUGFORCE("sphA");
              		}
				else {
             		/* h mean - using just hp probably ok */
					hav=0.5*(ph+sqrt(0.25*BALL2(q)));
					/* mu multiply by a to be consistent with physical c */
					absmu = -hav*dvdotdr*smf->a
						/(nnList[i].fDist2+0.01*hav*hav);
					/* mu terms for gas time step */
					if (absmu>pmumax) pmumax=absmu;
					if (absmu>q->mumax) q->mumax=absmu;
					/* viscosity term */

					visc = SWITCHCOMBINE(p,q)*
					  (smf->alpha*(pc + q->c) + smf->beta*2*absmu)
					  *absmu/(pDensity + q->fDensity);
					p->uDotPdV += rq*(PRES_PDV(pPoverRho2,q->PoverRho2))*dvdotdr;
					q->uDotPdV += rp*(PRES_PDV(q->PoverRho2,pPoverRho2))*dvdotdr;
					p->uDotAV += rq*(0.5*visc)*dvdotdr;
					q->uDotAV += rp*(0.5*visc)*dvdotdr;
					rq *= (PRES_ACC(pPoverRho2f,qPoverRho2f) + visc);
					rp *= (PRES_ACC(pPoverRho2f,qPoverRho2f) + visc);
					rp *= aFac; /* convert to comoving acceleration */
					rq *= aFac;

					pa[0] -= rq*dx;
					pa[1] -= rq*dy;
					pa[2] -= rq*dz;
					ACCEL(q,0) += rp*dx;
					ACCEL(q,1) += rp*dy;
					ACCEL(q,2) += rp*dz;
					DEBUGFORCE("sphAV");
              		}
				}
			else {
				/* q not active */
				if (dvdotdr>0.0) {

					p->uDotPdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
					rq *= (PRES_ACC(pPoverRho2f,qPoverRho2f));
					rq *= aFac; /* convert to comoving acceleration */

					pa[0] -= rq*dx;
					pa[1] -= rq*dy;
					pa[2] -= rq*dz;
					DEBUGFORCE("sphB ");
              		}
				else {
             		/* h mean */
					hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
					/* mu multiply by a to be consistent with physical c */
					absmu = -hav*dvdotdr*smf->a
						/(nnList[i].fDist2+0.01*hav*hav);
					/* mu terms for gas time step */
					if (absmu>pmumax) pmumax=absmu;
					/* viscosity term */

					visc = SWITCHCOMBINE(p,q)*
					  (smf->alpha*(pc + q->c) + smf->beta*2*absmu)
					  *absmu/(pDensity + q->fDensity);
					p->uDotPdV += rq*(PRES_PDV(pPoverRho2,q->PoverRho2))*dvdotdr;
					p->uDotAV += rq*(0.5*visc)*dvdotdr;
					rq *= (PRES_ACC(pPoverRho2f,qPoverRho2f) + visc);
					rq *= aFac; /* convert to comoving acceleration */

					pa[0] -= rq*dx;
					pa[1] -= rq*dy;
					pa[2] -= rq*dz;
					DEBUGFORCE("sphBV");
             		}
				}
	        }
		p->mumax = pmumax;
		ACCEL(p,0) += pa[0];
		ACCEL(p,1) += pa[1];
		ACCEL(p,2) += pa[2];
		}
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
	        q = nnList[i].pPart;
			if (!TYPEQueryACTIVE(q)) continue; /* neither active */

	        r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rp = rs1 * pMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) +
				nnList[i].fDist2*smf->H;

			qPoverRho2 = q->PoverRho2;
#ifdef PEXT
			{
			FLOAT qd2 = q->fDensity*q->fDensity;
			qPoverRho2f = (qPoverRho2*qd2-smf->Pext)/qd2;
			}
#else
			qPoverRho2f = qPoverRho2;
#endif

			if (dvdotdr>0.0) {
				q->uDotPdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
				rp *= (PRES_ACC(pPoverRho2f,qPoverRho2f));
				rp *= aFac; /* convert to comoving acceleration */

		        ACCEL(q,0) += rp*dx;
		        ACCEL(q,1) += rp*dy;
		        ACCEL(q,2) += rp*dz;
			DEBUGFORCE("sphC ");
				}
			else {
				/* h mean */
		        hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
			/* mu multiply by a to be consistent with physical c */
		        absmu = -hav*dvdotdr*smf->a
					/(nnList[i].fDist2+0.01*hav*hav);
				/* mu terms for gas time step */
				if (absmu>q->mumax) q->mumax=absmu;
				/* viscosity */

				visc = SWITCHCOMBINE(p,q)*
				  (smf->alpha*(pc + q->c) + smf->beta*2*absmu)
				  *absmu/(pDensity + q->fDensity);
				q->uDotPdV += rp*(PRES_PDV(q->PoverRho2,pPoverRho2))*dvdotdr;
				q->uDotAV += rp*(0.5*visc)*dvdotdr;
				rp *= (PRES_ACC(pPoverRho2f,qPoverRho2f) + visc);
				rp *= aFac; /* convert to comoving acceleration */

		        ACCEL(q,0) += rp*dx;
		        ACCEL(q,1) += rp*dy;
		        ACCEL(q,2) += rp*dz;
			DEBUGFORCE("sphCV");
				}
	        }
		}
	}

/* NB: ACCEL_PRES used here --
   with shock tracking disabled: #define NOSHOCKTRACK
   it is: a->aPres
   otherwise it is identical to p->a
   The postSphPressure function combines p->a and p->aPres
*/

/* Original Particle */
void initSphPressureParticle(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->uDotPdV = 0.0;
		}
	}

/* Cached copies of particle */
void initSphPressure(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->uDotPdV = 0.0;
		ACCEL_PRES(p,0) = 0.0;
		ACCEL_PRES(p,1) = 0.0;
		ACCEL_PRES(p,2) = 0.0;
		}
	}

void combSphPressure(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->uDotPdV += ((PARTICLE *)p2)->uDotPdV;
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL_PRES(p1,0) += ACCEL_PRES(p2,0);
		ACCEL_PRES(p1,1) += ACCEL_PRES(p2,1);
		ACCEL_PRES(p1,2) += ACCEL_PRES(p2,2);
		}
	}

void postSphPressure(PARTICLE *p, SMF *smf)
{
        if ( TYPEQuerySMOOTHACTIVE((PARTICLE *)p) ) {
            if ( TYPEQueryACTIVE((PARTICLE *)p) ) {
                    ACCEL_COMB_PRES(p,0);
                    ACCEL_COMB_PRES(p,1);
                    ACCEL_COMB_PRES(p,2);
                    }
            }
	}

/* Gather only version -- untested */
void SphPressure(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pa[3],pmumax;
	FLOAT ph,pc,pDensity,visc,hav,vFac,absmu;
	FLOAT fNorm,fNorm1,fNorm2;
	int i;

	assert(0);

	if (!TYPEQueryACTIVE(p)) return;

	pc = p->c;
	pDensity = p->fDensity;
	pPoverRho2 = p->PoverRho2;
	pmumax = p->mumax;
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (2./3.)/ph;
#else
	fNorm = M_1_PI*ih2/ph;
#endif
	fNorm1 = fNorm*ih2;
	fNorm2 = fNorm1*(smf->a);    /* Comoving accelerations */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	pa[0]=0.0;
	pa[1]=0.0;
	pa[2]=0.0;
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= q->fMass;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
	    if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
			dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
	    else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
	        dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;
		if (dvdotdr>0.0) {
			p->uDotPdV += rs1 * PRES_PDV(pPoverRho2,q->PoverRho2) * dvdotdr;
			rs1 *= (PRES_ACC(pPoverRho2,q->PoverRho2));
			pa[0] -= rs1 * dx;
			pa[1] -= rs1 * dy;
			pa[2] -= rs1 * dz;
			}
		else {
			hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
			/* mu
			   multiply by a to be consistent with physical c */
			absmu = -hav*dvdotdr*smf->a
			    / (nnList[i].fDist2+0.01*hav*hav);
			/* mu terms for gas time step */
			if (absmu>pmumax) pmumax=absmu;
			/* viscosity term */

			visc = SWITCHCOMBINE(p,q)*
			  (smf->alpha*(pc + q->c) + smf->beta*2*absmu)
			  *absmu/(pDensity + q->fDensity);
            p->uDotPdV += rs1 * (PRES_PDV(pPoverRho2,q->PoverRho2))*dvdotdr;
            p->uDotAV += rs1 * (0.5*visc)*dvdotdr;
			rs1 *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
			pa[0] -= rs1 * dx;
			pa[1] -= rs1 * dy;
			pa[2] -= rs1 * dz;
			}
 		}
	p->mumax = pmumax;
	ACCEL_PRES(p,0) += fNorm2*pa[0];
	ACCEL_PRES(p,0) += fNorm2*pa[1];
	ACCEL_PRES(p,0) += fNorm2*pa[2];
	}

void SphPressureSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pa[3],pMass,pmumax;
	FLOAT ph;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	pMass = p->fMass;
	pPoverRho2 = p->PoverRho2;
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (1./3.)/ph;
#else
	fNorm = 0.5*M_1_PI*ih2/ph;
#endif
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);    /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	if (TYPEQueryACTIVE(p)) {
		/* p active */
		pmumax = p->mumax;
		pa[0]=0.0;
		pa[1]=0.0;
		pa[2]=0.0;
		for (i=0;i<nSmooth;++i) {
			q = nnList[i].pPart;
			r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rq = rs1 * q->fMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
				+ nnList[i].fDist2*smf->H;

			if (TYPEQueryACTIVE(q)) {
				/* q active */
                rp = rs1 * pMass;
                p->uDotPdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
                q->uDotPdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
				rq *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rp *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rp *= aFac; /* convert to comoving acceleration */
				rq *= aFac;
				pa[0] -= rq * dx;
				pa[1] -= rq * dy;
				pa[2] -= rq * dz;
				ACCEL_PRES(q,0) += rp * dx;
				ACCEL_PRES(q,1) += rp * dy;
				ACCEL_PRES(q,2) += rp * dz;
              		        }
			else {
				/* q not active */
                p->uDotPdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
				rq *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rq *= aFac; /* convert to comoving acceleration */

				pa[0] -= rq*dx;
				pa[1] -= rq*dy;
				pa[2] -= rq*dz;
              		        }
             		}
		p->mumax = pmumax;
		ACCEL_PRES(p,0) += pa[0];
		ACCEL_PRES(p,1) += pa[1];
		ACCEL_PRES(p,2) += pa[2];
		}
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
	                q = nnList[i].pPart;
                        if (!TYPEQueryACTIVE(q)) continue; /* neither active */

                        r2 = nnList[i].fDist2*ih2;
                        DKERNEL(rs1,r2);
                        rs1 *= fNorm1;
			rp = rs1 * pMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) +
				nnList[i].fDist2*smf->H;
			q->uDotPdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
			rp *= (PRES_ACC(pPoverRho2,q->PoverRho2));
			rp *= aFac; /* convert to comoving acceleration */

		        ACCEL_PRES(q,0) += rp*dx;
		        ACCEL_PRES(q,1) += rp*dy;
		        ACCEL_PRES(q,2) += rp*dz;
	                }
		}
	}

/* Original Particle */
void initSphViscosityParticle(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
        ((PARTICLE *)p)->uDotAV = 0.0;
        }
}

/* Cached copies of particle */
void initSphViscosity(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->uDotAV = 0.0;
		ACCEL(p,0) = 0.0;
		ACCEL(p,1) = 0.0;
		ACCEL(p,2) = 0.0;
		}
	}

void combSphViscosity(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->uDotAV += ((PARTICLE *)p2)->uDotAV;
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL(p1,0) += ACCEL(p2,0);
		ACCEL(p1,1) += ACCEL(p2,1);
		ACCEL(p1,2) += ACCEL(p2,2);
		}
	}

/* Gather only */
void SphViscosity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
        assert(0);
	}

/* Symmetric Gather/Scatter version */
void SphViscositySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pa[3],pMass,pmumax;
	FLOAT ph,pc,pDensity,visc,hav,absmu;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (1./3.)/ph;
#else
	fNorm = 0.5*M_1_PI*ih2/ph;
#endif
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);    /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	if (TYPEQueryACTIVE(p)) {
		/* p active */
		pmumax = p->mumax;
		pa[0]=0.0;
		pa[1]=0.0;
		pa[2]=0.0;
		for (i=0;i<nSmooth;++i) {
			q = nnList[i].pPart;
			r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rq = rs1 * q->fMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
				+ nnList[i].fDist2*smf->H;
			if (dvdotdr > 0.0) continue;

			if (TYPEQueryACTIVE(q)) {
				/* q active */
				rp = rs1 * pMass;
				/* h mean - using just hp probably ok */
				hav=0.5*(ph+sqrt(0.25*BALL2(q)));
				/* mu multiply by a to be consistent with physical c */
				absmu = -hav*dvdotdr*smf->a
				  /(nnList[i].fDist2+0.01*hav*hav);
				/* mu terms for gas time step */
				if (absmu>pmumax) pmumax=absmu;
				if (absmu>q->mumax) q->mumax=absmu;
				/* viscosity term */

				visc =
				  (SWITCHCOMBINEA(p,q)*smf->alpha*(pc + q->c)
				   + SWITCHCOMBINEB(p,q)*smf->beta*2*absmu)
				  *absmu/(pDensity + q->fDensity);

				p->uDotAV += rq*0.5*visc*dvdotdr;
				q->uDotAV += rp*0.5*visc*dvdotdr;
				rq *= visc;
				rp *= visc;
				rp *= aFac; /* convert to comoving acceleration */
				rq *= aFac;

				pa[0] -= rq*dx;
				pa[1] -= rq*dy;
				pa[2] -= rq*dz;
				ACCEL(q,0) += rp*dx;
				ACCEL(q,1) += rp*dy;
				ACCEL(q,2) += rp*dz;
				}
			else {
				/* q not active */
			        /* h mean */
			        hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
				/* mu multiply by a to be consistent with physical c */
				absmu = -hav*dvdotdr*smf->a
				  /(nnList[i].fDist2+0.01*hav*hav);
				/* mu terms for gas time step */
				if (absmu>pmumax) pmumax=absmu;
				/* viscosity term */

				visc =
				  (SWITCHCOMBINEA(p,q)*smf->alpha*(pc + q->c)
				   + SWITCHCOMBINEB(p,q)*smf->beta*2*absmu)
				  *absmu/(pDensity + q->fDensity);

				p->uDotAV += rq*0.5*visc*dvdotdr;
				rq *= visc;
				rq *= aFac; /* convert to comoving acceleration */

				pa[0] -= rq*dx;
				pa[1] -= rq*dy;
				pa[2] -= rq*dz;
				}
	                }
		p->mumax = pmumax;
		ACCEL(p,0) += pa[0];
		ACCEL(p,1) += pa[1];
		ACCEL(p,2) += pa[2];
		}
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
           	        q = nnList[i].pPart;
			if (!TYPEQueryACTIVE(q)) continue; /* neither active */

	                r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rp = rs1 * pMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) +
				nnList[i].fDist2*smf->H;
			if (dvdotdr > 0.0) continue;

			/* h mean */
		        hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
			/* mu multiply by a to be consistent with physical c */
		        absmu = -hav*dvdotdr*smf->a
			  /(nnList[i].fDist2+0.01*hav*hav);
				/* mu terms for gas time step */
			if (absmu>q->mumax) q->mumax=absmu;
				/* viscosity */

			visc =
			  (SWITCHCOMBINEA(p,q)*smf->alpha*(pc + q->c)
			   + SWITCHCOMBINEB(p,q)*smf->beta*2*absmu)
			  *absmu/(pDensity + q->fDensity);

			q->uDotAV += rp*0.5*visc*dvdotdr;
			rp *= visc;
			rp *= aFac; /* convert to comoving acceleration */

			ACCEL(q,0) += rp*dx;
			ACCEL(q,1) += rp*dy;
			ACCEL(q,2) += rp*dz;
	                }
		}
	}



void initDivVort(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p )) {
		((PARTICLE *)p)->divv = 0.0;
		((PARTICLE *)p)->curlv[0] = 0.0;
		((PARTICLE *)p)->curlv[1] = 0.0;
		((PARTICLE *)p)->curlv[2] = 0.0;
#ifdef SHOCKTRACK
		((PARTICLE *)p)->divrhov = 0.0;
		((PARTICLE *)p)->gradrho[0] = 0.0;
		((PARTICLE *)p)->gradrho[1] = 0.0;
		((PARTICLE *)p)->gradrho[2] = 0.0;
#endif
		}
	}

void combDivVort(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1 )) {
		((PARTICLE *)p1)->divv += ((PARTICLE *)p2)->divv;
		((PARTICLE *)p1)->curlv[0] += ((PARTICLE *)p2)->curlv[0];
		((PARTICLE *)p1)->curlv[1] += ((PARTICLE *)p2)->curlv[1];
		((PARTICLE *)p1)->curlv[2] += ((PARTICLE *)p2)->curlv[2];
#ifdef SHOCKTRACK
		((PARTICLE *)p1)->divrhov += ((PARTICLE *)p2)->divrhov;
		((PARTICLE *)p1)->gradrho[0] += ((PARTICLE *)p2)->gradrho[0];
		((PARTICLE *)p1)->gradrho[1] += ((PARTICLE *)p2)->gradrho[1];
		((PARTICLE *)p1)->gradrho[2] += ((PARTICLE *)p2)->gradrho[2];
#endif
		}
	}

/* Gather only version -- untested */
void DivVort(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pcurlv[3],pdivv;
	FLOAT pDensity;
	FLOAT fNorm,vFac,a2;
	int i;

	if (!TYPEQueryACTIVE(p)) return;

	pDensity = p->fDensity;
	ih2 = 4.0/BALL2(p);
	a2 = (smf->a*smf->a);
#ifdef SPH1D
	fNorm = (2./3.)*ih2*sqrt(ih2);
#else
	fNorm = M_1_PI*ih2*ih2*sqrt(ih2);
#endif
	vFac = (smf->bCannonical ? 1./a2 : 1.0); /* converts v to xdot */

	pdivv=0.0;
	pcurlv[0]=0.0;
	pcurlv[1]=0.0;
	pcurlv[2]=0.0;
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= q->fMass/pDensity;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
	    if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
			dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
	    else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
	        dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;
		pdivv += rs1*dvdotdr;
		pcurlv[0] += rs1*(dvz*dy - dvy*dz);
		pcurlv[1] += rs1*(dvx*dz - dvz*dx);
		pcurlv[2] += rs1*(dvy*dx - dvx*dy);
 		}
	p->divv -=  fNorm*pdivv;  /* physical */
	p->curlv[0] += fNorm*vFac*pcurlv[0];
	p->curlv[1] += fNorm*vFac*pcurlv[1];
	p->curlv[2] += fNorm*vFac*pcurlv[2];
	}

/* Output is physical divv and curlv -- thus a*h_co*divv is physical */
void DivVortSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pMass;
	FLOAT fNorm,dv,vFac,a2;
	int i;

	mdlassert(smf->pkd->mdl, TYPETest(p,TYPE_GAS));

	pMass = p->fMass;
	ih2 = 4.0/BALL2(p);
	a2 = (smf->a*smf->a);
#ifdef SPH1D
	fNorm = (1./3.)*ih2*sqrt(ih2);
#else
	fNorm = 0.5*M_1_PI*ih2*ih2*sqrt(ih2);
#endif
	vFac = (smf->bCannonical ? 1./a2 : 1.0); /* converts v to xdot */

//#define GRADWONESIDED
#ifdef GRADWONESIDED
    fNorm *= 2;
#define QISACTIVE(q) (0)
#else
#define QISACTIVE(q) (TYPEQueryACTIVE(q))
#endif
#if (1)
#define DXFUNC(d) (d)
#define DYFUNC(d) (d)
#define DZFUNC(d) (d)
#else
#define DXFUNC(d) fabs(d)
#define DYFUNC(d) fabs(d)
#define DZFUNC(d) fabs(d)
#endif
#ifdef GRADW
#define DENSMULP(p,q) (1/q->fDensity)
#define DENSMULQ(p,q) (1/p->fDensity)
/* This gives grad 1_i = rhoi sum 0.5*(dWij+dWji)*m_j/(rhoi*rhoj)
                       =      sum 0.5*(dWij+dWji)*m_j/rhoj
   Same form as pressure acceleration (eqv to constant P case)
   Volume weighting should be really good for grad 1
   Ok for div v = div v - v . grad 1  */
    if (TYPETest(p,TYPE_FEEDBACK)) printf("GRADW %d: nS %d\n",p->iOrder,nSmooth);
#else
#define DENSMULP(p,q) (1/p->fDensity)
#define DENSMULQ(p,q) (1/q->fDensity)
/* This gives grad 1_i = 1/rhoi sum 0.5*(dWij+dWji)*m_j
   The volume weighting is really bad on this for grad 1
   Correct for div v = 1/rho div (rho v) - v . grad rho  (less noisy divv) */
#endif

	if (TYPEQueryACTIVE( p )) {
		/* p active */
		for (i=0;i<nSmooth;++i) {
	        q = nnList[i].pPart;
		mdlassert(smf->pkd->mdl, TYPETest(q,TYPE_GAS));
	        r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm;
			rq = rs1 * q->fMass*DENSMULP(p,q);

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) +
				nnList[i].fDist2*smf->H;

			if (QISACTIVE(q)) {
				/* q active */
				rp = rs1 * pMass*DENSMULQ(p,q);
				p->divv -= rq*dvdotdr;
				q->divv -= rp*dvdotdr;
#ifdef GRADW
				dv=DXFUNC(dx);
#else
				dv=vFac*(dvz*dy - dvy*dz);
#endif
				p->curlv[0] += rq*dv;
				q->curlv[0] += rp*dv;
#ifdef GRADW
				dv=DYFUNC(dy);
#else
				dv=vFac*(dvx*dz - dvz*dx);
#endif
				p->curlv[1] += rq*dv;
				q->curlv[1] += rp*dv;
#ifdef GRADW
				dv=DZFUNC(dz);
#else
				dv=vFac*(dvy*dx - dvx*dy);
#endif
				p->curlv[2] += rq*dv;
				q->curlv[2] += rp*dv;

#ifdef SHOCKTRACK
				p->divrhov -= rs1*dvdotdr*q->fMass;
				q->divrhov -= rs1*dvdotdr*pMass;
				p->gradrho[0] += rs1*q->fMass*dx;
				q->gradrho[0] -= rs1*pMass*dx;
				p->gradrho[1] += rs1*q->fMass*dy;
				q->gradrho[1] -= rs1*pMass*dy;
				p->gradrho[2] += rs1*q->fMass*dz;
				q->gradrho[2] -= rs1*pMass*dz;
#endif
		        }
			else {
		        /* q inactive */
				p->divv -= rq*dvdotdr;
#ifdef GRADW
				dv=DXFUNC(dx);
#else
				dv=vFac*(dvz*dy - dvy*dz);
#endif
				p->curlv[0] += rq*dv;
#ifdef GRADW
				dv=DYFUNC(dy);
#else
				dv=vFac*(dvx*dz - dvz*dx);
#endif
				p->curlv[1] += rq*dv;
#ifdef GRADW
				dv=DZFUNC(dz);
#else
				dv=vFac*(dvy*dx - dvx*dy);
#endif
				p->curlv[2] += rq*dv;

#ifdef SHOCKTRACK
				p->divrhov -= rs1*dvdotdr*q->fMass;
				p->gradrho[0] += rs1*q->fMass*dx;
				p->gradrho[1] += rs1*q->fMass*dy;
				p->gradrho[2] += rs1*q->fMass*dz;
#endif
		        }
	        }
		}
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
	        q = nnList[i].pPart;
            mdlassert(smf->pkd->mdl, TYPETest(q,TYPE_GAS));
            if (!QISACTIVE(q)) continue; /* neither active */

			r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *=fNorm;
			rp = rs1*pMass*DENSMULQ(p,q);

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
            if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
                dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
            else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
                dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
				+ nnList[i].fDist2*smf->H;
			/* q active */
			q->divv -= rp*dvdotdr;
#ifdef GRADW
            dv=DXFUNC(dx);
#else
            dv=vFac*(dvz*dy - dvy*dz);
#endif
            q->curlv[0] += rp*dv;
#ifdef GRADW
            dv=DYFUNC(dy);
#else
            dv=vFac*(dvx*dz - dvz*dx);
#endif
            q->curlv[1] += rp*dv;
#ifdef GRADW
            dv=DZFUNC(dz);
#else
            dv=vFac*(dvy*dx - dvx*dy);
#endif
            q->curlv[2] += rp*dv;

#ifdef SHOCKTRACK
			q->divrhov -= rs1*dvdotdr*pMass;
			q->gradrho[0] -= rs1*pMass*dx;
			q->gradrho[1] -= rs1*pMass*dy;
			q->gradrho[2] -= rs1*pMass*dz;
#endif
	        }
		}
	}

void initSurfaceNormal(void *p)
{
	}

void combSurfaceNormal(void *p1,void *p2)
{
	}

/* Gather only version */
void SurfaceNormal(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef SURFACEAREA
	FLOAT ih2,r2,rs,rs1,fNorm1;
	FLOAT dx,dy,dz;
	FLOAT grx=0,gry=0,grz=0, rcostheta;
	PARTICLE *q;
	int i;

	ih2 = 4.0/BALL2(p);

	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2*ih2;
	    q = nnList[i].pPart;
	    DKERNEL(rs1,r2);
	    rs1 *= q->fMass;
	    dx = nnList[i].dx;
	    dy = nnList[i].dy;
	    dz = nnList[i].dz;
	    grx += dx*rs1; /* Grad rho estimate */
	    gry += dy*rs1;
	    grz += dz*rs1;
	    }

	fNorm1 = -1/sqrt(grx*grx+gry*gry+grz*grz); /* to unit vector */
#ifdef NORMAL
        p->normal[0] = p->curlv[0] = grx *= fNorm1; /*REMOVE */
	p->normal[1] = p->curlv[1] = gry *= fNorm1;
	p->normal[2] = p->curlv[2] = grz *= fNorm1;
#else
        p->curlv[0] = grx *= fNorm1; /*REMOVE */
	p->curlv[1] = gry *= fNorm1;
	p->curlv[2] = grz *= fNorm1;
#endif

	p->fArea = 1;
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2;
	    rs = r2*ih2;
	    q = nnList[i].pPart;
	    DKERNEL(rs1,rs);
	    rs1 *= q->fMass;
	    dx = nnList[i].dx;
	    dy = nnList[i].dy;
	    dz = nnList[i].dz;
	    rcostheta = -(dx*grx + dy*gry + dz*grz);
	    /* cos^2 =  3/7. corresponds to an angle of about 49 degrees which is
	       the angle from vertical to the point from p x on hexagonal lattice
	       q  x  q  x  q
	          |   /
		  | /
                  p
	    */
/*
	    if ((p->iOrder % 10000)==0) {
		printf("Particle on Edge? %d %d, %g %g %g  %g %g %g   %g %g\n",p->iOrder,q->iOrder,grx,gry,grz,dx,dy,dz,rcostheta*fabs(rcostheta),r2);
		}
*/

	    if (rcostheta > 0 && rcostheta*rcostheta > (3/7.)*r2) {
		p->fArea = 0; /* You aren't on the edge */
		break;
		}
	    }
/*
	    if (p->fArea > 0 ) printf("Particle on Edge %d, %g %g %g\n",p->iOrder,p->r[0],p->r[1],p->r[2]);
*/

#endif
}

void initSurfaceArea(void *p)
{
	}

void combSurfaceArea(void *p1,void *p2)
{
	}

/* Gather only version */
void SurfaceArea(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef SURFACEAREA
	FLOAT ih2,ih,r2,rs,fColumnDensity,fNorm,fNorm1;
	FLOAT grx=0,gry=0,grz=0,dx,dy,dz;
	PARTICLE *q;
	int i;
	unsigned int qiActive;

	if (p->fArea > 0) {

	    ih2 = 4.0/BALL2(p);
	    ih = sqrt(ih2);
	    fNorm = 10./7.*M_1_PI*ih2;

	    for (i=0;i<nSmooth;++i) { /* Average surface normal */
		q = nnList[i].pPart;
		if (q->fArea > 0) {
		    r2 = nnList[i].fDist2;
		    dx = nnList[i].dx;
		    dy = nnList[i].dy;
		    dz = nnList[i].dz;
		    rs = dx*q->curlv[0]+dy*q->curlv[1]+dz*q->curlv[2];
		    assert(rs*rs <= r2); /* REMOVE */
		    r2 = (r2-rs*rs)*ih2;
		    KERNEL(rs,r2);
		    r2 *= q->fMass;
		    grx += rs*q->curlv[0];
		    gry += rs*q->curlv[1];
		    grz += rs*q->curlv[2];
		    }
		}

	    fNorm1 = 1/sqrt(grx*grx+gry*gry+grz*grz); /* to unit vector */
#ifdef NORMAL
	    p->normal[0] = p->curlv[0] = grx *= fNorm1; /*REMOVE */
	    p->normal[1] = p->curlv[1] = gry *= fNorm1;
	    p->normal[2] = p->curlv[2] = grz *= fNorm1;
#else
	    p->curlv[0] = grx *= fNorm1; /*REMOVE */
	    p->curlv[1] = gry *= fNorm1;
	    p->curlv[2] = grz *= fNorm1;
#endif
	    fColumnDensity = 0.0;
	    for (i=0;i<nSmooth;++i) { /* Calculate Column Density */
		q = nnList[i].pPart;
		if (q->fArea > 0) {
		    r2 = nnList[i].fDist2;
		    dx = nnList[i].dx;
		    dy = nnList[i].dy;
		    dz = nnList[i].dz;
		    rs = dx*grx+dy*gry+dz*grz;
		    assert(rs*rs <= r2); /* REMOVE */
		    r2 = (r2-rs*rs)*ih2;
		    KERNEL(rs,r2);
		    fColumnDensity += rs*q->fMass;
		    }
		}
	    p->fArea = p->fMass/(fColumnDensity*fNorm);
	    }
#endif
}



void initPreCorrDensPartitionParticle(void *p)
{
#ifdef CORRDENSPART
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
    ((PARTICLE *)p)->Eta_Factor=0.0;
  }
#endif
  }

  void initPreCorrDensPartition(void *p)
  {
#ifdef CORRDENSPART
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
    ((PARTICLE *)p)->Eta_Factor=0.0;
  }
#endif
    }

  void combPreCorrDensPartition(void *p1,void *p2)
  {
#ifdef CORRDENSPART
  if (TYPEQueryACTIVE((PARTICLE *) p1)) {
    ((PARTICLE *)p1)->Eta_Factor += ((PARTICLE *)p2)->Eta_Factor;
  }
#endif
    }

void PreCorrDensPartition(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef CORRDENSPART
  FLOAT ih2,ih,r2,rs,fDensity,fNorm,fNorm1,vFac,aFac;
  FLOAT dvx,dvy,dvz,dx,dy,dz,divvi,divvj,drs,rq,rp,qih2,Va;
  FLOAT corr=0.0;
  PARTICLE *q;
  int i;
  ih2 = 4.0/BALL2(p);
  ih = sqrt(ih2);
  fNorm = M_1_PI*ih2*ih*0.5;
  fNorm1 = fNorm*ih2;
  for (i=0;i<nSmooth;++i) {
    r2 = nnList[i].fDist2*ih2;
    q = nnList[i].pPart;
      dx = nnList[i].dx; /* NB: dx = px - qx */
      dy = nnList[i].dy;
      dz = nnList[i].dz;
    if (!TYPEQueryACTIVE(p) && !TYPEQueryACTIVE(q)) continue;
    KERNEL(rs,r2);
    DKERNEL(drs,r2);
    #ifdef COMBOKERNEL
      double rs2,rs3,rs4,rs5,rs6,rs7,rs8,rs9,rs10,rs11,rs12,rs13,rs14;
      double drs2,drs3,drs4,drs5,drs6,drs7,drs8,drs9,drs10,drs11,drs12,drs13,drs14;
      double alp=smf->alp; double alp2=smf->alp2; double alp3=smf->alp3;
      double alp4=smf->alp4; double alp5=smf->alp5;double alp6=smf->alp6;
      double alp7=smf->alp7; double alp8=smf->alp8;double alp9=smf->alp9;
      double alp10=smf->alp10; double alp11=smf->alp11; double alp12=smf->alp12;
      double alp13=smf->alp13;
      KERNEL2(rs2,r2);
      KERNEL3(rs3,r2);
      KERNEL4(rs4,r2);
      KERNEL5(rs5,r2);
      KERNEL6(rs6,r2);
      KERNEL7(rs7,r2);
      KERNEL8(rs8,r2);
      KERNEL9(rs9,r2);
      KERNEL10(rs10,r2);
      KERNEL11(rs11,r2);
      KERNEL12(rs12,r2);
      KERNEL13(rs13,r2);
      KERNEL14(rs14,r2);
      DKERNEL2(drs2,r2);
      DKERNEL3(drs3,r2);
      DKERNEL4(drs4,r2);
      DKERNEL5(drs5,r2);
      DKERNEL6(drs6,r2);
      DKERNEL7(drs7,r2);
      DKERNEL8(drs8,r2);
      DKERNEL9(drs9,r2);
      DKERNEL10(drs10,r2);
      DKERNEL11(drs11,r2);
      DKERNEL12(drs12,r2);
      DKERNEL13(drs13,r2);
      DKERNEL14(drs14,r2);

double rs15,rs16,rs17,rs18,rs19,rs20,rs21,rs22;
double drs15,drs16,drs17,drs18,drs19,drs20,drs21,drs22;
      double alp14=smf->alp14; double alp15=smf->alp15;double alp16=smf->alp16;
      double alp17=smf->alp17; double alp18=smf->alp18; double alp19=smf->alp19;
      double alp20=smf->alp20; double alp21=smf->alp21;
      KERNEL19(rs15,r2);
      KERNEL20(rs16,r2);
      KERNEL21(rs17,r2);
      KERNEL22(rs18,r2);
      KERNEL23(rs19,r2);
      KERNEL24(rs20,r2);
      KERNEL25(rs21,r2);
      KERNEL26(rs22,r2);
      DKERNEL19(drs15,r2);
      DKERNEL20(drs16,r2);
      DKERNEL21(drs17,r2);
      DKERNEL22(drs18,r2);
      DKERNEL23(drs19,r2);
      DKERNEL24(drs20,r2);
      DKERNEL25(drs21,r2);
      DKERNEL26(drs22,r2);
      rs = ((((((((((((((((((((((1.0-alp)*rs+alp*rs2)*(1.0-alp2)+alp2*rs3)*(1.0-alp3)+alp3*rs4)*(1.0-alp4)+alp4*rs5)*(1.0-alp5)+alp5*rs6)*(1.0-alp6)+alp6*rs7)*(1.0-alp7)+alp7*rs8)*(1.0-alp8)+alp8*rs9)*(1.0-alp9)+alp9*rs10)*(1.0-alp10)+alp10*rs11)*(1.0-alp11)+alp11*rs12)*(1.0-alp12)+alp12*rs13)*(1.0-alp13)+alp13*rs14)*(1.0-alp14)+alp14*rs15)*(1.0-alp15)+alp15*rs16)*(1.0-alp16)+alp16*rs17)*(1.0-alp17)+alp17*rs18)*(1.0-alp18)+alp18*rs19)*(1.0-alp19)+alp19*rs20)*(1.0-alp20)+alp20*rs21)*(1.0-alp21)+alp21*rs22);
      drs = ((((((((((((((((((((((1.0-alp)*drs+alp*drs2)*(1.0-alp2)+alp2*drs3)*(1.0-alp3)+alp3*drs4)*(1.0-alp4)+alp4*drs5)*(1.0-alp5)+alp5*drs6)*(1.0-alp6)+alp6*drs7)*(1.0-alp7)+alp7*drs8)*(1.0-alp8)+alp8*drs9)*(1.0-alp9)+alp9*drs10)*(1.0-alp10)+alp10*drs11)*(1.0-alp11)+alp11*drs12)*(1.0-alp12)+alp12*drs13)*(1.0-alp13)+alp13*drs14)*(1.0-alp14)+alp14*drs15)*(1.0-alp15)+alp15*drs16)*(1.0-alp16)+alp16*drs17)*(1.0-alp17)+alp17*drs18)*(1.0-alp18)+alp18*drs19)*(1.0-alp19)+alp19*drs20)*(1.0-alp20)+alp20*drs21)*(1.0-alp21)+alp21*drs22);
#endif
    double ab;
    rs *= fNorm;
    drs *= fNorm1;
    ab = 0.0;

        //rq = 1.0*rs*q->fMass/(pow(p->fDensity,ab)*pow(q->fDensity,1.0-ab));
    //rp = 1.0*rs*p->fMass/(pow(p->fDensity,1.0-ab)*pow(q->fDensity,ab));
    //one could also calculate p->Q4 before and use it here.
    //rq = -1.0*drs*q->fMass/(pow(p->fDensity,ab)*pow(q->fDensity,1.0-ab))*nnList[i].fDist2/3.0;
    //rp = -1.0*drs*p->fMass/(pow(p->fDensity,1.0-ab)*pow(q->fDensity,ab))*nnList[i].fDist2/3.0
    rq=p->fDensity*2.0;
    //rq=1.0;
      //rq = (q->fMass/q->fDensity*rs)*p->fDensity*2.0;
//rp = (p->fMass/p->fDensity*rs)*q->fDensity*2.0;

//eta_p = Q4_p * p->fDensity
//rq = -p->fDensity*(q->fMass/q->fDensity*drs)*nnList[i].fDist2/3.0;
//rp = -q->fDensity*(p->fMass/p->fDensity*drs)*nnList[i].fDist2/3.0;

//rq = p->fDensity;

if (TYPEQueryACTIVE(p)) {
    p->Eta_Factor += rq;
 }
 #ifdef CORRETASYM
 if (TYPEQueryACTIVE(q)) {
    q->Eta_Factor += rp;
  }
 #endif
  }
  
  #endif
  }


void initCorrDensPartitionParticle(void *p)
{
#ifdef CORRDENSPART
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
    ((PARTICLE *)p)->fDensity_Corrector=0.0;
  }
#endif
  }

  void initCorrDensPartition(void *p)
  {
#ifdef CORRDENSPART
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
    ((PARTICLE *)p)->fDensity_Corrector=0.0;
  }
#endif
    }

  void combCorrDensPartition(void *p1,void *p2)
  {
#ifdef CORRDENSPART
  if (TYPEQueryACTIVE((PARTICLE *) p1)) {
    ((PARTICLE *)p1)->fDensity_Corrector += ((PARTICLE *)p2)->fDensity_Corrector;
  }
#endif
    }

void CorrDensPartition(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef CORRDENSPART
  FLOAT ih2,ih,r2,rs,fDensity,fNorm,fNorm1,vFac,aFac;
  FLOAT dvx,dvy,dvz,dx,dy,dz,divvi,divvj,drs,rq,rp,qih2,Va;
  FLOAT corr=0.0;
  PARTICLE *q;
  int i;
  ih2 = 4.0/BALL2(p);
  ih = sqrt(ih2);
  #ifdef SINC
      double nloc = 9.0;
      double norm =0.027012593+0.020410827*sqrt(nloc)+0.0037451957*nloc+0.047013839*nloc*sqrt(nloc);
      fNorm = norm*ih2*ih;
  #else
  fNorm = M_1_PI*ih2*ih;
  #ifdef CORRDENSSYM
  fNorm = fNorm*0.5;
  #endif
  #endif
    fNorm1 = fNorm*ih2;
    //  printf(" %f ",p->Eta_Factor);
  for (i=0;i<nSmooth;++i) {
    r2 = nnList[i].fDist2*ih2;
    q = nnList[i].pPart;
      dx = nnList[i].dx; /* NB: dx = px - qx */
      dy = nnList[i].dy;
      dz = nnList[i].dz;
    if (!TYPEQueryACTIVE(p) && !TYPEQueryACTIVE(q)) continue;
    KERNEL(rs,r2);
    DKERNEL(drs,r2);
    #ifdef COMBOKERNEL
      double rs2,rs3,rs4,rs5,rs6,rs7,rs8,rs9,rs10,rs11,rs12,rs13,rs14;
      double drs2,drs3,drs4,drs5,drs6,drs7,drs8,drs9,drs10,drs11,drs12,drs13,drs14;
      double alp=smf->alp; double alp2=smf->alp2; double alp3=smf->alp3;
      double alp4=smf->alp4; double alp5=smf->alp5;double alp6=smf->alp6;
      double alp7=smf->alp7; double alp8=smf->alp8;double alp9=smf->alp9;
      double alp10=smf->alp10; double alp11=smf->alp11; double alp12=smf->alp12;
      double alp13=smf->alp13;
      KERNEL2(rs2,r2);
      KERNEL3(rs3,r2);
      KERNEL4(rs4,r2);
      KERNEL5(rs5,r2);
      KERNEL6(rs6,r2);
      KERNEL7(rs7,r2);
      KERNEL8(rs8,r2);
      KERNEL9(rs9,r2);
      KERNEL10(rs10,r2);
      KERNEL11(rs11,r2);
      KERNEL12(rs12,r2);
      KERNEL13(rs13,r2);
      KERNEL14(rs14,r2);
      DKERNEL2(drs2,r2);
      DKERNEL3(drs3,r2);
      DKERNEL4(drs4,r2);
      DKERNEL5(drs5,r2);
      DKERNEL6(drs6,r2);
      DKERNEL7(drs7,r2);
      DKERNEL8(drs8,r2);
      DKERNEL9(drs9,r2);
      DKERNEL10(drs10,r2);
      DKERNEL11(drs11,r2);
      DKERNEL12(drs12,r2);
      DKERNEL13(drs13,r2);
      DKERNEL14(drs14,r2);

double rs15,rs16,rs17,rs18,rs19,rs20,rs21,rs22;
double drs15,drs16,drs17,drs18,drs19,drs20,drs21,drs22;
      double alp14=smf->alp14; double alp15=smf->alp15;double alp16=smf->alp16;
      double alp17=smf->alp17; double alp18=smf->alp18; double alp19=smf->alp19;
      double alp20=smf->alp20; double alp21=smf->alp21;
      KERNEL19(rs15,r2);
      KERNEL20(rs16,r2);
      KERNEL21(rs17,r2);
      KERNEL22(rs18,r2);
      KERNEL23(rs19,r2);
      KERNEL24(rs20,r2);
      KERNEL25(rs21,r2);
      KERNEL26(rs22,r2);
      DKERNEL19(drs15,r2);
      DKERNEL20(drs16,r2);
      DKERNEL21(drs17,r2);
      DKERNEL22(drs18,r2);
      DKERNEL23(drs19,r2);
      DKERNEL24(drs20,r2);
      DKERNEL25(drs21,r2);
      DKERNEL26(drs22,r2);
      rs = ((((((((((((((((((((((1.0-alp)*rs+alp*rs2)*(1.0-alp2)+alp2*rs3)*(1.0-alp3)+alp3*rs4)*(1.0-alp4)+alp4*rs5)*(1.0-alp5)+alp5*rs6)*(1.0-alp6)+alp6*rs7)*(1.0-alp7)+alp7*rs8)*(1.0-alp8)+alp8*rs9)*(1.0-alp9)+alp9*rs10)*(1.0-alp10)+alp10*rs11)*(1.0-alp11)+alp11*rs12)*(1.0-alp12)+alp12*rs13)*(1.0-alp13)+alp13*rs14)*(1.0-alp14)+alp14*rs15)*(1.0-alp15)+alp15*rs16)*(1.0-alp16)+alp16*rs17)*(1.0-alp17)+alp17*rs18)*(1.0-alp18)+alp18*rs19)*(1.0-alp19)+alp19*rs20)*(1.0-alp20)+alp20*rs21)*(1.0-alp21)+alp21*rs22);
      drs = ((((((((((((((((((((((1.0-alp)*drs+alp*drs2)*(1.0-alp2)+alp2*drs3)*(1.0-alp3)+alp3*drs4)*(1.0-alp4)+alp4*drs5)*(1.0-alp5)+alp5*drs6)*(1.0-alp6)+alp6*drs7)*(1.0-alp7)+alp7*drs8)*(1.0-alp8)+alp8*drs9)*(1.0-alp9)+alp9*drs10)*(1.0-alp10)+alp10*drs11)*(1.0-alp11)+alp11*drs12)*(1.0-alp12)+alp12*drs13)*(1.0-alp13)+alp13*drs14)*(1.0-alp14)+alp14*drs15)*(1.0-alp15)+alp15*drs16)*(1.0-alp16)+alp16*drs17)*(1.0-alp17)+alp17*drs18)*(1.0-alp18)+alp18*drs19)*(1.0-alp19)+alp19*drs20)*(1.0-alp20)+alp20*drs21)*(1.0-alp21)+alp21*drs22);
#endif
    FLOAT etafacp,etafacq;
    rs *= fNorm;
    drs *= fNorm1;
    //etafacp = (q->Eta_Factor != 0 ? 1.0/q->Eta_Factor : 1.0);
    //etafacq = (p->Eta_Factor != 0 ? 1.0/p->Eta_Factor : 1.0);
    etafacp = p->Eta_Factor/q->Eta_Factor;
    etafacq = q->Eta_Factor/p->Eta_Factor;
    rq = rs*q->fMass*etafacp/p->fDensity;
    rp = rs*p->fMass*etafacq/q->fDensity;

if (TYPEQueryACTIVE(p)) {
    p->fDensity_Corrector += rq;
  }
#ifdef CORRDENSSYM
 if (TYPEQueryACTIVE(q)) {
    q->fDensity_Corrector += rp;
  }
#endif
  }
  #endif
  }


void initCorrEntropyParticle(void *p)
{
#ifdef CORRENTROPYPART
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
    ((PARTICLE *)p)->fDivv_Corrector=0.0;
    ((PARTICLE *)p)->fDivv_Corrector2=0.0;
  }
#endif
}

void initCorrEntropy(void *p)
{
#ifdef CORRENTROPYPART
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
    ((PARTICLE *)p)->fDivv_Corrector=0.0;
    ((PARTICLE *)p)->fDivv_Corrector2=0.0;
  }
#endif
}

void combCorrEntropy(void *p1,void *p2)
{
#ifdef CORRENTROPYPART
  if (TYPEQueryACTIVE((PARTICLE *) p1)) {
    ((PARTICLE *)p1)->fDivv_Corrector += ((PARTICLE *)p2)->fDivv_Corrector;
    ((PARTICLE *)p1)->fDivv_Corrector2 += ((PARTICLE *)p2)->fDivv_Corrector2;
  }
#endif
}

void CorrEntropy(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef CORRENTROPYPART
  FLOAT ih2,ih,r2,rs,fDensity,fNorm,fNorm1,vFac,aFac;
  FLOAT dvx,dvy,dvz,dx,dy,dz,divvi,divvj,rs1,rq,rp,qih2,Va;
  PARTICLE *q;
  int i;
  fDensity=p->fDensity;
  ih2 = 4.0/BALL2(p);
  ih = sqrt(ih2);
  #ifdef SINC
  double nloc = 9.0;
  double norm =0.027012593+0.020410827*sqrt(nloc)+0.0037451957*nloc+0.047013839*nloc*sqrt(nloc);
  fNorm = norm*ih2*ih;
  #else
  fNorm = M_1_PI*ih2*ih;
  #endif
  fNorm1 = fNorm*ih2*0.5;
  for (i=0;i<nSmooth;++i) {
    r2 = nnList[i].fDist2*ih2;
    q = nnList[i].pPart;
    if (!TYPEQueryACTIVE(p) && !TYPEQueryACTIVE(q)) continue;
    DKERNEL(rs1,r2);
#ifdef COMBOKERNEL
    double drs2,drs3,drs4,drs5,drs6,drs7,drs8,drs9,drs10,drs11,drs12,drs13,drs14;
    double alp=smf->alp; double alp2=smf->alp2; double alp3=smf->alp3;
    double alp4=smf->alp4; double alp5=smf->alp5;double alp6=smf->alp6;
    double alp7=smf->alp7; double alp8=smf->alp8;double alp9=smf->alp9;
    double alp10=smf->alp10; double alp11=smf->alp11; double alp12=smf->alp12;
    double alp13=smf->alp13;
    DKERNEL2(drs2,r2);
    DKERNEL3(drs3,r2);
    DKERNEL4(drs4,r2);
    DKERNEL5(drs5,r2);
    DKERNEL6(drs6,r2);
    DKERNEL7(drs7,r2);
    DKERNEL8(drs8,r2);
    DKERNEL9(drs9,r2);
    DKERNEL10(drs10,r2);
    DKERNEL11(drs11,r2);
    DKERNEL12(drs12,r2);
    DKERNEL13(drs13,r2);
    DKERNEL14(drs14,r2);

    double rs15,rs16,rs17,rs18,rs19,rs20,rs21,rs22;
      double drs15,drs16,drs17,drs18,drs19,drs20,drs21,drs22;
    double alp14=smf->alp14; double alp15=smf->alp15;double alp16=smf->alp16;
    double alp17=smf->alp17; double alp18=smf->alp18; double alp19=smf->alp19;
    double alp20=smf->alp20; double alp21=smf->alp21;
    DKERNEL19(drs15,r2);
    DKERNEL20(drs16,r2);
    DKERNEL21(drs17,r2);
    DKERNEL22(drs18,r2);
    DKERNEL23(drs19,r2);
    DKERNEL24(drs20,r2);
    DKERNEL25(drs21,r2);
    DKERNEL26(drs22,r2);
    rs1 = ((((((((((((((((((((((1.0-alp)*rs1+alp*drs2)*(1.0-alp2)+alp2*drs3)*(1.0-alp3)+alp3*drs4)*(1.0-alp4)+alp4*drs5)*(1.0-alp5)+alp5*drs6)*(1.0-alp6)+alp6*drs7)*(1.0-alp7)+alp7*drs8)*(1.0-alp8)+alp8*drs9)*(1.0-alp9)+alp9*drs10)*(1.0-alp10)+alp10*drs11)*(1.0-alp11)+alp11*drs12)*(1.0-alp12)+alp12*drs13)*(1.0-alp13)+alp13*drs14)*(1.0-alp14)+alp14*drs15)*(1.0-alp15)+alp15*drs16)*(1.0-alp16)+alp16*drs17)*(1.0-alp17)+alp17*drs18)*(1.0-alp18)+alp18*drs19)*(1.0-alp19)+alp19*drs20)*(1.0-alp20)+alp20*drs21)*(1.0-alp21)+alp21*drs22);
#endif
    rs1 *= nnList[i].fDist2*fNorm1;

    p->fDivv_Corrector2 -=  2*rs1*q->fMass/p->fDensity;
    rq = rs1*q->fMass/(q->fDensity);
    rp = rs1*p->fMass/(p->fDensity);
    if (TYPEQueryACTIVE(p)) {
      p->fDivv_Corrector -= rq;
    }
    if (TYPEQueryACTIVE(q)) {
      q->fDivv_Corrector -= rp;
    }
  }
  #endif
}

void initISphGradTermsParticle(void *p)
{
#ifdef ISPH
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
    ((PARTICLE *)p)->Igrxx=0.0;
    ((PARTICLE *)p)->Igrxy=0.0;
    ((PARTICLE *)p)->Igrxz=0.0;
    ((PARTICLE *)p)->Igryy=0.0;
    ((PARTICLE *)p)->Igryz=0.0;
    ((PARTICLE *)p)->Igrzz=0.0;
    }
#endif
  }

  void initISphGradTerms(void *p)
  {
#ifdef ISPH
    if (TYPEQueryACTIVE((PARTICLE *) p)) {
      ((PARTICLE *)p)->Igrxx=0.0;
      ((PARTICLE *)p)->Igrxy=0.0;
      ((PARTICLE *)p)->Igrxz=0.0;
      ((PARTICLE *)p)->Igryy=0.0;
      ((PARTICLE *)p)->Igryz=0.0;
      ((PARTICLE *)p)->Igrzz=0.0;
    	}
#endif
  	}

  void combISphGradTerms(void *p1,void *p2)
  {
#ifdef ISPH
    	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
((PARTICLE *)p1)->Igrxx += ((PARTICLE *)p2)->Igrxx;
((PARTICLE *)p1)->Igrxy += ((PARTICLE *)p2)->Igrxy;
((PARTICLE *)p1)->Igrxz += ((PARTICLE *)p2)->Igrxz;
((PARTICLE *)p1)->Igryy += ((PARTICLE *)p2)->Igryy;
((PARTICLE *)p1)->Igryz += ((PARTICLE *)p2)->Igryz;
((PARTICLE *)p1)->Igrzz += ((PARTICLE *)p2)->Igrzz;
  	}
#endif
    }

void ISphGradTerms(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef ISPH
  FLOAT ih2,ih,r2,rs,fDensity,fNorm,fNorm1,vFac,aFac;
  FLOAT dvx,dvy,dvz,dx,dy,dz,divvi,divvj,rs1,rq,rp,qih2;
  double Txx=0.0,Txy=0.0,Txz=0.0,Tyy=0.0,Tyz=0.0,Tzz=0.0;
  PARTICLE *q;
  int i;
  ih2 = 4.0/BALL2(p);
  ih = sqrt(ih2);
  #ifdef SINC
 double nloc = smf->nloc;
 double norm =0.027012593+0.020410827*sqrt(nloc)+0.0037451957*nloc+0.047013839*nloc*sqrt(nloc);
#ifdef COMBOKERNEL
 double nloc2 = smf->nloc2;
 double norm2 =0.027012593+0.020410827*sqrt(nloc2)+0.0037451957*nloc2+0.047013839*nloc2*sqrt(nloc2);
 double fNorm2 = 0.5*norm2*ih2*ih;
 double nloc3 = smf->nloc3;
 double norm3 =0.027012593+0.020410827*sqrt(nloc3)+0.0037451957*nloc3+0.047013839*nloc3*sqrt(nloc3);
 double fNorm3 = 0.5*norm3*ih2*ih;
#endif
 fNorm = 0.5*norm*ih2*ih;
  #else
 fNorm = 0.5*M_1_PI*ih2*ih;
  #endif

  #ifdef CORRENTROPYPART
 p->fDivv_Corrector = 1.0;
  #else
 p->fDivv_Corrector = 1.0;
  #endif

  for (i=0;i<nSmooth;++i) {
    r2 = nnList[i].fDist2*ih2;
    q = nnList[i].pPart;
    if (!TYPEQueryACTIVE(p) && !TYPEQueryACTIVE(q)) continue;
    KERNEL(rs,r2);
 #ifdef COMBOKERNEL
    double rs2,rs3,rs4,rs5,rs6,rs7,rs8,rs9,rs10,rs11,rs12,rs13,rs14;
    double alp=smf->alp; double alp2=smf->alp2; double alp3=smf->alp3;
    double alp4=smf->alp4; double alp5=smf->alp5;double alp6=smf->alp6;
    double alp7=smf->alp7; double alp8=smf->alp8;double alp9=smf->alp9;
    double alp10=smf->alp10; double alp11=smf->alp11; double alp12=smf->alp12;
    double alp13=smf->alp13;
    KERNEL2(rs2,r2);
    KERNEL3(rs3,r2);
    KERNEL4(rs4,r2);
    KERNEL5(rs5,r2);
    KERNEL6(rs6,r2);
    KERNEL7(rs7,r2);
    KERNEL8(rs8,r2);
    KERNEL9(rs9,r2);
    KERNEL10(rs10,r2);
    KERNEL11(rs11,r2);
    KERNEL12(rs12,r2);
    KERNEL13(rs13,r2);
    KERNEL14(rs14,r2);

double rs15,rs16,rs17,rs18,rs19,rs20,rs21,rs22;
      double alp14=smf->alp14; double alp15=smf->alp15;double alp16=smf->alp16;
      double alp17=smf->alp17; double alp18=smf->alp18; double alp19=smf->alp19;
      double alp20=smf->alp20; double alp21=smf->alp21;
      KERNEL19(rs15,r2);
      KERNEL20(rs16,r2);
      KERNEL21(rs17,r2);
      KERNEL22(rs18,r2);
      KERNEL23(rs19,r2);
      KERNEL24(rs20,r2);
      KERNEL25(rs21,r2);
      KERNEL26(rs22,r2);
      rs = ((((((((((((((((((((((1.0-alp)*rs+alp*rs2)*(1.0-alp2)+alp2*rs3)*(1.0-alp3)+alp3*rs4)*(1.0-alp4)+alp4*rs5)*(1.0-alp5)+alp5*rs6)*(1.0-alp6)+alp6*rs7)*(1.0-alp7)+alp7*rs8)*(1.0-alp8)+alp8*rs9)*(1.0-alp9)+alp9*rs10)*(1.0-alp10)+alp10*rs11)*(1.0-alp11)+alp11*rs12)*(1.0-alp12)+alp12*rs13)*(1.0-alp13)+alp13*rs14)*(1.0-alp14)+alp14*rs15)*(1.0-alp15)+alp15*rs16)*(1.0-alp16)+alp16*rs17)*(1.0-alp17)+alp17*rs18)*(1.0-alp18)+alp18*rs19)*(1.0-alp19)+alp19*rs20)*(1.0-alp20)+alp20*rs21)*(1.0-alp21)+alp21*rs22);
 #endif
    rs *= fNorm;
    //rs *= p->fDivv_Corrector;
    rq = rs*q->fMass/(q->fDensity);
    rp = rs*p->fMass/(p->fDensity);
    dx = nnList[i].dx; /* NB: dx = px - qx */
    //printf("1 %g 2 %g",dx,dx*dx*rs);
    dy = nnList[i].dy;
    dz = nnList[i].dz;
    if (TYPEQueryACTIVE(p)) {
   p->Igrxx +=  dx*dx*rq;
   p->Igrxy +=  dx*dy*rq;
   p->Igrxz +=  dx*dz*rq;
   p->Igryy +=  dy*dy*rq;
   p->Igryz +=  dy*dz*rq;
   p->Igrzz +=  dz*dz*rq;
 }
    if (TYPEQueryACTIVE(q)) {
   q->Igrxx +=  dx*dx*rp;
   q->Igrxy +=  dx*dy*rp;
   q->Igrxz +=  dx*dz*rp;
   q->Igryy +=  dy*dy*rp;
   q->Igryz +=  dy*dz*rp;
   q->Igrzz +=  dz*dz*rp;
 }
  }
  #endif
  }


void initICGeneratorParticle(void *p)
{
  #ifdef ICGEN
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
  int i;
  for(i=0;i<3;++i){
  ((PARTICLE *)p)->ICvel[i]=0.0;
  }
  #ifdef CHECKQUALITY
  ((PARTICLE *)p)->Q1=0.0;
  ((PARTICLE *)p)->Q2x=0.0;
  ((PARTICLE *)p)->Q2y=0.0;
  ((PARTICLE *)p)->Q2z=0.0;
  ((PARTICLE *)p)->E0x=0.0;
  ((PARTICLE *)p)->E0y=0.0;
  ((PARTICLE *)p)->E0z=0.0;
  ((PARTICLE *)p)->Q4xx=0.0;
  ((PARTICLE *)p)->Q4yy=0.0;
  ((PARTICLE *)p)->Q4zz=0.0;
  ((PARTICLE *)p)->Q4xy=0.0;
  ((PARTICLE *)p)->Q4yz=0.0;
  ((PARTICLE *)p)->Q4zx=0.0;

  //((PARTICLE *)p)->Veff=0.0;
  #endif
}
  #endif

	}

void initICGenerator(void *p)
{
#ifdef ICGEN
  if (TYPEQueryACTIVE((PARTICLE *) p)) {
  int i;
  for(i=0;i<3;++i){
    ((PARTICLE *)p)->ICvel[i]=0.0;
  }
        #ifdef CHECKQUALITY
  ((PARTICLE *)p)->Q1=0.0;
  ((PARTICLE *)p)->Q2x=0.0;
  ((PARTICLE *)p)->Q2y=0.0;
  ((PARTICLE *)p)->Q2z=0.0;
  ((PARTICLE *)p)->E0x=0.0;
  ((PARTICLE *)p)->E0y=0.0;
  ((PARTICLE *)p)->E0z=0.0;
  ((PARTICLE *)p)->Q4xx=0.0;
  ((PARTICLE *)p)->Q4yy=0.0;
  ((PARTICLE *)p)->Q4zz=0.0;
  ((PARTICLE *)p)->Q4xy=0.0;
  ((PARTICLE *)p)->Q4yz=0.0;
  ((PARTICLE *)p)->Q4zx=0.0;

  //  ((PARTICLE *)p)->Veff=0.0;
        #endif
}
#endif
}

void combICGenerator(void *p1,void *p2)
{
  #ifdef ICGEN
  if (TYPEQueryACTIVE((PARTICLE *) p1)) {
  int i;
  for(i=0;i<3;++i){
    ((PARTICLE *)p1)->ICvel[i] += ((PARTICLE *)p2)->ICvel[i];
  }
#ifdef CHECKQUALITY
	  ((PARTICLE *)p1)->Q1 += ((PARTICLE *)p2)->Q1;
	  ((PARTICLE *)p1)->Q2x += ((PARTICLE *)p2)->Q2x;
	  ((PARTICLE *)p1)->Q2y += ((PARTICLE *)p2)->Q2y;
	  ((PARTICLE *)p1)->Q2z += ((PARTICLE *)p2)->Q2z;
	  ((PARTICLE *)p1)->E0x += ((PARTICLE *)p2)->E0x;
	  ((PARTICLE *)p1)->E0y += ((PARTICLE *)p2)->E0y;
	  ((PARTICLE *)p1)->E0z += ((PARTICLE *)p2)->E0z;
	  ((PARTICLE *)p1)->Q4xx += ((PARTICLE *)p2)->Q4xx;
	  ((PARTICLE *)p1)->Q4yy += ((PARTICLE *)p2)->Q4yy;
	  ((PARTICLE *)p1)->Q4zz += ((PARTICLE *)p2)->Q4zz;
          ((PARTICLE *)p1)->Q4xy += ((PARTICLE *)p2)->Q4xy;
          ((PARTICLE *)p1)->Q4yz += ((PARTICLE *)p2)->Q4yz;
          ((PARTICLE *)p1)->Q4zx += ((PARTICLE *)p2)->Q4zx;
	  //((PARTICLE *)p1)->Veff += ((PARTICLE *)p2)->Veff;
}
#endif
#endif
	}


void ICGenerator(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
  
#ifdef ICGEN
  FLOAT ih2,ih,ph,r2,rs,drs,fDensity,fNorm,fNorm1,vFac,aFac,divvi,divvj,dxterm,dyterm,dzterm,nrs,pPeff,qPeff,pP,qP;
  FLOAT dx,dy,dz,rq,rp,igDensity2p,igDensity2q,ih2q,ihq,qh,rpartp,rpartq,rho0p,rho0q,igrho0p,igrho0q;
  FLOAT dxtermq,dytermq,dztermq,divvcorrq,divvcorrp,ttres,Rinside,smhrat;
  PARTICLE *q;
  int i;
  //  printf("\n iorder %d STARTING \n ",p->iOrder);
  ih2 = 4.0/BALL2(p);
  ih = sqrt(ih2);
  ph = 1/ih;
  #ifdef SINC
  double nloc = smf->nloc;
  double norm =0.027012593+0.020410827*sqrt(nloc)+0.0037451957*nloc+0.047013839*nloc*sqrt(nloc);
    #ifdef COMBOKERNEL
  double nloc2 = smf->nloc2;
  double norm2 =0.027012593+0.020410827*sqrt(nloc2)+0.0037451957*nloc2+0.047013839*nloc2*sqrt(nloc2);
  double fNorm2 = 0.5*norm2*ih2*ih;
  double fNorm12 = fNorm2*ih2;
  double nloc3 = smf->nloc3;
  double norm3 =0.027012593+0.020410827*sqrt(nloc3)+0.0037451957*nloc3+0.047013839*nloc3*sqrt(nloc3);
  double fNorm3 = 0.5*norm3*ih2*ih;
  double fNorm13 = fNorm3*ih2;
#endif
fNorm = 0.5*norm*ih2*ih;
  #else
  fNorm = 0.5*M_1_PI*ih2*ih;
  #endif
  fNorm1 = fNorm*ih2;

  //smf->densprofile 1 const dens 2 discontinuity dens (set by smf->densdiff) smf->densR sets the range of the high density region, smoothed dens profile (smooth region set by densRsmooth) 4 accretiondisk 5 ...
  //smf->densdir= 1 x, 2 y, 3 z, 4 r, 5 rcyl(x,y)
  if(smf->ICdensdir < 4){
      rpartp = sqrt(p->r[smf->ICdensdir-1]*p->r[smf->ICdensdir-1]);
    }
    else if(smf->ICdensdir==4){
      rpartp = sqrt(p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2]);
    }
    else if(smf->ICdensdir==5){
      rpartp = sqrt(p->r[0]*p->r[0]+p->r[1]*p->r[1]);
    }
    else if(smf->ICdensdir==6){
      rpartp = p->r[0];
    }
  
  if(smf->ICdensprofile==1){
    rho0p=smf->ICdensouter;
    }
  if(smf->ICdensprofile==2){
    if(rpartp > smf->ICdensR){rho0p=smf->ICdensouter;}
    else{rho0p=smf->ICdensinner;}
    }
  if(smf->ICdensprofile==3){
    if(smf->ICdensdir > 3){
      rho0p = smf->ICdensouter+(smf->ICdensinner-smf->ICdensouter)*1.0/(1.0+exp(-(smf->ICdensR-rpartp)/smf->ICdensRsmooth));
    }
    else{
      rho0p = smf->ICdensouter+(smf->ICdensinner-smf->ICdensouter)*0.5*(tanh((p->r[smf->ICdensdir-1]+smf->ICdensR)/smf->ICdensRsmooth)-tanh((p->r[smf->ICdensdir-1]-smf->ICdensR)/smf->ICdensRsmooth));
    }
  }
   if(smf->ICdensprofile >= 4){
     FLOAT zp = sqrt(p->r[2]*p->r[2]);
     FLOAT Hdisk = smf->ICdensRsmooth*rpartp;
     FLOAT qdisk = 0.5;
     FLOAT sigmadisk= pow(rpartp,-qdisk);
     //FLOAT mindens = 0.25*pow(smf->ICdensR,-qdisk);
     FLOAT mindens = smf->ICdensouter;
     FLOAT smooth = smf->ICdensRsmooth;
     smhrat = 1.0; //smoothing of unresolved density gradients, for example, if u want a 0.03 scaleheight disk and smoothing length is 1 it changes the scaleheight to be = smhrat*smoothing length
     Rinside = 0.5;
     if(smf->ICdensprofile==5 || smf->ICdensprofile==7){
       smooth = ( rpartp <= Rinside && ph*smhrat > smf->ICdensRsmooth ? ph*smhrat : smf->ICdensRsmooth );
       sigmadisk*=1.0/(1.0+exp(-(rpartp-Rinside)/smooth));
     }
     if(smf->ICdensprofile==6 || smf->ICdensprofile==7){
       smooth = ( rpartp >= smf->ICdensR && ph*smhrat > smf->ICdensRsmooth ? ph*smhrat : smf->ICdensRsmooth );
       sigmadisk*=1.0/(1.0+exp(-(smf->ICdensR-rpartp)/smooth));
	  }
     //resolved scaleheight, though leads to quite significant changing target density field, so probably bad idea here.
     //smooth = ( ph*smhrat*2.0 > Hdisk ? ph*smhrat*2.0 : Hdisk );
     smooth = Hdisk;
     sigmadisk *= 1.0/(smooth*sqrt(2*M_PI));
     //We dampen the density gradient in the accretion disk by changing it to a pow(-2) after 3 scale heights(representing a density difference of around 100 from centre)
     //There is probably a better way to do this... Tried to smooth similar to inner and outer disk but under estimates inner region by quite a lot in that case. 
     //sigmadisk *= 0.5*pow(10,-8)*exp(-zp*zp/(2*smooth*smooth));
     // Somehow ttres should depend on resolution though.... maybe needs to be an adjustable parameter
     //Though at ttres = 3.0 the integral_-3^3 rho dz = 0.9973, for ttres=4 0.999937, ttres=2 0.9545
     //ttres = 2.0 full integral with pow 2 decay 1.17046, pow 4 decay 1.02649 
     //ttres = 3.0 full integral with pow 2 decay 1.02389, pow 4 decay 1.00616
     //ttres = 4.0 full integral with pow 2 decay 1.00101, pow 4 decay 1.00029
     //ttres = 3.0 and pow 3 full integral 1.0106, seems enough
     ttres = 3.0;
     sigmadisk *= 0.75*pow(10,-8)*( zp > ttres*Hdisk ? exp(-ttres*ttres/2.0)*(pow(zp/(ttres*Hdisk),-3.0)) : exp(-zp*zp/(2*Hdisk*Hdisk)) );

     if (isnan(sigmadisk)) {
    printf("SIGMADISK IS NAN\n");
    sigmadisk = mindens;
} else if (sigmadisk <= mindens) {
       //printf("SIGMADISK is below threshold\n");
    sigmadisk = mindens;
}
     rho0p = sigmadisk; 
   }

   //  printf(" nsmooth %d ph %f u %f upred %f dens %g rho0p %g P %g \n", nSmooth,ph,p->u,p->uPred,p->fDensity,rho0p,pow(p->fDensity/(rho0p),4.0));
if (p->iOrder==0){printf(" cs %f P %f nsmooth %d ph %f u %f upred %f dens %f rho0p %f \n", p->PoverRho2*p->fDensity,p->PoverRho2*p->fDensity*p->fDensity,nSmooth,ph,p->u,p->uPred,p->fDensity,rho0p) \
;}


   
 // The power here determines the balance between getting the right density and the zero order errors (higher power weights density higher than getting a glassy/order distribution).
 // For stiff EOS a higher power is useful to fit the density more closely.
 p->fSoft = pow(p->fMass/rho0p,1./3.);
 p->fSoft0 = pow(p->fMass/rho0p,1./3.);
 //double alprho = max(2.0,1.0/(0.144*abs(log(p->fDensity/rho0p))+0.1));
 double alprho = smf->ICrhopow;
 
 pP = pow(p->fDensity/(rho0p),alprho);
 p->u=p->uPred=rho0p;

 if (p->iOrder==0){printf(" cs %f P %f nsmooth %d ph %f u %f upred %f dens %f rho0p %f \n", p->PoverRho2*p->fDensity,p->PoverRho2*p->fDensity*p->fDensity,nSmooth,ph,p->u,p->uPred,p->fDensity,rho0p) \
;}



#ifdef CORRENTROPYPART
  divvcorrp = (p->fDivv_Corrector != 0 ? p->fDivv_Corrector2/p->fDivv_Corrector : 1.0);
#else
  divvcorrp = 1.0;
#endif
  
  for (i=0;i<nSmooth;++i) {
    r2 = nnList[i].fDist2*ih2;
    q = nnList[i].pPart;
    ih2q = 4.0/BALL2(q);
    ihq = sqrt(ih2q);
    qh = 1/ihq;
    DKERNEL(drs,r2);
    KERNEL(rs,r2)
#ifdef COMBOKERNEL
      double rs2,rs3,rs4,rs5,rs6,rs7,rs8,rs9,rs10,rs11,rs12,rs13,rs14;
      double drs2,drs3,drs4,drs5,drs6,drs7,drs8,drs9,drs10,drs11,drs12,drs13,drs14;
      double alp=smf->alp; double alp2=smf->alp2; double alp3=smf->alp3;
      double alp4=smf->alp4; double alp5=smf->alp5;double alp6=smf->alp6;
      double alp7=smf->alp7; double alp8=smf->alp8;double alp9=smf->alp9;
      double alp10=smf->alp10; double alp11=smf->alp11; double alp12=smf->alp12;
      double alp13=smf->alp13;
      KERNEL2(rs2,r2);
      KERNEL3(rs3,r2);
      KERNEL4(rs4,r2);
      KERNEL5(rs5,r2);
      KERNEL6(rs6,r2);
      KERNEL7(rs7,r2);
      KERNEL8(rs8,r2);
      KERNEL9(rs9,r2);
      KERNEL10(rs10,r2);
      KERNEL11(rs11,r2);
      KERNEL12(rs12,r2);
      KERNEL13(rs13,r2);
      KERNEL14(rs14,r2);
      DKERNEL2(drs2,r2);
      DKERNEL3(drs3,r2);
      DKERNEL4(drs4,r2);
      DKERNEL5(drs5,r2);
      DKERNEL6(drs6,r2);
      DKERNEL7(drs7,r2);
      DKERNEL8(drs8,r2);
      DKERNEL9(drs9,r2);
      DKERNEL10(drs10,r2);
      DKERNEL11(drs11,r2);
      DKERNEL12(drs12,r2);
      DKERNEL13(drs13,r2);
      DKERNEL14(drs14,r2);

double rs15,rs16,rs17,rs18,rs19,rs20,rs21,rs22;
double drs15,drs16,drs17,drs18,drs19,drs20,drs21,drs22;
      double alp14=smf->alp14; double alp15=smf->alp15;double alp16=smf->alp16;
      double alp17=smf->alp17; double alp18=smf->alp18; double alp19=smf->alp19;
      double alp20=smf->alp20; double alp21=smf->alp21;
      KERNEL19(rs15,r2);
      KERNEL20(rs16,r2);
      KERNEL21(rs17,r2);
      KERNEL22(rs18,r2);
      KERNEL23(rs19,r2);
      KERNEL24(rs20,r2);
      KERNEL25(rs21,r2);
      KERNEL26(rs22,r2);
      DKERNEL19(drs15,r2);
      DKERNEL20(drs16,r2);
      DKERNEL21(drs17,r2);
      DKERNEL22(drs18,r2);
      DKERNEL23(drs19,r2);
      DKERNEL24(drs20,r2);
      DKERNEL25(drs21,r2);
      DKERNEL26(drs22,r2);
      rs = ((((((((((((((((((((((1.0-alp)*rs+alp*rs2)*(1.0-alp2)+alp2*rs3)*(1.0-alp3)+alp3*rs4)*(1.0-alp4)+alp4*rs5)*(1.0-alp5)+alp5*rs6)*(1.0-alp6)+alp6*rs7)*(1.0-alp7)+alp7*rs8)*(1.0-alp8)+alp8*rs9)*(1.0-alp9)+alp9*rs10)*(1.0-alp10)+alp10*rs11)*(1.0-alp11)+alp11*rs12)*(1.0-alp12)+alp12*rs13)*(1.0-alp13)+alp13*rs14)*(1.0-alp14)+alp14*rs15)*(1.0-alp15)+alp15*rs16)*(1.0-alp16)+alp16*rs17)*(1.0-alp17)+alp17*rs18)*(1.0-alp18)+alp18*rs19)*(1.0-alp19)+alp19*rs20)*(1.0-alp20)+alp20*rs21)*(1.0-alp21)+alp21*rs22);
      drs = ((((((((((((((((((((((1.0-alp)*drs+alp*drs2)*(1.0-alp2)+alp2*drs3)*(1.0-alp3)+alp3*drs4)*(1.0-alp4)+alp4*drs5)*(1.0-alp5)+alp5*drs6)*(1.0-alp6)+alp6*drs7)*(1.0-alp7)+alp7*drs8)*(1.0-alp8)+alp8*drs9)*(1.0-alp9)+alp9*drs10)*(1.0-alp10)+alp10*drs11)*(1.0-alp11)+alp11*drs12)*(1.0-alp12)+alp12*drs13)*(1.0-alp13)+alp13*drs14)*(1.0-alp14)+alp14*drs15)*(1.0-alp15)+alp15*drs16)*(1.0-alp16)+alp16*drs17)*(1.0-alp17)+alp17*drs18)*(1.0-alp18)+alp18*drs19)*(1.0-alp19)+alp19*drs20)*(1.0-alp20)+alp20*drs21)*(1.0-alp21)+alp21*drs22);
#endif
    rs *= fNorm;
    drs *= fNorm1;
    dx = nnList[i].dx;
    dy = nnList[i].dy;
    dz = nnList[i].dz;

    #ifdef ISPH
    dxterm = -(p->Igrxx*dx+p->Igrxy*dy+p->Igrxz*dz);
    dyterm = -(p->Igrxy*dx+p->Igryy*dy+p->Igryz*dz);
    dzterm = -(p->Igrxz*dx+p->Igryz*dy+p->Igrzz*dz);
    dxtermq = -(q->Igrxx*dx+q->Igrxy*dy+q->Igrxz*dz);
    dytermq = -(q->Igrxy*dx+q->Igryy*dy+q->Igryz*dz);
    dztermq = -(q->Igrxz*dx+q->Igryz*dy+q->Igrzz*dz);
    rq=rs*q->fMass;
    rp=rs*p->fMass;
    #else
    dxterm = dx; 
    dyterm = dy;
    dzterm = dz;
    dxtermq = dx; 
    dytermq = dy;
    dztermq = dz;
    rq=drs*q->fMass;
    rp=drs*p->fMass;
    #endif
    
#ifdef CORRENTROPYPART
      divvcorrq = (q->fDivv_Corrector != 0 ? q->fDivv_Corrector2/q->fDivv_Corrector : 1.0);
#else
      divvcorrq = 1.0;
#endif
    
    //if(p->iOrder==0) printf("q->iOrder %d rq %f rp %f dxterm %f dyterm %f dzterm %f \n",q->iOrder,rq,rp,dxterm,dyterm,dzterm);
#ifdef GDFORCE
      igDensity2p=1/(p->fDensity*q->fDensity);
      igDensity2q=igDensity2p;
#else
    igDensity2p=1/(p->fDensity*p->fDensity);
    igDensity2q=1/(q->fDensity*q->fDensity);    
#endif
    //V0 = p->fMass/(1.5-rpart*2.0);
    if(smf->ICdensdir < 4){
      rpartq = sqrt(q->r[smf->ICdensdir-1]*q->r[smf->ICdensdir-1]);
    }
    else if(smf->ICdensdir == 4){
      rpartq = sqrt(q->r[0]*q->r[0]+q->r[1]*q->r[1]+q->r[2]*q->r[2]);
    }
    else if(smf->ICdensdir == 5){
      rpartq = sqrt(q->r[0]*q->r[0]+q->r[1]*q->r[1]);
    }
    else if(smf->ICdensdir==6){
      rpartq = q->r[0];
    }
      
    if(smf->ICdensprofile==1){
      rho0q=smf->ICdensouter;
    }
    if(smf->ICdensprofile==2){
      if(rpartq > smf->ICdensR){rho0q=smf->ICdensouter;}
      else{rho0q=smf->ICdensinner;}
    }
    if(smf->ICdensprofile==3){
      if(smf->ICdensdir > 3){
	rho0q = smf->ICdensouter+(smf->ICdensinner-smf->ICdensouter)*1.0/(1.0+exp(-(smf->ICdensR-rpartq)/smf->ICdensRsmooth));
	}
    else{
      rho0q = smf->ICdensouter+(smf->ICdensinner-smf->ICdensouter)*0.5*(tanh((q->r[smf->ICdensdir-1]+smf->ICdensR)/smf->ICdensRsmooth)-tanh((q->r[smf->ICdensdir-1]-smf->ICdensR)/smf->ICdensRsmooth));
    }
    }
   if(smf->ICdensprofile >= 4){
     FLOAT zq = sqrt(q->r[2]*q->r[2]);
     FLOAT Hdisk = smf->ICdensRsmooth*rpartq;
     FLOAT qdisk = 0.5;
     FLOAT sigmadisk= pow(rpartq,-qdisk);
     //FLOAT mindens = 0.25*pow(smf->ICdensR,-qdisk);
     FLOAT mindens = smf->ICdensouter;
     FLOAT smooth = smf->ICdensRsmooth;
     if(smf->ICdensprofile==5 || smf->ICdensprofile==7){
       smooth = ( rpartq <= Rinside && qh*smhrat > smf->ICdensRsmooth ? qh*smhrat : smf->ICdensRsmooth );
       sigmadisk*=1.0/(1.0+exp(-(rpartq-Rinside)/smooth));
     }
     if(smf->ICdensprofile==6 || smf->ICdensprofile==7){
       smooth = ( rpartq >= smf->ICdensR && qh*smhrat > smf->ICdensRsmooth ? qh*smhrat : smf->ICdensRsmooth );
       sigmadisk*=1.0/(1.0+exp(-(smf->ICdensR-rpartq)/smooth));
          }
     //resolved scaleheight
     //smooth = ( qh*smhrat*2.0 > Hdisk ? qh*smhrat*2.0 : Hdisk );
     smooth = Hdisk;
     sigmadisk *= 1.0/(smooth*sqrt(2*M_PI));
     //sigmadisk *= 0.5*pow(10,-8)*exp(-zq*zq/(2*smooth*smooth));
     sigmadisk *= 0.75*pow(10,-8)*( zq > ttres*Hdisk ? exp(-ttres*ttres/2.0)*(pow(zq/(ttres*Hdisk),-3.0)) : exp(-zq*zq/(2*Hdisk*Hdisk)) );
     
     if (isnan(sigmadisk)) {
     printf("SIGMADISK IS NAN\n");
    sigmadisk = mindens;
} else if (sigmadisk <= mindens) {
       //printf("SIGMADISK is below threshold\n");
    sigmadisk = mindens;
}
     rho0q = sigmadisk;
   }


   //alprho = max(2.0,1.0/(0.144*abs(log(q->fDensity/rho0q))+0.1));
   qP = pow(q->fDensity/(rho0q),alprho);
   q->u=q->uPred=rho0q;
   pPeff=pP;
   qPeff=qP;
    
   igrho0p=pPeff*igDensity2p;
   igrho0q=qPeff*igDensity2q;

    
    if (TYPEQueryACTIVE(p)) {   
    p->ICvel[0] -= (igrho0p*dxterm*divvcorrp+igrho0q*dxtermq*divvcorrq)*rq;
    p->ICvel[1] -= (igrho0p*dyterm*divvcorrp+igrho0q*dytermq*divvcorrq)*rq;
    p->ICvel[2] -= (igrho0p*dzterm*divvcorrp+igrho0q*dztermq*divvcorrq)*rq;
    }
    if (TYPEQueryACTIVE(q)) {
    q->ICvel[0] += (igrho0p*dxterm*divvcorrp+igrho0q*dxtermq*divvcorrq)*rp;
    q->ICvel[1] += (igrho0p*dyterm*divvcorrp+igrho0q*dytermq*divvcorrq)*rp;
    q->ICvel[2] += (igrho0p*dzterm*divvcorrp+igrho0q*dztermq*divvcorrq)*rp;
    }
    
#ifdef CHECKQUALITY
    
    p->Q1  += q->fMass/q->fDensity*rs*2.0;
    p->Q2x -= q->fMass/q->fDensity*dx*rs*2.0;
    p->Q2y -= q->fMass/q->fDensity*dy*rs*2.0;
    p->Q2z -= q->fMass/q->fDensity*dz*rs*2.0;
    p->E0x += p->fDensity*ph*(igDensity2p*dxterm*divvcorrp+igDensity2q*dxtermq*divvcorrq)*rq;
    p->E0y += p->fDensity*ph*(igDensity2p*dyterm*divvcorrp+igDensity2q*dytermq*divvcorrq)*rq;
    p->E0z += p->fDensity*ph*(igDensity2p*dzterm*divvcorrp+igDensity2q*dztermq*divvcorrq)*rq;

    p->Q4xx -= rq*p->fDensity*igDensity2q*dxterm*dx*divvcorrp;
    p->Q4yy -= rq*p->fDensity*igDensity2q*dyterm*dy*divvcorrp;
    p->Q4zz -= rq*p->fDensity*igDensity2q*dzterm*dz*divvcorrp;
    p->Q4xy -= rq*p->fDensity*igDensity2q*dxterm*dy*divvcorrp;
    p->Q4yz -= rq*p->fDensity*igDensity2q*dyterm*dz*divvcorrp;
    p->Q4zx -= rq*p->fDensity*igDensity2q*dzterm*dx*divvcorrp;
    
    q->E0x -= q->fDensity*qh*(igDensity2p*dxterm*divvcorrp+igDensity2q*dxtermq*divvcorrq)*rp;
    q->E0y -= q->fDensity*qh*(igDensity2p*dyterm*divvcorrp+igDensity2q*dytermq*divvcorrq)*rp;
    q->E0z -= q->fDensity*qh*(igDensity2p*dzterm*divvcorrp+igDensity2q*dztermq*divvcorrq)*rp;
    q->Q4xx -= rp*q->fDensity*igDensity2p*dxtermq*dx*divvcorrq;
    q->Q4yy -= rp*q->fDensity*igDensity2p*dytermq*dy*divvcorrq;
    q->Q4zz -= rp*q->fDensity*igDensity2p*dztermq*dz*divvcorrq;
    q->Q4xy -= rp*q->fDensity*igDensity2p*dxtermq*dy*divvcorrq;
    q->Q4yz -= rp*q->fDensity*igDensity2p*dytermq*dz*divvcorrq;
    q->Q4zx -= rp*q->fDensity*igDensity2p*dztermq*dx*divvcorrq;
      
    //Calculate effecctive V
    //p->Veff -= (1.0/3.0)*rs*(dxterm*dx+dyterm*dy+dzterm*dz)*(p->fMass*q->fMass/(p->fDensity*q->fDensity));
    //q->Veff -= (1.0/3.0)*rs*(dxterm*dx+dyterm*dy+dzterm*dz)*(p->fMass*q->fMass/(p->fDensity*q->fDensity));

    //p->Veff -= (1.0/3.0)*drs*(dx*dx+dy*dy+dz*dz)*(p->fMass*q->fMass/(p->fDensity*q->fDensity));
    //q->Veff -= (1.0/3.0)*drs*(dx*dx+dy*dy+dz*dz)*(p->fMass*q->fMass/(p->fDensity*q->fDensity));
#endif
   
}

  if(p->iOrder==0) printf(" q1 %g q2x %g q2y %g q2z %g e0x %g e0y %g e0z %g q4x %g q4y %g q4z %g ",p->Q1,p->Q2x,p->Q2y,p->Q2z,p->E0x,p->E0y,p->E0z,p->Q4xx,p->Q4yy,p->Q4zz);

  //    printf("\n iorder %d ENDING \n",p->iOrder);
  if(smf->bStepZero==1){
    p->ICfmax = 0.0;
}
  #endif
  }

void initDenDVDX(void *p)
{
	}

void combDenDVDX(void *p1,void *p2)
{
    ((PARTICLE *)p1)->iActive |= (((PARTICLE *)p2)->iActive & TYPE_MASK);
	}

void postDenDVDX(PARTICLE *p, SMF *smf) {
#ifdef CULLENDEHNEN2
#ifdef CD_RDVDSONSFULL
    p->divv_old = p->dvdsonSFull;
#else
#ifdef CD_RDVDS
    p->divv_old = p->dvds;
#else
    p->divv_old = p->divv;
#endif
#endif
#endif
    }

#ifndef FITDVDX
/* Gather only version */
void DenDVDX(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	FLOAT ih2,ih,r2,rs,rs1,fDensity,fNorm,fNorm1,vFac;
    FLOAT dvxdx1 , dvxdy1 , dvxdz1, dvydx1 , dvydy1 , dvydz1, dvzdx1 , dvzdy1 , dvzdz1;
	FLOAT dvxdx , dvxdy , dvxdz, dvydx , dvydy , dvydz, dvzdx , dvzdy , dvzdz; /* comoving shear tensor */
	FLOAT dvx,dvy,dvz,dx,dy,dz,trace;
#ifdef DVDX_CD
    FLOAT dxdx=0,dxdy=0,dxdz=0,dydy=0,dydz=0,dzdz=0;
#else
    FLOAT divvnorm=0;
#endif
    double Hcorr = smf->H;
    FLOAT rgux,rguy,rguz;
	FLOAT grx,gry,grz,gnorm,dvds,dvdr,c;
    FLOAT R_CD=0,R_CDA=0,R_CDN=0,OneMinusR_CD,pdivv_old,pdvds_old,vSigMax=0,vmx=0,vmy=0,vmz=0,alphaNoise=0; // Cullen & Dehnen 2010
    FLOAT fDensity_old;
#if defined (DENSITYU) || defined(RTDENSITY) || defined(THERMALCOND)
	FLOAT fDensityU = 0;
#endif
	PARTICLE *q;
	int i;
	unsigned int qiActive;
#ifdef SFBOUND
    double fSigma2 = 0;
#endif
    int Rkern = 2.0;
    //Calculate effective volume
    //if(p->Veff>0.0){printf("Density %f new: %f",p->fDensity,p->fMass/p->Veff);}

	ih2 = 4.0/BALL2(p);
	ih = sqrt(ih2);
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */
  #ifdef RHOITERATIVE
    	qiActive = 0;
	double phmax,tolh,hdiff,hfact3,hfact3min,hdens,hstart,ph,phmin;
      double rhoit,phold,f,df,df2,dWdh,drs;
      int Rkern,conv,convmax,converged,nSmoothNew;
    	phmax = 1/ih; tolh = 0.0001; hdiff = 1.0; Rkern = 2;
    	hfact3 = 3.0*smf->nSmoothMean*M_1_PI/(Rkern*Rkern*Rkern*4.0);
    	hdens = pow(p->fMass*hfact3/p->fDensity,1.0/3.0);
    	//p->hPred = p->hPred == 0.0 ? hdens : p->hPred;
    	//double hstart = p->hPred > phmax || p->hPred <= 0.0 ? phmax : p->hPred;
  	hstart = hdens > phmax || hdens <= 0.0 ? phmax : hdens;
	hfact3min = 3.0*50*M_1_PI/(Rkern*Rkern*Rkern*4.0);
	phmin = pow(p->fMass*hfact3min/p->fDensity,1.0/3.0);
    	ph=hstart; conv = 0; convmax = 20; converged = 0;
    	ih=1/ph; ih2=ih*ih;
    	do {
  	  phold=ph;
    	  fNorm = M_1_PI*ih2*ih;
    	  fDensity = 0.0; dWdh = 0.0; nSmoothNew=0;
    	  rhoit = p->fMass*hfact3*ih2*ih;
    	  for(i=0;i<nSmooth;++i){
      	    r2 = nnList[i].fDist2*ih2;
          if(r2>=(Rkern*Rkern)){continue;}
      	    q = nnList[i].pPart;
          if (TYPETest(p,TYPE_ACTIVE)) TYPESet(q,TYPE_NbrOfACTIVE); /* important for SPH */
    	    qiActive |= q->iActive;
    	    DKERNEL(drs,r2);
    	    KERNEL(rs,r2);
    	    fDensity += rs*q->fMass*fNorm;
    	    dWdh -= M_1_PI*q->fMass*ih2*ih2*(drs*r2+3*rs);
    	    nSmoothNew=nSmoothNew+1;
    	  }
	  if(p->iOrder==0 || p->iOrder==6154){ printf("hstart %f ph %f phmax %f phmin %f f %f df %f nSmoothnew \
%d order %d nsmooth %d \n", hstart ,ph,phmax,phmin,f,df,nSmoothNew,p->iOrder,nSmooth);}
    	  f=fDensity-rhoit;
          df = 3*rhoit*ih*(1+ph/(3*rhoit)*dWdh);
    	  ph = ph-f/df;
          ph = ph > phmax ? phmax : ph;
          ph = ph < phmin ? phmax : ph;
	  if (nSmoothNew<40){printf("id %d nsmoothnew %d",p->iOrder,nSmoothNew);}
	  ih=1/ph; ih2=ih*ih;
    	  hdiff=sqrt((ph-phold)*(ph-phold))/hstart;
    	  conv += 1;
    	  converged = hdiff<tolh ? 1 : 0;
    	}while (conv<convmax && converged == 0);
    	if(converged==0){printf("NOT CONVERGED conv: %d h: %f hstart %f f: %f df: %f hdiff %.10f \n",conv,ph,hstart,f,df,hdiff);}
    	if(p->iOrder==0 || p->iOrder==6154){ printf("hstart %f ph %f phmax %f df %f nSmoothnew %d order %d \n", hstart ,ph,phmax,df,nSmoothNew,p->iOrder);}
    	dWdh=0.0; nSmoothNew=0; p->fBall2=4.0/ih2;
    #endif
#ifdef SETMAXH3
	if(p->fBall2 > 0.249*1.415*1.415){
	  p->fBall2=0.249*1.415*1.415;
	  ih2 = 4.0/BALL2(p);
	  ih = sqrt(ih2);
	  ph=1/ih;
          if(p->iOrder==0 || p->iOrder==6154){ printf("ph %f\n",ph);}}
#endif

#ifdef SPH1D
	fNorm = (2./3.)*ih;
#else
#ifdef SINC
	double nloc = smf->nloc;
	double norm =0.027012593+0.020410827*sqrt(nloc)+0.0037451957*nloc+0.047013839*nloc*sqrt(nloc);
#ifdef COMBOKERNEL
	double nloc2 = smf->nloc2;
	double norm2 =0.027012593+0.020410827*sqrt(nloc2)+0.0037451957*nloc2+0.047013839*nloc2*sqrt(nloc2);
	double fNorm2 = norm2*ih2*ih;
        double nloc3 = smf->nloc3;
        double norm3 =0.027012593+0.020410827*sqrt(nloc3)+0.0037451957*nloc3+0.047013839*nloc3*sqrt(nloc3);
        double fNorm3 = norm3*ih2*ih;
#endif	
fNorm = norm*ih2*ih;
  #else
	fNorm = M_1_PI*ih2*ih;
  #endif
#endif
	fNorm1 = fNorm*ih2;
	fDensity = 0.0;
	dvxdx1 = 0; dvxdy1 = 0; dvxdz1= 0;
	dvydx1 = 0; dvydy1 = 0; dvydz1= 0;
	dvzdx1 = 0; dvzdy1 = 0; dvzdz1= 0;

	grx = 0; gry = 0; grz = 0;
	rgux = 0; rguy = 0; rguz = 0;

	qiActive = 0;

#ifdef SHEARTEST2
    int nCount = 0; // SLIDING_PATCH_DEBUG
#endif
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		if(r2>=(Rkern*Rkern)){continue;}
	      
		q = nnList[i].pPart;
		if (TYPETest(p,TYPE_ACTIVE)) TYPESet(q,TYPE_NbrOfACTIVE); /* important for SPH */
		qiActive |= q->iActive;
#ifdef DKDENSITY
		DKERNEL(rs,r2);
		rs = -(1./3.)*r2*rs;
#else
		KERNEL(rs,r2);
#endif
#ifdef COMBOKERNEL
    double rs2,rs3,rs4,rs5,rs6,rs7,rs8,rs9,rs10,rs11,rs12,rs13,rs14;
    double alp=smf->alp; double alp2=smf->alp2; double alp3=smf->alp3;
    double alp4=smf->alp4; double alp5=smf->alp5;double alp6=smf->alp6;
    double alp7=smf->alp7; double alp8=smf->alp8;double alp9=smf->alp9;
    double alp10=smf->alp10; double alp11=smf->alp11; double alp12=smf->alp12;
    double alp13=smf->alp13;
    KERNEL2(rs2,r2);
    KERNEL3(rs3,r2);
    KERNEL4(rs4,r2);
    KERNEL5(rs5,r2);
    KERNEL6(rs6,r2);
    KERNEL7(rs7,r2);
    KERNEL8(rs8,r2);
    KERNEL9(rs9,r2);
    KERNEL10(rs10,r2);
    KERNEL11(rs11,r2);
    KERNEL12(rs12,r2);
    KERNEL13(rs13,r2);
    KERNEL14(rs14,r2);
double rs15,rs16,rs17,rs18,rs19,rs20,rs21,rs22;
      double alp14=smf->alp14; double alp15=smf->alp15;double alp16=smf->alp16;
      double alp17=smf->alp17; double alp18=smf->alp18; double alp19=smf->alp19;
      double alp20=smf->alp20; double alp21=smf->alp21;
      KERNEL19(rs15,r2);
      KERNEL20(rs16,r2);
      KERNEL21(rs17,r2);
      KERNEL22(rs18,r2);
      KERNEL23(rs19,r2);
      KERNEL24(rs20,r2);
      KERNEL25(rs21,r2);
      KERNEL26(rs22,r2);
      rs = ((((((((((((((((((((((1.0-alp)*rs+alp*rs2)*(1.0-alp2)+alp2*rs3)*(1.0-alp3)+alp3*rs4)*(1.0-alp4)+alp4*rs5)*(1.0-alp5)+alp5*rs6)*(1.0-alp6)+alp6*rs7)*(1.0-alp7)+alp7*rs8)*(1.0-alp8)+alp8*rs9)*(1.0-alp9)+alp9*rs10)*(1.0-alp10)+alp10*rs11)*(1.0-alp11)+alp11*rs12)*(1.0-alp12)+alp12*rs13)*(1.0-alp13)+alp13*rs14)*(1.0-alp14)+alp14*rs15)*(1.0-alp15)+alp15*rs16)*(1.0-alp16)+alp16*rs17)*(1.0-alp17)+alp17*rs18)*(1.0-alp18)+alp18*rs19)*(1.0-alp19)+alp19*rs20)*(1.0-alp20)+alp20*rs21)*(1.0-alp21)+alp21*rs22);
#endif

fDensity += rs*q->fMass;
#ifdef RHOITERATIVE
  	DKERNEL(drs,r2);
  	dWdh -= M_1_PI*q->fMass*ih2*ih2*(drs*r2+3*rs);
  	nSmoothNew += 1;
#endif
		DKERNEL(rs1,r2); /* rs1 is negative */
		rs1 *= q->fMass;
		dx = nnList[i].dx; /* NB: dx = px - qx */
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = (-p->vPred[0] + q->vPred[0])*vFac;
		dvy = (-p->vPred[1] + q->vPred[1])*vFac;
		dvz = (-p->vPred[2] + q->vPred[2])*vFac;

        // Convention here dvdx = vxq-vxp, dx = xp-xq  but rs < 0 so dvdx correct sign
		dvxdx1 += dvx*dx*rs1;
		dvxdy1 += dvx*dy*rs1;
		dvxdz1 += dvx*dz*rs1;
		dvydx1 += dvy*dx*rs1;
		dvydy1 += dvy*dy*rs1;
		dvydz1 += dvy*dz*rs1;
		dvzdx1 += dvz*dx*rs1;
		dvzdy1 += dvz*dy*rs1;
		dvzdz1 += dvz*dz*rs1;
		divvnorm += (dx*dx+dy*dy+dz*dz)*rs1;

		}
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);

    fDensity_old = p->fDensity;
    p->fDensity = fDensity*fNorm;

    divvnorm = (divvnorm != 0 ? 3/fabs(divvnorm) : 0); /* keep Norm positive consistent w/ std 1/rho norm */
    dvxdx = divvnorm*dvxdx1; dvydx = divvnorm*dvydx1; dvzdx = divvnorm*dvzdx1;
    dvxdy = divvnorm*dvxdy1; dvydy = divvnorm*dvydy1; dvzdy = divvnorm*dvzdy1;
    dvxdz = divvnorm*dvxdz1; dvydz = divvnorm*dvydz1; dvzdz = divvnorm*dvzdz1;

#ifdef CD_DEBUG
    p->dvxdx = dvxdx;
    p->dvydx = dvydx;
    p->dvzdx = dvzdx;
    p->dvxdy = dvxdy;
    p->dvydy = dvydy;
    p->dvzdy = dvzdy;
    p->dvxdz = dvxdz;
    p->dvydz = dvydz;
    p->dvzdz = dvzdz;
#endif

	trace = dvxdx+dvydy+dvzdz; /* same sign as divv */

	p->fDivv_t = trace; /* no H, comoving */
        pdivv_old = p->divv;
	p->divv =  trace + 3*Hcorr; /* include H, physical */
	p->curlv[0] = (dvzdy - dvydz); /* same in all coordinates */
	p->curlv[1] = (dvxdz - dvzdx);
	p->curlv[2] = (dvydx - dvxdy);
/* Prior: ALPHAMUL 10 on top -- make pre-factor for c instead then switch is limited to 1 or less */
	    {
	      gnorm = (grx*grx+gry*gry+grz*grz);
#ifdef CULLENDEHNEN
	      p->grx = grx/sqrt(gnorm);
	      p->gry = gry/sqrt(gnorm);
	      p->grz = grz/sqrt(gnorm);
#endif
	    if (gnorm > 0) gnorm=1/sqrt(gnorm);
	    grx *= gnorm;
	    gry *= gnorm;
	    grz *= gnorm;
	    dvdr = (((dvxdx+Hcorr)*grx+dvxdy*gry+dvxdz*grz)*grx
		+  (dvydx*grx+(dvydy+Hcorr)*gry+dvydz*grz)*gry
		+  (dvzdx*grx+dvzdy*gry+(dvzdz+Hcorr)*grz)*grz);


	    dvds = (p->divv < 0 ? 1.5*(dvdr -(1./3.)*p->divv) : dvdr );
	    }

#if defined(DIFFUSION) || defined(CULLENDEHNEN2)
        {
        double onethirdtrace = (1./3.)*trace;
        /* Build Traceless Strain Tensor (already normalized) */
        double sxx = dvxdx - onethirdtrace; /* pure compression/expansion doesn't diffuse */
        double syy = dvydy - onethirdtrace;
        double szz = dvzdz - onethirdtrace;
        double sxy = 0.5*(dvxdy + dvydx); /* pure rotation doesn't diffuse */
        double sxz = 0.5*(dvxdz + dvzdx);
        double syz = 0.5*(dvydz + dvzdy);
        double S2 = (sxx*sxx + syy*syy + szz*szz + 2*(sxy*sxy + sxz*sxz + syz*syz));
        // S2 Frobenius Norm^2 = trace(S ST) = Simpler as just sum_ij Sij^2 (used here)  Note: Sxy = Syx
/*	printf(" %g %g   %g %g %g  %g\n",p->fDensity,p->divv,p->curlv[0],p->curlv[1],p->curlv[2],fNorm1*sqrt(2*(sxx*sxx + syy*syy + szz*szz + 2*(sxy*sxy + sxz*sxz + syz*syz))) );*/
        }
#endif

#ifdef DENSITYU
    assert(p->fDensityU > 0);
#endif

	}

#else /* FITDVDX */

#include "stiff.h"
/* Gather only version -- don't use */
void DenDVDX(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	FLOAT ih2,ih,r2,rs,rs1,fDensity,fNorm,fNorm1,vFac,S;
	FLOAT dvxdx , dvxdy , dvxdz, dvydx , dvydy , dvydz, dvzdx , dvzdy , dvzdz;
	FLOAT Wdvxdx , Wdvxdy , Wdvxdz, Wdvydx , Wdvydy , Wdvydz, Wdvzdx , Wdvzdy , Wdvzdz;
	FLOAT Wdxdx, Wdxdy, Wdydy, Wdydz, Wdzdz, Wdxdz;
	FLOAT Wdvx,Wdvy,Wdvz,Wdx,Wdy,Wdz,W;
	FLOAT dvx,dvy,dvz,dx,dy,dz,trace;
	FLOAT grx,gry,grz,gnorm,dvds,dvdr,c;
	double **a, *x, d;
	int *indx;

	PARTICLE *q;
	int i,j;
	unsigned int qiActive;

	ih2 = 4.0/BALL2(p);
	ih = sqrt(ih2);
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */
#ifdef SPH1D
	fNorm = (2./3.)*ih;
#else
	fNorm = M_1_PI*ih2*ih;
#endif
	fNorm1 = fNorm*ih2;
	fDensity = 0.0;
	Wdvxdx = 0; Wdvxdy = 0; Wdvxdz= 0;
	Wdvydx = 0; Wdvydy = 0; Wdvydz= 0;
	Wdvzdx = 0; Wdvzdy = 0; Wdvzdz= 0;
	Wdxdx = 0; Wdxdy = 0; Wdxdz= 0;
	Wdydy = 0; Wdydz= 0;  Wdzdz= 0;
	Wdvx = 0; Wdvy = 0; Wdvz = 0;
	Wdx = 0; Wdy = 0; Wdz = 0;
	W = 0;

	grx = 0;
	gry = 0;
	grz = 0;

	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		q = nnList[i].pPart;
		if (TYPETest(p,TYPE_ACTIVE)) TYPESet(q,TYPE_NbrOfACTIVE); /* important for SPH */
		qiActive |= q->iActive;
		KERNEL(rs,r2);
		fDensity += rs*q->fMass;
		DKERNEL(rs1,r2);
		rs1 *= q->fMass;
		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = (-p->vPred[0] + q->vPred[0])*vFac - dx*smf->H; /* NB: dx = px - qx */
		dvy = (-p->vPred[1] + q->vPred[1])*vFac - dy*smf->H;
		dvz = (-p->vPred[2] + q->vPred[2])*vFac - dz*smf->H;
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
	    if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
			dvy += -1.5*smf->PP.dOrbFreq*smf->PP.dWidth*vFac;
	    else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
	        dvy -= -1.5*smf->PP.dOrbFreq*smf->PP.dWidth*vFac;
#endif

		grx += (q->uPred-p->uPred)*dx*rs1; /* Grad P estimate */
		gry += (q->uPred-p->uPred)*dy*rs1;
		grz += (q->uPred-p->uPred)*dz*rs1;

		rs = 1;
		Wdvxdx += dvx*dx*rs;
		Wdvxdy += dvx*dy*rs;
		Wdvxdz += dvx*dz*rs;
		Wdvydx += dvy*dx*rs;
		Wdvydy += dvy*dy*rs;
		Wdvydz += dvy*dz*rs;
		Wdvzdx += dvz*dx*rs;
		Wdvzdy += dvz*dy*rs;
		Wdvzdz += dvz*dz*rs;
		Wdxdx += dx*dx*rs;
		Wdxdy += dx*dy*rs;
		Wdxdz += dx*dz*rs;
		Wdydy += dy*dy*rs;
		Wdydz += dy*dz*rs;
		Wdzdz += dz*dz*rs;
		Wdvx += dvx*rs;
		Wdvy += dvy*rs;
		Wdvz += dvz*rs;
		Wdx += dx*rs;
		Wdy += dy*rs;
		Wdz += dz*rs;
		W += rs;
		}
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);

/*	printf("TEST %d  %g %g  %g %g %g",p->iOrder,p->fDensity,p->divv,p->curlv[0],p->curlv[1],p->curlv[2]);*/
	p->fDensity = fNorm*fDensity;
	fNorm1 /= p->fDensity;

/* new stuff */
	a = matrix(1,4,1,4);

	for (i=1;i<=4;i++) for (j=1;j<=4;j++) a[i][j] = 0;
	a[1][1] = W;
	a[1][2] = Wdx;
	a[1][3] = Wdy;
	a[1][4] = Wdy;
	a[2][1] = Wdx;
	a[2][2] = Wdxdx;
	a[2][3] = Wdxdy;
	a[2][4] = Wdxdz;
	a[3][1] = Wdy;
	a[3][2] = Wdxdy;
	a[3][3] = Wdydy;
	a[3][4] = Wdydz;
	a[4][1] = Wdz;
	a[4][2] = Wdxdz;
	a[4][3] = Wdydz;
	a[4][4] = Wdzdz;

	indx = ivector(1,4);
	ludcmp(a,4,indx,&d);

	x = vector(1,4);
	x[1] = Wdvx;
	x[2] = Wdvxdx;
	x[3] = Wdvxdy;
	x[4] = Wdvxdz;
	lubksb(a,4,indx,x);
	dvxdx =  -x[2];
	dvxdy =  -x[3];
	dvxdz =  -x[4];
	x[1] = Wdvy;
	x[2] = Wdvydx;
	x[3] = Wdvydy;
	x[4] = Wdvydz;
	lubksb(a,4,indx,x);
	dvydx =  -x[2];
	dvydy =  -x[3];
	dvydz =  -x[4];
	x[1] = Wdvz;
	x[2] = Wdvzdx;
	x[3] = Wdvzdy;
	x[4] = Wdvzdz;
	lubksb(a,4,indx,x);
	dvzdx =  -x[2];
	dvzdy =  -x[3];
	dvzdz =  -x[4];

/*
	if ((p->iOrder%400)==0) {
	    printf("%8d %7.3f %7.3f %7.3f  %f %f %f   %f %f %f   %f %f %f\n",p->iOrder,p->r[0],p->r[1],p->r[2],dvxdx,dvxdy,dvxdz,dvydz,dvydy,dvydz,dvzdx,dvzdy,dvzdz);
	    printf("%8d                          %f %f %f   %f %f %f   %f %f %f\n",p->iOrder,-2*3.14156*sin(2*3.14156*p->r[0]),0.,0.,0.,2*3.14156*cos(2*3.14156*p->r[1]),0., -2*3.14156*sin(2*3.14156*p->r[0]), 0., 0. );
	    }
*/

	free_ivector(indx,1,4);
	free_vector(x,1,4);
	free_matrix(a,1,4,1,4);

/* divv will be used to update density estimates in future */
	trace = dvxdx+dvydy+dvzdz;
	p->divv =  trace; /* physical */
/* Technically we could now dispense with curlv but it is used in several places as
   a work variable so we will retain it */
	p->curlv[0] = (dvzdy - dvydz);
	p->curlv[1] = (dvxdz - dvzdx);
	p->curlv[2] = (dvydx - dvxdy);

/* Prior: ALPHAMUL 10 on top -- make pre-factor for c instead then switch is limited to 1 or less */
#ifndef ALPHACMUL
#define ALPHACMUL 0.1
#endif
#ifndef DODVDS
	if (smf->iViscosityLimiter==2)
#endif
	    {
	    gnorm = (grx*grx+gry*gry+grz*grz);
	    if (gnorm > 0) gnorm=1/sqrt(gnorm);
	    grx *= gnorm;
	    gry *= gnorm;
	    grz *= gnorm;
	    dvdr = ((dvxdx*grx+dvxdy*gry+dvxdz*grz)*grx
		    +  (dvydx*grx+dvydy*gry+dvydz*grz)*gry
		    +  (dvzdx*grx+dvzdy*gry+dvzdz*grz)*grz);
#ifdef DODVDS
	    p->dvds =
#endif
            dvds = (p->divv < 0 ? 1.5*(dvdr -(1./3.)*p->divv) : dvdr );
	    }

	switch(smf->iViscosityLimiter) {
	case VISCOSITYLIMITER_NONE:
	    p->BalsaraSwitch=1;
	    break;
	case VISCOSITYLIMITER_BALSARA:
	    if (p->divv!=0.0) {
		p->BalsaraSwitch = fabs(p->divv)/
		    (fabs(p->divv)+sqrt(p->curlv[0]*p->curlv[0]+
					p->curlv[1]*p->curlv[1]+
					p->curlv[2]*p->curlv[2]));
		}
	    else {
		p->BalsaraSwitch = 0;
		}
	    break;
	case VISCOSITYLIMITER_JW:
	    c = sqrt(smf->gamma*p->uPred*(smf->gamma-1));
	    if (dvdr < 0 && dvds < 0 ) {
		p->BalsaraSwitch = -dvds/
		    (-dvds+sqrt(p->curlv[0]*p->curlv[0]+
				p->curlv[1]*p->curlv[1]+
				p->curlv[2]*p->curlv[2])+ALPHACMUL*c*ih)+ALPHAMIN;
		if (p->BalsaraSwitch > 1) p->BalsaraSwitch = 1;
		}
	    else {
		p->BalsaraSwitch = 0.01;
		}
	    break;
	    }
#ifdef DIFFUSION
        {
	double onethirdtrace = (1./3.)*trace;
	/* Build Traceless Strain Tensor (not yet normalized) */
	double sxx = dvxdx - onethirdtrace; /* pure compression/expansion doesn't diffuse */
	double syy = dvydy - onethirdtrace;
	double szz = dvzdz - onethirdtrace;
	double sxy = 0.5*(dvxdy + dvydx); /* pure rotation doesn't diffuse */
	double sxz = 0.5*(dvxdz + dvzdx);
	double syz = 0.5*(dvydz + dvzdy);
	/* diff coeff., nu ~ C L^2 S (add C via dMetalDiffusionConstant, assume L ~ h) */
	if (smf->bConstantDiffusion) p->diff = 1;
	else p->diff = fNorm1*0.25*BALL2(p)*sqrt(2*(sxx*sxx + syy*syy + szz*szz + 2*(sxy*sxy + sxz*sxz + syz*syz)));
/*	printf(" %g %g   %g %g %g  %g\n",p->fDensity,p->divv,p->curlv[0],p->curlv[1],p->curlv[2],fNorm1*sqrt(2*(sxx*sxx + syy*syy + szz*szz + 2*(sxy*sxy + sxz*sxz + syz*syz))) );*/
	}
#endif

	}
#endif /* FITDVDX */

void initSmoothBSw(void *p)
{
	}

void combSmoothBSw(void *p1,void *p2)
{
	}



void SmoothBSw(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{

#ifdef DODVDS
    FLOAT ih2,ih,r2,rs,c;
	FLOAT curlv[3],divv,fNorm,dvds;
	int i;
	PARTICLE *q;
	ih2 = 4.0/BALL2(p);
	ih = sqrt(ih2);
	curlv[0] = 0; curlv[1] = 0; curlv[2] = 0;
	divv = 0; dvds = 0;
	fNorm = 0;

	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		q = nnList[i].pPart;
		KERNEL(rs,r2);
		rs *= q->fMass;
		curlv[0] += q->curlv[0]*rs;
		curlv[1] += q->curlv[1]*rs;
		curlv[2] += q->curlv[2]*rs;
		divv += q->divv*rs;
		dvds += q->dvds*rs;
		fNorm += rs;
		}

	switch(smf->iViscosityLimiter) {
	case VISCOSITYLIMITER_NONE:
	    p->BalsaraSwitch=1;
	    break;
	case VISCOSITYLIMITER_BALSARA:
	    if (p->divv!=0.0) {
		p->BalsaraSwitch = fabs(divv)/ (fabs(divv)+sqrt(curlv[0]*curlv[0]+
					curlv[1]*curlv[1]+
					curlv[2]*curlv[2]));
		}
	    else {
		p->BalsaraSwitch = 0;
		}
	    break;
	case VISCOSITYLIMITER_JW:
/* Prior: ALPHAMUL 10 on top -- make pre-factor for c instead then switch is limited to 1 or less */
#ifndef ALPHACMUL
#define ALPHACMUL 0.1
#endif
	    c = sqrt(smf->gamma*p->uPred*(smf->gamma-1));
	    if (divv < 0 && dvds < 0 ) {
		p->BalsaraSwitch = -dvds/
		    (-dvds+sqrt(curlv[0]*curlv[0]+
				curlv[1]*curlv[1]+
				curlv[2]*curlv[2])+fNorm*ALPHACMUL*c*ih)+ALPHAMIN;
		if (p->BalsaraSwitch > 1) p->BalsaraSwitch = 1;
		}
	    else {
		p->BalsaraSwitch = 0.01;
		}
	    break;
	    }
#else
	assert(0);
#endif /* DODVDS */

	}

void initShockTrack(void *p)
{
#ifdef SHOCKTRACK
	if (TYPEQueryACTIVE((PARTICLE *) p )) {
		((PARTICLE *)p)->divrhov = 0.0;
		}
#endif
	}

void combShockTrack(void *p1,void *p2)
{
#ifdef SHOCKTRACK
	if (TYPEQueryACTIVE((PARTICLE *) p1 )) {
	        if (((PARTICLE *)p2)->divrhov > ((PARTICLE *)p2)->divrhov)
		  ((PARTICLE *)p1)->divrhov = ((PARTICLE *)p2)->divrhov;
		}
#endif
	}

/* Gather only version -- untested */
void ShockTrack(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	assert(0);
	}

/* Output is physical divv and curlv -- thus a*h_co*divv is physical */
void ShockTrackSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef SHOCKTRACK
	PARTICLE *q;
	double Mach,dv;
	int i,j;

	if (TYPEQueryACTIVE( p )) {
		/* p active */
		for (i=0;i<nSmooth;++i) {
	                q = nnList[i].pPart;
			Mach = 0;
			for (j=0;j<3;j++) {
			  dv = p->vPred[j] - q->gradrho[j];
			  Mach += dv*dv;
			}
			Mach = sqrt(Mach)/p->c;
			if (Mach > p->divrhov) p->divrhov = Mach;
			if (TYPEQueryACTIVE(q)) {
			  Mach = 0;
			  for (j=0;j<3;j++) {
			    dv = q->vPred[j] - p->gradrho[j];
			    Mach += dv*dv;
			  }
			  Mach = sqrt(Mach)/q->c;
			  if (Mach > q->divrhov) q->divrhov = Mach;
		        }
	        }
		}
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
 	                q = nnList[i].pPart;
			if (!TYPEQueryACTIVE(q)) continue; /* neither active */
			Mach = 0;
			for (j=0;j<3;j++) {
			  dv = q->vPred[j] - p->gradrho[j];
			  Mach += dv*dv;
			}
			Mach = sqrt(Mach)/q->c;
			if (Mach > q->divrhov) q->divrhov = Mach;
	        }
		}
#endif
	}

/* Original Particle */
void initHKPressureTermsParticle(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->uDotPdV = 0.0;
		((PARTICLE *)p)->uDotAV = 0.0;
#ifdef DEBUG
		((PARTICLE *)p)->PdVvisc = 0.0;
		((PARTICLE *)p)->PdVpres = 0.0;
#endif
		}
	}

/* Cached copies of particle */
void initHKPressureTerms(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->uDotPdV = 0.0;
		((PARTICLE *)p)->uDotAV = 0.0;
#ifdef DEBUG
		((PARTICLE *)p)->PdVvisc = 0.0;
		((PARTICLE *)p)->PdVpres = 0.0;
#endif
		ACCEL(p,0) = 0.0;
		ACCEL(p,1) = 0.0;
		ACCEL(p,2) = 0.0;
		}
	}

void combHKPressureTerms(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->uDotPdV += ((PARTICLE *)p2)->uDotPdV;
		((PARTICLE *)p1)->uDotAV += ((PARTICLE *)p2)->uDotAV;
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL(p1,0) += ACCEL(p2,0);
		ACCEL(p1,1) += ACCEL(p2,1);
		ACCEL(p1,2) += ACCEL(p2,2);
		}
	}

/* Gather only version -- (untested)  */
void HKPressureTerms(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pQonRho2,qQonRho2,qhdivv;
	FLOAT ph,pc,pDensity,visc,absmu,qh,hav;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	if (!TYPEQueryACTIVE(p)) return;

	pc = p->c;
	pDensity = p->fDensity;
	pPoverRho2 = p->PoverRho2;
	ph = sqrt(0.25*BALL2(p));
	/* QonRho2 given same scaling with a as PonRho2 */
	pQonRho2 = (p->divv>0.0 ? 0.0 : fabs(p->divv)*ph*smf->a
				*(smf->alpha*pc + smf->beta*fabs(p->divv)*ph*smf->a)/pDensity);
	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (2./3.)/ph;
#else
	fNorm = M_1_PI*ih2/ph;
#endif
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);        /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= fNorm1 * q->fMass;;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
	    if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
			dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
	    else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
	        dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;

		if (dvdotdr>0.0) {
			p->uDotPdV += rs1*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
			rs1 *= (PRES_ACC(pPoverRho2,q->PoverRho2));
			rs1 *= aFac;
			ACCEL(p,0) -= rs1*dx;
			ACCEL(p,1) -= rs1*dy;
			ACCEL(p,2) -= rs1*dz;
			}
		else {
			qh=sqrt(0.25*BALL2(q));
			qhdivv = qh*fabs(q->divv)*smf->a; /* units of physical velocity */
			qQonRho2 = (qhdivv>0.0 ? 0.0 :
						qhdivv*(smf->alpha*q->c + smf->beta*qhdivv)/q->fDensity);
			visc = pQonRho2 + qQonRho2;
			/* mu -- same timestep criteria as standard sph above (for now) */
			hav=0.5*(qh+ph);
			absmu = -hav*dvdotdr*smf->a
				/(nnList[i].fDist2+0.01*hav*hav);
			if (absmu>p->mumax) p->mumax=absmu;
			p->uDotPdV += rs1 * (PRES_PDV(pPoverRho2,q->PoverRho2)) * dvdotdr;
			p->uDotAV += rs1 * (0.5*visc) * dvdotdr;
			rs1 *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
			rs1 *= aFac; /* convert to comoving acceleration */
			ACCEL(p,0) -= rs1*dx;
			ACCEL(p,1) -= rs1*dy;
			ACCEL(p,2) -= rs1*dz;
			}
		}
	}

/* Bulk viscosity and standard pressure forces */
void HKPressureTermsSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pQonRho2,qQonRho2,qhdivv;
	FLOAT ph,pc,pDensity,visc,absmu,qh,pMass,hav;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	pPoverRho2 = p->PoverRho2;
	ph = sqrt(0.25*BALL2(p));
	/* QonRho2 given same scaling with a as PonRho2 */
	pQonRho2 = (p->divv>0.0 ? 0.0 : fabs(p->divv)*ph*smf->a
				*(smf->alpha*pc + smf->beta*fabs(p->divv)*ph*smf->a)/pDensity );
	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (1./3.)/ph;
#else
	fNorm = 0.5*M_1_PI*ih2/ph;
#endif
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);        /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= fNorm1;
		rq = rs1*q->fMass;
		rp = rs1*pMass;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
	    if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
			dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
	    else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
	        dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;

		if (dvdotdr>0.0) {
			if (TYPEQueryACTIVE(p)) {
                p->uDotPdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
				rq *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rq *= aFac;
				ACCEL(p,0) -= rq*dx;
				ACCEL(p,1) -= rq*dy;
				ACCEL(p,2) -= rq*dz;
		        }
			if (TYPEQueryACTIVE(q)) {
				q->uDotPdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
				rp *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rp *= aFac; /* convert to comoving acceleration */
				ACCEL(q,0) += rp*dx;
				ACCEL(q,1) += rp*dy;
				ACCEL(q,2) += rp*dz;
				}
			}
		else {
			qh=sqrt(0.25*BALL2(q));
			qhdivv = qh*fabs(q->divv)*smf->a; /* units of physical velocity */
			qQonRho2 = (qhdivv>0.0 ? 0.0 :
						qhdivv*(smf->alpha*q->c + smf->beta*qhdivv)/q->fDensity);
			visc = pQonRho2 + qQonRho2;
			/* mu -- same timestep criteria as standard sph above (for now) */
			hav=0.5*(qh + ph);
			absmu = -hav*dvdotdr*smf->a/(nnList[i].fDist2+0.01*hav*hav);
			if (TYPEQueryACTIVE(p)) {
				if (absmu>p->mumax) p->mumax=absmu;
				p->uDotPdV += rq*(PRES_PDV(pPoverRho2,q->PoverRho2) + 0.5*visc)*dvdotdr;
				rq *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
				rq *= aFac; /* convert to comoving acceleration */
				ACCEL(p,0) -= rq*dx;
				ACCEL(p,1) -= rq*dy;
				ACCEL(p,2) -= rq*dz;
		        }
			if (TYPEQueryACTIVE(q)) {
				if (absmu>q->mumax) q->mumax=absmu;
				q->uDotPdV += rp*(PRES_PDV(q->PoverRho2,pPoverRho2) + 0.5*visc)*dvdotdr;
				rp *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
				rp *= aFac; /* convert to comoving acceleration */
				ACCEL(q,0) += rp*dx;
				ACCEL(q,1) += rp*dy;
				ACCEL(q,2) += rp*dz;
				}
			}
		}
	}

/* Original Particle */
void initHKViscosityParticle(void *p)
{
	ACCEL(p,0) += ACCEL_PRES(p,0);
	ACCEL(p,1) += ACCEL_PRES(p,1);
	ACCEL(p,2) += ACCEL_PRES(p,2);
	}

/* Cached copies of particle */
void initHKViscosity(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->uDotAV = 0.0;
		ACCEL(p,0) = 0.0;
		ACCEL(p,1) = 0.0;
		ACCEL(p,2) = 0.0;
		}
	}

void combHKViscosity(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->uDotAV += ((PARTICLE *)p2)->uDotAV;
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL(p1,0) += ACCEL(p2,0);
		ACCEL(p1,1) += ACCEL(p2,1);
		ACCEL(p1,2) += ACCEL(p2,2);
		}
	}

/* Gather only version */
void HKViscosity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
        assert(0);
	}

/* Bulk viscosity */
void HKViscositySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pQonRho2,qQonRho2,qhdivv;
	FLOAT ph,pc,pDensity,visc,absmu,qh,pMass,hav;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	ph = sqrt(0.25*BALL2(p));
	/* QonRho2 given same scaling with a as PonRho2 */
	pQonRho2 = (p->divv>0.0 ? 0.0 : fabs(p->divv)*ph*smf->a
				*(smf->alpha*pc + smf->beta*fabs(p->divv)*ph*smf->a)/pDensity );
	ih2 = 4.0/BALL2(p);
#ifdef SPH1D
	fNorm = (1./3.)/ph;
#else
	fNorm = 0.5*M_1_PI*ih2/ph;
#endif
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);        /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= fNorm1;
		rq = rs1*q->fMass;
		rp = rs1*pMass;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
	    if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
			dvy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
	    else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
	        dvy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;

		if (dvdotdr<0.0) {
			qh=sqrt(0.25*BALL2(q));
			qhdivv = qh*fabs(q->divv)*smf->a; /* units of physical velocity */
			qQonRho2 = (qhdivv>0.0 ? 0.0 :
						qhdivv*(smf->alpha*q->c + smf->beta*qhdivv)/q->fDensity);

			visc = SWITCHCOMBINE(p,q)*(pQonRho2 + qQonRho2);
			/* mu -- same timestep criteria as standard sph above (for now) */
			hav=0.5*(qh + ph);
			absmu = -hav*dvdotdr*smf->a/(nnList[i].fDist2+0.01*hav*hav);
			if (TYPEQueryACTIVE(p)) {
				if (absmu>p->mumax) p->mumax=absmu;
		        p->uDotPdV += rq*0.5*visc*dvdotdr;
				rq *= visc;
				rq *= aFac; /* convert to comoving acceleration */
		        ACCEL(p,0) -= rq*dx;
		        ACCEL(p,1) -= rq*dy;
		        ACCEL(p,2) -= rq*dz;
		        }
			if (TYPEQueryACTIVE(q)) {
				if (absmu>q->mumax) q->mumax=absmu;
				q->uDotPdV += rp*0.5*visc*dvdotdr;
				rp *= visc;
				rp *= aFac; /* convert to comoving acceleration */
		        ACCEL(q,0) += rp*dx;
		        ACCEL(q,1) += rp*dy;
		        ACCEL(q,2) += rp*dz;
				}
			}
		}
	}


int CompISORT(const void * a, const void * b) {
    return ( (((ISORT *) a)->r2 < ((ISORT *) b)->r2) ? -1 : 1 );
    }

#ifdef STARFORM
void initDistDeletedGas(void *p1)
{
	if(!TYPETest(((PARTICLE *)p1), TYPE_DELETED)) {
    /*
     * Zero out accumulated quantities.
     */
		((PARTICLE *)p1)->fMass = 0;
		((PARTICLE *)p1)->v[0] = 0;
		((PARTICLE *)p1)->v[1] = 0;
		((PARTICLE *)p1)->v[2] = 0;
		((PARTICLE *)p1)->u = 0;
		((PARTICLE *)p1)->uPred = 0;
		((PARTICLE *)p1)->uDot = 0.0;
		((PARTICLE *)p1)->fMetals = 0.0;
		((PARTICLE *)p1)->fMFracOxygen = 0.0;
		((PARTICLE *)p1)->fMFracIron = 0.0;
#ifdef TWOPHASE
		((PARTICLE *)p1)->fMassHot = 0;
		((PARTICLE *)p1)->uHot = 0;
		((PARTICLE *)p1)->uHotPred = 0;
#endif
		}
    }

void combDistDeletedGas(void *vp1,void *vp2)
{
    /*
     * Distribute u, v, and fMetals for particles returning from cache
     * so that everything is conserved nicely.
     */
	PARTICLE *p1 = vp1;
	PARTICLE *p2 = vp2;

	if(!TYPETest((p1), TYPE_DELETED)) {
		FLOAT delta_m = p2->fMass;
		FLOAT m_new,f1,f2;
		FLOAT fTCool; /* time to cool to zero */

		m_new = p1->fMass + delta_m;
		if (delta_m > 0) {
			f1 = p1->fMass /m_new;
			f2 = delta_m  /m_new;
			if(p1->uDot < 0.0) /* margin of 1% to avoid roundoff
					    * problems */
				fTCool = 1.01*p1->uPred/p1->uDot;
			p1->v[0] = f1*p1->v[0] + f2*p2->v[0];
			p1->v[1] = f1*p1->v[1] + f2*p2->v[1];
			p1->v[2] = f1*p1->v[2] + f2*p2->v[2];
			p1->fMetals = f1*p1->fMetals + f2*p2->fMetals;
			p1->fMFracOxygen = f1*p1->fMFracOxygen + f2*p2->fMFracOxygen;
			p1->fMFracIron = f1*p1->fMFracIron + f2*p2->fMFracIron;
#ifdef COOLDEBUG
			assert(p1->u >= 0.0);
#endif
			if(p1->uDot < 0.0)
				p1->uDot = p1->uPred/fTCool;

#ifdef TWOPHASE
		    FLOAT mHot_new = p1->fMassHot + p2->fMassHot;
            if (mHot_new > 0) {
                FLOAT f1_hot = p1->fMassHot/mHot_new;
                FLOAT f2_hot = p2->fMassHot/mHot_new;
                FLOAT mCold_new = m_new-mHot_new;
                if(!(mCold_new > 0)) printf("mCold_new: %e m1: %e m2: %e m1_hot: %e m2_hot: %e\n", mCold_new, p1->fMass, p2->fMass, p1->fMassHot, p2->fMassHot);
                assert(mCold_new > 0);
                FLOAT f1_cold = (p1->fMass-p1->fMassHot)/mCold_new;
                FLOAT f2_cold = (delta_m-p2->fMassHot)/mCold_new;
                p1->uHot =     f1_hot*p1->uHot     + f2_hot*p2->uHot;
                p1->uHotPred = f1_hot*p1->uHotPred + f2_hot*p2->uHotPred;
                p1->u =        f1_cold*p1->u       + f2_cold*p2->u;
                p1->uPred =    f1_cold*p1->uPred   + f2_cold*p2->uPred;

                p1->fMassHot = mHot_new;
                }
            else
#endif
                {
                p1->u = f1*p1->u + f2*p2->u;
                p1->uPred = f1*p1->uPred + f2*p2->uPred;
                }
            assert(p1->u > 0);
            assert(p1->uPred > 0);

			p1->fMass = m_new;
#ifdef TWOPHASE
            assert(p1->fMassHot < p1->fMass);
#endif
            }
		}
    }

void DistDeletedGas(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
    PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs,rstot,delta_m,m_new,f1,f2;
	FLOAT fTCool; /* time to cool to zero */
	int i;

	assert(TYPETest(p, TYPE_GAS));

    if (p->fMass == 0.0) return;   /* the particle to be deleted has NOTHING */

	ih2 = 4.0/BALL2(p);
        rstot = 0;
	for (i=0;i<nSmooth;++i) {
            q = nnList[i].pPart;
	    if(TYPETest(q, TYPE_DELETED)) continue;
	    assert(TYPETest(q, TYPE_GAS));
            r2 = nnList[i].fDist2*ih2;
            KERNEL(rs,r2);
            rstot += rs;
        }
	if(rstot <= 0.0) {
	    /* we have a particle to delete and nowhere to put its mass
	     * => we will keep it around */
	    pkdUnDeleteParticle(smf->pkd, p);
	    return;
	    }
	assert(rstot > 0.0);
	fNorm = 1./rstot;
	assert(p->fMass >= 0.0);
#ifdef TWOPHASE
    if (p->fMassHot > 0)
    {
        ISORT *isort;
        isort = (ISORT *) malloc(sizeof(ISORT)*nSmooth);
        for (i=0;i<nSmooth;++i) {
            isort[i].pNN = &(nnList[i]);
            isort[i].r2=-1.*nnList[i].pPart->fMassHot;
            }
        qsort( isort, nSmooth, sizeof(ISORT), CompISORT );
        for (i=0;i<nSmooth;++i) {
            q = isort[i].pNN->pPart;
            m_new = q->fMass + p->fMass;
            double mHot_new = q->fMassHot + p->fMassHot;
            /* Cached copies can have zero mass:skip them */
            if (m_new == 0) continue;
            f1 = q->fMass/m_new;
            f2 = p->fMass/m_new;
            double mCold_new = m_new-mHot_new;
            if(!(mCold_new > 0)) printf("mCold_new: %e m1: %e m2: %e m1_hot: %e m2_hot: %e\n", mCold_new, p->fMass, q->fMass, p->fMassHot, q->fMassHot);
            assert(mCold_new > 0);
            double f1_cold = (q->fMass-q->fMassHot)/mCold_new;
            double f2_cold = (p->fMass-p->fMassHot)/mCold_new;
            q->fMass = m_new;
            q->fMassHot = mHot_new;
            assert(p->fMassHot < p->fMass);
            if(q->uDot < 0.0) /* margin of 1% to avoid roundoff error */
            fTCool = 1.01*q->uPred/q->uDot;

                /* Only distribute the properties
                 * to the other particles on the "home" machine.
                 * u, v, and fMetals will be distributed to particles
                 * that come through the cache in the comb function.
                 */
            q->u = f1_cold*q->u+f2_cold*p->u;
            q->uPred = f1_cold*q->uPred+f2_cold*p->uPred;
            if (mHot_new > 0) {
                double f1_hot = q->fMassHot/mHot_new;
                double f2_hot = p->fMassHot/mHot_new;
                q->uHot = f1_hot*q->uHot+f2_hot*p->uHot;
                q->uHotPred = f1_hot*q->uHotPred+f2_hot*p->uHotPred;
                }

#ifdef COOLDEBUG
            assert(q->u >= 0.0);
#endif
#ifdef SLIDING_PATCH
#error Cannot do dist deleted gas with SLIDING_PATCH
#endif
            q->v[0] = f1*q->v[0]+f2*p->v[0];
            q->v[1] = f1*q->v[1]+f2*p->v[1];
            q->v[2] = f1*q->v[2]+f2*p->v[2];
            q->fMetals = f1*q->fMetals + f2*p->fMetals;
            q->fMFracOxygen = f1*q->fMFracOxygen + f2*p->fMFracOxygen;
            q->fMFracIron = f1*q->fMFracIron + f2*p->fMFracIron;
            if(q->uDot < 0.0) /* make sure we don't shorten cooling time */
                q->uDot = q->uPred/fTCool;
            return;
        }
    }
#endif
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		if(TYPETest(q, TYPE_DELETED)) continue;

		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
	    /*
	     * All these quantities are per unit mass.
	     * Exact if only one gas particle being distributed or in serial
	     * Approximate in parallel (small error).
	     */
		delta_m = rs*fNorm*p->fMass;
		m_new = q->fMass + delta_m;
		/* Cached copies can have zero mass: skip them */
		if (m_new == 0) continue;
		f1 = q->fMass /m_new;
		f2 = delta_m  /m_new;
		q->fMass = m_new;
		if(q->uDot < 0.0) /* margin of 1% to avoid roundoff error */
			fTCool = 1.01*q->uPred/q->uDot;

                /* Only distribute the properties
                 * to the other particles on the "home" machine.
                 * u, v, and fMetals will be distributed to particles
                 * that come through the cache in the comb function.
                 */
		q->u = f1*q->u+f2*p->u;
		q->uPred = f1*q->uPred+f2*p->uPred;
#if defined(UNONCOOL) && !defined(TWOPHASE)
		q->uHot = f1*q->uHot+f2*p->uHot;
		q->uHotPred = f1*q->uHotPred+f2*p->uHotPred;
#endif
		assert(q->u > 0.0);
		assert(q->uPred > 0.0);
#ifdef SLIDING_PATCH
#error Cannot do dist deleted gas with SLIDING_PATCH
#endif
		q->v[0] = f1*q->v[0]+f2*p->v[0];
		q->v[1] = f1*q->v[1]+f2*p->v[1];
		q->v[2] = f1*q->v[2]+f2*p->v[2];
		q->fMetals = f1*q->fMetals + f2*p->fMetals;
                q->fMFracOxygen = f1*q->fMFracOxygen + f2*p->fMFracOxygen;
                q->fMFracIron = f1*q->fMFracIron + f2*p->fMFracIron;
		if(q->uDot < 0.0) /* make sure we don't shorten cooling time */
			q->uDot = q->uPred/fTCool;
        }
}

#define PROMOTE_SUMWEIGHT(p_) (((PARTICLE *) (p_))->curlv[0])
#define PROMOTE_SUMUPREDWEIGHT(p_) (((PARTICLE *) (p_))->curlv[1])
#define PROMOTE_UPREDINIT(p_) (((PARTICLE *) (p_))->curlv[2])

void initPromoteToHotGas(void *p1)
    {
    TYPEReset(((PARTICLE *) p1),TYPE_PROMOTED);
    PROMOTE_SUMWEIGHT(p1) = 0; /* store weight total */
    PROMOTE_SUMUPREDWEIGHT(p1) = 0; /* store u x weight total */
    PROMOTE_UPREDINIT(p1) = ((PARTICLE *) p1)->uPred; /* store uPred */
    }

void combPromoteToHotGas(void *p1,void *p2)
    {
    if(TYPETest(((PARTICLE *) p2), TYPE_PROMOTED)) {
        TYPESet(((PARTICLE *) p1),TYPE_PROMOTED);
        if (((PARTICLE *) p2)->fTimeCoolIsOffUntil > ((PARTICLE *) p1)->fTimeCoolIsOffUntil) ((PARTICLE *) p1)->fTimeCoolIsOffUntil = ((PARTICLE *) p2)->fTimeCoolIsOffUntil;
        }
    PROMOTE_SUMWEIGHT(p1) += PROMOTE_SUMWEIGHT(p2);
    PROMOTE_SUMUPREDWEIGHT(p1) += PROMOTE_SUMUPREDWEIGHT(p2);
    }

void PromoteToHotGas(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
#ifdef NOCOOLING
    return;
#else
    PARTICLE *q;
    FLOAT fFactor,ph,ih2,r2,rs,rstot;
    FLOAT Tp,Tq,up52,Prob,mPromoted;
    double xc,yc,zc,dotcut2,dot;
	int i,nCold,nHot;

	assert(TYPETest(p, TYPE_GAS));
	assert(TYPETest(p, TYPE_FEEDBACK));
	assert(!TYPETest(p, TYPE_PROMOTED));
	assert(!TYPETest(p, TYPE_DELETED));
    ph = sqrt(BALL2(p)*0.25);
    ih2 = 1/(ph*ph);
    /* Exclude cool particles */
    Tp = CoolCodeEnergyToTemperature( smf->pkd->Cool, &p->CoolParticle, p->uPred, p->fDensity, p->fMetals );
    if (Tp <= smf->dEvapMinTemp) return;

    up52 = pow(p->uPred,2.5);
    rstot = 0;
    xc = 0; yc = 0; zc = 0;
    nCold = 0;
    assert(nSmooth > 0);
	for (i=0;i<nSmooth;++i) {
        q = nnList[i].pPart;
        if (p->iOrder == q->iOrder) continue;
	    if (TYPETest(q, TYPE_DELETED) || (TYPETest(q, TYPE_FEEDBACK) && !TYPETest(q, TYPE_PROMOTED))) continue;
        Tq = CoolCodeEnergyToTemperature( smf->pkd->Cool, &q->CoolParticle, q->uPred, q->fDensity, q->fMetals );
        if (Tq >= smf->dEvapMinTemp) continue;  /* Exclude hot particles */
	    assert(TYPETest(q, TYPE_GAS));
        assert(!TYPETest(p, TYPE_STAR));
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
        rstot += rs;
  		xc += rs*nnList[i].dx;
		yc += rs*nnList[i].dy;
		zc += rs*nnList[i].dz;
        nCold++;
        }
    if (nSmooth <= nCold) {
        printf("PROMOTION ERROR: nSmooth %d nCold %d\n", nSmooth, nCold);
        assert(0);
    }

    if (rstot == 0) return;

    /* Check for non-edge hot particle  theta = 45 deg, cos^2 = 0.5 */
    dotcut2 = (xc*xc+yc*yc+zc*zc)*0.5;

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		if (p->iOrder == q->iOrder) continue;
		if (TYPETest(q, TYPE_DELETED)) continue;
		Tq = CoolCodeEnergyToTemperature( smf->pkd->Cool, &q->CoolParticle, q->uPred, q->fDensity, q->fMetals );
#ifdef TWOPHASE
		if (q->uHot == 0 && Tq <= smf->dEvapMinTemp) continue;
#else
		if (Tq <= smf->dEvapMinTemp) continue;
#endif /* TWOPHASE */
		dot = xc*nnList[i].dx + yc*nnList[i].dy + zc*nnList[i].dz;
		if (dot > 0 && dot*dot > dotcut2*nnList[i].fDist2) {
            //printf("promote (hot excluded): %d %d  %g %g  (%g %g %g) (%g %g %g)\n",p->iOrder,q->iOrder,Tp, Tq,xc,yc,zc,nnList[i].dx,nnList[i].dy,nnList[i].dz);
            return;
            }
        }

    /* Area = h^2 4 pi nCold/nSmooth */
	nHot=nSmooth-nCold;
	assert(nHot > 0);
    fFactor = smf->dDeltaStarForm*smf->dEvapCoeffCode*ph*12.5664*1.5/(nHot)/rstot;
	//printf("CHECKAREA2: %e %d %e %d %d %e %e %e\n", smf->dTime, p->iOrder, 12.5664*ph*ph*1.5/(nHot), nSmooth, nCold, xc, yc, zc);

    mPromoted = 0;
	for (i=0;i<nSmooth;++i) {
        q = nnList[i].pPart;
        if (p->iOrder == q->iOrder) continue;
	    if(TYPETest(q, TYPE_DELETED) || (TYPETest(q, TYPE_FEEDBACK) && !TYPETest(q, TYPE_PROMOTED))) continue;
        Tq = CoolCodeEnergyToTemperature( smf->pkd->Cool, &q->CoolParticle, q->uPred, q->fDensity, q->fMetals );
        if (Tq >= smf->dEvapMinTemp ) continue;  /* Exclude hot particles */
	    assert(TYPETest(q, TYPE_GAS));
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
        PROMOTE_SUMWEIGHT(q) += p->fMass;
        PROMOTE_SUMUPREDWEIGHT(q) += p->fMass*p->uPred;

        /* cf. Weaver etal'77 mdot = 4.13d-14 * (dx^2/4 !pi) (Thot^2.5-Tcold^2.5)/dx - 2 udot mHot/(k T/mu)
           Kernel sets total probability to 1 */
        Prob = fFactor*(up52-pow(q->uPred,2.5))*rs/q->fMass;
		//printf("promote?: %d %d %g %g %g  %g %g %g\n",p->iOrder,q->iOrder,Tp, Tq, ph, fFactor*(up52-pow(q->uPred,2.5))*rs, q->fMass, Prob);
        if ( (rand()/((double) RAND_MAX)) < Prob) {
            mPromoted += q->fMass;
            //printf("promote? MASS: %d %d %g %g %g  %g + %g %g\n",p->iOrder,q->iOrder,Tp, Tq, ph, fFactor*(up52-pow(q->uPred,2.5))*rs, q->fMass, Prob);
            }
        }

    if (mPromoted > 0) {
        double dTimeCool = smf->dTime + 0.9999*smf->dDeltaStarForm;
        ISORT *isort;

        isort = (ISORT *) malloc(sizeof(ISORT)*nSmooth);
        for (i=0;i<nSmooth;++i) {
            isort[i].pNN = &(nnList[i]);
            isort[i].r2=nnList[i].fDist2;
            }
        qsort( isort, nSmooth, sizeof(ISORT), CompISORT );

        for (i=0;i<nSmooth;++i) {
            q = isort[i].pNN->pPart;
            if (p->iOrder == q->iOrder) continue;
            if (TYPETest(q, TYPE_DELETED) || TYPETest(q, TYPE_FEEDBACK) || TYPETest(q, TYPE_PROMOTED)) continue;
            Tq = CoolCodeEnergyToTemperature( smf->pkd->Cool, &q->CoolParticle, q->uPred, q->fDensity, q->fMetals );
            if (Tq >= smf->dEvapMinTemp ) continue;  /* Exclude hot particles */
            assert(TYPETest(q, TYPE_GAS));

            if (dTimeCool > q->fTimeCoolIsOffUntil) q->fTimeCoolIsOffUntil = dTimeCool;
            TYPESet(q, TYPE_PROMOTED|TYPE_FEEDBACK);
            mPromoted -= q->fMass;
            //printf("promote? YES: %d %d %g %g %g  %g - %g %g\n",p->iOrder,q->iOrder,Tp, Tq, ph, fFactor*(up52-pow(q->uPred,2.5))*rs, q->fMass, Prob);
            if (mPromoted < q->fMass*0.1) break;
            }
        free(isort);
        }

#endif
}

void initShareWithHotGas(void *p1)
    {
	if(!TYPETest(((PARTICLE *)p1), TYPE_DELETED)) {
		((PARTICLE *)p1)->u = 0;
		((PARTICLE *)p1)->uPred = 0;
		}
    }

void combShareWithHotGas(void *vp1,void *vp2)
    {
	PARTICLE *p1 = vp1;
	PARTICLE *p2 = vp2;

	if(!TYPETest((p1), TYPE_DELETED)) {
        p1->u = p1->u + p2->u;
        p1->uPred = p1->uPred + p2->uPred;
		}
    }

void ShareWithHotGas(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
	PARTICLE *q;
	FLOAT rsmax,uavg,umin;
	FLOAT dE,Eadd,factor,Tp;
	int i,nPromoted;

	assert(TYPETest(p, TYPE_GAS));
	assert(TYPETest(p, TYPE_FEEDBACK));
	assert(!TYPETest(p, TYPE_PROMOTED));
    Tp = CoolCodeEnergyToTemperature( smf->pkd->Cool, &p->CoolParticle, p->uPred, p->fDensity, p->fMetals );
    if (Tp <= smf->dEvapMinTemp) return;

    rsmax = 1.0;
    nPromoted = 0;

    dE = 0;
    umin = FLT_MAX;
	for (i=0;i<nSmooth;++i) {
        q = nnList[i].pPart;
	    if (TYPETest(q, TYPE_PROMOTED)) {
            nPromoted++;
            uavg = (rsmax*q->fMass*PROMOTE_UPREDINIT(q) + PROMOTE_SUMUPREDWEIGHT(q))/
                (rsmax*q->fMass + PROMOTE_SUMWEIGHT(q));
            if (uavg < umin) umin=uavg;
            Eadd = (uavg-PROMOTE_UPREDINIT(q))*q->fMass;
            if (Eadd < 0) Eadd=0;
            dE += p->fMass/PROMOTE_SUMWEIGHT(q)*Eadd;
            }
        }

    if (!nPromoted || dE == 0 || p->uPred <= umin) return;
    factor = ((p->uPred-umin)*p->fMass)/dE;
    if (factor > 1) factor=1;

	for (i=0;i<nSmooth;++i) {
        q = nnList[i].pPart;
	    if (TYPETest(q, TYPE_PROMOTED)) {
            nPromoted++;
            uavg = (rsmax*q->fMass*PROMOTE_UPREDINIT(q) + PROMOTE_SUMUPREDWEIGHT(q))/
                (rsmax*q->fMass + PROMOTE_SUMWEIGHT(q));
            if (uavg < umin) umin=uavg;
            Eadd = (uavg-PROMOTE_UPREDINIT(q))*q->fMass;
            if (Eadd < 0) Eadd=0;
            dE = factor*p->fMass/PROMOTE_SUMWEIGHT(q)*Eadd;
            q->uPred += dE/q->fMass;
            q->u += dE/q->fMass;
            p->uPred -=  dE/p->fMass;
            p->u -=  dE/p->fMass;
            assert(q->uPred > 0);
            assert(q->u > 0);
            assert(p->uPred > 0);
            assert(p->u > 0);
            }
        }
    }

void DeleteGas(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
/* flag low mass gas particles for deletion */
	PARTICLE *q;
	FLOAT fMasstot,fMassavg;
	int i;

	assert(TYPETest(p, TYPE_GAS));
	fMasstot = 0;
#ifdef COOLDEBUG
	assert(p->fMass >= 0.0);
#endif

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
	    assert(TYPETest(q, TYPE_GAS));
		fMasstot += q->fMass;
        }

	fMassavg = fMasstot/(FLOAT) nSmooth;

	if (p->fMass < smf->dMinMassFrac*fMassavg) {
		pkdDeleteParticle(smf->pkd, p);
        }
	else {
		assert (p->fMass > 0.0);
		}

}

void initTreeParticleDistFBEnergy(void *p1)
{
    /* Convert energy and metals to non-specific quantities (not per mass)
     * to make it easier to divvy up SN energy and metals.
     */

    if(TYPETest((PARTICLE *)p1, TYPE_GAS)){
#ifdef TWOPHASE
		if(((PARTICLE *)p1)->fMassHot > 0) {
			((PARTICLE *)p1)->uDotFB *= ((PARTICLE *)p1)->fMassHot;
		}
		else {
			((PARTICLE *)p1)->uDotFB *= ((PARTICLE *)p1)->fMass;
		}
#else
        ((PARTICLE *)p1)->uDotFB *= ((PARTICLE *)p1)->fMass;
#endif /* TWOPHASE*/
        ((PARTICLE *)p1)->uDotESF *= ((PARTICLE *)p1)->fMass;
        ((PARTICLE *)p1)->fMetals *= ((PARTICLE *)p1)->fMass;
        ((PARTICLE *)p1)->fMFracOxygen *= ((PARTICLE *)p1)->fMass;
        ((PARTICLE *)p1)->fMFracIron *= ((PARTICLE *)p1)->fMass;
        }

    }

void initDistFBEnergy(void *p1)
    {
    /*
     * Warning: kludgery.  We need to accumulate mass in the cached
     * particle, but we also need to keep the original mass around.
     * Let's use the curlv field in the cached particle copy to hold the original
     * mass.  Note: original particle curlv's never modified.
     */
    ((PARTICLE *)p1)->curlv[0] = ((PARTICLE *)p1)->fMass;
#ifdef TWOPHASE
    ((PARTICLE *)p1)->curlv[1] = ((PARTICLE *)p1)->fMassHot;
#endif

    /*
     * Zero out accumulated quantities.
     */
    ((PARTICLE *)p1)->uDotFB = 0.0;
    ((PARTICLE *)p1)->uDotESF = 0.0;
    ((PARTICLE *)p1)->fMetals = 0.0;
    ((PARTICLE *)p1)->fMFracOxygen = 0.0;
    ((PARTICLE *)p1)->fMFracIron = 0.0;
    }

void combDistFBEnergy(void *p1,void *p2)
    {
    /*
     * See kludgery notice above.
     */
    FLOAT fAddedMass = ((PARTICLE *)p2)->fMass - ((PARTICLE *)p2)->curlv[0];
#ifdef TWOPHASE
    FLOAT fAddedMassHot = ((PARTICLE *)p2)->fMassHot - ((PARTICLE *)p2)->curlv[1];
    ((PARTICLE *)p1)->fMassHot += fAddedMassHot;
#endif

    ((PARTICLE *)p1)->fMass += fAddedMass;
#ifdef TWOPHASE
    assert(((PARTICLE *)p1)->fMassHot < ((PARTICLE *)p1)->fMass);
#endif
    ((PARTICLE *)p1)->uDotFB += ((PARTICLE *)p2)->uDotFB;
    ((PARTICLE *)p1)->uDotESF += ((PARTICLE *)p2)->uDotESF;
    ((PARTICLE *)p1)->fMetals += ((PARTICLE *)p2)->fMetals;
    ((PARTICLE *)p1)->fMFracOxygen += ((PARTICLE *)p2)->fMFracOxygen;
    ((PARTICLE *)p1)->fMFracIron += ((PARTICLE *)p2)->fMFracIron;
    ((PARTICLE *)p1)->fTimeCoolIsOffUntil = max( ((PARTICLE *)p1)->fTimeCoolIsOffUntil,
                ((PARTICLE *)p2)->fTimeCoolIsOffUntil );
    ((PARTICLE *)p1)->fTimeForm = max( ((PARTICLE *)p1)->fTimeForm,
                ((PARTICLE *)p2)->fTimeForm ); /* propagate FB time JMB 2/24/10 */
    }

#ifdef TOPHATFEEDBACK
#define FEEDBACKKERNEL(_rs,_r2)  {		\
	_rs = 1.0;				\
	}
#else
#define FEEDBACKKERNEL(_rs,_r2)  KERNEL(_rs,_r2)
#endif

void DistESF(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
    assert(TYPETest(p, TYPE_STAR));
	int i;
	PARTICLE *q;
	double fNorm, fNorm_u, ih2, r2, rs;
	double ESFRate = p->fMass*smf->dESFEnergy/smf->dESFTime;

	ih2 = 4.0/BALL2(p);
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
        FEEDBACKKERNEL(rs,r2);
		q = nnList[i].pPart;
		if(q->fMass > smf->dMaxGasMass) {
			continue;		/* Skip heavy particles */
		}
		fNorm_u += q->fMass*rs;
		rs *= fNorm;
		assert(TYPETest(q, TYPE_GAS));
	}
  assert(fNorm_u > 0.0);  /* be sure we have at least one neighbor */
  fNorm_u = 1./fNorm_u;
  for (i=0;i<nSmooth;++i) {
	FLOAT weight;
	q = nnList[i].pPart;
	if(q->fMass > smf->dMaxGasMass)
        continue;		/* Skip heavy particles */
	r2 = nnList[i].fDist2*ih2;
	FEEDBACKKERNEL(rs,r2);
	/* Remember: We are dealing with total energy rate and total metal
	 * mass, not energy/gram or metals per gram.
	 * q->fMass is in product to make units work for fNorm_u.
	 */
#ifdef VOLUMEFEEDBACK
	weight = rs*fNorm_u*q->fMass/q->fDensity;
#else
	weight = rs*fNorm_u*q->fMass;
#endif
	q->uDotESF += weight*ESFRate;
  }
}

void DistIonize(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
  PARTICLE *q;
  FLOAT mgot,mwanted,tCoolAgain;
  double T;
  int i,ngot;
  ISORT *isort;

  if ( p->fTimeForm < smf->dTime ) return; /* Is this my very first feedback? */
  assert(TYPETest(p, TYPE_STAR));

  isort = (ISORT *) malloc(sizeof(ISORT)*nSmooth);
  for (i=0;i<nSmooth;++i) {
      assert(TYPETest(nnList[i].pPart, TYPE_GAS));
      isort[i].pNN = &(nnList[i]);
      isort[i].r2=nnList[i].fDist2;
      }

  qsort( isort, nSmooth, sizeof(ISORT), CompISORT );

  mwanted = p->fMass*smf->dIonizeMultiple;
  mgot = 0;
  ngot = 0;

  for (i=0;i<nSmooth;++i) {
      q = isort[i].pNN->pPart;
	  T = CoolCodeEnergyToTemperature(smf->pkd->Cool,&q->CoolParticle, q->u, q->fDensity, q->fMetals/q->fMass);
      if (T > smf->dIonizeTMin) continue;
      /* Stop once we have enough cold mass ionized
	 -- currently only checks nSmooth neighbours max so may truncate if in a hot region
         Could use Stromgren sphere type calculation to select the amount of mass */
      if (mgot+q->fMass*(rand()/((double) RAND_MAX)) > mwanted) break;

      mgot += q->fMass;
      ngot++;
      tCoolAgain = smf->dIonizeTime+smf->dTime;
      if (tCoolAgain > q->fTimeCoolIsOffUntil) q->fTimeCoolIsOffUntil=tCoolAgain;
      if (T < smf->dIonizeT) {
          CoolInitEnergyAndParticleData( smf->pkd->Cool, &q->CoolParticle, &q->u, q->fDensity, smf->dIonizeT, q->fMetals );
          }
      }
//  printf("Ionize: Star %d: %g %g %d\n",p->iOrder,mgot,mwanted,ngot);
}

void DistFBMME(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
  PARTICLE *q;
  FLOAT fNorm,ih2,r2,rs,fNorm_u,fNorm_Pres,fAveDens;
  int i;

  if ( p->fMSN == 0.0 ){return;} /* Is there any feedback mass? */
  assert(TYPETest(p, TYPE_STAR));
  ih2 = 4.0/BALL2(p);
  fNorm_u = 0.0;
  fNorm_Pres = 0.0;
  fAveDens = 0.0;

  fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
  for (i=0;i<nSmooth;++i) {
    r2 = nnList[i].fDist2*ih2;
    FEEDBACKKERNEL(rs,r2);
    q = nnList[i].pPart;
    if(q->fMass > smf->dMaxGasMass) {
	continue;		/* Skip heavy particles */
	}
    fNorm_u += q->fMass*rs;
    rs *= fNorm;
    fAveDens += q->fMass*rs;
    fNorm_Pres += q->fMass*q->uPred*rs;
    assert(TYPETest(q, TYPE_GAS));
  }
  assert(fNorm_u > 0.0);  /* be sure we have at least one neighbor */
  fNorm_Pres *= (smf->gamma-1.0);

  assert(fNorm_u != 0.0);
  fNorm_u = 1./fNorm_u;
  for (i=0;i<nSmooth;++i) {
    FLOAT weight;
    q = nnList[i].pPart;
    if(q->fMass > smf->dMaxGasMass)
	  continue;		/* Skip heavy particles */
    r2 = nnList[i].fDist2*ih2;
    FEEDBACKKERNEL(rs,r2);
    /* Remember: We are dealing with total energy rate and total metal
     * mass, not energy/gram or metals per gram.
     * q->fMass is in product to make units work for fNorm_u.
     */
#ifdef VOLUMEFEEDBACK
    weight = rs*fNorm_u*q->fMass/q->fDensity;
#else
    weight = rs*fNorm_u*q->fMass;
#endif
    if (p->fNSN == 0.0) q->uDotFB += weight*p->uDotFB;  /* uDot is erg/s not erg/g/s here */
    q->fMetals += weight*p->fSNMetals;
    q->fMFracOxygen += weight*p->fMOxygenOut;
    q->fMFracIron += weight*p->fMIronOut;
    q->fMass += weight*p->fMSN;
	if(weight > 0) {
		TYPESet(q, TYPE_FEEDBACK);
	}
#ifdef TWOPHASE
	FLOAT Tq = CoolCodeEnergyToTemperature( smf->pkd->Cool, &q->CoolParticle, q->uPred, q->fDensity, q->fMetals );
	if(Tq < smf->dMultiPhaseMinTemp && weight > 0 && p->fNSN > 0.0) {
		double fMassHot = q->fMassHot + weight*p->fMSN;
		double deltaMassLoad = weight*p->fMSN*smf->dFBInitialMassLoad;
		if (fMassHot+deltaMassLoad >= q->fMass) {
			deltaMassLoad = q->fMass - fMassHot;
			fMassHot = q->fMass;
			}
		else {
			fMassHot += deltaMassLoad;
		}
		q->fMassHot = fMassHot;
		assert(q->fMassHot < q->fMass);
		assert(q->fMassHot >= 0);
	}
#endif /* TWOPHASE */
  }
}

void DistFBEnergy(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
  PARTICLE *q;
  FLOAT fNorm,ih2,r2,rs,fNorm_u,fNorm_Pres,fAveDens,f2h2;
  FLOAT fBlastRadius,fShutoffTime,fmind;
  double dAge,  aFac, dCosmoDenFac;
  int i,imind;

  if (smf->bIonize && smf->dTime-p->fTimeForm >= 0 && smf->dTime-p->fTimeForm < smf->dIonizeTime)
      DistIonize(p,nSmooth,nnList,smf);

  if (smf->bESF && smf->dTime-p->fTimeForm >= 0 && smf->dTime-p->fTimeForm < smf->dESFTime)
      DistESF(p, nSmooth, nnList, smf);

  if ( p->fMSN == 0.0 ){return;}

  /* "Simple" ejecta distribution (see function above) */
  DistFBMME(p,nSmooth,nnList,smf);

  if (p->fNSN == 0) return;
  if (p->fTimeForm < 0.0 && smf->bBHTurnOffCooling == 1) return; /* don't let BHs have Supernovae  JMB 4/28/09 */

  /* The following ONLY deals with SNII Energy distribution */
  assert(TYPETest(p, TYPE_STAR));
  ih2 = 4.0/BALL2(p);
    aFac = smf->a;
    dCosmoDenFac = aFac*aFac*aFac;
  f2h2=BALL2(p);
  fNorm_u = 0.0;
  fNorm_Pres = 0.0;
  fAveDens = 0.0;
  dAge = smf->dTime - p->fTimeForm;
  if (dAge == 0.0) return;

  fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
  for (i=0;i<nSmooth;++i) {
    r2 = nnList[i].fDist2*ih2;
    FEEDBACKKERNEL(rs,r2);
    q = nnList[i].pPart;
    fNorm_u += q->fMass*rs;
    rs *= fNorm;
    fAveDens += q->fMass*rs;
    fNorm_Pres += q->fMass*q->uPred*rs;
    assert(TYPETest(q, TYPE_GAS));
  }
  fNorm_Pres *= (smf->gamma-1.0)/dCosmoDenFac;
    fAveDens /= dCosmoDenFac;
  if (smf->sn.iNSNIIQuantum > 0) {
    /* McCray + Kafatos (1987) ApJ 317 190*/
    fBlastRadius = smf->sn.dRadPreFactor*pow(p->fNSN / fAveDens, 0.2) *
	    pow(dAge,0.6)/aFac; /* eq 3 */
    /* TOO LONG    fShutoffTime = smf->sn.dTimePreFactor*pow(p->fMetals, -1.5)*
       pow(p->fNSN,0.3) / pow(fAveDens,0.7);*/
  }
  else {
    /* from McKee and Ostriker (1977) ApJ 218 148 */
    fBlastRadius = smf->sn.dRadPreFactor*pow(p->fNSN,0.32)*
	    pow(fAveDens,-0.16)*pow(fNorm_Pres,-0.2)/aFac;
  }
  if (smf->bShortCoolShutoff){
    /* End of snowplow phase */
    fShutoffTime = smf->sn.dTimePreFactor*pow(p->fNSN,0.31)*
      pow(fAveDens,0.27)*pow(fNorm_Pres,-0.64);
  } else{        /* McKee + Ostriker 1977 t_{max} */
    fShutoffTime = smf->sn.dTimePreFactor*pow(p->fNSN,0.32)*
      pow(fAveDens,0.34)*pow(fNorm_Pres,-0.70);
  }
  /* Shut off cooling for 3 Myr for stellar wind */
  if (p->fNSN < smf->sn.iNSNIIQuantum)
    fShutoffTime= 3e6 * SECONDSPERYEAR / smf->dSecUnit;

  fmind = BALL2(p);
  imind = 0;
  if ( p->uDotFB > 0.0 ) {
    if(smf->bSmallSNSmooth) {
      /* Change smoothing radius to blast radius
       * so that we only distribute mass, metals, and energy
       * over that range.
       */
      f2h2 = fBlastRadius*fBlastRadius;
      ih2 = 4.0/f2h2;
    }

    fNorm_u = 0.0;

    for (i=0;i<nSmooth;++i) {
      if ( nnList[i].fDist2 < fmind ){imind = i; fmind = nnList[i].fDist2;}
      if ( nnList[i].fDist2 < f2h2 || !smf->bSmallSNSmooth) {
	r2 = nnList[i].fDist2*ih2;
	FEEDBACKKERNEL(rs,r2);
	q = nnList[i].pPart;
#ifdef VOLUMEFEEDBACK
	fNorm_u += q->fMass/q->fDensity*rs;
#else
	fNorm_u += q->fMass*rs;
#endif
	assert(TYPETest(q, TYPE_GAS));
      }
    }
  }

  /* If there's no gas particle within blast radius,
     give mass and energy to nearest gas particle. */
  if (fNorm_u ==0.0){
    r2 = nnList[imind].fDist2*ih2;
    FEEDBACKKERNEL(rs,r2);
    /*
     * N.B. This will be NEGATIVE, but that's OK since it will
     * cancel out down below.
     */
#ifdef VOLUMEFEEDBACK
    fNorm_u = nnList[imind].pPart->fMass/nnList[imind].pPart->fDensity*rs;
#else
    fNorm_u = nnList[imind].pPart->fMass*rs;
#endif
  }

  assert(fNorm_u != 0.0);
  fNorm_u = 1./fNorm_u;
  for (i=0;i<nSmooth;++i) {
    FLOAT weight;
    q = nnList[i].pPart;
    if (smf->bSmallSNSmooth) {
      if ( (nnList[i].fDist2 <= f2h2) || (i == imind) ) {
	if( smf->bSNTurnOffCooling &&
	    (fBlastRadius*fBlastRadius >= nnList[i].fDist2)) {
	  q->fTimeCoolIsOffUntil = max(q->fTimeCoolIsOffUntil,
				       smf->dTime + fShutoffTime);
	  q->fTimeForm = smf->dTime; /* store SN FB time here JMB 2/24/10 */
	}

	r2 = nnList[i].fDist2*ih2;
	FEEDBACKKERNEL(rs,r2);
	/* Remember: We are dealing with total energy rate and total metal
	 * mass, not energy/gram or metals per gram.
	 * q->fMass is in product to make units work for fNorm_u.
	 */
#ifdef VOLUMEFEEDBACK
	weight = rs*fNorm_u*q->fMass/q->fDensity;
#else
	weight = rs*fNorm_u*q->fMass;
#endif
	q->uDotFB += weight*p->uDotFB;
      }
    } else {
      r2 = nnList[i].fDist2*ih2;
      FEEDBACKKERNEL(rs,r2);
      /* Remember: We are dealing with total energy rate and total metal
       * mass, not energy/gram or metals per gram.
       * q->fMass is in product to make units work for fNorm_u.
       */
#ifdef VOLUMEFEEDBACK
      weight = rs*fNorm_u*q->fMass/q->fDensity;
#else
      weight = rs*fNorm_u*q->fMass;
#endif
      q->uDotFB += weight*p->uDotFB;
      /*		printf("SNTEST: %d %g %g %g %g\n",q->iOrder,weight,sqrt(q->r[0]*q->r[0]+q->r[1]*q->r[1]+q->r[2]*q->r[2]),q->uDotFB,q->fDensity);*/

      if ( p->uDotFB > 0.0 && smf->bSNTurnOffCooling &&
	   (fBlastRadius*fBlastRadius >= nnList[i].fDist2)){
	q->fTimeCoolIsOffUntil = max(q->fTimeCoolIsOffUntil,
				     smf->dTime + fShutoffTime);
	q->fTimeForm = smf->dTime;  /* store SN FB time here JMB 2/24/10 */
      }
      /*	update mass after everything else so that distribution
		is based entirely upon initial mass of gas particle */
    }
  }
}

void postDistFBEnergy(PARTICLE *p1, SMF *smf)
{
    /* Convert energy and metals back to specific quantities (per mass)
       because we are done with our conservative calculations */

    if(TYPETest(p1, TYPE_GAS)){
#ifdef TWOPHASE
		if (p1->fMassHot > 0) {
			p1->uDotFB /= p1->fMassHot;
            if (!TYPETest(p1,TYPE_TWOPHASE)) {
                double E = 0;
                CoolInitEnergyAndParticleData(smf->pkd->Cool, &p1->CoolParticleHot, &E,
                    smf->dHotInitCodeDensity, smf->dHotInitTemp, p1->fMetals);
                TYPESet(p1,TYPE_TWOPHASE);
                }
            }
		else {
			p1->uDotFB /= p1->fMass;
		}
#else
        p1->uDotFB /= p1->fMass;
#endif
        p1->uDotESF /= p1->fMass;
        p1->fMetals /= p1->fMass;
        p1->fMFracIron /= p1->fMass;
        p1->fMFracOxygen /= p1->fMass;
        }

    }

void initStarClusterForm(void *p)
{
    TYPEReset( ((PARTICLE *) p), TYPE_STARFORM );
	}

void combStarClusterForm(void *p1,void *p2)
    {
    if (TYPETest( ((PARTICLE *) p2), TYPE_STARFORM ) && !TYPETest( ((PARTICLE *) p1), TYPE_STARFORM ) ) {
		TYPESet( ((PARTICLE *) p1), TYPE_STARFORM );
        StarClusterFormiOrder(((PARTICLE *) p1)) = StarClusterFormiOrder(((PARTICLE *) p1));
        }
    }

void StarClusterForm(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
	int i;
	PARTICLE *q;
    double rs;
    double dv2,Etherm,Mass;
	double vFac,dvx,dvy,dvz;

    vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

    dv2 = 0;
    Etherm = 0;
    Mass = 0;
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
#if SCFSMOOTHED
		double r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rs *= q->fMass;
#else
        rs = q->fMass;
#endif

		dvx = (-p->vPred[0] + q->vPred[0])*vFac;
		dvy = (-p->vPred[1] + q->vPred[1])*vFac;
		dvz = (-p->vPred[2] + q->vPred[2])*vFac;
#ifdef SLIDING_PATCH						/* Gas patch: boost dvy for radial wrap */
	    if (dx < 0.0 && (p->r[0] - q->r[0]) > 0.0)		/* if q particle is outside +x boundary */
			dvy += -1.5*smf->PP.dOrbFreq*smf->PP.dWidth*vFac;
	    else if (dx > 0.0 && (p->r[0] - q->r[0]) < 0.0)	/* if q particle is outside -x boundary */
	        dvy -= -1.5*smf->PP.dOrbFreq*smf->PP.dWidth*vFac;
#endif
        dv2 +=  (dvx*dvx+dvy*dvy*dvz*dvz)*rs;
        Etherm += q->uPred*rs;  /* Note: only first 3 dof count, not rotations, not P_floor etc... */
        Mass += rs;
	    }

    if (0.5*dv2+Etherm <= smf->dStarClusterRatio*Mass*Mass/sqrt(p->fBall2)) { // G=1, default ratio 0.5*3/5
        printf("SCF: %8d %g %g %g Cluster BOUND ENOUGH: FORM\n",p->iOrder,0.5*dv2,Etherm,smf->dStarClusterRatio*Mass*Mass/sqrt(p->fBall2));
        for (i=0;i<nSmooth;++i) {
            q = nnList[i].pPart;
            TYPESet(q,TYPE_STARFORM);
            StarClusterFormiOrder(q) = p->iOrder;
            }
        }
    else {
        printf("SCF: %8d %g %g %g Cluster NOT BOUND ENOUGH\n",p->iOrder,0.5*dv2,Etherm,smf->dStarClusterRatio*Mass*Mass/sqrt(p->fBall2));
        if (StarClusterFormiOrder(p) == -1) TYPEReset(p,TYPE_STARFORM); /* Allow for other cluster to use this particle */
        }

    p->fBall2 = StarClusterFormfBall2Save(p); /* Restore fBall2 */
#endif
}

#endif /* STARFORM */

#ifdef SIMPLESF
void initSimpleSF_Feedback(void *p1)
{
    /*
     * Zero out accumulated quantities.
     */
    ((PARTICLE *)p1)->u = 0.0;
    ((PARTICLE *)p1)->fMetals = 0.0;
    ((PARTICLE *)p1)->fTimeForm= 0.0;
    }

void combSimpleSF_Feedback(void *p1,void *p2)
{
    ((PARTICLE *)p1)->u += ((PARTICLE *)p2)->u;
    ((PARTICLE *)p1)->fMetals += ((PARTICLE *)p2)->fMetals;
    ((PARTICLE *)p1)->fTimeForm += ((PARTICLE *)p2)->fTimeForm;
    }

void SimpleSF_Feedback(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs,rstot,fNorm_u,fNorm_t;
	int i;

	assert(TYPETest(p, TYPE_STAR));
	ih2 = 4.0/BALL2(p);
	rstot = 0;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rstot += rs;
        }

	fNorm = 1./rstot;
	fNorm_u = fNorm*p->fMass*p->fESN;
	assert(fNorm_u > 0.0);

	fNorm_t = fNorm*p->uDotPdV; /* p->PdV store the cooling delay dtCoolingShutoff */
	assert(fNorm_t > 0.0);

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
	    assert(TYPETest(q, TYPE_GAS));
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		q->u += rs*fNorm_u/q->fMass;
		q->fMetals += rs*fNorm;
		q->fTimeForm += rs*fNorm_t;
        }
	}

#endif /* SIMPLESF */

#endif /* GASOLINE */

#ifdef COLLISIONS

void initFindRejects(void *p)
{
	((PARTICLE *)p)->dtCol = 0.0;
	}

void combFindRejects(void *p1,void *p2)
{
	if (((PARTICLE *)p2)->dtCol < 0.0) ((PARTICLE *)p1)->dtCol = -1.0;
	}

void FindRejects(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 ** Checks "nSmooth" neighbours of "p" for overlap (physical, Hill
	 ** radius, or a combination). When the particles in question are
	 ** of different "size", the one with the smallest size is rejected
	 ** first, otherwise the one with the higher iOrder is rejected.
	 ** Note that only planetesimal neighbours not already rejected
	 ** are considered. This procedure uses a combiner cache so that
	 ** the neighbours of "p" can be flagged with maximum efficiency.
	 */

	PARTICLE *pn;
	double r,rn,sr;
	int i;

#ifndef SLIDING_PATCH
	double a=0.0,r2,v2,an,rh;
#endif

	if (p->dtCol < 0.0) return;

	r = RADIUS(p); /* radius = 2 * softening */

#ifndef SLIDING_PATCH /*DEBUG really should handle this more generally*/
	if (smf->dCentMass > 0.0) {
		r2 = p->r[0]*p->r[0] + p->r[1]*p->r[1] + p->r[2]*p->r[2];
		v2 = p->v[0]*p->v[0] + p->v[1]*p->v[1] + p->v[2]*p->v[2];
		assert(r2 > 0.0); /* particle must not be at origin */
		a = 2/sqrt(r2) - v2/(smf->dCentMass + p->fMass);
		assert(a != 0.0); /* can't handle parabolic orbits */
		a = 1/a;
		}
#endif

	for (i=0;i<nSmooth;i++) {
		pn = nnList[i].pPart;
		if (pn->iOrder == p->iOrder || pn->dtCol < 0.0) continue;
		rn = RADIUS(pn);
#ifndef SLIDING_PATCH /*DEBUG as above*/
		if (smf->dCentMass > 0.0) {
			r2 = pn->r[0]*pn->r[0] + pn->r[1]*pn->r[1] + pn->r[2]*pn->r[2];
			v2 = pn->v[0]*pn->v[0] + pn->v[1]*pn->v[1] + pn->v[2]*pn->v[2];
			assert(r2 > 0.0);
			an = 2/sqrt(r2) - v2/(smf->dCentMass + pn->fMass);
			assert(an != 0.0);
			an = 1/an;
			rh = pow((p->fMass + pn->fMass)/(3*smf->dCentMass),1.0/3)*(a+an)/2;
			if (rh > r) r = rh;
			if (rh > rn) rn = rh;
			}
#endif
		if (rn > r || (rn == r && pn->iOrder < p->iOrder)) continue;
		sr = r + rn;
		if (nnList[i].fDist2 <= sr*sr) pn->dtCol = -1.0; /* cf REJECT() macro */
	}
}

void
_CheckForCollapse(PARTICLE *p,double dt,double rdotv,double r2,SMF *smf)
{
	/*
	 ** Sets bTinyStep flag of particle "p" to 1 if collision time "dt"
	 ** represents a fractional motion of less than "smf->dCollapseLimit" X
	 ** the current separation distance ("rdotv" is the dot product of the
	 ** relative position and relative velocity; "r2" is the square of the
	 ** distance). This routine should only be called by CheckForCollision().
	 */

	double dRatio;

	dRatio = rdotv*(p->dtPrevCol - dt)/r2;
	if (!smf->bFixCollapse)
		assert(dRatio > 0.0);
	if (dRatio < smf->dCollapseLimit) {
#if (INTERNAL_WARNINGS)
		static int bGiveWarning = 1;
		if (bGiveWarning) {
			(void) fprintf(stderr,"WARNING [T=%e]: Tiny step %i & %i "
						   "(dt=%.16e, dRatio=%.16e)\n",smf->dTime,
						   p->iOrder,p->iOrderCol,dt,dRatio);
#if (INTERNAL_WARNINGS_ONCE)
			bGiveWarning = 0;
#endif
			}
#endif /* INTERNAL_WARNINGS */
		p->bTinyStep = 1;
		}
}


void FindBinary(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/* This subroutine looks for binaries. Search only those particle
	 * on rungs higher than iMinBinaryRung. The particles have been
	 * predetermined to be on high rungs which may result from being
	 * bound. The particles are considered bound if they have a total
	 * energy < 0.
	 */

	PARTICLE *pn;
	int i,j;
	FLOAT *x,*v,v_circ2,ke=0,pe,r=0,vsq=0;

	x=malloc(3*sizeof(FLOAT));
	v=malloc(3*sizeof(FLOAT));

	for (i=0;i<nSmooth;i++) {
	  pn=nnList[i].pPart;
	  if (p == pn) continue;
	  /* Now transform to p's rest frame */
	  for (j=0;j<3;j++) {
	    v[j]=p->v[j] - pn->v[j];
	    x[j]=p->r[j] - pn->r[j];
	    vsq+=v[j]*v[j];
	    ke+=0.5*vsq;
	    r+=x[j]*x[j];
	  }
	  r=sqrt(r);

#ifdef SLIDING_PATCH
	  if (p->r[0] > pn->r[0] && nnList[i].dx < 0)
	    v[1] += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
	  else if (p->r[0] < pn->r[0] && nnList[i].dx > 0)
	    v[1] -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif

	  /* First cut: Is the particle bound? */
	  pe=-p->fMass/r;
	  /* WARNING! Here I am overloading the dtCol field. In this
	   ** situation I am going to fill it with the binding energy of
	   ** the binary. This is used later in pst/pkdFindTightestBinary.
	   */
	  p->dtCol = FLOAT_MAXVAL;
	  if ((ke+pe) < 0 ) {
	    v_circ2=ke-pe;
	    /* This is quick, but not optimal. We need to be sure we don't
	       automatically merge any bound particles, as those on highly
	       eccentric orbits may be distrupted at apocenter. Therefore
	       we assume that these particles will only reach this point of
	       the code near pericenter (due to iMinBinaryRung), and we can
	       safely exclude particles that do not meet the following
	       criterion. Some fiddling with iMinBinaryRung and dMaxBinaryEcc
	       may be necessary to acheive an appropriate balance of merging.
	    */
	    if (vsq < sqrt((1+smf->dMaxBinaryEcc)/(1-smf->dMaxBinaryEcc))*v_circ2) {
	      p->dtCol = ke+pe;
	      p->iOrderCol = pn->iOrder;
	    }
	  }
	}
	free((void *) x);
	free((void *) v);
}

void CheckForCollision(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 ** Checks whether particle "p" will collide with any of its "nSmooth"
	 ** nearest neighbours in "nnList" during drift interval smf->dStart to
	 ** smf->dEnd, relative to current time.  If any collisions are found,
	 ** the relative time to the one that will occur first is noted in
	 ** p->dtCol along with the iOrder of the collider in p->iOrderCol.
	 ** Note that it can happen that only one particle of a colliding pair
	 ** will "know" about the collision, since the other may be inactive.
	 */

	PARTICLE *pn;
	FLOAT vx,vy,vz,rdotv,v2,sr,dr2,D,dt;
	int i;

#ifdef AGGS
	FLOAT qx,qy,qz;
#endif

    assert(TYPEQueryACTIVE(p)); /* just to be sure */

	p->dtCol = DBL_MAX; /* initialize */
	p->iOrderCol = -1;
	p->bTinyStep = 0;
	if (smf->dStart == 0.0) { /* these are set in PutColliderInfo() */
		p->dtPrevCol = 0.0;
		p->iPrevCol = INT_MAX;
		}

	for (i=0;i<nSmooth;i++) {
		pn = nnList[i].pPart;
		/*
		 ** Consider all valid particles, even if inactive.  We can skip
		 ** this neighbour if it's the last particle we collided with in
		 ** this search interval, eliminating the potential problem of
		 ** "ghost" particles colliding due to roundoff error.  Note that
		 ** iPrevCol must essentially be reset after each KICK.
		 */
		if (pn == p || pn->iOrder < 0 || pn->iOrder == p->iPrevCol) continue;
#ifdef AGGS
		/*
		 ** Do not consider collisions between particles in the same
		 ** aggregate.
		 */
		if (IS_AGG(p) && IS_AGG(pn) && AGG_IDX(p) == AGG_IDX(pn)) continue;
#endif
		vx = p->v[0] - pn->v[0];
		vy = p->v[1] - pn->v[1];
		vz = p->v[2] - pn->v[2];
#ifdef AGGS
		/* centripetal terms (vx,vy,vz contain omega x r term for aggs) */
		qx = qy = qz = 0.0;
		if (IS_AGG(p)) {
			qx = p->a[0];
			qy = p->a[1];
			qz = p->a[2];
			}
		if (IS_AGG(pn)) {
			qx -= p->a[0];
			qy -= p->a[1];
			qz -= p->a[2];
			}
#endif
#ifdef SLIDING_PATCH
		if (p->r[0] > pn->r[0] && nnList[i].dx < 0.0)
			vy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
		else if (p->r[0] < pn->r[0] && nnList[i].dx > 0.0)
			vy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
		rdotv = nnList[i].dx*vx + nnList[i].dy*vy + nnList[i].dz*vz;
		if (rdotv >= 0) /*DEBUG not guaranteed for aggregates?*/
			continue; /* skip if particles not approaching */
		v2 = vx*vx + vy*vy + vz*vz;
#ifdef AGGS
		v2 += nnList[i].dx*qx + nnList[i].dy*qy + nnList[i].dz*qz;
#endif
#ifdef RUBBLE_ZML
		/*
		 ** If both particles are planetesimals increase sr by a factor
		 ** of 2 to redusce missed collisions. Inflate radius so time step
		 ** is dropped to lowest rung. The factor 2*1.4 comes from assuming
		 ** v_rms = v_esc - fg = 1+(v_esc/v_rms)^2 = 2.
		 */
		if (p->iColor == PLANETESIMAL && pn->iColor == PLANETESIMAL &&
			p->iRung == 0 && pn->iRung == 0)
			sr = 5*(p->fSoft + pn->fSoft); /*experimenting with expansion*/
/*			sr = 2.8*(p->fSoft + pn->fSoft);*/ /* softening = 0.5 * particle radius */
		else
			sr = 2*(p->fSoft + pn->fSoft); /* softening = 0.5 * particle radius */
#else
		sr = RADIUS(p) + RADIUS(pn); /* radius = twice softening */
#endif
		dr2 = nnList[i].fDist2 - sr*sr; /* negative ==> overlap... */
		D = rdotv*rdotv - dr2*v2;
		if (D <= 0.0)
			continue; /* no real solutions (or graze) ==> no collision */
		D = sqrt(D);
		/* if rdotv < 0, only one valid root possible */
		dt = (- rdotv - D)/v2; /* minimum time to surface contact */
		/*
		 ** Normally there should be no touching or overlapping particles
		 ** at the start of the step.  But inelastic collapse and other
		 ** numerical problems may make it necessary to relax this...
		 */
		if (smf->dStart == 0.0 && dt <= 0.0) {
			if (smf->bFixCollapse) {
#if (INTERNAL_WARNINGS)
				static int bGiveWarning1 = 1,bGiveWarning2 = 1;
				FLOAT fOverlap = 1.0 - sqrt(nnList[i].fDist2)/sr;
				if (bGiveWarning1 && p->iOrder < pn->iOrder) {
					(void) fprintf(stderr,"WARNING [T=%e]: "
								   "POSITION FIX %i & %i D=%e dt=%e\n",
								   smf->dTime,p->iOrder,pn->iOrder,D,dt);
#if (INTERNAL_WARNINGS_ONCE)
					bGiveWarning1 = 0;
#endif
					}
				if (bGiveWarning2 && p->iOrder < pn->iOrder && fOverlap > 0.01) {
					(void) fprintf(stderr,"WARNING [T=%e]: "
								   "LARGE OVERLAP %i & %i (%g%%)\n",
								   smf->dTime,p->iOrder,pn->iOrder,100*fOverlap);
#if (INTERNAL_WARNINGS_ONCE)
					bGiveWarning2 = 0;
#endif
					}
#endif /* INTERNAL_WARNINGS */
				if (dt < p->dtCol) { /* take most negative */
					p->dtCol = dt;
					p->iOrderCol = pn->iOrder;
					continue;
					}
				}
			else {
				(void) fprintf(stderr,"OVERLAP [T=%e]:\n"
							   "%i (r=%g,%g,%g,iRung=%i) &\n"
							   "%i (rn=%g,%g,%g,iRung=%i)\n"
							   "fDist=%g v=%g,%g,%g v_mag=%g\n"
							   "rv=%g sr=%g sep'n=%g D=%g dt=%g\n",
							   smf->dTime,
							   p->iOrder,p->r[0],p->r[1],p->r[2],p->iRung,
							   pn->iOrder,pn->r[0],pn->r[1],pn->r[2],pn->iRung,
							   sqrt(nnList[i].fDist2),vx,vy,vz,sqrt(vx*vx+vy*vy+vz*vz),
							   rdotv,sr,sqrt(nnList[i].fDist2) - sr,D,dt);
				assert(0); /* particle not allowed to touch or overlap initially */
				}
			} /* if overlap */
		/* finally, decide if this collision should be stored */
		if (dt > smf->dStart && dt <= smf->dEnd) {
			if (dt > p->dtCol) continue; /* skip if this one happens later */
			assert(dt < p->dtCol); /* can't handle simultaneous collisions */
			p->dtCol = dt;
			p->iOrderCol = pn->iOrder;
			/*DEBUG rdotv slightly different for aggs in following---ok?*/
			if (smf->dCollapseLimit > 0.0)
				_CheckForCollapse(p,dt,rdotv,nnList[i].fDist2,smf);
#ifdef RUBBLE_ZML
			/*
			 ** At the start of the top step we need to know if any
			 ** *planetesimals* (i.e., not rubble pieces) are predicted
			 ** to collide during the top step so that we can demote
			 ** them to the lowest timestep rungs (since the rubble
			 ** pieces need to be treated much more carefully, and the
			 ** collision will generally occur sometime in the middle
			 ** of the step).  Note that we only check at the beginning
			 ** of the step (i.e. when dStart is exactly zero), since
			 ** there are no resmoothing circumstances that can lead
			 ** two planetesimals to collide that have *not* already
			 ** undergone a collision themselves during the step.
			 */
			if (smf->dStart == 0 && p->iColor == PLANETESIMAL &&
				pn->iColor == PLANETESIMAL)	p->bMayCollide = 1;
			/* note: flag is reset with call to pkdRubbleResetColFlag() */
#endif
			}
		}

#ifdef SAND_PILE
#ifdef TUMBLER
	{
	WALLS *w = &smf->walls;
	double R,ndotr,ndotv,target,dt,r2;
	int j;
	R = RADIUS(p); /* this is the particle radius (twice the softening) */

	rdotv = r2 = 0; /* to keep compiler happy */
	for (i=0;i<w->nWalls;i++) {
		/* skip check if previous collision was with this wall, *unless*
		   this is a cylindrical wall. In that case we have to check for
		   collision with a different part of the wall */
		if ( (p->iPrevCol == -i - 1) && (w->wall[i].type != 1) ) continue;
		ndotr = ndotv = 0;
		for (j=0;j<3;j++) {
			ndotr += w->wall[i].n[j] * p->r[j];
			ndotv += w->wall[i].n[j] * p->v[j];
			}
		if (w->wall[i].type == 1) {     /* cylinder wall */
			double rprime[3],vprime[3],rr=0,rv=0,vv=0,dsc;
			for (j=0;j<3;j++) {
				rprime[j] = p->r[j] - ndotr * w->wall[i].n[j];
				vprime[j] = p->v[j] - ndotv * w->wall[i].n[j];
				rr += rprime[j]*rprime[j];
				rv += rprime[j]*vprime[j];
				vv += vprime[j]*vprime[j];
				}
			if (vv == 0) continue;
			target = w->wall[i].radius - R;
			dsc = rv*rv - rr*vv + target*target*vv;
			if (dsc < 0) continue;	/* no intersection with cylinder */
			/* if dsc > 0, there are two values of dt that satisfy the
			   equations. Whether we start inside or outside, we always
			   want the larger so we're colliding from the inside. */
			dt = ( sqrt(dsc) - rv ) / vv;
			if (dt <= 0) continue; /* should do an overlap check here */
			if (smf->dCollapseLimit) {
				if (rv > 0) {
					rdotv = (1 - target/sqrt(rr))*rv;
					r2 = (target - sqrt(rr));
					}
				else {
					rdotv = (1 + target/sqrt(rr))*rv;
					r2 = (target + sqrt(rr));
					}
				r2 *= r2;
				}
			}
		else { /* infinite flat wall */
			double distance;
			if (ndotv == 0) continue;
			target = w->wall[i].ndotp;
			/* first check for overlap */
			if ( (smf->dStart == 0) && (fabs(target-ndotr) <= R) ) {
				(void) fprintf(stderr,"OVERLAP [T=%g]: %i & wall %i dt %g\n",
							   smf->dTime,p->iOrder,i,dt);
				(void) fprintf(stderr,"x=%f y=%f z=%f vx=%f vy=%f vz=%f\n",
							   p->r[0],p->r[1],p->r[2],p->v[0],p->v[1],p->v[2]);
				}
			/* now check for collision */
			target += (ndotr < w->wall[i].ndotp) ? -R : R;
			distance = target - ndotr;
			dt = distance / ndotv;
			if (dt <= 0) continue; /* should do an overlap check here */
			if (smf->dCollapseLimit) {
				rdotv = - distance * ndotv;
				r2 = distance*distance;
				}
			}
		assert(smf->dStart > 0.0 || dt > 0.0);
		if (dt > smf->dStart && dt <= smf->dEnd && dt < p->dtCol) {
			p->dtCol = dt;
			p->iOrderCol = -i -1;
			if (smf->dCollapseLimit)
				_CheckForCollapse(p,dt,rdotv,r2,smf);
			}
		}
	}
#else /* TUMBLER */
	{
	WALLS *w = &smf->walls;
	double R,lx,lz,r2=0,x0,z0,d,m,b,dp,ldotv,l2,st,nx,nz,l;
	int approaching;

	R = RADIUS(p);
	v2 = p->v[0]*p->v[0] + p->v[2]*p->v[2]; /* velocity relative to frame */
	for (i=0;i<w->nWalls;i++) {
		/* check endpoints first */
		approaching = 0;
		/* first endpoint */
		lx = p->r[0] - w->wall[i].x1;
		lz = p->r[2] - w->wall[i].z1;
		rdotv = lx*p->v[0] + lz*p->v[2];
		if (p->iPrevCol != -w->nWalls - i*2 - 1 && rdotv < 0) {
			approaching = 1;
			r2 = lx*lx + lz*lz;
			D = 1 - v2*(r2 - R*R)/(rdotv*rdotv);
			if (D >= 0) {
				dt = rdotv*(sqrt(D) - 1)/v2;
				assert(smf->dStart > 0.0 || dt > 0.0);
				if (dt > smf->dStart && dt <= smf->dEnd && dt < p->dtCol) {
					p->dtCol = dt;
					p->iOrderCol = -w->nWalls - i*2 - 1; /* endpt 1 encoding */
					if (smf->dCollapseLimit)
						_CheckForCollapse(p,dt,rdotv,r2,smf);
					}
				}
			}
		/* second endpoint */
		lx = p->r[0] - w->wall[i].x2;
		lz = p->r[2] - w->wall[i].z2;
		rdotv = lx*p->v[0] + lz*p->v[2];
		if (p->iPrevCol != -w->nWalls - i*2 - 2 && rdotv < 0) {
			approaching = 1;
			r2 = lx*lx + lz*lz;
			D = 1 - v2*(r2 - R*R)/(rdotv*rdotv);
			if (D >= 0) {
				dt = rdotv*(sqrt(D) - 1)/v2;
				assert(smf->dStart > 0.0 || dt > 0.0);
				if (dt > smf->dStart && dt <= smf->dEnd && dt < p->dtCol) {
					p->dtCol = dt;
					p->iOrderCol = -w->nWalls - i*2 - 2; /* endpt 2 encoding */
					if (smf->dCollapseLimit)
						_CheckForCollapse(p,dt,rdotv,r2,smf);
					}
				}
			}
		if (!approaching) continue; /* must be approaching at least 1 endpt */
		/* now check wall */
		if (p->iPrevCol == -i - 1) continue;
		lx = w->wall[i].x2 - w->wall[i].x1;
		lz = w->wall[i].z2 - w->wall[i].z1;
		if (lx == 0) { /* vertical wall */
			if (w->wall[i].x1 < p->r[0] && p->v[0] < 0) {
				d = p->r[0] - w->wall[i].x1 - R;
				dt = -d/p->v[0];
				if (smf->dCollapseLimit) rdotv = d*p->v[0];
				}
			else if (w->wall[i].x1 > p->r[0] && p->v[0] > 0) {
				d = w->wall[i].x1 - p->r[0] - R;
				dt = d/p->v[0];
				if (smf->dCollapseLimit) rdotv = -d*p->v[0];
				}
			else continue;
			x0 = 0; /* to satisfy compiler */
			z0 = p->r[2] + p->v[2]*dt;
			if (smf->dCollapseLimit) r2 = d*d;
			}
		else if (lz == 0) { /* horizontal wall */
			if (w->wall[i].z1 < p->r[2] && p->v[2] < 0) {
				d = p->r[2] - w->wall[i].z1 - R;
				dt = -d/p->v[2];
				if (smf->dCollapseLimit) rdotv = d*p->v[2];
				}
			else if (w->wall[i].z1 > p->r[2] && p->v[2] > 0) {
				d = w->wall[i].z1 - p->r[2] - R;
				dt = d/p->v[2];
				if (smf->dCollapseLimit) rdotv = -d*p->v[2];
				}
			else continue;
			x0 = p->r[0] + p->v[0]*dt;
			z0 = 0;
			if (smf->dCollapseLimit) r2 = d*d;
			}
		else { /* oblique wall */
			m = lz/lx;
			b = w->wall[i].z1 - m*w->wall[i].x1;
			dp = (m*p->r[0] - p->r[2] + b)/sqrt(1 + m*m); /* perp. distance */
			if (dp < 0) dp = -dp;
			ldotv = lx*p->v[0] + lz*p->v[2];
			l2 = lx*lx + lz*lz;
			st = sqrt(1 - ldotv*ldotv/(l2*v2)); /* sin theta */
			d = (dp - R)/st; /* travel distance until contact */
			dt = d/sqrt(v2); /* travel time */
			if (dt <= smf->dStart) continue; /* to save time */
			x0 = p->r[0] + p->v[0]*dt;
			z0 = p->r[2] + p->v[2]*dt;
			nx = lz;
			nz = -lx;
			if ((w->wall[i].x1 - x0)*nx + (w->wall[i].z1 - z0)*nz < 0) {
				nx = -nx;
				nz = -nz;
				}
			rdotv = -(nx*p->v[0] + nz*p->v[2]); /* wrong magnitude... */
			if (rdotv >= 0) continue; /* moving away */
			l = R/sqrt(l2);
			x0 += nx*l;
			z0 += nz*l;
			if (smf->dCollapseLimit) {
				rdotv = (dp - R)*rdotv*l/R; /* roundabout... */
				r2 = (dp - R)*(dp - R);
				}
			}
		/* check perpendicular contact point lies on or between endpoints */
		if ((lx < 0 && (x0 > w->wall[i].x1 || x0 < w->wall[i].x2)) ||
			(lx > 0 && (x0 < w->wall[i].x1 || x0 > w->wall[i].x2)) ||
			(lz < 0 && (z0 > w->wall[i].z1 || z0 < w->wall[i].z2)) ||
			(lz > 0 && (z0 < w->wall[i].z1 || z0 > w->wall[i].z2)))
			continue;
		if (smf->dStart == 0 && dt <= 0) {
			(void) fprintf(stderr,"OVERLAP [T=%g]: %i & wall %i dt %g\n",
					smf->dTime,p->iOrder,i,dt);
			assert(0); /* no backsteps allowed */
			}
		if (dt >= smf->dStart && dt <= smf->dEnd && dt <= p->dtCol) {
			assert(dt < p->dtCol); /* no simultaneous collisions allowed */
			p->dtCol = dt;
			p->iOrderCol = -1 - i; /* wall index encoding */
			if (smf->dCollapseLimit) _CheckForCollapse(p,dt,rdotv,r2,smf);
			}
		}
	}
#endif /* !TUMBLER */
#endif /* SAND_PILE */

	}

#endif /* COLLISIONS */

#ifdef COLLISIONS
#ifdef SLIDING_PATCH

void initFindOverlaps(void *p)
{
	((PARTICLE *)p)->dtCol = 0.0;
	}

void combFindOverlaps(void *p1,void *p2)
{
	if (((PARTICLE *)p2)->dtCol < 0.0) ((PARTICLE *)p1)->dtCol = -1.0;
	}

void FindOverlaps(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 ** Streamlined version of FindRejects() designed specifically
	 ** for the sliding patch when we want to randomize particle
	 ** data following an azimuthal boundary wrap.  As part of the
	 ** randomization procedure, we need to make sure we don't
	 ** overlap any particles.  That's what's checked for here.
	 */

	PARTICLE *pn = NULL;
	double r,rn,sr;
	int i;

	assert(nSmooth > 1); /* for now */

	assert(p->dtCol >= 0.0); /* can't already be rejected */

	assert(p->bAzWrap == 1); /* must have wrapped */

	r = RADIUS(p); /* radius = 2 * softening */

	for (i=0;i<nSmooth;i++) {
		pn = nnList[i].pPart;
		if (pn->iOrder == p->iOrder) continue;
		rn = RADIUS(pn);
		sr = r + rn;
		if (nnList[i].fDist2 <= sr*sr) p->dtCol = -1.0; /* cf REJECT() macro */
	}

	if (p->dtCol >= 0.0) p->bAzWrap = 0; /* if not rejected, do not need to regenerate */
	/*DEBUG*/if (p->bAzWrap) printf("FindOverlaps(): particle %i overlaps particle %i\n",p->iOrder,pn->iOrder);
}

#endif /* SLIDING_PATCH */
#endif /* COLLISIONS */
