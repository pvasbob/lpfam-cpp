#ifndef UNEDF_H
#define UNEDF_H

#include <vector>
#include <complex>
#include <string>

class UNEDF
{
public:
    void read_UNEDF_NAMELIST();
    void printUNEDF();

protected:
    //!===================================================================================================================================
    // !Module UNEDF
    // !===================================================================
    //   !------------------------------------------------
    // static bool use_charge_density, use_cm_cor;
    static bool use_charge_density;
    static bool use_cm_cor;
    // static double Urhorho[4][8], Urhotau[4][8], UrhoDrho[4][8], Unablarho[4][8]; //! ph DME amplitudes
    static double Urhorho[4][8];   //! ph DME amplitudes
    static double Urhotau[4][8];   //! ph DME amplitudes
    static double UrhoDrho[4][8];  //! ph DME amplitudes
    static double Unablarho[4][8]; //! ph DME amplitudes

    // static double UJnablarho[4][8], UrhonablaJ[4][8], UJJ[4][8];
    static double UJnablarho[4][8];
    static double UrhonablaJ[4][8];
    static double UJJ[4][8];

    static double Urhorhopr[4][8]; //! pp amplitudes

    // static double UEnonstdr[2], UFnonstdr[2], URnonstdr[2];                              //! Other amplitudes
    static double UEnonstdr[2]; //! Other amplitudes
    static double UFnonstdr[2]; //! Other amplitudes
    static double URnonstdr[2]; //! Other amplitudes

    // static double hbzero, sigma, e2charg, CExPar;                                        //! hbr^2/2m, DD sigma, e^2 charge, coul. exchange parameter
    static double hbzero;
    static double sigma;
    static double e2charg;
    static double CExPar; //! hbr^2/2m; DD sigma; e^2 charge; coul. exchange parameter

    // static double Crho[2], Cdrho[2], Ctau[2], CrDr[2], CrdJ[2], CJ[2], CpV0[2], CpV1[2]; //! basic coupling constants
    static double Crho[2];
    static double Cdrho[2];
    static double Ctau[2];
    static double CrDr[2];
    static double CrdJ[2];
    static double CJ[2];
    static double CpV0[2];
    static double CpV1[2]; //! basic coupling constants

    // static double E_NM, K_NM, SMASS_NM, RHO_NM, ASS_NM, LASS_NM, VMASS_NM, P_NM, KA_NM;
    static double E_NM;
    static double K_NM;
    static double SMASS_NM;
    static double RHO_NM;
    static double ASS_NM;
    static double LASS_NM;
    static double VMASS_NM;
    static double P_NM;
    static double KA_NM;

    static double CHrho; //! Crho(0) from the Hartree term in NM
    // static bool use_DME3N_terms, use_j2terms;
    static bool use_DME3N_terms;
    static bool use_j2terms;

    // static int DMEorder, DMElda, use_TMR_pairing;
    static int DMEorder;
    static int DMElda;
    static int use_TMR_pairing;

    // static double mpi, gA, fpi, c1, c3, c4, cd, ce, LambdaX;
    static double mpi;
    static double gA;
    static double fpi;
    static double c1;
    static double c3;
    static double c4;
    static double cd;
    static double ce;
    static double LambdaX;

    // static bool use_INM, use_Namelist, Print_Namelist;
    // static bool use_INM, use_Namelist, Print_Namelist;
    static bool use_INM;
    static bool use_Namelist;
    static bool Print_Namelist;

    static std::string FunctionalName;
    //!
    //! === PRIVATE VARIABLES ===
    //!
protected:
    // static double nuCrho[2], nuCdrho[2], nuCtau[2], nuCrDr[2]; //! basic coupling constants in natural units
    static double nuCrho[2];
    static double nuCdrho[2];
    static double nuCtau[2];
    static double nuCrDr[2]; //! basic coupling constants in natural units

    // static double nuCrdJ[2], nuCJ[2], nuCpV0[2], nuCpV1[2];    //!
    static double nuCrdJ[2];
    static double nuCJ[2];
    static double nuCpV0[2];
    static double nuCpV1[2]; //!

    // static double t0, t1, t2, t3, x0, x1, x2, x3, b4, b4p, te, to;
    static double t0;
    static double t1;
    static double t2;
    static double t3;
    static double x0;
    static double x1;
    static double x2;
    static double x3;
    static double b4;
    static double b4p;
    static double te;
    static double to;

    // static double nuLambda, nufpi;   //! parameters associated to natural units
    static double nuLambda;
    static double nufpi; //! parameters associated to natural units

    // static double Cnrho[2], CJdr[2]; //! hidden and always zero
    static double Cnrho[2];
    static double CJdr[2]; //! hidden and always zero

    static int i_cut;  //! dmeorder: -1=Standard Skyrme, 0=LO, 1=NLO, 2=N2LO
    static double eps; //! dmelda: 0=Kf-LDA, 1=CB-LDA
    // static double kfconst, CK;       //! (3Pi^2/2)^(1/3)
    static double kfconst;
    static double CK; //! (3Pi^2/2)^(1/3)

    static double mevfm;
    // static double qqrho[2], qqtau[2], nrho2[2], lrho[2];
    static double qqrho[2];
    static double qqtau[2];
    static double nrho2[2];
    static double lrho[2];

    // static double mpi2, fpi2, fpi4, gA2, gA4, gA6, CHartree;
    static double mpi2;
    static double fpi2;
    static double fpi4;
    static double gA2;
    static double gA4;
    static double gA6;
    static double CHartree;

    // static double arhorho, brhorho, arhodrho, brhodrho, arhotau, brhotau, ajj, bjj, adrdr, bdrdr;
    static double arhorho;
    static double brhorho;
    static double arhodrho;
    static double brhodrho;
    static double arhotau;
    static double brhotau;
    static double ajj;
    static double bjj;
    static double adrdr;
    static double bdrdr;

    // static double darhorho, dbrhorho, darhodrho, dbrhodrho, darhotau, dbrhotau, dajj, dbjj, dadrdr, dbdrdr;
    static double darhorho;
    static double dbrhorho;
    static double darhodrho;
    static double dbrhodrho;
    static double darhotau;
    static double dbrhotau;
    static double dajj;
    static double dbjj;
    static double dadrdr;
    static double dbdrdr;

    // static double ddarhodrho, ddbrhodrho, ddarhotau, ddbrhotau, ddarhorho, ddbrhorho;
    static double ddarhodrho;
    static double ddbrhodrho;
    static double ddarhotau;
    static double ddbrhotau;
    static double ddarhorho;
    static double ddbrhorho;

    // static double hrho0rho0, hrho1rho1, hdr0dr0, hdr1dr1, hrho0Drho0, hrho1Drho0,
    //     hrho1Drho1, hrho0tau0, hrho1tau0, hrho1tau1, hJ0dr0, hrho0DJ0, hJ1dr1, hrho1DJ1,
    //     hJ0dr1, hrho1DJ0, hJ1dr0, hJ0J0, hJ0J1, hJ1J1;
    static double hrho0rho0;
    static double hrho1rho1;
    static double hdr0dr0;
    static double hdr1dr1;
    static double hrho0Drho0;
    static double hrho1Drho0;
    static double hrho1Drho1;
    static double hrho0tau0;
    static double hrho1tau0;
    static double hrho1tau1;
    static double hJ0dr0;
    static double hrho0DJ0;
    static double hJ1dr1;
    static double hrho1DJ1;
    static double hJ0dr1;
    static double hrho1DJ0;
    static double hJ1dr0;
    static double hJ0J0;
    static double hJ0J1;
    static double hJ1J1;

    // static double dhrho0rho0, dhrho1rho1, dhdr0dr0, dhdr1dr1, dhrho0Drho0,
    //     dhrho1Drho0, dhrho1Drho1, dhrho0tau0, dhrho1tau0, dhrho1tau1, dhJ0dr0, dhrho0DJ0,
    //     dhJ1dr1, dhrho1DJ1, dhJ0dr1, dhrho1DJ0, dhJ1dr0, dhJ0J0, dhJ0J1, dhJ1J1;
    static double dhrho0rho0;
    static double dhrho1rho1;
    static double dhdr0dr0;
    static double dhdr1dr1;
    static double dhrho0Drho0;
    static double dhrho1Drho0;
    static double dhrho1Drho1;
    static double dhrho0tau0;
    static double dhrho1tau0;
    static double dhrho1tau1;
    static double dhJ0dr0;
    static double dhrho0DJ0;
    static double dhJ1dr1;
    static double dhrho1DJ1;
    static double dhJ0dr1;
    static double dhrho1DJ0;
    static double dhJ1dr0;
    static double dhJ0J0;
    static double dhJ0J1;
    static double dhJ1J1;

    // static double ddhrho0rho0, ddhrho1rho1, ddhrho0Drho0, ddhrho1Drho0,
    //     ddhrho1Drho1, ddhrho0tau0, ddhrho1tau0, ddhrho1tau1;
    static double ddhrho0rho0;
    static double ddhrho1rho1;
    static double ddhrho0Drho0;
    static double ddhrho1Drho0;
    static double ddhrho1Drho1;
    static double ddhrho0tau0;
    static double ddhrho1tau0;
    static double ddhrho1tau1;

    // static double ctr0r0[3][3][33], ctr1r1[3][3][33], ctdr0dr0[3][3][33], ctdr1dr1[3][3][33], //! coefficients for 3N part
    //     ctr0Dr0[3][3][33], ctr1Dr0[3][3][33], ctr1Dr1[3][3][33], ctr0t0[3][3][33], ctr1t0[3][3][33], ctr1t1[3][3][33], ctJ0dr0[3][3][33], ctr0dJ0[3][3][33], ctJ1dr1[3][3][33],
    //     ctr1dJ1[3][3][33], ctJ0dr1[3][3][33], ctr1dJ0[3][3][33], ctJ1dr0[3][3][33], ctJ0J0[3][3][33], ctJ0J1[3][3][33], ctJ1J1[3][3][33];
    static double ctr0r0[3][3][33];
    static double ctr1r1[3][3][33];
    static double ctdr0dr0[3][3][33];
    static double ctdr1dr1[3][3][33]; //! coefficients for 3N part
    static double ctr0Dr0[3][3][33];
    static double ctr1Dr0[3][3][33];
    static double ctr1Dr1[3][3][33];
    static double ctr0t0[3][3][33];
    static double ctr1t0[3][3][33];
    static double ctr1t1[3][3][33];
    static double ctJ0dr0[3][3][33];
    static double ctr0dJ0[3][3][33];
    static double ctJ1dr1[3][3][33];
    static double ctr1dJ1[3][3][33];
    static double ctJ0dr1[3][3][33];
    static double ctr1dJ0[3][3][33];
    static double ctJ1dr0[3][3][33];
    static double ctJ0J0[3][3][33];
    static double ctJ0J1[3][3][33];
    static double ctJ1J1[3][3][33];

    // static double u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12;
    static double u2;
    static double u3;
    static double u4;
    static double u5;
    static double u6;
    static double u7;
    static double u8;
    static double u9;
    static double u10;
    static double u11;
    static double u12;

    // static double ual, lual, atu, asqu, asqu4;
    static double ual;
    static double lual;
    static double atu;
    static double asqu;
    static double asqu4;

    static double acoord;
    static double ac2;
    static double ac3;
    //!
    // static double A1_1, A1_2, A1_3, A1_4, A1_5, b1_1, b1_2, b1_3, b1_4, b1_5;
    static double A1_1;
    static double A1_2;
    static double A1_3;
    static double A1_4;
    static double A1_5;
    static double b1_1;
    static double b1_2;
    static double b1_3;
    static double b1_4;
    static double b1_5;

    // static double A3_1, A3_2, A3_3, A3_4, A3_5, b3_1, b3_2, b3_3, b3_4, b3_5;
    static double A3_1;
    static double A3_2;
    static double A3_3;
    static double A3_4;
    static double A3_5;
    static double b3_1;
    static double b3_2;
    static double b3_3;
    static double b3_4;
    static double b3_5;

    // static double h0mpi6, h0mpi6c1, h0mpi6c3, h0mpi6NM, h0mpi6c1NM, h0mpi6c3NM;
    static double h0mpi6;
    static double h0mpi6c1;
    static double h0mpi6c3;
    static double h0mpi6NM;
    static double h0mpi6c1NM;
    static double h0mpi6c3NM;

    //!
    //   Namelist /UNEDF_NAMELIST/ FunctionalName,DMEorder,DMElda,use_INM,hbzero,use_TMR_pairing, &
    //        Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,sigma,CpV0,CpV1,e2charg, &
    //        E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM, &
    //        mpi,gA,fpi,c1,c3,c4,cd,ce,LambdaX, &
    //        use_cm_cor,use_charge_density,use_DME3N_terms,use_j2terms,CExPar, &
    //        Print_Namelist
    //
};

#endif
