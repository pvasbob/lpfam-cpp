#include <iostream>
#include <string>

#include "UNEDF.h"

bool UNEDF::use_charge_density;
bool UNEDF::use_cm_cor;
// static double Urhorho[4][8], Urhotau[4][8], UrhoDrho[4][8], Unablarho[4][8]; //! ph DME amplitudes
double UNEDF::Urhorho[4][8];   //! ph DME amplitudes
double UNEDF::Urhotau[4][8];   //! ph DME amplitudes
double UNEDF::UrhoDrho[4][8];  //! ph DME amplitudes
double UNEDF::Unablarho[4][8]; //! ph DME amplitudes

// static double UJnablarho[4][8], UrhonablaJ[4][8], UJJ[4][8];
double UNEDF::UJnablarho[4][8];
double UNEDF::UrhonablaJ[4][8];
double UNEDF::UJJ[4][8];
double UNEDF::Urhorhopr[4][8]; //! pp amplitudes

// static double UEnonstdr[2], UFnonstdr[2], URnonstdr[2];                              //! Other amplitudes
double UNEDF::UEnonstdr[2]; //! Other amplitudes
double UNEDF::UFnonstdr[2]; //! Other amplitudes
double UNEDF::URnonstdr[2]; //! Other amplitudes

// static double hbzero, sigma, e2charg, CExPar;                                        //! hbr^2/2m, DD sigma, e^2 charge, coul. exchange parameter
double UNEDF::hbzero;
double UNEDF::sigma;
double UNEDF::e2charg;
double UNEDF::CExPar; //! hbr^2/2m; DD sigma; e^2 charge; coul. exchange parameter

// static double Crho[2], Cdrho[2], Ctau[2], CrDr[2], CrdJ[2], CJ[2], CpV0[2], CpV1[2]; //! basic coupling constants
double UNEDF::Crho[2];
double UNEDF::Cdrho[2];
double UNEDF::Ctau[2];
double UNEDF::CrDr[2];
double UNEDF::CrdJ[2];
double UNEDF::CJ[2];
double UNEDF::CpV0[2];
double UNEDF::CpV1[2]; //! basic coupling constants

// static double E_NM, K_NM, SMASS_NM, RHO_NM, ASS_NM, LASS_NM, VMASS_NM, P_NM, KA_NM;
double UNEDF::E_NM;
double UNEDF::K_NM;
double UNEDF::SMASS_NM;
double UNEDF::RHO_NM;
double UNEDF::ASS_NM;
double UNEDF::LASS_NM;
double UNEDF::VMASS_NM;
double UNEDF::P_NM;
double UNEDF::KA_NM;

double UNEDF::CHrho; //! Crho(0) from the Hartree term in NM
// static bool use_DME3N_terms, use_j2terms;
bool UNEDF::use_DME3N_terms;
bool UNEDF::use_j2terms;

// static int DMEorder, DMElda, use_TMR_pairing;
int UNEDF::DMEorder;
int UNEDF::DMElda;
int UNEDF::use_TMR_pairing;

// static double mpi, gA, fpi, c1, c3, c4, cd, ce, LambdaX;
double UNEDF::mpi;
double UNEDF::gA;
double UNEDF::fpi;
double UNEDF::c1;
double UNEDF::c3;
double UNEDF::c4;
double UNEDF::cd;
double UNEDF::ce;
double UNEDF::LambdaX;

// static bool use_INM, use_Namelist, Print_Namelist;
// static bool use_INM, use_Namelist, Print_Namelist;
bool UNEDF::use_INM;
bool UNEDF::use_Namelist;
bool UNEDF::Print_Namelist;
std::string UNEDF::FunctionalName;

// Private in Orign

// static double nuCrho[2], nuCdrho[2], nuCtau[2], nuCrDr[2]; //! basic coupling constants in natural units
double UNEDF::nuCrho[2];
double UNEDF::nuCdrho[2];
double UNEDF::nuCtau[2];
double UNEDF::nuCrDr[2]; //! basic coupling constants in natural units

// static double nuCrdJ[2], nuCJ[2], nuCpV0[2], nuCpV1[2];    //!
double UNEDF::nuCrdJ[2];
double UNEDF::nuCJ[2];
double UNEDF::nuCpV0[2];
double UNEDF::nuCpV1[2]; //!

// static double t0, t1, t2, t3, x0, x1, x2, x3, b4, b4p, te, to;
double UNEDF::t0;
double UNEDF::t1;
double UNEDF::t2;
double UNEDF::t3;
double UNEDF::x0;
double UNEDF::x1;
double UNEDF::x2;
double UNEDF::x3;
double UNEDF::b4;
double UNEDF::b4p;
double UNEDF::te;
double UNEDF::to;

// static double nuLambda, nufpi;   //! parameters associated to natural units
double UNEDF::nuLambda;
double UNEDF::nufpi; //! parameters associated to natural units

// static double Cnrho[2], CJdr[2]; //! hidden and always zero
double UNEDF::Cnrho[2];
double UNEDF::CJdr[2]; //! hidden and always zero

int UNEDF::i_cut;  //! dmeorder: -1=Standard Skyrme, 0=LO, 1=NLO, 2=N2LO
double UNEDF::eps; //! dmelda: 0=Kf-LDA, 1=CB-LDA
// static double kfconst, CK;       //! (3Pi^2/2)^(1/3)
double UNEDF::kfconst;
double UNEDF::CK; //! (3Pi^2/2)^(1/3)

double UNEDF::mevfm;
// static double qqrho[2], qqtau[2], nrho2[2], lrho[2];
double UNEDF::qqrho[2];
double UNEDF::qqtau[2];
double UNEDF::nrho2[2];
double UNEDF::lrho[2];

// static double mpi2, fpi2, fpi4, gA2, gA4, gA6, CHartree;
double UNEDF::mpi2;
double UNEDF::fpi2;
double UNEDF::fpi4;
double UNEDF::gA2;
double UNEDF::gA4;
double UNEDF::gA6;
double UNEDF::CHartree;

// static double arhorho, brhorho, arhodrho, brhodrho, arhotau, brhotau, ajj, bjj, adrdr, bdrdr;
double UNEDF::arhorho;
double UNEDF::brhorho;
double UNEDF::arhodrho;
double UNEDF::brhodrho;
double UNEDF::arhotau;
double UNEDF::brhotau;
double UNEDF::ajj;
double UNEDF::bjj;
double UNEDF::adrdr;
double UNEDF::bdrdr;

// static double darhorho, dbrhorho, darhodrho, dbrhodrho, darhotau, dbrhotau, dajj, dbjj, dadrdr, dbdrdr;
double UNEDF::darhorho;
double UNEDF::dbrhorho;
double UNEDF::darhodrho;
double UNEDF::dbrhodrho;
double UNEDF::darhotau;
double UNEDF::dbrhotau;
double UNEDF::dajj;
double UNEDF::dbjj;
double UNEDF::dadrdr;
double UNEDF::dbdrdr;

// static double ddarhodrho, ddbrhodrho, ddarhotau, ddbrhotau, ddarhorho, ddbrhorho;
double UNEDF::ddarhodrho;
double UNEDF::ddbrhodrho;
double UNEDF::ddarhotau;
double UNEDF::ddbrhotau;
double UNEDF::ddarhorho;
double UNEDF::ddbrhorho;

// static double hrho0rho0, hrho1rho1, hdr0dr0, hdr1dr1, hrho0Drho0, hrho1Drho0,
//     hrho1Drho1, hrho0tau0, hrho1tau0, hrho1tau1, hJ0dr0, hrho0DJ0, hJ1dr1, hrho1DJ1,
//     hJ0dr1, hrho1DJ0, hJ1dr0, hJ0J0, hJ0J1, hJ1J1;
double UNEDF::hrho0rho0;
double UNEDF::hrho1rho1;
double UNEDF::hdr0dr0;
double UNEDF::hdr1dr1;
double UNEDF::hrho0Drho0;
double UNEDF::hrho1Drho0;
double UNEDF::hrho1Drho1;
double UNEDF::hrho0tau0;
double UNEDF::hrho1tau0;
double UNEDF::hrho1tau1;
double UNEDF::hJ0dr0;
double UNEDF::hrho0DJ0;
double UNEDF::hJ1dr1;
double UNEDF::hrho1DJ1;
double UNEDF::hJ0dr1;
double UNEDF::hrho1DJ0;
double UNEDF::hJ1dr0;
double UNEDF::hJ0J0;
double UNEDF::hJ0J1;
double UNEDF::hJ1J1;

// static double dhrho0rho0, dhrho1rho1, dhdr0dr0, dhdr1dr1, dhrho0Drho0,
//     dhrho1Drho0, dhrho1Drho1, dhrho0tau0, dhrho1tau0, dhrho1tau1, dhJ0dr0, dhrho0DJ0,
//     dhJ1dr1, dhrho1DJ1, dhJ0dr1, dhrho1DJ0, dhJ1dr0, dhJ0J0, dhJ0J1, dhJ1J1;
double UNEDF::dhrho0rho0;
double UNEDF::dhrho1rho1;
double UNEDF::dhdr0dr0;
double UNEDF::dhdr1dr1;
double UNEDF::dhrho0Drho0;
double UNEDF::dhrho1Drho0;
double UNEDF::dhrho1Drho1;
double UNEDF::dhrho0tau0;
double UNEDF::dhrho1tau0;
double UNEDF::dhrho1tau1;
double UNEDF::dhJ0dr0;
double UNEDF::dhrho0DJ0;
double UNEDF::dhJ1dr1;
double UNEDF::dhrho1DJ1;
double UNEDF::dhJ0dr1;
double UNEDF::dhrho1DJ0;
double UNEDF::dhJ1dr0;
double UNEDF::dhJ0J0;
double UNEDF::dhJ0J1;
double UNEDF::dhJ1J1;

// static double ddhrho0rho0, ddhrho1rho1, ddhrho0Drho0, ddhrho1Drho0,
//     ddhrho1Drho1, ddhrho0tau0, ddhrho1tau0, ddhrho1tau1;
double UNEDF::ddhrho0rho0;
double UNEDF::ddhrho1rho1;
double UNEDF::ddhrho0Drho0;
double UNEDF::ddhrho1Drho0;
double UNEDF::ddhrho1Drho1;
double UNEDF::ddhrho0tau0;
double UNEDF::ddhrho1tau0;
double UNEDF::ddhrho1tau1;

// static double ctr0r0[3][3][33], ctr1r1[3][3][33], ctdr0dr0[3][3][33], ctdr1dr1[3][3][33], //! coefficients for 3N part
//     ctr0Dr0[3][3][33], ctr1Dr0[3][3][33], ctr1Dr1[3][3][33], ctr0t0[3][3][33], ctr1t0[3][3][33], ctr1t1[3][3][33], ctJ0dr0[3][3][33], ctr0dJ0[3][3][33], ctJ1dr1[3][3][33],
//     ctr1dJ1[3][3][33], ctJ0dr1[3][3][33], ctr1dJ0[3][3][33], ctJ1dr0[3][3][33], ctJ0J0[3][3][33], ctJ0J1[3][3][33], ctJ1J1[3][3][33];
double UNEDF::ctr0r0[3][3][33];
double UNEDF::ctr1r1[3][3][33];
double UNEDF::ctdr0dr0[3][3][33];
double UNEDF::ctdr1dr1[3][3][33]; //! coefficients for 3N part
double UNEDF::ctr0Dr0[3][3][33];
double UNEDF::ctr1Dr0[3][3][33];
double UNEDF::ctr1Dr1[3][3][33];
double UNEDF::ctr0t0[3][3][33];
double UNEDF::ctr1t0[3][3][33];
double UNEDF::ctr1t1[3][3][33];
double UNEDF::ctJ0dr0[3][3][33];
double UNEDF::ctr0dJ0[3][3][33];
double UNEDF::ctJ1dr1[3][3][33];
double UNEDF::ctr1dJ1[3][3][33];
double UNEDF::ctJ0dr1[3][3][33];
double UNEDF::ctr1dJ0[3][3][33];
double UNEDF::ctJ1dr0[3][3][33];
double UNEDF::ctJ0J0[3][3][33];
double UNEDF::ctJ0J1[3][3][33];
double UNEDF::ctJ1J1[3][3][33];

// static double u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12;
double UNEDF::u2;
double UNEDF::u3;
double UNEDF::u4;
double UNEDF::u5;
double UNEDF::u6;
double UNEDF::u7;
double UNEDF::u8;
double UNEDF::u9;
double UNEDF::u10;
double UNEDF::u11;
double UNEDF::u12;

// static double ual, lual, atu, asqu, asqu4;
double UNEDF::ual;
double UNEDF::lual;
double UNEDF::atu;
double UNEDF::asqu;
double UNEDF::asqu4;

double UNEDF::acoord = 0.50;
double UNEDF::ac2 = 4.0 * (pow(acoord, 2) - acoord + 0.50);
double UNEDF::ac3 = 2.0 * (pow(acoord, 2) - acoord + 0.50);
//!
// static double A1_1, A1_2, A1_3, A1_4, A1_5, b1_1, b1_2, b1_3, b1_4, b1_5;
double UNEDF::A1_1;
double UNEDF::A1_2;
double UNEDF::A1_3;
double UNEDF::A1_4;
double UNEDF::A1_5;
double UNEDF::b1_1;
double UNEDF::b1_2;
double UNEDF::b1_3;
double UNEDF::b1_4;
double UNEDF::b1_5;

// static double A3_1, A3_2, A3_3, A3_4, A3_5, b3_1, b3_2, b3_3, b3_4, b3_5;
double UNEDF::A3_1;
double UNEDF::A3_2;
double UNEDF::A3_3;
double UNEDF::A3_4;
double UNEDF::A3_5;
double UNEDF::b3_1;
double UNEDF::b3_2;
double UNEDF::b3_3;
double UNEDF::b3_4;
double UNEDF::b3_5;

// static double h0mpi6, h0mpi6c1, h0mpi6c3, h0mpi6NM, h0mpi6c1NM, h0mpi6c3NM;
double UNEDF::h0mpi6;
double UNEDF::h0mpi6c1;
double UNEDF::h0mpi6c3;
double UNEDF::h0mpi6NM;
double UNEDF::h0mpi6c1NM;
double UNEDF::h0mpi6c3NM;

void UNEDF::read_UNEDF_NAMELIST()
{

    //   std::string FunctionalName = "ORIGINAL_SKMS";
    //   FunctionalName,DMEorder,DMElda,use_INM,hbzero,use_TMR_pairing
    //   Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,sigma,CpV0,CpV1,e2charg, &
    //   E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM, &
    //   mpi,gA,fpi,c1,c3,c4,cd,ce,LambdaX, &
    //   use_cm_cor,use_charge_density,use_DME3N_terms,use_j2terms,CExPar, &
    //   Print_Namelist
    //

    FunctionalName = "ORIGINAL_SKMS";
    use_cm_cor = true;
    use_j2terms = true;
    DMEorder = -1;
    DMElda = 0;
    use_INM = false;
    hbzero = 20.7339830;
    use_TMR_pairing = 0;
    Crho[0] = -991.875000000000000;
    Crho[1] = 390.137499999999989;
    Cdrho[0] = 974.687500000000000;
    Cdrho[1] = -324.895833333333314;
    Ctau[0] = 34.6875000000000000;
    Ctau[1] = -34.0625000000000000;
    CrDr[0] = -68.2031250000000000;
    CrDr[1] = 17.1093750000000000;
    CrdJ[0] = -97.5000000000000000;
    CrdJ[1] = -32.5000000000000000;
    CJ[0] = 0.00000000000000000;
    CJ[1] = 0.00000000000000000;
    CpV0[0] = -160.000000000000000;
    CpV0[1] = -180.000000000000000;
    CpV1[0] = 0.00000000000000000;
    CpV1[1] = 0.00000000000000000;
    sigma = 0.166666666666666713;
    hbzero = 20.7339829999999985;
    e2charg = 1.43997839999999999;
    CExPar = 1.00000000000000000;

    mpi = 138.03 / 197.3;
    fpi = 92.4 / 197.3;
    gA = 1.29;
    c1 = -0.81 / 1000.0 * 197.3;
    c3 = -3.4 / 1000.0 * 197.3;
    c4 = 3.4 / 1000.0 * 197.3;
    cd = -2062.0 / 1000.0;
    ce = -625.0 / 1000.0;
    LambdaX = 700.0 / 197.3;
}

void UNEDF::printUNEDF()
{
    // std::cout << CExPar << " " << Crho[0] << " " << Cdrho[0] << std::endl;
    std::cout << UNEDF::FunctionalName << std::endl;
    std::cout << "CpV0[0]" << CpV0[0] << "CpV0[1]" << CpV0[1] << std::endl;
}
