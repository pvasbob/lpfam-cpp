#ifndef HFBTHO_H
#define HFBTHO_H

#include <vector>
#include <complex>
#include <string>

#include "UNEDF.h"

class HFBTHO : public UNEDF
{
public:
    void read_HFBTHO_NAMELIST();
    void printHFBTHO();

public:
    // int n00_INI, iLST_INI, inin_INI, icou_INI;
    static int n00_INI;
    static int iLST_INI;
    static int inin_INI;
    static int icou_INI;

    // int npr_INI[3], kindhfb_INI;
    static int npr_INI[3];
    static int kindhfb_INI;

    // int keyblo1_INI, keyblo2_INI, IDEBUG_INI;
    static int keyblo1_INI;
    static int keyblo2_INI;
    static int IDEBUG_INI;

    // double b0_INI, q_INI, cdef_INI, cqad_INI;
    static double b0_INI;
    static double q_INI;
    static double cdef_INI;
    static double cqad_INI;

    static char skyrme_INI[30];
    static double epsi_INI;
    // bool Add_Pairing_INI, Print_HFBTHO_Namelist_INI, DO_FITT_INI;
    static bool Add_Pairing_INI;
    static bool Print_HFBTHO_Namelist_INI;
    static bool DO_FITT_INI;

    static int nkblo_INI[2][5];
    static bool Parity_INI;

    static int MAX_ITER_INI;
    static int keypj_INI;
    static int iproj_INI;
    static int npr1pj_INI;
    static int npr2pj_INI;

    // protected:
public:
    static char Version[6];
    static char pairing_Version[6];

    // int NNN0, mv, ivpair = 0;
    static int NNN0;
    static int mv;
    static int ivpair;

    // double gl0, gal, gl[2], vfac[2], e_pair[2], d_pair[2];
    static double gl0;
    static double gal;
    static double gl[2];
    static double vfac[2];
    static double e_pair[2];
    static double d_pair[2];

    // std::vector<std::vector<double>> wnn, bin, rk_pair, ak_pair;
    static std::vector<std::vector<double>> wnn;
    static std::vector<std::vector<double>> bin;
    static std::vector<std::vector<double>> rk_pair;
    static std::vector<std::vector<double>> ak_pair;

    static std::vector<int> jsort;

    //! Output for regression optimization
    static double efit_0;
    // double efit_rhorho[2], efit_rhorhoD[2], efit_rhotau[2], efit_rhoDrho[2];
    static double efit_rhorho[2];
    static double efit_rhorhoD[2];
    static double efit_rhotau[2];
    static double efit_rhoDrho[2];

    // double efit_rhonablaJ[2], efit_JJ[2], efitV0[2], dfitV0[2], efV_0[2];
    static double efit_rhonablaJ[2];
    static double efit_JJ[2];
    static double efitV0[2];
    static double dfitV0[2];
    static double efV_0[2];
    //! serial output (1:on/0:off)
    static int IDEBUG;
    static bool DO_FITT;
    //! For loop over used particle types. For normal nuclei min=1, max=2. For n droplets min=max=1.
    // int itmin, itmax, irestart;
    static int itmin;
    static int itmax;
    static int irestart;
    //! Global numbers
    // double zero = 0.00, half = 0.50, one = 1.00, two = 2.0, three = 3.0, four = 4.0, five = 5.0, six = 6.0, seven = 7.0, eight = 8.0, nine = 9.0, ten = 10.0;
    static double zero;
    static double half;
    static double one;
    static double two;
    static double three;
    static double four;
    static double five;
    static double six;
    static double seven;
    static double eight;
    static double nine;
    static double ten;
    //! Whole global numbers pp#
    // double pp12 = 12.0, pp16 = 16.0, pp15 = 15.0, pp20 = 20.0, pp24 = 24.0, pp27 = 27.0, pp32 = 32.0, pp64 = 64.0, pp40 = 40.0;
    static double pp12;
    static double pp16;
    static double pp15;
    static double pp20;
    static double pp24;
    static double pp27;
    static double pp32;
    static double pp64;
    static double pp40;
    //! Fractional global numbers p#
    // double p12 = one / two, p13 = one / three, p14 = 0.250, p23 = two / three, p43 = four / three, p32 = 1.50, p34 = three / four, p53 = five / three, p18 = one / eight, p38 = three / eight, p59 = five / nine, p52 = 2.50, p54 = five / four, p74 = seven / four;
    // double p12 = one / two, p13 = one / three, p14 = 0.250, p23 = two / three, p43 = four / three, p32 = 1.50, p34 = three / four, p53 = five / three, p18 = one / eight, p38 = three / eight, p59 = five / nine, p52 = 2.50, p54 = five / four, p74 = seven / four;
    static double p12;
    static double p13;
    static double p14;
    static double p23;
    static double p43;
    static double p32;
    static double p34;
    static double p53;
    static double p18;
    static double p38;
    static double p59;
    static double p52;
    static double p54;
    static double p74;
    //! Frequent Constants
    // double PI, ffdef3, ffdef4, ffdef5, ffdef6, ffdef7;
    static double PI;
    static double ffdef3;
    static double ffdef4;
    static double ffdef5;
    static double ffdef6;
    static double ffdef7;
    //! Single constants
    // double bet, beta0, q, bp, bpp, bz, hom, hb0, b0, etot, coex, t0s, t0a, drs, dra, ts, ta, t3alp, t3al0, t3alm, alp, alm, wla0, wla1, cex, cdef, cqad, ty20, vin, rin, ain, qin, pwi, si, siold, epsi, xmix, xmix0, xmax, alst, clst, sklst, alphi, amas, skass, varmas, v0ws, akv, hqc, amu, r0, r00, r02, r04, decay, rmm3, amm3, bmm3, cmm3, chargee2, EBASECUT;
    static double bet;
    static double beta0;
    static double q;
    static double bp;
    static double bpp;
    static double bz;
    static double hom;
    static double hb0;
    static double b0;
    static double etot;
    static double coex;
    static double t0s;
    static double t0a;
    static double drs;
    static double dra;
    static double ts;
    static double ta;
    static double t3alp;
    static double t3al0;
    static double t3alm;
    static double alp;
    static double alm;
    static double wla0;
    static double wla1;
    static double cex;
    static double cdef;
    static double cqad;
    static double ty20;
    static double vin;
    static double rin;
    static double ain;
    static double qin;
    static double pwi;
    static double si;
    static double siold;
    static double epsi;
    static double xmix;
    static double xmix0;
    static double xmax;
    static double alst;
    static double clst;
    static double sklst;
    static double alphi;
    static double amas;
    static double skass;
    static double varmas;
    static double v0ws;
    static double akv;
    static double hqc;
    static double amu;
    static double r0;
    static double r00;
    static double r02;
    static double r04;
    static double decay;
    static double rmm3;
    static double amm3;
    static double bmm3;
    static double cmm3;
    static double chargee2;
    static double EBASECUT;

    static double rho_c;
    // int lfile, lin, lout, lwin, lwou, lplo, lwel, lres, icstr, icou, ncut, iLST1, iLST, maxi, iiter, inin, nzm, nrm, icacou, iqrpa, icacoupj, icahartree, nlm, nb, nt, n00, itass, kindhfb, iappend, iError_in_HO, iError_in_THO, ierest;
    static int lfile;
    static int lin;
    static int lout;
    static int lwin;
    static int lwou;
    static int lplo;
    static int lwel;
    static int lres;
    static int icstr;
    static int icou;
    static int ncut;
    static int iLST1;
    static int iLST;
    static int maxi;
    static int iiter;
    static int inin;
    static int nzm;
    static int nrm;
    static int icacou;
    static int iqrpa;
    static int icacoupj;
    static int icahartree;
    static int nlm;
    static int nb;
    static int nt;
    static int n00;
    static int itass;
    static int kindhfb;
    static int iappend;
    static int iError_in_HO;
    static int iError_in_THO;
    static int ierest;

    // int n00max = 50, esu;
    static int n00max;
    static int esu;
    //! Results
    // const int ieresu = 50, ieresl = 20, ieresj = 50, ieresbl = 6;
    static const int ieresu;
    static const int ieresl;
    static const int ieresj;
    static const int ieresbl;

    // const int ieres = ieresu + ieresl + ieresj + ieresbl;
    static const int ieres;
    // double eres[ieres];
    static double eres[96];
    static char ereslbl[2][14];
    static char nucname[2];
    // double eresu[ieresu], eresl[ieresl], eresbl[ieresbl], eresj[ieresj];
    // double eresu[50], eresl[20], eresbl[6], eresj[50];
    static double eresu[50];
    static double eresl[20];
    static double eresbl[6];
    static double eresj[50];

    // char hlabels[ieres + 3][13];
    static char hlabels[99][13];
    //! Common small arrays
    // double alast[2], ala[2], ala1[2], tz[2], ass[2], drhoi[2], del[2], vso[2], r0v[2], av[2], rso[2], aso[2], Sumnz[2], Dispersion[2], v2min[2], v2minv[2], rms[3], ept[3], q2[3], Dnfactor[3], varmasnz[2], pjmassnz[2];
    static double alast[2];
    static double ala[2];
    static double ala1[2];
    static double tz[2];
    static double ass[2];
    static double drhoi[2];
    static double del[2];
    static double vso[2];
    static double r0v[2];
    static double av[2];
    static double rso[2];
    static double aso[2];
    static double Sumnz[2];
    static double Dispersion[2];
    static double v2min[2];
    static double v2minv[2];
    static double rms[3];
    static double ept[3];
    static double q2[3];
    static double Dnfactor[3];
    static double varmasnz[2];
    static double pjmassnz[2];

    // int npr[3], inz[2], ldel[2], nk[2], itbl[2], kbl[2], tpar[2], ipbl[2], nbl[2], ibbl[2], klmax[2], inner[2], iasswrong[3], lcc; //! remove
    static int npr[3];
    static int inz[2];
    static int ldel[2];
    static int nk[2];
    static int itbl[2];
    static int kbl[2];
    static int tpar[2];
    static int ipbl[2];
    static int nbl[2];
    static int ibbl[2];
    static int klmax[2];
    static int inner[2];
    static int iasswrong[3];
    static int lcc; //! remove

    //! Lipkin-Nogami
    // double ala2[2], etr[3], ssln[3][2], Geff[2];
    static double ala2[2];
    static double etr[3];
    static double ssln[3][2];
    static double Geff[2];
    //! Blocking
    // double pwiblo = 2.0, eqpmin[2] = {0.0};
    static double pwiblo;
    static double eqpmin[2];

    static const int bloall = 200;
    // int bloblo[201][2], blo123[201][2] = {0}, blok1k2[201][2] = {0};
    static int bloblo[201][2];
    static int blo123[201][2];
    static int blok1k2[201][2];

    static double bloqpdif[201][2];

    // int iparenti[2], keyblo[3], nkblo_INI[2][5], nkblo[2][5] = {0};
    static int iparenti[2];
    static int keyblo[3];
    // static int nkblo_INI[2][5];
    static int nkblo[2][5];

    // int blocross[2], blomax[2], blo123d[2], blok1k2d[2], blocanon[2];
    static int blocross[2];
    static int blomax[2];
    static int blo123d[2];
    static int blok1k2d[2];
    static int blocanon[2];

    //! manualBlocking
    static int manualBlocking;
    //! bool and character variables
    // char tq, tp[2], tl[21], tis[2];
    static char tq;
    static char tp[2];
    static char tl[21];
    static char tis[2];

    static char skyrme[30];
    static std::string tit[2];
    static char protn[2][8];
    //! Allocatable arrays
    static std::vector<std::vector<char>> tb;
    static std::vector<std::vector<char>> txb;
    // std::vector<std::vector<double>> rk, ak, hh0, de0, qh, qh1, ek, dk, vk, vk1, uk, hfb1, vkmax;
    static std::vector<std::vector<double>> rk;
    static std::vector<std::vector<double>> ak;
    static std::vector<std::vector<double>> hh0;
    static std::vector<std::vector<double>> de0;
    static std::vector<std::vector<double>> qh;
    static std::vector<std::vector<double>> qh1;
    static std::vector<std::vector<double>> ek;
    static std::vector<std::vector<double>> dk;
    static std::vector<std::vector<double>> vk;
    static std::vector<std::vector<double>> vk1;
    static std::vector<std::vector<double>> uk;
    static std::vector<std::vector<double>> hfb1;
    static std::vector<std::vector<double>> vkmax;

    // std::vector<std::vector<std::vector<double>>> ddc, ddc1, ql, ql1;
    static std::vector<std::vector<std::vector<double>>> ddc;
    static std::vector<std::vector<std::vector<double>>> ddc1;
    static std::vector<std::vector<std::vector<double>>> ql;
    static std::vector<std::vector<std::vector<double>>> ql1;

    // std::vector<double> fdsx, fdsy, fdsy1, fdsy2, fdsy3, fspb0, fspc0, fspd0, fspb1, fspc1, fspd1, fspb2, fspc2, fspd2, fspb3, fspc3, fspd3, fak, fi, sq, sqi, wf, wfi;
    static std::vector<double> fdsx;
    static std::vector<double> fdsy;
    static std::vector<double> fdsy1;
    static std::vector<double> fdsy2;
    static std::vector<double> fdsy3;
    static std::vector<double> fspb0;
    static std::vector<double> fspc0;
    static std::vector<double> fspd0;
    static std::vector<double> fspb1;
    static std::vector<double> fspc1;
    static std::vector<double> fspd1;
    static std::vector<double> fspb2;
    static std::vector<double> fspc2;
    static std::vector<double> fspd2;
    static std::vector<double> fspb3;
    static std::vector<double> fspc3;
    static std::vector<double> fspd3;
    static std::vector<double> fak;
    static std::vector<double> fi;
    static std::vector<double> sq;
    static std::vector<double> sqi;
    static std::vector<double> wf;
    static std::vector<double> wfi;

    static std::vector<std::vector<double>> rkass;
    // std::vector<int> id, ia, ikb, ipb, nz, nr, nl, ns, npar, iv;
    static std::vector<int> id;
    static std::vector<int> ia;
    static std::vector<int> ikb;
    static std::vector<int> ipb;
    static std::vector<int> nz;
    static std::vector<int> nr;
    static std::vector<int> nl;
    static std::vector<int> ns;
    static std::vector<int> npar;
    static std::vector<int> iv;

    // std::vector<std::vector<int>> ka, kd, numax, lcanon;
    static std::vector<std::vector<int>> ka;
    static std::vector<std::vector<int>> kd;
    static std::vector<std::vector<int>> numax;
    static std::vector<std::vector<int>> lcanon;

    // std::vector<double> AN, ANK, PFIU, PFID;
    static std::vector<double> AN;
    static std::vector<double> ANK;
    static std::vector<double> PFIU;
    static std::vector<double> PFID;

    // std::vector<double> FIU, FID, FIUR, FIDR;
    static std::vector<double> FIU;
    static std::vector<double> FID;
    static std::vector<double> FIUR;
    static std::vector<double> FIDR;

    // std::vector<double> FIUD2N, FIDD2N, FIUZ, FIDZ;
    static std::vector<double> FIUD2N;
    static std::vector<double> FIDD2N;
    static std::vector<double> FIUZ;
    static std::vector<double> FIDZ;

    //! optimization arrays
    // std::vector<std::vector<double>> QHLA_opt, FI1R_opt, FI1Z_opt, FI2D_opt;
    static std::vector<std::vector<double>> QHLA_opt;
    static std::vector<std::vector<double>> FI1R_opt;
    static std::vector<std::vector<double>> FI1Z_opt;
    static std::vector<std::vector<double>> FI2D_opt;

    static std::vector<double> y_opt;
    //! Arrays depending on mesh points
    // int ngh, ngl, nleg, nghl, nbx, ntx, nzx, nrx, nlx, ndx, ndx2, ndxs, nqx;
    static int ngh;
    static int ngl;
    static int nleg;
    static int nghl;
    static int nbx;
    static int ntx;
    static int nzx;
    static int nrx;
    static int nlx;
    static int ndx;
    static int ndx2;
    static int ndxs;
    static int nqx;

    // int nhfbqx, nb2x, nhfbx, nkx, nzrlx, iqqmax;
    static int nhfbqx;
    static int nb2x;
    static int nhfbx;
    static int nkx;
    static int nzrlx;
    static int iqqmax;

    // std::vector<double> xh, wh, xl, sxl, wl, xleg, wleg;
    static std::vector<double> xh;
    static std::vector<double> wh;
    static std::vector<double> xl;
    static std::vector<double> sxl;
    static std::vector<double> wl;
    static std::vector<double> xleg;
    static std::vector<double> wleg;

    // std::vector<double> vhbn, vn, vrn, vzn, vdn, vsn, dvn;
    static std::vector<double> vhbn;
    static std::vector<double> vn;
    static std::vector<double> vrn;
    static std::vector<double> vzn;
    static std::vector<double> vdn;
    static std::vector<double> vsn;
    static std::vector<double> dvn;

    // std::vector<double> vhbp, vp, vrp, vzp, vdp, vsp, dvp;
    static std::vector<double> vhbp;
    static std::vector<double> vp;
    static std::vector<double> vrp;
    static std::vector<double> vzp;
    static std::vector<double> vdp;
    static std::vector<double> vsp;
    static std::vector<double> dvp;

    static std::vector<std::vector<double>> vc;
    // std::vector<double> vSZFIn, vSFIZn, vSRFIn, vSFIRn;
    static std::vector<double> vSZFIn;
    static std::vector<double> vSFIZn;
    static std::vector<double> vSRFIn;
    static std::vector<double> vSFIRn;

    // std::vector<double> vSZFIp, vSFIZp, vSRFIp, vSFIRp;
    static std::vector<double> vSZFIp;
    static std::vector<double> vSFIZp;
    static std::vector<double> vSRFIp;
    static std::vector<double> vSFIRp;

    // std::vector<std::vector<double>> aka, ro, tau, dro, dj, SZFI, SFIZ, SRFI, SFIR, NABLAR, NABLAZ;
    static std::vector<std::vector<double>> aka;
    static std::vector<std::vector<double>> ro;
    static std::vector<std::vector<double>> tau;
    static std::vector<std::vector<double>> dro;
    static std::vector<std::vector<double>> dj;
    static std::vector<std::vector<double>> SZFI;
    static std::vector<std::vector<double>> SFIZ;
    static std::vector<std::vector<double>> SRFI;
    static std::vector<std::vector<double>> SFIR;
    static std::vector<std::vector<double>> NABLAR;
    static std::vector<std::vector<double>> NABLAZ;

    // std::vector<double> fl, fli, fh, fd, fp1, fp2, fp3, fp4, fp5, fp6, fs1, fs2, fs3, fs4, fs5, fs6, wdcor, wdcori, cou;
    static std::vector<double> fl;
    static std::vector<double> fli;
    static std::vector<double> fh;
    static std::vector<double> fd;
    static std::vector<double> fp1;
    static std::vector<double> fp2;
    static std::vector<double> fp3;
    static std::vector<double> fp4;
    static std::vector<double> fp5;
    static std::vector<double> fp6;
    static std::vector<double> fs1;
    static std::vector<double> fs2;
    static std::vector<double> fs3;
    static std::vector<double> fs4;
    static std::vector<double> fs5;
    static std::vector<double> fs6;
    static std::vector<double> wdcor;
    static std::vector<double> wdcori;
    static std::vector<double> cou;

    // std::vector<std::vector<double>> vDHartree, vhart00, vhart01, vhart11;
    static std::vector<std::vector<double>> vDHartree;
    static std::vector<std::vector<double>> vhart00;
    static std::vector<std::vector<double>> vhart01;
    static std::vector<std::vector<double>> vhart11;

    //! PAV Projection
    // int keypj, ilpj, ilpj2, ilnqx, ilnghl;
    static int keypj;
    static int ilpj;
    static int ilpj2;
    static int ilnqx;
    static int ilnghl;

    // double rehfbcan, ehfb, retotpj, depnp, iproj, npr1pj, npr2pj;
    static double rehfbcan;
    static double ehfb;
    static double retotpj;
    static double depnp;
    static double iproj;
    static double npr1pj;
    static double npr2pj;

    static std::complex<double> onei;
    // std::vector<std::complex<double>> phypj, sinphy, exp1iphy, exp1iphym, exp2iphy, exp2iphym;
    static std::vector<std::complex<double>> phypj;
    static std::vector<std::complex<double>> sinphy;
    static std::vector<std::complex<double>> exp1iphy;
    static std::vector<std::complex<double>> exp1iphym;
    static std::vector<std::complex<double>> exp2iphy;
    static std::vector<std::complex<double>> exp2iphym;

    // std::vector<std::vector<std::complex<double>>> coupj, pjk, epj;
    static std::vector<std::vector<std::complex<double>>> coupj;
    static std::vector<std::vector<std::complex<double>>> pjk;
    static std::vector<std::vector<std::complex<double>>> epj;

    // std::vector<std::vector<std::vector<std::complex<double>>>> ropj, taupj, dropj, djpj, akapj, SZFIpj, SFIZpj, SRFIpj, SFIRpj, ddepj, cpj, ypj, rpj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> ropj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> taupj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> dropj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> djpj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> akapj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> SZFIpj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> SFIZpj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> SRFIpj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> SFIRpj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> ddepj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> cpj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> ypj;
    static std::vector<std::vector<std::vector<std::complex<double>>>> rpj;

    // double polem[2], polep[2];
    static double polem[2];
    static double polep[2];

    //! CMC
    static int ICMinput;
    // double ECMHFB[3], ECMPAV[3];
    static double ECMHFB[3];
    static double ECMPAV[3];

    //! CRC
    static int ICRinput;
    // double DEROT[3], SQUJ[3], CRAN[3], ERIGHFB[3];
    static double DEROT[3];
    static double SQUJ[3];
    static double CRAN[3];
    static double ERIGHFB[3];

    //! hfbdiagonal
    // std::vector<double> erhfb, drhfb, erhfb1, drhfb1, evvk, evvkcan, zhfb;
    static std::vector<double> erhfb;
    static std::vector<double> drhfb;
    static std::vector<double> erhfb1;
    static std::vector<double> drhfb1;
    static std::vector<double> evvk;
    static std::vector<double> evvkcan;
    static std::vector<double> zhfb;

    // std::vector<std::vector<double>> hfb, hfbcan;
    // static std::vector<std::vector<double>> hfb;
    static std::vector<double> hfb;
    static std::vector<std::vector<double>> hfbcan;

    //! Broyden
    static char bbroyden;
    static int nbroyden;
    static double alphamix;
    // int nhhdim, nhhdim2, nhhdim3, nhhdim4, ialwork, ilwork;
    static int nhhdim;
    static int nhhdim2;
    static int nhhdim3;
    static int nhhdim4;
    static int ialwork;
    static int ilwork;

    // std::vector<double> brout, brin;
    static std::vector<double> brout;
    static std::vector<double> brin;

    static std::vector<double> alwork;
    static std::vector<int> lwork;
    //! new keys
    // bool Parity, Parity_INI;
    static bool Parity;
    // static bool Parity_INI;

    static bool Print_Screen;
    // bool Add_Pairing, Print_HFBTHO_Namelist;
    static bool Add_Pairing;
    static bool Print_HFBTHO_Namelist;

    // int MAX_ITER_INI, keypj_INI, iproj_INI, npr1pj_INI, npr2pj_INI;
    // static int MAX_ITER_INI;
    // static int keypj_INI;
    // static int iproj_INI;
    // static int npr1pj_INI;
    // static int npr2pj_INI;

    //! Eqp U,V
    // int nuv, nqp;
    static int nuv;
    static int nqp;

    // std::vector<double> RVqpN, RVqpP, RUqpN, RUqpP, REqpN, REqpP;
    static std::vector<double> RVqpN;
    static std::vector<double> RVqpP;
    static std::vector<double> RUqpN;
    static std::vector<double> RUqpP;
    static std::vector<double> REqpN;
    static std::vector<double> REqpP;

    // std::vector<int> KpwiN, KpwiP, KqpN, KqpP;
    static std::vector<int> KpwiN;
    static std::vector<int> KpwiP;
    static std::vector<int> KqpN;
    static std::vector<int> KqpP;

    //! error indicator
    static int ierror_flag;
    static std::string ierror_info[11];
    //! namelist
    //   Namelist /HFBTHO_NAMELIST/ MAX_ITER_INI,epsi_INI,Add_Pairing_INI &
    //        ,icou_INI,iLST_INI,keypj_INI,iproj_INI,npr1pj_INI,npr2pj_INI &
    //        ,DO_FITT_INI,IDEBUG_INI,Parity_INI,Print_HFBTHO_Namelist_INI

    //!
    //! mpi setup
    static int iam_mpi; //!! mpi id number of this process
    // int num_processors_mpi, num_nodes_mpi, icom_err_mpi; //!! for mpi
    static int num_processors_mpi;
    static int num_nodes_mpi;
    static int icom_err_mpi; //!! for mpi
};

#endif
