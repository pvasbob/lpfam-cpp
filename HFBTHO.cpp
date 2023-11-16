#include <fstream>
#include <iostream>

#include "HFBTHO.h"

// define static variables declared in HFBTHO.h
char HFBTHO::Version[6] = {"133"};
char HFBTHO::pairing_Version[6] = {"9"};

// int NNN0, mv, ivpair = 0;
int HFBTHO::NNN0;
int HFBTHO::mv;
int HFBTHO::ivpair = 0;

// double gl0, gal, gl[2], vfac[2], e_pair[2], d_pair[2];
double HFBTHO::gl0;
double HFBTHO::gal;
double HFBTHO::gl[2];
double HFBTHO::vfac[2];
double HFBTHO::e_pair[2];
double HFBTHO::d_pair[2];

// std::vector<std::vector<double>> wnn, bin, rk_pair, ak_pair;
std::vector<std::vector<double>> HFBTHO::wnn;
std::vector<std::vector<double>> HFBTHO::bin;
std::vector<std::vector<double>> HFBTHO::rk_pair;
std::vector<std::vector<double>> HFBTHO::ak_pair;

std::vector<int> HFBTHO::jsort;

// int n00_INI, iLST_INI, inin_INI, icou_INI;
int HFBTHO::n00_INI;
int HFBTHO::iLST_INI;
int HFBTHO::inin_INI;
int HFBTHO::icou_INI;

// int npr_INI[3], kindhfb_INI;
int HFBTHO::npr_INI[3];
int HFBTHO::kindhfb_INI;

// int keyblo1_INI, keyblo2_INI, IDEBUG_INI;
int HFBTHO::keyblo1_INI = 0;
int HFBTHO::keyblo2_INI = 0;
int HFBTHO::IDEBUG_INI;

// double b0_INI, q_INI, cdef_INI, cqad_INI;
double HFBTHO::b0_INI;
double HFBTHO::q_INI;
double HFBTHO::cdef_INI;
double HFBTHO::cqad_INI;

char HFBTHO::skyrme_INI[30];
double HFBTHO::epsi_INI;
// bool Add_Pairing_INI, Print_HFBTHO_Namelist_INI, DO_FITT_INI;
bool HFBTHO::Add_Pairing_INI;
bool HFBTHO::Print_HFBTHO_Namelist_INI;
bool HFBTHO::DO_FITT_INI;

//! Output for regression optimization
double HFBTHO::efit_0;
// double efit_rhorho[2], efit_rhorhoD[2], efit_rhotau[2], efit_rhoDrho[2];
double HFBTHO::efit_rhorho[2];
double HFBTHO::efit_rhorhoD[2];
double HFBTHO::efit_rhotau[2];
double HFBTHO::efit_rhoDrho[2];

// double efit_rhonablaJ[2], efit_JJ[2], efitV0[2], dfitV0[2], efV_0[2];
double HFBTHO::efit_rhonablaJ[2];
double HFBTHO::efit_JJ[2];
double HFBTHO::efitV0[2];
double HFBTHO::dfitV0[2];
double HFBTHO::efV_0[2];
//! serial output (1:on/0:off)
int HFBTHO::IDEBUG;
bool HFBTHO::DO_FITT;
//! For loop over used particle types. For normal nuclei min=1, max=2. For n droplets min=max=1.
// int itmin, itmax, irestart;
int HFBTHO::itmin;
int HFBTHO::itmax;
int HFBTHO::irestart;
//! Global numbers
// double zero = 0.00, half = 0.50, one = 1.00, two = 2.0, three = 3.0, four = 4.0, five = 5.0, six = 6.0, seven = 7.0, eight = 8.0, nine = 9.0, ten = 10.0;
double HFBTHO::zero = 0.00;
double HFBTHO::half = 0.50;
double HFBTHO::one = 1.00;
double HFBTHO::two = 2.0;
double HFBTHO::three = 3.0;
double HFBTHO::four = 4.0;
double HFBTHO::five = 5.0;
double HFBTHO::six = 6.0;
double HFBTHO::seven = 7.0;
double HFBTHO::eight = 8.0;
double HFBTHO::nine = 9.0;
double HFBTHO::ten = 10.0;
//! Whole global numbers pp#
// double pp12 = 12.0, pp16 = 16.0, pp15 = 15.0, pp20 = 20.0, pp24 = 24.0, pp27 = 27.0, pp32 = 32.0, pp64 = 64.0, pp40 = 40.0;
double HFBTHO::pp12 = 12.0;
double HFBTHO::pp16 = 16.0;
double HFBTHO::pp15 = 15.0;
double HFBTHO::pp20 = 20.0;
double HFBTHO::pp24 = 24.0;
double HFBTHO::pp27 = 27.0;
double HFBTHO::pp32 = 32.0;
double HFBTHO::pp64 = 64.0;
double HFBTHO::pp40 = 40.0;
//! Fractional global numbers p#
// double p12 = one / two, p13 = one / three, p14 = 0.250, p23 = two / three, p43 = four / three, p32 = 1.50, p34 = three / four, p53 = five / three, p18 = one / eight, p38 = three / eight, p59 = five / nine, p52 = 2.50, p54 = five / four, p74 = seven / four;
// double p12 = one / two, p13 = one / three, p14 = 0.250, p23 = two / three, p43 = four / three, p32 = 1.50, p34 = three / four, p53 = five / three, p18 = one / eight, p38 = three / eight, p59 = five / nine, p52 = 2.50, p54 = five / four, p74 = seven / four;
double HFBTHO::p12 = one / two;
double HFBTHO::p13 = one / three;
double HFBTHO::p14 = 0.250;
double HFBTHO::p23 = two / three;
double HFBTHO::p43 = four / three;
double HFBTHO::p32 = 1.50;
double HFBTHO::p34 = three / four;
double HFBTHO::p53 = five / three;
double HFBTHO::p18 = one / eight;
double HFBTHO::p38 = three / eight;
double HFBTHO::p59 = five / nine;
double HFBTHO::p52 = 2.50;
double HFBTHO::p54 = five / four;
double HFBTHO::p74 = seven / four;
//! Frequent Constants
// double PI, ffdef3, ffdef4, ffdef5, ffdef6, ffdef7;
double HFBTHO::PI = 3.14159265357;
double HFBTHO::ffdef3;
double HFBTHO::ffdef4;
double HFBTHO::ffdef5;
double HFBTHO::ffdef6;
double HFBTHO::ffdef7;
//! Single constants
// double bet, beta0, q, bp, bpp, bz, hom, hb0, b0, etot, coex, t0s, t0a, drs, dra, ts, ta, t3alp, t3al0, t3alm, alp, alm, wla0, wla1, cex, cdef, cqad, ty20, vin, rin, ain, qin, pwi, si, siold, epsi, xmix, xmix0, xmax, alst, clst, sklst, alphi, amas, skass, varmas, v0ws, akv, hqc, amu, r0, r00, r02, r04, decay, rmm3, amm3, bmm3, cmm3, chargee2, EBASECUT;
double HFBTHO::bet;
double HFBTHO::beta0;
double HFBTHO::q;
double HFBTHO::bp;
double HFBTHO::bpp;
double HFBTHO::bz;
double HFBTHO::hom;
double HFBTHO::hb0;
double HFBTHO::b0;
double HFBTHO::etot;
double HFBTHO::coex;
double HFBTHO::t0s;
double HFBTHO::t0a;
double HFBTHO::drs;
double HFBTHO::dra;
double HFBTHO::ts;
double HFBTHO::ta;
double HFBTHO::t3alp;
double HFBTHO::t3al0;
double HFBTHO::t3alm;
double HFBTHO::alp;
double HFBTHO::alm;
double HFBTHO::wla0;
double HFBTHO::wla1;
double HFBTHO::cex;
double HFBTHO::cdef;
double HFBTHO::cqad;
double HFBTHO::ty20;
double HFBTHO::vin;
double HFBTHO::rin;
double HFBTHO::ain;
double HFBTHO::qin;
double HFBTHO::pwi;
double HFBTHO::si;
double HFBTHO::siold;
double HFBTHO::epsi;
double HFBTHO::xmix;
double HFBTHO::xmix0;
double HFBTHO::xmax;
double HFBTHO::alst;
double HFBTHO::clst;
double HFBTHO::sklst;
double HFBTHO::alphi;
double HFBTHO::amas;
double HFBTHO::skass;
double HFBTHO::varmas;
double HFBTHO::v0ws;
double HFBTHO::akv;
double HFBTHO::hqc;
double HFBTHO::amu;
double HFBTHO::r0;
double HFBTHO::r00;
double HFBTHO::r02;
double HFBTHO::r04;
double HFBTHO::decay;
double HFBTHO::rmm3;
double HFBTHO::amm3;
double HFBTHO::bmm3;
double HFBTHO::cmm3;
double HFBTHO::chargee2;
double HFBTHO::EBASECUT;

double HFBTHO::rho_c;
// int lfile, lin, lout, lwin, lwou, lplo, lwel, lres, icstr, icou, ncut, iLST1, iLST, maxi, iiter, inin, nzm, nrm, icacou, iqrpa, icacoupj, icahartree, nlm, nb, nt, n00, itass, kindhfb, iappend, iError_in_HO, iError_in_THO, ierest;
int HFBTHO::lfile;
int HFBTHO::lin;
int HFBTHO::lout;
int HFBTHO::lwin;
int HFBTHO::lwou;
int HFBTHO::lplo;
int HFBTHO::lwel;
int HFBTHO::lres;
int HFBTHO::icstr;
int HFBTHO::icou;
int HFBTHO::ncut;
int HFBTHO::iLST1;
int HFBTHO::iLST;
int HFBTHO::maxi;
int HFBTHO::iiter;
int HFBTHO::inin;
int HFBTHO::nzm;
int HFBTHO::nrm;
int HFBTHO::icacou;
int HFBTHO::iqrpa;
int HFBTHO::icacoupj;
int HFBTHO::icahartree;
int HFBTHO::nlm;
int HFBTHO::nb;
int HFBTHO::nt;
int HFBTHO::n00;
int HFBTHO::itass;
int HFBTHO::kindhfb;
int HFBTHO::iappend;
int HFBTHO::iError_in_HO;
int HFBTHO::iError_in_THO;
int HFBTHO::ierest;

// int n00max = 50, esu;
int HFBTHO::n00max = 50;
int HFBTHO::esu;
//! Results
// const int ieresu = 50, ieresl = 20, ieresj = 50, ieresbl = 6;
const int HFBTHO::ieresu = 50;
const int HFBTHO::ieresl = 20;
const int HFBTHO::ieresj = 50;
const int HFBTHO::ieresbl = 6;

// const int ieres = ieresu + ieresl + ieresj + ieresbl;
const int HFBTHO::ieres = 96;
// double eres[ieres];
double HFBTHO::eres[96];
char HFBTHO::ereslbl[2][14] = {" 00[00,00,00]"};
char HFBTHO::nucname[2];
// double eresu[ieresu], eresl[ieresl], eresbl[ieresbl], eresj[ieresj];
// double eresu[50], eresl[20], eresbl[6], eresj[50];
double HFBTHO::eresu[50] = {0};
double HFBTHO::eresl[20] = {0};
double HFBTHO::eresbl[6] = {0};
double HFBTHO::eresj[50] = {0};

// char hlabels[ieres + 3][13];
char HFBTHO::hlabels[99][13];
//! Common small arrays
// double alast[2], ala[2], ala1[2], tz[2], ass[2], drhoi[2], del[2], vso[2], r0v[2], av[2], rso[2], aso[2], Sumnz[2], Dispersion[2], v2min[2], v2minv[2], rms[3], ept[3], q2[3], Dnfactor[3], varmasnz[2], pjmassnz[2];
double HFBTHO::alast[2];
double HFBTHO::ala[2];
double HFBTHO::ala1[2];
double HFBTHO::tz[2];
double HFBTHO::ass[2];
double HFBTHO::drhoi[2];
double HFBTHO::del[2];
double HFBTHO::vso[2];
double HFBTHO::r0v[2];
double HFBTHO::av[2];
double HFBTHO::rso[2];
double HFBTHO::aso[2];
double HFBTHO::Sumnz[2];
double HFBTHO::Dispersion[2];
double HFBTHO::v2min[2];
double HFBTHO::v2minv[2];
double HFBTHO::rms[3];
double HFBTHO::ept[3];
double HFBTHO::q2[3];
double HFBTHO::Dnfactor[3];
double HFBTHO::varmasnz[2];
double HFBTHO::pjmassnz[2];

// int npr[3], inz[2], ldel[2], nk[2], itbl[2], kbl[2], tpar[2], ipbl[2], nbl[2], ibbl[2], klmax[2], inner[2], iasswrong[3], lcc; //! remove
int HFBTHO::npr[3];
int HFBTHO::inz[2];
int HFBTHO::ldel[2];
int HFBTHO::nk[2];
int HFBTHO::itbl[2];
int HFBTHO::kbl[2];
int HFBTHO::tpar[2];
int HFBTHO::ipbl[2];
int HFBTHO::nbl[2];
int HFBTHO::ibbl[2];
int HFBTHO::klmax[2];
int HFBTHO::inner[2];
int HFBTHO::iasswrong[3];
int HFBTHO::lcc; //! remove

//! Lipkin-Nogami
// double ala2[2], etr[3], ssln[3][2], Geff[2];
double HFBTHO::ala2[2];
double HFBTHO::etr[3];
double HFBTHO::ssln[3][2];
double HFBTHO::Geff[2];
//! Blocking
// double pwiblo = 2.0, eqpmin[2] = {0.0};
double HFBTHO::pwiblo = 2.0;
double HFBTHO::eqpmin[2] = {0.0};

const int bloall = 200;
// int bloblo[201][2], blo123[201][2] = {0}, blok1k2[201][2] = {0};
int HFBTHO::bloblo[201][2] = {0};
int HFBTHO::blo123[201][2] = {0};
int HFBTHO::blok1k2[201][2] = {0};

double bloqpdif[201][2] = {0};

// int iparenti[2], keyblo[3], nkblo_INI[2][5], nkblo[2][5] = {0};
int HFBTHO::iparenti[2];
int HFBTHO::keyblo[3] = {0};
int HFBTHO::nkblo_INI[2][5] = {0};
int HFBTHO::nkblo[2][5] = {0};

// int blocross[2], blomax[2], blo123d[2], blok1k2d[2], blocanon[2];
int HFBTHO::blocross[2] = {0};
int HFBTHO::blomax[2] = {0};
int HFBTHO::blo123d[2] = {0};
int HFBTHO::blok1k2d[2] = {0};
int HFBTHO::blocanon[2] = {0};

//! manualBlocking
int HFBTHO::manualBlocking = 0;
//! bool and character variables
// char tq, tp[2], tl[21], tis[2];
char HFBTHO::tq;
char HFBTHO::tp[2];
char HFBTHO::tl[21];
char HFBTHO::tis[2];

char HFBTHO::skyrme[30];
std::string HFBTHO::tit[2];
char HFBTHO::protn[2][8] = {"neutron", "proton"};
//! Allocatable arrays
std::vector<std::vector<char>> HFBTHO::tb;
std::vector<std::vector<char>> HFBTHO::txb;
// std::vector<std::vector<double>> rk, ak, hh0, de0, qh, qh1, ek, dk, vk, vk1, uk, hfb1, vkmax;
std::vector<std::vector<double>> HFBTHO::rk;
std::vector<std::vector<double>> HFBTHO::ak;
std::vector<std::vector<double>> HFBTHO::hh0;
std::vector<std::vector<double>> HFBTHO::de0;
std::vector<std::vector<double>> HFBTHO::qh;
std::vector<std::vector<double>> HFBTHO::qh1;
std::vector<std::vector<double>> HFBTHO::ek;
std::vector<std::vector<double>> HFBTHO::dk;
std::vector<std::vector<double>> HFBTHO::vk;
std::vector<std::vector<double>> HFBTHO::vk1;
std::vector<std::vector<double>> HFBTHO::uk;
std::vector<std::vector<double>> HFBTHO::hfb1;
std::vector<std::vector<double>> HFBTHO::vkmax;

// std::vector<std::vector<std::vector<double>>> ddc, ddc1, ql, ql1;
std::vector<std::vector<std::vector<double>>> HFBTHO::ddc;
std::vector<std::vector<std::vector<double>>> HFBTHO::ddc1;
std::vector<std::vector<std::vector<double>>> HFBTHO::ql;
std::vector<std::vector<std::vector<double>>> HFBTHO::ql1;

// std::vector<double> fdsx, fdsy, fdsy1, fdsy2, fdsy3, fspb0, fspc0, fspd0, fspb1, fspc1, fspd1, fspb2, fspc2, fspd2, fspb3, fspc3, fspd3, fak, fi, sq, sqi, wf, wfi;
std::vector<double> HFBTHO::fdsx;
std::vector<double> HFBTHO::fdsy;
std::vector<double> HFBTHO::fdsy1;
std::vector<double> HFBTHO::fdsy2;
std::vector<double> HFBTHO::fdsy3;
std::vector<double> HFBTHO::fspb0;
std::vector<double> HFBTHO::fspc0;
std::vector<double> HFBTHO::fspd0;
std::vector<double> HFBTHO::fspb1;
std::vector<double> HFBTHO::fspc1;
std::vector<double> HFBTHO::fspd1;
std::vector<double> HFBTHO::fspb2;
std::vector<double> HFBTHO::fspc2;
std::vector<double> HFBTHO::fspd2;
std::vector<double> HFBTHO::fspb3;
std::vector<double> HFBTHO::fspc3;
std::vector<double> HFBTHO::fspd3;
std::vector<double> HFBTHO::fak;
std::vector<double> HFBTHO::fi;
std::vector<double> HFBTHO::sq;
std::vector<double> HFBTHO::sqi;
std::vector<double> HFBTHO::wf;
std::vector<double> HFBTHO::wfi;

std::vector<std::vector<double>> HFBTHO::rkass;
// std::vector<int> id, ia, ikb, ipb, nz, nr, nl, ns, npar, iv;
std::vector<int> HFBTHO::id;
std::vector<int> HFBTHO::ia;
std::vector<int> HFBTHO::ikb;
std::vector<int> HFBTHO::ipb;
std::vector<int> HFBTHO::nz;
std::vector<int> HFBTHO::nr;
std::vector<int> HFBTHO::nl;
std::vector<int> HFBTHO::ns;
std::vector<int> HFBTHO::npar;
std::vector<int> HFBTHO::iv;

// std::vector<std::vector<int>> ka, kd, numax, lcanon;
std::vector<std::vector<int>> HFBTHO::ka;
std::vector<std::vector<int>> HFBTHO::kd;
std::vector<std::vector<int>> HFBTHO::numax;
std::vector<std::vector<int>> HFBTHO::lcanon;

// std::vector<double> AN, ANK, PFIU, PFID;
std::vector<double> HFBTHO::AN;
std::vector<double> HFBTHO::ANK;
std::vector<double> HFBTHO::PFIU;
std::vector<double> HFBTHO::PFID;

// std::vector<double> FIU, FID, FIUR, FIDR;
std::vector<double> HFBTHO::FIU;
std::vector<double> HFBTHO::FID;
std::vector<double> HFBTHO::FIUR;
std::vector<double> HFBTHO::FIDR;

// std::vector<double> FIUD2N, FIDD2N, FIUZ, FIDZ;
std::vector<double> HFBTHO::FIUD2N;
std::vector<double> HFBTHO::FIDD2N;
std::vector<double> HFBTHO::FIUZ;
std::vector<double> HFBTHO::FIDZ;

//! optimization arrays
// std::vector<std::vector<double>> QHLA_opt, FI1R_opt, FI1Z_opt, FI2D_opt;
std::vector<std::vector<double>> HFBTHO::QHLA_opt;
std::vector<std::vector<double>> HFBTHO::FI1R_opt;
std::vector<std::vector<double>> HFBTHO::FI1Z_opt;
std::vector<std::vector<double>> HFBTHO::FI2D_opt;

std::vector<double> HFBTHO::y_opt;
//! Arrays depending on mesh points
// int ngh, ngl, nleg, nghl, nbx, ntx, nzx, nrx, nlx, ndx, ndx2, ndxs, nqx;
int HFBTHO::ngh;
int HFBTHO::ngl;
int HFBTHO::nleg;
int HFBTHO::nghl;
int HFBTHO::nbx;
int HFBTHO::ntx;
int HFBTHO::nzx;
int HFBTHO::nrx;
int HFBTHO::nlx;
int HFBTHO::ndx;
int HFBTHO::ndx2;
int HFBTHO::ndxs;
int HFBTHO::nqx;

// int nhfbqx, nb2x, nhfbx, nkx, nzrlx, iqqmax;
int HFBTHO::nhfbqx;
int HFBTHO::nb2x;
int HFBTHO::nhfbx;
int HFBTHO::nkx;
int HFBTHO::nzrlx;
int HFBTHO::iqqmax;

// std::vector<double> xh, wh, xl, sxl, wl, xleg, wleg;
std::vector<double> HFBTHO::xh;
std::vector<double> HFBTHO::wh;
std::vector<double> HFBTHO::xl;
std::vector<double> HFBTHO::sxl;
std::vector<double> HFBTHO::wl;
std::vector<double> HFBTHO::xleg;
std::vector<double> HFBTHO::wleg;

// std::vector<double> vhbn, vn, vrn, vzn, vdn, vsn, dvn;
std::vector<double> HFBTHO::vhbn;
std::vector<double> HFBTHO::vn;
std::vector<double> HFBTHO::vrn;
std::vector<double> HFBTHO::vzn;
std::vector<double> HFBTHO::vdn;
std::vector<double> HFBTHO::vsn;
std::vector<double> HFBTHO::dvn;

// std::vector<double> vhbp, vp, vrp, vzp, vdp, vsp, dvp;
std::vector<double> HFBTHO::vhbp;
std::vector<double> HFBTHO::vp;
std::vector<double> HFBTHO::vrp;
std::vector<double> HFBTHO::vzp;
std::vector<double> HFBTHO::vdp;
std::vector<double> HFBTHO::vsp;
std::vector<double> HFBTHO::dvp;

std::vector<std::vector<double>> HFBTHO::vc;
// std::vector<double> vSZFIn, vSFIZn, vSRFIn, vSFIRn;
std::vector<double> HFBTHO::vSZFIn;
std::vector<double> HFBTHO::vSFIZn;
std::vector<double> HFBTHO::vSRFIn;
std::vector<double> HFBTHO::vSFIRn;

// std::vector<double> vSZFIp, vSFIZp, vSRFIp, vSFIRp;
std::vector<double> HFBTHO::vSZFIp;
std::vector<double> HFBTHO::vSFIZp;
std::vector<double> HFBTHO::vSRFIp;
std::vector<double> HFBTHO::vSFIRp;

// std::vector<std::vector<double>> aka, ro, tau, dro, dj, SZFI, SFIZ, SRFI, SFIR, NABLAR, NABLAZ;
std::vector<std::vector<double>> HFBTHO::aka;
std::vector<std::vector<double>> HFBTHO::ro;
std::vector<std::vector<double>> HFBTHO::tau;
std::vector<std::vector<double>> HFBTHO::dro;
std::vector<std::vector<double>> HFBTHO::dj;
std::vector<std::vector<double>> HFBTHO::SZFI;
std::vector<std::vector<double>> HFBTHO::SFIZ;
std::vector<std::vector<double>> HFBTHO::SRFI;
std::vector<std::vector<double>> HFBTHO::SFIR;
std::vector<std::vector<double>> HFBTHO::NABLAR;
std::vector<std::vector<double>> HFBTHO::NABLAZ;

// std::vector<double> fl, fli, fh, fd, fp1, fp2, fp3, fp4, fp5, fp6, fs1, fs2, fs3, fs4, fs5, fs6, wdcor, wdcori, cou;
std::vector<double> HFBTHO::fl;
std::vector<double> HFBTHO::fli;
std::vector<double> HFBTHO::fh;
std::vector<double> HFBTHO::fd;
std::vector<double> HFBTHO::fp1;
std::vector<double> HFBTHO::fp2;
std::vector<double> HFBTHO::fp3;
std::vector<double> HFBTHO::fp4;
std::vector<double> HFBTHO::fp5;
std::vector<double> HFBTHO::fp6;
std::vector<double> HFBTHO::fs1;
std::vector<double> HFBTHO::fs2;
std::vector<double> HFBTHO::fs3;
std::vector<double> HFBTHO::fs4;
std::vector<double> HFBTHO::fs5;
std::vector<double> HFBTHO::fs6;
std::vector<double> HFBTHO::wdcor;
std::vector<double> HFBTHO::wdcori;
std::vector<double> HFBTHO::cou;

// std::vector<std::vector<double>> vDHartree, vhart00, vhart01, vhart11;
std::vector<std::vector<double>> HFBTHO::vDHartree;
std::vector<std::vector<double>> HFBTHO::vhart00;
std::vector<std::vector<double>> HFBTHO::vhart01;
std::vector<std::vector<double>> HFBTHO::vhart11;

//! PAV Projection
// int keypj, ilpj, ilpj2, ilnqx, ilnghl;
int HFBTHO::keypj;
int HFBTHO::ilpj;
int HFBTHO::ilpj2;
int HFBTHO::ilnqx;
int HFBTHO::ilnghl;

// double rehfbcan, ehfb, retotpj, depnp, iproj, npr1pj, npr2pj;
double HFBTHO::rehfbcan;
double HFBTHO::ehfb;
double HFBTHO::retotpj;
double HFBTHO::depnp;
double HFBTHO::iproj;
double HFBTHO::npr1pj;
double HFBTHO::npr2pj;

std::complex<double> HFBTHO::onei = (0.0, 1.0);
// std::vector<std::complex<double>> phypj, sinphy, exp1iphy, exp1iphym, exp2iphy, exp2iphym;
std::vector<std::complex<double>> HFBTHO::phypj;
std::vector<std::complex<double>> HFBTHO::sinphy;
std::vector<std::complex<double>> HFBTHO::exp1iphy;
std::vector<std::complex<double>> HFBTHO::exp1iphym;
std::vector<std::complex<double>> HFBTHO::exp2iphy;
std::vector<std::complex<double>> HFBTHO::exp2iphym;

// std::vector<std::vector<std::complex<double>>> coupj, pjk, epj;
std::vector<std::vector<std::complex<double>>> HFBTHO::coupj;
std::vector<std::vector<std::complex<double>>> HFBTHO::pjk;
std::vector<std::vector<std::complex<double>>> HFBTHO::epj;

// std::vector<std::vector<std::vector<std::complex<double>>>> ropj, taupj, dropj, djpj, akapj, SZFIpj, SFIZpj, SRFIpj, SFIRpj, ddepj, cpj, ypj, rpj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::ropj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::taupj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::dropj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::djpj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::akapj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::SZFIpj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::SFIZpj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::SRFIpj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::SFIRpj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::ddepj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::cpj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::ypj;
std::vector<std::vector<std::vector<std::complex<double>>>> HFBTHO::rpj;

// double polem[2], polep[2];
double HFBTHO::polem[2];
double HFBTHO::polep[2];

//! CMC
int HFBTHO::ICMinput;
// double ECMHFB[3], ECMPAV[3];
double HFBTHO::ECMHFB[3];
double HFBTHO::ECMPAV[3];

//! CRC
int HFBTHO::ICRinput;
// double DEROT[3], SQUJ[3], CRAN[3], ERIGHFB[3];
double HFBTHO::DEROT[3];
double HFBTHO::SQUJ[3];
double HFBTHO::CRAN[3];
double HFBTHO::ERIGHFB[3];

//! hfbdiagonal
// std::vector<double> erhfb, drhfb, erhfb1, drhfb1, evvk, evvkcan, zhfb;
std::vector<double> HFBTHO::erhfb;
std::vector<double> HFBTHO::drhfb;
std::vector<double> HFBTHO::erhfb1;
std::vector<double> HFBTHO::drhfb1;
std::vector<double> HFBTHO::evvk;
std::vector<double> HFBTHO::evvkcan;
std::vector<double> HFBTHO::zhfb;

// std::vector<std::vector<double>> hfb, hfbcan;
// std::vector<std::vector<double>> HFBTHO::hfb;
std::vector<double> HFBTHO::hfb;
std::vector<std::vector<double>> HFBTHO::hfbcan;

//! Broyden
char HFBTHO::bbroyden;
int HFBTHO::nbroyden = 7;
double HFBTHO::alphamix = 0.70;
// int nhhdim, nhhdim2, nhhdim3, nhhdim4, ialwork, ilwork;
int HFBTHO::nhhdim;
int HFBTHO::nhhdim2;
int HFBTHO::nhhdim3;
int HFBTHO::nhhdim4;
int HFBTHO::ialwork;
int HFBTHO::ilwork;

// std::vector<double> brout, brin;
std::vector<double> HFBTHO::brout;
std::vector<double> HFBTHO::brin;

std::vector<double> HFBTHO::alwork;
std::vector<int> HFBTHO::lwork;
//! new keys
// bool Parity, Parity_INI;
bool HFBTHO::Parity;
bool HFBTHO::Parity_INI;

bool HFBTHO::Print_Screen = true;
// bool Add_Pairing, Print_HFBTHO_Namelist;
bool HFBTHO::Add_Pairing;
bool HFBTHO::Print_HFBTHO_Namelist;

// int MAX_ITER_INI, keypj_INI, iproj_INI, npr1pj_INI, npr2pj_INI;
int HFBTHO::MAX_ITER_INI;
int HFBTHO::keypj_INI;
int HFBTHO::iproj_INI;
int HFBTHO::npr1pj_INI;
int HFBTHO::npr2pj_INI;

//! Eqp U,V
// int nuv, nqp;
int HFBTHO::nuv;
int HFBTHO::nqp;

// std::vector<double> RVqpN, RVqpP, RUqpN, RUqpP, REqpN, REqpP;
std::vector<double> HFBTHO::RVqpN;
std::vector<double> HFBTHO::RVqpP;
std::vector<double> HFBTHO::RUqpN;
std::vector<double> HFBTHO::RUqpP;
std::vector<double> HFBTHO::REqpN;
std::vector<double> HFBTHO::REqpP;

// std::vector<int> KpwiN, KpwiP, KqpN, KqpP;
std::vector<int> HFBTHO::KpwiN;
std::vector<int> HFBTHO::KpwiP;
std::vector<int> HFBTHO::KqpN;
std::vector<int> HFBTHO::KqpP;

//! error indicator
int HFBTHO::ierror_flag = 0;
std::string HFBTHO::ierror_info[11];
//! namelist
//   Namelist /HFBTHO_NAMELIST/ MAX_ITER_INI,epsi_INI,Add_Pairing_INI &
//        ,icou_INI,iLST_INI,keypj_INI,iproj_INI,npr1pj_INI,npr2pj_INI &
//        ,DO_FITT_INI,IDEBUG_INI,Parity_INI,Print_HFBTHO_Namelist_INI

//!
//! mpi setup
int HFBTHO::iam_mpi; //!! mpi id number of this process
// int num_processors_mpi, num_nodes_mpi, icom_err_mpi; //!! for mpi
int HFBTHO::num_processors_mpi;
int HFBTHO::num_nodes_mpi;
int HFBTHO::icom_err_mpi; //!! for mpi

void HFBTHO::read_HFBTHO_NAMELIST()
{
  //!------------------------------------
  iLST_INI = 0;                     //! 0:HO, -1:HO->THO, 1:THO
  MAX_ITER_INI = 201;               //! max number of iterations
  keypj_INI = 1;                    // ! number of gaouge points
  iproj_INI = 0;                    // ! projecting on different nucleus
  npr1pj_INI = 0;                   // !  its neutron number
  npr2pj_INI = 0;                   // !  its proton number
  epsi_INI = 0.0000001;             // ! relaxing tests
  Add_Pairing_INI = true;           // ! add pairing starting from file
  icou_INI = 2;                     // ! coul: no-(0), dir.only-(1), plus exchange-(2)
  DO_FITT_INI = false;              // ! calculates quantities for reg.optimization
  IDEBUG_INI = 0;                   // ! debug
  Parity_INI = true;                // ! reflection symmetry
  Print_HFBTHO_Namelist_INI = true; //! Print Namelist
  //!------------------------------------
  // iniialize_HFBTHO_SOLVER();
  //!
}

void HFBTHO::printHFBTHO()
{
  std::cout << "MAX_ITER_INI: " << MAX_ITER_INI << "Parity_INI: " << Parity_INI
            << "Do_FITT_INI: " << DO_FITT_INI << std::endl;
}
