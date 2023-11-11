#include <iostream>
#include <cstdlib>

#include "HFBTHO_solver.h"

void HFBTHO_solver::solver()
{
  //!-------------------------------------------------------------
  //! inatializing all according to *_INI values
  //!-------------------------------------------------------------
  std::cout << inin_INI << inin  << icstr << std::endl;
  iniialize_HFBTHO_SOLVER();
  // If(ierror_flag.Ne.0) Return
  // If(lout.Lt.lfile) Open(lfile,file='thoout.dat',status='unknown')
  std::cout << inin_INI << inin  << icstr << std::endl;
  Constraint_or_not(inin_INI, inin, icstr);
  std::cout << inin_INI << inin  << icstr << std::endl;
}

void HFBTHO_solver::iniialize_HFBTHO_SOLVER()
{
  //!------------------------------------
  //! tapes
  //!------------------------------------

  lwin = 21;
  lwou = 22;
  lwel = 12;
  lres = 17;
  lin = 3;
  lout = 6;

  //!------------------------------------
  //! From Namelist or default values
  //!------------------------------------

  epsi = epsi_INI;                                   //! stop criteria
  Add_Pairing = Add_Pairing_INI;                     //! add pairing starting from file
  icou = icou_INI;                                   //! coul: no-(0), dir.only-(1), plus exchange-(2)
  DO_FITT = DO_FITT_INI;                             //! calculates quantities for reg.optimization
  IDEBUG = IDEBUG_INI;                               //! debug
  Parity = Parity_INI;                               //! reflection symmetry
  Print_HFBTHO_Namelist = Print_HFBTHO_Namelist_INI; //! Print Namelist

  //!------------------------------------
  //! output control
  //!------------------------------------
  if (n00_INI > 0)
  {
    Print_Screen = true;
    lfile = lout + 1;
  } //! output to screen & thoout.dat
  else
  {
    Print_Screen = false;
    lfile = lout - 1; //! no output to screen & thoout.dat
  }

  //!------------------------------------
  //! Pi
  //!------------------------------------
  // PI=four*Atan(one)  // initalized in HFBTHO.cpp
  //!------------------------------------
  //! blocking
  //!------------------------------------
  // bloblo=0; blo123=0; blok1k2=0;  keyblo=0
  // blomax=0; nkblo=0;  iparenti=0; irestart=0
  // blocanon=zero;      eqpmin=zero;     // initialized in HFBTHO.cpp

  // !------------------------------------
  // ! buffers
  // !------------------------------------
  // eres = zero;
  // eresu = zero;
  // eresl = zero;
  // eresj = zero;
  // eresbl = zero;
  // ereslbl = ' 00[00,00,00]';

  // !------------------------------------
  // ! def parameters
  // !------------------------------------
  ffdef3 = sqrt(five / (four * PI)) / two;
  ffdef4 = sqrt(117.0) / (four * PI);
  ffdef5 = sqrt(nine / (four * PI)) / eight;
  ffdef6 = sqrt(five * PI) / three;
  ffdef7 = sqrt(PI) / four;

  // !------------------------------------
  // ! former linear mixing
  // !------------------------------------
  xmix0 = 0.30; //! lowest mixing parameter  (redefined later)
  xmix = 0.3;   //! initial mixing parameter (changes every iteration)
  xmax = 1.0;   //! mario
                // !------------------------------------
  // ! misck (redefined later)
  // !------------------------------------
  rehfbcan = 0.0;
  depnp = 0.0;
  // ala2 = 0.00;
  for (auto &m : ala2)
    m = 0.00;
  // ept = -2.0;
  for (auto &m : ept)
    m = -2.0;
  // del = 1.0;
  for (auto &m : del)
    m = 1.00;
  // ala = -7.0;
  for (auto &m : ala)
    m = -7.00;
  ala1[0] = -14.6851;
  ala1[1] = -3.7522;
  si = 1.0;
  iqrpa = 0;
  icacou = 0;
  icacoupj = 0;
  icahartree = 0;
  // iasswrong = 0 ;
  for (auto &m : iasswrong)
    m = 0;
  iError_in_HO = 0;
  iError_in_THO = 0;
  // ECMHFB = 0.0;
  for (auto &m : ECMHFB)
    m = 0.00;
  // ECMPAV = 0.0
  for (auto &m : ECMPAV)
    m = 0.00;

  // !------------------------------------
  // ! Saxon-Woods: von koepf und ring, z.phys. (1991)
  // !------------------------------------
  v0ws = -71.28;
  akv = 0.4616;
  // r0v = 1.2334;
  for (auto &m : r0v)
    m = 1.2334;
  // av = 0.6150;
  for (auto &m : av)
    m = 0.6150;
  // vso = 11.1175;
  for (auto &m : vso)
    m = 11.1175;
  // rso = 1.1443;
  for (auto &m : rso)
    m = 1.1443;
  // aso = 0.6476;
  for (auto &m : aso)
    m = 0.6476;

  //!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  //! fixed text
  //!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  tp[0] = '+';
  tp[1] = '-';
  tis[0] = 'n';
  tis[1] = 'p';
  tit[0] = "neutrons";
  tit[1] = "protons";
  tl[0] = 's';
  tl[1] = 'p';
  tl[2] = 'd';
  tl[3] = 'f';
  tl[4] = 'g';
  tl[5] = 'h';
  tl[6] = 'i';
  tl[7] = 'j';
  tl[8] = 'k';
  tl[9] = 'l';
  tl[10] = 'm';
  tl[11] = 'n';
  tl[12] = 'o';
  tl[13] = 'p';
  tl[14] = 'q';
  tl[15] = 'r';
  tl[16] = 's';
  tl[17] = 't';
  tl[18] = 'u';
  tl[19] = 'v';
  tl[20] = 'w';

  // !------------------------------------
  // ! fixed parity sign
  // !------------------------------------
  tpar[0] = +1;
  tpar[1] = -1;
  // !------------------------------------
  // ! physical constants
  // !------------------------------------
  amu = 938.90590;
  r0 = 1.20;
  alphi = 137.036020;
  hqc = 197.328910; // !hqc = 197.3282849;

  // !------------------------------------
  // ! e2 for protons (set now in elsewhere)
  // !------------------------------------
  // !chargee2=hqc/alphi
  // !chargee2=1.43997840_pr
  // !-----------------------------------
  // ! set the loops over particle types
  // !-----------------------------------
  itmin = 1;
  itmax = 2;
  if ((npr_INI[0] < 1) && (npr_INI[1] < 1))
  {
    ierror_flag = ierror_flag + 1;
    ierror_info[ierror_flag] = "STOP: Nonsense in particle numbers N<0 and Z<0.";
    exit;
  }
  if (npr_INI[0] == 0)
    itmin = 2;
  if (npr_INI[1] == 0)
    itmax = 1;

  // !-----------------------------------
  // ! error flag and info
  // !-----------------------------------
  ierror_flag = 0;
  ierror_info[ierror_flag] = "No errors in the solver!";

  set_functional_parameters(skyrme_INI, false);
};

void HFBTHO_solver::set_functional_parameters(std::string fname, bool lpr)
{
  // !--------------------------------------------------------------------------------
  // ! set functional parameters
  // !--------------------------------------------------------------------------------
  // Implicit None
  // Logical, Intent(in) :: lpr
  // Character (30), Intent(inout) :: fname
  bool regularization;
  const int lin = 15, lout = 6;
  //!
  //! parameters
  FunctionalName = fname;
  // eps=Spacing(1.0_pr);
  eps = 0.0000001;
  // Pi=4.0_pr*Atan(1.0_pr)
  kfconst = pow((1.50 * pow(PI, 2)), (1.0 / 3.0));
  //!(3Pi ^ 2 / 2) ^ (1 / 3)
  CK = 3.0 / 5.0 * pow(kfconst, 2);
  nuLambda = 700.0;
  nufpi = 93.0;

  //!
  // Make_Parameter_Free_Useful_Combintions();  // DEMORDER = -1, this function not executed.
  //!
  //! exact Hartree CHrho from INM
  CHrho = 0.0; //!!!!If (dmeorder.eq.3) Call CHrho_from_NM()

  // If(use_INM) Then
  //  Call calculate_C_form_NM(E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM)
  // Else
  Crho[0] = Crho[0] + CHrho; //!*(0.0_pr) !mario
  // End If
  // calculate_NM_properties();
  //!
  Crho[0] = Crho[0] - CHrho; //!*(0.00_pr) !mario
                             //!
  calculate_natural_units();
  //  !
  //  ! Print output
  //  if (lpr)
  print_functional_parameters(lout);
}

void HFBTHO_solver::calculate_natural_units()
{

  nuCrho[0] = Crho[0] * pow(nufpi, 2) / pow(mevfm, 3);
  nuCrho[1] = Crho[1] * pow(nufpi, 2) / pow(mevfm, 3);

  nuCdrho[0] = Cdrho[0] * pow(nufpi, 2) * pow((nuLambda * nufpi * nufpi), sigma) / pow(mevfm, (3.0 * (1.0 + sigma)));
  nuCdrho[1] = Cdrho[1] * pow(nufpi, 2) * pow((nuLambda * nufpi * nufpi), sigma) / pow(mevfm, (3.0 * (1.0 + sigma)));

  nuCtau[0] = Ctau[0] * pow((nufpi * nuLambda), 2) / pow(mevfm, 5);
  nuCtau[1] = Ctau[1] * pow((nufpi * nuLambda), 2) / pow(mevfm, 5);

  nuCrDr[0] = CrDr[0] * pow((nufpi * nuLambda), 2) / pow(mevfm, 5);
  nuCrDr[1] = CrDr[1] * pow((nufpi * nuLambda), 2) / pow(mevfm, 5);

  nuCrdJ[0] = CrdJ[0] * pow((nufpi * nuLambda), 2) / pow(mevfm, 5);
  nuCrdJ[1] = CrdJ[1] * pow((nufpi * nuLambda), 2) / pow(mevfm, 5);

  nuCJ[0] = CJ[0] * pow((nufpi * nuLambda), 2) / pow(mevfm, 5);
  nuCJ[1] = CJ[1] * pow((nufpi * nuLambda), 2) / pow(mevfm, 5);

  nuCpV0[0] = CpV0[0] * pow(nufpi, 2) / pow(mevfm, 3);
  nuCpV0[1] = CpV0[1] * pow(nufpi, 2) / pow(mevfm, 3);

  nuCpV1[0] = CpV1[0] * pow(nufpi, 4) * nuLambda / pow(mevfm, 6);
  nuCpV1[1] = CpV1[1] * pow(nufpi, 4) * nuLambda / pow(mevfm, 6);
};

void HFBTHO_solver::print_functional_parameters(int lout)
{
  t_from_C();
  std::cout << "print_functional_parameters()" << std::endl;
};

void HFBTHO_solver::t_from_C()
{
  t0 = (8.0 / 3) * Crho[0];
  t1 = 4.0 / 3.0 * (Ctau[0] - 4.0 * CrDr[0]);
  t2 = 4.0 / 3.0 * (3.0 * Ctau[0] - 6.0 * Ctau[1] + 4.0 * CrDr[0] - 8.0 * CrDr[1]);
  t3 = 16.0 * Cdrho[0];
  x0 = -0.50 * (3.0 * Crho[1] / Crho[0] + 1.0);
  x1 = 2.0 * (-Ctau[0] - 3.0 * Ctau[1] + 4.0 * CrDr[0] + 12.0 * CrDr[1]) / t1 / 3.0;
  x2 = -2.0 * (3.0 * Ctau[0] - 15.0 * Ctau[1] + 4.0 * CrDr[0] - 20.0 * CrDr[1]) / t2 / 3.0;
  x3 = -0.50 * (3.0 * Cdrho[1] / Cdrho[0] + 1.0);
  b4 = CrdJ[1] - CrdJ[0];
  b4p = -2.0 * CrdJ[1];
  te = (4.0 / 15.0) * (3.0 * CJ[0] - 9.0 * CJ[1] - 4.0 * CrDr[0] + 12.0 * CrDr[1] - 2.0 * Ctau[0] + 6.0 * Ctau[1]);
  to = (4.0 / 15.0) * (3.0 * CJ[0] + 3.0 * CJ[1] + 4.0 * CrDr[0] + 4.0 * CrDr[1]);
};

void HFBTHO_solver::Constraint_or_not(int &inin_INI0, int &inin0, int &icstr0)
{
  if (std::abs(inin_INI0) > 100)
  {
    icstr0 = 1;
    inin0 = inin_INI0 / 100;
  }
  else
  {
    icstr0 = 0;
    inin0 = inin_INI0;
  }
};
