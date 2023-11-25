#include <iostream>
#include <cstdlib>
#include <cstring>
#include <iomanip>

#include <algorithm>
#include <bits/stdc++.h>

#include "HFBTHO_solver.h"
#include "lapack.h"

void HFBTHO_solver::solver()
{
  //!-------------------------------------------------------------
  //! inatializing all according to *_INI values
  //!-------------------------------------------------------------
  std::cout << inin_INI << inin << icstr << std::endl;
  iniialize_HFBTHO_SOLVER();
  // If(ierror_flag.Ne.0) Return
  // If(lout.Lt.lfile) Open(lfile,file='thoout.dat',status='unknown')
  std::cout << inin_INI << inin << icstr << std::endl;
  std::cout << "Inside solver" << std::endl;
  Constraint_or_not(inin_INI, inin, icstr);
  std::cout << inin_INI << inin << icstr << std::endl;
  //-------------------------------------------------------------------------
  // Loop recalculating eventually the even-even solution for an odd nucleus
  //-------------------------------------------------------------------------
  do
  {
    irestart = 0;
    n00 = std::abs(n00_INI);
    b0 = b0_INI;
    q = q_INI;
    iLST = iLST_INI;
    maxi = MAX_ITER_INI;
    npr[0] = npr_INI[0];
    npr[1] = npr_INI[1];
    npr[2] = npr[0] + npr[1];
    strncpy(skyrme, skyrme_INI, 30);
    kindhfb = kindhfb_INI;
    cdef = cdef_INI;
    cqad = cqad_INI;
    keypj = keypj_INI;
    iproj = iproj_INI;
    npr1pj = npr1pj_INI;
    npr2pj = npr2pj_INI;
    // nkblo=nkblo_INI;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 5; j++)
        nkblo[i][j] = nkblo_INI[i][j];

    // blocking
    // Do .. end do loop drop for now. Dont care blocking case.

    //-------------------------------------------------------------
    // HFB+HO calculations
    //-------------------------------------------------------------
    if (iLST <= 0)
    {
      icacou = 0;
      icahartree = 0;
      preparer(true);
      // inout(1);
    }
  } while (irestart != 0);
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
  if (std::abs(inin_INI0) >= 100)
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

void HFBTHO_solver::preparer(bool lpr)
{
  // If n00 skip
  // call nucleus skip
  // If lpr skip

  //-----------------------------------------
  // pairing parameters (NB! modify later)
  //-----------------------------------------
  // rho_c=0.160; pwi=2000.0 !60.0;
  rho_c = 0.160;
  pwi = 60.0;
  //-----------------------------------------
  // particle number as real variable
  //-----------------------------------------
  tz[0] = static_cast<double>(npr[0]);
  tz[1] = static_cast<double>(npr[1]);
  amas = tz[0] + tz[1];
  //-----------------------------------------
  // default combinations
  //-----------------------------------------
  r00 = r0 * pow(amas, p13);
  r02 = pow(r00, 2);
  r04 = pow(r02, 2);
  chargee2 = e2charg;
  coex = -chargee2 * pow((3.0 / PI), p13);
  cex = -0.750 * coex;
  hom = 41.0 * pow(amas, (-p13)) * r0;
  //-----------------------------------------
  // hbzero from forces [hqc**2/(two*amu)]
  //-----------------------------------------
  hb0 = hbzero;
  if (use_cm_cor)
    hb0 = hb0 * (one - one / amas);
  //-----------------------------------------
  // basis parameter q
  //-----------------------------------------
  beta0 = q;
  q = exp((3.0 * sqrt(5.0 / (16.0 * PI))) * beta0);
  //-----------------------------------------
  // basis parameters b0,bp,bz
  //-----------------------------------------
  if (b0 <= zero)
    b0 = sqrt(two * hbzero / hom);
  if (b0 <= zero)
    b0 = sqrt(pow(hqc, 2) / (hom * amu));
  bp = b0 * pow(q, (-one / 6.0));
  bz = b0 * pow(q, (one / 3.0));
  bpp = bp * bp;
  //-----------------------------------------
  // constraint in terms of beta
  //-----------------------------------------
  ty20 = sqrt(5.0 / PI) * hom / pow(b0, 2) / two;
  //-----------------------------------------
  // projection: number of grid points
  //-----------------------------------------
  keypj = std::max(1, keypj);
  ilpj = keypj;
  ilpj2 = pow(ilpj, 2);
  //-----------------------------------------
  // projecting on different nucleus
  //-----------------------------------------
  if (iproj == 0)
  {
    npr1pj = npr[0];
    npr2pj = npr[1];
  }
  else
  {
    npr1pj = npr[0] + npr1pj;
    npr2pj = npr[2] + npr2pj;
  }
  //-----------------------------------------
  // blocking window
  //-----------------------------------------
  pwiblo = std::min(std::max(25.0 / sqrt(static_cast<float>(npr[0] + npr[1])), 2.0), 8.0);
  //---------------------------------------------------------
  // statistics to screen('lout')/file('lfile')
  //---------------------------------------------------------
  //-----------------------------------------
  // BASIS, GAUSS POINTS, HOWF
  //-----------------------------------------
  gfv(); // factorials
         // if lpr skip
  base0(lpr);
  // std::cout << "Outside base0: " << std::endl;
  // std::cout << std::setw(8) << std::left << ndx
  //           << std::setw(8) << std::left << nqp
  //           << std::setw(8) << std::left << nuv
  //           << std::setw(8) << std::left << nbx
  //           << std::setw(8) << std::left << ntx;

  // std::cout << "Outside gfv" << std::endl;
  // for (int i = 0; i < iv.size() ; i++)
  //{
  //   std::cout << std::setw(4) << std::right << static_cast<double>(iv[i]) << " " << std::setw(15) << std::left << sq[i] << " " << std::setw(15) << std::left << fak[i] << " " << std::setw(15) << wf[i] << std::endl;
  // }
  thoalloc();
  gausspoints();
}

void HFBTHO_solver::gfv()
{
  //---------------------------------------------------------------------
  // Calculates sign, Sqrt, factorials, etc. of integers and half int.
  // iv(n)=(-1)**n, sq(n)=Sqrt(n), sqi(n)=1/Sqrt(n)
  // fak(n)=n!; wf(n)=Sqrt(n!); wfi(n)=1/Sqrt(n!)
  //---------------------------------------------------------------------
  // Use HFBTHO
  // Implicit None
  int igfv = 170; // maximal number for GFV
  // If(Allocated(iv)) Deallocate(iv,fak,fi,sq,sqi,wf,wfi)
  // Allocate(iv(-igfv:igfv),fak(0:igfv),fi(0:igfv),sq(0:igfv),sqi(0:igfv))
  // Allocate(wf(0:igfv),wfi(0:igfv))
  iv.resize(igfv + 1); // In Fortran version size of iv is iv(-igfv:igfv), in c++ we can't do this. We use the relation iv(-i) = iv(i) to get around.
  fak.resize(igfv + 1);
  fi.resize(igfv + 1);
  sq.resize(igfv + 1);
  sqi.resize(igfv + 1);
  wf.resize(igfv + 1);
  wfi.resize(igfv + 1);

  iv[0] = 1;
  sq[0] = zero;
  sqi[0] = 1.0;
  fak[0] = one;
  fi[0] = one;
  wf[0] = one;
  wfi[0] = one;

  for (int i = 1; i < igfv + 1; i++)
  {
    iv[i] = -iv[i - 1];
    // iv(-i) = iv(i)
    sq[i] = sqrt(static_cast<double>(i));
    sqi[i] = one / sq[i];
    fak[i] = static_cast<double>(i) * fak[i - 1];
    fi[i] = one / fak[i];
    wf[i] = sq[i] * wf[i - 1];
    wfi[i] = one / wf[i];
  }
  std::cout << "Inside gfv" << std::endl;
  std::cout << "Inside gfv" << iv.size() << std::endl;
  // for (int i = 0; i <= iv.size(); i++)
  //{
  //   std::cout << std::setw(4) << std::left <<  i << std::setw(4) << std::right << static_cast<double>(iv[i]) << " " << std::setw(15) << std::left << sq[i] << " " << std::setw(15) << std::left << fak[i] << " " << std::setw(15) << wf[i] << std::endl;
  // }
}

//=======================================================================
void HFBTHO_solver::base0(bool lpr)
{
  //---------------------------------------------------------------------
  // selects HO basis configurations in cylindrical coordinates
  //---------------------------------------------------------------------
  // Use HFBTHO
  // Use HFBTHO_ord
  // Implicit None
  // Logical :: lpr
  int iw, k, nre, nze, la, le, ip, ir, iz, il, is, Iall, ilauf, jlauf, ib, nd;
  int NOSCIL;
  std::vector<double> e;
  double hbz, hbp, ee;
  std::cout << "base0 start" << std::endl;
  //
  // if n00.Gt.n00max skip
  //-----------------------------------------------
  // MAXIMUM NUMBER OF THE HO SHELLS (n00,NOSCIL)
  // (7,120),(8,165),(9,220),(10,286),(11,364)
  // (12,455),(14,680),(16,969),(18,1330),(20,1771)
  //-----------------------------------------------
  NOSCIL = (n00 + 1) * (n00 + 2) * (n00 + 3) / 6;
  //-----------------------------------------------
  // count all states for n00max
  //-----------------------------------------------
  nze = n00max;
  nre = n00max / 2;
  Iall = 0;
  for (int k = 1; k <= n00max + 1; k++)
  {
    la = k - 1;
    le = std::min(n00max, k);
    for (int ip = 1; ip <= 2; ip++)
    {
      for (int ir = 0; ir <= nre; ir++)
      {
        for (int iz = 0; iz <= nze; iz++)
        {
          for (int il = la; il <= le; il++)
          {
            for (int is = +1; is >= -1; is = is - 2)
            {
              if (iz + 2 * ir + il > n00max)
                continue;
              if (il + (is + 1) / 2 != k)
                continue;
              if ((iz + il) % 2 != (ip - 1))
                continue;
              Iall = Iall + 1;
              // std::cout << std::setw(4) << std::left << k << std::setw(4) << std::left << ip << std::setw(4) << std::left << ir << std::setw(4) << std::left << iz << std::setw(4) << std::left << il << std::setw(4) << std::left << is << std::setw(4) << std::left << Iall << std::endl;
            }
          }
        }
      }
    }
  };
  std::cout << "1qq" << std::endl;
  //-----------------------------------------------
  //  charge all energies for n00max
  //-----------------------------------------------
  //  Allocate(e(Iall))
  e.resize(Iall);
  std::cout << "hbzero: " << hbzero << "two: " << two << std::endl;
  std::cout << "bz: " << bz << "bp: " << bp << std::endl;
  hbz = two * hbzero / pow(bz, 2);
  hbp = two * hbzero / pow(bp, 2);
  std::cout << "hbz: " << hbz << "hbp: " << hbp << std::endl;
  Iall = 0;
  for (int k = 1; k <= n00max + 1; k++)
  {
    la = k - 1;
    le = std::min(n00max, k);
    for (int ip = 1; ip <= 2; ip++)
    {
      for (int ir = 0; ir <= nre; ir++)
      {
        for (int iz = 0; iz <= nze; iz++)
        {
          for (int il = la; il <= le; il++)
          {
            for (int is = +1; is >= -1; is = is - 2)
            {
              if (iz + 2 * ir + il > n00max)
                continue;
              if (il + (is + 1) / 2 != k)
                continue;
              if ((iz + il) % 2 != ip - 1)
                continue;
              Iall = Iall + 1;
              // std::cout << std::setw(8) << std::left << Iall << std::endl;
              // Iall -1 because e in c++ start from zero not one
              e[Iall - 1] = hbz * (iz + half) + hbp * (two * ir + il + one);
              // std::cout << std::setw(8) << std::left << Iall << std::setw(8) << std::left << e[Iall] << std::endl;
            }
          }
        }
      }
    }
  }
  //-----------------------------------------------
  // sort energies  and derive base cut-off energy
  //-----------------------------------------------
  // std::cout << "e before order" << std::endl;
  // for (int i = 0; i < e.size(); i++)
  //{
  //  std::cout << std::setw(8) << std::left << i << e[i] << std::endl;
  //}

  ord(Iall, e);
  // std::cout << "e after order" << std::endl;
  // for (int i = 0; i < e.size(); i++)
  //{
  //   std::cout << std::setw(8) << std::left << i << std::setprecision(8) <<  e[i] << std::endl;
  // }
  EBASECUT = e[NOSCIL - 1] + 1.0 * exp(-5);
  std::cout << "NOSCIL: " << NOSCIL << std::endl;
  std::cout << "EBASECUT: " << EBASECUT << std::endl;
  //-----------------------------------------------
  // calculate the actual states
  //-----------------------------------------------
  nze = n00max;
  nre = n00max / 2;
  ib = 0;
  ilauf = 0;
  nzx = 0;
  nrx = 0;
  nlx = 0;
  nqp = 0;
  nuv = 0;
  // loop over k-quantum number
  for (int k = 1; k <= n00max + 1; k++)
  {
    la = k - 1;
    le = std::min(n00max, k);
    // std::cout << "Parity: " << Parity << std::endl;
    //  loop over parity
    if (!Parity)
      jlauf = ilauf; // Nop
    for (int ip = 1; ip <= 2; ip++)
    {
      if (Parity)
        jlauf = ilauf; // Yesp
      // std::cout << "jlauf: " << jlauf << std::endl;
      for (int ir = 0; ir <= nre; ir++)
      {
        for (int iz = 0; iz <= nze; iz++)
        {
          for (int il = la; il <= le; il++)
          {
            for (int is = +1; is >= -1; is = is - 2)
            {
              if (iz + 2 * ir + il > n00max)
                continue;
              if (il + (is + 1) / 2 != k)
                continue;
              if ((iz + il) % 2 != ip - 1)
                continue;
              ee = hbz * (iz + half) + hbp * (two * ir + il + one);
              // std::cout << "ee: " << ee << std::endl;
              if (ee <= EBASECUT)
              {
                std::cout << "ee: " << ee << std::endl;
                ilauf = ilauf + 1;
                nzx = std::max(nzx, iz);
                nrx = std::max(nrx, ir);
                nlx = std::max(nlx, il);
              }
            }
          }
        }
      }
      if (Parity)
      { // Yesp
        if (ilauf > jlauf)
        {
          ib = ib + 1;
          nd = ilauf - jlauf;
          ndx = std::max(ndx, nd);
          nqp = nqp + nd;
          nuv = nuv + nd * nd;
        }
      }
    }
    if (!Parity)
    {
      if (ilauf >= jlauf)
      {
        ib = ib + 1;
        nd = ilauf - jlauf;
        ndx = std::max(ndx, nd);
        nqp = nqp + nd;
        nuv = nuv + nd * nd;
      }
    }
  }
  nbx = ib;
  ntx = ilauf;
  std::cout << std::setw(4) << std::left << ib
            << std::setw(8) << std::left << nd
            << std::setw(8) << std::left << ndx
            << std::setw(8) << std::left << nqp
            << std::setw(8) << std::left << nuv
            << std::setw(8) << std::left << nbx
            << std::setw(8) << std::left << ntx;
  //-----------------------------------------------
  // print statistics
  //-----------------------------------------------
  // if lpr  skip
}

void HFBTHO_solver::ord(int n, std::vector<double> &e)
{
  //---------------------------------------------------------------------
  // orders a set of numbers according to their size
  //---------------------------------------------------------------------
  // Use HFBTHO
  // Implicit None
  // Integer(ipr) :: n,i,k,j
  // Real(pr) :: e(n)
  int k;
  double p;
  for (int i = 0; i < n; i++)
  {
    k = i;
    p = e[i];
    if (i < n)
    {
      for (int j = i + 1; j < n; j++)
      {
        if (e[j] <= p)
        {
          k = j;
          p = e[j];
        }
      }
      if (k != i)
      {
        e[k] = e[i];
        e[i] = p;
      }
    }
  }
}

void HFBTHO_solver::thoalloc()
{
  //! Use HFBTHO
  //  Implicit None
  // int ier, ib, ND;
  int ier;
  //  !
  //  ! number of int.points
  if (Parity)
  {
    //!     ngh=30; ngl=30; nleg=30                                    !Yesp
    //     !value from CPC(2013)
    ngh = 40;
    ngl = 40;
    nleg = 80; //! Nop
  }
  else
  {
    //!     ngh=60; ngl=30; nleg=30                                    !Nop
    //     !value from CPC(2013)
    ngh = 80;
    ngl = 40;
    nleg = 80; //! Nop
  }

  std::cout << "ilpj" << std::endl;
  std::cout << std::setw(8) << std::left << ilpj;

  nzrlx = (nzx + 1) * (nrx + 1) * (nlx + 1); //! phy(:,:,nzrlx)
  nghl = ngh * ngl;                          //! nghl=ngh*ngl
  nqx = ndx * ndx;
  nb2x = nbx + nbx;
  ndx2 = ndx + ndx;
  ilnqx = ilpj * nqx;
  ilnghl = ilpj * nghl;
  nhfbx = ndx + ndx;
  nhfbqx = nhfbx * nhfbx;
  nkx = ntx;
  ndxs = ndx * (ndx + 1) / 2;
  std::cout << "thoalloc" << std::endl;
  std::cout << std::setw(8) << std::left << nzrlx
            << std::setw(8) << std::left << nghl
            << std::setw(8) << std::left << nqx
            << std::setw(8) << std::left << nb2x
            << std::setw(8) << std::left << ndx2
            << std::setw(8) << std::left << ilnqx
            << std::setw(8) << std::left << ilnghl
            << std::setw(8) << std::left << nhfbx
            << std::setw(8) << std::left << nhfbqx
            << std::setw(8) << std::left << nkx
            << std::setw(8) << std::left << ndxs;

  //!-----------------------------------------
  //! Arrays depending on gauss points
  //!-----------------------------------------
  // If(Allocated(xh)) Deallocate(xh,wh,xl,sxl,wl,xleg,wleg,vc &
  //      ,vhbn,vn,vrn,vzn,vdn,vsn,dvn,vhbp,vp,vrp,vzp,vdp,vsp,dvp  &
  //      ,vSZFIn,vSFIZn,vSRFIn,vSFIRn,vSZFIp,vSFIZp,vSRFIp,vSFIRp &
  //      ,fl,fli,fh,fd,fp1,fp2,fp3,fp4,fp5,fp6  &
  //      ,fs1,fs2,fs3,fs4,fs5,fs6,wdcor,wdcori,cou,vDHartree,vhart00,vhart01,vhart11)
  // Allocate(xh(ngh),wh(ngh),xl(ngl),sxl(ngl),wl(ngl),xleg(nleg),wleg(nleg),vc(nghl,nghl))
  xh.resize(ngh);
  wh.resize(ngh);
  xl.resize(ngl);
  sxl.resize(ngl);
  wl.resize(ngl);
  xleg.resize(nleg);
  wleg.resize(nleg);
  vc.resize(nghl, std::vector<double>(nghl));

  // Allocate(vhbn(nghl),vn(nghl),vrn(nghl),vzn(nghl),vdn(nghl),vsn(nghl),dvn(nghl)  &
  //      ,vhbp(nghl),vp(nghl),vrp(nghl),vzp(nghl),vdp(nghl),vsp(nghl),dvp(nghl)  &
  //      ,vSZFIn(nghl),vSFIZn(nghl),vSRFIn(nghl),vSFIRn(nghl)  &
  //      ,vSZFIp(nghl),vSFIZp(nghl),vSRFIp(nghl),vSFIRp(nghl))

  vhbn.resize(nghl);
  vn.resize(nghl);
  vrn.resize(nghl);
  vzn.resize(nghl);
  vdn.resize(nghl);
  vsn.resize(nghl);
  dvn.resize(nghl);
  ;
  vhbp.resize(nghl);
  vp.resize(nghl);
  vrp.resize(nghl);
  vzp.resize(nghl);
  vdp.resize(nghl);
  vsp.resize(nghl);
  dvp.resize(nghl);
  ;
  vSZFIn.resize(nghl);
  vSFIZn.resize(nghl);
  vSRFIn.resize(nghl);
  vSFIRn.resize(nghl);
  ;
  vSZFIp.resize(nghl);
  vSFIZp.resize(nghl);
  vSRFIp.resize(nghl);
  vSFIRp.resize(nghl);

  // Allocate(fl(nghl),fli(nghl),fh(nghl),fd(nghl),fp1(nghl),fp2(nghl),fp3(nghl)  &
  //      ,fp4(nghl),fp5(nghl),fp6(nghl),fs1(nghl),fs2(nghl),fs3(nghl),fs4(nghl)  &
  //      ,fs5(nghl),fs6(nghl),wdcor(nghl),wdcori(nghl),cou(nghl),vDHartree(nghl,2) &
  //      ,vhart00(nghl,nghl),vhart01(nghl,nghl),vhart11(nghl,nghl))

  fl.resize(nghl);
  fli.resize(nghl);
  fh.resize(nghl);
  fd.resize(nghl);
  fp1.resize(nghl);
  fp2.resize(nghl);
  fp3.resize(nghl);
  fp4.resize(nghl);
  fp5.resize(nghl);
  fp6.resize(nghl);
  fs1.resize(nghl);
  fs2.resize(nghl);
  fs3.resize(nghl);
  fs4.resize(nghl);
  fs5.resize(nghl);
  fs6.resize(nghl);
  wdcor.resize(nghl);
  wdcori.resize(nghl);
  cou.resize(nghl);
  vDHartree.resize(nghl, std::vector<double>(2));
  vhart00.resize(nghl, std::vector<double>(nghl));
  vhart01.resize(nghl, std::vector<double>(nghl));
  vhart11.resize(nghl, std::vector<double>(nghl));

  // If(Allocated(aka)) Deallocate(aka,ro,tau,dro,dj,NABLAR,NABLAZ,SZFI,SFIZ,SRFI,SFIR)
  // Allocate(aka(nghl,2),ro(nghl,2),tau(nghl,2),dro(nghl,2),dj(nghl,2)  &
  //      ,SZFI(nghl,2),SFIZ(nghl,2),SRFI(nghl,2),SFIR(nghl,2)  &
  //      ,NABLAR(nghl,2),NABLAZ(nghl,2))

  // If(Allocated(aka)) Deallocate(aka,ro,tau,dro,dj,NABLAR,NABLAZ,SZFI,SFIZ,SRFI,SFIR)
  aka.resize(nghl, std::vector<double>(2));
  ro.resize(nghl, std::vector<double>(2));
  tau.resize(nghl, std::vector<double>(2));
  dro.resize(nghl, std::vector<double>(2));
  dj.resize(nghl, std::vector<double>(2));
  SZFI.resize(nghl, std::vector<double>(2));
  SFIZ.resize(nghl, std::vector<double>(2));
  SRFI.resize(nghl, std::vector<double>(2));
  SFIR.resize(nghl, std::vector<double>(2));
  NABLAR.resize(nghl, std::vector<double>(2));
  NABLAZ.resize(nghl, std::vector<double>(2));

  //-----------------------------------------
  // Arrays depending on configurations
  //-----------------------------------------
  // If(Allocated(rk)) Deallocate(rk,ak,qh,qh1,ql,ql1,nz,nr,nl,ns,npar,id  &
  //     ,ia,ikb,ipb,ka,kd,tb,txb,numax,ek,dk,vk,vk1,uk,vkmax,ddc,ddc1,hfb1,lcanon)
  // Allocate(rk(nqx,nb2x),ak(nqx,nb2x),qh(0:nzx,1:ngh+1)  &
  //     ,qh1(0:nzx,1:ngh+1),ql(0:nrx,0:nlx,1:ngl+1),ql1(0:nrx,0:nlx,1:ngl+1)  &
  //     ,nz(ntx),nr(ntx),nl(ntx),ns(ntx),npar(ntx),id(nbx),ia(nbx),ikb(nbx),lcanon(0:nbx,2)  &
  //     ,ipb(nbx),ka(nbx,2),kd(nbx,2),tb(ntx),txb(nbx),numax(0:nkx,2)  &
  //     ,ek(nkx,2),dk(nkx,2),vk(nkx,2),vk1(nkx,2),uk(nkx,2),vkmax(nkx,2)  &
  //     ,ddc(ndx,nkx,2),ddc1(ndx,nkx,2),hfb1(nhfbx,2))

  //-----------------------------------------
  // Arrays depending on configurations
  //-----------------------------------------
  // If(Allocated(rk)) Deallocate(rk,ak,qh,qh1,ql,ql1,nz,nr,nl,ns,npar,id  &
  //     ,ia,ikb,ipb,ka,kd,tb,txb,numax,ek,dk,vk,vk1,uk,vkmax,ddc,ddc1,hfb1,lcanon)
  rk.resize(nqx, std::vector<double>(nb2x));
  ak.resize(nqx, std::vector<double>(nb2x));
  qh.resize(nzx + 1, std::vector<double>(ngh + 1));
  qh1.resize(nzx + 1, std::vector<double>(ngh + 1));
  ql.resize(nrx + 1, std::vector<std::vector<double>>(nlx + 1, std::vector<double>(ngl + 1)));
  ql1.resize(nrx + 1, std::vector<std::vector<double>>(nlx + 1, std::vector<double>(ngl + 1)));
  nz.resize(ntx);
  nr.resize(ntx);
  nl.resize(ntx);
  ns.resize(ntx);
  npar.resize(ntx);
  id.resize(nbx);
  ia.resize(nbx);
  ikb.resize(nbx);
  lcanon.resize(nbx + 1, std::vector<int>(2));
  ipb.resize(nbx);
  ka.resize(nbx, std::vector<int>(2));
  kd.resize(nbx, std::vector<int>(2));
  tb.resize(ntx);
  txb.resize(nbx);
  numax.resize(nkx + 1, std::vector<int>(2));
  ek.resize(nkx, std::vector<double>(2));
  dk.resize(nkx, std::vector<double>(2));
  vk.resize(nkx, std::vector<double>(2));
  vk1.resize(nkx, std::vector<double>(2));
  uk.resize(nkx, std::vector<double>(2));
  vkmax.resize(nkx, std::vector<double>(2));
  ddc.resize(ndx, std::vector<std::vector<double>>(nkx, std::vector<double>(2)));
  ddc1.resize(ndx, std::vector<std::vector<double>>(nkx, std::vector<double>(2)));
  hfb1.resize(nhfbx, std::vector<double>(2));

  //-----------------------------------------
  // HFB Arrays
  //-----------------------------------------
  // If(Allocated(erhfb)) Deallocate(erhfb,drhfb,erhfb1,drhfb1)
  // Allocate(erhfb(nkx),drhfb(nkx),erhfb1(nkx),drhfb1(nkx))
  // If(Allocated(hfb)) Deallocate(hfb,zhfb,evvk,hfbcan,evvkcan)
  // Allocate(hfb(ndx2,ndx2),zhfb(ndx2),evvk(ndx2),hfbcan(ndx,ndx),evvkcan(ndx))
  // If(Allocated(AN)) Deallocate(AN,ANk,PFIU,PFID,FIU,FID,FIUR,FIDR,FIUD2N,FIDD2N,FIUZ,FIDZ)
  // Allocate(AN(nqx),ANk(nqx),PFIU(ndx),PFID(ndx),FIU(ndx),FID(ndx)  &
  //     ,FIUR(ndx),FIDR(ndx),FIUD2N(ndx),FIDD2N(ndx),FIUZ(ndx),FIDZ(ndx))

  //-----------------------------------------
  // HFB Arrays
  //-----------------------------------------
  // If(Allocated(erhfb)) Deallocate(erhfb,drhfb,erhfb1,drhfb1)
  erhfb.resize(nkx);
  drhfb.resize(nkx);
  erhfb1.resize(nkx);
  drhfb1.resize(nkx);
  // If(Allocated(hfb)) Deallocate(hfb,zhfb,evvk,hfbcan,evvkcan)
  // hfb.resize(ndx2, std::vector<double>(ndx2));
  hfb.resize(ndx2 * ndx2);
  zhfb.resize(ndx2);
  evvk.resize(ndx2);
  hfbcan.resize(ndx, std::vector<double>(ndx));
  evvkcan.resize(ndx);
  // If(Allocated(AN)) Deallocate(AN,ANk,PFIU,PFID,FIU,FID,FIUR,FIDR,FIUD2N,FIDD2N,FIUZ,FIDZ)
  AN.resize(nqx);
  ANK.resize(nqx);
  PFIU.resize(ndx);
  PFID.resize(ndx);
  FIU.resize(ndx);
  FID.resize(ndx);
  FIUR.resize(ndx);
  FIDR.resize(ndx);
  FIUD2N.resize(ndx);
  FIDD2N.resize(ndx);
  FIUZ.resize(ndx);
  FIDZ.resize(ndx);

  //!-----------------------------------------
  //! Optimal LAPACK storage
  //!-----------------------------------------
  // ialwork=1; ilwork=1;
  // If(Allocated(alwork)) Deallocate(alwork,lwork)
  // Allocate(alwork(ialwork),lwork(ilwork))
  // ier=0; Call DSYEVD('V','L',ndx2,hfb,ndx2,evvk,ALWORK,-1,LWORK,-1,ier)
  // If(ier.Ne.0) Then
  //    ierror_flag=ierror_flag+1
  //    ierror_info(ierror_flag)='STOP: FATAL ERROR CONDITION IN DSYEVD'
  //    Return
  // Endif
  // ialwork=Int(alwork(1)); ilwork=lwork(1)
  // If(Allocated(alwork)) Deallocate(alwork,lwork)
  // Allocate(alwork(ialwork),lwork(ilwork))

  //!-----------------------------------------
  //! Optimal LAPACK storage
  //!-----------------------------------------
  ialwork = 1;
  ilwork = 1;
  // If(Allocated(alwork)) Deallocate(alwork,lwork)
  alwork.resize(ialwork);
  lwork.resize(ilwork);
  ier = 0;

  const int a = -1;
  const int b = -1;
  // LAPACK_dsyevd("V", "L", &ndx2, hfb, &ndx2, evvk, alwork, &a, lwork, &b, &ier);
  LAPACK_dsyevd("V", "L", &ndx2, &hfb[0], &ndx2, &evvk[0], &alwork[0], &a, &lwork[0], &b, &ier);
  // If(ier.Ne.0) Then
  //    ierror_flag=ierror_flag+1
  //    ierror_info(ierror_flag)='STOP: FATAL ERROR CONDITION IN DSYEVD'
  //    Return
  // Endif
  ialwork = static_cast<int>(alwork[0]);
  ilwork = lwork[0];
  // If(Allocated(alwork)) Deallocate(alwork,lwork)
  alwork.resize(ialwork);
  lwork.resize(ilwork);

  //!-----------------------------------------
  //! Eqp, U,V
  //!-----------------------------------------
  // If(Allocated(RVqpN)) Deallocate(RVqpN,RVqpP,RUqpN,RUqpP,REqpN,REqpP)
  // Allocate(RVqpN(nuv),RVqpP(nuv),RUqpN(nuv),RUqpP(nuv),REqpN(nqp),REqpP(nqp))
  // If(Allocated(KpwiP)) Deallocate(KpwiP,KpwiN,KqpN,KqpP)
  // Allocate(KpwiN(nqp),KpwiP(nqp),KqpN(nqp),KqpP(nqp))

  //!-----------------------------------------
  //! Eqp, U,V
  //!-----------------------------------------
  // If(Allocated(RVqpN)) Deallocate(RVqpN,RVqpP,RUqpN,RUqpP,REqpN,REqpP)
  RVqpN.resize(nuv);
  RVqpP.resize(nuv);
  RUqpN.resize(nuv);
  RUqpP.resize(nuv);
  REqpN.resize(nqp);
  REqpP.resize(nqp);
  // If(Allocated(KpwiP)) Deallocate(KpwiP,KpwiN,KqpN,KqpP)
  KpwiN.resize(nqp);
  KpwiP.resize(nqp);
  KqpN.resize(nqp);
  KqpP.resize(nqp);

  //!-----------------------------------------
  //! PNP ARRAYS: CONF. AND GAUGE ANGLE
  //!-----------------------------------------
  // If(Allocated(exp1iphy))Deallocate(ropj,taupj,dropj,djpj,akapj,coupj,pjk  &
  //      ,SZFIpj,SFIZpj,SRFIpj,SFIRpj,epj,cpj,ypj,rpj,ddepj,phypj,sinphy  &
  //      ,exp1iphy,exp2iphy,exp1iphym,exp2iphym)
  // ropj(nghl,ilpj,2);taupj(nghl,ilpj,2);dropj(nghl,ilpj,2)
  //      ;djpj(nghl,ilpj,2);akapj(nghl,ilpj,2);coupj(nghl,ilpj);pjk(ilpj,2)
  //      ;SZFIpj(nghl,ilpj,2);SFIZpj(nghl,ilpj,2);SRFIpj(nghl,ilpj,2)
  //      ;SFIRpj(nghl,ilpj,2);epj(ilpj,2);cpj(nkx,ilpj,2);ypj(nkx,ilpj,2)
  //      ;rpj(nkx,ilpj,2);ddepj(nqx,ilpj,nb2x);phypj(ilpj);sinphy(ilpj);
  //      exp1iphy(ilpj);exp2iphy(ilpj);exp1iphym(ilpj);exp2iphym(ilpj);

  //!-----------------------------------------
  //! PNP ARRAYS: CONF. AND GAUGE ANGLE
  //!-----------------------------------------
  // If(Allocated(exp1iphy))Deallocate(ropj,taupj,dropj,djpj,akapj,coupj,pjk  &
  //      ,SZFIpj,SFIZpj,SRFIpj,SFIRpj,epj,cpj,ypj,rpj,ddepj,phypj,sinphy  &
  //      ,exp1iphy,exp2iphy,exp1iphym,exp2iphym)
  ropj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  taupj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  dropj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  djpj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  akapj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  coupj.resize(nghl, std::vector<std::complex<double>>(ilpj));
  pjk.resize(ilpj, std::vector<std::complex<double>>(2));
  SZFIpj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  SFIZpj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  SRFIpj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  SFIRpj.resize(nghl, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  epj.resize(ilpj, std::vector<std::complex<double>>(2));
  cpj.resize(nkx, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  ypj.resize(nkx, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  rpj.resize(nkx, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(2)));
  ddepj.resize(nqx, std::vector<std::vector<std::complex<double>>>(ilpj, std::vector<std::complex<double>>(nb2x)));
  phypj.resize(ilpj);
  sinphy.resize(ilpj);
  exp1iphy.resize(ilpj);
  exp2iphy.resize(ilpj);
  exp1iphym.resize(ilpj);
  exp2iphym.resize(ilpj);

  //!-----------------------------------------
  //! FIELDS INITIALIZATION (NB! optimize)
  //!-----------------------------------------
  // ro=zero;     tau=zero;    dro=zero;    dj=zero;  aka=zero; rk=zero;
  // vn=zero;     vsn=zero;    vhbn=zero;   vrn=zero; vzn=zero; vdn=zero;
  // vp=zero;     vsp=zero;    vhbp=zero;   vrp=zero; vzp=zero; vdp=zero;
  // vSFIZn=zero; vSZFIn=zero; vSFIRn=zero; vSRFIn=zero;  vDHartree=zero;
  // vSFIZp=zero; vSZFIp=zero; vSFIRp=zero; vSRFIp=zero;

  //!-----------------------------------------
  //! FIELDS INITIALIZATION (NB! optimize)
  //!-----------------------------------------
  //  std::fill(ro.begin(), ro.end(), zero);
  //  std::fill(tau.begin(), tau.end(), zero);
  //  std::fill(dro.begin(), dro.end(), zero);
  //  std::fill(dj.begin(), dj.end(), zero);
  //  std::fill(aka.begin(), aka.end(), zero);
  //  std::fill(rk.begin(), rk.end(), zero);
  //  std::fill(rk.begin(), rk.end(), zero);
  //  std::fill(vsn.begin(), vsn.end(), zero);
  //  std::fill(vhbn.begin(), vhbn.end(), zero);
  //  std::fill(vrn.begin(), vrn.end(), zero);
  //  std::fill(vzn.begin(), vzn.end(), zero);
  //  std::fill(vdn.begin(), vdn.end(), zero);
  //  std::fill(vp.begin(), vp.end(), zero);
  //  std::fill(vsp.begin(), vsp.end(), zero);
  //  std::fill(vhbp.begin(), vhbp.end(), zero);
  //  std::fill(vrp.begin(), vrp.end(), zero);
  //  std::fill(vzp.begin(), vzp.end(), zero);
  //  std::fill(vdp.begin(), vdp.end(), zero);
  //  std::fill(vSFIZn.begin(), vSFIZn.end(), zero);
  //  std::fill(vSZFIn.begin(), vSZFIn.end(), zero);
  //  std::fill(vSFIRn.begin(), vSFIRn.end(), zero);
  //  std::fill(vSRFIn.begin(), vSRFIn.end(), zero);
  //  std::fill(vDHartree.begin(), vDHartree.end(), zero);
  //  std::fill(vSFIZp.begin(), vSFIZp.end(), zero);
  //  std::fill(vSZFIp.begin(), vSZFIp.end(), zero);
  //  std::fill(vSFIRp.begin(), vSFIRp.end(), zero);
  //  std::fill(vSRFIp.begin(), vSRFIp.end(), zero);
}

void HFBTHO_solver::gausspoints()
{

  //  !---------------------------------------------------------------------
  //  !PROGRAM DETERMINES POINTS AND WEIGHTS FOR THE GAUSS INTEGRATION.
  //  !GAUSS-LEGENDRE,LAGUERRE AND HERMITTE INTEGRATION.
  //  !---------------------------------------------------------------------
  //!  Use HFBTHO
  //  Implicit None
  double Work[1000], al, be, sparity;
  int N, N1, N2, N3, i, j, j1, j2, KINDI, kpts, nparity;
  //!
  al = 0.0;
  be = 0.0;
  kpts = 0;
  //!
  //!--------------------------------------------------------------------
  //!------------------>> Gauss-Hermite (positive nodes) <<--------------
  //!--------------------------------------------------------------------
  if (Parity)
  {

    KINDI = 4;
    N = 2 * ngh;
    N1 = 3;
    N2 = N1 + N;
    N3 = N2 + N; //! Yesp
  }
  else
  {

    KINDI = 4;
    N = ngh;
    N1 = 3;
    N2 = N1 + N;
    N3 = N2 + N; //! Nop
  }
  //
  gaussq(KINDI, N, al, be, kpts, &Work[1 - 1], &Work[N1 - 1], &Work[N2 - 1], &Work[N3 - 1]);
  // If(ierror_flag.Ne.0) Return
  //! WRITE(11,*) ' !'
  //! WRITE(11,*) ' !GAUSS-HERMITTE INTEGRATION POINTS AND WEIGHTS'
  if (Parity)
  {
    nparity = ngh; //! Yesp
    sparity = 2.0;
  }
  else
  {
    nparity = 0; //! Nop
    sparity = 1.0;
  }
  //
  for (int i = nparity + 1; i <= N; i++)
  {

    j = i - nparity;
    j1 = N2 + i - 1;
    j2 = N3 + i - 1;
    xh[j - 1] = Work[j1 - 1];
    wh[j - 1] = sparity * exp(Work[j1 - 1] * Work[j1 - 1] + log(Work[j2 - 1]));
    //! write(11,'(a,i2,a,D24.18,a,I2,a,D24.18)')  '  xh(',J,')=',xh(J),';   wh(',J,')=',wh(J)
  }
}

//!=======================================================================
void HFBTHO_solver::gaussq(const int &kindi, const int &N, const double &ALPHA, const double &BETA, const int &KPTS, double *ENDPTS,
                           double *B, double *T, double *W)
{
  //! Use HFBTHO, Only: pr,ipr
  // Implicit None
  // Integer(ipr) :: N,kindi
  double MUZERO, GAM, T1;
  int j1, J2, IERR;
  // Real(pr):: B(N),T(N),W(N),ENDPTS(2)
  // double GBSLVE;
  // !
  std::cout << "Inside gaussq, kindi: " << kindi << std::endl;
  std::cout << "Before call Class: " << MUZERO << std::endl;
  Class(kindi, N, ALPHA, BETA, B, T, MUZERO);
  std::cout << "After  call Class: " << MUZERO << std::endl;
  if (KPTS == 0)
  {
    // W = 0.0;
    W[-1 + 1] = 1.0;
    GBTQL2(N, T, B, W, IERR);
    for (int i = 0; i < N; i++)
      W[i] = MUZERO * W[i] * W[i];
    return;
  }
  if (KPTS == 2)
  {
    GAM = GBSLVE(ENDPTS[1], N, T, B);
    T1 = ((ENDPTS[1] - ENDPTS[2]) / (GBSLVE(ENDPTS[2], N, T, B) - GAM));
    B[N - 1] = sqrt(T1);
    T[N] = ENDPTS[1] + GAM * T1;
    // W = 0.0;
    W[1] = 1.0;
    GBTQL2(N, T, B, W, IERR);
    // W = MUZERO * W * W;
    for (int i = 0; i < N; i++)
      W[i] = MUZERO * W[i] * W[i];
    return;
  }
  T[-1 + N] = GBSLVE(ENDPTS[1], N, T, B) * pow(B[N - 1], 2) + ENDPTS[1];
  // W = 0.0;
  W[1] = 1.0;
  GBTQL2(N, T, B, W, IERR);
  // W = MUZERO * W * W
  for (int i = 0; i < N; i++)
    W[i] = MUZERO * W[i] * W[i];
  // End Subroutine GAUSSQ
  //!=======================================================================
  //!
}

//!=======================================================================
void HFBTHO_solver::Class(const int &kindi, const int &N, const double &ALPHA, const double &BETA, double *B, double *A, double &MUZERO)
{
  //! Use HFBTHO, Only: pr,ipr
  // Implicit None
  int i, NM1;
  // Real(pr) :: MUZERO,ALPHA,BETA,A(N),B(N)
  double ABI, DI20, AB, A2B2, FI;
  // double DGAMMA;
  const double PI = 3.1415926535897930;
  NM1 = N - 1;
  switch (kindi)
  {
  case 1:
    MUZERO = 2.0;
    for (int I = 1 - 1; i <= NM1 - 1; i++)
    // minus 1 here instead of block cause block depend on I.
    {
      A[I] = 0;
      ABI = I;
      B[I] = ABI / sqrt(4.0 * ABI * ABI - 1.0);
    }
    A[N - 1] = 0.0;
  case 2:
    MUZERO = PI;
    for (int I = 1 - 1; I <= NM1 - 1; I++)
    {
      A[I] = 0.0;
      B[I] = 0.50;
    };
    B[1] = sqrt(0.50);
    A[N] = 0.0;
  case 3:
    MUZERO = PI / 2.0;
    for (int I = 1 - 1; I <= NM1 - 1; I++)
    {
      A[I] = 0.0;
      B[I] = 0.50;
    };
    A[N] = 0.0;
  case 4:
    MUZERO = sqrt(PI);
    for (int I = 1 - 1; I <= NM1 - 1; I++)
    {
      A[I] = 0.0;
      DI20 = I / 2.0;
      B[I] = sqrt(DI20);
    }
    A[N] = 0.0;
  case 5:
    AB = ALPHA + BETA;
    ABI = 2.0 + AB;
    MUZERO = pow(2.0, (AB + 1.0)) * DGAMMA(ALPHA + 1.0) * DGAMMA(BETA + 1.0) / DGAMMA(ABI);
    A[1 - 1] = (BETA - ALPHA) / ABI;
    B[1 - 1] = sqrt(4.0 * (1.0 + ALPHA) * (1.0 + BETA) / ((ABI + 1.0) * ABI * ABI));
    A2B2 = BETA * BETA - ALPHA * ALPHA;
    for (int I = 2 - 1; I <= NM1 - 1; I++)
    {
      ABI = 2.0 * I + AB;
      A[I] = A2B2 / ((ABI - 2.0) * ABI);
      FI = I;
      B[I] = sqrt(4.0 * FI * (FI + ALPHA) * (FI + BETA) * (FI + AB) / ((ABI * ABI - 1.0) * ABI * ABI));
    }
    ABI = 2.0 * N + AB;
    A[N - 1] = A2B2 / ((ABI - 2.0) * ABI);
  case 6:
    MUZERO = DGAMMA(ALPHA + 1.0);
    for (int I = 1 - 1; I <= NM1 - 1; I++)
    {
      FI = I;
      A[I] = 2.0 * FI - 1.0 + ALPHA;
      B[I] = sqrt(FI * (FI + ALPHA));
    }
    A[N - 1] = 2.0 * N - 1.0 + ALPHA;
  default:
    break;
  }
  // End Subroutine Class
  //!=======================================================================
}

//!=======================================================================
double HFBTHO_solver::DGAMMA(double Z)
{
  //! Use HFBTHO
  //  Implicit None
  // Real(pr)::Z
  int K;
  double A[18], T, P;
  double DGAMMA;
  A[-1 + 1] = 1.00;
  A[-1 + 2] = .42278433509846780;
  A[-1 + 3] = .41184033042636720;
  A[-1 + 4] = .08157691925026090;
  A[-1 + 5] = .07424901068009040;
  A[-1 + 6] = -.00026698103334840;
  A[-1 + 7] = .01115403602403440;
  A[-1 + 8] = -.00285258214461970;
  A[-1 + 9] = .00210362870245980;
  A[-1 + 10] = -.00091848436909910;
  A[-1 + 11] = .00048742279447680;
  A[-1 + 12] = -.00023472040189190;
  A[-1 + 13] = .00011153395196660;
  A[-1 + 14] = -.00004787479838340;
  A[-1 + 15] = .00001751027271790;
  A[-1 + 16] = -.00000492037509040;
  A[-1 + 17] = .00000091991564070;
  A[-1 + 18] = -.00000008399404960;
  if (Z <= 1.00)
  {
    T = Z;
    P = A[-1 + 18];
    for (int K1 = -1 + 1; K1 <= -1 + 17; K1++)
    {
      K = 18 - K1;
      P = T * P + A[K];
    };
    DGAMMA = P / (Z * (Z + 1.00));
    return DGAMMA;
  }
  else if (Z > 1.0)
  {
    DGAMMA = P / Z;
    return DGAMMA;
  };

  if (Z <= 2.00)
  {
    T = Z - 1.00;
    P = A[-1 + 18];
    for (int K1 = -1 + 1; K1 <= -1 + 17; K1++)
    {
      K = 18 - K1;
      P = T * P + A[K];
    };
    DGAMMA = P / (Z * (Z + 1.00));
    return DGAMMA;
  }
  else if (Z > 2.0)
  {
    DGAMMA = P;
    return DGAMMA;
  }
  T = Z - 2.00;
  P = A[-1 + 18];
  for (int K1 = -1 + 1; K1 <= -1 + 17; K1++)
  {
    K = 18 - K1;
    P = T * P + A[K];
  }
  DGAMMA = P / (Z * (Z + 1.00));

  return DGAMMA;
  //!=======================================================================
  //!====================END library gauss points=========================
  //!=======================================================================
}

//!=======================================================================
void HFBTHO_solver::GBTQL2(const int &N, double *D, double *E, double *Z, int &IERR)
{
  //! Use HFBTHO
  // Implicit None
  // Integer(ipr) :: N,IERR
  // Real(pr) :: D(N),E(N),Z(N)
  int I, J, K, L, M, II, MML;
  double MACHEP, P, G, R, S, C, F, B;
  MACHEP = pow(16.0, -14);
  IERR = 0;
  if (N == 1)
    return;
  E[N] = 0.0;
  for (int L = 1; L <= N; L++)
  {
    J = 0;
    while (true)
    {
      for (int M = L; M <= N; M++)
      {
        if (M == N)
          exit;
        if (abs(E[M]) <= MACHEP * (abs(D[M]) + abs(D[M + 1])))
          exit;
        continue;
      }
      P = D[L];
      if (M == L)
        exit;
      if (J == 30)
      {
        IERR = L;
        return;
      }
      J = J + 1;
      G = (D[L + 1] - P) / (2.0 * E[L]);
      R = sqrt(G * G + 1.0);
      // G = D(M) - P + E(L) / (G + Sign(R, G));
      G = D[M] - P + E[L] / (G + pow(-1, std::signbit(G)) * R);
      S = 1.0;
      C = 1.0;
      P = 0.0;
      MML = M - L;
      for (int II = 1; II <= MML; II++)
      {
        I = M - II;
        F = S * E[I];
        B = C * E[I];
        if (abs(F) >= abs(G))
        {
          C = G / F;
          R = sqrt(C * C + 1.0);
          E[I + 1] = F * R;
          S = 1.0 / R;
          C = C * S;
        }
        else
        {
          S = F / G;
          R = sqrt(S * S + 1.0);
          E[I + 1] = G * R;
          C = 1.0 / R;
          S = S * C;
        };
        G = D[I + 1] - P;
        R = (D[I] - G) * S + 2.0 * C * B;
        P = S * R;
        D[I + 1] = G + P;
        G = C * R - B;
        F = Z[I + 1];
        Z[I + 1] = S * Z[I] + C * F;
        Z[I] = C * Z[I] - S * F;
      };
      D[L] = D[L] - P;
      E[L] = G;
      E[M] = 0.0;
    }
  }
  for (int II = 2; II <= N; II++)
  {
    I = II - 1;
    K = I;
    P = D[I];
    for (int J = II; J <= N; J++)
    {
      if (D[J] >= P)
        continue;
      K = J;
      P = D[J];
    }
    if (K == I)
      continue;
    D[K] = D[I];
    D[I] = P;
    P = Z[I];
    Z[I] = Z[K];
    Z[K] = P;
  }
  // End Subroutine GBTQL2
  //!=======================================================================
  //!
}

//!=======================================================================
double HFBTHO_solver::GBSLVE(const double &SHIFT, const int &N, double *A, double *B)
{
  //! Use HFBTHO
  //  Implicit None
  // double GBSLVE;
  int NM1, i;
  double ALPHA;
  double GBSLVE;
  ALPHA = A[1] - SHIFT;
  NM1 = N - 1;
  for (int I = 2; I <= NM1; I++)
  {
    ALPHA = A[I] - SHIFT - pow(B[I - 1], 2) / ALPHA;
  };
  GBSLVE = 1.0 / ALPHA;
  // End Function GBSLVE
  //!=======================================================================
  //!
  return GBSLVE;
}
