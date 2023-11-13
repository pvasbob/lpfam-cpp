#include <iostream>
#include <cstdlib>
#include <cstring>
#include <iomanip>

#include "HFBTHO_solver.h"

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
