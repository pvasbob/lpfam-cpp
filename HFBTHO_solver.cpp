#include <iostream>
#include <cstdlib>
#include <cstring>
#include <iomanip>

#include <algorithm>
#include <bits/stdc++.h>
#include <ctime>

#include "HFBTHO_solver.h"
#include <lapack.h>
#include <cblas.h>

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
    skyrme = skyrme_INI;
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
      inout(1);
      // If(ierror_flag.Ne.0) Return
      //!-------------------------------------------------------------
      //! Preliminary constrained calculations
      //!-------------------------------------------------------------
      if (inin > 0 && icstr == 0)
      {
        epsi0 = epsi; //! remember accuracu
        icstr = 1;    //! constraint true
        epsi = 1.0;   //! small accuracu
        // Do iw=lout,lfile
        //    If(Parity) Then
        //       Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, reflection symmetry used)'
        //    Else
        //       Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, no reflection symmetry used)'
        //    Endif
        // Enddo
        iter(true); //! small constraint iterations
        // If(ierror_flag.Ne.0) Return
        icstr = 0;    //! requested unconstraint calculation
        epsi = epsi0; //! return to requested requested accuracy
      }
      //!-------------------------------------------------------------
      //! REGULAR HFB+HO ITERATIONS
      //!-------------------------------------------------------------
      // Do iw=lout,lfile
      //    If(Parity) Then
      //       Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (reflection symmetry imposed)'
      //    Else
      //       Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (no reflection symmetry imposed)'
      //    Endif
      // Enddo
      iter(true);
      // If(ierror_flag.Ne.0) Return
      // resu(1);
      // If(ierror_flag.Ne.0) Return
    }
    //

    //
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
  base(lpr);
  gaupol(lpr);
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
  hfb.resize(ndx2, std::vector<double>(ndx2));
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
  LAPACK_dsyevd("V", "L", &ndx2, &hfb[0][0], &ndx2, &evvk[0], &alwork[0], &a, &lwork[0], &b, &ier);
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

  // for (int qqi = 0; qqi < ngh; qqi++)
  //{
  //   std::cout << xh[qqi] << "** " << wh[qqi] << std::endl;
  // }

  //!--------------------------------------------------------------------
  //!---------------------------->> Gauss-Laguerre <<--------------------|
  //!--------------------------------------------------------------------
  KINDI = 6;
  // std::cout << "gauss_ngl " << ngl << std::endl;
  N = ngl;
  N1 = 3;
  N2 = N1 + N;
  N3 = N2 + N;
  gaussq(KINDI, N, al, be, kpts, &Work[-1 + 1], &Work[-1 + N1], &Work[-1 + N2], &Work[-1 + N3]);
  // if (ierror_flag.Ne .0)
  //   return;
  //! WRITE(11,*) ' !'
  //! WRITE(11,*) ' !GAUSS-LAGUERRE INTEGRATION POINTS AND WEIGHTS'
  for (int J = 1; J <= ngl; J++)
  {
    j1 = N2 + J - 1;
    j2 = N3 + J - 1;
    xl[-1 + J] = Work[-1 + j1];
    wl[-1 + J] = exp(Work[-1 + j1] + log(Work[-1 + j2]));
    sxl[-1 + J] = sqrt(xl[-1 + J]);
    //! write(11,'(a,i2,a,D24.18,a,I2,a,D24.18)')  '  xl(',J,')=',xl(J),';   wl(',J,')=',wl(J)
  }

  //!--------------------------------------------------------------------
  //!----------------->> Gauss-Legendre (positive nodes) <<--------------
  //!--------------------------------------------------------------------
  KINDI = 1;
  N = 2 * nleg;
  N1 = 3;
  N2 = N1 + N;
  N3 = N2 + N;
  gaussq(KINDI, N, al, be, kpts, &Work[-1 + 1], &Work[-1 + N1], &Work[-1 + N2], &Work[-1 + N3]);
  // If(ierror_flag.Ne.0) Return
  //! WRITE(11,*) ' !'
  //! WRITE(11,*) ' !GAUSS-LEGENDRE INTEGRATION POINTS AND WEIGHTS'
  for (int J = 1; J <= nleg; J++)
  {
    j1 = N2 + nleg + J - 1;
    j2 = N3 + nleg + J - 1;
    xleg[J] = Work[j1];
    wleg[J] = Work[j2];
    //! write(11,'(a,I2,a,D24.18,a,I2,a,D24.18)') ' xleg(',J,')=',xleg(J),'; wleg(',J,')=',wleg(J)
  }
};

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
  // std::cout << "kindi: " << kindi << std::endl;
  // std::cout << "Before call Class: " << MUZERO << std::endl;
  Class(kindi, N, ALPHA, BETA, B, T, MUZERO);
  // std::cout << "After  call Class: " << MUZERO << std::endl;

  // std::cout << "B: " << B << std::endl;
  // std::cout << "*B: " << *B << std::endl;
  // std::cout << "*(B+1): " << *(B + 1) << std::endl;
  // for (int qqi = 0; qqi < N; qqi++)
  //  std::cout << B[qqi] << " ** " << T[qqi] << std::endl;

  if (KPTS == 0)
  {
    // W = 0.0;
    for (int qqi = 0; qqi < N; qqi++)
      W[qqi] = 0.0;

    W[-1 + 1] = 1.0;

    // std::cout << std::endl;
    // for (int qqi = 0; qqi < N; qqi++)
    //   std::cout << "Before GBTQL2, T[qqi], B[qqi], W[qqi]: "
    //             << T[qqi] << "**" << B[qqi] << "**" << W[qqi] << "**" << std::endl;

    GBTQL2(N, T, B, W, IERR);

    // for (int qqi = 0; qqi < N; qqi++)
    //   std::cout << "After  GBTQL2, T[qqi], B[qqi], W[qqi]: "
    //             << T[qqi] << "**" << B[qqi] << "**" << W[qqi] << "**" << std::endl;

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
    for (int qqi = 0; qqi < N; qqi++)
      W[qqi] = 0.0;
    W[1] = 1.0;
    GBTQL2(N, T, B, W, IERR);
    // W = MUZERO * W * W;
    for (int i = 0; i < N; i++)
      W[i] = MUZERO * W[i] * W[i];
    return;
  }
  // std::cout << "Ever executed" << std::endl;
  T[-1 + N] = GBSLVE(ENDPTS[1], N, T, B) * pow(B[N - 1], 2) + ENDPTS[1];
  // W = 0.0;
  for (int qqi = 0; qqi < N; qqi++)
    W[qqi] = 0.0;
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
    for (int I = 1; I <= NM1; I++)
    // minus 1 here instead of block cause block depend on I.
    {
      A[-1 + I] = 0;
      ABI = I;
      B[-1 + I] = ABI / sqrt(4.0 * ABI * ABI - 1.0);

      // std::cout << "Inside Class, kindi, B(I), A(I): " << kindi << "*  " << B[-1 + I] << "*  " << A[-1 + I] << std::endl;
    }
    A[N - 1] = 0.0;
    break;
  case 2:
    MUZERO = PI;
    for (int I = 1; I <= NM1; I++)
    {
      A[-1 + I] = 0.0;
      B[-1 + I] = 0.50;
    };
    B[-1 + 1] = sqrt(0.50);
    A[-1 + N] = 0.0;
    break;
  case 3:
    MUZERO = PI / 2.0;
    for (int I = 1; I <= NM1; I++)
    {
      A[-1 + I] = 0.0;
      B[-1 + I] = 0.50;
    };
    A[-1 + N] = 0.0;
    break;
  case 4:
    MUZERO = sqrt(PI);
    for (int I = 1; I <= NM1; I++)
    {
      A[-1 + I] = 0.0;
      DI20 = I / 2.0;
      B[-1 + I] = sqrt(DI20);
    }
    A[-1 + N] = 0.0;

    // std::cout << "Inside Class: " << std::endl;
    // std::cout << "B: " << B << std::endl;
    // std::cout << "&B[0]: " << &B[0] << std::endl;
    // std::cout << "*B: " << *B << std::endl;

    break;
  case 5:
    AB = ALPHA + BETA;
    ABI = 2.0 + AB;
    MUZERO = pow(2.0, (AB + 1.0)) * DGAMMA(ALPHA + 1.0) * DGAMMA(BETA + 1.0) / DGAMMA(ABI);
    A[1 - 1] = (BETA - ALPHA) / ABI;
    B[1 - 1] = sqrt(4.0 * (1.0 + ALPHA) * (1.0 + BETA) / ((ABI + 1.0) * ABI * ABI));
    A2B2 = BETA * BETA - ALPHA * ALPHA;
    for (int I = 2; I <= NM1; I++)
    {
      ABI = 2.0 * I + AB;
      A[-1 + I] = A2B2 / ((ABI - 2.0) * ABI);
      FI = I;
      B[-1 + I] = sqrt(4.0 * FI * (FI + ALPHA) * (FI + BETA) * (FI + AB) / ((ABI * ABI - 1.0) * ABI * ABI));
    }
    ABI = 2.0 * N + AB;
    A[-1 + N] = A2B2 / ((ABI - 2.0) * ABI);
    break;
  case 6:
    // std::cout << "kindi, MUZERO, before: " << kindi << "* " << MUZERO << std::endl;
    MUZERO = DGAMMA(ALPHA + 1.0);
    // std::cout << "kindi, MUZERO, after : " << kindi << "* " << MUZERO << std::endl;
    for (int I = 1; I <= NM1; I++)
    {
      FI = I;
      A[-1 + I] = 2.0 * FI - 1.0 + ALPHA;
      B[-1 + I] = sqrt(FI * (FI + ALPHA));
    }
    A[-1 + N] = 2.0 * N - 1.0 + ALPHA;
    break;
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
    for (int K1 = 1; K1 <= 17; K1++)
    {
      K = 18 - K1;
      P = T * P + A[-1 + K];
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
    for (int K1 = 1; K1 <= 17; K1++)
    {
      K = 18 - K1;
      P = T * P + A[-1 + K];
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
  for (int K1 = 1; K1 <= 17; K1++)
  {
    K = 18 - K1;
    P = T * P + A[-1 + K];
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
  // std::cout << "MACHEP: " << MACHEP << std::endl;
  IERR = 0;
  if (N == 1)
    return;

  // std::cout << "E[N-1]: " << E[N - 1] << std::endl;
  // std::cout << "E[N-2]: " << E[N - 2] << std::endl;
  E[N - 1] = 0.0;

  // for (int qqi = 0; qqi < N; qqi++)
  //{
  //   std::cout << "qqi: " << qqi << "** " << D[qqi] << "** " << E[qqi] << "** " << Z[qqi] << std::endl;
  // }

  for (int L = 1; L <= N; L++)
  {
    J = 0;
    while (true)
    {
      for (M = L; M <= N; M++)
      {
        if (M == N)
          break;
        if (abs(E[-1 + M]) <= MACHEP * (abs(D[-1 + M]) + abs(D[-1 + M + 1])))
          break;
        continue;
      }
      // std::cout << L << "*  " << J << "*  " << M << std::endl;
      P = D[-1 + L];
      if (M == L)
        break;
      if (J == 30)
      {
        IERR = L;
        return;
      }
      J = J + 1;
      G = (D[-1 + L + 1] - P) / (2.0 * E[-1 + L]);
      R = sqrt(G * G + 1.0);

      // G = D(M) - P + E(L) / (G + Sign(R, G));

      // std::cout << "M: " << M << std::endl;
      G = D[-1 + M] - P + E[-1 + L] / (G + copysign(R, G));

      S = 1.0;
      C = 1.0;
      P = 0.0;
      MML = M - L;

      for (II = 1; II <= MML; II++)
      {
        I = M - II;
        F = S * E[-1 + I];
        B = C * E[-1 + I];
        if (abs(F) >= abs(G))
        {
          C = G / F;
          R = sqrt(C * C + 1.0);
          E[-1 + I + 1] = F * R;
          S = 1.0 / R;
          C = C * S;
        }
        else
        {
          S = F / G;
          R = sqrt(S * S + 1.0);
          E[-1 + I + 1] = G * R;
          C = 1.0 / R;
          S = S * C;
        };
        G = D[-1 + I + 1] - P;
        R = (D[-1 + I] - G) * S + 2.0 * C * B;
        P = S * R;
        D[-1 + I + 1] = G + P;
        G = C * R - B;
        F = Z[-1 + I + 1];
        Z[-1 + I + 1] = S * Z[-1 + I] + C * F;
        Z[-1 + I] = C * Z[-1 + I] - S * F;
      };
      D[-1 + L] = D[-1 + L] - P;
      E[-1 + L] = G;
      E[-1 + M] = 0.0;
    }
  }
  for (II = 2; II <= N; II++)
  {
    I = II - 1;
    K = I;
    P = D[-1 + I];
    for (J = II; J <= N; J++)
    {
      if (D[-1 + J] >= P)
        continue;
      K = J;
      P = D[-1 + J];
    }
    if (K == I)
      continue;
    D[-1 + K] = D[-1 + I];
    D[-1 + I] = P;
    P = Z[-1 + I];
    Z[-1 + I] = Z[-1 + K];
    Z[-1 + K] = P;
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

void HFBTHO_solver::base(bool lpr)
{
  //!---------------------------------------------------------------------
  //! set HO basis configurations in cylindrical coordinates
  //!---------------------------------------------------------------------
  // Use HFBTHO
  // Implicit None
  // Logical :: lpr
  int nze, nre, ib, ilauf, jlauf, nom, nnm, k, la, le, ip, ir, iz, il, is, nn, ND, IBX, N1, N2, iw;
  double hbz, hbp, ee;
  //!
  hbz = two * hbzero / pow(bz, 2);
  hbp = two * hbzero / pow(bp, 2);
  //!
  nze = n00max;
  nre = n00max / 2;
  ib = 0;
  ilauf = 0;
  nzm = 0;
  nrm = 0;
  nlm = 0;
  nom = 0;
  nnm = 0;
  //!-----------------------------------------------
  //! loop over k-quantum number
  //!-----------------------------------------------

  for (k = 1; k <= n00max + 1; k++)
  {
    la = k - 1;
    le = std::min(n00max, k);
    //! loop over parity
    if (!Parity)
      jlauf = ilauf; //! Nop
    for (ip = 1; ip <= 2; ip++)
    {
      if (Parity)
        jlauf = ilauf; //! Yesp
      for (ir = 0; ir <= nre; ir++)
      {
        for (iz = 0; iz <= nze; iz++)
        {
          for (il = la; il <= le; il++)
          {
            for (is = +1; is >= -1; is = is - 2)
            {
              if ((iz + 2 * ir + il) > n00max)
                continue;
              if (il + (is + 1) / 2 != k)
                continue;
              if ((iz + il) % 2 != ip - 1)
                continue;
              ee = hbz * (iz + half) + hbp * (two * ir + il + one);
              // std::cout << k << ip << ir << iz << il << is << "ee: " << ee << std::endl;
              // std::cout << k << ip << ir << iz << il << is  << "ee: " << ee << std::endl;
              if (ee <= EBASECUT)
              {
                // std::cout << "ee, EBASECUT: " << ee << "   " << EBASECUT << std::endl;
                ilauf = ilauf + 1;
                // if (ilauf.Gt.ntx) {
                //    ierror_flag=ierror_flag+1
                //    ierror_info(ierror_flag)='STOP: in base: ntx too small'
                //    Return
                // }
                nz[-1 + ilauf] = iz;
                nr[-1 + ilauf] = ir;
                nl[-1 + ilauf] = il;
                ns[-1 + ilauf] = is;
                npar[-1 + ilauf] = ip;
                nn = iz + 2 * ir + il;

                // std::cout << "nz: "
                //           << nz[-1 + ilauf]
                //           << nr[-1 + ilauf]
                //           << nl[-1 + ilauf]
                //           << ns[-1 + ilauf]
                //           << npar[-1 + ilauf]
                //           << nn
                //           << std::endl;
                //  Write(tb(ilauf),100) 2*k-1,tp(ip),nn,iz,il
                //  100                    //Format(i2,a1,'[',i2,',',i2,',',i2,']')
                //  Do iw=lout,lfile
                //     If(lpr.And.IDEBUG.Gt.10) Write(iw,'(i4,a,i2,a,i2,a,i2,a,i2,a,i2,a,2x,a,1x,a,f14.8)')  &
                //          ilauf,'   nn=',nn,'   nz=',iz,'   nr=',ir,  &
                //          '   ml=',il,'  ms=',is,' /2',tb(ilauf),'e=',ee
                //  Enddo
                nzm = std::max(nzm, iz);
                nrm = std::max(nrm, ir);
                nlm = std::max(nlm, il);
                nnm = std::max(nnm, iz + 2 * ir + il);
              }
            }
          }
        }
      }
      // std::cout << nz.size() << nr.size() << nl.size() << ns.size() << npar.size() << std::endl;
      //!-----------------------------------------------
      //! Block memory
      //!-----------------------------------------------
      if (Parity)
      { //! Yesp
        if (ilauf > jlauf)
        {
          ib = ib + 1;
          ia[-1 + ib] = jlauf;
          id[-1 + ib] = ilauf - jlauf;
          ikb[-1 + ib] = k;
          ipb[-1 + ib] = ip;
          std::cout << ilauf
                    << " " << jlauf
                    << " " << ib
                    << " " << ikb[-1 + ib]
                    << " " << ipb[-1 + ib]
                    << " " << std::endl;
          // Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
          //! ir=(ib+1)/2
          //! write(*,*)  ib,2*k-1,'2*Omega=',2*ir - 1
          // Do iw=lout,lfile
          //    If(lpr.And.IDEBUG.Gt.10) Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
          // End Do
        }
        // if(id(ib) == 0) {
        //    ierror_flag=ierror_flag+1
        //    ierror_info(ierror_flag)='STOP: in base Block Memory(1)'
        //    Return
        // }
      }
    }
    //!-----------------------------------------------
    //! Block memory
    //!-----------------------------------------------
    if (!Parity)
    { //! Nop
      if (ilauf > jlauf)
      {
        ib = ib + 1;
        ia[-1 + ib] = jlauf;
        id[-1 + ib] = ilauf - jlauf;
        ikb[-1 + ib] = k;
        ipb[-1 + ib] = ip;
        // Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
        // Write(txb(ib),'(i3,a,i2,a)') ib,'. block:  k=',k+k-1,'/2'
        // Do iw=lout,lfile
        //    If(lpr.And.IDEBUG.Gt.10) Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
        // End Do
      }
      // If(id(ib).Eq.0) Then
      //    ierror_flag=ierror_flag+1
      //    ierror_info(ierror_flag)='STOP: in base Block Memory(2)'
      //    Return
      // Endif
    }
  }
  nb = ib;
  nt = ilauf;
  //!-----------------------------------------------
  //! broyden/linear mixing (storage)
  //!-----------------------------------------------
  nhhdim = 0;
  for (ib = 1; ib <= nb; ib++)
  {
    ND = id[-1 + ib];
    for (N1 = 1; N1 <= ND; N1++)
    {
      for (N2 = 1; N2 <= N1; N2++)
      {
        nhhdim = nhhdim + 1;
      }
    }
  }
  nhhdim2 = 2 * nhhdim;
  nhhdim3 = 3 * nhhdim;
  nhhdim4 = 4 * nhhdim;
  std::cout << "nhhdimx: "
            << nhhdim2
            << nhhdim3
            << nhhdim4
            << std::endl;
  // If(Allocated(brin)) Deallocate(brin,brout)
  // Allocate(brin(nhhdim4),brout(nhhdim4))
  brin.resize(nhhdim4), brout.resize(nhhdim4);
  //!-----------------------------------------------
  //! Print statistics
  //!-----------------------------------------------
  // If(lpr) Then
  //    Do iw=lout,lfile
  //       Write(iw,'(a,i4)')   '  actual basis used:'
  //       Write(iw,'(a,i4)')   '  number of blocks: nb=',nb
  //       Write(iw,'(a,i4)')   '  number of levels: nt=',nt
  //       Write(iw,'(a,i4)')   '  maximal 2*omega : nom=',nom
  //       Write(iw,'(a,i4)')   '  maximal nz:       nzm=',nzm
  //       Write(iw,'(a,i4)')   '  maximal nr:       nrm=',nrm
  //       Write(iw,'(a,i4)')   '  maximal ml:       nlm=',nlm
  //       Write(iw,'(a,i4)')   '  maximal nn=nz+2*nr+nl=',nnm
  //       Write(iw,'(a,i4)')   '  2 x bigest block dim.=',ndx2
  //       Write(iw,'(a,i8)')   '  Nonzero elemets of HH=',nhhdim
  //       Write(iw,'(a,i8)')   '  Number Broyden elemens=',nhhdim4
  //       Write(iw,'(a,i4)')
  //    Enddo
  // Endif
  // If(nzm.Ge.n00max.Or.(nom-1)/2.Eq.n00max) Then
  //    Write(*,*) 'nzm=',nzm,'  (nom-1)/2=',(nom-1)/2,'  n00max=',n00max
  //    ierror_flag=ierror_flag+1
  //    ierror_info(ierror_flag)='STOP: Please increase n00max to have correct basis'
  // Endif
}

void HFBTHO_solver::gaupol(bool lpr)
{
  //!---------------------------------------------------------------------
  //! HO wave functions in cylindrical coordinates
  //!---------------------------------------------------------------------
  // Use HFBTHO
  // Implicit None
  // Logical :: lpr
  double w0, z, x, s, s0, s1, w00, w4pii, dsq, d1, d2, d3, d4, hs0, hs1;
  int ih, il, iw, ix, n, l, n1, n2;
  //!-----------------------------------------------
  //! Hermit
  //!-----------------------------------------------
  w4pii = pow(PI, (-0.250));
  for (ih = 1; ih <= ngh; ih++)
  {
    z = xh[-1 + ih];
    w0 = w4pii * exp(-half * z * z);
    // std::cout << "z, w0: " << -1 + ih << "   " << z << w0 << std::endl;
    //! functions qh, qh1; norm: \sum_{ih} qh(n1,ih)*qh(n2,ih)=\delta_{n1,n2}
    w0 = w0 * sqrt(wh[-1 + ih]);
    qh[0][-1 + ih] = w0;
    // std::cout << "sq[-1+2]" << sq[-1 + 2] << std::endl;
    qh[1][-1 + ih] = sq[2] * w0 * z;
    qh1[0][-1 + ih] = -w0 * z;
    qh1[1][-1 + ih] = sq[2] * w0 * (one - z * z);

    // std::cout << "z: " << z << "w0: " << w0 << "qh"
    //           << "1" << -1 + ih << "  " << qh[1][-1 + ih] << std::endl;

    for (n = 2; n <= nzm; n++)
    {
      qh[n][-1 + ih] = sqi[n] * (sq[2] * z * qh[n - 1][-1 + ih] - sq[n - 1] * qh[n - 2][-1 + ih]);
      qh1[n][-1 + ih] = sq[n + n] * qh[n - 1][-1 + ih] - z * qh[n][-1 + ih];
    }
  };

  // for (int i = 0; i < sq.size(); i++)
  //{
  //   std::cout << sq.size() << "sq[i]: " << sq[i] << std::endl;
  // }

  //!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
  // std::cout << "dimension: " << qh.size() << " " << qh[0].size() << std::endl;
  // std::cout << "xh size: " << xh.size() << std::endl;
  // for (int i = 0; i < xh.size(); i++)
  //{
  //  std::cout << i << " " << xh[i] << std::endl;
  //}

  // for (int i = 0; i < qh.size(); i++)
  //{
  //   for (int j = 0; j < qh[0].size() - 1; j++)
  //   {
  //     std::cout << qh.size() << " " << qh[0].size() << " " << i << " " << j << " " << qh[i][j] << std::endl;
  //   }
  // }

  //!-----------------------------------------------
  //! Laguerre
  //!-----------------------------------------------
  for (il = 1; il <= ngl; il++)
  {
    x = xl[-1 + il];
    w00 = sq[2] * exp(-half * x);
    // std::cout << "-1+il"
    //           << " " << -1 + il << " " << x << " " << w00 << std::endl;
    for (l = 0; l <= nlm; l++)
    {
      //! functions ql, ql1; norm: \sum_{il} ql(n1,l,il)*ql(n2,l,il)=\delta_{n1,n2}
      w0 = w00 * sqrt(half * wl[-1 + il] * pow(x, l));
      //
      // std::cout << "-1+il, l" << -1 + il << " " << l << " " << w0 << std::endl;
      //
      ql[0][l][-1 + il] = wfi[l] * w0;
      ql[1][l][-1 + il] = (l + 1 - x) * wfi[l + 1] * w0;
      ql1[0][l][-1 + il] = (l - x) * wfi[l] * w0;
      ql1[1][l][-1 + il] = ((l * l + l) - x * (l + l + 3) + x * x) * wfi[l + 1] * w0;
      // std::cout << "l, -1+il" << l << " " << -1 + il << " " << ql1[0][l][-1 + il] << "  " << ql1[1][l][-1 + il] << std::endl;
      for (n = 2; n <= nrm; n++)
      {
        dsq = sq[n] * sq[n + l];
        d1 = (n + n + l - 1) - x;
        d2 = sq[n - 1] * sq[n - 1 + l];
        d3 = n + n + l - x;
        d4 = two * dsq;
        ql[n][l][-1 + il] = (d1 * ql[n - 1][l][-1 + il] - d2 * ql[n - 2][l][-1 + il]) / dsq;
        ql1[n][l][-1 + il] = d3 * ql[n][l][-1 + il] - d4 * ql[n - 1][l][-1 + il];
        // std::cout << "n, l, -1+il:    " << n << " " << l << " " << -1 + il << " " << ql[n][l][-1 + il] << "  " << ql1[n][l][-1 + il] << std::endl;
      }
    }
  }

  // for (int i = 0; i < ql.size(); i++)
  //{
  //   for (int j = 0; j < ql[0].size(); j++)
  //   {
  //     for (int k = 0; k < ql[0][0].size() - 1; k++)
  //     {
  //       std::cout << "i, j, k: " << i << " " << j << " " << k << " " << ql[i][j][k] << "  " << ql1[i][j][k] << std::endl;
  //     }
  //   }
  // }
  //!
  //!-----------------------------------------------
  //! test accuracy for hermit orthonormalization
  //!-----------------------------------------------
  hs0 = 0.0;
  hs1 = 2.0;
  for (n1 = 0; n1 <= nzm; n1++)
  {
    for (n2 = 0; n2 <= n1; n2++)
    {
      if ((n1 - n2) % 2 == 0)
      {
        s = zero;
        for (ih = 1; ih <= ngh; ih++)
        {
          s = s + qh[n1][-1 + ih] * qh[n2][-1 + ih];
        }
        if (n1 != n2)
          hs0 = std::max(s, hs0);
        else
          hs1 = std::min(s, hs1);
      }
    }
  }
  //!-----------------------------------------------
  //! test accuracy for laguerre orthonormalization
  //!-----------------------------------------------
  s0 = 0.0;
  s1 = 2.0;
  for (l = 0; l <= nlm; l++)
  {
    for (n1 = 0; n1 <= nrm; n1++)
    {
      for (n2 = 0; n2 <= n1; n2++)
      {
        s = zero;
        for (il = 1; il <= ngl; il++)
        {
          s = s + ql[n1][l][-1 + il] * ql[n2][l][-1 + il];
        }
        if (n1 != n2)
          s0 = std::max(s, s0);
        else
          s1 = std::min(s, s1);
      }
    }
  }

  std::cout << "hs0, hs1: " << hs0 << "  " << hs1 << std::endl;
  std::cout << "s0, s1: " << s0 << "  " << s1 << std::endl;
  //
  coordinateLST(false);
}

void HFBTHO_solver::coordinateLST(bool lpr)
{
  //  !------------------------------------------------------------------
  //  ! HO/THO
  //  !------------------------------------------------------------------
  //  Use HFBTHO
  // #ifndef hide_tho
  //  Use THO_MODULE, Only: f01234
  // #endif
  //  Implicit None
  //  Logical :: lpr
  int i, il, ih;

  std::cout << "iLST1: " << iLST1 << std::endl;

  if (iLST1 == 0)
  {
    //! HO-basis
    for (il = 1; il <= ngl; il++)
    {
      for (ih = 1; ih <= ngh; ih++)
      {
        i = ih + (il - 1) * ngh;
        fh[-1 + i] = bz * xh[-1 + ih];
        fl[-1 + i] = bp * sqrt(xl[-1 + il]);
        wdcor[-1 + i] = PI * wh[-1 + ih] * wl[-1 + il] * bz * bp * bp;
        wdcori[-1 + i] = one / wdcor[-1 + i];
        //
        // std::cout << "il, ih: " << il << " " << ih << " " << fh[-1 + i] << " " << fl[-1 + i] << std::endl;
      }
    }
  }
  // #ifndef hide_tho
  else
  {
    //! THO basis
    // Call f01234(.False.)
    // If(ierror_flag.Ne.0) Return
    // #endif
  }            //!
  optHFBTHO(); //! optimal HO/THO combinations
  //  If(ierror_flag.Ne.0) Return
  //!
}

void HFBTHO_solver::optHFBTHO()
{
  // Subroutine optHFBTHO
  //   !---------------------------------------------------------------------
  //   ! oprimization arrays
  //   ! NB FI2D_opt(JA,ihil) == Laplassisan(r,z) HOwf
  //   !    FID2D-xlamy2*FID  == Laplassisan(r,z,phy) FID
  //   !    FIU2D-xlapy2*FIU  == Laplassisan(r,z,phy) FIU
  //   !---------------------------------------------------------------------
  //   Use HFBTHO
  //   Implicit None
  int i, ih, il, ib, ibx, nd, nza, nra, nla, nsa;
  int ihil, laplus, im, JA, N1, N2, ndnd, n12, n21;
  double qla, v2, v4, yi, y, y2, qha, qhla, xmi, u, u2, un, up, xxx;
  double sml2, cnzaa, cnraa, a, b;
  // double FITW1,FITW2,FITW3,FITW4;
  double fitw1, fitw2, fitw3, fitw4;
  // double fi1r,fi1z,fi2d,QHL1A,QH1LA,vh,vdh,vsh,hbh;
  double fi1r, fi1z, fi2d, qhl1a, qh1la, vh, vdh, vsh, hbh;
  // double SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh;
  double srfih, sfirh, sfizh, szfih, snablarh, snablazh;
  // double xlam, xlam2, xlamy, xlamy2, xlap, xlap2, xlapy, xlapy2, XLAMPY;
  double xlam, xlam2, xlamy, xlamy2, xlap, xlap2, xlapy, xlapy2, xlampy;
  double bpi, bpi2, bzi, bzi2, xh2;
  //!
  bpi = one / bp;
  bpi2 = bpi * bpi;
  bzi = one / bz;
  bzi2 = bzi * bzi;
  //!
  //!-----------------------------------------
  //! Allocate the optimization arrays
  //!-----------------------------------------
  // If(Allocated(QHLA_opt)) Deallocate(QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt)
  QHLA_opt.resize(ntx, std::vector<double>(nghl));
  FI1R_opt.resize(ntx, std::vector<double>(nghl));
  FI1Z_opt.resize(ntx, std::vector<double>(nghl));
  FI2D_opt.resize(ntx, std::vector<double>(nghl));
  y_opt.resize(nghl);
  //!----------------------------------------------
  //! START BLOCKS
  //!----------------------------------------------
  for (ib = 1; ib <= nb; ib++)
  {
    nd = id[-1 + ib];
    im = ia[-1 + ib];
    if (Parity)
      laplus = (ib + 1) / 2; //! Yesp
    else
      laplus = ib; //! Nop
    // Endif
    xlap = laplus;
    xlam = xlap - one;
    xlap2 = xlap * xlap;
    xlam2 = xlam * xlam;
    //!----------------------------------------------
    //! SUM OVER GAUSS INTEGRATION POINTS
    //!----------------------------------------------
    for (il = 1; il <= ngl; il++)
    {
      v2 = half / xl[-1 + il];
      v4 = v2 * v2;
      for (ih = 1; ih <= ngh; ih++)
      {
        ihil = ih + (il - 1) * ngh;
        xh2 = pow(xh[-1 + ih], 2);
        if (iLST1 == 0)
        {
          //! HO-basis
          yi = sqrt(xl[-1 + il]) * bp;
          y = one / yi;
          y2 = y * y;
          xlamy = xlam * y;
          xlamy2 = xlam2 * y2;
          xlapy = xlap * y;
          xlapy2 = xlap2 * y2;
          // XLAMPY = XLAMY + XLAPY;
          xlampy = xlamy + xlapy;
        }
        // #ifndef hide_tho
        else
        {
          ////! THO - basis
          // y = fli(ihil);
          // y2 = y * y;
          // xlamy = xlam * y;
          // u = xh(ih);
          // u2 = u * u;
          // xlamy2 = xlam2 * y2;
          // xlapy = xlap * y;
          // xlapy2 = xlap2 * y2;
          // XLAMPY = XLAMY + XLAPY;
        }
        // #endif
        // End If
        y_opt[-1 + ihil] = y;
        //!----------------------------------------------
        //! SCAN OVER BASIS STATES
        //!----------------------------------------------
        for (N1 = 1; N1 <= nd; N1++)
        {
          JA = N1 + im;
          nla = nl[-1 + JA];
          nra = nr[-1 + JA];
          nza = nz[-1 + JA];
          nsa = ns[-1 + JA];
          // SML2 = NLA * NLA;
          sml2 = nla * nla;
          // CNZAA = NZA + NZA + 1;
          cnzaa = nza + nza + 1;
          // CNRAA = NRA + NRA + NLA + 1;
          cnraa = nra + nra + nla + 1;
          qha = qh[nza][-1 + ih];
          qla = ql[nra][nla][-1 + il];
          qhla = qha * qla;
          qhl1a = qha * ql1[nra][nla][-1 + il] * v2;
          qh1la = qh1[nza][-1 + ih] * qla;
          if (iLST1 == 0)
          {
            //! HO-basis
            fi1r = (two * sqrt(xl[-1 + il]) * bpi) * qhl1a;
            fi1z = bzi * qh1la;
            fi2d = ((xh2 - cnzaa) * bzi2 + four * (p14 - cnraa * v2 + sml2 * v4) * xl[-1 + il] * bpi2) * qhla;
          }
          // #ifndef hide_tho
          else
          {
            //  //! THO-basis
            //  u = xh(ih);
            //  u2 = u * u;
            //  FI1R = FP4(IHIL) * QHLA + FP5(IHIL) * QH1LA + FP6(IHIL) * QHL1A;
            //  FI1Z = FP1(IHIL) * QHLA + FP2(IHIL) * QH1LA + FP3(IHIL) * QHL1A;
            //  FI2D = (FS1(IHIL) * QH1LA * QH1LA + FS2(IHIL) * QHL1A * QHL1A  \ 
          //            +
            //          FOUR * FS4(IHIL) * QH1LA * QHL1A + TWO * (FS5(IHIL) * QH1LA + FS6(IHIL) * QHL1A) * QHLA \ 
          //            +
            //          ((U2 - CNZAA) * FS1(IHIL) + (p14 - CNRAA * V2 + SML2 * V4) * FS2(IHIL) + FS3(IHIL)) * QHLA * QHLA - TWO * (FI1R * FI1R + FI1Z * FI1Z)) /
            //         (TWO * QHLA);
          }
          // #endif
          // Endif
          QHLA_opt[-1 + JA][-1 + ihil] = qhla;
          FI2D_opt[-1 + JA][-1 + ihil] = fi2d;
          FI1R_opt[-1 + JA][-1 + ihil] = fi1r;
          FI1Z_opt[-1 + JA][-1 + ihil] = fi1z;
          // std::cout << "ib: " << QHLA_opt[-1 + JA][-1 + ihil] << " " << FI2D_opt[-1 + JA][-1 + ihil] << " " << FI1R_opt[-1 + JA][-1 + ihil] << " " << FI1Z_opt[-1 + JA][-1 + ihil] << std::endl;
        } //! N1
        //!
      } //! IH
    }   //! IL
  }     //! IB
}

void HFBTHO_solver::iter(bool lpr)
{
  //!------------------------------------------------------------------
  //! Iterations through successive diagonalisation
  //!------------------------------------------------------------------
  // Use HFBTHO
  // Implicit None
  // Logical :: lpr
  double assprn, delln[2];
  clock_t time;
  int iw, it, ite;
  clock_t time1, time2, time3;
  //!---------------------------------------------------
  //! print to screen('lout')/thoout.dat('lfile')
  //!---------------------------------------------------
  // Do iw=lout,lfile
  //    If(iLST.Eq.0) Then
  //       Write(iw,'(a,f7.3,4(a,i3),a)')  &
  //            '  |HFB+HO> iterations(b0=',b0,', Nsh=',n00,  &
  //            ', inin=',inin,', N=',npr(1),', Z=',npr(2),')...'
  //    Else
  //       If(iLST1.Eq.0.Or.iasswrong(3).Ne.0) Then
  //          If(iasswrong(3).Ne.0) Then
  //             Write(iw,'(a,f7.3,a,i3,a)')  &
  //                  '  |HFB+THO substituted by HFB+HO> iterations (b0=',  &
  //                  b0,', Nsh=',n00,')...'
  //          Else
  //             Write(iw,'(a,f7.3,a)')'  towards |hfb+tho> iterations...'
  //             Write(iw,'(a,f7.3,a)')
  //             Write(iw,'(a,f7.3,a,i3,a)')  &
  //                  '  |Preliminary HFB+HO> iterations (b0=',b0,', Nsh=',n00,')...'
  //          End If
  //       Else
  //          If(itass.Eq.1) Then
  //             Write(iw,'(2(a,f7.3),a,i3,a)')  &
  //                  '  |HFB+THO> iterations(b0=',b0,', neutron density decay=',  &
  //                  decay,', Nsh=',n00,')...'
  //          Else
  //             Write(iw,'(2(a,f7.3),a,i3,a)')  &
  //                  '  |HFB+THO> iterations(b0=',b0,', proton density decay=',  &
  //                  decay,', Nsh=',n00,')...'
  //          End If
  //       End If
  //    End If
  //    Write(iw,1)
  //    Write(iw,'(20(a))')'  i','          si ','    mix ','  beta', '  &
  //         &  Etot ','    A ','      rn','      rp ','        En', '   &
  //         &    Dn','      Ep','      Dp','        Ln  ','   Lp ', '   &
  //         &    time' !Idro '
  //    Write(iw,1)
  // End Do
  //
  std::cout << "iter() execute: " << std::endl;

  for (ite = 1; ite <= maxi; ite++)
  {

    // Cpu_time(time1);
    time1 = std::clock();
    // std::cout << "ite, lpr, time1 " << ite << " " << lpr << " " << time1 << std::endl;

    //!
    iiter = ite;
    // std::cout << "lpr, iiter " << lpr << " " << iiter << std::endl;
    if (lpr && iiter == 1)
    {
      // std::cout << "Inside if. lpr, iiter " << lpr << " " << iiter << std::endl;
      assprn = ass[-1 + 1];
      if (assprn > ass[-1 + 2])
        assprn = -ass[-1 + 2]; //! protons come with '-'
                               // delLN = del;
                               // std::cout << "sizeof(del): " << sizeof(del) / sizeof(del[0]) << std::endl;
      for (int i = 0; i < sizeof(del) / sizeof(del[0]); i++)
        delln[i] = del[i];
      //
      std::cout << "kindhfb: " << kindhfb << std::endl;
      //
      if (kindhfb < 0)
      {
        // delLN = del + ala2; //! LN case
        for (int i = 0; i < sizeof(del) / sizeof(del[0]); i++)
        {
          delln[i] = del[i] + ala2[i];
          std::cout << "delln[i]: " << delln[i] << std::endl;
        }
      }
      //! during iterations print
      // Do iw=lout,lfile
      //    If(Max(Abs(drhoi(1)),Abs(drhoi(2))).Gt.1.0e-10_pr) Then
      //       Write(*,*) '  WARNING! Int(Dro)=',Max(Abs(drhoi(1)),Abs(drhoi(2)))
      //    Endif
      //    Write(iw,2) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
      //         ept(2),delLN(2),alast(1),alast(2),time
      // End Do
    }
    //!-------------------------------------------------
    //! HFBDIAG
    //!-------------------------------------------------
    for (it = itmin; it <= itmax; it++)
      hfbdiag(it, 0); //! hfb diagonalization with minimal canonical
                      // If(ierror_flag.Ne.0) Return
    //     !-------------------------------------------------
    //     ! EXPECT, DENSIT, COULOMB, FILED, GAMDEL
    //     !-------------------------------------------------
    //     Call expect(.False.)    ! expectation values
    //     If(ierror_flag.Ne.0) Return
    //     Call field              ! new fields
    //     If(ierror_flag.Ne.0) Return
    //     Call gamdel             ! hf-matrix
    //     If(ierror_flag.Ne.0) Return
    //     !-------------------------------------------------
    //     ! Dumping control (old linear mixing)
    //     !-------------------------------------------------
    //     xmix0=0.10 !original 0.1
    //     If(si.Lt.siold) Then
    //        xmix=Min(xmax,xmix * 1.130_pr);  !old value 1.13
    //     Else
    //        xmix=xmix0
    //     End If
    //     siold=si
    //     !-------------------------------------------------
    //     ! time per iteration
    //     !-------------------------------------------------
    //     Call Cpu_time(time2)
    //     time=time2-time1; time3=time3+time
    //     !-------------------------------------------------
    //     ! Solution is OK within the iteration limit
    //     !-------------------------------------------------
    //     If(iiter.Ge.2.And.si.Lt.epsi) Then
    //        If(iLST1.Eq.0) Then
    //           iError_in_HO=0
    //        Else
    //           iError_in_THO=0
    //        End If
    //        ! iteration interrupted print
    //        If(.Not.lpr) Then
    //           delLN=del; If(kindhfb.Lt.0) delLN=del+ala2
    //           Do iw=lout,lfile
    //              Write(iw,3) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
    //                   ept(2),delLN(2),alast(1),alast(2),time !Max(Abs(drhoi(1)),Abs(drhoi(2)))
    //              Write(iw,'(a,f8.3,a)') '  Total CPU time=',time3/60.0_pr,' minutes'
    //           Enddo
    //        End If
    //        ! converged print
    //        Do iw=lout,lfile
    //           Write(iw,4) iiter,si,iError_in_HO,iError_in_THO
    //           Write(iw,'(a,f8.3,a)') '  Total CPU time=',time3/60.0_pr,' minutes'
    //        Enddo
    //        iiter=iiter+1
    //        Return
    //     End If
    //     !-------------------------------------------------
    //     ! Slow convergence and lambda >0 (stop iterations)
    //     !-------------------------------------------------
    //     If(iiter.Ge.1000.And.(alast(1).Gt.zero.Or.alast(2).Gt.zero)) Exit
    //     !
    //
    //  }
  }
}

void HFBTHO_solver::hfbdiag(int it, int icanon)
{
  // std::cout << "it, icacon: " << it << " " << icancon << std::endl;
  //!------------------------------------------------------------------
  //! Skyrme-HFB diagonalization in axial HO/THO basis
  //!------------------------------------------------------------------
  // Use HFBTHO
  // Implicit None
  bool lpr_pwi, norm_to_improve;
  double al, al2, emin, hla, dla, pn, eqpe, ela, enb, enb1, ekb, s1, s2, s3, alnorm, sitest;
  int iw, i0, ibiblo, ier, i, j, k, k0, kl, lc, ib, nd;
  int nhfb, n1, n2, kaib, m, ndk, nd1, nd2, kdib, k1, k2, id1, id2;
  int n12, n21, ntz, nhhph, nhhpp, ibro, ibroib, i_uv, i_eqp;
  std::vector<double> *EqpPo;
  std::vector<double> *VqpPo;
  std::vector<double> *UqpPo;
  std::vector<int> *KpwiPo;
  std::vector<int> *KqpPo;
  //!
  // Call get_CPU_time('hfbdiag',0)
  //!
  if (it == 1)
  {
    EqpPo = &REqpN;
    VqpPo = &RVqpN;
    UqpPo = &RUqpN;
    KpwiPo = &KpwiN;
    KqpPo = &KqpN;
  }
  else
  {
    EqpPo = &REqpP;
    VqpPo = &RVqpP;
    UqpPo = &RUqpP;
    KpwiPo = &KpwiP;
    KqpPo = &KqpP;
  }
  //
  std::cout << "it, icanon: " << it << " " << icanon << std::endl;
  std::cout << "(*KpwiPo).size(): " << (*EqpPo).size() << " " << (*VqpPo).size() << " " << (*UqpPo).size() << " " << (*KpwiPo).size() << " " << KpwiN.size() << std::endl;

  for (int i = 0; i < (*KpwiPo).size(); i++)
  {
    (*KpwiPo)[i] = 0;
    (*KqpPo)[i] = 0;
  }
  //!
  nhhph = (it - 1) * nhhdim;
  nhhpp = (it + 1) * nhhdim;
  //
  std::cout << "nhhph, nhhpp, nhhdim: " << nhhph << " " << nhhpp << " " << nhhdim << std::endl;
  //!
  std::cout << "brin.size(), brout.size(): " << brin.size() << " " << brout.size() << std::endl;
  for (int qqi = 1; qqi < brin.size(); qqi++)
  {
    std::cout << "qqi, brin(qqi), brout(qqi): " << qqi << " " << brin[qqi] << " " << brout[qqi] << std::endl;
  }
  //!
  //!------------------------------------------------------------------
  //! Loop the internal normalization
  //!------------------------------------------------------------------
  sitest = std::max(std::min(0.10, si * 0.010), 0.000010);
  norm_to_improve = true;
  inner[-1 + it] = -1;
  sumnz[-1 + it] = one;

  while (norm_to_improve)
  {
    std::cout << "sumnz: " << sumnz[-1 + it] << std::endl;
    //!
    inner[-1 + it] = inner[-1 + it] + 1;
    //!
    if (abs(sumnz[-1 + it]) < sitest && inner[-1 + it] == 20)
      norm_to_improve = false;
    //!
    sumnz[-1 + it] = zero;
    v2min[-1 + it] = one;
    Dispersion[-1 + it] = zero;
    std::cout << "it, sumnz, v2min, Dispersion: " << it << " " << sumnz[-1 + it] << " " << v2min[-1 + it] << " " << Dispersion[-1 + it] << std::endl;
    //!
    kl = 0;
    emin = 1000.0;
    al = ala[-1 + it];
    std::cout << "it, al: " << it << " " << al << std::endl;
    //!
    //! blocking
    // If(iparenti(it).Eq.0) blomax(it)=0
    // blo123d(it)=0; blok1k2d(it)=0; blocanon(it)=0;
    // ibiblo=bloblo(keyblo(it),it)
    //!------------------------------------------------------------------
    //! Runs over blocks
    //!------------------------------------------------------------------
    i_uv = 0;
    i_eqp = 0;
    lc = 0;
    lcanon[0][-1 + it] = 0;
    klmax[0] = 0;
    klmax[1] = 0;
    ibro = 0;
    for (ib = 1; ib <= nb; ib++)
    {
      nd = id[-1 + ib];
      nhfb = nd + nd;
      i0 = ia[-1 + ib];
      m = ib + (it - 1) * nbx;
      ibroib = ibro;
      std::cout << "nd, nhfb, i0, m, ibroib: " << nd << " " << nhfb << " " << i0 << " " << m << " " << ibroib << std::endl;
      //!------------------------------------------------------------------
      //!  hfb-matrix
      //!------------------------------------------------------------------
      for (n1 = 1; n1 <= nd; n1++)
      {
        nd1 = n1 + nd;
        std::cout << "nd1: " << nd1 << std::endl;
        for (n2 = 1; n2 <= n1; n2++)
        {
          nd2 = n2 + nd;
          ibro = ibro + 1;
          std::cout << "nd1, n2, nd2, ibro: " << nd1 << " " << n2 << " " << nd2 << " " << ibro << std::endl;
          hla = brin[-1 + nhhph + ibro];
          dla = brin[-1 + nhhpp + ibro];
          hfb[-1 + n1][-1 + n2] = hla;
          hfb[-1 + nd2][-1 + n1] = dla;
          hfb[-1 + nd1][-1 + n2] = dla;
          hfb[-1 + nd1][-1 + nd2] = -hla;
          std::cout << "hla, dla: " << hla << " " << dla << std::endl;
        }
        hfb[-1 + n1][-1 + n1] = hfb[-1 + n1][-1 + n1] - al;
        hfb[-1 + nd1][-1 + nd1] = hfb[-1 + nd1][-1 + nd1] + al;
      }
      ier = 0;
      // DSYEVD('V', 'L', nhfb, hfb, ndx2, evvk, ALWORK, ialwork, LWORK, ilwork, ier);
      LAPACK_dsyevd("V", "L", &nhfb, &hfb[0][0], &ndx2, &evvk[0], &alwork[0], &ialwork, &lwork[0], &ilwork, &ier);
      // std::cout << "LAPACK_dsyevd call: " << std::endl;
      //  LAPACK_dsyevd("V", "L", &ndx2, &hfb[0], &ndx2, &evvk[0], &alwork[0], &a, &lwork[0], &b, &ier);
      //! Call dsyev('V','L',nhfb,hfb,ndx2,evvk,ALWORK,ialwork,ier)
      //!------------------------------------------------------------------
      //!  NB! Diagonalization bug in LAPACK
      //!------------------------------------------------------------------
      // no need back If(ier.Gt.0) Then
      // no need back    Do iw=lout,lfile
      // no need back       Write(iw,*) 'FATAL ERROR CONDITION IN HFBDIAG DSYEVD, ier=',ier,'(RECOVERED)'
      // no need back    End Do
      // no need back    ibro=ibroib
      // no need back    Do n1=1,nd
      // no need back       nd1=n1+nd
      // no need back       Do n2=1,n1
      // no need back          nd2=n2+nd; ibro=ibro+1
      // no need back          hla=brin(nhhph+ibro);     dla=brin(nhhpp+ibro)
      // no need back          hfb(n1,n2)=hla;           hfb(nd2,n1)=dla
      // no need back          hfb(nd1,n2)=dla;          hfb(nd1,nd2)=-hla
      // no need back       End Do
      // no need back       hfb(n1,n1)=hfb(n1,n1)-al; hfb(nd1,nd1)=hfb(nd1,nd1)+al
      // no need back    End Do
      // no need back    Call sdiag(ndx2,nhfb,hfb,evvk,hfb,zhfb,+1)
      // no need back End If
      // wait!------------------------------------------------------------------
      // wait! Blocking
      // wait!------------------------------------------------------------------
      // wait! external blocking
      // waitIf(iiter.Eq.1.And.inner(it).Eq.0) Then
      // wait   If(iparenti(it).Ne.0.And.keyblo(it).Eq.0) Then
      // wait      ! eventually charging
      // wait      !   keyblo(it)=1
      // wait      !   bloblo(keyblo(it),it)=ib
      // wait      !   blo123(keyblo(it),it)=requested level (k0)
      // wait      Call requested_blocked_level(ib,it)
      // wait      If(ierror_flag.Ne.0) Return
      // wait  ibiblo=bloblo(keyblo(it),it)
      // wait   Endif
      // waitEndif
      // wait! general blocking
      // waitk0=0
      // waitIf(ibiblo.Eq.ib) Then
      // wait   If(iiter.Eq.1.And.inner(it).Eq.0) Then
      // wait      ! blocked level as in the even-even nucleus
      // wait      k0=blo123(keyblo(it),it); ndk=k0+nd
      // wait      Do n2=1,nd
      // wait         nd2=n2+nd
      // wait         hfb1(n2,it)=hfb(n2,ndk)    !U
      // wait         hfb1(nd2,it)=hfb(nd2,ndk)  !V
      // wait      Enddo
      // wait      ! number of states in the block to be tested
      // wait      blocross(it)=Min(blomax(it)+10,nd)
      // wait   Endif
      // wait   ! overlap between new and old blocked levels
      // wait   s3=zero
      // wait   Do n1=1,blocross(it)
      // wait      ndk=n1+nd; s1=zero
      // wait      Do n2=1,nd
      // wait         nd2=n2+nd
      // wait         s1=s1+Abs(hfb1(nd2,it)*hfb(nd2,ndk)) !VV
      // wait         s1=s1+Abs(hfb1(n2,it)*hfb(n2,ndk))   !UU
      // wait      Enddo
      // wait      If(s1.Gt.s3) Then
      // wait         s3=s1; k0=n1
      // wait      Endif
      // wait   Enddo
      // wait   blo123d(it)=k0
      // wait   If(.Not.norm_to_improve) Then
      // wait      ! find maximal HO component
      // wait      ndk=k0+nd
      // wait      s1=zero
      // wait      Do n1=1,nd
      // wait         nd1=n1+nd
      // wait         hfb1(n1,it)=hfb(n1,ndk); hfb1(nd1,it)=hfb(nd1,ndk)
      // wait         s2=Max(s1,Abs(hfb(n1,ndk)),Abs(hfb(nd1,ndk)))
      // wait         If(s2.Gt.s1) Then
      // wait            s1=s2; i=n1+i0  ! labels in k[k1,k2] numbering
      // wait         End If
      // wait      End Do
      // wait      ! print blocked state
      // wait      Do iw=lout,lfile
      // wait         Write(iw,'(4x,a,2(a,i3),2x,3(a,1x,f12.8,1x),(i3,a,i3,1x),a)')  &
      // wait              protn(it),' Blocking: block=',ib,  &
      // wait              ' state=',k0,  &
      // wait              ' Eqp=',evvk(k0+nd),  &
      // wait              ' Dqpe=',evvk(k0+nd)-eqpmin(it),  &
      // wait              ' Ovlp=',s3  &
      // wait              , keyblo(it),'/',blomax(it)  &
      // wait              , tb(i)
      // wait      Enddo
      // wait      ! ieresbl=6, 'BLKN','BLKZ'
      // wait      ereslbl(it)=tb(i)
      // wait      If(it.Eq.1) Then
      // wait         ! 'BlEqpN','BlDEqpN','BlOvrN'
      // wait         eresbl(1)=evvk(k0+nd); eresbl(2)=evvk(k0+nd)-eqpmin(it); eresbl(3)=s1
      // wait      Else
      // wait         ! 'BlEqpZ','BlDEqpZ','BlOvrZ'
      // wait         eresbl(4)=evvk(k0+nd); eresbl(5)=evvk(k0+nd)-eqpmin(it); eresbl(6)=s1
      // wait      Endif
      // wait   End If
      // waitEnd If
      //!------------------------------------------------------------------
      //! Run over all qp states k in the block
      //!------------------------------------------------------------------
      // std::cout << "kaib: " << std::endl;
      kaib = kl;
      for (k = 1; k <= nd; k++)
      {
        ndk = k + nd;
        //! referent spectra
        pn = zero;
        for (i = 1; i <= nd; i++)
        {
          hla = pow(hfb[-1 + i + nd][-1 + ndk], 2);
          pn = pn + hla;
          // no need back If(i.Eq.1) Then
          // no need back    !write(*,*) hfb(i+nd,ndk),hfb(i,ndk),evvk(nd+k)         ! Vak, Uak, Sign( V*U)  Ek
          // no need back    !write(*,*) hfb(i,nd-k+1),hfb(i+nd,nd-k+1),evvk(nd-k+1) ! Vak, Uak, Sign(-V*U) -Ek
          // no need back Endif
        }
        // blocking ? ! Blocking
        // blocking ? If(k.Eq.k0) Then
        // blocking ?    n1=k0+nd
        // blocking ?    Do i=1,nd
        // blocking ?       hla=hfb(i+nd,n1)**2; dla=hfb(i,n1)**2; pn=pn-half*(hla-dla)
        // blocking ?    Enddo
        // blocking ? Endif
        eqpe = evvk[-1 + nd + k];
        ela = eqpe * (one - two * pn);
        enb = ela + al;
        ekb = sqrt(abs(eqpe * eqpe - pow(ela, 2)));
        //!------------------------------------------------------------------
        //! cut-off condition: energy pwi + Fermi cut-off function
        //!------------------------------------------------------------------
        lpr_pwi = enb <= pwi && abs(one / (one + exp(100.0 * (enb - pwi)))) > 0.000001;
        // #ifndef hide_qrpa
        lpr_pwi = enb <= pwi; //! cristina sharp cut off for qrpa
        // #endif
        if (use_TMR_pairing == 0)
          lpr_pwi = true; //! no pairing window with TMR pairing
        //!------------------------------------------------------------------
        //! Remember the whole qp solution
        //!------------------------------------------------------------------
        if (!norm_to_improve)
        {
          i_eqp = i_eqp + 1;
          (*EqpPo)[-1 + i_eqp] = evvk[-1 + nd + k]; //! Eqp_k
          if (lpr_pwi)
            (*KqpPo)[-1 + kl + 1] = i_eqp; //! below pwi otherwise zero
          if (lpr_pwi)
            (*KpwiPo)[-1 + kl + 1] = i_uv; //! below pwi otherwise zero
          for (n2 = 1; n2 <= nd; n2++)
          {
            nd2 = n2 + nd;
            i_uv = i_uv + 1;
            (*UqpPo)[-1 + i_uv] = hfb[-1 + n2][-1 + ndk];  //! U_ak
            (*VqpPo)[-1 + i_uv] = hfb[-1 + nd2][-1 + ndk]; //! V_ak
          }
        }
        //!------------------------------------------------------------------
        //! Pairing window
        //!------------------------------------------------------------------
        if (lpr_pwi)
        {
          kl = kl + 1; //! number of active states
          // if(k0.Eq.k) blok1k2d(it)=kl                  //!blocking: dynamic #: k[k1,k2] numbering
          if ((eqpe <= emin) && (pn > 0.0001))
          { //! to avoid unocc at magic numbers
            emin = eqpe;
            alnorm = pn; //! min qpe and its occupation
          }
          erhfb[-1 + kl] = enb;
          drhfb[-1 + kl] = ekb;
          uk[-1 + kl][-1 + it] = pn; //! ref.s.p. energies, deltas, occupancies
          std::cout << "pn: " << pn << std::endl;
          sumnz[-1 + it] = sumnz[-1 + it] + two * pn; //! internal normalization
        }
      }
      if (norm_to_improve)
        continue;
      //!------------------------------------------------------------------
      //!  Density matrices
      //!------------------------------------------------------------------
      kdib = kl - kaib;
      ka[-1 + ib][-1 + it] = kaib;
      kd[-1 + ib][-1 + it] = kdib;
      k1 = kaib + 1;
      k2 = kaib + kdib;
      eqpe = 0.0;
      for (n2 = 1; n2 <= nd; n2++)
      {
        for (n1 = n2; n1 <= nd; n1)
        {
          s1 = zero;
          s2 = zero;
          if (k1 <= k2)
          {
            for (k = k1; k <= k2; k++)
            {
              nd1 = (*KpwiPo)[-1 + k] + n1;
              nd2 = (*KpwiPo)[-1 + k] + n2;
              s1 = s1 + (*VqpPo)[-1 + nd1] * (*VqpPo)[-1 + nd2];
              s2 = s2 + (*UqpPo)[-1 + nd1] * (*VqpPo)[-1 + nd2] + (*VqpPo)[-1 + nd1] * (*UqpPo)[-1 + nd2];
            }
            //! if(it.eq.1) write(200,*) n1,n2,s1-s2
            s1 = two * s1;
            s2 = half * s2; //! two:due to m-projection, half:due to symmetrization
                            //! blocking
                            // If(ibiblo.Eq.ib) Then
            //    i=blok1k2d(it); id1=KpwiPo(i)+n1; id2=KpwiPo(i)+n2
            //    s1=s1-VqpPo(id1)*VqpPo(id2)+UqpPo(id1)*UqpPo(id2)
            //    s2=s2-half*(UqpPo(id1)*VqpPo(id2)+VqpPo(id1)*UqpPo(id2))
            // Endif
          }
          n12 = n1 + (n2 - 1) * nd;
          n21 = n2 + (n1 - 1) * nd;
          rk[-1 + n12][-1 + m] = s1;
          rk[-1 + n21][-1 + m] = s1; //!  V V'
          ak[-1 + n12][-1 + m] = -s2;
          ak[-1 + n21][-1 + m] = -s2; //!- U V', ak=half*(pairing density)
          hfbcan[-1 + n1][-1 + n2] = s1;
          hfb[-1 + n1][-1 + n2] = s1;
        } //! n1
      }   //! n2
      //!------------------------------------------------------------------
      //! Canonical basis
      //!------------------------------------------------------------------
      // If(k1.Le.k2) Then
      //    Call Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
      //    If(ierror_flag.Ne.0) Return
      // Endif
      // lcanon(ib,it)=lc
    }
    //! ib
    // If(kl.Eq.0) Then
    //    ierror_flag=ierror_flag+1
    //    ierror_info(ierror_flag)=' STOP: kl=zero, no states below pwi!!!'
    //    Return
    // Endif
    // If(iparenti(it).Ne.0.And.ibiblo.Eq.0) Then
    //    ierror_flag=ierror_flag+1
    //    ierror_info(ierror_flag)='STOP: No blocking candidate found!!!'
    //    Return
    // Endif
    eqpmin[-1 + it] = emin;
    klmax[-1 + it] = kl;
    sumnz[-1 + it] = sumnz[-1 + it] - tz[-1 + it];
    //!------------------------------------------------------------------
    //! Lambda search
    //!------------------------------------------------------------------
    Alambda(al, it, kl);
    // If(ierror_flag.Ne.0) Return
    if (keyblo[-1 + it] == 0)
      ala[-1 + it] = al;
    else
      ala[-1 + it] = ala[-1 + it] + 0.50 * (al - ala[-1 + it]);
    //
    //! NB! 'alast' instead of 'al' at small pairing
    alast[-1 + it] = al;
    if (abs(ept[-1 + it]) < 0.0001)
    {
      ntz = tz[-1 + it] + 0.1;
      ntz = ntz / 2;
      for (k = 1; k <= kl; k++)
        drhfb[-1 + k] = erhfb[-1 + k];
      //
      ord(kl, drhfb);
      alast[-1 + it] = drhfb[-1 + ntz]; //! last bound s.p. energy
    }
    //!------------------------------------------------------------------
    //! THO asymptotic decay
    //!------------------------------------------------------------------
    //! density asymptotic decay \rho(r)->Exp(-ass(it)*r)
    //! ass(it)=2*Sqrt((E_min-\lambda)/((A-1)/A)*hbar**2/(2*m)))
    al2 = zero;
    if (kindhfb < 0)
      al2 = al + two * ala2[-1 + it] * (one - two * alnorm); //! al=al+two*ala2(it)
    //
    al2 = (emin - al2) / hb0;
    //! wrong asymptotic
    iasswrong[-1 + it] = 0;
    if (al2 <= zero)
      iasswrong[-1 + it] = 1;
    ass[-1 + it] = two * sqrt(abs(al2));
    //!
  } //! While(norm_to_improve)
}

void HFBTHO_solver::Alambda(double al, int it, int kl)
{
  //!---------------------------------------------------------------------
  //! Adjusting Fermi energy
  //!---------------------------------------------------------------------
  // Use HFBTHO
  // Implicit None
  double fm7 = 0.0000001, fm10 = 0.0000000001;
  double vh, xinf, xsup, esup, ez, dez, dvh, y, a, b, einf, absez, sn;
  int i, k, icze, lit, ntz;
  //
  // std::cout << "Alambda execute: " << std::endl;
  //!-------------------------------------------------
  //! Chemical potential without pairing
  //!-------------------------------------------------
  if (CpV0[it - 1] == zero)
  {
    ntz = tz[-1 + it] + 0.1;
    ntz = ntz / 2;
    for (k = 1; k <= kl; k++)
      drhfb[-1 + k] = erhfb[-1 + k];
    ord(kl, drhfb);
    if (ntz < kl)
      al = half * (drhfb[-1 + ntz] + drhfb[-1 + ntz + 1]);
    else
      al = drhfb[-1 + ntz] + 0.001;
    //
    return;
  }
  //!-------------------------------------------------
  //! Chemical potential with pairing
  //!-------------------------------------------------
  xinf = -100.0;
  xsup = 80.0;
  esup = one;
  icze = 0;
  for (lit = 1; lit <= 500; lit++)
  {
    sn = zero;
    dez = zero;
    for (i = 1; i <= kl; i++)
    {
      vh = zero;
      dvh = zero;
      y = erhfb[-1 + i] - al;
      a = y * y + pow(drhfb[-1 + i], 2);
      b = sqrt(a);
      if (b > zero)
        vh = half * (one - y / b);
      if (b < fm7 && icze == 1)
        vh = -einf / (esup - einf); //! no pairing
      if (vh < zero)
        vh = zero;
      if (vh > one)
        vh = one;
      if (b > zero)
        dvh = pow(drhfb[-1 + i], 2) / (a * b); //! D[ez,al](i)
      ////! blocking
      // if (i.Eq.blok1k2d(it))
      //{
      //   vh = half;
      //   dvh = zero;
      // }
      ////
      sn = sn + two * vh;
      dez = dez + dvh; //! D[ez,al]
    }
    ez = sn - tz[-1 + it];
    absez = abs(ez);
    //!-------------------------------------------------
    //! Correcting bounds
    //!-------------------------------------------------
    if (ez < zero)
    {
      xinf = std::max(xinf, al);
      einf = ez;
    }
    else
    {
      xsup = std::min(xsup, al);
      esup = ez;
    }
    //
    if (lit == 1)
    {
      if (absez <= 0.10)
        al = al - ez;
      else
        al = al - 0.10 * copysign(one, ez);
      //
    }
    else
      al = al - ez / (dez + pow(10, -20)); //! newton method
    //
    if (xsup - xinf < fm7)
      icze = 1; //! low/upp close
    if (al < xinf || al > xsup)
      al = half * (xinf + xsup); //! mean upp/low
    if (absez <= fm10)
      return;
  }
  //!-------------------------------------------------
  //! Low accuracy warning
  //!-------------------------------------------------
  // Write(lout,'(a,2(e12.5,2x),a,2(2x,f8.4),a,i2)') ' Low accuracy=',sn,ez,' for N,Z=',tz,' it=',it
}

int HFBTHO_solver::inout(int is)
{

  //!---------------------------------------------------------------------
  //! is=1: reads matrix elements from tape and exit
  //! is=2: writes matrix elements to tape and exit
  //! NB! if the welfile is missing or corrupt call start
  //!     to restart calculations from scratch
  //!---------------------------------------------------------------------
  // Use HFBTHO
  // Implicit None
  int iw, n1, n2, nd, ib, bloall1;
  std::string nucname1;
  std::string filelabel;
  std::string welfile;
  double tz1[2], b01, beta1, v0r[2], v1r[2], pwir;
  int npr1, npr11, ngh1, ngl1, n001, nt1;
  int ntx1, nb1, nhhdim1, NLANSA0, NLANSA1, NZA2NRA, NZA1, NLA1;
  int ID1[nbx];
  //!==== HFBODD interface
  int i, nza, nra, nla, nsa, ibasis;
  int ibro;
  //!====
  //! label organization
  // not used yet FileLabels(NPR, ININ, FILELABEL);
  // If(ierror_flag.Ne.0) Return
  // If(iLST1.Le.0) Write(welfile,'(a8,a4)')  FILELABEL,'.hel'
  // If(iLST1.Gt.0) Write(welfile,'(a8,a4)')  FILELABEL,'.tel'
  //!
  if (is == 1)
  {
    //!---------------------------------------------------------------------
    //! read matrix elements from 'welfile' file or start from scratch
    //!---------------------------------------------------------------------
    if (inin > 0)
    {
      start();
      return 0;
    }
    // not consider now     !read matrix elements
    // not consider now     Open(lwin,file=welfile,status='unknown',form='unformatted', ERR=100)
    // not consider now     !hel
    // not consider now     Read(lwin,Err=100,End=100) nucname1,npr11,npr1,ngh1,ngl1,n001,nb1,nt1
    // not consider now     If(Abs(n001).Ne.Abs(n00).And.nb1.Ne.nb) go to 100
    // not consider now     Read(lwin,Err=100,End=100) b01,beta1,si,etot,rms,bet,xmix,v0r,v1r,pwir  &
    // not consider now          ,del,ept,ala,ala2,alast,tz1,varmas,varmasNZ,pjmassNZ,ass,skass
    // not consider now     brin=zero; si=one; bbroyden='L'
    // not consider now     Read(lwin,Err=100,End=100) ntx1,nb1,nhhdim1
    // not consider now     Read(lwin,Err=100,End=100) id1
    // not consider now     Read(lwin,Err=100,End=100) brin
    // not consider now     !
    // not consider now     ! Add small pairing de=de+0.1 in the no-LN
    // not consider now     ! case to prevent pairing collaps
    // not consider now     If(kindhfb.Eq.1.And.Add_Pairing) Then
    // not consider now        ibro=0
    // not consider now        Do ib=1,NB
    // not consider now           ND=ID1(ib)
    // not consider now           I=ibro
    // not consider now           Do N1=1,ND
    // not consider now              Do N2=1,N1
    // not consider now                 I=I+1
    // not consider now                 brin(i+nhhdim2)=brin(i+nhhdim2)+0.10_pr
    // not consider now                 brin(i+nhhdim3)=brin(i+nhhdim3)+0.10_pr
    // not consider now              End Do !N2
    // not consider now           End Do !N1
    // not consider now           ibro=i
    // not consider now        End Do !IB
    // not consider now     Endif
    // not consider now     Do ib=1,NB
    // not consider now        ND=ID1(ib)
    // not consider now        Do N1=1,ND
    // not consider now           Read(lwin,Err=100,End=100) NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1
    // not consider now        End Do
    // not consider now     End Do
    // not consider now     ! blocking
    // not consider now     Read(lwin,Err=100,End=100)  bloall1
    // not consider now     Read(lwin,Err=100,End=100)  bloblo,blo123,blok1k2,blomax,bloqpdif
    // not consider now     If(bloall1.Ne.bloall) go to 100
    // not consider now     !tel
    // not consider now     If(iLST.Gt.0) Then
    // not consider now        Read(lwin,Err=100,End=100) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
    // not consider now        If(Allocated(fdsx)) Deallocate(fdsx,fdsy,fdsy1,fdsy2,fdsy3,  &
    // not consider now             fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,  &
    // not consider now             fspb3,fspc3,fspd3)
    // not consider now        Allocate(fdsx(iqqmax),fdsy(iqqmax),fdsy1(iqqmax),  &
    // not consider now             fdsy2(iqqmax),fdsy3(iqqmax),fspb0(iqqmax),fspc0(iqqmax),  &
    // not consider now             fspd0(iqqmax),fspb1(iqqmax),fspc1(iqqmax),fspd1(iqqmax),  &
    // not consider now             fspb2(iqqmax),fspc2(iqqmax),fspd2(iqqmax),fspb3(iqqmax),  &
    // not consider now             fspc3(iqqmax),fspd3(iqqmax))
    // not consider now        Read(lwin,Err=100,End=100) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
    // not consider now             ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
    // not consider now     End If
    // not consider now     Do iw=lout,lfile
    // not consider now        Write(iw,*)
    // not consider now        Write(iw,*) ' Reading from wel_file: ',welfile
    // not consider now        Write(iw,*)
    // not consider now     Enddo
    // not consider now     Close(lwin)
    // not consider now     !---------------------------------------------------------------------
    // not consider now     !for HFODD interface
    // not consider now     !ib=abs(inin); npr_temp(1)=npr11; npr_temp(2)=npr1;
    // not consider now     !Call HFBTHO_HFODD(1616,npr_temp,ib)
    // not consider now     !---------------------------------------------------------------------
    // not consider now     Return
    // not consider now     !
    // not consider now100  Continue
    // not consider now     !---------------------------------------------------------------------
    // not consider now     ! missing or corrupt 'welfile' file
    // not consider now     !---------------------------------------------------------------------
    // not consider now     Close(lwin)
    // not consider now     Do iw=lout,lfile
    // not consider now        Write(iw,'(1x,a,a,a)')
    // not consider now        Write(iw,'(1x,a,a,a)')   ' The file ',welfile,' is corrupted, missing, or wrong!'
    // not consider now        Write(iw,'(1x,a,a,a)')   ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!'
    // not consider now        Write(iw,'(1x,a,a,a)')
    // not consider now     Enddo
    // not consider now     Call start
    // not consider now     Return
  }
  // not consider now  !---------------------------------------------------------------------
  // not consider now  ! write matrix elements to 'welfile' file
  // not consider now  !---------------------------------------------------------------------
  // not consider now  If (is.Eq.2.And.iasswrong(3).Eq.0) Then
  // not consider now  ! don't write with mpi qrpa
  // not consider now#ifdef hide_mpi_qrpa
  // not consider now     Open(lwou,file=welfile,status='unknown',form='unformatted')
  // not consider now     !hel
  // not consider now     npr11=npr(1); npr1=npr(2)
  // not consider now     Write(lwou) nucname,npr11,npr1,ngh,ngl,n00,nb,nt
  // not consider now     Write(lwou) b0,beta0,si,etot,rms,bet,xmix,CpV0,CpV1,pwi,del,ept  &
  // not consider now          ,ala,ala2,alast,tz,varmas,varmasNZ,pjmassNZ,ass,skass
  // not consider now     Write(lwou) ntx,nb,nhhdim
  // not consider now     Write(lwou) id
  // not consider now     Write(lwou) brin
  // not consider now     ibasis=0
  // not consider now     Do ib=1,NB
  // not consider now        ND=ID(ib)
  // not consider now        Do N1=1,ND
  // not consider now           ibasis=ibasis+1
  // not consider now           NLA=NL(ibasis); NRA=NR(ibasis); NZA=NZ(ibasis); NSA=NS(ibasis); NLANSA1=(-1)**(NZA+NLA)
  // not consider now           Write(lwou) 2*NLA+NSA,NLANSA1,NZA+2*NRA+NLA,NZA,NLA
  // not consider now        End Do
  // not consider now     End Do
  // not consider now     !---------------------------------------------------------------------
  // not consider now     ! blocking: sort blocking candidates first
  // not consider now     !---------------------------------------------------------------------
  // not consider now     Do ib=1,2
  // not consider now        Call blosort(ib,blomax(ib))
  // not consider now     Enddo
  // not consider now     Write(lwou) bloall
  // not consider now     Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
  // not consider now     !tel
  // not consider now     If(Allocated(fdsx)) Then
  // not consider now        Write(lwou) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
  // not consider now        Write(lwou) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
  // not consider now             ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
  // not consider now     End If
  // not consider now     Close(lwou)
  // not consider now     Do iw=lout,lfile
  // not consider now        Write(iw,'(a,a,a)')
  // not consider now        Write(iw,'(a,a,a)') '  Writing to wel_file: ',welfile
  // not consider now        Write(iw,'(a,a,a)') ' __________________________________  '
  // not consider now        Write(iw,'(a,a,a)') '  The tape ',welfile,' recorded:     '
  // not consider now        Write(iw,'(a,a,a)') '  nucname,npr,ngh,ngl,n00,nb,nt      '
  // not consider now        Write(iw,'(a,a,a)') '  b0,beta0,si,etot,rms,bet,xmix      '
  // not consider now        Write(iw,'(a,a,a)') '  pairing:     CpV0,CpV1,pwi         '
  // not consider now        Write(iw,'(a,a,a)') '  delta:       del,ept               '
  // not consider now        Write(iw,'(a,a,a)') '  lambda:      ala,ala2,alast,tz     '
  // not consider now        Write(iw,'(a,a,a)') '  asymptotic:  varmas,ass,skass      '
  // not consider now        Write(iw,'(a,a,a)') '  ntx,nb,nhhdim,id,N_rz,n_r,n_z      '
  // not consider now        Write(iw,'(a,a,a)') '  Omega2,Sigma2,Parirty,Lambda       '
  // not consider now        Write(iw,'(a,a,a)') '  matrices(inbro):    hh,de          '
  // not consider now        Write(iw,'(a,a,a)') '  *all blocking candidates           '
  // not consider now        If(Allocated(fdsx)) Write(iw,'(a,a,a)') '  *all THO arrays                    '
  // not consider now        Write(iw,'(a,a,a)') ' __________________________________  '
  // not consider now        Write(iw,'(a,a,a)')
  // not consider now     Enddo
  // not consider now#else
  // not consider now   ! with qrpa no hel file: do nothing here
  // not consider now#endif
  // not consider now  Endif
  // not consider now  !
  // not consider now  !ib=abs(inin); npr_temp(1)=npr11; npr_temp(2)=npr1;
  // not consider now  !Call HFBTHO_HFODD(1717,npr_temp,ib)
  // not consider now  !
  // not consider now
  return 0;
}

void HFBTHO_solver::start()
{
  std::cout << "start() execute" << std::endl;
  // Subroutine start
  //   !---------------------------------------------------------------------
  //   ! initializes scratch Saxon-Woods potentials
  //   !---------------------------------------------------------------------
  //   Use HFBTHO
  //   Implicit None
  int iw, i, ih, il, ihl, it, ita;
  double zb[ngh], rrb[ngl], rb[ngl], rav, rao, vpws, vls, betas, gamma;
  double fac, facb, zz, rr, r, ctet, cphi, p2, p20, p22, s, u, w, f, rc, c, beta00;
  //!--------------------------------------------------------------------
  //! initializing
  //!--------------------------------------------------------------------
  //!----------------------------------------------------------------------------
  //! reinatializing all again since scratch calculation
  //!----------------------------------------------------------------------------
  iniialize_HFBTHO_SOLVER();
  // If(ierror_flag.Ne.0) Return
  //! mario
  for (it = itmin; it < itmax; it++)
  {
    if (npr[-1 + it] != 2 * (npr[-1 + it] / 2))
    {
      irestart = irestart + 1;
      npr[-1 + it] = npr_INI[-1 + it];
    }
    //
  }
  // If(irestart.Ne.0) Then
  //    ! odd nucleus requested but no even-even solution
  //    ! recalculate the even-even nucleus from scratch
  //    Do iw=lout,lfile
  //       Write(iw,'(1x,a,2i4)')
  //       Write(iw,'(1x,a,2i4)') ' Initialization for the even-even core (N,Z)=: ',npr(1:2)
  //    Enddo
  // Else
  //    ! scratch for the even-even nucleus requested
  //    Do iw=lout,lfile
  //       Write(iw,'(1x,a,2i4)')
  //       Write(iw,'(a,a,3i4)')    '  Scratch initialization for the nucleus: ',nucname,npr(1:2)
  //       Write(iw,'(1x,a,2i4)')
  //    Enddo
  // Endif
  n00 = abs(n00_INI);
  b0 = b0_INI;
  q = q_INI;
  iLST = iLST_INI;
  maxi = MAX_ITER_INI;
  inin = inin_INI;
  skyrme = skyrme_INI;
  kindhfb = kindhfb_INI;
  cdef = cdef_INI;
  cqad = cqad_INI;
  iproj = iproj_INI;
  npr1pj = npr1pj_INI;
  npr2pj = npr2pj_INI;
  Constraint_or_not(inin_INI, inin, icstr);
  // If(ierror_flag.Ne.0) Return
  preparer(false);
  // If(ierror_flag.Ne.0) Return
  inin = abs(inin); //! positive even if inin_INI is not
  bet = cdef;
  //!-----------------------------------
  //! Saxon-Woods potentials
  //!-----------------------------------
  // Do iw=lout,lfile
  //    Write(iw,'(/,a)') '  Initial potentials of Saxon-Woods shape '
  // End Do
  beta00 = bet; //! wf to requested deformation
  // Do iw=lout,lfile
  //    Write(iw,'(a,f10.4)')  '  v0ws   =',v0ws
  //    Write(iw,'(a,f10.4)')  '  kappa  =',akv
  //    Write(iw,'(a,2f10.4)') '  vs0    =',vso
  //    Write(iw,'(a,2f10.4)') '  r0     =',r0v
  //    Write(iw,'(a,2f10.4)') '  a      =',av
  //    Write(iw,'(a,2f10.4)') '  r0-so  =',rso
  //    Write(iw,'(a,2f10.4)') '  a-so   =',aso
  //    Write(iw,'(a,f10.4)')  '  beta00 =',beta00
  // Enddo
  //!-----------------------------------
  //! Densities
  //!-----------------------------------
  for (it = itmin; it < itmax; it++)
  {
    ita = 3 - it;
    rav = r0v[-1 + it] * pow(amas, p13);
    rao = rso[-1 + it] * pow(amas, p13);
    vpws = v0ws * (one - akv * (npr[-1 + it] - npr[-1 + ita]) / amas);
    vls = half * pow((hqc / amu), 2) * vpws * vso[-1 + it];
    betas = beta00 * sqrt(5.0 / (16.0 * PI));
    gamma = zero;
    fac = one + betas * cos(gamma * PI / 180.0);
    fac = (one + betas * cos((gamma + 120.0) * PI / 180.0)) * fac;
    fac = (one + betas * cos((gamma - 120.0) * PI / 180.0)) * fac;
    fac = pow(fac, (-p13));
    //! z,r-coordinates in fm
    // zb = xh * bz;
    for (int i = 0; i < xh.size(); i++)
    {
      zb[i] = xh[i] * bz;
    }
    //
    // rrb = xl * bp * bp;
    for (int i = 0; i < xl.size(); i++)
    {
      rrb[i] = xl[i] * bp * bp;
    }
    //
    // rb = Sqrt(rrb);
    //
    for (int i = 0; i < sizeof(rrb) / sizeof(rrb[0]); i++)
    {
      rb[i] = sqrt(rrb[i]);
    }

    for (ih = 1; ih <= ngh; ih++)
    {
      zz = pow(zb[-1 + ih], 2);
      for (il = 1; il <= ngl; il++)
      {
        rr = rrb[-1 + il] + zz;
        r = sqrt(rr);
        //! woods saxon
        ctet = zz / rr;
        cphi = zero;
        p20 = 3.0 * ctet - one;
        p22 = sqrt(3.0) * cphi;
        p2 = p20 * cos(gamma * PI / 180) + p22 * sin(gamma * PI / 180.0);
        facb = fac * (one + betas * p2);
        u = vpws / (one + exp((r - rav * facb) / av[-1 + it]));
        w = -vls / (one + exp((r - rao * facb) / aso[-1 + it]));
        ihl = ih + (il - 1) * ngh;
        if (it == 1)
        {
          vhbn[-1 + ihl] = hb0;
          vn[-1 + ihl] = u;
          vsn[-1 + ihl] = w;
          vrn[-1 + ihl] = zero;
          vzn[-1 + ihl] = zero;
          vdn[-1 + ihl] = zero;
          vSFIZn[-1 + ihl] = zero;
          vSZFIn[-1 + ihl] = zero;
          vSFIRn[-1 + ihl] = zero;
          vSRFIn[-1 + ihl] = zero;
        }
        else
        {
          vhbp[-1 + ihl] = hb0;
          vp[-1 + ihl] = u;
          vsp[-1 + ihl] = w;
          vrp[-1 + ihl] = zero;
          vzp[-1 + ihl] = zero;
          vdp[-1 + ihl] = zero;
          vSFIZp[-1 + ihl] = zero;
          vSZFIp[-1 + ihl] = zero;
          vSFIRp[-1 + ihl] = zero;
          vSRFIp[-1 + ihl] = zero;
        }
        ro[-1 + ihl][-1 + it] = u;
        aka[-1 + ihl][-1 + it] = pow(10, -3) * exp((r - rav * facb) / 2.0);
      }
    }
    // s = npr[-1 + it] / Sum(ro( :, it));
    double sum_ro_it = 0.0;
    for (int i = 0; i < sizeof(ro) / sizeof(ro[0]); i++)
      sum_ro_it = sum_ro_it + ro[i][-1 + it];
    //
    for (il = 1; il <= ngl; il++)
    {
      for (ih = 1; ih <= ngh; ih++)
      {
        ihl = ih + (il - 1) * ngh;
        f = s / (PI * wh[-1 + ih] * wl[-1 + il] * bz * bp * bp);
        ro[-1 + ihl][-1 + it] = f * ro[-1 + ihl][-1 + it];
      }
    }
    //!-----------------------------------
    //! pairing
    //!-----------------------------------
    for (il = 1; il <= nghl; il++)
    {
      f = (ro[-1 + il][-1 + 1] + ro[-1 + il][-1 + 2]) / rho_c;
      if (it == 1)
        dvn[-1 + il] = CpV0[it - 1] * (one - CpV1[it - 1] * f) * aka[-1 + il][-1 + it];
      else
        dvp[-1 + il] = CpV0[it - 1] * (one - CpV1[it - 1] * f) * aka[-1 + il][-1 + it];
      //
    }
  }
  //!-----------------------------------
  //! coulomb
  //!-----------------------------------
  if (icou == 0)
  // cou = zero;
  {
    for (int i = 0; i < sizeof(cou) / sizeof(cou[0]); i++)
      cou[i] = zero;
  }
  else
  {
    rc = r0v[-1 + 2] * pow(amas, p13);
    for (il = 1; il <= ngl; il++)
    {
      for (ih = 1; ih <= ngh; ih++)
      {
        r = sqrt(pow(zb[-1 + ih], 2) + rrb[-1 + il]);
        if (r < rc)
          c = half * (3 / rc - r * r / pow(rc, 3));
        else
          c = one / r;
        //
        cou[-1 + ih + (il - 1) * ngh] = c * npr[-1 + 2] / alphi;
      }
    }
  }
  //!-----------------------------------
  //! initial ph+pp matrix elements
  //!-----------------------------------
  // ak = 0.1;
  for (auto &m : ak)
    for (auto &n : m)
      n = 0.1;
  //
  // rk = 0.1; //! initial density matrix elements (improve later)
  for (auto &m : rk)
    for (auto &n : m)
      n = 0.1;
  //
  // brin = zero; //! initial matrix elements to zero
  for (auto &m : brin)
    m = zero;
  iiter = 0; //! iteration number iiter to zero
  std::cout << "before execute gamdel(): " << std::endl;
  gamdel();
  //!
  // End Subroutine start
}

void HFBTHO_solver::gamdel()
{
  std::cout << "Inside gamdel: " << std::endl;
  //!---------------------------------------------------------------------
  //! ph- and pp- matrices in configurational space
  //!---------------------------------------------------------------------
  // Use HFBTHO
  // Use pairing_HFBTHO !mario
  // Implicit None
  int i, ih, il, ib, ibx, nd, nd2, nza, nra, nla, nsa, nsb, nsab;
  // int ihil, laplus, im, JA, N1, N2, ndnd, n12, n21;
  int ihil, laplus, IM, JA, N1, N2, ndnd, n12, n21;
  int i1, i2, i3;
  double qla, yi, y, y2, qha, qhla, xmi, u2, un, up, xxx;
  double sml2, cnzaa, cnraa, SSU, SSD;
  double FITW1, FITW2, FITW3, FITW4;
  double fi1r, fi1z, fi2d, QHL1A, QH1LA;
  double vh, vdh, vsh, hbh, vsum;
  double SRFIh, SFIRh, SFIZh, SZFIh, SNABLARh, SNABLAZh;
  // double xlam, xlam2, xlamy, xlamy2, xlap, xlap2, xlapy, xlapy2, XLAMPY;
  double xlam, xlam2, xlamy, xlamy2, xlap, xlap2, xlapy, xlapy2, xlampy;
  double FIUN1, FIDN1, FIURN1, FIDRN1, FIUZN1, FIDZN1, FIUD2N1, FIDD2N1;
  double FIUN2, FIDN2, FIURN2, FIDRN2, FIUZN2, FIDZN2, FIUD2N2, FIDD2N2;
  double FIUN12, FIDN12, FIURN12, FIDRN12, FIUZN12, FIDZN12;
  double vnhl, vrnhl, vznhl, vdnhl, vsnhl, vhbnhl, vSRFInhl, vSFIRnhl;
  double vSFIZnhl, vSZFInhl, vphl, vrphl, vzphl, vdphl, vsphl, vhbphl;
  double vSRFIphl, vSFIRphl, vSFIZphl, vSZFIphl, dvnhl, dvphl;
  int ibro;
  //!
  // Call get_CPU_time('gamdel',0)
  //!
  //!----------------------------------------------
  //! START BLOCKS
  //!----------------------------------------------
  for (int qqi = 0; qqi < sizeof(brout) / sizeof(brout[0]); qqi++)
    brout[qqi] = 0;
  // brout=zero;
  ibro = 0;
  for (ib = 1; ib <= nb; ib++)
  {
    nd = id[-1 + ib];
    IM = ia[-1 + ib];
    ibx = ib + nbx;
    if (Parity)
      laplus = (ib + 1) / 2; //! Yesp
    else
      laplus = ib; //! Nop
    // Endif
    // XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
    xlap = laplus;
    xlam = xlap - one;
    xlap2 = xlap * xlap;
    xlam2 = xlam * xlam;
    //!----------------------------------------------
    //! SUM OVER GAUSS INTEGRATION POINTS
    //!----------------------------------------------
    for (ihil = 1; ihil <= nghl; ihil++)
    {
      y = y_opt[-1 + ihil];
      xlamy = xlam * y;
      xlapy = xlap * y;
      // XLAMPY=XLAMY+XLAPY;
      xlampy = xlamy + xlapy;
      y2 = y * y;
      xlamy2 = xlam2 * y2;
      xlapy2 = xlap2 * y2;
      //!
      vnhl = vn[-1 + ihil];
      vrnhl = vrn[-1 + ihil];
      vznhl = vzn[-1 + ihil];
      vdnhl = vdn[-1 + ihil];
      vsnhl = vsn[-1 + ihil];
      vhbnhl = vhbn[-1 + ihil];
      vSRFInhl = vSRFIn[-1 + ihil];
      vSFIRnhl = vSFIRn[-1 + ihil];
      vSFIZnhl = vSFIZn[-1 + ihil];
      vSZFInhl = vSZFIn[-1 + ihil];
      vphl = vp[-1 + ihil];
      vrphl = vrp[-1 + ihil];
      vzphl = vzp[-1 + ihil];
      vdphl = vdp[-1 + ihil];
      vsphl = vsp[-1 + ihil];
      vhbphl = vhbp[-1 + ihil];
      vSRFIphl = vSRFIp[-1 + ihil];
      vSFIRphl = vSFIRp[-1 + ihil];
      vSFIZphl = vSFIZp[-1 + ihil];
      vSZFIphl = vSZFIp[-1 + ihil];
      dvnhl = dvn[-1 + ihil];
      dvphl = dvp[-1 + ihil];
      //!
      for (N1 = 1; N1 <= nd; N1++)
      {
        JA = IM + N1;
        nsa = ns[-1 + JA];
        SSU = std::max(nsa, 0);
        SSD = std::max(-nsa, 0);
        qhla = QHLA_opt[-1 + JA][-1 + ihil];
        fi1r = FI1R_opt[-1 + JA][-1 + ihil];
        fi1z = FI1Z_opt[-1 + JA][-1 + ihil];
        fi2d = FI2D_opt[-1 + JA][-1 + ihil];
        FIU[-1 + N1] = qhla * SSU;
        FIUR[-1 + N1] = fi1r * SSU;
        FIUZ[-1 + N1] = fi1z * SSU;
        FIUD2N[-1 + N1] = (fi2d - xlamy2 * qhla) * SSU;
        FID[-1 + N1] = qhla * SSD;
        FIDR[-1 + N1] = fi1r * SSD;
        FIDZ[-1 + N1] = fi1z * SSD;
        FIDD2N[-1 + N1] = (fi2d - xlapy2 * qhla) * SSD;
      } // End Do
      //!
      i = ibro;
      for (N1 = 1; N1 <= nd; N1++)
      {
        JA = IM + N1;
        nsa = ns[-1 + JA];
        FIUN1 = FIU[-1 + N1];
        FIURN1 = FIUR[-1 + N1];
        FIUZN1 = FIUZ[-1 + N1];
        FIUD2N1 = FIUD2N[-1 + N1];
        FIDN1 = FID[-1 + N1];
        FIDRN1 = FIDR[-1 + N1];
        FIDZN1 = FIDZ[-1 + N1];
        FIDD2N1 = FIDD2N[-1 + N1];
        for (N2 = 1; N2 <= N1; N2++)
        {
          i = i + 1;
          i1 = i + nhhdim;
          i2 = i + nhhdim2;
          i3 = i + nhhdim3;
          nsb = ns[-1 + N2 + IM];
          nsab = nsa + nsb;
          if (nsab != 0)
          {
            if (nsb > 0)
            {
              //! spin : UpUp
              //!  SRFIh=ZERO; SFIRh=ZERO; SZFIh=ZERO
              FIUN2 = FIU[-1 + N2];
              FIURN2 = FIUR[-1 + N2];
              FIUD2N2 = FIUD2N[-1 + N2];
              FIUZN2 = FIUZ[-1 + N2];
              vh = FIUN1 * FIUN2;
              hbh = vh * xlamy2 + FIURN1 * FIURN2 + FIUZN1 * FIUZN2;
              vdh = hbh + hbh + FIUN1 * FIUD2N2 + FIUN2 * FIUD2N1;
              SNABLARh = FIURN1 * FIUN2 + FIURN2 * FIUN1;
              SNABLAZh = FIUZN1 * FIUN2 + FIUZN2 * FIUN1;
              vsh = SNABLARh * xlamy;
              SFIZh = (vh + vh) * xlamy;
            }
            else
            //! spin : DoDo
            {
              //! SRFIh=ZERO; SFIRh=ZERO; SZFIh=ZERO
              FIDN2 = FID[-1 + N2];
              FIDRN2 = FIDR[-1 + N2];
              FIDZN2 = FIDZ[-1 + N2];
              FIDD2N2 = FIDD2N[-1 + N2];
              vh = FIDN1 * FIDN2;
              hbh = vh * xlapy2 + FIDRN1 * FIDRN2 + FIDZN1 * FIDZN2;
              vdh = hbh + hbh + FIDN1 * FIDD2N2 + FIDN2 * FIDD2N1;
              SNABLARh = FIDRN1 * FIDN2 + FIDRN2 * FIDN1;
              SNABLAZh = FIDZN1 * FIDN2 + FIDZN2 * FIDN1;
              vsh = -SNABLARh * xlapy;
              SFIZh = -(vh + vh) * xlapy;
            } // End If
            brout[-1 + i] = brout[-1 + i] + vSFIZnhl * SFIZh + vh * vnhl + SNABLARh * vrnhl + SNABLAZh * vznhl + vdh * vdnhl + vsh * vsnhl + hbh * vhbnhl;
            brout[-1 + i1] = brout[-1 + i1] + vSFIZphl * SFIZh + vh * vphl + SNABLARh * vrphl + SNABLAZh * vzphl + vdh * vdphl + vsh * vsphl + hbh * vhbphl;
            brout[-1 + i2] = brout[-1 + i2] + vh * dvnhl;
            brout[-1 + i3] = brout[-1 + i3] + vh * dvphl;
          }
          else
          {
            if (nsb > 0)
            { //! spin:DoUp
              // vh = ZERO hbh = zero vdh = zero SNABLARh = zero SNABLAZh = zero SFIZh = zero;
              FIUN2 = FIU[-1 + N2];
              FIURN2 = FIUR[-1 + N2];
              FIUD2N2 = FIUD2N[-1 + N2];
              FIUZN2 = FIUZ[-1 + N2];
              FITW3 = -FIDZN1 * FIUN2;
              FITW4 = FIUZN2 * FIDN1;
              vsh = -FIDRN1 * FIUZN2 + FIURN2 * FIDZN1 + FITW3 * xlamy - FITW4 * xlapy;
              SRFIh = -FIDRN1 * FIUN2 + FIURN2 * FIDN1;
              SFIRh = FIDN1 * FIUN2 * xlampy;
              SZFIh = FITW3 + FITW4;
            }
            else // !spin:UpDo
            {
              //! vh = ZERO;
              hbh = zero;
              vdh = zero;
              SNABLARh = zero;
              SNABLAZh = zero;
              SFIZh = zero;
              FIDN2 = FID[-1 + N2];
              FIDRN2 = FIDR[-1 + N2];
              FIDZN2 = FIDZ[-1 + N2];
              FIDD2N2 = FIDD2N[-1 + N2];
              FITW3 = -FIDZN2 * FIUN1;
              FITW4 = FIUZN1 * FIDN2;
              vsh = FIURN1 * FIDZN2 - FIDRN2 * FIUZN1 - FITW4 * xlapy + FITW3 * xlamy;
              SRFIh = FIURN1 * FIDN2 - FIDRN2 * FIUN1;
              SFIRh = FIUN1 * FIDN2 * xlampy;
              SZFIh = FITW3 + FITW4;
            } // Endif
            brout[-1 + i] = brout[-1 + i] + vsh * vsnhl + vSRFInhl * SRFIh + vSFIRnhl * SFIRh + vSZFInhl * SZFIh;
            brout[-1 + i1] = brout[-1 + i1] + vsh * vsphl + vSRFIphl * SRFIh + vSFIRphl * SFIRh + vSZFIphl * SZFIh;
          } // End If
          //!----------------------------------------------
          //! LN PH PART
          //!----------------------------------------------
          if (kindhfb < 0)
          {
            if (ihil == 1)
            {
              un = zero;
              up = zero;
              if (N1 == N2)
              {
                un = -ala2[-1 + 1];
                up = -ala2[-1 + 2];
              } // End If
              n12 = N1 + (N2 - 1) * nd;
              brout[-1 + i] = brout[-1 + i] + two * (ala2[-1 + 1] * rk[-1 + n12][-1 + ib] + un);
              brout[-1 + i1] = brout[-1 + i1] + two * (ala2[-1 + 2] * rk[-1 + n12][-1 + ibx] + up);
            } // End If
          }   // End If
        }     // End Do !N2
      }       // End Do !N1
    }         // End Do !ihil
    ibro = i;
  } // End Do !IB
  //
  std::cout << "End of gamdel()" << std::endl;
  // not consider for now ! TMR pairing matrix elements
  // not consider for now If(use_TMR_pairing.Ne.0) Then
  // not consider for now  !if(inin.lt.0.or.(inin.gt.0.and.iiter.gt.0))
  // not consider for now  Call delta
  // not consider for now Endif
  // not consider for now !
  //!----------------------------------------------
  //! BROYDEN/LINEAR MIXING
  //!----------------------------------------------
  //!
  // Call get_CPU_time('broyden',0)
  //!
  broyden_min(nhhdim4, brout, brin, alphamix, si, iiter, nbroyden, bbroyden);
  //!
  // Call get_CPU_time('broyden',1)
  //!
}

//!=======================================================================
void HFBTHO_solver::broyden_min(int &N, std::vector<double> &vout, std::vector<double> &vin, double &alpha,
                                double &si, int &iter, int &M, char &bbroyden)
{
  //!---------------------------------------------------------------------
  //! Modified Broyden's method: D.D.Johnson, PRB 38, 12807 (1988)
  //! Adopted from: (C) 2001 PWSCF group
  //! Input :
  //!  N      dimension of arrays vin,vout
  //!  vin    outpu at previous iteration
  //!  vout   output at current iteration
  //!  alpha  mixing factor (0 < alpha <= 1)
  //!  iter   current iteration number
  //!  M      number of iterations in Broyden history
  //!  M=0    Linear mixing
  //! Output:
  //!  si     MaxVal(|vout-vin|)
  //!  vin    Broyden/Linear mixing result
  //!  vout   vout-vin
  //!  bbroyden='B' Broyden mixing, curvature>0
  //!  bbroyden='L' Linear mixing,  curvature<0
  //!---------------------------------------------------------------------
  // Use HFBTHO, Only: pr,ipr,ierror_flag,ierror_info
  // Implicit None
  // Integer(ipr),     Intent(In)    :: N,iter,M
  // Real(pr),        Intent(In)     :: alpha
  // Real(pr),        Intent(Out)    :: si
  // Character(1),      Intent(Out)  :: bbroyden
  // Real(pr),        Intent(InOut)  :: vout(N),vin(N)
  int i, j, iter_used, ipos, inext;
  std::vector<int> iwork;
  std::vector<double> work, curv;
  std::vector<std::vector<double>> beta, df, dv;
  // Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
  // Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
  double w0;
  double DDOT, DNRM2, normi, gamma, curvature, sf;
  int qqone;
  //!
  sf = -1.0;
  qqone = 1.0;

  // DAXPY(N, sf, vin, 1, vout, 1);
  // qqdaxpy(N, sf, vin, qqone, vout, qqone);
  // need to check storage pattern for c++ and fortran
  cblas_daxpy(N, sf, &vin[0], qqone, &vout[0], qqone);
  // for (int qqi = 0; qqi < sizeof(vin) / sizeof(vin[0]); qqi++)
  //{
  //  vout[qqi] = vout[i] + sf * vin[qqi];
  //}
  // si = Maxval(Abs(vout));
  si = std::max(*std::max_element(vout.begin(), vout.end()), *std::min_element(vout.begin(), vout.end()));
  //! Linear mixing
  if (M == 0 || iter == 0)
  {
    bbroyden = 'L';
    // need to check storage pattern for c++ and fortran
    cblas_daxpy(N, alpha, &vout[0], qqone, &vin[0], qqone);
    //! If(iter.Eq.0) Write(6,*) '  Linear mixing (alpha) : ',alpha
    return;
  }
  //! Broyden mixing
  iter_used = std::min(iter - 1, M);
  ipos = iter - 1 - ((iter - 2) / M) * M;
  inext = iter - ((iter - 1) / M) * M;
  if (iter == 1)
  {
    w0 = 0.010;
    // if (Allocated(df))
    //   Deallocate(curv, df, dv, beta, work, iwork);
    // Allocate(curv(N), df(N, M), dv(N, M), beta(M, M), work(M), iwork(M));
    //! write(6,'(a,i3,3(2x,f18.8),a)') '   Broyden mixing (M,alpha,w0,mem) : '  &
    //!      ,M,alpha,w0,(2*N*M+N)*8._pr/1.e6,' MB'
    curv.resize(N);
    work.resize(N);
    iwork.resize(N);
    df.resize(N, std ::vector<double>(M));
    dv.resize(N, std ::vector<double>(M));
    beta.resize(M, std ::vector<double>(M));
  }
  else
  {
    // df( :, ipos) = vout( :) - df( :, ipos);
    // dv( :, ipos) = vin( :) - dv( :, ipos);
    for (int qqi = 1; qqi <= N; qqi++)
    {
      df[-1 + qqi][-1 + ipos] = vout[-1 + qqi] - df[-1 + qqi][-1 + ipos];
      dv[-1 + qqi][-1 + ipos] = vin[-1 + qqi] - dv[-1 + qqi][-1 + ipos];
    }

    // need to check storage pattern for c++ and fortran
    normi = 1.0 / std::sqrt(pow((cblas_dnrm2(N, &df[-1 + 1][-1 + ipos], 1)), 2));
    // need to check storage pattern for c++ and fortran
    cblas_dscal(N, normi, &df[-1 + 1][-1 + ipos], 1);
    // need to check storage pattern for c++ and fortran
    cblas_dscal(N, normi, &dv[-1 + 1][-1 + ipos], 1);
  }
  // Endif
  for (i = 1; i <= iter_used; i++)
  {
    for (j = i + 1; j <= iter_used; j++)
    {
      // need to check the storage pattern for c++ and fortran
      beta[-1 + i][-1 + j] = cblas_ddot(N, &df[-1 + 1][-1 + j], 1, &df[-1 + 1][-1 + i], 1);
    };
    beta[-1 + i][-1 + i] = 1.0 + w0 * w0;
  };
  // DSYTRF('U', iter_used, beta, M, iwork, work, M, i);
  LAPACK_dsytrf("U", &iter_used, &beta[0][0], &M, &iwork[0], &work[0], &M, &i);
  if (i != 0)
  {
    // ierror_flag=ierror_flag+1
    // ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRF '
    return;
  }
  // DSYTRI('U', iter_used, beta, M, iwork, work, i);
  LAPACK_dsytri("U", &iter_used, &beta[0][0], &M, &iwork[0], &work[0], &i);
  if (i != 0)
  {
    // ierror_flag=ierror_flag+1
    // ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRI '
    return;
  }
  for (i = 1; i <= iter_used; i++)
  {
    for (j = i + 1; iter_used; j++)
    {
      beta[-1 + j][-1 + i] = beta[-1 + i][-1 + j];
    }
    work[-1 + i] = cblas_ddot(N, &df[-1 + 1][-1 + i], 1, &vout[0], 1);
  }
  // curv = alpha * vout;
  for (int qqi = 1; qqi <= sizeof(vout) / sizeof(vout[0]); qqi++)
  {
    curv[qqi] = alpha * vout[qqi];
  }
  //
  for (i = 1; i <= iter_used; i++)
  {
    gamma = 0.0;
    for (j = 1; j <= iter_used; j++)
    {
      gamma = gamma + beta[-1 + j][-1 + i] * work[-1 + j];
    }
    // curv = curv - gamma * (dv( :, i) + alpha * df( :, i));
    for (int qqi = 1; qqi <= N; qqi++)
    {
      curv[-1 + qqi] = curv[-1 + qqi] - gamma * (dv[-1 + qqi][-1 + i] + alpha * df[-1 + qqi][-1 + i]);
    }
  }
  cblas_dcopy(N, &vout[0], 1, &df[-1 + 1][-1 + inext], 1);
  cblas_dcopy(N, &vin[0], 1, &dv[-1 + 1][-1 + inext], 1);
  curvature = cblas_ddot(N, &vout[0], 1, &curv[0], 1);
  if (curvature > -1.0)
  {
    bbroyden = 'B';
    sf = +1.0;
    cblas_daxpy(N, sf, &curv[0], 1, &vin[0], 1);
  }
  else
  {
    bbroyden = 'L';
    sf = alpha * 0.50;
    cblas_daxpy(N, sf, &vout[0], 1, &vin[0], 1);
  } // End If
  // End Subroutine broyden_min
  //!=======================================================================
  //!
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// self made tools.  can't find blas header so has to do it self.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HFBTHO_solver::qqdaxpy(int &N, double &sf, std::vector<double> &x, int &incx, std::vector<double> &y, int &incy)
{
  // incx and incy not used here, default as 1
  for (int qqi = 0; qqi < N; qqi++)
  {
    y[qqi] = y[qqi] + sf * x[qqi];
  }
}