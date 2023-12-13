#include <iostream>
#include <cstdlib>
#include <cstring>
#include <iomanip>

#include <algorithm>
#include <bits/stdc++.h>
#include <ctime>

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
  //!
  //!------------------------------------------------------------------
  //! Loop the internal normalization
  //!------------------------------------------------------------------
  sitest = std::max(std::min(0.10, si * 0.010), 0.000010);
  norm_to_improve = true;
  inner[-1 + it] = -1;
  sumnz[-1 + it] = one;
}
