#include <fstream>

#include "HFBTHO_solver.h"
#include "printArray.h"

int main()
{
    HFBTHO_solver hfbtho_solver;

    // qq defined variables to interface with corresponding ones in HFBTHO.h
    int qqn00_INI;
    int qqiLST_INI;
    int qqinin_INI;
    int qqicou_INI;

    // int npr_INI[3], kindhfb_INI;
    int qqnpr_INI[3] = {0};
    int qqkindhfb_INI;

    // int keyblo1_INI, keyblo2_INI, IDEBUG_INI;
    int qqkeyblo1_INI;
    int qqkeyblo2_INI;
    int qqIDEBUG_INI;

    // double b0_INI, q_INI, cdef_INI, cqad_INI;
    double qqb0_INI;
    double qqq_INI;
    double qqcdef_INI;
    double qqcqad_INI;

    char qqskyrme_INI[30];
    double qqepsi_INI;
    // bool Add_Pairing_INI, Print_HFBTHO_Namelist_INI, DO_FITT_INI;
    bool qqAdd_Pairing_INI;
    bool qqPrint_HFBTHO_Namelist_INI;
    bool qqDO_FITT_INI;

    bool lpr;
    int iblocase[2] = {0}, nkblocase[2][5] = {0};
    int i, it, iodd, ios, envstat, envlen;
    int qrpa_points_thiscore, qrpa_start_thiscore, idivm;
    int arrayTaskID, arrayTaskCount, arrayTaskMin, arrayTaskMax;
    char envRead[20];
    int testarray[3][7][8] = {0};

    std::fstream thoRead("tho.dat");
    if (thoRead.is_open())
    {
        thoRead >> qqn00_INI >> qqb0_INI >> qqq_INI;
        thoRead >> qqnpr_INI[0] >> qqnpr_INI[1];
        thoRead >> qqskyrme_INI;
        thoRead >> qqkindhfb_INI >> qqinin_INI >> qqcdef_INI >> qqcqad_INI;
        thoRead >> nkblocase[0][0] >> nkblocase[0][1] >> nkblocase[0][2] >> nkblocase[0][3] >> nkblocase[0][4];
        thoRead >> nkblocase[1][0] >> nkblocase[1][1] >> nkblocase[1][2] >> nkblocase[1][3] >> nkblocase[1][4];

        HFBTHO::n00_INI = qqn00_INI;
        HFBTHO::iLST_INI = qqiLST_INI;
        HFBTHO::inin_INI = qqinin_INI;
        HFBTHO::icou_INI = qqicou_INI;

        // int npr_INI[3], kindhfb_INI;
        // int npr_INI[3] = qqnpr_INI[3];
        for (int i = 0; i < 3; i++)
            HFBTHO::npr_INI[i] = qqnpr_INI[i];

        HFBTHO::kindhfb_INI = qqkindhfb_INI;

        // int keyblo1_INI, keyblo2_INI, IDEBUG_INI;
        HFBTHO::keyblo1_INI = qqkeyblo1_INI;
        HFBTHO::keyblo2_INI = qqkeyblo2_INI;
        HFBTHO::IDEBUG_INI = qqIDEBUG_INI;

        // double b0_INI, q_INI, cdef_INI, cqad_INI;
        HFBTHO::b0_INI = qqb0_INI;
        HFBTHO::q_INI = qqq_INI;
        HFBTHO::cdef_INI = qqcdef_INI;
        HFBTHO::cqad_INI = qqcqad_INI;

        // char skyrme_INI[30] = qqskyrme_INI[30];
        for (int i = 0; i < 30; i++)
            HFBTHO::skyrme_INI[i] = qqskyrme_INI[i];

        HFBTHO::epsi_INI = qqepsi_INI;
        // bool Add_Pairing_INI, Print_HFBTHO_Namelist_INI, DO_FITT_INI;
        HFBTHO::Add_Pairing_INI = qqAdd_Pairing_INI;
        HFBTHO::Print_HFBTHO_Namelist_INI = qqPrint_HFBTHO_Namelist_INI;
        HFBTHO::DO_FITT_INI = qqDO_FITT_INI;

        hfbtho_solver.read_UNEDF_NAMELIST();
        hfbtho_solver.read_HFBTHO_NAMELIST();
        hfbtho_solver.printUNEDF();

        // printArray1D<int>(qqnpr_INI, 3);
        // printArray1D<int>(iblocase, 2);
        // printArray2D<int>(nkblocase, 2, 5);
        // printArray2D<int>(testarray[0], 7, 8);

        std::cout << HFBTHO::n00_INI << HFBTHO::b0_INI << HFBTHO::q_INI << std::endl;
        std::cout << HFBTHO::npr_INI[0] << HFBTHO::npr_INI[1] << HFBTHO::npr_INI[2] << std::endl;
        std::cout << HFBTHO::skyrme_INI << qqkindhfb_INI << qqinin_INI << qqcdef_INI << qqcqad_INI << std::endl;
        std::cout << nkblocase[0][0] << nkblocase[0][1] << nkblocase[0][2] << nkblocase[0][3] << nkblocase[0][4] << std::endl;
        std::cout << nkblocase[1][0] << nkblocase[1][1] << nkblocase[1][2] << nkblocase[1][3] << nkblocase[1][4] << std::endl;

        std::cout << "print inin_INI in main" << HFBTHO::inin_INI << std::endl;
        hfbtho_solver.printHFBTHO();
        // nkblo_INI initialized to 0;
    }

    thoRead.close();
    // use brute way to hard implement unedf_namelist.dat first
    //!================================= HFBTHO_NAMELIST =========================================
    //&HFBTHO_NAMELIST
    HFBTHO::MAX_ITER_INI = 201;
    HFBTHO::epsi_INI = 1.000E-08;
    HFBTHO::Add_Pairing_INI = true;
    HFBTHO::icou_INI = 2;
    //! HFBTHO:: ILST_INI =         -1,
    HFBTHO::iLST_INI = 0;
    HFBTHO::keypj_INI = 11;
    HFBTHO::iproj_INI = 1;
    HFBTHO::npr1pj_INI = 0;
    HFBTHO::npr2pj_INI = 0;
    HFBTHO::DO_FITT_INI = false;
    HFBTHO::IDEBUG_INI = 0;
    HFBTHO::Parity_INI = true;
    HFBTHO::Print_HFBTHO_Namelist_INI = true;
    //
    // use brute way to hard implement skyrme parameters
    //
    //    ! JASON 15 September
    HFBTHO::DMEorder = -1;
    HFBTHO::DMElda = 0;
    HFBTHO::use_cm_cor = false;
    HFBTHO::use_TMR_pairing = 0; // HBZERO = 20.735530D0;
    HFBTHO::e2charg = 1.439978400000000;
    HFBTHO::Crho[0] = -731.222785829509800000;
    HFBTHO::Cdrho[0] = 855.690051584978500000;
    HFBTHO::Ctau[0] = -0.543988860905982100;
    HFBTHO::Crho[1] = 263.710305524676100000;
    HFBTHO::Cdrho[1] = -176.864195604041100000;
    HFBTHO::Ctau[1] = -33.361881866521340000;
    HFBTHO::CrDr[0] = -43.29008975531540;
    HFBTHO::CJ[0] = 0.000000000000000000;
    HFBTHO::CrdJ[0] = -75.260870048289410000;
    HFBTHO::CrDr[1] = -164.1379857135440;
    HFBTHO::CJ[1] = 0.000000000000000000;
    HFBTHO::CrdJ[1] = -22.652819964871300000;
    HFBTHO::CpV0[0] = -186.192246596249;
    HFBTHO::CpV1[0] = 0.5000000000000000;
    HFBTHO::sigma = 0.298783982778235700;
    HFBTHO::CpV0[1] = -206.746416898386;
    HFBTHO::CpV1[1] = 0.500000000000000;
    HFBTHO::CExPar = 0.639129523762364;
    // RHO_NM = 0.157327162963946800;
    // E_NM = -15.800000484870580000;
    // K_NM = 225.943394883389600000;
    // SMASS_NM = 0.995872580822852000;
    // ASS_NM = 28.348338556986590000;
    // LASS_NM = 40.001962979089330000;
    // VMASS_NM = 1.248999953269958000;

    hfbtho_solver.solver();
}
