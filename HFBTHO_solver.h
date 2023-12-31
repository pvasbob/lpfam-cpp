#ifndef SOLVEHFBTHO_H
#define SOLVEHFBTHO_H

#include <string>

#include "HFBTHO.h"

class HFBTHO_solver : public HFBTHO
{
public:
  void solver();
  void iniialize_HFBTHO_SOLVER();
  void set_functional_parameters(std::string fname, bool lpr);
  void calculate_natural_units();
  void print_functional_parameters(int lout);
  void t_from_C();
  void Constraint_or_not(int &inin_INI, int &inin, int &icstr);
  void preparer(bool lpr);
  void gfv();
  void base0(bool lpr);
  void ord(int, std::vector<double> &);
  void thoalloc();
  void gausspoints();
  double DGAMMA(double);
  void gaussq(const int &KINDI, const int &N, const double &al, const double &be, const int &kpts, double *, double *, double *, double *);
  void Class(const int &kindi, const int &N, const double &ALPHA, const double &BETA, double *B, double *A, double &MUZERO);
  void GBTQL2(const int &, double *, double *, double *, int &);
  double GBSLVE(const double &, const int &, double *, double *);
  void base(bool);
  void gaupol(bool);
  void coordinateLST(bool);
  void optHFBTHO();
  void iter(bool);
  void hfbdiag(int, int);
  void Alambda(double, int, int);
  int inout(int);
  void start();
  void gamdel();
  void broyden_min(int &nhhdim4, std::vector<double> &brout, std::vector<double> &brin, double &alphamix,
                   double &si, int &iiter, int &nbroyden, char &bbroyden);
  void qqdaxpy(int &, double &, std::vector<double> &, int &, std::vector<double> &, int &);

private:
  int iw, ib, j, i, it;
  double epsi0;
};

#endif
