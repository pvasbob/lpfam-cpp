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
  void Constraint_or_not(int& inin_INI,int& inin,int& icstr);  

private:
  int iw, ib, j, i, it;
  double epsi0;
};

#endif
