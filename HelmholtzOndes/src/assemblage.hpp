#ifndef ASSEMBLAGE_HPP
#define ASSEMBLAGE_HPP

#include "mat_cr.hpp"
#include "maillage.hpp"
#include "vect.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <stdlib.h>

void build_mat(matrice&, matrice&, const vector< point >&, const vector< edge >&, const vector< triangle >&, int, char**);
void pseudo_elim(matrice&, const vector< edge >&);
void build_f(vect&, const vector< edge >&, const vector< point >&, double, int=1);
double g(double, double, double, vect&);

#endif
