// MatriceElementaire.hpp --- 
// 
// Filename: MatriceElementaire.hpp
// Description: 
// Author: E.NAYIR O.TISSOT
// Created: ven. janv. 24 17:13:23 2014 (+0100)
// Version: 
// Last-Updated: dim. mars 23 19:16:24 2014 (+0100)
//           By: Olivier
//     Update #: 242
// Compatibility: 
// 
// RNM.hpp, R1.hpp, R2.hpp

// Commentary: 
// 
// Calcul des matrices elementaires associees
// 
// 

// Code:

#ifndef MATRICELEMENTAIRE_HPP
#define MATRICELEMENTAIRE_HPP

#define CHEK_KN
#include "Maillage.hpp"    // pour avoir les classes Sommet et Triangle
#include "R1.hpp"
#include "R2.hpp"
#include "RNM.hpp"         // pour faciliter l'ecriture : classe de tableaux
#include "QuadratureFormular.hpp"
#include "P1P2Lagrange.hpp"
#include <cmath>
#include <cassert>

// La transformation affine
class Bk{
public:
  void static init(const Triangle * T, R b[2][2]) {
    const R2& s0 = *(T->s[0]);
    const R2& s1 = *(T->s[1]);
    const R2& s2 = *(T->s[2]);
    b[0][0] = s2.y - s0.y;
    b[0][1] = s0.x - s2.x;
    b[1][0] = s0.y - s1.y;
    b[1][1] = s1.x - s0.x;
  }
  double static det(const Triangle * T) {
    const R2& s0 = *(T->s[0]);
    const R2& s1 = *(T->s[1]);
    const R2& s2 = *(T->s[2]);
    return (s2.y - s0.y)*(s1.x - s0.x) - (s0.y - s1.y)*(s0.x - s2.x);
  }
  double static sgnDet(const Triangle * T) {
    return copysign(1,Bk::det(T));
  }
};

// Quelques operations sur les tableaux utiles pour faciliter l'ecriture
KNM<R> operator*(const double& alpha, const KNM<R>& A);

// On va definir les matrices elementaires sur le Kchapeau
void MatricesElementairesKchap(KNM<R>& Mchap,
			       KNM<R>& Bchap1, 
			       KNM<R>& Bchap2, 
			       KNM<R>& Wchap11, 
			       KNM<R>& Wchap12, 
			       KNM<R>& Wchap22);

void MatricesElementaires(const Triangle& K, 
			  KNM<R>& A, 
			  KNM<R>& B,
			  const KNM<R>& Bchap1,
			  const KNM<R>& Bchap2,
			  const KNM<R>& Wchap11,
			  KNM<R>& Wchap12,
			  const KNM<R>& Wchap22);

// Masse P2 avec Quadrature
double MelemQ(int, int);
double Wchap11Q(int, int);
double Wchap12Q(int, int);
double Wchap22Q(int, int);
double Bchap1Q(int, int);
double Bchap2Q(int, int);
#endif


// 
// MatriceElementaire.hpp ends here
