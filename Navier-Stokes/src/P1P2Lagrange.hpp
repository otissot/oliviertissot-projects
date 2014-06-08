// P1P2Lagrange.hpp --- 
// 
// Filename: P1P2Lagrange.hpp
// Description: 
// Author: O.TISSOT 
// Created: dim. mars 23 00:28:38 2014 (+0100)
// Version: 
// Last-Updated: dim. mars 23 08:57:48 2014 (+0100)
//           By: Olivier
//     Update #: 26
// Compatibility: 
// 
// 

// Commentary: 
// 
// Petites classes d'enrobage pour implementer les fonctions de bases associees
// aux elements P1, P2 Lagrange
// 


// Code:

#ifndef P1P2LAGRANGE_HPP
#define P1P2LAGRANGE_HPP

#include "Maillage.hpp"

// Fonctions de base P1
class P1 {
public:
  // lambda
  void static init(const R2 & P, R * l) {
    l[0] = 1 - P.x - P.y;
    l[1] = P.x;
    l[2] = P.y;
  }
  void static init(const R2 & P, const Triangle * T, R * l) {  
    const R2 Q[3]={ R2((*T)[0]),R2((*T)[1]),R2((*T)[2]) };
    l[0] = det(P,Q[1],Q[2]);
    l[1] = det(Q[0],P,Q[2]);
    l[2] = det(Q[0],Q[1],P);
    R Det = l[0]+l[1]+l[2];
    l[0] /= Det;
    l[1] /= Det;
    l[2] /= Det;
  }
  // grad(lambda)
  void static init(R2 * dl) {
    dl[0] = R2(-1,-1);
    dl[1] = R2(1,0);
    dl[2] = R2(0,1);
  }
  void static init(const Triangle * T, R2 * dl) {
    const R2& s0 = *(T->s[0]);
    const R2& s1 = *(T->s[1]);
    const R2& s2 = *(T->s[2]);
    
    dl[0] = R2(s1,s2).perp()/(2.*(T->mesure()));
    dl[1] = R2(s2,s0).perp()/(2.*(T->mesure()));
    dl[2] = -dl[0]-dl[1];
  }
};

// Fonctions de base P2
class P2 {
public:
  // phi
  // Une alloc-desalloc pour alleger l'appel (mais legerement plus lent du coup)
  // Ne pas utiliser cette methode...
  void static init(const R2& P, R * phi) {
    R l[3];
    P1::init(P,l);
    P2::init(l,phi);
  }  
  void static init(const R2& P, const Triangle * T, R * phi) {
    R l[3];
    P1::init(P,T,l);
    for(int ii=0;ii<3;++ii){
      phi[ii] = l[ii]*(2*l[ii]-1);
      phi[ii+3] = 4*l[(ii+1)%3]*l[(ii+2)%3];
    }
    //if(l) delete[] l;
  }
  // A partir des coordonnees barycentriques (fonctions P1)
  void static init(R* l, R* phi){
    for(int ii=0;ii<3;++ii){
      phi[ii] = l[ii]*(2*l[ii]-1);
      phi[ii+3] = 4*l[(ii+1)%3]*l[(ii+2)%3];
    }
  }
  // grad(phi)
  void static init(const R2& P, R2 * dphi) {
    dphi[0].x = dphi[0].y = 1-4*(1-P.x-P.y);
    dphi[1] = R2(4*P.x-1,0);
    dphi[2] = R2(0,4*P.y-1);
    dphi[3] = R2(4*P.y,4*P.x);
    dphi[4] = R2(-4*P.y,4-4*P.x-8*P.y);
    dphi[5] = R2(4-4*P.y-8*P.x,-4*P.x);
  }
  // bug
  /*
  void static init(const R2& P, const Triangle * T, R2 * dphi) {
    R l[3];
    R2 dl[3];
    P1::init(P,T,l);
    P1::init(T,dl);
    // derivee du produit
    for(int ii=0;ii<3;++ii){
      dphi[ii] = l[ii]*2*dl[ii] + dl[ii]*(2*l[ii]-1);
      dphi[ii+3] = 4*dl[(ii+1)%3]*l[(ii+2)%3] + 4*l[(ii+1)%3]*dl[(ii+2)%3];
    }
  }
  */
};

#endif

// 
// P1P2Lagrange.hpp ends here
