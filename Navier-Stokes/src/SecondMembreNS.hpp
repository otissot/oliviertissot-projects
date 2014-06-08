// SecondMembreNS.hpp --- 
// 
// Filename: SecondMembreNS.hpp
// Description: 
// Author: O.TISSOT 
// Created: lun. mars  3 19:02:46 2014 (+0100)
// Version: 
// Last-Updated: dim. mars 23 01:05:08 2014 (+0100)
//           By: Olivier
//     Update #: 34
// Compatibility: 
// 
// 

// Commentary: 
// 
// Les fonctions necessaires a la construction du terme
// u_n°X_n
// 

// Code:

#ifndef SECONDMEMBRENS_HPP
#define SECONDMEMBRENS_HPP

#include "Maillage.hpp"
#include "P1P2Lagrange.hpp"

// Projection P2 pour calculer une valeur de la vitesse a l'interieur du triangle
R2 Projection(const Maillage &,
	      map< pair<int,int>, int > &,
	      const Triangle *,
	      R *,
	      const KN_<R> &,
	      const KN_<R> &);

// Calcule le transporte du point
R2 Transport(const Maillage & Th,
	     int it, 
	     double *lambda,
             const  KN_<R> & U,
	     const  KN_<R> & V, 
	     const double & dt);


// Assemble le second membre
void buildSecondMembreFind(const Maillage &,
		       map< pair<int,int>, int > &,
		       const KN_<R> &,
		       const KN_<R> &,
		       const KN_<R> &,
		       const double &);

double g(const R2&);
double g(double);

int onBorder(R *, const Triangle *);

void buildSecondMembre(const Maillage &,
		       const KN_<R> &,
		       const KN_<R> &,
		       const KN_<R> &,
		       const double,
		       int,
		       map< pair<int,int>, int > &);


// Des methodes pour projeter sur le maillage
// - on se promene pour trouver les points proches
// - on interpole grace aux coordonnees barycentriques
int  WalkInTriangle(const Maillage &,
		    int, 
		    double*,
		    R, 
		    R, 
		    R &);       
int  WalkInTriangle(const Maillage &,
		    int, 
		    double*,
                    const  KN_<R> &,
		    const  KN_<R> &, 
		    R &);
int Walk(const Maillage &,
	 int&, 
	 R*,
         const KN_<R>  &,
	 const KN_<R>  &, 
	 R);


#endif


// 
// SecondMembreNS.hpp ends here
