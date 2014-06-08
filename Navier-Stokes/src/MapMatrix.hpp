// MapMatrix.hpp --- 
// 
// Filename: MapMatrix.hpp
// Description: 
// Author: O.TISSOT 
// Created: mar. déc.  3 09:34:47 2013 (+0100)
// Version: 
// Last-Updated: mar. févr. 11 11:12:05 2014 (+0100)
//           By: Olivier
//     Update #: 103
// Compatibility: 
// 
// Compatible avec le format CSR grace a la methode de conversion -> compatibilite UMFPACK (grace a la transposee !)

// Commentary: 
// 
// Implementation d'une matrice creuse sous la forme de map< pair<int,int> , matrixresult >. Pas le plus rapide
// mais l'ecriture est commode. A ne pas utiliser pour une matrice pleine.
// 

// Code:

#ifndef MAPMATRIX_HPP
#define MAPMATRIX_HPP

#include <iostream>
#include <map>
#include <stdexcept>
#include "R2.hpp"
#include "Maillage.hpp"

using namespace std;

// Une petite classe d'enrobage pour essayer de regler le pb de l'operator()(int,int)
class matrixresult{
public:
  // Attribut
  double value;
  // Constructeur vide
  //matrixresult(){cout << "Constructeur vide matrixresult !" << endl;}
  // Copie
  matrixresult(const double& val=0.):value(val){/*if(value==0){cout<<"Destruction instantannee !"<<endl; this->~matrixresult();}*/}
  matrixresult& operator=(const double& val){
    if(val!=0){
      value = val;
      //cout << "Affectation matrixresult from double : " << val << endl; 
    }
    return *this;
  }
  matrixresult(const matrixresult& mr):value(mr.value){}
  matrixresult& operator=(matrixresult& mr){
    if(this!=&mr && mr.value!=0){
      value = mr.value;
      //cout << "Affectation matrixresult : " << mr.value << endl; 
    }
    return *this;
  }
  // Cast a partir d'un double
  operator double() const {
    //cout << "Cast matrixresult !" << endl;
    return value;
  }
  // Operateur +=
  matrixresult& operator+=(const double& inc){
    value += inc;
    return *this;
  }
  // Operateur *=
  matrixresult& operator*=(const double& inc){
    value *= inc;
    return *this;
  }
};

class MapMatrix : map< pair<int,int> , double > {
private:
  // La taille de la matrice
  unsigned int dimension;
public:
  using map< pair<int,int> , double >::size;
  using map< pair<int,int> , double >::erase;
  // Constructeur de la matrice nulle + surcharge du constructeur vide :
  MapMatrix(unsigned int n=0) : dimension(n) {};
  void Identite();
  int getDimension() { return dimension; };
  /*
  // Operateur d'acces en lecture/ecriture dans un element de la matrice : complexite en log(n) (a ne pas utiliser de maniere trop intensive)
  // /!\ les indices commencent à 0 : A(0,0) est tout en haut à gauche de la matrice /!\
  */
  double& operator()(int,int);
  double operator()(int,int) const;
  // Converti la map vers le format CSR ligne (et pas colonne attention) avec une compatibilite vers le C (on utilise des pointeurs)
  // Pas d'allocation memoire !!
  int toCSR(int*,int*,double*);

  // Pour affichage des coefs non nuls et de la dimension directement
  friend ostream& operator<<(ostream&, const MapMatrix&);
};


#endif

// 
// MapMatrix.hpp ends here
