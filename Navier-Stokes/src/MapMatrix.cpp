// MapMatrix.cpp --- 
// 
// Filename: MapMatrix.cpp
// Description: 
// Author: O.TISSOT 
// Created: mar. déc.  3 10:11:14 2013 (+0100)
// Version: 
// Last-Updated: jeu. févr.  6 11:19:57 2014 (+0100)
//           By: Olivier
//     Update #: 238
// Compatibility: 
// 
// 

// Commentary: 
// 
// Implementation du .hpp

// Code:

#include "MapMatrix.hpp"

/*
// /!\ Si on fait s'afficher une valeur nulle de la matrice alors la taille de la map est augmentee /!\
*/

double& MapMatrix::operator()(int i, int j){
  return (*this)[make_pair(i,j)];
}

/*
// Idem ci-dessus mais appelle une exception si l'element n'existe pas donc on renvoit zero dans ce cas-la.
// La difference entre les 2 c'est qu'ici la taille de la map n'est pas modifiee contrairement a ce qui se
// passe si on utilise les []. Donc la methode est bien const et on peut l'utiliser avec une reference
// constante sur l'objet.
*/

double MapMatrix::operator()(int i, int j) const{
  try{
    return this->at(make_pair(i,j));
  }
  catch(const std::out_of_range &){
    return 0;
  }
}

void MapMatrix::Identite(){
  for(unsigned int i=0;i<dimension;++i)
    operator()(i,i) = 1.;
}

int MapMatrix::toCSR(int* row, int* col, double* val){
  unsigned int kk =0; 
  row[0] = 0;
  //  set unset value to zero ...
  for(int l=1;l<getDimension()+1;l++)
    row[l] = -1;
  double aij;
  unsigned int i;
  unsigned int j;
  //int countn = 0;
  for(map< pair<int,int>, double >::const_iterator k = begin(); k != end(); ++k){
    i = k->first.first;
    j = k->first.second;
    aij = k->second;
    //cout << "(" << i << "," << j << ") : " << aij << endl;
    col[kk] = j;
    val[kk] = aij;
    /*if(aij == 0.){ 
      cout << "(i,j) = (" << i << "," << j << ")" << endl;
      countn ++;
      cout << countn << "eme valeur nulle..." << endl;
      }*/
    row[i+1] = ++kk;
  }
  int p = 0; 
  for(i = 1; i < dimension+1; ++i)
    if(row[i] == -1) row[i] = p; 
    else p = row[i];
  return size();
}

ostream& operator<<(ostream& os, const MapMatrix& A){
  unsigned int i,j;
  for(i=0;i<A.dimension;i++){
    for(j=0;j<A.dimension;j++)
      os << A(i,j) << '\t';
    if(i!=A.dimension-1) os << endl;    // pour etre consistant avec le standard C++
  }
  return os;
}


// 
// MapMatrix.cpp ends here
