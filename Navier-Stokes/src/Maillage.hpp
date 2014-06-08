// Maillage.hpp --- 
// 
// Filename: Maillage.hpp
// Description: 
// Author: O.TISSOT E.NAYIR 
// Created: ven. déc. 13 10:00:43 2013 (+0100)
// Version: 
// Last-Updated: dim. mars 23 01:05:37 2014 (+0100)
//           By: Olivier
//     Update #: 241
// Compatibility: 
// 
// GQuadTree.hpp, R1.hpp, R2.hpp

// Commentary: 
// 
// Une classe qui definit notre maillage triangulaire 2D.
// Nous avons utilise comme modele la classe qui realise le meme
// travail pour Freefem++ : Mesh2dn.
// Notre classe est moins generique mais aussi moins compliquee.

// Code:

#ifndef MAILLAGE_HPP
#define MAILLAGE_HPP

#include "R1.hpp"
#include "R2.hpp"
#include "HashTable.hpp"
#include "RNM.hpp"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <string>
#include <map>

using namespace std;

inline int BinaryRand()
{
#ifdef RAND_MAX
  const long HalfRandMax = RAND_MAX/2;
  return rand()<HalfRandMax;
#else
  return rand() & 16384; // 2^14
#endif
};


class Sommet :public R2   
{
public:
  // Attributs :
  typedef R2 Rd;
  int numero;
  int label;         // pour savoir si c'est un point du bord et de quel bord
  // Constructeur :
  Sommet() : numero(-1) {};
  //Methodes :
  bool onGamma(int gammai=1) const {return label==gammai;}
  operator int () const {return numero;} // conversion sommet en entier  (**)
private: // Pas de copie 
  Sommet(const Sommet &);
  Sommet & operator=(Sommet &); 
};

class Segment {
public:
  // Attributs :
  Sommet * s[2];
  double mes;              // meme remarque que pour Triangle
  const static int nv = 2; // nombre de sommets
  // Constructeur :
  Segment() { s[0]=s[1]=0; }
  //~Segment() { if(s) delete[] s;}
  // Methode :
  Sommet & operator[](int i) {return *(s[i]);}
  const Sommet & operator[](int i) const {return *(s[i]);}
  double mesure() const { return mes; }
  bool contient(const R2& P) const{
    const R2& A = *s[1] - *s[0];
    const R2& B = P - *s[0];
    if((A^B)!=0) return false;
    else{
      return ( 0 < (A,B) && (A,B) < A.norme2() );
    }
  }
  void init(Sommet *v0, int n[2])
  {
    for(int i=0;i<=1;++i)
      s[i] = & v0[n[i]]; 
    mes = R2( *(s[0]),*(s[1]) ).norme();
    assert( mes > 0) ; 
  }
private: // Pas de copie
  Segment(const Segment &); 
  Segment & operator=(Segment &); 
};

class Triangle {
public:
  // Attributs :
  Sommet * s[3];// tableau de 3 pointeurs sur des Sommets    
  double mes;   // pour eviter de faire le calcul de la mesure plusieurs fois : on le fait juste dans init et on sauve la valeur
  const static int nva = 2;             // nb pt sur hypersurface
  const static int nea = 3;             // nb triangles voisins
  static const int (* const nvadj)[2];  // les segments du triangle
  // Constructeur :
  Triangle() {s[0]=s[1]=s[2]=0;}
  //~Triangle() {if(s) delete[] s;} 
  // Methodes :
  Sommet & operator[](int i) {return *(s[i]);}
  const Sommet & operator[](int i) const {return *(s[i]);}
  // Application affine : Khat -> K (exemple : milieu "bas" de K = R2 M(K(R2(1/.2,0))) et on obtient facilement les autres milieux de la même manière)
  R2 operator()(const R2& Shat) const {
    const R2& s0 = *s[0];
    const R2& s1 = *s[1];
    const R2& s2 = *s[2];
    return (1 - Shat.x - Shat.y)*s0 + Shat.x*s1 + Shat.y*s2;
  }
  // Application inverse : K -> Khat
  R2 fK(const R2& S) const{
    R2 res;
    const R2& s0 = *s[0];
    const R2& s1 = *s[1];
    const R2& s2 = *s[2];
    res.x = (s2.y-s0.y)*(S-s0).x - (s2.x-s0.x)*(S-s0).y;
    res.y = (s0.y-s1.y)*(S-s0).x - (s1.x-s0.x)*(S-s0).y;
    res *= 2*mes;
    return res;
  }
  // Si le triangle contient le point ou non
  bool contient(const R2& P) const{
    const R2& PA = P - *s[0];
    const R2& PB = P - *s[1];
    const R2& PC = P - *s[2];
    return ( ((PB^PC) > 0 && (PC^PA) > 0 && (PA^PB) > 0) || ((PB^PC) < 0 && (PC^PA) < 0 && (PA^PB) < 0) );
  }
  double mesure() const { return mes; }
  double mesure(Sommet *v[3]) {
    return det(*v[0],*v[1],*v[2])/2.;  
  }
  // Allocation de la memoire
  void init(Sommet *v0,int n[3])
  {
    for(int i=0;i<=2;++i)
      s[i] = & v0[n[i]]  ; 
    mes = mesure(s);
    assert( mes >0) ; 
  }
private: // Pas de copie 
  Triangle(const Triangle &); 
  Triangle & operator=(Triangle &); 
};

class Maillage 
{
public:
  // Attributs :
  int nv;                                  // nombre de sommets
  int nt;                                  // nombre de triangles
  int ne;                                  // nombre de segments au bord
  Triangle * t;                            // tableau des triangles
  Segment * be;                            // tableau des segments au bord
  Sommet * s;                              // tableau des sommets
  int *TheAdjacencesLink;                  // pour connaitre les triangles adjacents
  int *BoundaryElementHeadLink;            // idem sur le bord
  static const int nea = Triangle::nea;    // NbOfAdjElem (pour un triangle)
  static const int nva = Triangle::nva;    // NbOfVertexOnHyperFace (pour un triangle)
  static int ksearch, kthrough;            // Nb re recherches et nb d'elements traverses
  // Constructeur/Destructeur :
  Maillage(const string& filename);
  ~Maillage(){
    if(t) delete[] t;
    if(be) delete[] be;
    if(s) delete[] s;
    if(TheAdjacencesLink) delete[] TheAdjacencesLink;
  }
  // Methodes :
  int CheckT(int k) const { assert(k>=0 && k<nt); return k;}
  int CheckS(int k) const { assert(k>=0 && k<nv); return k;}
  int CheckBE(int k) const { assert(k>=0 && k<ne); return k;} 
  Triangle & operator[](int k) const { return t[CheckT(k)];}
  Segment & borderEdge(int k) const { return be[CheckBE(k)];}
  Sommet & operator()(int k) const { return s[CheckS(k)];}     // (***)
  int operator()(int k,int j) const { return CheckS(t[k][j]);} // Attention on utilise le cast du sommet en un entier voir ligne (**)

  // Pour avoir les numerotations associees (petit trick : on fait des soustractions sur des pointeurs -> renvoie un entier !!)
  int operator()(const Triangle & tri) const { return CheckT(&tri - t); }
  int operator()(const Triangle * tri) const { return CheckT(tri - t); }
  int operator()(const Sommet & so) const { return CheckS(&so - s); }
  int operator()(const Sommet * so) const { return CheckS(so - s); }
  int operator()(const Segment & seg) const { return CheckBE(&seg - be); }
  int operator()(const Segment * seg) const { return CheckBE(seg - be); }

  // Pour construire les segments qui ne sont pas au bord
  void buildEdge(map< pair<int,int>, int >&) const;
  
  // Pour contruire les triangles adjacents
  void BuildAdj();
  
  // Pour recuperer l'element adjacent avec numerotation induite
  int ElementAdj(int k,int &j) const  {
    int p=TheAdjacencesLink[nea*k+j];
    j=p%nea;
    return p>=0 ? p/nea: -1;
  }
  
  // Pour retrouver le triangle qui contient (ou le plus proche) de P
  const Triangle* Find(R2, R2&, bool&, const Triangle* =0) const;
  
  template<int N,int M>
  SortArray<int,N> iteme(const int (* const  nu )[N],int k,int i){
    int nnv[N];
    Triangle & K(t[CheckT(k)]);
    assert(i>=0 && i <M);
    for (int j=0;j<N;++j){
      nnv[j] = operator()(K[nu[i][j]]);
    }
    return SortArray<int,N>(nnv);
  }
  
  SortArray<int,Segment::nv> itemadj(int k,int i){
    return iteme<Segment::nv,Triangle::nea>(Triangle::nvadj,k,i);
  }
 
  SortArray<int,Segment::nv> itembe(int k){
    int nnv[Segment::nv];
    Segment & K(be[CheckBE(k)]);
    for (int j=0;j<Segment::nv;++j){
      nnv[j] = operator()(K[j]);
    }
    return SortArray<int,Segment::nv>(nnv);
  }

  
  // Pour avoir une correspondance entre numerotation globale et locale
  int numedge(int, map< pair<int,int>, int >&, int) const;
  int numglob(int, map< pair<int,int>, int >&, int) const;

private:                                                 // Pas de copie  (dur a programme, y reflechir).
  Maillage(const Maillage &); 
  Maillage & operator=(Maillage &); 
}; 


#endif


// 
// Maillage.hpp ends here
