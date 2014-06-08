// MatriceElementaire.cpp --- 
// 
// Filename: MatriceElementaire.cpp
// Description: 
// Author: E.NAYIR O.TISSOT 
// Created: ven. janv. 31 15:07:42 2014 (+0100)
// Version: 
// Last-Updated: dim. mars 23 20:43:39 2014 (+0100)
//           By: Olivier
//     Update #: 333
// Compatibility: 
// 
// 

// Commentary: 
// 
// Implementation du .hpp du meme nom
// 
// 

// Code:

#include "MatriceElementaire.hpp"
// Pour forcer la precision du cout
//#include <iomanip>
//#include <string> 

// Quelques operations sur les tableaux utiles pour faciliter l'ecriture
KNM<R> operator*(const double& alpha, const KNM<R>& A){
  KNM<R> res(A);
  for(int i=0;i<res.N();++i)
    for(int j=0;j<res.M();++j)
      res(i,j) *= alpha;
  return res;
}

// On va definir les matrices elementaires sur le Kchapeau
void MatricesElementairesKchap(KNM<R>& Mchap,
			       KNM<R>& Bchap1, 
			       KNM<R>& Bchap2, 
			       KNM<R>& Wchap11, 
			       KNM<R>& Wchap12, 
			       KNM<R>& Wchap22){  
  // Initialisation des tableaux
  for(int i=0;i<6;++i){
    for(int j=0;j<6;++j){
      Mchap(i,j) = 0.;
      Wchap11(i,j) = 0.;
      Wchap12(i,j) = 0.;
      Wchap22(i,j) = 0.;
      if(j<3){ Bchap1(i,j) = 0.; Bchap2(i,j) = 0.; }
    }
  }

  // Matrice de masse P2
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      Mchap(i,j) = -1 + 7*(i==j);
      Mchap(i+3,j) = Mchap(i,j+3) = -4*(i==j);
      Mchap(i+3,j+3) = 16*(1+(i==j));
    }
  }
  Mchap *= 1/180.;
  
  // On va d abord definir les matrices B chapeau (matrice de divergence) sur le Kchapeau
  Bchap1(0,0)=-1; 
  Bchap1(1,1)=1; 
  Bchap1(3,0)=1; Bchap1(3,1)=1; Bchap1(3,2)=2; 
  Bchap1(4,0)=-1; Bchap1(4,1)=-1; Bchap1(4,2)=-2;
  Bchap1(5,0)=1; Bchap1(5,1)=-1;
  Bchap1 *= 1/6.;

  Bchap2(0,0)=-1;
  Bchap2(2,2)=1;
  Bchap2(3,0)=1; Bchap2(3,1)=2; Bchap2(3,2)=1;
  Bchap2(4,0)=1; Bchap2(4,2)=-1; 
  Bchap2(5,0)=-1; Bchap2(5,1)=-2; Bchap2(5,2)=-1;
  Bchap2 *= 1/6.;
  
  // On definit ensuite les matrices de rigidite Wchapeau sur le Kchapeau
  Wchap11(0,0)=3; Wchap11(0,1)=1; Wchap11(0,5)=-4; 
  Wchap11(1,0)=1; Wchap11(1,1)=3; Wchap11(1,5)=-4;
  Wchap11(3,3)=8; Wchap11(3,4)=-8; 
  Wchap11(4,3)=-8; Wchap11(4,4)=8; 
  Wchap11(5,0)=-4; Wchap11(5,1)=-4; Wchap11(5,5)=8;
  Wchap11 *= 1/6.;
  
  Wchap22(0,0)=3; Wchap22(0,2)=1; Wchap22(0,4)=-4;
  Wchap22(2,0)=1; Wchap22(2,2)=3; Wchap22(2,4)=-4; 
  Wchap22(3,3)=8; Wchap22(3,5)=-8;
  Wchap22(4,0)=-4; Wchap22(4,2)=-4; Wchap22(4,4)=8;
  Wchap22(5,3)=-8; Wchap22(5,5)=8;
  Wchap22 *= 1/6.;

  Wchap12(0,0)=3; Wchap12(0,2)=1; Wchap12(0,4)=-4; 
  Wchap12(1,0)=1; Wchap12(1,2)=-1; Wchap12(1,3)=4; Wchap12(1,5)=-4;
  Wchap12(3,2)=4; Wchap12(3,3)=4; Wchap12(3,4)=-4; Wchap12(3,5)=-4;
  Wchap12(4,2)=-4; Wchap12(4,3)=-4; Wchap12(4,4)=4; Wchap12(4,5)=4;
  Wchap12(5,0)=-4; Wchap12(5,3)=-4; Wchap12(5,4)=4; Wchap12(5,5)=4;
  Wchap12 *= 1/6.;

}

void MatricesElementaires(const Triangle& K, 
			  KNM<R>& A, 
			  KNM<R>& B,
			  const KNM<R>& Bchap1,
			  const KNM<R>& Bchap2,
			  const KNM<R>& Wchap11,
			  KNM<R>& Wchap12,
			  const KNM<R>& Wchap22){  
  
  double aire = K.mesure();
  double valeur = 0.5*(1./aire);
  double x0 = K[0].x, y0 = K[0].y;
  double x1 = K[1].x, y1 = K[1].y;
  double x2 = K[2].x, y2 = K[2].y;
  KNM<R> W11(6,6), W22(6,6), B1(6,3), B2(6,3);
  for(int i=0;i<6;++i){
    for(int j=0;j<6;++j){
      W11(i,j) =  W22(i,j) = 0.;
      if(j<3) B1(i,j) = B2(i,j) = 0.;
    }
  }
  
  // Legere transformation de la formule pour eviter d'avoir a ecrire des soustractions de matrices (les | indiquent les endroits modifies)
  //                                           |    |      |
  //W11 = valeur*(pow(K[2].y-K[0].y,2)*Wchap11 + (K[0].y-K[1].y)*(K[2].y-K[0].y)*(Wchap12.t()+Wchap12) + pow(K[1].y-K[0].y,2)*Wchap22);
  // Decomposition de la ligne sinon probleme a la compilation : operations non reconnues.
  // Plus precisement les operations non reconnues sont : operator*(double, KNM_<R>) (le produit de la transposee par un reel)
  //                                                      affectation d'une somme de matrice   
  W11  = Wchap12.t();                    // transposee
  W11 += Wchap12;
  //W11 *= (K[0].y-K[1].y)*(K[2].y-K[0].y);
  W11 *= (y0-y1)*(y2-y0);
  //W11 += pow(K[2].y-K[0].y,2)*Wchap11;
  W11 += pow(y2-y0,2)*Wchap11;
  //W11 += pow(K[1].y-K[0].y,2)*Wchap22;
  W11 += pow(y1-y0,2)*Wchap22;
  W11 *= valeur;
  //W22 = valeur*(pow(K[2].x-K[0].x,2)*Wchap11+(K[0].x-K[1].x)*(K[2].x-K[0].x)*(Wchap12+Wchap12.t())+pow(K[1].x-K[0].x,2)*Wchap22);
  W22  = Wchap12.t();
  W22 += Wchap12;
  //W22 *= (K[0].x-K[1].x)*(K[2].x-K[0].x);
  W22 *= (x0-x1)*(x2-x0);
  //W22 += pow(K[2].x-K[0].x,2)*Wchap11;
  W22 += pow(x2-x0,2)*Wchap11;
  //W22 += pow(K[1].x-K[0].x,2)*Wchap22;
  W22 += pow(x1-x0,2)*Wchap22;
  W22 *= valeur;
  // Determinant de la matrice Bk (voir polycopie) Fk=Bkx+bk
  //double detBk = (K[1].x-K[0].x)*(K[2].y-K[0].y)-(K[2].x-K[0].x)*(K[1].y-K[0].y);
  double detBk = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
  // Le signe du determinant
  detBk = copysign(1,detBk);
  // Attention : (Bk^-1).transposee() = 1./det(Bk)*J !!
  // J(0,0) = y2-y0; J(0,1) = y0-y1;
  // J(1,0) = x0-x2; J(1,1) = x1-x0;
  // On definit les matrices de divergence B sur K
  //B1 = sgn(detBk)*( J(0,0)*Bchap1 + J(0,1)*Bchap2 );
  B1  = (y2-y0)*Bchap1;
  B1 += (y0-y1)*Bchap2;
  B1 *= detBk;
  //B2 = sgn(detBk)*( J(1,0)*Bchap1 + J(1,1)*Bchap2 );
  B2  = (x0-x2)*Bchap1;
  B2 += (x1-x0)*Bchap2;
  // Idem ci-dessus
  B2 *= detBk;

  // Initialisation A et B
  for(int i=0;i<12;++i){
    for(int j=0;j<12;++j){
      A(i,j) = 0.;
      if(j<3) B(i,j) = 0.;
    }
  }
  // Definition de la matrice A élémentaire (le petit a de la formulation variationelle)
  // Voir polycop pour la definition de A
  for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){
      A(i,j) = A(i+6,j+6) = W11(i,j) + W22(i,j);
    }
  }
  // Definition de la matrice B élémentaire (le petit b de la formulation variationelle)
  // B = (B1,B2)
  for(int i=0;i<=5;i++){
    for(int j=0;j<=2;j++){
      B(i,j) = B1(i,j);
      B(i+6,j) = B2(i,j);
    }
  }
}


double MelemQ(int i, int j)
{
  assert(i<6 && j<6);
  R phi[6];
  int sQP = QuadratureFormular_T_5.n;
  double res = 0;
  for(int iqp=0;iqp<sQP;++iqp){
    P2::init(QuadratureFormular_T_5(iqp),phi);
    res += QuadratureFormular_T_5(iqp).a*phi[j]*phi[i];
  };
  res *= 0.5;
  return res;
}

double Wchap11Q(int i, int j)
{
  assert(i<6 && j<6);
  R2 dphi[6];
  int sQP = QuadratureFormular_T_5.n;
  R2 QP[sQP];
  for(int ii=0;ii<sQP;++ii)
    QP[ii] = (QuadratureFormular_T_5(ii));
  double res = 0;
  for(int iqp=0;iqp<sQP;++iqp){
    P2::init(QP[iqp],dphi);
    res += QuadratureFormular_T_5(iqp).a*dphi[i].x*dphi[j].x;
  }
  res *= 0.5;
  return res;
}

double Wchap12Q(int i, int j)
{
  assert(i<6 && j<6);
  R2 dphi[6];
  int sQP = QuadratureFormular_T_5.n;
  R2 QP[sQP];
  for(int ii=0;ii<sQP;++ii)
    QP[ii] = (QuadratureFormular_T_5(ii));
  double res = 0;
  for(int iqp=0;iqp<sQP;++iqp){
    P2::init(QP[iqp],dphi);
    res += QuadratureFormular_T_5(iqp).a*dphi[i].x*dphi[j].y;
  }
  res *= 0.5;
  return res;
}

double Wchap22Q(int i, int j)
{
  assert(i<6 && j<6);
  R2 dphi[6];
  int sQP = QuadratureFormular_T_5.n;
  R2 QP[sQP];
  for(int ii=0;ii<sQP;++ii)
    QP[ii] = (QuadratureFormular_T_5(ii));
  double res = 0;
  for(int iqp=0;iqp<sQP;++iqp){
    P2::init(QP[iqp],dphi);
    res += QuadratureFormular_T_5(iqp).a*dphi[i].y*dphi[j].y;
  }
  res *= 0.5;
  return res;
}

double Bchap1Q(int i, int j)
{
  assert(i<6 && j<3);
  R2 dphi[6];
  R lambda[3];
  int sQP = QuadratureFormular_T_5.n;
  R2 QP[sQP];
  for(int ii=0;ii<sQP;++ii)
    QP[ii] = (QuadratureFormular_T_5(ii));
  double res = 0;
  for(int iqp=0;iqp<sQP;++iqp){
    P1::init(QP[iqp],lambda);
    P2::init(QP[iqp],dphi);
    res += QuadratureFormular_T_5(iqp).a*lambda[j]*dphi[i%6].x;
  }
  res *= 0.5;
  return res;
}

double Bchap2Q(int i, int j)
{
  assert(i<6 && j<3);
  R2 dphi[6];
  R lambda[3];
  int sQP = QuadratureFormular_T_5.n;
  R2 QP[sQP];
  for(int ii=0;ii<sQP;++ii)
    QP[ii] = (QuadratureFormular_T_5(ii));
  double res = 0;
  for(int iqp=0;iqp<sQP;++iqp){
    P1::init(QP[iqp],lambda);
    P2::init(QP[iqp],dphi);
    res += QuadratureFormular_T_5(iqp).a*lambda[j]*dphi[i].y;
  }
  res *= 0.5;
  return res;
}


// 
// MatriceElementaire.cpp ends here
