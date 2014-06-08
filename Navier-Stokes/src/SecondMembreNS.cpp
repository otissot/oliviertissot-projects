// SecondMembreNS.cpp --- 
// 
// Filename: SecondMembreNS.cpp
// Description: 
// Author: O.TISSOT 
// Created: lun. mars  3 19:07:27 2014 (+0100)
// Version: 
// Last-Updated: ven. mai 30 07:46:11 2014 (+0200)
//           By: Olivier
//     Update #: 553
// Compatibility: 
// 
// 

// Commentary: 
// 
// Implemenation du .hpp
// 
// 

// Code:

#include "SecondMembreNS.hpp"
#include "QuadratureFormular.hpp"


// Condition au bord (avec surcharge pour faciliter l'ecriture)
double g(const R2& p){
  return 16*(p.y -.5)*(1.-p.y);
}
double g(double y){
  return 16*(y -.5)*(1.-y);
}

// Calcule le transporte du point
R2 Transport(const Maillage & Th,
	     int it, 
	     double *lambda,
             const  KN_<R> & U,
	     const  KN_<R> & V, 
	     const double & dt)
{
  const Triangle & T(Th[it]);
  const R2 Q[3]={ R2(T[0]),R2(T[1]),R2(T[2]) };
  R2 P  = lambda[0]*Q[0]  + lambda[1]*Q[1]  + lambda[2]*Q[2];
  int i0=T[0];
  int i1=T[1];
  int i2=T[2];
  R u = lambda[0]*U[i0] + lambda[1]*U[i1] + lambda[2]*U[i2];
  R v = lambda[0]*V[i0] + lambda[1]*V[i1] + lambda[2]*V[i2];
  return P - R2(u,v)*dt;
};

// Projection P2 : on utilise les 6 degres de libertes connus
R2 Projection(const Maillage & Th,
	      map< pair<int,int>, int > & edge,
	      const Triangle * T,
	      R * phi,
	      const KN_<R> & U,
	      const KN_<R> & V)
{
  R2 res;
  // Le numero du triangle (car pas de reallocation de la memoire dans Find() )
  int nt = T - Th.t;
  int i[6] = { (*T)[0],(*T)[1],(*T)[2],Th.numedge(nt,edge,0),Th.numedge(nt,edge,1),Th.numedge(nt,edge,2) };
  for(int ii=0;ii<6;++ii){
    res.x += phi[ii]*U[i[ii]];
    res.y += phi[ii]*V[i[ii]];
  }
  return res;
};

// Version Find() bugguee
void buildSecondMembreFind(const Maillage & Th,
		       map< pair<int,int>, int > & edge,
		       const KN_<R> & Fh,
		       const KN_<R> & u1,
		       const KN_<R> & u2,
		       const double & dt)
{
  int iedge=0;
  R2 PF;
  R2 Phat;
  R2 projres;
  bool outside;
  R phi[6];
  // Remarque : iK est un pointeur sur une constante MAIS pas un pointeur constant, c'est la difference entre const Triangle * et const * Triangle
  // -> il peut pointer sur different objets (mais qui doivent rester constants)
  const Triangle * iK=0;
  double lambda[6][3] = { {.1,.0,.0},{.0,.1,.0 },{.0,.0,.1},{0.,.5,.5},{.5,0.,.5},{.5,.5,0.} }; 
  for(int itria=0;itria<Th.nt;++itria){
    // Les degres de libertes P2
    for(int DoF=0;DoF<6;++DoF){
      //index = itria;
      PF = Transport(Th,itria,lambda[DoF],u1,u2,dt);
      // Peut-etre une meilleure facon de choisir tstart... Moyenne de plusieurs tirages aleatoires? (Quelle loi? Combien?)
      // Pour recuperer le numero du point
      //cout << "Transporte ----> " << PF << endl;
      iK = Th.Find(PF,Phat,outside,Th.t);
      //cout << "Apres le transport 1... -> " << PF << endl;
      // Si on est sur la frontiere
      if(outside && !iK){
	// On cherche le segment auquel appartient le point Phat
	// La recherche n'est pas optimisee MAIS il n'y a pas beaucoup de points sur la frontiere
	cout << "Projete sur la frontiere -> " << Phat << endl;
	while(!Th.be[iedge].contient(Phat))
	  ++iedge;
	double val = (PF-Th.be[iedge][0]).norme();
	// On projette la valeur a partir des valeurs aux sommets
	projres.x = val*u1[Th.be[iedge][0]] + (1-val)*u1[Th.be[iedge][1]];
	projres.x = val*u2[Th.be[iedge][0]] + (1-val)*u2[Th.be[iedge][1]];
	// Integration numerique
	projres *= .5*Th.be[iedge].mesure();
      }
      else if(outside && iK){
	/*cout << "Que faire?" << endl;
	cout << PF << endl;
	cout << (*iK)[0] << " ; " << (*iK)[1] << " ; " << (*iK)[2] << endl;
	cout << Phat << endl;*/
      }
      // Si on est PAS sur la frontiere
      else if (!outside) {
	P2::init(PF,iK,phi);
	projres = Projection(Th,edge,iK,phi,u1,u2);
	// Integration numerique
	projres *= .333333333333*iK->mesure();
	}
      // On met au bon endroit dans Fh
      Fh[Th.numglob(itria,edge,DoF)] = projres.x;
      Fh[Th.numglob(itria,edge,DoF+6)] = projres.y;
    }
  }
  // On desalloue la memoire
  //if(phi) delete[] phi;
  //if(iK) delete(iK);
};

int onBorder(R* l, const Triangle* T){
  if(l[0]>0 && l[1]>0 && l[2]>0) return 0;
  else{
    int is;
    for(is=0;is<2;++is){
      if(l[is] == 0)
	break;
      else
	continue;
    }
    if( (*T)[(is+1)%3].label!=0 && (*T)[(is+2)%3].label!=0 )
      return (*T)[(is+1)%3].label;
    else return 0;
  }
};

// Version Walk()
void buildSecondMembre(const Maillage & Th,
		       const KN_<R> & Fh,
		       const KN_<R> & U,
		       const KN_<R> & V,
		       const double dt,
		       int step_time,
		       map< pair<int,int>, int> & edge)
{
  R lambda[3];
  R phi[6];
  // QuadratureFormular_T_5 : formule de quadrature exacte a l'ordre 5
  int iit, sQP = QuadratureFormular_T_5.n;
  KN<R> UConvect(sQP), VConvect(sQP);
  double Ubuf, Vbuf;
  // Reinitialisation de Fh pour U1 et U2
  for(int i=0;i<2*(Th.nv+edge.size());++i)
    Fh[i] = .0;
  R2 QP[sQP];
  for(int it=0;it<Th.nt;++it){
    // Nos points de quadrature
    for(int i=0;i<sQP;++i)
      QP[i] = Th[it](QuadratureFormular_T_5(i));
    // On calcule u_n°X_n pour chacun de ces points et ensuite on utilise la formule de quadrature
    for(int iqp=0;iqp<sQP;++iqp){
      // On calcule les coordonnees barycentriques du pt a deplacer 
      P1::init(QP[iqp],Th.t+it,lambda);
      iit = it;
      // -dt car on remonte le temps
      Walk(Th,iit,lambda,U,V,-dt);

      // On projete la valeur grace aux coordonnees barycentriques
      P2::init(lambda,phi);
      UConvect[iqp] = U[Th.numglob(iit,edge,0)]*phi[0];
      VConvect[iqp] = V[Th.numglob(iit,edge,0)]*phi[0];
      for(int iDoF=1;iDoF<6;++iDoF){
	// Calcul de "u_n°X_n.x(QP[.])" et "u_n°X_n.y(QP[.])"
	UConvect[iqp] += U[Th.numglob(iit,edge,iDoF)]*phi[iDoF];
	VConvect[iqp] += V[Th.numglob(iit,edge,iDoF)]*phi[iDoF];
      }
    }
    // On utilise une formule de quadrature pour calculer int_K( phi[iDoF]*u_n°X_n )
    for(int DoF=0;DoF<6;++DoF){
      Ubuf = Vbuf = 0.;
      for(int iqp=0;iqp<sQP;++iqp){
	P2::init(QP[iqp],Th.t+it,phi);
	Ubuf += QuadratureFormular_T_5(iqp).a*UConvect[iqp]*phi[DoF];
	Vbuf += QuadratureFormular_T_5(iqp).a*VConvect[iqp]*phi[DoF];
      }
      double tau = Th[it].mesure()/dt;
      Ubuf *= tau;
      Vbuf *= tau;
      Fh[Th.numglob(it,edge,DoF)]   += Ubuf;
      Fh[Th.numglob(it,edge,DoF+6)] += Vbuf;
    }
  }
};

// 
// SecondMembreNS.cpp ends here
