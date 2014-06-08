// poisson2d.cpp --- 
// 
// Filename: poisson2d.cpp
// Description: 
// Author: O.TISSOT 
// Created: sam. déc. 14 00:59:12 2013 (+0100)
// Version: 
// Last-Updated: mar. févr. 11 15:56:28 2014 (+0100)
//           By: Olivier
//     Update #: 291
// Compatibility: 
// 
// UMFPACK pour la resolution du systeme lineaire

// Commentary: 
// 
// Juste un petit test sur une equation simple de la
// matrice de rigidite A.
// On resoud 2 pbs de Poisson independants avec condition 
// de Dirichlet sur tous les bords sauf le numero 4.
// On a choisi de prendre nu = 1. mais changer sa valeur
// ne change pas l'allure de la solution.
//

// Code:

#include <ctime>                        // Timer
#include "MapMatrix.hpp"                // Matrices creuses
#include "Maillage.hpp"                 // Maillage
#include "MatriceElementaire.hpp"       // Matrices Elementaires associees
#include "Param.hpp"                    // Parametres du calculs
#include "umfpack.h"                    // UMFPACK

#define CHECK_MATRIX 0                  // Sauvegarde de la matrice dans un fichier

typedef KN<double> Vector;

// Condition au bord (avec surcharge pour faciliter l'ecriture)
double g(R2 p){
  return 16*(p.y -.5)*(1.-p.y);
}
double g(double y){
  return 16*(y -.5)*(1.-y);
}


int main(int argc, char ** argv){
  
  /*****************************************/
  // RECUPERATION DES PARAMETRES DU CALCUL //
  /*****************************************/
  
  Param paramstokes;
  string filename;
  clock_t tstart;
  tstart = clock();
  cout << "*********************************************" << endl << endl;
  cout << "Importation des parametres du calcul..." << endl << endl;
  if(argc==1)
    paramstokes.init("terminal");        // On devra rentrer les parametres dans le terminal
  else{
    filename = argv[1];
    if(argc>2){
      cout << "/!\\ Attention trop d'arguments donnes /!\\ " << endl;               // "\\" est une sequence d'echapement pour pouvoir afficher "\"
      cout << "Essai de construction des parametres..." << endl;
    }
    paramstokes.init(filename);
  }
  cout << "Parametres du calcul recuperes !" << endl << endl;

  cout << "*********************************************" << endl << endl;

  /**************************************/
  // CHARGEMENT DES DONNEES DU MAILLAGE //
  /**************************************/
 
  cout << "Chargement des donnes du maillage..." << endl;
  Maillage Th(paramstokes.Mesh);
  cout << "Donnees du maillage recuperees !" << endl;
  map< pair<int,int>, int > edge;
  cout << "Construction des arretes a l'interieur du domaine..." << endl;
  Th.buildEdge(edge);
  cout << "Arretes construites !" << endl << endl; 
  cout << "*********************************************" << endl << endl;

  /************************************/
  // ASSEMBLAGE DE LA MATRICE GLOBALE //
  /************************************/
  
  cout << "Assemblage de la matrice..." << endl;
  int nbedge = edge.size();
  int dimtot = 2*(Th.nv + nbedge);
  MapMatrix A(dimtot);
  KNM<R> Mchap(6,6), Wchap11(6,6), Wchap12(6,6), Wchap22(6,6), Bchap1(6,3), Bchap2(6,3);
  KNM<R> Aelem(12,12), Belem(12,3);
  int I, J;
  // Initialisation des matrices elementaires associees a Kchap
  MatricesElementairesKchap(Mchap,Bchap1,Bchap2,Wchap11,Wchap12,Wchap22);
  for(int k=0;k<Th.nt;++k){
    //cout << "triangle n°" << k << endl;
    MatricesElementaires(Th[k],Aelem,Belem,Bchap1,Bchap2,Wchap11,Wchap12,Wchap22);
    for(int i=0;i<12;++i){
      I = Th.numglob(k,edge,i);                          // Numero global associe
      for(int j=0;j<12;++j){
	if(Aelem(i,j)!=0){                               // Pour conserver une structure creuse
	  J = Th.numglob(k,edge,j);
	  A(I,J) += Aelem(i,j);
	}
      }
    }  
  }
  cout << "Matrice assemblee !" << endl << endl;
  cout << "*********************************************" << endl << endl;

  /*******************************/
  // ASSEMBLAGE DU SECOND MEMBRE //
  /*******************************/
  
  Vector Fh(dimtot);
  cout << "Assemblage du second membre..." << endl;
  map< pair<int,int>, int>::iterator it;
  const double tgv = 1e30;
  double midx, midy;
  int si, sim;
  // On initialise le second membre
  /*for(int i=0;i<dimtot;++i)
    Fh[i] = 0.;*/
  for(int e=0;e<Th.ne;++e){
    for(int i=0;i<1;++i){
      if(Th.borderEdge(e)[i].onGamma(6)){
	si = Th.borderEdge(e)[i];
	it = edge.find(make_pair(min((int)Th.borderEdge(e)[i],(int)Th.borderEdge(e)[(i+1)%2]),max((int)Th.borderEdge(e)[i],(int)Th.borderEdge(e)[(i+1)%2])));
	sim = it->second + Th.nv;
      	midx = (Th.borderEdge(e)[i].x+Th.borderEdge(e)[(i+1)%2].x)/2.;
	midy = (Th.borderEdge(e)[i].y+Th.borderEdge(e)[(i+1)%2].y)/2.;
	//cout << "(" << Th.borderEdge(e)[0].x << "," << Th.borderEdge(e)[0].y << ")" << endl;
	//cout << "(" << Th.borderEdge(e)[1].x << "," << Th.borderEdge(e)[1].y << ")" << endl;
	//cout << "(" << midx << "," << midy << ")" << endl;
	//cout << "(" << (Th(it->first.first).x+Th(it->first.second).x)/2. << "," << (Th(it->first.first).y+Th(it->first.second).y)/2. << ")" << endl;
	// U1
	Fh[si] = tgv*g(Th.borderEdge(e)[i]);
	if(midx == 0.) Fh[sim] = tgv*g(midy);
	// U2
	Fh[ si + nbedge + Th.nv ] = tgv;
	if(midx == 0.) Fh[ sim + Th.nv + nbedge ] = tgv;
      }
    }
  }
  // On fixe un point a pression egale a 1 (par exemple) :
  //Fh[ (int)(Th.borderEdge(0)[0]) + 2*(nbedge+Th.nv) ] = tgv;
  cout << "Second membre assemble !" << endl << endl;
  cout << "*********************************************" << endl << endl;
  
  /********************************/
  // AJOUT DES CONDITIONS DE BORD //
  /********************************/
  
  Vector u(dimtot);
  cout << "Ajout des conditions de bord..." << endl;
  // Penalisation pour modifier A sur les points du bord (methode TGV)
  for(int e=0;e<Th.ne;++e){
    for(int i=0;i<1;++i){
      // Les bords ou on met des conditions de Dirichlet :
      if(Th.borderEdge(e)[i].onGamma(1)||Th.borderEdge(e)[i].onGamma(2)||Th.borderEdge(e)[i].onGamma(3)||Th.borderEdge(e)[i].onGamma(5)||Th.borderEdge(e)[i].onGamma(6)){
      //if(!Th.borderEdge(e)[i].label){
	si = Th.borderEdge(e)[i];
	it = edge.find(make_pair(min((int)Th.borderEdge(e)[i],(int)Th.borderEdge(e)[(i+1)%2]),max((int)Th.borderEdge(e)[i],(int)Th.borderEdge(e)[(i+1)%2])));
	sim = it->second + Th.nv;
	/*
	  cout << "(" << Th.borderEdge(e)[0].x << "," << Th.borderEdge(e)[0].y << ")" << endl;
	  cout << "(" << (Th.borderEdge(e)[0].x+Th.borderEdge(e)[1].x)/2. << "," << (Th.borderEdge(e)[0].y+Th.borderEdge(e)[1].y)/2. << ")" << endl;
	*/
	// U1
	A(si,si) = tgv;
	A(sim,sim) = tgv;
	// U2
	A(si + nbedge + Th.nv,si + nbedge + Th.nv) = tgv;
	A(sim + nbedge + Th.nv,sim + nbedge + Th.nv ) = tgv;
      }
    }
  }
  cout << "Conditions de bord ajoutees !" << endl << endl;
  cout << "*********************************************" << endl << endl;

  /**********************************/
  // RESOLUTION DU SYSTEME LINEAIRE //
  /**********************************/
 
  cout << "Resolution du systeme lineaire..." << endl << endl;
  // Interfacage avec UMFPACK
  cout << "Preparation des donnes pour pouvoir utiliser UMFPACK..." << endl;
  int nnz = A.size();
  int * Ai = new int[A.getDimension()+1];
  int * Aj = new int[nnz];
  double * Ax = new double[nnz];
  A.toCSR(Ai,Aj,Ax);
  assert( nnz == Ai[dimtot] );
  cout << "Donnees preparees !" << endl;
  if(CHECK_MATRIX){
    cout << "Sauvegarde de la matrice dans un fichier pour check..." << endl;
    ofstream m("matriceglob.dat");
    m << dimtot << endl;
    for(int i=0;i<dimtot;++i){
      m << "ligne " << i << " : ";
      int kk=0;
      int j0=-1;
      for (int k=Ai[i];k<Ai[i+1];++k)
	{
	  kk++;
	  int j = Aj[k];
	  m << "(" << i << "," << j << ") : " << Ax[k] << " ; ";
	  // Check sur les colonnes
	  assert( j>=0 && j<dimtot);
	  assert(j0 <j);
	  j0 =j;
	}
      m << endl;
      if(kk==0) m << " ligne " << i << " vide"<< endl;
    }
    cout << "Matrice sauvegardee !" << endl;
  }
  // Resolution du systeme lineaire avec UMFPACK
  double *null = (double *) NULL ;
  void *Symbolic=0, *Numeric=0 ; // modif ..
  int status;
  cout << "Resolution du systeme lineaire avec UMFPACK..." << endl;
  status =  umfpack_di_symbolic(dimtot, dimtot, (const int *) Ai, Aj, Ax, &Symbolic, null, null);
  switch(status){
  case UMFPACK_ERROR_out_of_memory: cout << "/!\\ Error umpfack -> umfpack_di_symbolic, status : UMFPACK_ERROR_out_of_memory" << endl;break;
  case UMFPACK_ERROR_argument_missing: cout << "/!\\ Error umpfack -> umfpack_di_symbolic, status : UMFPACK_ERROR_argument_missing" << endl;break;
  case UMFPACK_ERROR_invalid_matrix: cout << "/!\\ Error umpfack -> umfpack_di_symbolic, status : UMFPACK_ERROR_invalid_matrix" << endl;break;
  case UMFPACK_ERROR_n_nonpositive: cout << "/!\\ Error umpfack -> umfpack_di_symbolic, status : UMFPACK_ERROR_n_nonpositive" << endl;break;
  case UMFPACK_ERROR_internal_error: cout << "/!\\ Error umpfack -> umfpack_di_symbolic, status : UMFPACK_ERROR_internal_error" << endl;break;  
  }
  assert(!status);
  status = umfpack_di_numeric((const int *) Ai, Aj, Ax, Symbolic, &Numeric, null, null);
  switch(status){
  case UMFPACK_WARNING_singular_matrix: cout << "/!\\ Error umpfack -> umfpack_di_numeric, status : UMFPACK_WARNING_singular_matrix" << endl;break;  
  case UMFPACK_ERROR_out_of_memory: cout << "/!\\ Error umpfack -> umfpack_di_numeric, status : UMFPACK_ERROR_out_of_memory" << endl;break;
  case UMFPACK_ERROR_argument_missing: cout << "/!\\ Error umpfack -> umfpack_di_numeric, status : UMFPACK_ERROR_argument_missing" << endl;break;
  case UMFPACK_ERROR_invalid_Symbolic_object: cout << "/!\\ Error umpfack -> umfpack_di_numeric, status : UMFPACK_ERROR_invalid_Symbolic_object" << endl;break;
  case UMFPACK_ERROR_different_pattern: cout << "/!\\ Error umpfack -> umfpack_di_numeric, status : UMFPACK_ERROR_different_pattern" << endl;break;
  }
  umfpack_di_free_symbolic(&Symbolic) ;
  assert(!status);
  status = umfpack_di_solve(UMFPACK_At,(const int *) Ai, Aj, Ax, u, Fh, Numeric, null, null); 
  switch(status){
  case UMFPACK_WARNING_singular_matrix: cout << "/!\\ Error umpfack -> umfpack_di_solve, status : UMFPACK_WARNING_singular_matrix" << endl;break;  
  case UMFPACK_ERROR_out_of_memory: cout << "/!\\ Error umpfack -> umfpack_di_solve, status : UMFPACK_ERROR_out_of_memory" << endl;break;
  case UMFPACK_ERROR_argument_missing: cout << "/!\\ Error umpfack -> umfpack_di_solve, status : UMFPACK_ERROR_argument_missing" << endl;break;
  case UMFPACK_ERROR_invalid_Numeric_object: cout << "/!\\ Error umpfack -> umfpack_di_solve, status : UMFPACK_ERROR_invalid_Numeric_object" << endl;break;
  case UMFPACK_ERROR_invalid_system: cout << "/!\\ Error umpfack -> umfpack_di_solve, status : UMFPACK_ERROR_invalid_system" << endl;break;
  }
  umfpack_di_free_numeric(&Numeric) ;
  assert(!status);
  cout << "Systeme lineaire resolu !" << endl << endl;
  
  cout << "Desallocation de la memoire..." << endl;
  if(Ai) delete[] Ai;
  if(Aj) delete[] Aj;
  if(Ax) delete[] Ax;
  cout << "Memoire liberee !" << endl << endl;
  
  cout << "*********************************************" << endl << endl;
  
  /*******************************************/
  // ECRITURE DE LA SOLUTION DANS UN FICHIER //
  /*******************************************/
  
  cout << "Sauvegarde de la solution..." << endl;
  Vector u1(Th.nv + nbedge), u2(Th.nv + nbedge);
  for(int k=0;k<Th.nv+nbedge;++k){
    u1[k] = u[k];
    u2[k] = u[k + Th.nv+nbedge];
  }
  cout << "U1max = " << u1.max() << " U1min = " << u1.min() << endl;
  cout << "U2max = " << u2.max() << " U2min = " << u2.min() << endl;
  {
    ofstream gp("results/plotU1.sol");
    gp << Th.nv+nbedge << endl;
    for(int k=0;k<Th.nv+nbedge;++k)
	  gp << u1[k] << endl;
  }
  {  
    ofstream gp("results/plotU2.sol");
    gp << Th.nv+nbedge << endl;
    for(int k=0;k<Th.nv+nbedge;++k)
	  gp << u2[k] << endl;
  }
  
  cout << "Solution sauvegardee !" << endl << endl;
  
  cout << "*********************************************" << endl << endl;
  cout << "Temps d'execution : " << (clock()-tstart)/(double)CLOCKS_PER_SEC << "s" << endl << endl;
  cout << "*********************************************" << endl << endl;
  
  string command = "./glplotisoP2 " + paramstokes.Mesh + " results/";
  string choix;
  cout << "Que souhaitez-vous visualiser? (U1 ou U2)" << endl;
  cin >> choix;
  if(choix == "U1")
    command += "plotU1.sol";
  else if(choix == "U2")
    command += "plotU2.sol";
  else{
    cout << "Choix non reconnu -> choix possibles : U1 ou U2 (et attention a la casse...)." << endl;
    exit(1);
  }
  assert(!system(command.c_str()));
}

// 
// poisson2d.cpp ends here
