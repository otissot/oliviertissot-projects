// stokes2d.cpp --- 
// 
// Filename: stokes2d.cpp
// Description: 
// Author: O.TISSOT 
// Created: lun. janv. 27 09:31:10 2014 (+0100)
// Version: 
// Last-Updated: sam. mai 31 13:26:26 2014 (+0200)
//           By: otissot
//     Update #: 1003
// Compatibility: 
// 
// UMFPACK

// Commentary: 
// 
// Resolution du probleme de Stokes en 2d
// 
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
    paramstokes.init("terminal");                                                   // On devra rentrer les parametres dans le terminal
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
  int dimtot = 3*Th.nv + 2*nbedge;
  const double eps = 1e-6;
  const double nu = paramstokes.Viscosite;
  const double dt = paramstokes.dt;
  const double cM = 1/12.;
  double cKM;
  MapMatrix A(dimtot);
  KNM<R> Melem(6,6), Wchap11(6,6), Wchap12(6,6), Wchap22(6,6), Bchap1(6,3), Bchap2(6,3);
  KNM<R> Belem(12,3), Aelem(12,12);
  int I, J;
  // Initialisation des matrices elementaires associees a Kchap
  MatricesElementairesKchap(Melem,Bchap1,Bchap2,Wchap11,Wchap12,Wchap22);
  /*
  for(int i=0;i<6;++i){
    for(int j=0;j<6;++j){
      Wchap11(i,j)  = Wchap11Q(i,j);
      Wchap12(i,j)  = Wchap12Q(i,j);
      Wchap22(i,j)  = Wchap22Q(i,j);
      Melem(i,j)    = MelemQ(i,j);
      Bchap1(i,j%3) = Bchap1Q(i,j%3);
      Bchap2(i,j%3) = Bchap2Q(i,j%3);
    }
  }
  */
  // Rappel : structure de la matrice elementaire (15x15)
  // Aelem       Belem
  // Belem.t()   -eps*MasseP1
  for(int k=0;k<Th.nt;++k){
    MatricesElementaires(Th[k],Aelem,Belem,Bchap1,Bchap2,Wchap11,Wchap12,Wchap22);
    // Terme de rigidite
    for(int i=0;i<12;++i){
      I = Th.numglob(k,edge,i);                          // Numero global associe : soit U1, soit U2
      for(int j=0;j<12;++j){
	if(Aelem(i,j)!=0){                               // Pour conserver une structure creuse
	  J = Th.numglob(k,edge,j);
	  A(I,J) += nu*Aelem(i,j);
	}
      }
    }
    // Terme de divergence
    for(int i=0;i<12;++i){
      I = Th.numglob(k,edge,i);                          // Numero global associe : soit U1, soit U2
      for(int j=12;j<15;++j){
	if(Belem(i,j-12)!=0){                            // Pour conserver une structure creuse
	  J = Th.numglob(k,edge,j);
	  A(I,J) -= Belem(i,j-12);
	  A(J,I) -= Belem(i,j-12);
	}
      }
    }
    // Un terme pour stabiliser le systeme et rendre la matrice inversible
    cKM = -eps*Th[k].mesure()*cM;
    for(int i=12;i<15;++i){
      I = Th.numglob(k,edge,i);
      for(int j=12;j<15;++j){
	J = Th.numglob(k,edge,j);
	A(I,J) += cKM*(1.+(i==j));
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
  double midy;
  int si, sim;
  for(int e=0;e<Th.ne;++e){
    if(Th.borderEdge(e)[0].onGamma(6)){
      si = Th.borderEdge(e)[0];
      it = edge.find(make_pair(min((int)Th.borderEdge(e)[0],(int)Th.borderEdge(e)[1]),max((int)Th.borderEdge(e)[0],(int)Th.borderEdge(e)[1])));
      sim = it->second + Th.nv;
      midy = (Th.borderEdge(e)[0].y+Th.borderEdge(e)[1].y)/2.;
      // U1
      Fh[si] = tgv*g(Th.borderEdge(e)[0]);
      Fh[sim] = tgv*g(midy);
      // U2
      Fh[ si + nbedge + Th.nv ] = 0.;
      Fh[ sim + Th.nv + nbedge ] = 0.;
    }
  }
  cout << "Second membre assemble !" << endl << endl;
  cout << "*********************************************" << endl << endl;
  
  /********************************/
  // AJOUT DES CONDITIONS DE BORD //
  /********************************/
  
  Vector u(dimtot);
  cout << "Ajout des conditions de bord..." << endl;
  // Penalisation pour modifier A sur les points du bord (methode TGV)
  for(int e=0;e<Th.ne;++e){
    // Les bords ou on met des conditions de Dirichlet :
    // Par defaut le point dans le coin en bas a droite a un label qui vaut 4, or on veut mettre du Dirichlet dessus donc on ajoute 
    // la condition supplementaire dans le if et on teste avant de mettre tgv sur les points milieux
    if(Th.borderEdge(e)[0].onGamma(1)||Th.borderEdge(e)[0].onGamma(2)||Th.borderEdge(e)[0].onGamma(3)||Th.borderEdge(e)[0].onGamma(5)||Th.borderEdge(e)[0].onGamma(6)||(Th.borderEdge(e)[0].onGamma(4)&&Th.borderEdge(e)[0].y==0.)){
      si = Th.borderEdge(e)[0];
      it = edge.find(make_pair(min((int)Th.borderEdge(e)[0],(int)Th.borderEdge(e)[1]),max((int)Th.borderEdge(e)[0],(int)Th.borderEdge(e)[1])));
      sim = it->second + Th.nv;
      // U1
      A(si,si) = tgv;
      if(!Th.borderEdge(e)[0].onGamma(4)) A(sim,sim) = tgv;
      // U2
      A(si + nbedge + Th.nv,si + nbedge + Th.nv) = tgv;
      if(!Th.borderEdge(e)[0].onGamma(4)) A(sim + nbedge + Th.nv,sim + nbedge + Th.nv) = tgv;
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
      if(kk==0) m << " ligne " << i << " vide" << endl;
    }
    cout << "Matrice sauvegardee !" << endl;
  }
  // Resolution du systeme lineaire avec UMFPACK
  double *null = (double *) NULL ;
  void *Symbolic=0, *Numeric=0 ;
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
  Vector p(Th.nv), u1(Th.nv + nbedge), u2(Th.nv + nbedge);
  for(int k=0;k<Th.nv;++k){
    p[k] = u[k+2*(Th.nv+nbedge)];
  }
  for(int k=0.;k<Th.nv+nbedge;++k){
    u1[k] = u[k];
    u2[k] = u[k + Th.nv+nbedge];
  }
  cout << "Pmax = " << p.max() << " Pmin = " << p.min() << endl;   
  cout << "U1max = " << u1.max() << " U1min = " << u1.min() << endl;
  cout << "U2max = " << u2.max() << " U2min = " << u2.min() << endl;
  {  
    ofstream gp("results/plotPression.sol");
    gp << Th.nv << endl;
    for(int k=0;k<Th.nv;++k)
      gp << p[k] << endl;
  }
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
  
  string command;
  string choix;
  
  cout << "Que souhaitez-vous visualiser? (U1, U2 ou P)" << endl;
  cin >> choix;
  if(choix == "P")
    command = "./glplotiso " + paramstokes.Mesh + " results/plotPression.sol";
  else if(choix == "U1")
    command = "./glplotisoP2 " + paramstokes.Mesh + " results/plotU1.sol";
  else if(choix == "U2")
    command = "./glplotisoP2 " + paramstokes.Mesh + " results/plotU2.sol";
  else{
    cout << "Choix non reconnu -> choix possibles : P, U1, U2 (et attention a la casse...)." << endl;
    exit(1);
  }
  assert(!system(command.c_str()));
  
  /*
  // Little postprocessing
  cout << "-- POST-PROCESSING --" << endl;
  string solPffname  = "results/plotPFreefem" + paramstokes.Mesh + ".sol";
  solPffname.replace(solPffname.find(".msh"),4,"");
  solPffname.replace(solPffname.find("mesh_data/"),10,"");
  fstream solPff(solPffname.c_str(), ifstream::in);
  Vector pff(Th.nv);
  int bufferff = 0;
  double l2resp = 0.0;
  solPff >> bufferff;
  for(int i=0;i<Th.nv;i++)
    solPff >> pff(i);
  solPff.close();
  double pffbary = 0.0, pbary = 0.0;
  R lambda[3];
  for(int it=0;it<Th.nt;it++){
    const R2 bary = Th[it](R2(1/3.,1/3.));
    P1::init(bary,Th.t+it,lambda);
    pbary   = p(Th[it][0])*lambda[0] + p(Th[it][1])*lambda[1] + p(Th[it][2])*lambda[2];
    pffbary = pff(Th[it][0])*lambda[0] + pff(Th[it][1])*lambda[1] + pff(Th[it][2])*lambda[2];
    l2resp  += Th[it].mesure()*(pffbary - pbary)*(pffbary - pbary);
  }
  l2resp = sqrt(l2resp);
  fstream resl2f("results/resl2Stokes01.dat", ostream::out | ostream::app);
  resl2f << paramstokes.Mesh << " " << l2resp << endl;
  cout << "Residu L^2 de P  (ff++ et P calculee ) : " << l2resp << endl << endl;
  
  cout << "*********************************************" << endl << endl;
  cout << "Fin du programme..." << endl;
  */
  /*
  string solU1ffname = "results/plotU1Freefem" + paramstokes.Mesh + ".sol";
  solU1ffname.replace(solU1ffname.find(".msh"),4,"");
  solU1ffname.replace(solU1ffname.find("mesh_data/"),10,"");  
  string solU2ffname = "results/plotU2Freefem" + paramstokes.Mesh + ".sol";
  solU2ffname.replace(solU2ffname.find(".msh"),4,"");
  solU2ffname.replace(solU2ffname.find("mesh_data/"),10,"");
  Vector u1ff(Th.nv+nbedge);
  Vector u2ff(Th.nv+nbedge);
  double l2resu1 = 0.0, l2resu2 = 0.0
  fstream solU1ff(solU1ffname.c_str(), ifstream::in);
  solU1ff >> bufferff;
  for(int i=0;i<Th.nv+nbedge;i++)
    solU1ff >> u1ff(i);
  solU1ff.close();
  fstream solU2ff(solU2ffname.c_str(), ifstream::in);
  for(int i=0;i<Th.nv+nbedge;i++)
    solU2ff >> u2ff(i);
  solU2ff.close();
  l2resu1 *= Th[0].mesure();
  l2resu2 *= Th[0].mesure();
  cout << "Residu L^2 de U1 (ff++ et U1 calculee) : " << l2resu1 << endl << endl;
  cout << "Residu L^2 de U2 (ff++ et U2 calculee) : " << l2resu2 << endl << endl;
  */  


}

// 
// stokes2d.cpp ends here
