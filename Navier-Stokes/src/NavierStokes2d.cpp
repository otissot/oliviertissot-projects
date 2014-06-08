// NavierStokes.cpp --- 
// 
// Filename: stokes2d.cpp
// Description: 
// Author: O.TISSOT 
// Created: jeu. févr. 13 13:20:10 2014 (+0100)
// Version: 
// Last-Updated: sam. mai 31 15:39:29 2014 (+0200)
//           By: otissot
//     Update #: 828
// Compatibility: 
// 
// UMFPACK

// Commentary: 
// 
// Resolution du probleme de Navier Stokes en 2d
// avec la methode des caracteristiques et des 
// Elements Finis P1/P2 (pression/vitesse)

// Code:

#include <ctime>                        // Timer
#include "MapMatrix.hpp"                // Matrices creuses
#include "Maillage.hpp"                 // Maillage
#include "MatriceElementaire.hpp"       // Matrices Elementaires associees
#include "SecondMembreNS.hpp"           // Le second membre de Navier-Stokes
#include "Param.hpp"                    // Parametres du calculs
#include "umfpack.h"                    // UMFPACK

typedef KN<double> Vector;


int main(int argc, char ** argv){
  
  /*****************************************/
  // RECUPERATION DES PARAMETRES DU CALCUL //
  /*****************************************/
  
  Param paramNS;
  string paramfilename;
  clock_t tstart;
  tstart = clock();
  cout << "*********************************************" << endl << endl;
  cout << "Importation des parametres du calcul..." << endl << endl;
  if(argc==1)
    paramNS.init("terminal");                                                   // On devra rentrer les parametres dans le terminal
  else{
    paramfilename = argv[1];
    if(argc>2){
      cout << "/!\\ Attention trop d'arguments donnes /!\\ " << endl;               // "\\" est une sequence d'echapement pour pouvoir afficher "\"
      cout << "Essai de construction des parametres..." << endl;
    }
    paramNS.init(paramfilename);
  }
  cout << "Parametres du calcul recuperes !" << endl << endl;

  cout << "*********************************************" << endl << endl;

  /**************************************/
  // CHARGEMENT DES DONNEES DU MAILLAGE //
  /**************************************/
 
  cout << "Chargement des donnes du maillage..." << endl;
  Maillage Th(paramNS.Mesh);
  cout << "Donnees du maillage recuperees !" << endl;
  cout << "Construction du tableau des triangles adjacents..." << endl;
  Th.BuildAdj();
  cout << "Triangles adjacents connus !" << endl;
  map< pair<int,int>, int > edge;
  cout << "Construction des arretes a l'interieur du domaine..." << endl;
  Th.buildEdge(edge);
  cout << "Arretes construites !" << endl << endl; 
  cout << "*********************************************" << endl << endl;
  
  /********************************/
  // INITIALISATION DES VARIABLES //
  /********************************/
  
  const int nbedge = edge.size();
  const int dimtot = 3*Th.nv + 2*nbedge;
  const double Tmax = paramNS.Tmax;
  const double TSnapShot = paramNS.TSnapShot;
  const double eps = 1e-9;
  const double nu = paramNS.Viscosite;
  const double cM = 1/12.;
  const double dt = paramNS.dt;
  const double tau = 1./dt;
  const double tgv = 1e30;
  map< pair<int,int>, int>::iterator it;
  double midy;
  int I, J, si, sim;
  string filename, timename;
  char * buf = (char *) malloc(sizeof(*buf) * 10);
  double cKM;
  Vector u(dimtot), Fh(dimtot), p(Th.nv), u1(Th.nv + nbedge), u2(Th.nv + nbedge);
  MapMatrix A(dimtot), ANS(dimtot);
  KNM<R> Melem(6,6), Wchap11(6,6), Wchap12(6,6), Wchap22(6,6), Bchap1(6,3), Bchap2(6,3);
  KNM<R> Belem(12,3), Aelem(12,12);
  int nnz;
  int * Ai;
  int * Aj;
  double * Ax;
  
  /***********************************/
  // INITIALISATION DU SECOND MEMBRE //
  /***********************************/
  
  cout << "Initialisation second membre..." << endl;
  for(int e=0;e<Th.ne;++e){
    if(Th.be[e][0].onGamma(6)){
      si = Th.be[e][0];
      it = edge.find(make_pair(min((int)Th.be[e][0],(int)Th.be[e][1]),max((int)Th.be[e][0],(int)Th.be[e][1])));
      sim = it->second + Th.nv;
      midy = (Th.be[e][0].y+Th.be[e][1].y)/2.;
      // U1
      Fh[si]  = tgv*g(Th.be[e][0]);
      Fh[sim] = tgv*g(midy);
      // U2
      Fh[ si + nbedge + Th.nv ] = 0.;
      Fh[ sim + Th.nv + nbedge ] = 0.;
    }
  }
  cout << "Second membre initialise !" << endl << endl;
  cout << "*********************************************" << endl << endl;

  /****************************/
  // ASSEMBLAGE NAVIER-STOKES //
  /****************************/
  
  cout << "Assemblage de la matrice de Navier-Stokes..." << endl;
  // Initialisation des matrices elementaires associees a Kchap
  MatricesElementairesKchap(Melem,Bchap1,Bchap2,Wchap11,Wchap12,Wchap22);
  
  // Rappel : structure de la matrice elementaire (15x15)
  // Aelem       Belem
  // Belem.t()   -eps*MasseP1

  for(int k=0;k<Th.nt;++k){
    MatricesElementaires(Th[k],Aelem,Belem,Bchap1,Bchap2,Wchap11,Wchap12,Wchap22);
    // Terme de rigidite
    for(int i=0;i<12;++i){
      I = Th.numglob(k,edge,i);                                           // Numero global associe : soit U1, soit U2
      for(int j=0;j<12;++j){
	if(Aelem(i,j)!=0){                                                // Pour conserver une structure creuse
	  J = Th.numglob(k,edge,j);
	  A(I,J) += nu*Aelem(i,j);
	}
      }
    }
    // Terme de masse P2
    // U1
    for(int i=0;i<6;++i){
      I = Th.numglob(k,edge,i);                                           // Numero global associe : soit U1, soit U2
      for(int j=0;j<6;++j){
	if(Melem(i,j)!=0){                                                // Pour conserver une structure creuse
	  J = Th.numglob(k,edge,j);
	  A(I,J) += tau*Melem(i,j)*Th[k].mesure();
	}
      }
    }
    // U2
    for(int i=0;i<6;++i){
      I = Th.numglob(k,edge,i+6);                                         // Numero global associe : soit U1, soit U2
      for(int j=0;j<6;++j){
	if(Melem(i,j)!=0){                                                // Pour conserver une structure creuse
	  J = Th.numglob(k,edge,j+6);
	  A(I,J) += tau*Melem(i,j)*Th[k].mesure();
	}
      }
    }
    // Terme de divergence
    for(int i=0;i<12;++i){
      I = Th.numglob(k,edge,i);                          // Numero global associe : soit U1, soit U2
      for(int j=12;j<15;++j){
	if(Belem(i,j-12)!=0){                            // Pour conserver une structure creuse
	  J = Th.numglob(k,edge,j);
	  A(I,J) = A(J,I) -= Belem(i,j-12);
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
  cout << "Matrice Assemblee !" << endl;
  cout << "Ajout des conditions de bord..." << endl;
    // Penalisation pour modifier A sur les points du bord (methode TGV)
  for(int e=0;e<Th.ne;++e){
    // Les bords ou on met des conditions de Dirichlet :
    // Par defaut le point dans le coin en bas a droite a un label qui vaut 4, or on veut mettre du Dirichlet dessus donc on ajoute 
    // la condition supplementaire dans le if et on teste avant de mettre tgv sur les points milieux
    if(Th.be[e][0].onGamma(1)||Th.be[e][0].onGamma(2)||Th.be[e][0].onGamma(3)||Th.be[e][0].onGamma(5)||Th.be[e][0].onGamma(6)||(Th.be[e][0].onGamma(4)&&Th.be[e][0].y==0.)){
      si = Th.be[e][0];
      it = edge.find(make_pair(min((int)Th.be[e][0],(int)Th.be[e][1]),max((int)Th.be[e][0],(int)Th.be[e][1])));
      sim = it->second + Th.nv;
      // U1
      A(si,si) = tgv;
      if(!Th.be[e][0].onGamma(4))
	A(sim,sim) = tgv;
      // U2
      A(si + nbedge + Th.nv,si + nbedge + Th.nv) = tgv;
      if(!Th.be[e][0].onGamma(4))
	A(sim + nbedge + Th.nv,sim + nbedge + Th.nv) = tgv;
    }
  }
  cout << "Conditions de bords ajoutees !" << endl << endl;
  
  cout << "*********************************************" << endl << endl;
  
  // Interfacage UMFPACK
  nnz = A.size();
  Ai = new int[A.getDimension()+1];
  Aj = new int[nnz];
  Ax = new double[nnz];
  A.toCSR(Ai,Aj,Ax);
  assert( nnz == Ai[dimtot] );

  /*********************************/
  // DEBUT DES ITERATIONS EN TEMPS //
  /*********************************/

  cout << "Debut du calcul..." << endl << endl;
  
  for(int l=1;l<Tmax/dt;++l){  


    for(int k=0;k<Th.nv+nbedge;++k){
      u1[k] = u[k];
      u2[k] = u[k + Th.nv + nbedge];
    }
    for(int k=0;k<Th.nv;++k){
      p[k] = u[k+2*(Th.nv+nbedge)];
    }
    
    /********************************/
    // MISE A JOUR DU SECOND MEMBRE //
    /********************************/
    
    // Nouveau second Membre
    buildSecondMembre(Th,Fh,u1,u2,dt,l,edge);

    /********************************/
    // AJOUT DES CONDITIONS DE BORD //
    /********************************/
    
    // On remet les conditons de bord pour simuler l'entree du fluide
    
    for(int e=0;e<Th.ne;++e){
      if(Th.be[e][0].onGamma(6)){
	si = Th.be[e][0];
	it = edge.find(make_pair(min((int)Th.be[e][0],(int)Th.be[e][1]),max((int)Th.be[e][0],(int)Th.be[e][1])));
	sim = it->second + Th.nv;
	midy = (Th.be[e][0].y+Th.be[e][1].y)/2.;
	// U1
	Fh[si]  = g(Th.be[e][0]);
	Fh[sim] = g(midy);
	// U2
	Fh[ si + nbedge + Th.nv ]  = 0.;
	Fh[ sim + Th.nv + nbedge ] = 0.;
      }
    }
    
    // Penalisation pour modifier Fh sur les points du bord (methode TGV)
    for(int e=0;e<Th.ne;++e){
      // Les bords ou on met des conditions de Dirichlet :
      // Par defaut le point dans le coin en bas a droite a un label qui vaut 4, or on veut mettre du Dirichlet dessus donc on ajoute 
      // la condition supplementaire dans le if et on teste avant de mettre tgv sur les points milieux
      if(Th.be[e][0].onGamma(1)||Th.be[e][0].onGamma(2)||Th.be[e][0].onGamma(3)||Th.be[e][0].onGamma(5)||Th.be[e][0].onGamma(6)||(Th.be[e][0].onGamma(4)&&Th.be[e][0].y==0.)){
	si = Th.be[e][0];
	it = edge.find(make_pair(min((int)Th.be[e][0],(int)Th.be[e][1]),max((int)Th.be[e][0],(int)Th.be[e][1])));
	sim = it->second + Th.nv;
	// U1
	Fh(si) *= tgv;
	if(!Th.be[e][0].onGamma(4))
	  Fh(sim) *= tgv;
	// U2
	Fh(si + nbedge + Th.nv) *= tgv;
	if(!Th.be[e][0].onGamma(4))
	  Fh(sim + nbedge + Th.nv) *= tgv;
      }
    }
    
    /**********************************/
    // RESOLUTION DU SYSTEME LINEAIRE //
    /**********************************/
    
    double *null = (double *) NULL ;
    void *Symbolic=0, *Numeric=0 ;
    int status;
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
    
    
    /*******************************************/
    // ECRITURE DE LA SOLUTION DANS UN FICHIER //
    /*******************************************/
      for(int k=0;k<Th.nv+nbedge;++k){
      u1[k] = u[k];
      u2[k] = u[k + Th.nv+nbedge];
    }
    
    //if( l%(int)(TSnapShot/dt) == 0){
      for(int k=0;k<Th.nv;++k){
	p[k] = u[k+2*(Th.nv+nbedge)];
      }
      
      cout << " Iteration " << l << " -> u1 : " << u1.min() << " " << u1.max() << endl;
      cout << "                 u2 : " << u2.min() << " " << u2.max() << endl;
      cout << "                 p  : " << p.min() << " " << p.max() << endl << endl;
      
      sprintf(buf,"%d",l);
      timename = buf;
      { 
	filename = "results/P/plotPression"+timename+".sol";
	ofstream gp(filename.c_str());
	gp << Th.nv << endl;
	for(int k=0;k<Th.nv;++k)
	  gp << p[k] << endl;
      }
      {  
	filename = "results/U1/plotU1"+timename+".sol";
	ofstream gp(filename.c_str());
	gp << Th.nv+nbedge << endl;
	for(int k=0;k<Th.nv+nbedge;++k)
	  gp << u1[k] << endl;
      }
      {  
	filename = "results/U2/plotU2"+timename+".sol";
	ofstream gp(filename.c_str());
	gp << Th.nv+nbedge << endl;
	for(int k=0;k<Th.nv+nbedge;++k)
	  gp << u2[k] << endl;
      }
      //cout << "Solution sauvegardee !" << endl;
      // }
    
  }  
  // Fin de la boucle en temps
  
  cout << endl << "*********************************************" << endl << endl;
  cout << "Temps d'execution : " << (clock()-tstart)/(double)CLOCKS_PER_SEC << "s" << endl << endl;
  cout << "*********************************************" << endl << endl;
  cout << "Desallocation de la memoire..." << endl;
  if(Ai) delete[] Ai;
  if(Aj) delete[] Aj;
  if(Ax) delete[] Ax;
  cout << "Memoire liberee !" << endl << endl;

  string command;
  string choix;
  string stp_time;
  int step_time;
  
  cout << "Que souhaitez-vous visualiser? (U1, U2 ou P)" << endl;
  cin >> choix;
  cout << "Quel pas de temps? (entier <" << Tmax/dt << ")" << endl; 
  cin >> step_time;
  char * buft = (char *) malloc(sizeof(*buft) * 10);
  sprintf(buft,"%d",step_time);
  stp_time = buft;
  if(choix == "P")
    command = "./glplotiso " + paramNS.Mesh + " results/P/plotPression"+stp_time+".sol";
  else if(choix == "U1")
    command = "./glplotisoP2 " + paramNS.Mesh + " results/U1/plotU1"+stp_time+".sol";
  else if(choix == "U2")
    command = "./glplotisoP2 " + paramNS.Mesh + " results/U2/plotU2"+stp_time+".sol";
  else{
    cout << "Choix non reconnu -> choix possibles : P, U1, U2 (et attention a la casse...)." << endl;
    exit(1);
  }
  assert(!system(command.c_str()));
  
  /*
  cout << "-- POST-PROCESSING --" << endl;
  string solPffname  = "results/P/plotPFreefem" + paramNS.Mesh + "t99.sol";
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
  
  string solPffname1  = "results/P/plotPFreefem" + paramNS.Mesh + "t50.sol";
  solPffname1.replace(solPffname1.find(".msh"),4,"");
  solPffname1.replace(solPffname1.find("mesh_data/"),10,"");
  fstream solPff1(solPffname1.c_str(), ifstream::in);
  Vector pff1(Th.nv);
  double l2resp1 = 0.0;
  solPff1 >> bufferff;
  for(int i=0;i<Th.nv;i++)
    solPff1 >> pff1(i);
  solPff1.close();
  string solPname1  = "results/P/plotPression50.sol";
  fstream solP1(solPname1.c_str(), ifstream::in);
  Vector p1(Th.nv);
  solP1 >> bufferff;
  for(int i=0;i<Th.nv;i++)
    solP1 >> p1(i);

  for(int it=0;it<Th.nt;it++){
    const R2 bary = Th[it](R2(1/3.,1/3.));
    P1::init(bary,Th.t+it,lambda);
    pbary   = p1(Th[it][0])*lambda[0] + p1(Th[it][1])*lambda[1] + p1(Th[it][2])*lambda[2];
    pffbary = pff1(Th[it][0])*lambda[0] + pff1(Th[it][1])*lambda[1] + pff1(Th[it][2])*lambda[2];
    l2resp1  += Th[it].mesure()*(pffbary - pbary)*(pffbary - pbary);
  }

  cout << "Residu L^2 de P(Tf)   (ff++ et P calculee ) : " << l2resp << endl << endl;
  cout << "Residu L^2 de P(Tf/2)  (ff++ et P calculee ) : " << l2resp1 << endl << endl;

  fstream resl2f("results/resl2NS0025.dat", ostream::out | ostream::app);
  resl2f << paramNS.Mesh << " " << l2resp1 << " " << l2resp << endl;

  cout << "*********************************************" << endl << endl;
  cout << "Fin du programme..." << endl;
  */
}

// 
// stokes2d.cpp ends here
