#include <GL/glut.h>
#include <GL/gl.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "maillage.hpp"
#include "mat_cr.hpp"
#include "assemblage.hpp"		//routines d'assemblage

using namespace std;

//====================//
// Variables globales //
//====================//

int Height,Width;
double Rho;
double Phi;
double Theta;

double X0,Y0, Z0;
double Xcam,Ycam, Zcam;

int xold;
int yold;

vect u;
vector< triangle > t;
vector< point > p;
vector< edge > e;

//=======================//
// Fonction pour gerer   //
// la carte de couleur   //
//=======================//

void hsvToRgb (float h, float s, float v, float & r, float & g, float & b)
{
  int i;
  float aa, bb, cc, f;
  
  if (s == 0) /* Grayscale */
    r = g = b = v;
  else {
    if (h == 1.0) h = 0;
    h *= 6.0;
    i =  int(h);
    f = h - i;
    aa = v * (1 - s);
    bb = v * (1 - (s * f));
    cc = v * (1 - (s * (1 - f)));
    switch (i) {
    case 0: r = v;  g = cc; b = aa; break;
    case 1: r = bb; g = v;  b = aa; break;
    case 2: r = aa; g = v;  b = cc; break;
    case 3: r = aa; g = bb; b = v;  break;
    case 4: r = cc; g = aa; b = v;  break;
    case 5: r = v;  g = aa; b = bb; break;
    }
  }
}


void SetColor(double f)
{
  float r=0,g=0,b=0;
  double fmax=1; // bornes de la fonction
  double fmin=-1;
  
  hsvToRgb(.66*(fmax-f)/(fmax-fmin),1,1,r,g,b);		//.66 car rouge correspond a 0% et bleu a 66% ensuite on revient vers le rouge par le violet et de 0% a 66% on passe par le vert
  glColor3f(r,g,b);
}

void Clean() {
  glClearColor(1.0, 1.0, 1.0, 0.0);  // couleur du fond =  blanc en RGB
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}


// Cette procedure recupere, lorsque l'on tape une touche du clavier,
// la valeur de la touche tapee (key) et la position de la souris
// (x,y), donnee en pixels.
static void Key( unsigned char key, int x, int y ) { 

  switch (key) {

  case 27: 
    // Touche "escape" pour quitter
    exit(0);
    break;
  case 'a':
    Rho = Rho/1.1;
    break;
  case 'e':
    Rho = Rho*1.1;
    break;
  case 'z':
    Y0   = Y0-0.1;
    break;
  case 's':
    Y0   = Y0+0.1;
    break;
  case 'd':
    X0   = X0-0.1;
    break;
  case 'q':
    X0   = X0+0.1;
    break;
    
  }
  // redessine la fenetre
  glutPostRedisplay();
}



static void Mouse( int button,int state,int x,int y )
{
  // state up or down 
  if (state == GLUT_DOWN) { xold=x; yold=y;return;}
}



static void MotionMouse(int x,int y )
{
  Theta -= (y-yold)/(2.*180.);
  Phi   -= (x-xold)/(2*180.);
  xold = x;
  yold = y;
  glutPostRedisplay();
}


void SetView() {
  
  // Orientation de la  camera
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  Xcam  = X0 + Rho*sin(Theta)*cos(Phi);
  Ycam  = Y0 + Rho*sin(Theta)*sin(Phi);
  Zcam =  Z0 + Rho*cos(Theta);
  gluLookAt(Xcam,Ycam,Zcam,X0,Y0,Z0,0.,0.,1.);  

  // On choisit le mode de visualisation
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 
  gluPerspective(60., Width/Height,0,1000);

}

// fonction que l'on veut 
// representer sur le maillage
double s(double X, double Y, int k){
  vector< double > v(2,1/sqrt(2));
  return sin(k*X*v[0]+k*Y*v[1]);

}

// Cette procedure permet de dessiner dans la fenetre.
void Draw() {

  // Choix de couleur
  glColor3f(0., 0., 1.);

//======================//  
//  TRACE DU MAILLAGE   //
//======================// 
  
  int I; double X, Y;
  for(unsigned j=0; j<t.size(); j++){
    //glBegin(GL_LINES);
    glBegin(GL_TRIANGLES);

    // 1er sommet
    I = t[j].tri[0];
    X = p[I].coord[0], Y = p[I].coord[1];
    //SetColor(s(X,Y,10));
    //glColor3f(0.,0.,0.);
    SetColor(u[I]);
    glVertex3f(X,Y,0.);       

    // 2eme sommet
    I = t[j].tri[1];
    X = p[I].coord[0], Y = p[I].coord[1];
    //SetColor(s(X,Y,10));
    //glColor3f(0.,0.,0.);
    SetColor(u[I]);
    glVertex3f(X,Y,0.);

    // 3eme sommet
    I = t[j].tri[2];
    X = p[I].coord[0], Y = p[I].coord[1];
    //glColor3f(0.,0.,0.);
    //SetColor(s(X,Y,10));
    SetColor(u[I]);
    glVertex3f(X,Y,0.);    
    
    glEnd();  
  }
 
    //Affichage de l'echelle de couleur :

    glBegin(GL_QUADS);
    SetColor(-1.);
    glVertex3f(1.5,0.2,0.0);
    SetColor(-1.);
    glVertex3f(1.55,0.2,0.0);
    SetColor(-.33);
    glVertex3f(1.55,0.4,0.0);
    SetColor(-.33);
    glVertex3f(1.5,0.4,0.0);
    glEnd();
    glBegin(GL_QUADS);
    SetColor(-.33);
    glVertex3f(1.5,0.4,0.0);
    SetColor(-.33);
    glVertex3f(1.55,0.4,0.0);
    SetColor(.33);
    glVertex3f(1.55,0.6,0.0);
    SetColor(.33);
    glVertex3f(1.5,0.6,0.0);
    glEnd();
    glBegin(GL_QUADS);
    SetColor(.33);
    glVertex3f(1.5,0.6,0.0);
    SetColor(.33);
    glVertex3f(1.55,0.6,0.0);
    SetColor(1.);
    glVertex3f(1.55,0.8,0.0);
    SetColor(1.);
    glVertex3f(1.50,0.8,0.0);
    glEnd();

    //affichage des unites sur le cote :
    //ce morceau est assez long car de base la bibliotheque Glut n'a pas la fonction glutBitmapString implementee...
    SetColor(2.0);		//pour afficher en noir
    char * buf = (char *) malloc( sizeof(*buf) * 5);
    for(int i=1;i<4;i++){
    	glRasterPos2f(1.55,.2+i*.8/10.);
    	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'_');
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,' ');
	sprintf(buf,"%f",i/4.-1);
	for(int k=0;k<5;k++)
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,buf[k]);
    }
    for(int i=1;i<4;i++){
    	glRasterPos2f(1.55,.5+i*.8/10.);
    	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'_');
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,' ');
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'+');
	sprintf(buf,"%f",i/4.);
	for(int k=0;k<4;k++)
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,buf[k]);
    }
    glRasterPos2f(1.55,0.8);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'_');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,' ');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'+');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'1');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'.');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'0');
    glRasterPos2f(1.55,0.2);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'_');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,' ');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'-');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'1');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'.');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'0');
    glRasterPos2f(1.55,0.5);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'_');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,' ');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'0');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'.');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'0');

    free(buf);
}


static void Reshape( int width, int height ) {
  Width  = width;
  Height = height;
  SetView();
  glutPostRedisplay();  //  pour redessiner
}


void Display(void) {
  Clean();              // nettoie l'ecran
  SetView();
  Draw();
  glFlush();            // vide les tampons GL
  glutSwapBuffers();    // affiche le buffer du dessin
}

int main(int argc, char* argv[]){


//=====================//  
//  Initialisation des //
//  variables globales //
//=====================//
  
  X0 = 3/4.; 
  Y0 = 1/2.; 
  Z0 = 0;
  Width  = 500; 
  Rho    = 2.5; 
  Height = 500; 
  Phi   = M_PI/2.;
  Theta = -0.000001;		//si on met 0 on ne voit rien
  
  clock_t start;
  clock_t endtime;

//=============================//  
//  IMPORTATION DU MAILLAGE    //
//=============================//
 
  start = clock();
  vector< edge > e;
  cout<<"************************************************"<<endl;
  cout<<"Importation du maillage ..."<<endl;
  if(argc==1){
	  cout<<"Importation ratee : fichier non specifie !"<<endl;
	  exit(1);
  }	
  maillage m(argv[1]);
  p = m.get_p();
  e = m.get_e();
  t = m.get_t();
  endtime= clock();
  cout<<"Temps d'importation : "<<(endtime-start)/(double) CLOCKS_PER_SEC<<endl;
  cout<<"************************************************"<<endl<<endl;

//=========================//
//    ASSEMBLAGE MATRICE   //
//=========================//
  
  matrice a;
  matrice A;
  cout<<"************************************************"<<endl;
  cout<<"Debut assemblage matrice"<<endl;
  start = clock();
  build_mat(a, A, p, e, t, argc, argv);
  endtime = clock();
  cout<<"Fin assemblage matrice"<<endl;
  cout<<"Temps d'assemblage : "<<(endtime-start)/(double) CLOCKS_PER_SEC<<endl;
  cout<<"************************************************"<<endl<<endl;

//==============================//
//  CONSTRUCTION SECOND MEMBRE  //
//==============================//

  vect f, init;
  cout<<"************************************************"<<endl;
  cout<<"Construction second membre ..."<<endl;
  start = clock();
  build_f(f, e, p, atof(argv[3]));
  endtime= clock();
  cout<<"Temps de construction : "<<(endtime-start)/(double) CLOCKS_PER_SEC<<endl;
  cout<<"************************************************"<<endl<<endl;

//==============================//  
//  RESOLUTION SYSTEME LINAIRE  //
//==============================//

  init.resize(a.size(),0);
  cout<<"************************************************"<<endl;
  cout<<"Resolution systeme de taille "<<p.size()<<" ..."<<endl;
  if(argc==2){
	cout<<"Impossible de resoudre le systeme : methode de resolution non specifiee !"<<endl;
	exit(1);
  }
  matrice lu;
  int max_iter=1000;
  const double eps = 1e-12;
  start = clock();
  switch(*argv[2]){
	case 'g' : u = a.solveur_gc(f,init,eps,max_iter);
		   break;
	case 'j' : u = a.solveur_jacobi(f,init,eps,max_iter);
		   break;
	case 'l' : lu = a.facto_lu();
		   u = a.solveur_lu(f,lu);
		   break;
	default  : cout<<"Methode inconnue ..."<<endl;
		   exit(1);
		   break;
  }

  //creation d'un nouveau repertoire de sauvegarde des resultats :
  char * COMMAND = (char *) malloc(sizeof(*COMMAND) * 256);
  char * directory_name = (char *) malloc(sizeof(*directory_name) * 128);
  strcpy(COMMAND,"mkdir results/");	//concatenation des 2 chaines dans COMMAND
  strtok(argv[1], "/");			//separation en 2 sous-chaine separees par le caratere '/'
  directory_name = strtok(NULL, "/"); 	//acces a la deuxieme sous-chaine : le nom du maillage
  for(int i=2;i<argc;i++)
  	strcat(directory_name, argv[i]);
  strcat(COMMAND, directory_name);
  system(COMMAND);
  free(COMMAND);

  fstream fichier;
		
	//sauvegarde du fichier dans le repertoire cree :
	char * file_name = (char *) malloc(sizeof(*file_name) * 128);
	char * buffer = (char *) malloc(sizeof(*buffer) * 128);
	strcpy(file_name, "results/");			//copie "results/" dans file_name
	strcat(file_name,directory_name);
	strcat(file_name, "/");
	if(*argv[2]!='l')
		sprintf(buffer, "%f.dat", eps);
	else
		strcpy(buffer, "exact.dat");
	strcat(file_name, buffer);
	fichier.open(file_name,fstream::out|fstream::trunc);  
	for(unsigned i=0; i<p.size(); i++)
		fichier <<u[i]<<endl;
	fichier.close();
	free(file_name);
	free(buffer);

  endtime = clock();
  cout<<"Temps de calcul : "<<(endtime-start)/(double) CLOCKS_PER_SEC<<endl;
  cout<<"************************************************"<<endl<<endl;

//======================//  
//  CALCUL DE L'ERREUR  //
//======================//

  cout<<"************************************************"<<endl;
  cout<<"Calcul de l'erreur..."<<endl;  

  double max_seg=0;
  double erreur=0;
  double x,y;
  int j;
  for(unsigned int i=0; i<t.size(); i++){  

	//on fait une boucle sur tous les triangles et on
	//utilise la 2eme formule de quadrature du cours

	j = t[i].tri[0];
        x = p[j].coord[0];
	y = p[j].coord[1];
	erreur+=(1/3.)*aire(t[i],p)*pow(s(x,y,10)-u[j],2);   //a la puissance 2 car norme L2 donc g est qqch au carre

	j = t[i].tri[1];
        x = p[j].coord[0];
	y = p[j].coord[1];
	erreur+=(1/3.)*aire(t[i],p)*pow(s(x,y,10)-u[j],2);

	j = t[i].tri[2];
        x = p[j].coord[0];
	y = p[j].coord[1];
	erreur+=(1/3.)*aire(t[i],p)*pow(s(x,y,10)-u[j],2);

	//on calcule la longueur max de tous les segments pour pouvoir calculer la vitesse de convergence apres
	double seg1=sqrt(pow(p[t[i].tri[1]].coord[0]-p[t[i].tri[0]].coord[0],2)+pow(p[t[i].tri[1]].coord[1]-p[t[i].tri[0]].coord[1],2));
	double seg2=sqrt(pow(p[t[i].tri[1]].coord[0]-p[t[i].tri[2]].coord[0],2)+pow(p[t[i].tri[1]].coord[1]-p[t[i].tri[2]].coord[1],2));
	double seg3=sqrt(pow(p[t[i].tri[0]].coord[0]-p[t[i].tri[2]].coord[0],2)+pow(p[t[i].tri[0]].coord[1]-p[t[i].tri[2]].coord[1],2));
	if(max_seg<seg1)max_seg=seg1;
	if(max_seg<seg2)max_seg=seg2;
	if(max_seg<seg3)max_seg=seg3;
  }

  //vitesse de convergence
  double c;
  int alp = 1;
  erreur = sqrt(erreur);
  c = erreur/pow(max_seg,alp);  //il faut tester avec les 4 maillages et voir quand c est constant

  cout<<"erreur : "<< erreur <<endl;
  cout<<"alpha : "<< alp <<endl;
  cout<<"C : "<< c <<endl;
   cout<<"************************************************"<<endl<<endl;

//=============//  
//  AFFICHAGE  //
//=============//

	cout<<"************************************************"<<endl;
	cout<<"Affichage ..."<<endl;
	cout<<"************************************************"<<endl<<endl;


//=====================//  
//  Initialisation de  //
//  l'affichage        //
//=====================//

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowSize(Width, Height);
  glutInitWindowPosition(500, 150);
  glutCreateWindow("Helmholtz");
  glutPopWindow();

//=========================//  
// Affichage et evenements // 
//=========================//

  SetView();
  glutReshapeFunc  ( Reshape     );  // si la fenetre graphique change 
  glutKeyboardFunc ( Key         );  // pour les evenements clavier
  glutMouseFunc    ( Mouse       );  // pour les evenements souris
  glutMotionFunc   ( MotionMouse );
  glutDisplayFunc  ( Display     );  // pour l'affichage
  glutMainLoop();               
  
}
