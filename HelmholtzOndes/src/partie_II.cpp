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

//======================//
// Fonction pour gerer  //
// la carte de couleur  //
//======================//

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
  double fmin=0;
  
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
    //glColor3f(0.,0.,0.);
    SetColor(u[I]);
    glVertex3f(X,Y,0.);       

    // 2eme sommet
    I = t[j].tri[1];
    X = p[I].coord[0], Y = p[I].coord[1];
    //glColor3f(0.,0.,0.);
    SetColor(u[I]);
    glVertex3f(X,Y,0.);

    // 3eme sommet
    I = t[j].tri[2];
    X = p[I].coord[0], Y = p[I].coord[1];
    //glColor3f(0.,0.,0.);
    SetColor(u[I]);
    glVertex3f(X,Y,0.);    
    
    glEnd();  
  }
 
    //Affichage de l'echelle de couleur :

    glBegin(GL_QUADS);
    SetColor(0.);
    glVertex3f(1.2,-0.4,0.0);
    SetColor(0.);
    glVertex3f(1.25,-0.4,0.0);
    SetColor(.33);
    glVertex3f(1.25,-0.2,0.0);
    SetColor(.33);
    glVertex3f(1.2,-0.2,0.0);
    glEnd();
    glBegin(GL_QUADS);
    SetColor(.33);
    glVertex3f(1.2,-0.2,0.0);
    SetColor(.33);
    glVertex3f(1.25,-0.2,0.0);
    SetColor(.66);
    glVertex3f(1.25,0.2,0.0);
    SetColor(.66);
    glVertex3f(1.2,0.2,0.0);
    glEnd();
    glBegin(GL_QUADS);
    SetColor(.66);
    glVertex3f(1.2,0.2,0.0);
    SetColor(.66);
    glVertex3f(1.25,0.2,0.0);
    SetColor(1.);
    glVertex3f(1.25,0.4,0.0);
    SetColor(1.);
    glVertex3f(1.20,0.4,0.0);
    glEnd();

    //affichage des unites sur le cote :
    SetColor(2.0);		//pour afficher en noir
    char * buf = (char *) malloc( sizeof(*buf) * 2);
    for(int i=0;i<10;i++){
    	glRasterPos2f(1.25,-.4+i*.8/10.);
    	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'_');
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,' ');
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'0');
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'.');
	sprintf(buf,"%d",i);
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*buf);
    }
    glRasterPos2f(1.25,.4);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'_');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,' ');
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'1');
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
// Initialisation des  //
// variables globales  //
//=====================//
  
  X0 = 1/4.; 
  Y0 = 0; 
  Z0 = 0;
  Width  = 500; 
  Rho    = 2.5; 
  Height = 500; 
  Phi   = M_PI/2.;
  Theta = -0.000001;		//si on met 0 on ne voit rien
  
  clock_t start;
  clock_t endtime;

  const double T = 2.25;
  double dt;
  int N;
  int t_snap_shot;
  double min_seg=1;

//===========================//  
//  IMPORTATION DU MAILLAGE  //
//===========================//
 
  start = clock();
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

	//on calcule la longueur minimale de tous les segments pour la CFL
  for(unsigned i=0;i<t.size();i++){
	double seg1=sqrt(pow(p[t[i].tri[1]].coord[0]-p[t[i].tri[0]].coord[0],2)+pow(p[t[i].tri[1]].coord[1]-p[t[i].tri[0]].coord[1],2));
	double seg2=sqrt(pow(p[t[i].tri[1]].coord[0]-p[t[i].tri[2]].coord[0],2)+pow(p[t[i].tri[1]].coord[1]-p[t[i].tri[2]].coord[1],2));
	double seg3=sqrt(pow(p[t[i].tri[0]].coord[0]-p[t[i].tri[2]].coord[0],2)+pow(p[t[i].tri[0]].coord[1]-p[t[i].tri[2]].coord[1],2));
	if(min_seg>seg1)min_seg=seg1;
	if(min_seg>seg2)min_seg=seg2;
	if(min_seg>seg3)min_seg=seg3;
  }

//=======================//
//  ASSEMBLAGE MATRICE   //
//=======================//
  
  matrice M;
  matrice AB;
  cout<<"************************************************"<<endl;
  cout<<"Debut assemblage matrice"<<endl;
  start = clock();
  build_mat(M, AB, p, e, t, argc, argv);
  endtime = clock();
  cout<<"Fin assemblage matrice"<<endl;
  cout<<"Temps d'assemblage : "<<(endtime-start)/(double) CLOCKS_PER_SEC<<endl;

  cout<<"************************************************"<<endl<<endl;

//===================//  
//  SCHEMA EN TEMPS  //
//===================//

  int max_iter=1000;
  u.resize(p.size(),0);
  vect init;
  vect u_0;
  vect v;
  vect v_0;
  vect bu;
  vect bv;
  vect buf;
  vect buf2;
  vect u2, v2;

  u_0.resize(p.size());
  v.resize(p.size());
  v_0.resize(p.size());
  bu.resize(p.size());
  bv.resize(p.size());
  buf.resize(p.size());
  buf2.resize(p.size());
  u2.resize(p.size());
  v2.resize(p.size());

  dt = atof(argv[2]);
  t_snap_shot = atoi(argv[3]);	//k*dt ou on sauve avec k entre en argument du programme
  N = T/dt;
  cout<<"************************************************"<<endl;
  cout<<"Debut du schema en temps :"<<endl;
  start = clock();
	
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

  //recuperation du numero du point associe a l'origine :
  int index_or;
  double k = 5;
  for(unsigned i=0;i<p.size();i++){
	if(k != min(k, pow(p[i].coord[0],2)+pow(p[i].coord[1],2)))
		   index_or=i;
        k = min(k, pow(p[i].coord[0],2)+pow(p[i].coord[1],2));
  }

  //creation d'un fichier ou sauvegarder uh(0,0) :
  fstream fichier_uh;
  char * file_uh_name = (char *) malloc(sizeof(*file_uh_name) * 128);
  char * aux = (char *) malloc(sizeof(*aux) * 128);
  strcpy(file_uh_name,"results/uh(0,0)/");
  sprintf(aux,"%f%c.dat",atof(argv[2]), *argv[4]);
  strcat(file_uh_name,aux);
  fichier_uh.open(file_uh_name,fstream::out|fstream::trunc);		//on efface le contenu deja existant
	
  //initialisation :		
  for(unsigned i=0;i<p.size();i++)
	u_0[i] = exp(-100*(pow(p[i].coord[0]-.1,2)+pow(p[i].coord[1]-.1,2)));
  
  init.resize(p.size(),0);
  //schema general :
  for(int n=1;n<N;n++){

	buf = M*u_0;
	v2 = v_0*dt;
	buf += M*v2;
	u2 = u_0*(-(dt*dt)/2.);
	buf += AB*u2;
	bu = buf;
	buf.clear();
	u = M.solveur_gc(bu,init,1e-12, max_iter);

	buf = u+u_0;
	buf = buf*(-dt/2.);
	buf2 = AB*buf;
	v2 = M*v_0;
	buf2 += v2;
	bv = buf2;
	v = M.solveur_gc(bv,init,1e-12, max_iter);

	
	u_0=u;
	v_0=v;	
	buf.clear();
	buf2.clear();
	
	//Sauvegarde uh(0,0) :
	fichier_uh << n*dt << "\t" << u[index_or] << endl;

	//Sauvegarde dans un fichier :
	if( (n%t_snap_shot) == 0 ){
		fstream fichier;
		
		//sauvegarde du fichier dans le repertoire cree :
		char * file_name = (char *) malloc(sizeof(*file_name) * 128);
		char * buffer = (char *) malloc(sizeof(*buffer) * 128);
		strcpy(file_name, "results/");			//copie "results/" dans file_name
		strcat(file_name,directory_name);
		strcat(file_name, "/");
		sprintf(buffer, "%f.dat", n*dt);
		strcat(file_name, buffer);
		fichier.open(file_name,fstream::out|fstream::trunc);  
		for(unsigned i=0; i<p.size(); i++)
			fichier <<u[i]<<"\t"<<v[i]<<endl;
		fichier.close();
		free(file_name);
		free(buffer);	
	}
  }
  fichier_uh.close();
  free(aux);
  free(file_uh_name);
  endtime = clock();
  cout<<"Temps de calcul : "<<(endtime-start)/(double) CLOCKS_PER_SEC<<endl;
  cout<<"CFL dt<"<<min_seg/2<<endl;
  cout<<"************************************************"<<endl<<endl;

//=============//  
//  AFFICHAGE  //
//=============//

	cout<<"************************************************"<<endl;
	cout<<"Affichage ..."<<endl;
	cout<<"************************************************"<<endl<<endl;


//=====================//  
// Initialisation de   //
// l'affichage         //
//=====================//

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowSize(Width, Height);
  glutInitWindowPosition(500, 150);
  glutCreateWindow("Helmotz");
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
	
  free(directory_name);
}
