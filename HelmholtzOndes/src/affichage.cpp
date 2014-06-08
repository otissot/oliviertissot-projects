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
#include <gmm/gmm.h>

using namespace std;
using namespace gmm;

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

vector< double > u;
vector< double > v;
vector< triangle > t;
vector< point > p;
vector< edge > e;

//=====================//
// Fonction pour gerer //
// la carte de couleur //
//=====================//

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

//=====================//  
//  TRACE DU MAILLAGE  //
//=====================//  
  
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


//====================//  
// Initialisation des //
// variables globales //
//====================//
  
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

  double t_aff;			//l'instant affiche
  ifstream fichier;

  cout << "************************************************" << endl;
  cout << "Bienvenue sur la routine d'affichage des resultats !" << endl << endl;
  cout << "Parametres de la simulation :"<<endl;
  cout << "Maillage : "<< argv[1] << endl;
  cout << "dt : "<< argv[2] << endl;
  cout << "t_snap_shot : "<<argv[3]<<endl;
  cout << "Mass-lumping : "<<argv[4]<<endl;
  cout << "Quel instant voulez-vous voir? (ce doit etre un multiple de t_snap_shot)"<<endl;
  cin >> t_aff;
  cout << "************************************************" << endl << endl;


  cout<<"************************************************"<<endl;
  cout<<"Importation du maillage ..."<<endl;
  start = clock();
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
  
  u.resize(p.size());
  v.resize(p.size());

  cout<<"************************************************"<<endl<<endl;
  cout<<"Recuperation des donnes ..."<<endl;
  start = clock();
  char * directory_name = (char *)malloc( sizeof(*directory_name) * 128);
  strtok(argv[1], "/");					//separation en 2 sous-chaine separees par le caratere '/'
  directory_name = strtok(NULL, "/"); 			//acces a la deuxieme sous-chaine : le nom du maillage
  for(int i=2;i<argc;i++)
  	strcat(directory_name, argv[i]);
  char * file_name = (char *) malloc(sizeof(*file_name) * 128);
  char * buffer = (char *) malloc(sizeof(*buffer) * 128);
  strcpy(file_name, "results/");			//copie "results/" dans file_name
  strcat(file_name,directory_name);
  strcat(file_name, "/");
  sprintf(buffer, "%f.dat", t_aff);
  strcat(file_name, buffer);
  cout << endl;
  fichier.open(file_name,ifstream::in);  
  for(unsigned i=0; i<p.size(); i++){
	fichier>>u[i];
	fichier>>v[i];
  }
  fichier.close();
  free(file_name);
  free(buffer);
  endtime = clock();
  cout<<"Temps de recuperation : "<<(endtime-start)/(double) CLOCKS_PER_SEC<<endl;
  cout<<"************************************************"<<endl<<endl;

//=============//  
//  AFFICHAGE  //
//=============//

	cout<<"************************************************"<<endl;
	cout<<"Affichage ..."<<endl;
	cout<<"************************************************"<<endl<<endl;


//===================//  
// Initialisation de //
// l'affichage       //
//===================//

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowSize(Width, Height);
  glutInitWindowPosition(500, 150);
  glutCreateWindow("Equation des ondes");
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
