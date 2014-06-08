#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

#include <fstream>
#include <iostream>
#include <cstring>

#include <vector>

#include "Maillage.hpp"
//#include "RNM.hpp"

const R pi=M_PI; 
using namespace std;

typedef  vector<R> Rn;
const int TheDrawList=1;  // numero de la list d'affichage 
bool   TheDrawListCreate=false;

class Global {    
// une petite classe pour stoker toutes les variables globales
public:
  int Width , Height;                  // taille de l'ecran en pixel
  Maillage &Th;                        // maillage courant
  bool stereo;
  bool windowdump;
  R rapz;
  int nbiso;
  R * viso;                            // les valeurs des iso values
  Rn &f ;                              // solution P1 a afficher
  R theta,phi,coef_dist;               // coordonnee polaire de la camera
  R dtheta;                            // vitesse de rotation de la camera 
  R xmin,xmax,ymin,ymax,zmin,zmax;     // borne de la scene
  R xm,ym,zm;                          // point regarde
  R ceyes;
  R focal;
  int xold,yold;
  Global(Maillage & TTh,Rn &ff,int height,int width,R rpz,int nbisovalue,bool st) ;
  void SetView(int eyes=0) const;      // defini le point de vue
  void DefaultView();
  void MoveXView(R dx,R dy);
  //  void MakeListDraw() const  ;     // construit la list d'affichage
  void Affiche();
} *global ;                            // la variable globale


Global::Global(Maillage & TTh,Rn &ff,int height,int width,R rpz,int nbisovalue,bool st) 
  : Th(TTh), stereo(st),f(ff)
{
  nbiso = nbisovalue;
  Width= width ;
  Height=height;
  rapz=rpz;
  windowdump=false;
  // first compute the mesh bound
  const  Sommet & v0=Th(0);
  xmin=xmax=v0.x;
  ymin=ymax=v0.y;
  zmin = f[0], zmax=f[0];
  for (int i=0;i<Th.nv;i++)
    {
      const  Sommet  & v=Th(i);
      xmin= min(xmin,v.x);
      ymin= min(ymin,v.y);
      xmax= max(xmax,v.x);
      ymax= max(ymax,v.y); 
      zmin= min(zmin,f[i]);
      zmax= max(zmax,f[i]);                   
    }
  if(nbiso>2)
    {
      viso = new R[nbiso];
      R diso=(zmax-zmin)/(nbiso-1);
      for (int i=0;i<nbiso;i++)
	viso[i]=zmin+i*diso;
    }
  else 
    viso=0;
  
   DefaultView();
  
}

void Global::DefaultView() 
{
  coef_dist = 1.;
  theta = -90*pi/180.;
  dtheta = 0.; // 0.1 degree  si trop rapide
  focal = 30; // en degree
  ceyes = 0.05; // de l'objets
  coef_dist=1;
  phi = 90*pi/180.;
  xm = (xmin+xmax)*0.5;
  ym = (ymin+ymax)*0.5;
  zm = rapz*(zmin+zmax)*0.5;    
}

void Global::SetView(int eye) const
{
  glViewport( 0, 0, Width, Height );
  
  R ratio= (double) Width / (double)  Height; 
  R dx =(xmax-xmin), dy= (ymax-ymin), dz=(zmax-zmin)*rapz;
  R dmax= sqrt(dx*dx+dy*dy+dz*dz);
  R dist = dmax/2/asin(focal/2*pi/180.)*coef_dist;
  R camx=xm+cos(phi)*cos(theta)*dist;
  R camy=ym+cos(phi)*sin(theta)*dist;
  R camz=zm+dist*sin(phi);  
  R znear=max(dist-dmax/2,1e-5);
  R zfare=dist+dmax/2;
  R aspect=ratio;
  if (eye)
    {
      R dmm = -dmax*ceyes;
      R dx = -dmm*sin(theta);
      R dy = dmm*cos(theta);
      camx += dx*eye;
      camy += dy*eye;
    }   
  
  //  R hx= (  ratio*dy < dx  ) ? dx : dy*ratio ;
  //  R hy= hx/ratio ;
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); 
  gluLookAt(camx,camy,camz,xm,ym,zm,0.,0.,1.);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); 
  gluPerspective(focal,aspect,znear,zfare);
  
  
}

void Global::MoveXView(R dx,R dy)
{
  // cout << xm << " " << ym << " " << zm << " apres " ;
  zm -= dy*(zmax-zmin)*rapz/50.;
  xm += dx*(xmax-xmin)*sin(theta)/50;
  ym -= dx*(ymax-ymin)*cos(theta)/50;   
  // cout << xm << " " << ym << " " << zm << endl;
}

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

void   SetColor(R f)
{
  float r=0,g=0,b=0; 
  assert(global);
  R fmin=global->zmin; // borne de la fonction
  R fmax=global->zmax;
  hsvToRgb(0.99*(f-fmin)/(fmax-fmin),1,1,r,g,b);
  glColor3f(r,g,b);
}

void DrawVertex(const  R2 & v,R z=0,R rapz=1) 
{  
  SetColor(z);             // la couleur
  glVertex2f(v.x, v.y);/*, z*rapz);*/ // le sommet
}

void DrawIsoTfill(const R2 Pt[3],const R ff[3],const R * Viso,int NbIso, R rapz=1)
{
  R2 PQ[10];
  R z[10];
  
  R eps= (Viso[NbIso-1]-Viso[0])*1e-6;
  for(int l=1;l< NbIso;l++)  //   loop on the level curves 
    {
      R xfb = Viso[l-1];
      R xfh = Viso[l];
      assert(xfb < xfh);
      int im=0;
      for(int i=0;i<3;i++) // for the  3 edges 
	{
          int j=(i+1)%3;
          R fi=(ff[i]);
          R fj=(ff[j]);
          R xxfb =  xfb;
          R xxfh =  xfh;
	  if (fj<fi ) swap(xxfb,xxfh);
          R xf  = xxfb;
	  if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
	    {
	      if (fabs(fi-fj)>=0.1e-20)
		{
		  R  xlam=(fi-xf)/(fi-fj);
		  z[im] =  ff[i] * (1.F-xlam)  +  ff[j]* xlam;
		  PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
		}
	    }
          xf = xxfh;  
	  if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
	    {
	      if (fabs(fi-fj)>=0.1e-20)
		{
		  R  xlam=(fi-xf)/(fi-fj);
		  z[im] =  ff[i] * (1.F-xlam)  +  ff[j]* xlam;
		  PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
		}
	    }
	  if (  xfb-eps <=fj  && fj <= xfh+eps) 
	    z[im]=ff[j],PQ[im++] = Pt[j];
	  
	}
      if (im>2) 
	{
          glBegin(GL_POLYGON);
          SetColor((xfb+xfh)/2); 
	  for (int i=0;i<im;i++)
	    {// cout << i << " \t : " << PQ[i].x << " " <<  PQ[i].y << " " << z[i]*rapz << endl;
	      glVertex2f(PQ[i].x, PQ[i].y);/*,z[i]*rapz);*/
	    }
          glEnd();
	  
	}
    }
} 

bool WindowDump(int width,int height)
{
	int i,j;
	FILE *fptr;
	static int counter = 0;
	char fname[32];
	unsigned char *image;
	
	/* Allocate our buffer for the image */
	if ((image = new unsigned char[3*width*height]) == NULL) {
		fprintf(stderr,"WindowDump - Failed to allocate memory for image\n");
		return(false);
	}
	
	/* Open the file */
	sprintf(fname,"L_%04d.ppm",counter);
	if ((fptr = fopen(fname,"w")) == NULL) {
		fprintf(stderr,"WindowDump - Failed to open file for window dump\n");
		return(false);
	}
	cout << " WindowDump in " << fname << endl;
	/* Copy the image into our buffer */
	glReadBuffer(GL_FRONT);
	glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);
	
	/* Write the PPM file */
	fprintf(fptr,"P6\n%d %d\n255\n",width,height);
	for (j=height-1;j>=0;j--) {
		for (i=0;i<width;i++) {
			fputc(image[3*j*width+3*i+0],fptr);
			fputc(image[3*j*width+3*i+1],fptr);
			fputc(image[3*j*width+3*i+2],fptr);
		}
	}
	fclose(fptr);
	
	delete [] image;
	counter++;
	return(true);
}

void Clean() 
{
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

static void Reshape( int width, int height )
{   
  global->Width  = width;
  global->Height = height;  
  glutPostRedisplay();
}

static void Key( unsigned char key, int x, int y )
{
  switch (key) {
  case 27: // esc char
    exit(0);
    break;
  case '+':  
    global->coef_dist /= 1.2;
    break;
  case '-':  
    global->coef_dist *= 1.2;
    break;
  case 'g':  
    global->theta += pi/180.;
    break;
  case 'd':  
    global->theta -= pi/180.;
    break;
  case 'h':  
    global->phi += pi/180.;
    break;
  case 'b':  
    global->phi -= pi/180.;
    break;
  case 'a':
    global->dtheta = global->dtheta ? 0 : pi/1800.;
    break;	
  case 'A':
    global->dtheta *=1.2;
    break;	
  case 's':
    global->dtheta /= 1.2;
    break;	
  case 'f':
    global->focal /= 1.2 ;
    break;	
  case 'F':
    global->focal *= 1.2 ;
    break;	
  case 'e':
    global->ceyes /= 1.2 ;
    break;	
  case 'E':
    global->ceyes *= 1.2 ;
    break;
  case 'p':
    cout << " ceyes = " << global->ceyes << " " << " focal = " << global->focal 
	 << " coef_dist " << global->coef_dist  << endl;
    break;	
   case 'w':
	   global->windowdump=true;
	   break;
  case '=':
    global->DefaultView();
    break;
  default:
    cout << " Key Character " << (int) key << " " << key << endl;  
    cout << 
        "=  defautl \n"
        "Ee  +/- ecart en les oeils\n"
        "fF  +/ focal \n"
        "s   slow \n"
        "A    accelere \n"
        "a   animation ou arret\n"
        "bhgd  bas haut gauche droite \n"
	"w     window dump \n"
        "+-   zoom \n"
        "(ESC) FIN \n"
        ;
        
    
  }
  glutPostRedisplay();
}

void Global::Affiche()
{
  if(TheDrawListCreate)
    glCallList(1);
  else 
    {
      glNewList(1,GL_COMPILE_AND_EXECUTE);
      glPolygonMode(GL_FRONT,GL_FILL); // mode affichage des polygones  
      // constructions des triangles colore
      if (nbiso)
	{
	  for (int i=0;i<Th.nt;i++)
	    {
	      const Triangle & K(Th[i]);
	      R2 Pt[3]={K[0],K[1],K[2]};
	      R ff[3]={f[Th(K[0])],f[Th(K[1])],f[Th(K[2])]};
	      DrawIsoTfill(Pt,ff,viso,nbiso,rapz);
	    }
	}
      else
	{
	  for (int i=0;i<Th.nt;i++)
	    {
	      const Triangle & K(Th[i]); 
	      int i0= Th(K[0]),  i1= Th(K[1]),   i2= Th(K[2]) ;    
	      glBegin(GL_TRIANGLES);
	      DrawVertex(K[0],f[i0],rapz);
	      DrawVertex(K[1],f[i1],rapz);
	      DrawVertex(K[2],f[i2],rapz);
	      glEnd();
	    }
	}
      glEndList();
      TheDrawListCreate=true;
    }
}

void Display(void)
{ 

  if (global->stereo)
    {
      glClearColor(1.0, 1.0, 1.0, 0.0);
      glDrawBuffer(GL_BACK_RIGHT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      global->SetView(-1);
      //glCallList(TheDrawList);    
      global->Affiche();
      glClearColor(1.0, 1.0, 1.0, 0.0);
      glDrawBuffer(GL_BACK_LEFT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      global->SetView(+1);
      //glCallList(TheDrawList);
      global->Affiche();
      glFlush();
      glutSwapBuffers();
    }
  else 
    {
      Clean();
      global->SetView();
      //glCallList(TheDrawList);    
      global->Affiche();
      glFlush();
      glutSwapBuffers();
	  if ( global->windowdump)
		  WindowDump(global->Width,global->Height);
	  global->windowdump=false;
    }
}

static void Idle( void )
{
  if (global->dtheta)
    {
      global->theta += global->dtheta;
      glutPostRedisplay();
    }
  
}
static void Mouse( int button,int state,int x,int y )
{
  // state up or down 
  if (state == GLUT_DOWN) { global->xold=x,global->yold=y;return;}
  // cout << "Mouse " << button<< " " << state << " " << x-global->xold << " " << y-global->yold << endl;
  //  x gauche -> droite
  //  y  haut -> bas`
  global->phi   += (y-global->yold)/(2.*180.);
  global->theta -= (x-global->xold)/(2*180.);
  glutPostRedisplay();
  
}

static void MotionMouse(int x,int y )
{
  // cout << " MotionMouse " << " " << x << " " << y << endl;
  // GLuint gtime = glutGet(GLUT_ELAPSED_TIME); //   
  global->phi   += (y-global->yold)/(2.*180.);
  global->theta -= (x-global->xold)/(2*180.);
  global->xold=x;
  global->yold=y;
  glutPostRedisplay();
}

void SpecialKey(int key, int x, int y)
{
 // cout << " SpecialKey " << key << " " << x << " " << y << " : ";
  R dx(0),dy(0);
  switch (key) 
    {
    case  GLUT_KEY_LEFT:   dx = -1; break;
    case  GLUT_KEY_RIGHT:  dx = +1; break;
    case  GLUT_KEY_DOWN:   dy = -1; break;
    case  GLUT_KEY_UP:     dy = +1; break;
    }
  // calcul du deplacement de xm,ym,zm;
  // cout << " " << dx << " " << dy << endl;
  global->MoveXView(dx,dy);
  glutPostRedisplay();
}


int main(int argc, char** argv)
{
  if (argc <3)
   {
      cerr << " utilisation : " << argv[0] << " meshfile solfile  [rap in z ] [nbisovalue] " << endl;
      return 1;
   }
  int imesh = 1, isol=2;
  bool stereo=false;
  bool fullscreen = false;
  assert(argc>2);
  int karg=1;
  for (int k=0;karg!=k;)
   {
     k=karg;
     cout << "arg [" << karg << "] = " << argv[karg] ;
       cout << " stereo =" << stereo << " fullscreen =" << fullscreen << endl;

    if( strcmp("-s",argv[karg])==0) 
     { stereo=1; karg++;}
    if (strcmp("-f",argv[karg])==0) 
     { fullscreen=1; karg++;}
     cout << " stereo =" << stereo << " fullscreen =" << fullscreen << endl;
   }
  cout << " stereo =" << stereo << " fullscreen =" << fullscreen << endl;
  imesh = karg++;
  isol = karg++;
  cout << argc << " "<< karg << endl;
  assert(argc>=karg);
  R rapz=1; 
  int nbiso=25;
  if (argc>karg) rapz=atof(argv[karg++]);
  if (argc>karg) nbiso=atoi(argv[karg++]);
  cout << " Rap z " << rapz << endl;
 
  
  // On lit la fonction a afficher
  Maillage Th(argv[imesh]);
  Rn f(Th.nv); 
  {
    ifstream fdat(argv[isol]);
    assert(fdat.good());
    int n; 
    fdat >> n ; 
    f.resize(n); 
    for(int i=0;i <n; ++i)
    fdat >> f[i];
  } // pour ferme le fichier (la variable fdat est detruite)
  glutInit(&argc, argv);
  
  if(stereo)
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);
  else  
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  
  int Height = 512;
  int Width = 512; 
  glutInitWindowSize(Width , Height);
  glutInitWindowPosition(100, 100);

  string titre = "vue de ";
  titre += argv[imesh] ;
  titre += ", ";
  titre += argv[isol];
  glutCreateWindow(titre.c_str());
  glutPushWindow();
  if (fullscreen)
    glutFullScreen();

  global=new Global(Th,f,Height,Width,rapz,nbiso,stereo);
  //global->MakeListDraw();    
  
  glEnable(GL_DEPTH_TEST); 
  glutReshapeFunc( Reshape ); // pour changement de fenetre 
  glutKeyboardFunc( Key );    // pour les evenements clavier
  glutSpecialFunc(SpecialKey);
  glutMouseFunc(Mouse);       // pour les evenements sourie
  glutMotionFunc(MotionMouse); // les mouvements  de la sourie 
  glutDisplayFunc( Display ); // l'affichage
  glutIdleFunc( Idle );       // l'animation automatique
  
  glutMainLoop(); 
  return 0;
}


