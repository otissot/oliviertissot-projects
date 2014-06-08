// Maillage.cpp --- 
// 
// Filename: Maillage.cpp
// Description: 
// Author: O.TISSOT E.NAYIR
// Created: ven. déc. 13 09:59:53 2013 (+0100)
// Version: 
// Last-Updated: ven. mars 21 09:36:55 2014 (+0100)
//           By: Olivier
//     Update #: 329
// Compatibility: 
// 
// 

// Commentary: 
// 
// 
// 
// 

// Code:
 
#include "Maillage.hpp"

// Initialisation des membres statiques
static const int  nvSegTria[3][2] = { {1,2},{2,0},{0,1} };
const int (* const Triangle::nvadj)[2] = nvSegTria ;
int Maillage::ksearch = 0;
int Maillage::kthrough = 0;


Maillage::Maillage(const string& filename) 
  : TheAdjacencesLink(0)
{
  // Check sur le format du fichier de maillage (p-e a mettre dans le main) :
  unsigned pos = filename.find(".");
  string meshFormat = filename.substr(pos);   // on recupere la partie a droite du . dans le fichier de maillage
  cout << "Format maillage : " << meshFormat << endl;
  if(meshFormat.compare(".msh")){
    cout << "Mauvais format de maillage... Seuls les .msh sont acceptes ;-)" << endl;
    exit(1);         // on jette l'user si l'extension n'est pas .msh (eventuellement on pourrait ajouter un autre lecteur dans ce cas)
  }
  ifstream f(filename.c_str());
  assert(f);
  int notused;
  int isommet[3];
  f >> nv >> nt >> ne;
  t = new Triangle[nt]; // utilise le constructeur par defaut de la class Triangle
  be = new Segment[ne]; // ---------------------------------------------- Segment
  s = new Sommet[nv];   // ---------------------------------------------- Sommet
  assert(t && s);
  // On lit les sommmets :
  for( int i=0;i<nv;++i)
    {
      f >> s[i].x >> s[i].y >> s[i].label;
      s[i].numero = i;
    }
  double mes=0;
  // On lit les triangles :
  for( int i=0;i<nt;++i)
    {
      for(int l=0;l<3;++l){
	f >> isommet[l];                           // on lit les 3 numeros des sommets du triangle
        isommet[l]--;                              // on saute le label (la 4eme valeur) car on ne l'utilise pas et sinon on renormalise les numero de
      }	                                           // sommet a cause des convention de numerotation des tableaux en C++
      f >> notused;
      // On "construit" le triangle i
      t[i].init(s,isommet);
      mes += t[i].mesure();
    }
  double mesgamma = 0;
  // On lit les segments du bord :
  for(int i=0;i<ne;++i){
    for(int l=0;l<2;++l){
	f >> isommet[l];                           // on lit les 2 numeros des sommets du segment
        isommet[l]--;                              // on saute le label (la 3eme valeur) car on ne l'utilise pas et sinon on renormalise les numero de
      }	                                           // sommet a cause des convention de numerotation des tableaux en C++
    f >> notused;
    // On "construit" le segment i
    be[i].init(s,isommet);
    mesgamma += be[i].mesure();
  }
  // Petit check si tout est ok :
  cout << "Fin de la lecture du maillage !" << endl;
  cout << "nv = " << nv << " nt = " << nt << " ne = " << ne << endl;
  cout << "mesure d'omega = " << mes << endl;
  cout << "mesure de gamma = " << mesgamma << endl;
}

void Maillage::buildEdge(map< pair<int,int>, int>& edge) const{
  int neinside = -1;
  for(int k=0;k<nt;++k){
    for(int a=0;a<3;++a)
      {
        int i0 = (a+1)%3;
        int i1 = (a+2)%3;
        int j0 = t[k][i0], j1 = t[k][i1];
        if( j0 > j1) swap(j0,j1);
	pair<int,int> j01(j0,j1); 
	if( edge.find(j01) == edge.end() )
	  edge[j01] = ++neinside; 
      }
  }
  cout << " neinside =" << neinside+1 << endl;
  cout << " nb of  connex comp.  - nb of holes :  " << nv - neinside-1 + nt << endl;
}


// Pour associer un numero global aux milieux des segments
int Maillage::numedge(int k, map< pair<int,int>, int >& edge, int i) const{
  int res = -1;
  map< pair<int,int>, int >::iterator iter;
  // Le segment que l'on associe a i : point milieu oppose au sommet i%3 mais les segments sont entres de maniere ordonnee pour garantir que (s1,s2) = (s2, s1)
  pair<int,int> segment( min(t[k][(i+2)%3],t[k][(i+1)%3]), max(t[k][(i+2)%3],t[k][(i+1)%3]) );
  iter = edge.find(segment);
  res = iter->second + nv;
  return res;
};

// Pour faire la correspondance entre numerotation locale et numerotation globale
int Maillage::numglob(int k, map< pair<int,int>, int >& edge, int i) const{
  int res = -1;
  // U1 :
  if(i<3)
    res = t[k][i];                   // Cast d'un sommet en entier (son numero)
  else if(i>=3 && i<6)
    res = numedge(k,edge,i);
  // U2 :
  else if(i>=6 && i<9)
    res = t[k][i%3] + nv + edge.size();
  else if(i>=9 && i<12)
    res = numedge(k,edge,i) + nv + edge.size();
  // P :
  else if(i>=12 && i<15)
    res = t[k][i%3] + 2*(nv + edge.size()); 
  return res;
};



void Maillage::BuildAdj()
{
  // const int nva   = T::nva;
  // const int nea   = T::nea;
  if(TheAdjacencesLink!=0) return ;//  already build ...
  TheAdjacencesLink = new int[nea*nt];
  BoundaryElementHeadLink = new int[ne];
  HashTable<SortArray<int,nva>,int> h(nea*nt,nv);
  int nk=0,nba=0;
  int err=0;
  //if(verbosity>5)
  //  cout << "   -- BuildAdj:nva=// nea=" << nva << " " << nea << " "<< ne << endl;
  for (int k=0;k<nt;++k)
    for (int i=0;i<nea;++i)
      {
        SortArray<int,nva> a(itemadj(k,i));
        //cout << " ### "   << " item(k,i)= " << itemadj(k,i) << " a= " << a << " k " << k << " i " << i << endl;
        typename HashTable<SortArray<int,nva>,int>::iterator p= h.find(a);
        if(!p)
          {
            h.add(a,nk);
            TheAdjacencesLink[nk]=-1;
            nba++;
          }
        else
          {      
            assert(p->v>=0);
            TheAdjacencesLink[nk]=p->v;
            TheAdjacencesLink[p->v]=nk;
            p->v=-1-nk;
            nba--;
          }
        ++nk;
      }
   
  for (int k=0;k<ne;++k)
     {
        SortArray<int,nva> a(itembe(k));

        typename HashTable<SortArray<int,nva>,int>::iterator p= h.find(a);
        //cout << k << " ### "   << " item(k,i)= " << itembe(k) << " a= " << a << endl;
        if(!p) { err++;
          if(err==1) cerr << "Err  Border element not in mesh \n";
          if (err<10)  cerr << " \t " << k << " " << a << endl;
        }
	else
	  {
	    BoundaryElementHeadLink[k] = p->v <0 ? -p->v-1 : p->v;
          #ifndef NDEBUG
	    int tt=BoundaryElementHeadLink[k]/nea;
	    int ee=BoundaryElementHeadLink[k]%nea;
	    //cout << k << " ### "   << a << " = " << itemadj(t,e) << " t " << t << " e " << e << endl;
	    assert(itemadj(tt,ee)==a);
          #endif
	  }
     }

  assert(err==0);
  int na= h.n;
  cout << "  -- BuildAdj: nb Elememt " << nt << " nb vertices " << nv << endl;
  cout << "             : nb adj  = "<< na << " on border " << nba << " nea = " << nea << " nva = " << nva << endl;  
};



const Triangle *  Maillage::Find( R2 P, R2 & Phat,bool & outside,const Triangle * tstart) const
{
  int it,j;
  if ( tstart )
    it =  (*this)(tstart);
  else  {  
    assert(0); 
    // bug not implementa
    /*
    const Vertex * v=quadtree->NearestVertexWithNormal(P);
    if (!v) {
      v=quadtree->NearestVertex(P);
      assert(v);
    }
   
    it=Contening(v);
    */
  }
 
  int iit=-1;
  R delta=-1;
  R2 Phatt;
  int k=0;    
  ksearch++;
  while (1)
    {
      //loop:
      kthrough++;
      if (k++>=1000)
        {
          assert(k++<1000);
        }
      int kk,n=0,nl[3];
     
      const Triangle & K(t[it]);
     
      const R2 & A(K[0]), & B(K[1]), & C(K[2]);
      R l[3] = {0,0,0};
      R area2 = K.mesure()*2;
      R eps = -area2*1e-6;
      l[0] = det(P,B,C);
      l[1] = det(A,P,C);
      l[2] = area2-l[0]-l[1];
      if (l[0] < eps) nl[n++]=0;
      if (l[1] < eps) nl[n++]=1;
      if (l[2] < eps) nl[n++]=2;
      if (n==0)
        {
          outside=false;
          Phat=R2(l[1]/area2,l[2]/area2);
          return &K;
        }
      else if (n==1)
        j=nl[0];
      else  
        {
          kk=BinaryRand() ? 1 : 0;
          j= nl[ kk ];
        }
     
      //int jj  = j;
      int itt =  ElementAdj(it,j);
     
      if(itt==it || itt <0)  
        {
          if ( n==2 )
            {
              /*jj = */j = nl[ 1-kk ];
              itt =  ElementAdj(it,j);                  
              if (itt && itt != it)
                {
                  it=itt;
                  continue;
                }
            }
          // projection du point sur la frontiere
          l[nl[0]]=0;
          if(n==2) l[nl[1]]=0;
          R ll=l[0]+l[1]+l[2];
          Phat=R2(l[1]/ll,l[2]/ll);
          R2 PQ(K(Phat),P);
          R dd=(PQ,PQ);
          /*if (dd>delta && iit>=0)          
            {
              Phat=Phatt;
              return t+iit;
	      }*/
          /*
          int j0=(j+1)%3,i0= &K[j0]-s;
          int j1=(j+2)%3,i1= &K[j1]-s;
          int ii=-1,jj;
          if ( l[j0]> ll/2 ) ii=i0,jj=i1;
          if ( l[j1]> ll/2 ) ii=i1,jj=i0;
	  // cout << ii << " " << jj << " it = " << it << " " << delta << " " << dd <<  endl;
          //  pour la gestion de la frontiere
          //  on se promene sur la frontiere
          
          if (ii>0 && iib != ii )
            for (int p=BoundaryAdjacencesHead[ii];p>=0;p=BoundaryAdjacencesLink[p])
              { int e=p/2, ie=p%2, je=2-ie;
                // cout << number(bedges[e][0]) << " " << number(bedges[e][1]) << endl;
                if (! bedges[e].in( s+jj))
                  {  
                    iib = ii;
                    iit=it;
                    delta=dd;
                    Phatt=Phat;
                    it= BoundaryElement(e,ie);                      
                    // cout << "  ------ " << it << " " << Phatt <<  endl;
                    goto loop;
                  }
              }
          */
          outside=true;
          if (dd>delta && iit>=0)          
            {
	      Phat=Phatt;
	      return t+iit;
	    }
          else
            return t+it;
        }
      it=itt;
    }
};

int  WalkInTriangle(const Maillage & Th,
		    int it, 
		    double *lambda,
		    R u, 
		    R v, 
		    R & dt)
{
  const Triangle & T(Th[it]);
  const R2 Q[3]={ R2(T[0]),R2(T[1]),R2(T[2]) };
  
  R2 P  = lambda[0]*Q[0]  + lambda[1]*Q[1]  + lambda[2]*Q[2];
  
  //  cout << " " << u << " " << v ;
  R2 PF = P + R2(u,v)*dt;
  
  R l[3];
  l[0] = det(PF,Q[1],Q[2]);
  l[1] = det(Q[0],PF,Q[2]);
  l[2] = det(Q[0],Q[1],PF);
  R Det = l[0]+l[1]+l[2];
  l[0] /= Det;
  l[1] /= Det;
  l[2] /= Det;
  const R eps = 1e-5;
  int neg[3],k=0;
  int kk=-1;
  if (l[0]>-eps && l[1]>-eps && l[2]>-eps) 
    {
      dt =0;
      lambda[0] = l[0];
      lambda[1] = l[1];
      lambda[2] = l[2];
    }
  else 
    {
      if (l[0]<eps && lambda[0] != l[0]) neg[k++]=0;
      if (l[1]<eps && lambda[1] != l[1]) neg[k++]=1;
      if (l[2]<eps && lambda[2] != l[2]) neg[k++]=2;
      R eps1 = T.mesure()     * 1.e-5;
      if (k==2) // 2  
	{
	  // let j be the vertex beetween the 2 edges 
	  int j = 3-neg[0]-neg[1];
	  R S = det(P,PF,Q[j]);
	  
	  if (S>eps1)
	    kk = (j+1)%3;
	  else if (S<-eps1)
	    kk = (j+2)%3;
	  else if (BinaryRand())
	    kk = (j+1)%3;
	  else
	    kk = (j+2)%3;
	  
	} 
      else if (k==1)
	kk = neg[0];
      if(kk>=0)
	{
	  R d=lambda[kk]-l[kk];
	  
	  assert(d);
	  R coef =  lambda[kk]/d;
	  R coef1 = 1-coef;
	  dt        = dt*coef1;
	  lambda[0] = lambda[0]*coef1 + coef *l[0];
	  lambda[1] = lambda[1]*coef1 + coef *l[1];
	  lambda[2] = lambda[2]*coef1 + coef *l[2];
	  lambda[kk] = 0;
	}
    }
  int jj=0;
  R lmx=lambda[0];
  if (lmx<lambda[1])  jj=1, lmx=lambda[1];
  if (lmx<lambda[2])  jj=2, lmx=lambda[2];
  if(lambda[0]<0) lambda[jj] += lambda[0],lambda[0]=0;
  if(lambda[1]<0) lambda[jj] += lambda[1],lambda[1]=0;
  if(lambda[2]<0) lambda[jj] += lambda[2],lambda[2]=0;
  return kk;
};       

int  WalkInTriangle(const Maillage & Th,
		    int it, 
		    double *lambda,
                    const  KN_<R> & U,
		    const  KN_<R> & V, 
		    R & dt)
{
  const Triangle & T(Th[it]);
  int i0=T[0];
  int i1=T[1];
  int i2=T[2];
  R u   = lambda[0]*U[i0] + lambda[1]*U[i1] + lambda[2]*U[i2];
  R v   = lambda[0]*V[i0] + lambda[1]*V[i1] + lambda[2]*V[i2];
  return WalkInTriangle( Th,it,lambda,u,v,dt);
};

int Walk(const Maillage & Th,
	 int& it, 
	 R *l,
         const KN_<R>  & U,
	 const KN_<R>  & V, 
	 R dt) 
{
  int k=0;
  int j; 
  while ( (j=WalkInTriangle(Th,it,l,U,V,dt))>=0) 
    { 
      assert( l[j] == 0);
      R a= l[(j+1)%3], b= l[(j+2)%3];
      int itt = Th.ElementAdj(it,j);
      if(itt==it || itt<0)  return -1;
      it = itt;
      l[j]=0;
      l[(j+1)%3] = b;
      l[(j+2)%3] = a;
      assert(k++<1000);
    }
  return it;
};


// 
// Maillage.cpp ends here
