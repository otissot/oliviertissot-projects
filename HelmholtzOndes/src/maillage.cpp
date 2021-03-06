#include"maillage.hpp"

maillage::maillage(const char* nom){
	ifstream f;
	f.open(nom, ifstream::in);
		int nbr_p;
		int nbr_t;
		int nbr_e;
		f>>nbr_p;
		f>>nbr_t;
		f>>nbr_e;

		p.resize(nbr_p);
		t.resize(nbr_t);
		e.resize(nbr_e);


			for(int i=0; i<nbr_p; i++){
				for(int j=0; j<2; j++){
					f>>p[i].coord[j];
				}
				f>>p[i].ref;
			}
			
			for(int i=0; i<nbr_t; i++){
				for(int j=0; j<3; j++){
					f>>t[i].tri[j];
					t[i].tri[j]--;
				}
				f>>t[i].ref;
			}

			for(int i=0; i<nbr_e; i++){
				for(int j=0; j<2; j++){
					f>>e[i].ext[j];
					e[i].ext[j]--;
				}
				f>>e[i].ref;
			}
			
		f.close();	
}
 
maillage::maillage(int l){
	double dx = 1./double(l);
  	vector<double> bord;

  	// Nombre de noeuds
  	int Nv = (l+1)*(l+1);
 	p.resize(Nv);

 	int num = 0;
  	for(int j=0; j<l+1; j++){
 		for(int k=0; k<l+1; k++){
      
     			p[num].coord[0] = j*dx;
      			p[num].coord[1] = 1-k*dx;
      			p[num].ref=0;
     			if(k==0){
				p[num].ref=1;		//On met des references non nulles aux bords du domaine
				bord.push_back(num);
			}
      			if(k==l){
				p[num].ref=3;
				bord.push_back(num);
      			}
      			if(j==0){
				p[num].ref=2;
				bord.push_back(num);
      			}
      			if(j==l){
				p[num].ref=4;
				bord.push_back(num);
      			}
      			num++;
    		} 
  	}
	e.resize(bord.size());
	for(unsigned i=0;i<bord.size()-1;i++){
		e[i].ext[0]=bord[i];
		e[i].ext[1]=bord[i+1];
		e[i].ref=p[bord[i+1]].ref;
	}
  	// Nombre de triangles
  	int Nt = 2*l*l;
  	t.resize(Nt);  
 	num = 0;
	for(int j=0; j<l; j++){
    		for(int k=0; k<l; k++){            
      
      		t[num].tri[0] = j+ k*(l+1);
      		t[num].tri[1] = j+ (k+1)*(l+1);
      		t[num].tri[2] = j+1+ k*(l+1);
      		t[num].ref=0;
      		num++;
      
      		t[num].tri[0] = j+1+ k*(l+1);
     		t[num].tri[1] = j+ (k+1)*(l+1);
      		t[num].tri[2] = j+1+ (k+1)*(l+1);
      		t[num].ref=0;
      		num++;
      
    		}
  	}
};

int maillage::voisinage(int k, set< int > & voisin ){	//k designe le numero du point	
	vector< int > v;
	int s=t.size();
	for(int i=0; i<s; i++){
		for(int j=0; j<3; j++){
			if (t[i].tri[j]==k){
				for(int l=0; l<3; l++)
					voisin.insert(t[i].tri[l]);
			}
		}
	}
	voisin.erase(voisin.find(k));
	return voisin.size();
}


void maillage::build_edge(){
	vector< edge >e2;
	set< int >w;
	edge ed;
	for(unsigned i=1; i<p.size()+1; i++){
		w.clear();
		voisinage(i,w);	
		for(set<int>::iterator it=w.begin();it!=w.end();++it){	
			ed.ext[0]=i;
			ed.ext[1]=*it;
			ed.ref=1;
			e.push_back(ed);		
		}
	}
	int l=e.size();
	e2.resize(l);
	for(int i=0; i<l; i++)
	e2[i]=e[i];

	//Suppression des doublons (ab ab) :

	for(int k=l-1; k>-1; k--){						
		for(int n=0; n<k; n++){
			if ( ((e2[k].ext[0]==e2[n].ext[0]) && (e2[k].ext[1]==e2[n].ext[1])) )
					e.erase (e.begin()+k);
		}
	}
	int p=e.size();
	e2.clear();
	e2.resize(p);
	for(int i=0; i<p; i++)
	e2[i]=e[i];
	
	//Suppression des doublons inverse (ab ba) :

	for(int h=p-1; h>-1; h--){						
		for(int n=0; n<h; n++){
			if ( ((e2[h].ext[0]==e2[n].ext[1]) && (e2[h].ext[1]==e2[n].ext[0])) )
					e.erase (e.begin()+h);
		}
	}
}

double aire(const triangle& t, const vector< point >& p){
	double res=(p[t.tri[1]].coord[0]-p[t.tri[0]].coord[0])*(p[t.tri[2]].coord[1]-p[t.tri[0]].coord[1]);
	res-=(p[t.tri[2]].coord[0]-p[t.tri[0]].coord[0])*(p[t.tri[1]].coord[1]-p[t.tri[0]].coord[1]);
	res=(1/2.)*(abs(res));
	return res;
};

double longueur(const edge& e, const vector< point >& p){
	double res = pow(p[e.ext[1]].coord[0]-p[e.ext[0]].coord[0],2);
	res += pow(p[e.ext[1]].coord[1]-p[e.ext[0]].coord[1],2);
	return sqrt(res);
};
