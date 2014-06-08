#include "assemblage.hpp"

void build_mat(matrice& a, matrice& A, const vector< point >& p, const vector< edge >& e,  const vector< triangle >& t, int argc, char** argv){
	a.nulle(p.size());
	A.nulle(p.size());
	int sij=0;
  	int sik=0;
  	double cont=0;
  	double d=0;
  	vect temp;
  	vect ds21; 
  	vect ds31;					//accroissement des coordonnees sur (s2,s1) et (s3,s1)
  	temp.resize(2);
  	ds21.resize(2);
  	ds31.resize(2);
  	vector< vect > b(2);				//matrice liee au changement de coordonnees
 	b[0].resize(2);
  	b[1].resize(2);
  	vector< vect > grad(3);				//les gradients elementaires
  	grad[0].push_back(-1);
  	grad[0].push_back(-1);
  	grad[1].push_back(1);
  	grad[1].push_back(0);
  	grad[2].push_back(0);
  	grad[2].push_back(1);
	double mass_elem[3][3];
	double bord[2][2];
	if(argc==3)	
		cout<<"Construction pour equation de Poisson ..."<<endl;
	if(argc==4){
		cout<<"Construction pour equation d'Helmotz de parametre "<<argv[3]<<" ..."<<endl;
		float k = atof(argv[3]);
		cout<<"k : "<<k<<endl;
		//La masse elementaire :
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				if (j==i)
					mass_elem[i][j] = -k*k* (1/6.);
				else
					mass_elem[i][j] = -k*k* (1/12.);
			}
		}
		//Terme de bord elementaire :
		bord[0][0] = 1/3.;
		bord[0][1] = 1/6.;
		bord[1][0] = 1/6.;
		bord[1][1] = 1/3.;
		
	}

	if(argc==5){
		cout<<"Construction pour equation des ondes..."<<endl;
		//La masse elementaire :
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				if (j==i)
					mass_elem[i][j] = 1/6.;
				else
					mass_elem[i][j] = 1/12.;
			}
		}
		//Terme de bord elementaire :
		bord[0][0] = 1/3.;
		bord[0][1] = 1/6.;
		bord[1][0] = 1/6.;
		bord[1][1] = 1/3.;
		
	}

  	for(unsigned i=0;i<t.size();i++){

		ds21[0] = p[t[i].tri[1]].coord[0]-p[t[i].tri[0]].coord[0];
		ds21[1] = p[t[i].tri[1]].coord[1]-p[t[i].tri[0]].coord[1];
		ds31[0] = p[t[i].tri[2]].coord[0]-p[t[i].tri[0]].coord[0];
		ds31[1] = p[t[i].tri[2]].coord[1]-p[t[i].tri[0]].coord[1];

		b[0][0] = ds31*ds31;
		b[0][1] = ds31*ds21*(-1);
		b[1][0] = ds31*ds21*(-1);
		b[1][1] = ds21*ds21;

		d=4*aire(t[i],p);

		//Double boucle sur les sommets :
		for(int j=0;j<3;j++){
			sij=t[i].tri[j];
			temp[0] = b[0]*grad[j];
			temp[1] = b[1]*grad[j];
			for(int k=0;k<3;k++){
				sik = t[i].tri[k];
				cont = grad[k]*temp;
				cont /= d;
				//Si Helmotz :
				cont += aire(t[i],p)*mass_elem[j][k];
				if(argc == 4)
					a.set_val(sik,sij,a.valeur(sik,sij)+cont);

				//Rigidite si equation des ondes :
				if(argc==5){
					A.set_val(sik,sij,A.valeur(sik,sij)+cont);
					if(*argv[4] == 'o')
						//Masse avec mass-lumping :
						a.set_val(sik,sij,a.valeur(sik,sij)+((sij == sik) * aire(t[i],p) / 3.));
					if(*argv[4] == 'n')
						//Masse sans mass-lumping :
						a.set_val(sik,sij,a.valeur(sik,sij)+aire(t[i],p)*mass_elem[j][k]);
				}
			}
		}
  	}
	//On ajoute le terme de bord si Helmotz car condition de Robin :
	double h = 0;
	if(argc==4){
		for(unsigned i=0;i<e.size();i++){
			for(int j=0;j<2;j++){
				sij = e[i].ext[j];
				h = longueur(e[i],p);
				for(int k=0;k<2;k++){
					sik=e[i].ext[k];
					a.set_val(sik,sij,a.valeur(sik,sij)+h*bord[j][k]);
				}
			}
		}	
	}

	if(argc==5){
		for(unsigned i=0;i<e.size();i++){
			for(int j=0;j<2;j++){
				sij = e[i].ext[j];
				h = longueur(e[i],p);
				for(int k=0;k<2;k++){
					sik=e[i].ext[k];
					A.set_val(sik,sij,A.valeur(sik,sij)+h*bord[j][k]);
				}
			}
		}	
	}
}

void pseudo_elim(matrice& A, const vector< edge >& e){
   	int l=0;
	int si, sj;
  	for(unsigned i=0;i<e.size();i++){
		si = e[i].ext[0];
		for(unsigned j=0;j<e.size();j++){
			sj = e[j].ext[0];
			A.set_val(si,sj,(sj==si));
			l++;
		}
  	}
	cout<<"Nombre d'iteration lors de la pseudo-elimination : "<<l<<endl<<endl;
}

double g(double X, double Y, double k, vect& w){
	double res;
	vect n;
	n.resize(2);
	if(X==0 && Y>0 && Y<1){
		n[0] = -1;
		n[1] = 0;
	}
	if(X==1 && Y>0 && Y<1){
		n[0] = 1;
		n[1] = 0;     
	}
	if(Y==0){
		n[0] = 0;
		n[1] = -1;
	}
	if(Y==1){
		n[0] = 0;
		n[1] = 1;
	}
	res = sin(k*X*w[0]+k*Y*w[1]) + k*cos(k*X*w[0]+k*Y*w[1])*w*n;
	return res;
}

void build_f(vect& f, const vector< edge >& e, const vector< point >& p, double k, int flag){
  	f.resize(p.size(),0);
	double buf = 5;
	int index_mil = 0;
	vect w;
	switch(flag){
		//Source au centre avec Dirichlet homogene a 0 au bord :
		case 0   : for(unsigned i=0;i<p.size();i++){
				   if(buf != min(buf, pow(p[i].coord[0]-0.5,2)+pow(p[i].coord[1]-0.5,2)))
					   index_mil=i;
				   buf = min(buf, pow(p[i].coord[0]-0.5,2)+pow(p[i].coord[1]-0.5,2));
			   }
			   f[index_mil] = 1;
		   	   break;
		//Second membre projet :
		case 1   : w.resize(2, 1./1.41421356);
			   for(unsigned i=0;i<e.size();i++)
				   for(int j=0;j<2;j++)
				 	   f[e[i].ext[j]] += (1/2.)* longueur(e[i],p) * g(p[e[i].ext[1-j]].coord[0], p[e[i].ext[1-j]].coord[1], k, w);
			   break;
		default  : cout<<"Construction second membre non specifiee ..."<<endl;
		   	   exit(1);
		   	   break;
  	}
	
}




















