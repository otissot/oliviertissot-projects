#include "mat_cr.hpp"

matrice::matrice(vector<int> r, vector<int> co, vector<double> v){
	int L, M;
	
	L=r.size();
	M=co.size();		//row et col on la meme taille
	row.resize(L);
	col.resize(M);
	val.resize(M);

	for(int i=0; i<L; i++){
		row[i]=r[i];
	}

	for(int i=0; i<M; i++){
		col[i]=co[i];
		val[i]=v[i];
	}
}

matrice::matrice(vector<double> ann, vector<double> a1n, vector<double> an1){
	int N = an1.size();

	row.push_back(0);
	
	//initialisation ligne 1 :
	
	if(ann[0]!=0){
		val.push_back(ann[0]);
		col.push_back(0);
		row.push_back(1);
	}
	if(an1[0]!=0){
		val.push_back(an1[0]);
		col.push_back(1);
		if(row.size()==1)
			row.push_back(1);
		else
			row[1]++;
	}
	
	//Remplissage lignes du milieu :

	unsigned R;
	int j=0;
	for(int i =2; i<3*N-1; i=i+3){
		R=row.size();
		if(a1n[j]!=0){
			val.push_back(a1n[j]);
			col.push_back(j);
			row.push_back(row[R-1]+1);
		}
		if(ann[j+1]!=0){
			val.push_back(ann[j+1]);
			col.push_back(j+1);
			if(R==row.size())
				row.push_back(row[R-1]+1);
			else
				row[row.size()-1]++;
		}
		if(an1[j+1]!=0){
			val.push_back(an1[j+1]);
			col.push_back(j+2);
			if(R==row.size())
				row.push_back(row[R-1]+1);
			else
				row[row.size()-1]++;
		}
		j++;		
	}

	//Remplissage derniere ligne :

	R=row.size();
	if(a1n[N-1]!=0){
		val.push_back(a1n[N-1]);
		col.push_back(N-1);
		row.push_back(row[R-1]+1);
	}
	if(ann[N]!=0){
		val.push_back(ann[N]);
		col.push_back(N);
		if(R==row.size())
			row.push_back(row[R-1]+1);
		else
			row[row.size()-1]++;
	}
}

matrice::matrice(int n){
	row.resize(n+1);
}

void matrice::identite(int d){
	val.resize(d);
	row.resize(d+1);
	col.resize(d);
	row[0]=0;
	val[0]=1;
	col[0]=0;
	for(int i=1; i<d; i++){
		val[i]=1;
		col[i]=i;
		row[i]=row[i-1]+1;
	}
	row[d]=row[d-1]+1;
}

void matrice::nulle(int d){
	row.resize(d+1);
	col.clear();
	val.clear();
}

vect matrice::operator*(vect & co){
	//les matrices sont supposees carrees
	vect s;
	s.resize(co.size());
	int M=row.size();
	int R=0;

	for(int i=0; i<M-1; i++){
		R=row[i+1]-row[i];
		s[i]=0;
		 for(int j=0; j<R; j++){
			s[i]+=co[col[row[i]+j]]*val[row[i]+j];
		 }
	
	}
	return s;
}

vect matrice::get_line(int i){
	vect res;
	res.resize(row.size()-1);
	for(unsigned j=0;j<row.size()-1;j++)
		res[j]=valeur(i,j);
	return res;
}

vect matrice::get_col(int j){
	vect res;
	res.resize(row.size()-1);
	for(unsigned i=0;i<row.size()-1;i++)
		res[i]=valeur(i,j);
	return res;
}

matrice matrice::operator*(matrice & m){
	//les matrices sont supposees carrees
	matrice res;
	double buf;
	vect l;
	vect col;
	res.row.resize(row.size(),0);
	res.row[0]=0;

	for(unsigned i=0;i<row.size()-1;i++){
		res.row[i+1]=res.row[i];
		for(unsigned j=0;j<row.size()-1;j++){
			l=get_line(i); 
			col=m.get_col(j);
			buf=l*col;
			if(buf!=0){
				res.val.push_back(buf);
				res.col.push_back(j);
				res.row[i+1]++;	
			}
		}
	}
	
	return res;	
}

double matrice::valeur(int a,int b){
	double s=0;
	double R;
	R=row[a+1]-row[a];

	if (R!=0){
		for(int i=row[a]; i<row[a]+R ;i++){
		  	if (col[i]==b)
			  	s=val[i];
		}
	}

	return s;
}

void matrice::set_val(const int & i, const int & j, const double & a){
	int R = row[i+1] - row[i];		//nombre de valeurs non nulles sur la ligne i
	int index_val = row[i];
	bool find_val = false;
	bool find_index = false;
	for (int k=0;(k<R) && (find_val==false) && (find_index==false);k++){
		if ( col[row[i]+k] == j){ 
			find_val = true;
			if(a!=0)
				val[row[i]+k] = a;
			else {
				val.erase(val.begin()+row[i]+k);
				col.erase(col.begin()+row[i]+k);
				for (unsigned k=i;k<row.size()-1;k++)
					row[k+1]--;
			}
		}
		else {
			if (j>col[row[i]+k]){
				index_val = row[i] + k + 1;
			}
			else{
				index_val = row[i] + k;
				find_index=true;
			}
		}
	}
		
	if (!find_val){
		for (unsigned k=i;k<row.size()-1;k++)
			row[k+1]++;
		col.insert (col.begin()+index_val, j);
		val.insert (val.begin()+index_val, a);
	}
}

matrice matrice::facto_lu(){
	int n = size();
	matrice LU(n);
	double s;
	LU.set_val(0,0,valeur(0,0));		//initialisation des premieres valeurs sinon probleme dans la boucle k quand j ou i =0
	for(int j=1; j<n; j++){
		LU.set_val(0,j,valeur(0,j));
		LU.set_val(j,0,(valeur(j,0)/LU.valeur(0,0)));
	}
	
	for(int i=1; i<n; i++){
		for(int j=1; j<n; j++){
			s=valeur(i,j);
			for(int k=0; k<min(i,j); k++)
					s-=(LU.valeur(i,k)*LU.valeur(k,j));
			if(j<i){
				double c=s/LU.valeur(j,j);
				LU.set_val(i,j,c);
			}
			else LU.set_val(i,j,s);
		}
	}
	return LU;
}

vect matrice::solveur_lu(vect & b, matrice & LU){
	int n=row.size()-1;
	vect x;

	//On veut resoudre a*x=b <=> LU*x=b
	//-Resolution de LY=b avec Y=U*x :

	vector< double > y;
	y.resize(n);
	y[0]=b[0];
	double res=0;
	for(int i=1; i<n; i++){
		for(int j=0; j<i; j++){
			res+=(LU.valeur(i,j))*y[j];
		}
	y[i]=b[i]-res;
	res=0;				
	}

	//-Resolution de U*x=Y :

	x.resize(n);

	x[n-1]=y[n-1]/LU.valeur(n-1,n-1);
	for(int i=n-2; i>-1; i--){
		for(int j=i+1; j<n; j++){
			res+=LU.valeur(i,j)*x[j];
		}
	x[i]=(1./LU.valeur(i,i))*(y[i]-res);
	res=0;
	}

	return x;
}

double matrice::norme_c(vect & u){
	double res;
	vect v = (*this)*u; 
	res=u*v;

	return res;
}

vect matrice::solveur_gc(vect& b, vect& init, const double eps, int maxIteration){
	vect res;
	vect r;
	vect p;
	double alp,bet;
	vect buf;		//pour pouvoir calculer b
	
	res = init;
	r = (*this) * init - b;
	p=r;
	int i=0;
	if(maxIteration == 0)
		maxIteration = INT_MAX;
	do{
		alp = r.norme_2c()/norme_c(p);
		res -= alp*p;
		buf = r;
		r -= alp*((*this)*p);
		bet = r.norme_2c()/buf.norme_2c();
		p = r + bet*p;
		alp=sqrt(r.norme_2c());
		i++;
		if(sqrt(r.norme_2c())<eps)
			i=maxIteration;
	}while(i<maxIteration);

	//cout<<"Nombre iterration : "<<i<<endl;
	
	return res;

}

vect matrice::solveur_jacobi(vect& b, vect& init, const double eps, int maxIteration){
	vect res = init;
	int n = row.size()-1;
	matrice d(n), ef(n);
	vect ann, buf;
	ann.resize(n);
	buf.resize(n);

	//Construction des 2 matrices D et E+F :

	for(int i=0; i<n; i++){
		ann[i]=1/valeur(i,i);  //en fait on cree directement D-1 et on utilise un vect pour gagner de la place
		d.set_val(i,i,1/valeur(i,i));
	}
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i!=j) ef.set_val(i,j,-(valeur(i,j)));
		}
	}
	cout<<"Produit en cours ..."<<endl;
	//Pour optimiser un peu le calcul :

	ef=d*ef;
	buf=d*b;

	//Construction de la solution approchee :
	int i=0;
	do{
		res = ef*res + buf;
		ann = (*this)*res - b;								//on utilise ann pour gagner de la memoire
		i++;
	}while( (sqrt(ann.norme_2c())/sqrt(b.norme_2c()) > eps) && (i<maxIteration) );		//test d'arret du cours
	
	cout<<"Nombre iterration : "<<i<<endl;
	cout<<"Erreur : "<<(sqrt(ann.norme_2c())/sqrt(b.norme_2c()))<<endl;
	cout<<"Eps : "<<eps<<endl;
	return res;
}

ostream& operator<<(ostream& o, const matrice& co){
	
	int M=co.row.size();			//taille du vecteur row de c
	int R=0;				//nombre de valeurs de la ligne i

	//on affiche pour toutes les valeurs non nulles pour chaque ligne (c, valeur non nulle)
	//on passe a la ligne une fois la ligne de la matrice terminee
	for(int i=0; i<M-1; i++){
		R=co.row[i+1]-co.row[i];
		o<<"ligne "<<i<<"\t";
			for(int j=0; j<R; j++){
				o<<"("<<co.col[co.row[i]+j]<<";"<<co.val[co.row[i]+j]<<")\t";
			}
		o<<endl;
	}
	return o;
}
