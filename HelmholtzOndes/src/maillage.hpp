#ifndef MAILLAGE_HPP
#define MAILLAGE_HPP

#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<set>
#include<cmath>
#include<stdlib.h>

using namespace std;

typedef struct{
	double coord[2];
	int ref;
} point;

typedef struct{
	int ext[2];
	int ref;

} edge;

typedef struct{
	int tri[3];
	int ref;
} triangle;



class maillage {
private :
	vector< point > p;
	vector< edge > e;
	vector< triangle > t;

public :
	maillage(){};
	maillage(const char*);
	maillage(int l);				//maillage cartesien carre avec l points sur le cote
	int voisinage(int, set< int > &);		 
	void build_edge();						 
	int get_nb_point(){return p.size();};
	int get_nb_edge(){return e.size();};
	int get_nb_tria(){return t.size();};
	const vector< point >& get_p(){return p;};
	const vector< edge >& get_e(){return e;};
	const vector< triangle >& get_t(){return t;};
};

double aire(const triangle&, const vector< point >&);
double longueur(const edge&, const vector< point >&);

#endif
