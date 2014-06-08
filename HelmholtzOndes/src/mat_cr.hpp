#ifndef MAT_CR_HPP
#define MAT_CR_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <climits>
#include "vect.hpp"

using namespace std;

class matrice{
private :
	vector<int> row; 						//nombre de valeurs non nulles par ligne
	vector<int> col; 						//contient les colones des valeurs non nulles
	vector<double> val; 						//contient les valeurs non nulles de la matrice

public :
	matrice(int=0);																
	matrice(vector<int>, vector<int>, vector<double>);
	matrice(vector<double>, vector<double>, vector<double>);	//prise en compte des 0
	int size(){return row.size()-1;};
	void nulle(int);			
	void identite(int);										
	vect operator*(vect&);
	matrice operator*(matrice&);		
	double valeur(int,int);
	void set_val(const int &,const int &,const double &);
	vect get_line(int);
	vect get_col(int);
	matrice facto_lu();
	vect solveur_lu(vect&, matrice&);				//on donne en parametre la matrice LU deja construite
	double norme_c(vect&);						//calcule v~Av (v~ = v transpose)
	vect solveur_gc(vect&, vect&, const double, int=0);
	vect solveur_jacobi(vect&, vect&, const double, int);
	friend ostream& operator<<(ostream&, const matrice&);

};

#endif
