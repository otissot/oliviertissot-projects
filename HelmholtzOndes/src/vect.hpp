#ifndef VECT_HPP
#define VECT_HPP

#include <vector>
#include <iostream>
using namespace std;


class vect : public vector<double>{
	public:
		using vector<double>::push_back;
		using vector<double>::begin;
		using vector<double>::end;
		using vector<double>::size;
		using vector<double>::operator[];
		vect operator+(vect);
		vect operator-(vect);
		void operator+=(vect);
		void operator-=(vect);
		vect operator*(const double&);
		double operator*(vect& a);    			//produit scalaire
		double norme_2c();				//norme 2 du vecteur AU CARRE
		friend ostream& operator<<(ostream&, const vect&);
		
};

vect operator*(const double&, vect);				//pour pouvoir ecrire le produit dans le sens "naturel"

#endif
