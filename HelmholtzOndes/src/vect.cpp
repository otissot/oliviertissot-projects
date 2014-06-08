#include "vect.hpp"

vect vect::operator+(vect v){
	vect res;
	for(unsigned i=0; i<size(); i++)
		res.push_back(v[i]+(*this)[i]);
	return res;
}

vect vect::operator-(vect v){
	return (*this)+v*(-1);
}

void vect::operator+=(vect v){
	(*this)=(*this)+v;
}

void vect::operator-=(vect v){
	(*this)=(*this)-v;
}

vect vect::operator*(const double& a){
	vect res;
	for(unsigned i=0; i<size(); i++)
		res.push_back(a*(*this)[i]);
	return res;
}

double vect::operator*(vect& a){
	double res=0;
	for(unsigned i=0; i<size(); i++)
		res+=a[i]*(*this)[i];
	return res;
}

double vect::norme_2c(){
	double res=0;
	for(unsigned i=0; i<size();i++)
		res+=(*this)[i]*(*this)[i];
	return res;
}

ostream& operator<<(ostream& o, const vect& v){
	for(unsigned i=0;i<v.size();i++)
		o<<v[i]<<endl;
	return o;
}

vect operator*(const double& a, vect v){
	return v*a;
}
