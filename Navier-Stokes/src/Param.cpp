// Param.cpp --- 
// 
// Filename: Param.cpp
// Description: 
// Author: O.TISSOT 
// Created: ven. déc. 20 20:25:24 2013 (+0100)
// Version: 
// Last-Updated: lun. févr. 10 10:09:05 2014 (+0100)
//           By: Olivier
//     Update #: 25
// Compatibility: 
// 
// 

// Commentary: 
// 
// Implementation de Param.hpp
// 
// 

// Code:

#include "Param.hpp"

// Allocation de memoire
void Param::init(const string& filename){
  // Les parametres vont etre recuperes grace au terminal
  if(filename == "terminal"){
    cout << "Tmax? ";
    cin >> Tmax;
    cout << "dt? ";
    cin >> dt;
    cout << "TSnapShot? ";
    cin >> TSnapShot;
    cout << "Fichier de maillage (avec le chemin)? ";
    cin >> Mesh;
    cout << "Viscosite? ";
    cin >> Viscosite;
  }
  // Les parametres sont dans un fichier
  else{
    fstream paramfile(filename.c_str(), ifstream::in);
    // On verifie que le fichier existe bien
    if(paramfile){
      string buffer;
      // On lit les donnees : l'ordre d'apparition n'est pas important
      do{
	paramfile >> buffer;
	if(buffer == "Tmax")
	  paramfile >> Tmax;
	else if(buffer == "TSnapShot")
	  paramfile >> TSnapShot;
	else if(buffer == "Mesh")
	  paramfile >> Mesh;
	else if(buffer == "Viscosite")
	  paramfile >> Viscosite;
	else if(buffer == "dt")
	  paramfile >> dt;
	// On check si la structure du fichier est conforme sinon on jette l'utilisateur
	else{
	  cerr << "Param::Param(const string&) : Le mot cle : " << buffer << ", est inconnu..." << endl;
	  exit(1);
	}
      }while(!paramfile.eof());
      paramfile.close();
    }
    else{
      cerr << "Param::Param(const string&) : fichier de parametre non trouve." << endl;
      exit(1);
    }
  }
}
// 
// Param.cpp ends here
