// Param.hpp --- 
// 
// Filename: Param.hpp
// Description: 
// Author: O.TISSOT 
// Created: ven. déc. 20 20:10:33 2013 (+0100)
// Version: 
// Last-Updated: lun. févr. 10 10:07:53 2014 (+0100)
//           By: Olivier
//     Update #: 43
// Compatibility: 
// 
// 

// Commentary: 
// 
// Une classe plus ou moins generique pour rassembler les parametres utilisateur du calcul.
// On la met static pour eviter d'avoir a l'instancier.

// Code:

#ifndef PARAM_HPP
#define PARAM_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

class Param{
public:
  // Attributs
  double Tmax;             // La borne en temps
  int TSnapShot;           // Pour definir les instants que l'on sauvegarde
  string Mesh;             // Le nom du fichier de maillage
  double Viscosite;        // Valeur de la viscosite
  double dt;               // Le pas de temps 
  // Constructeur :
  Param():Tmax(0.),TSnapShot(0),Mesh(),Viscosite(0.),dt(0.){};
  
  // Initialisation via un fichier.
  // Le format du fichier (l'ordre des donnees n'est pas important) :
  // Mesh
  // Path_to_Mesh/Mesh_filename
  // Tmax
  // Tmax_value
  // TSnapShot
  // TSnapShot_value
  // dt
  // dt_value
  // Viscosite
  // Viscosite_value
  // ...
  //
  // Si on entre en nom de fichier (de parametres) "terminal" (attention a la casse !) alors
  // la construction se fera via le terminal.
  void init(const string&);

};

#endif

// 
// Param.hpp ends here
