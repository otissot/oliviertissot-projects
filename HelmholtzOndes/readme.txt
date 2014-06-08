Ce document contient des informations concernant l'utilisation du code source des parties I et II du projet


Partie I : Equation de Helmholtz
''''''''''''''''''''''''''''''

Notre code contient :
_Un dossier "mesh_data" dans lequel se trouvent les triangulations

_Un dossier "results" dans lequel les resultats seront sauvegardes

_Un Makefile qui permet de compiler les differentes parties du code sans avoir a etre modifie

_Deux fichiers mat_cr.cpp et mat_cr.hpp qui sont relatifs a l'utilisation de matrices creuses.
 C'est dans ce fichier que sont définis les solveurs lineaires.

_Deux fichiers vect.cpp et vect.hpp qui permettent d'utiliser la classe vect.

_Deux fichiers maillage.cpp et maillage.hpp qui contiennent les outils necessaires a l'importation des maillages a partir des fichiers fournis.

_Deux fichiers assemblage.cpp et assemblage.hpp qui servent a assembler la matrice principale ainsi que le second membre.

_Un fichier partieI.cpp qui contient le main


Concernant l'utilisation du code,
voici les arguments a entrer lors de l'execution du programme :
_Le nom de l'executable (./PartieI)
_Le nom complet (avec chemin d'acces si necessaire) du fichier ou est stocke le maillage
_Le nom de la methode de resolution (g, j ou l)
_Le coefficient pour l'equation d'Helmotz

Ainsi pour resoudre et afficher le resultat de l'equation de Helmholtz de coeficient k=10 avec la methode du gradient conjugue et la triangulation carre4, il faut entrer :

./PartieI mesh_data/carre4.msh g 10

Une fois le programme lance :

_A                  : zoom
_E                  : dezoom
_Z                  : bouge l'image vers le haut
_Q                  : bouge l'image vers le bas
_S                  : bouge l'image vers la gauche
_D                  : bouge l'image vers la droite
_clic gauche souris : change le point de vue 3d de la scene (non conseille)
_echap              : ferme la fenetre et quitte le programme




Partie II : Equation des ondes
''''''''''''''''''''''''''''''

Notre code contient les même fichiers que pour la partie I. Mais le fichier partie_I.cpp est remplacé par le fichier partie_II.cpp.
Il y a en plus un fichier affichage.cpp qui permet de visualiser des solutions deja calculees auparavant et qui ont ete enregistrees dans le sous-repertoire results.


Pour cette partie, les arguments a entrer lors de l'execution du programme sont :
_le nom de l'executable (./PartieII)
_Le nom complet (avec chemin d'acces si necessaire) du fichier ou est stocke le maillage
_le pas dt de discretisation souhaité
_un entier qui servira a faire les snapshots de la solution calculee pour l'enregistrer dans results
_o ou n pour utiliser ou non mass lumping

Donc pour afficher le resultat de l'equation des ondes avec dt=0.0225, la triangulation cercle2, faire des snapshot de la solution tous les 10 pas et ne pas utiliser mass lumping mais le schema classique, il faut entrer :

./PartieII mesh_data/cercle2.msh 0.0225 10 n


Pour utiliser le programme affichage :
il faut rentrer en ligne de commande exactement la meme chose que lors du calcul de la solution mais en specifiant que l'executable est affichage. Si on a calcule la solution precedente, on peut afficher de nouveau les resultats en entrant :

./affichage mesh_data/cercle2.msh 0.0225 10 n

Le programme vous demandera ensuite quel moment vous souhaitez visualiser.

Une fois les programmes lances, les touches qui permettaient precedement de modifier l'affichage restent les memes.
