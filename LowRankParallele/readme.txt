Contenu de chaque dossier :

- include : header pour pouvoir appeler les wrapers C des fonctions ScaLAPACK sans avoir de warnings à la compilations.
- tools : routines utilitaires (écriture dans un fichier, génération des matrices...)
- prototype : 
	* le programme principal : prototype.c
	* les Makefile : il faut modifier Makefile.opts en mettant les bons chemins vers les bibliothèques
	* le dossier data qui contient les fichiers d'input/output
	* répertoire build pour les fichiers *.o
