/* 
 * Code to move the doc++ comments /// from the end of a line to 
 * the beginning of previous line, in order that doc++ processes the
 * source file correctly. The name of the output reorganized file 
 * has a suffix _r. 
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef OBSOLETE_HEADERS

#include <iostream.h>

#else

#include <iostream>
using namespace std ;
#endif 

int main(int argc, char** argv)
{

    // Ouverture fichiers entree et sortie
    if(argc<2) {
    	fprintf(stderr,"usage: reorg_comments filename(s) \n") ;
    	abort() ;
    }

    char nom_sortie[1000] ;
    char nom_ori[1000] ;
    char ligne[300] ;

// Boucle sur tous les fichiers

for (int ific=1 ; ific < argc ; ific++) {

    strcpy( nom_ori, argv[ific]) ;
    cout << "            reorganization of the /// comments of " << nom_ori 
	 << endl ;
 
    strcpy( nom_sortie, nom_ori ) ; 
    strcat( nom_sortie, "_r" ) ; 

    FILE* inf = fopen(nom_ori,"r");
    if(!inf) {
    	fprintf(stderr,"reorg_comments: cannot open %s ! \n", nom_ori);
    	abort() ;
    }


    FILE* inf2 = fopen(nom_sortie,"w");
    if(!inf2) {
    	fprintf(stderr,"reorg_comments: cannot open %s ! \n", nom_sortie);
    	abort() ;
    }

    // Lecture ligne a ligne
    while (fgets(ligne, 300, inf)) {
    
    	// Recherche de /// dans la ligne
    	char* triple_com = strstr(ligne, "///") ;
    	if (triple_com != NULL) {   	    	    // oui un triple com
	    // est-ce un commentaire isole ?
	    char ligne_tmp[80] ;
	    sscanf(ligne,"%s",ligne_tmp) ;  	    // vire les blancs...
	    ligne_tmp[3] = 0 ;	    	    	    // fin de ligne_tmp
	    char* occurence = strstr(ligne_tmp, "///") ;
	    if (occurence != NULL) {	    	    // c'est un comm. isole
    	    	fprintf(inf2, "%s", ligne) ; 	    // la ligne entiere
	    }
	    else {  	    	    	    	    // il faut travailler
		fprintf(inf2, "%s", triple_com) ;   // Le commentaire d'abord
		triple_com[0] = 0 ;	    	    // raccourci la ligne
		fprintf(inf2, "%s\n", ligne) ; 	    // la ligne raccourcie
	    }
    	}
	else {	    	    	    	    	    // non pas de triple com
    	    fprintf(inf2, "%s", ligne) ; 	    // la ligne entiere
    	}
    }

} // Fin de boucle sur les fichiers

    return 0 ;
      
}
