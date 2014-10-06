/* 
 * Code to move the doc++ comments /// from the end of a line to 
 * the beginning of previous line, in order that doc++ processes the
 * source file correctly. The name of the output reorganized file 
 * has a suffix _r. 
 *
 */

/*
 *   Copyright (c) 1999 Jean-Alain Marck
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:05:46  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2003/10/15 08:43:32  e_gourgoulhon
 * Modif. comments.
 *
 *
 *
 * $Header$
 *
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <iostream>
using namespace std ;

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
