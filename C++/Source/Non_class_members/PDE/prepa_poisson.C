/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char prepa_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/01/31 08:49:58  e_gourgoulhon
 * Increased the number nmax of stored matrices from 100 to 200.
 *
 * Revision 1.2  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.17  2000/05/22  13:40:27  phil
 * ajout du cas dzpuis == 3
 *
 * Revision 2.16  2000/01/18  14:15:31  phil
 * enleve assert sur nobre de points min en r
 *
 * Revision 2.15  2000/01/04  19:00:21  phil
 * Double nmax
 *
 * Revision 2.14  1999/10/12  09:38:16  phil
 * passage en const Matrice &
 *
 * Revision 2.13  1999/10/11  14:29:25  phil
 * & -> &&
 *
 * Revision 2.12  1999/09/30  09:20:19  phil
 * remplacement des && en &
 *
 * Revision 2.11  1999/09/17  15:24:46  phil
 * correction definition de NMAX
 *
 * Revision 2.10  1999/06/23  12:34:44  phil
 * ajout de dzpuis = 2
 *
 * Revision 2.9  1999/04/28  10:47:12  phil
 * augmentation de NMAX a 50
 *
 * Revision 2.8  1999/04/19  14:28:47  phil
 * *** empty log message ***
 *
 * Revision 2.7  1999/04/19  14:05:32  phil
 * *** empty log message ***
 *
 * Revision 2.6  1999/04/16  13:18:37  phil
 * *** empty log message ***
 *
 * Revision 2.5  1999/04/16  13:16:24  phil
 * *** empty log message ***
 *
 * Revision 2.4  1999/04/14  13:56:42  phil
 * Sauvegarde des Matrices deja calculees
 *
 * Revision 2.3  1999/04/07  14:38:05  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/04/07  14:26:22  phil
 * *** empty log message ***
 *
 * Revision 2.1  1999/04/07  14:24:14  phil
 * Les matrices sont passees par reference
 *
 * Revision 2.0  1999/04/07  14:11:07  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//fichiers includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrice.h"
#include "tbl.h"
#include "type_parite.h"
#include "indent.h"
#include "proto.h"

/*
 * Fonctions supprimant le nombre de colonnes (les premieres)
  et de lignes (les dernieres) a l'operateur renvoye par laplacien_mat, de facon
  a ce qu'il ne soit plus degenere. Ceci doit etre fait apres les combinaisons
  lineaires. La mise a bandes et la decomposition LU sont egalement effectuees ici
  
  
  Entree : lap :  resultat de laplacien_mat
	    l : associe a lap
	    puis : puissance dans la ZEC
	    base_r : base de developpement
	    
  Sortie : renvoie un operateur non degenere ....
 */
 
 


		//------------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _prepa_nondege_pas_prevu(const Matrice &lap, int l, double echelle, int puis) {
    cout << "Construction non degeneree pas prevue..." << endl ;
    cout << "l : " << l << endl ;
    cout << "lap : " << lap << endl ;
    cout << "echelle : " << echelle << endl ;
    cout << " puis : " << puis << endl ;
    abort() ;
    exit(-1) ;
    Matrice res(1, 1) ;
    return res;
}



	     	//-------------------
	       //--  R_CHEB   -------
	      //--------------------

Matrice _prepa_nondege_r_cheb (const Matrice &lap, int l, double echelle, int) {
    
    
    int n = lap.get_dim(0) ;
    
   const int nmax = 200 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
   static double vieux_echelle = 0;
   
   // Si on a change l'echelle : on detruit tout :
   if (vieux_echelle != echelle) {
       for (int i=0 ; i<nb_dejafait ; i++) {
	   l_dejafait[i] = -1 ;
	   nr_dejafait[i] = -1 ;
	   delete tab[i] ;	 
       }
	vieux_echelle = echelle ;
	 nb_dejafait = 0 ;
   }
      
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_prepa_nondege_r_cheb : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
    //assert (l<n) ;
    
    Matrice res(n-2, n-2) ;
    res.set_etat_qcq() ;
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n-2 ; j++)
	    res.set(i, j) = lap(i, j+2) ;
	
    res.set_band(2, 2) ;
    res.set_lu() ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;  
}




	     	//------------------
	       //--  R_CHEBP   ----
	      //------------------
	      
Matrice _prepa_nondege_r_chebp (const Matrice &lap, int l, double, int) {
    
    int n = lap.get_dim(0) ;
    
       
   const int nmax = 200 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
    
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_prepa_nondege_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    assert (div(l, 2).rem == 0) ;
  //  assert (l<=2*n-2) ;
    
    if (l==0) {
	Matrice res(n-1, n-1) ;
	res.set_etat_qcq() ;
	for (int i=0 ; i<n-1 ; i++)
	    for (int j=0 ; j<n-1 ; j++)
		res.set(i, j) = lap(i, j+1) ;
	res.set_band(3, 0) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
    }
    else {
	Matrice res(n-2, n-2) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-2 ; i++)
	    for (int j=0 ; j<n-2 ; j++)
		res.set(i, j) = lap(i, j+2) ;
	
	res.set_band(2, 1) ;	
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
     }
    }
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}




	     	//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
	      
Matrice _prepa_nondege_r_chebi (const Matrice &lap, int l, double, int) {
    
    int n = lap.get_dim(0) ;
       
   const int nmax = 200 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
    
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_prepa_nondege_r_chebi : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
    assert (div(l, 2).rem == 1) ;
  //  assert (l<=2*n-1) ;
    
    if (l==1) {
	Matrice res(n-1, n-1) ;
	res.set_etat_qcq() ;
	for (int i=0 ; i<n-1 ; i++)
	    for (int j=0 ; j<n-1 ; j++)
		res.set(i, j) = lap(i, j+1) ;
	res.set_band(3, 0) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
    }
    else {
	Matrice res(n-2, n-2) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-2 ; i++)
	    for (int j=0 ; j<n-2 ; j++)
		res.set(i, j) = lap(i, j+2) ;
	
	res.set_band(2, 1) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
	} 
    }
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}




	     	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
	      
	      
Matrice _prepa_nondege_r_chebu (const Matrice &lap, int l, double, int puis) {

    switch (puis) {
	case 4 :
	    return _prepa_nondege_r_chebu_quatre (lap, l) ;
	case 3 :
	    return _prepa_nondege_r_chebu_trois (lap, l) ;
	case 2 :
	    return _prepa_nondege_r_chebu_deux (lap, l) ;
	default :
	    abort() ;
	    exit(-1) ;
	    return Matrice(0, 0) ;
    }
}

// Cas dzpuis = 4 ;
Matrice _prepa_nondege_r_chebu_quatre (const Matrice &lap, int l) {
    
    int n = lap.get_dim(0) ;
    
   const int nmax = 200; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
    
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_prepa_nondege_r_chebu : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
  //  assert (l<=n-2) ;
    
    if (l==0) {
	Matrice res(n-2, n-2) ;
	res.set_etat_qcq() ;
	for (int i=0 ; i<n-2 ; i++)
	    for (int j=0 ; j<n-2 ; j++)
		res.set(i, j) = lap(i, j+2) ;
	res.set_band(3, 0) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
    }
    else {
	Matrice res(n-3, n-3) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-3 ; i++)
	    for (int j=0 ; j<n-3 ; j++)
		res.set(i, j) = lap(i, j+3) ;
	
	res.set_band(2, 1) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
	} 
    }
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}
// Cas dzpuis = 3 ;
Matrice _prepa_nondege_r_chebu_trois (const Matrice &lap, int l) {
    
    int n = lap.get_dim(0) ;
    
   const int nmax = 200; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
    
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_prepa_nondege_r_chebu_trois : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice res(n-2, n-2) ;
    res.set_etat_qcq() ;
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n-2 ; j++)
	    res.set(i, j) = lap(i, j+2) ;
    res.set_band(2, 1) ;
    res.set_lu() ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    }
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

// Cas dzpuis = 2 ;
Matrice _prepa_nondege_r_chebu_deux (const Matrice &lap, int l) {
    
   int n = lap.get_dim(0) ;
    
   const int nmax = 200; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
    
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	indice = conte ;
    
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_prepa_nondege_r_chebu : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
  //  assert (l<=n-2) ;
    
    if (l==0) {
	Matrice res(n-2, n-2) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-2 ; i++)
	    for (int j=0 ; j<n-2 ; j++)
		res.set(i, j) = lap(i, j+2) ;
	res.set_band(3, 2) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
    }
    else {
	Matrice res(n-1, n-1) ;
	res.set_etat_qcq() ;
	for (int i=0 ;i<n-1 ; i++)
	    for (int j=0 ; j<n-1 ; j++)
		res.set(i, j) = lap(i, j+1) ;
	res.set_band(4, 1) ;
	res.set_lu() ;
	tab[nb_dejafait] = new Matrice(res) ;
	nb_dejafait ++ ;
	return res ;
	}
    }
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

	     	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      

Matrice prepa_nondege(const Matrice &lap, int l, double echelle, int puis, int base_r)
{

		// Routines de derivation
    static Matrice (*prepa_nondege[MAX_BASE])(const Matrice&, int, double, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    prepa_nondege[i] = _prepa_nondege_pas_prevu ;
	}
		// Les routines existantes
	prepa_nondege[R_CHEB >> TRA_R] = _prepa_nondege_r_cheb ;
	prepa_nondege[R_CHEBU >> TRA_R] = _prepa_nondege_r_chebu ;
	prepa_nondege[R_CHEBP >> TRA_R] = _prepa_nondege_r_chebp ;
	prepa_nondege[R_CHEBI >> TRA_R] = _prepa_nondege_r_chebi ;
    }
    
    Matrice res(prepa_nondege[base_r](lap, l, echelle, puis)) ;
    return res ;
}

