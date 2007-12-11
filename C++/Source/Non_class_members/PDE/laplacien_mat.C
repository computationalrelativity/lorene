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


char laplacien_mat_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2007/12/11 15:28:22  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.8  2007/06/06 07:43:28  p_grandclement
 * nmax increased to 200
 *
 * Revision 1.7  2005/01/27 10:19:43  j_novak
 * Now using Diff operators to build the matrices.
 *
 * Revision 1.6  2004/10/05 15:44:21  j_novak
 * Minor speed enhancements.
 *
 * Revision 1.5  2004/02/06 10:53:54  j_novak
 * New dzpuis = 5 -> dzpuis = 3 case (not ready yet).
 *
 * Revision 1.4  2003/01/31 08:49:58  e_gourgoulhon
 * Increased the number nmax of stored matrices from 100 to 200.
 *
 * Revision 1.3  2002/10/16 14:37:11  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/10/09 12:47:32  j_novak
 * Execution speed improved
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.15  2000/05/22  13:36:33  phil
 * ajout du cas dzpuis == 3
 *
 * Revision 2.14  2000/01/04  19:00:09  phil
 * Double nmax
 *
 * Revision 2.13  1999/10/11  14:27:26  phil
 * & -> &&
 *
 * Revision 2.12  1999/09/30  09:17:11  phil
 * remplacement && en & et initialisation des variables statiques.
 *
 * Revision 2.11  1999/09/17  15:19:30  phil
 * correction de definition de nmax
 *
 * Revision 2.10  1999/09/03  14:07:15  phil
 * pas de modif
 *
 * Revision 2.9  1999/07/08  09:54:20  phil
 * *** empty log message ***
 *
 * Revision 2.8  1999/07/07  10:02:39  phil
 * Passage aux operateurs 1d
 *
 * Revision 2.7  1999/06/23  12:34:07  phil
 * ajout de dzpuis = 2
 *
 * Revision 2.6  1999/04/28  10:45:54  phil
 * augmentation de NMAX a 50
 *
 * Revision 2.5  1999/04/19  14:03:42  phil
 * *** empty log message ***
 *
 * Revision 2.4  1999/04/16  13:15:52  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/04/14  13:57:26  phil
 * Sauvegarde des Matrices deja calculees
 *
 * Revision 2.2  1999/04/13  13:58:30  phil
 * ajout proto.h
 *
 * Revision 2.1  1999/04/07  14:22:17  phil
 * *** empty log message ***
 *
 * Revision 2.0  1999/04/07  14:09:41  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//fichiers includes
#include <stdio.h>
#include <stdlib.h>

#include "diff.h"
#include "proto.h"

/*
 * Routine calculant l'operateur suivant :
 * 
 * -R_CHEB : r^2d2sdx2+2*rdsdx-l*(l+1)Id
 * 
 * -R_CHEBP et R_CHEBI : d2sdx2+2/r dsdx-l(l+1)/r^2
 * 	
 * -R_CHEBU : d2sdx2-l(l+1)/x^2
 * 
 * -R_JACO02 : d2sdx2 + 2/r dsdx - l(l+1)/r^2
 * 
 * Entree :
 *	-n nbre de points en r
	-l voire operateur.
	-echelle utile uniquement pour R_CHEB : represente beta/alpha 
						(cf mapping)
	
	- puis : exposant de multiplication dans la ZEC
	- base_r : base de developpement
	
    Sortie :
	La fonction renvoie la matrice.
	
 */
		//-----------------------------------
		// Routine pour les cas non prevus --
		//-----------------------------------

Matrice _laplacien_mat_pas_prevu(int n, int l, double echelle, int puis) {
    cout << "laplacien pas prevu..." << endl ;
    cout << "n : " << n << endl ;
    cout << "l : " << l << endl ;
    cout << "puissance : " << puis << endl ;
    cout << "echelle : " << echelle << endl ;
    abort() ;
    exit(-1) ;
    Matrice res(1, 1) ;
    return res;
}


		   //-------------------------
		   //--   CAS R_JACO02   -----
		   //-------------------------
		    

Matrice _laplacien_mat_r_jaco02 (int n, int l, double, int) {
   
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_mat_r_jaco02 : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       

    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Diff_dsdx2 d2(R_JACO02, n) ;
    Diff_sxdsdx sxd(R_JACO02, n) ;
    Diff_sx2 sx2(R_JACO02, n) ;
    
    tab[nb_dejafait] = new Matrice(d2 + 2.*sxd -(l*(l+1))*sx2) ;

    indice = nb_dejafait ;
    nb_dejafait ++ ;	
   }
    
   return *tab[indice] ;
}


		   //-------------------------
		   //--   CAS R_CHEBP    -----
		   //--------------------------
		    

Matrice _laplacien_mat_r_chebp (int n, int l, double, int) {
   
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_mat_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       

    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Diff_dsdx2 d2(R_CHEBP, n) ;
    Diff_sxdsdx sxd(R_CHEBP, n) ;
    Diff_sx2 sx2(R_CHEBP, n) ;
    
    tab[nb_dejafait] = new Matrice(d2 + 2.*sxd -(l*(l+1))*sx2) ;

    indice = nb_dejafait ;
    nb_dejafait ++ ;	
   }
    
   return *tab[indice] ;
}



		   //------------------------
		   //--   CAS R_CHEBI    ----
		   //------------------------
		    

Matrice _laplacien_mat_r_chebi (int n, int l, double, int) {
   
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_mat_r_chebi : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Diff_dsdx2 d2(R_CHEBI, n) ;
    Diff_sxdsdx sxd(R_CHEBI, n) ;
    Diff_sx2 sx2(R_CHEBI, n) ;
    
    tab[nb_dejafait] = new Matrice(d2 + 2.*sxd - (l*(l+1))*sx2) ;
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   } 
    
   return *tab[indice] ;
}




		   //-------------------------
		   //--   CAS R_CHEBU    -----
		   //-------------------------

Matrice _laplacien_mat_r_chebu( int n, int l, double, int puis) {
    Matrice res(n, n) ;
    res.set_etat_qcq() ;
    switch (puis) {
	case 4 :
	    res = _laplacien_mat_r_chebu_quatre (n, l) ;
	    break ;
	case 3 :
	    res = _laplacien_mat_r_chebu_trois (n, l) ;
	    break ;
	case 2 :
	    res = _laplacien_mat_r_chebu_deux (n, l) ;
	    break ;
	case 5 :
	    res = _laplacien_mat_r_chebu_cinq (n, l) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
    }
    return res ;
}
    
    // Cas ou dzpuis = 4
Matrice _laplacien_mat_r_chebu_quatre (int n, int l) {
        
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_mat_r_chebu_quatre : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Diff_dsdx2 dd(R_CHEBU, n) ;
    Diff_sx2 xx(R_CHEBU, n) ;
    
    tab[nb_dejafait] = new Matrice(dd-(l*(l+1))*xx) ;
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   } 
    
   return *tab[indice] ;
}

// Cas ou dzpuis =3 
Matrice _laplacien_mat_r_chebu_trois (int n, int l) {
        
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_mat_r_chebu_trois : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Diff_xdsdx2 xd2(R_CHEBU, n) ;
    Diff_sx sx(R_CHEBU, n) ;

    tab[nb_dejafait] = new Matrice(xd2 -(l*(l+1))*sx) ;
    
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   } 
    
   return *tab[indice] ;
}


    //Cas ou dzpuis = 2
Matrice _laplacien_mat_r_chebu_deux (int n, int l) {
        
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_mat_r_chebu_deux : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;

    Diff_x2dsdx2 x2dd(R_CHEBU, n) ;
    Diff_id id(R_CHEBU, n) ;
    
    tab[nb_dejafait] = new Matrice(x2dd - (l*(l+1))*id) ;
    
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   } 

   return *tab[indice] ;
}

    //Cas ou dzpuis = 5
Matrice _laplacien_mat_r_chebu_cinq (int n, int l) {
        
   const int nmax = 200 ;// Nombre de Matrices stockees
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
	   cout << "_laplacien_mat_r_chebu_cinq : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Diff_x2dsdx2 x2dd(R_CHEBU, n) ;
    Diff_xdsdx xd1(R_CHEBU, n) ;
    Diff_id id(R_CHEBU, n) ;

    tab[nb_dejafait] = new Matrice( x2dd + 6.*xd1 + (6-l*(l+1))*id ) ;

    indice = nb_dejafait ;
    nb_dejafait ++ ;
   } 
    
   return *tab[indice] ;
}

		   //-------------------------
		   //--   CAS R_CHEB    -----
		   //-----------------------
		    

Matrice _laplacien_mat_r_cheb (int n, int l, double echelle, int) {
            
   const int nmax = 200 ;// Nombre de Matrices stockees
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
       
       nb_dejafait = 0 ;
       vieux_echelle = echelle ;
   }
   
   int indice = -1 ;
   
   // On determine si la matrice a deja ete calculee :
   for (int conte=0 ; conte<nb_dejafait ; conte ++)
       if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
	   indice = conte ;
   
   // Calcul a faire : 
   if (indice  == -1) {
       if (nb_dejafait >= nmax) {
	   cout << "_laplacien_mat_r_cheb : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
       
       l_dejafait[nb_dejafait] = l ;
       nr_dejafait[nb_dejafait] = n ;
       
       Diff_dsdx2 d2(R_CHEB, n) ;
       Diff_xdsdx2 xd2(R_CHEB, n) ;
       Diff_x2dsdx2 x2d2(R_CHEB, n) ;
       Diff_dsdx d1(R_CHEB, n) ;
       Diff_xdsdx xd1(R_CHEB, n) ;
       Diff_id id(R_CHEB, n) ;

       tab[nb_dejafait] = 
	   new Matrice(x2d2 + (2*echelle)*xd2 + (echelle*echelle)*d2
		       + 2*xd1 + (2*echelle)*d1 - (l*(l+1))*id) ;
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   } 
   
   return *tab[indice] ;  
}


		 //--------------------------
		//- La routine a appeler  ---
	       //----------------------------
Matrice laplacien_mat(int n, int l, double echelle, int puis, int base_r)
{

		// Routines de derivation
    static Matrice (*laplacien_mat[MAX_BASE])(int, int, double, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    laplacien_mat[i] = _laplacien_mat_pas_prevu ;
	}
		// Les routines existantes
	laplacien_mat[R_CHEB >> TRA_R] = _laplacien_mat_r_cheb ;
	laplacien_mat[R_CHEBU >> TRA_R] = _laplacien_mat_r_chebu ;
	laplacien_mat[R_CHEBP >> TRA_R] = _laplacien_mat_r_chebp ;
	laplacien_mat[R_CHEBI >> TRA_R] = _laplacien_mat_r_chebi ;
	laplacien_mat[R_JACO02 >> TRA_R] = _laplacien_mat_r_jaco02 ;
    }
    
    return laplacien_mat[base_r](n, l, echelle, puis) ;
}

