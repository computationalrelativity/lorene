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
#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

/*
 * Routine caluclant l'operateur suivant :
 * 
 * -R_CHEB : r^2d2sdx2+2*rdsdx-l*(l+1)Id
 * 
 * -R_CHEBP et R_CHEBI : d2sdx2+2/r dsdx-l(l+1)/r^2
 * 
 * -R_CHEBU : d2sdx2-l(l+1)/x^2
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
		   //--   CAS R_CHEBP    -----
		   //--------------------------
		    

Matrice _laplacien_mat_r_chebp (int n, int l, double, int) {
   
   const int nmax = 100 ;// Nombre de Matrices stockees
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
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

   double* vect  = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBP) ;
	
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j] ; 
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxdsdx_1d (n, &vect, R_CHEBP) ;
	    for (int j=0 ; j<n ; j++)
		xd.set(j, i) = vect[j] ;
	
	}
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEBP) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
   
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+2*xd-l*(l+1)*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}



		   //------------------------
		   //--   CAS R_CHEBI    ----
		   //------------------------
		    

Matrice _laplacien_mat_r_chebi (int n, int l, double, int) {
   
   const int nmax = 100 ;// Nombre de Matrices stockees
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
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBI) ;  // appel dans le cas impair
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxdsdx_1d (n, &vect, R_CHEBI) ;
	    for (int j=0 ; j<n ; j++)
		xd.set(j, i) = vect[j] ;
	}
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEBI) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+2*xd-l*(l+1)*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
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
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEBU) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd-l*(l+1)*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
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
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    double* auxi = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
	mult_xm1_1d_cheb (n, vect, auxi) ;
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = auxi[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxm1_1d_cheb (n, vect) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    delete [] auxi ;
    
    Matrice res(n, n) ;
    res = dd-l*(l+1)*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
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
    
    Matrice res(n, n) ;
    res.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    double* x2vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
	mult2_xm1_1d_cheb (n, vect, x2vect) ; // multiplication par (x-1)^2
	for (int j=0 ; j<n ; j++)
	    res.set(j, i) = x2vect[j] ;
    }
    
    delete [] vect ;
    delete [] x2vect ;
    
    for (int i=0 ; i<n ; i++)
	res.set(i, i) = res.set(i, i) - l*(l+1) ;
    
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

		   //-------------------------
		   //--   CAS R_CHEB    -----
		   //-----------------------
		    

Matrice _laplacien_mat_r_cheb (int n, int l, double echelle, int) {
            
   const int nmax = 100 ;// Nombre de Matrices stockees
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
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	for (int j=0 ; j<n ; j++)
	    dd.set(j, i) = vect[j]*echelle*echelle ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	multx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	    dd.set(j, i) += 2*echelle*vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
	multx_1d (n, &vect, R_CHEB) ;
	multx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	    dd.set(j, i) += vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxdsdx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<n ; j++)
	    xd.set(j, i) = vect[j]*echelle ;
	}
    
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sxdsdx_1d (n, &vect, R_CHEB) ;
	multx_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	    xd.set(j, i) += vect[j] ;
	}
	   
    for (int i=0 ; i<n ; i++) {
	for (int j=0 ; j<n ; j++)
	    vect[j] = 0 ;
	vect[i] = 1 ;
	sx2_1d (n, &vect, R_CHEB) ;
	for (int j=0 ; j<n ; j++)
	    xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+2*xd-l*(l+1)*xx ;   
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
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
    }
    
    Matrice res(laplacien_mat[base_r](n, l, echelle, puis)) ;
    return res ;
}

