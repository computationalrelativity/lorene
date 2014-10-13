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


char comb_lin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2014/10/13 08:53:28  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/10/06 15:16:08  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2008/02/18 13:53:42  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.7  2007/12/12 12:30:48  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.6  2007/06/06 07:43:28  p_grandclement
 * nmax increased to 200
 *
 * Revision 1.5  2004/02/09 09:33:56  j_novak
 * Minor modif.
 *
 * Revision 1.4  2004/02/06 10:53:54  j_novak
 * New dzpuis = 5 -> dzpuis = 3 case (not ready yet).
 *
 * Revision 1.3  2002/10/16 14:37:11  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/10/09 12:47:31  j_novak
 * Execution speed improved
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.15  2000/09/07  12:52:04  phil
 * *** empty log message ***
 *
 * Revision 2.14  2000/05/22  13:39:01  phil
 * ajout du cas dzpuis == 3
 *
 * Revision 2.13  2000/01/04  18:59:58  phil
 * Double nmax
 *
 * Revision 2.12  1999/10/12  09:38:52  phil
 * passage en const Matrice&
 *
 * Revision 2.11  1999/10/11  14:26:07  phil
 * & _> &&
 *
 * Revision 2.10  1999/09/30  09:25:36  phil
 * ajout condition sur dirac=0
 *
 * Revision 2.9  1999/09/30  09:21:51  phil
 * remplacement des && en &
 *
 * Revision 2.8  1999/09/17  15:22:47  phil
 * correction definition NMAX
 *
 * Revision 2.7  1999/06/23  12:34:29  phil
 * ajout de dzpuis = 2
 *
 * Revision 2.6  1999/04/28  10:48:19  phil
 * augmentation de NMAX a 50
 *
 * Revision 2.5  1999/04/19  14:01:45  phil
 * *** empty log message ***
 *
 * Revision 2.4  1999/04/16  13:15:45  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/04/14  13:56:17  phil
 * Sauvegarde des Matrices deja calculees
 *
 * Revision 2.2  1999/04/07  14:37:34  phil
 * *** empty log message ***
 *
 * Revision 2.1  1999/04/07  14:29:51  phil
 * passage par reference
 *
 * Revision 2.0  1999/04/07  14:09:11  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//fichiers includes
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

/*   FONCTIONS FAISANT DES COMBINAISONS LINEAIRES DES LIGNES.
 * Existe en deux versions Tbl ou Matrice.
 * 
 * Entree :
 *  La Matrice ou le Tbl ( a une dimension)
 *  l : nbre quantique
 *  puis = puissance dans la ZEC
 *  La base de developpement en R
 * 
 * Sortie :
 *  Renvoie la matrice ou le tableau apres Combinaison lineaire
 * 
 */

// Version Matrice --> Matrice
namespace Lorene {
Matrice _cl_pas_prevu (const Matrice &source, int l, double echelle, int puis) {
    cout << "Combinaison lineaire pas prevu..." << endl ;
    cout << "Source : " << source << endl ;
    cout << "l : " << l << endl ;
    cout << "dzpuis : " << puis << endl ;
    cout << "Echelle : " << echelle << endl ;
    abort() ;
    exit(-1) ;
    return source;
}


		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Matrice _cl_r_cheb (const Matrice &source, int l, double echelle, int) {
    int n = source.get_dim(0) ;assert (n == source.get_dim(1)) ;
    
                
   const int nmax = 200 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
   static double vieux_echelle = 0 ;
   
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
	   cout << "_cl_r_cheb : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=i ; j<(n>(i+7)? i+7 : n) ; j++)
	    barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j))
				/(i+1) ;
	if (i==0) dirac = 0 ;
    }
    
    Matrice res(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=i ; j<(n>(i+5)? i+5 : n) ; j++)
	    res.set(i, j) = barre(i, j)-barre(i+2, j) ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;  
}


		//-------------------
	       //--  R_JACO02 ------
	      //-------------------

Matrice _cl_r_jaco02 (const Matrice &source, int l, double echelle, int) {
    int n = source.get_dim(0) ;assert (n == source.get_dim(1)) ;
    
                
   const int nmax = 200 ; // Nombre de Matrices stockees
   static Matrice* tab[nmax] ;  // les matrices calculees
   static int nb_dejafait = 0 ; // nbre de matrices calculees
   static int l_dejafait[nmax] ;
   static int nr_dejafait[nmax] ;
   static double vieux_echelle = 0 ;
   
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
	   cout << "_cl_r_jaco02 : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
    for (int i=0 ; i<n ; i++) {
	for (int j=i ; j<n ; j++)
	    barre.set(i, j) = source(i, j) ;
    }
    
    Matrice res(barre) ;
    for (int i=0 ; i<n ; i++)
	for (int j=i ; j<n ; j++)
	    res.set(i, j) = barre(i, j);
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;  
}


		//-------------------
	       //--  R_CHEBP   -----
	      //-------------------


Matrice _cl_r_chebp (const Matrice &source, int l, double, int) {
    
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
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
	   cout << "_cl_r_chebp : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	    barre.set(i, j) = (1+dirac)*source(i, j)-source(i+2, j) ;
	if (i==0) dirac = 0 ;
    }

    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = barre(i, j)-barre(i+2, j) ;

    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = tilde(i, j)-tilde(i+1, j) ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    }
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

		//-------------------
	       //--  R_CHEBI   -----
	      //-------------------


Matrice _cl_r_chebi (const Matrice &source, int l, double, int) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
       
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
	   cout << "_cl_r_chebi : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
   
    for (int i=0 ; i<n-2 ; i++)
	for (int j=0 ; j<n ; j++)
	    barre.set(i, j) = source(i, j)-source(i+2, j) ;

    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = barre(i, j)-barre(i+2, j) ;    

    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = tilde(i, j)-tilde(i+1, j) ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}
		//-------------------
	       //--  R_CHEBU   -----
	      //-------------------

Matrice _cl_r_chebu (const Matrice &source, int l, double, int puis) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
    Matrice res(n, n) ;
    res.set_etat_qcq() ;
    
    switch (puis) {
	case 5 :
	    res = _cl_r_chebu_cinq(source, l) ;
	    break ;
	case 4 :
	    res = _cl_r_chebu_quatre(source, l) ;
	    break ;
	case 3 :
	    res = _cl_r_chebu_trois (source, l) ;
	    break ;
	case 2 :
	    res = _cl_r_chebu_deux(source, l) ;
	    break ;
	default :
	    abort() ;
	    exit(-1) ;
    }
    
    return res ;
}


// Cas dzpuis = 4
Matrice _cl_r_chebu_quatre (const Matrice &source, int l) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
            
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
	   cout << "_cl_r_chebu_quatre : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	     barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j)) ;
	if (i==0) dirac = 0 ;
	}
    
    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = (barre(i, j)-barre(i+2, j)) ;
	    
    Matrice prime(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    prime.set(i, j) = (tilde(i, j)-tilde(i+1, j)) ;
    
    Matrice res(prime) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = (prime(i, j)-prime(i+2, j)) ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

// Cas dzpuis == 3
Matrice _cl_r_chebu_trois (const Matrice &source, int l) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
            
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
	   cout << "_cl_r_chebu_trois : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	     barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j)) ;
	if (i==0) dirac = 0 ;
	}
    
    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = (barre(i, j)-barre(i+2, j)) ;
    
    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = (tilde(i, j)+tilde(i+1, j)) ;
    
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}


//Cas dzpuis == 2
Matrice _cl_r_chebu_deux (const Matrice &source, int l) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
            
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
	   cout << "_cl_r_chebu_deux : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	     barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j)) ;
	if (i==0) dirac = 0 ;
	}
   
    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = (barre(i, j)-barre(i+2, j)) ;
	    
    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = (tilde(i, j)+tilde(i+1, j)) ;

    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}


//Cas dzpuis == 5
Matrice _cl_r_chebu_cinq (const Matrice &source, int l) {
    int n = source.get_dim(0) ;
    assert (n == source.get_dim(1)) ;
    
            
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
	   cout << "_cl_r_chebu_cinq : trop de matrices" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice barre(source) ;
  
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	for (int j=0 ; j<n ; j++)
	     barre.set(i, j) = ((1+dirac)*source(i, j)-source(i+2, j)) ;
	if (i==0) dirac = 0 ;
	}
   
    Matrice tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    tilde.set(i, j) = (barre(i, j)-barre(i+2, j)) ;
	    
    Matrice res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	for (int j=0 ; j<n ; j++)
	    res.set(i, j) = (tilde(i, j)+tilde(i+1, j)) ;

//     cout << "Apres comb. lin. : " << endl ;
//     cout << res ;
//     int fg ; cin >> fg ;

    return res ;
    } 
    
    // Cas ou le calcul a deja ete effectue :
    else
	return *tab[indice] ;
}

		//-------------------------
	       //- La routine a appeler ---
	      //---------------------------

Matrice combinaison (const Matrice &source, int l, double echelle, int puis, int base_r) {
    
		// Routines de derivation
    static Matrice (*combinaison[MAX_BASE])(const Matrice &, int, double, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    combinaison[i] = _cl_pas_prevu ;
	}
		// Les routines existantes
	combinaison[R_CHEB >> TRA_R] = _cl_r_cheb ;
	combinaison[R_CHEBU >> TRA_R] = _cl_r_chebu ;
	combinaison[R_CHEBP >> TRA_R] = _cl_r_chebp ;
	combinaison[R_CHEBI >> TRA_R] = _cl_r_chebi ;
	combinaison[R_JACO02 >> TRA_R] = _cl_r_jaco02 ;
    }
    
    Matrice res(combinaison[base_r](source, l, echelle, puis)) ;
    return res ;
}

//--------------------------------------------------
// Version Tbl --> Tbl a 1D pour la source 
//--------------------------------------------------


Tbl _cl_pas_prevu (const Tbl &source, int puis) {
     cout << "Combinaison lineaire pas prevue..." << endl ;
    cout << "source : " << &source << endl ;
    cout << "dzpuis : " << puis << endl ;
    abort() ;
    exit(-1) ;
    return source;
}



		//-------------------
	       //--  R_CHEB  -------
	      //--------------------

Tbl _cl_r_cheb (const Tbl &source, int) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	    barre.set(i) = ((1+dirac)*source(i)-source(i+2))
				/(i+1) ;
	if (i==0) dirac = 0 ;
    }
    
    Tbl res(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = barre(i)-barre(i+2) ;
   return res ;        
}


		//-------------------
	       //--  R_JACO02 ------
	      //-------------------

Tbl _cl_r_jaco02 (const Tbl &source, int) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    for (int i=0 ; i<n ; i++) {
	    barre.set(i) = source(i) ;
    }
    
    Tbl res(barre) ;
    for (int i=0 ; i<n ; i++)
	    res.set(i) = barre(i);
   return res ;        
}


		//-------------------
	       //--  R_CHEBP   -----
	      //-------------------

Tbl _cl_r_chebp (const Tbl &source, int) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	    barre.set(i) = (1+dirac)*source(i)-source(i+2) ;
	if (i==0) dirac = 0 ;
    }

    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = barre(i)-barre(i+2) ;

    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = tilde(i)-tilde(i+1) ;
	
   return res ;
}


		//-------------------
	       //--  R_CHEBI   -----
	      //-------------------

Tbl _cl_r_chebi (const Tbl &source, int) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    for (int i=0 ; i<n-2 ; i++)
	    barre.set(i) = source(i)-source(i+2) ;

    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = barre(i)-barre(i+2) ;    

    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = tilde(i)-tilde(i+1) ;
	
   return res ;
}


		//-------------------
	       //--  R_CHEBU   -----
	      //-------------------

Tbl _cl_r_chebu (const Tbl &source, int puis) {
    
    int n=source.get_dim(0) ;
    Tbl res(n) ;
    res.set_etat_qcq() ;
    
    switch(puis) {
	case 5 :
	    res = _cl_r_chebu_cinq(source) ;
	    break ;
	case 4 :
	    res = _cl_r_chebu_quatre(source) ;
	    break ;
	case 3 :
	    res = _cl_r_chebu_trois (source) ;
	    break ;
	case 2 :
	    res = _cl_r_chebu_deux(source) ;
	    break ;
	
	default :
	    abort() ;
	    exit(-1) ;    
    }
   return res ;
}

// Cas dzpuis = 4 ;
Tbl _cl_r_chebu_quatre (const Tbl &source) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	     barre.set(i) = ((1+dirac)*source(i)-source(i+2)) ;
	if (i==0) dirac = 0 ;
	}
    
    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = (barre(i)-barre(i+2)) ;
	    
    Tbl prime(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    prime.set(i) = (tilde(i)-tilde(i+1)) ;
    
     Tbl res(prime) ;
	for (int i=0 ; i<n-4 ; i++)
		res.set(i) = (prime(i)-prime(i+2)) ;
 
   return res ;
}
// cas dzpuis = 3
Tbl _cl_r_chebu_trois (const Tbl &source) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	     barre.set(i) = ((1+dirac)*source(i)-source(i+2)) ;
	if (i==0) dirac = 0 ;
	}
    
    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = (barre(i)-barre(i+2)) ;
	    
    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = (tilde(i)+tilde(i+1)) ;
 
   return res ;
}

// Cas dzpuis = 2 ;
Tbl _cl_r_chebu_deux (const Tbl &source) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	     barre.set(i) = ((1+dirac)*source(i)-source(i+2)) ;
	if (i==0) dirac = 0 ;
	}
    
    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = (barre(i)-barre(i+2)) ;
	    
    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = (tilde(i)+tilde(i+1)) ;
   return res ;
}

// Cas dzpuis = 5 ;
Tbl _cl_r_chebu_cinq (const Tbl &source) {
    Tbl barre(source) ;
    int n = source.get_dim(0) ;
    
    int dirac = 1 ;
    for (int i=0 ; i<n-2 ; i++) {
	     barre.set(i) = ((1+dirac)*source(i)-source(i+2)) ;
	if (i==0) dirac = 0 ;
	}
    
    Tbl tilde(barre) ;
    for (int i=0 ; i<n-4 ; i++)
	    tilde.set(i) = (barre(i)-barre(i+2)) ;
	    
    Tbl res(tilde) ;
    for (int i=0 ; i<n-4 ; i++)
	    res.set(i) = (tilde(i)+tilde(i+1)) ;
   return res ;
}


		//----------------------------
	       //- Routine a appeler        ---
	      //------------------------------

Tbl combinaison (const Tbl &source, int puis, int base_r) {
    
		// Routines de derivation
    static Tbl (*combinaison[MAX_BASE])(const Tbl &, int) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    combinaison[i] = _cl_pas_prevu ;
	}
		// Les routines existantes
	combinaison[R_CHEB >> TRA_R] = _cl_r_cheb ;
	combinaison[R_CHEBU >> TRA_R] = _cl_r_chebu ;
	combinaison[R_CHEBP >> TRA_R] = _cl_r_chebp ;
	combinaison[R_CHEBI >> TRA_R] = _cl_r_chebi ;
	combinaison[R_JACO02 >> TRA_R] = _cl_r_jaco02 ;
    }
    
    Tbl res(combinaison[base_r](source, puis)) ;
    return res ;
}
}
