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


char solh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2014/10/13 08:53:30  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/10/06 15:16:10  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2008/02/18 13:53:43  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.7  2007/12/20 09:11:09  jl_cornou
 * Correction of an error in op_sxpun about Jacobi(0,2) polynomials
 *
 * Revision 1.6  2007/12/13 15:48:46  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.5  2007/12/12 12:30:49  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.4  2004/10/05 15:44:21  j_novak
 * Minor speed enhancements.
 *
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
 * Revision 2.13  2000/01/18  14:15:50  phil
 * enleve assert sur nobre de points min en r
 *
 * Revision 2.12  2000/01/04  18:59:39  phil
 * Double nmax
 *
 * Revision 2.11  1999/10/11  14:29:37  phil
 * &-> &&
 *
 * Revision 2.10  1999/09/30  09:23:25  phil
 * remplacement des && en &
 *
 * Revision 2.9  1999/09/17  15:26:25  phil
 * correction definition de NMAX
 *
 * Revision 2.8  1999/06/23  12:35:00  phil
 * *** empty log message ***
 *
 * Revision 2.7  1999/04/28  10:49:19  phil
 * augmentation de nmax a 50
 *
 * Revision 2.6  1999/04/27  13:12:34  phil
 * *** empty log message ***
 *
 * Revision 2.5  1999/04/19  14:07:02  phil
 * *** empty log message ***
 *
 * Revision 2.4  1999/04/16  13:19:07  phil
 * *** empty log message ***
 *
 * Revision 2.3  1999/04/16  13:17:02  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/04/14  13:57:05  phil
 * Sauvegarde des solutions deja calculees
 *
 * Revision 2.1  1999/04/07  14:39:23  phil
 * esthetique
 *
 * Revision 2.0  1999/04/07  14:11:18  phil
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

#include "proto.h"
#include "matrice.h"
#include "type_parite.h"


/*
 * 
 * Renvoie une ou 2 solutions homogenes
 *  Si base_r = R_CHEB deux solutions (x+echelle)^l dans (0, *) et
 *				1/(x+echelle)^(l+1) dans (1, *)
 *  Si base_r = R_CHEBU 1 solution (x-1)^l+1 dans (*) 
 *  Si base_r = R_CHEBP ou R_CHEBI x^l dans (*)
 * 
 * Entree : 
 *	n : nbre de points en r
 *	l : nbre quantique associe
 *	echelle : cf ci-dessus, utile que dans le cas R_CHEB
 *	base_r : base de decomposition
 * 
 * Sortie :
 *	Tbl contenant les coefficient de la ou des solution homogenes
 * 
 */

		//------------------------------------
		// Routine pour les cas non prevus --
		//------------------------------------
namespace Lorene {
Tbl _solh_pas_prevu (int n, int l, double echelle) {

    cout << " Solution homogene pas prevue ..... : "<< endl ;
    cout << " N : " << n << endl ;
    cout << " l : " << l << endl ;
    cout << " echelle : " << echelle << endl ;
    abort() ;
    exit(-1) ;
    Tbl res(1) ;
    return res;
}
	
	
		//-------------------
	       //--  R_CHEB   ------
	      //-------------------

Tbl _solh_r_cheb (int n, int l, double echelle) {
                
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_solh_r_cheb : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    	
  //  assert (l < n) ;
    
    tab[nb_dejafait] = new Tbl(2, n) ;
    Tbl* pres = tab[nb_dejafait] ;
    pres->set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    //Construction de la premiere solution homogene :
    // cad celle polynomiale.
    
    if (l==0) {
	pres->set(0, 0) = 1 ;
	for (int i=1 ; i<n ; i++)
	    pres->set(0, i) = 0 ;
	    }
    else {
	for (int i=0 ; i<n ; i++)
	    coloc[i] = pow(echelle-cos(M_PI*i/(n-1)), double(l)) ;
	
	cfrcheb(deg, deg, coloc, deg, coloc) ;
	for (int i=0 ; i<n ;i++)
	    pres->set(0, i) = coloc[i] ;
	}
    
    
    // construction de la seconde solution homogene :
    // cad celle fractionnelle.
    for (int i=0 ; i<n ; i++)
	coloc[i] = 1/pow(echelle-cos(M_PI*i/(n-1)), double(l+1)) ;
	
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
	pres->set(1, i) = coloc[i] ;	
    
    
    delete [] coloc ;
    delete [] deg ;
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   }
   
   return *tab[indice] ;
}	
	

		//-------------------
	       //--  R_JACO02 ------
	      //-------------------

Tbl _solh_r_jaco02 (int n, int l, double echelle) {
                
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_solh_r_jaco02 : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    	
  //  assert (l < n) ;
    
    tab[nb_dejafait] = new Tbl(2, n) ;
    Tbl* pres = tab[nb_dejafait] ;
    pres->set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    

    double* zeta = pointsgausslobatto(n-1) ;
    //Construction de la premiere solution homogene :
    // cad celle polynomiale.
    
    if (l==0) {
	pres->set(0, 0) = 1 ;
	for (int i=1 ; i<n ; i++)
	    pres->set(0, i) = 0 ;
	    }
    else {
	for (int i=0 ; i<n ; i++)
	    coloc[i] = pow((echelle + zeta[i]), double(l)) ;
	
	cfrjaco02(deg, deg, coloc, deg, coloc) ;
	for (int i=0 ; i<n ;i++)
	    pres->set(0, i) = coloc[i] ;
	}
    
    
    // construction de la seconde solution homogene :
    // cad celle fractionnelle.
    for (int i=0 ; i<n ; i++)
	coloc[i] = 1/pow((echelle + zeta[i]), double(l+1)) ;
	
    cfrjaco02(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
	pres->set(1, i) = coloc[i] ;	
    
    
    delete [] coloc ;
    delete [] deg ;
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   }
   
   return *tab[indice] ;
}	



		//-------------------
	       //--  R_CHEBP  ------
	      //-------------------

Tbl _solh_r_chebp (int n, int l, double) {
   
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_solh_r_chebp : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
   	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    assert (div(l, 2).rem ==0) ;
    
  //  assert (l < 2*n-1) ;
    
    tab[nb_dejafait] = new Tbl(n) ;
    Tbl* pres = tab[nb_dejafait] ;
    pres->set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    if (l==0) {
	pres->set(0) = 1 ;
	for (int i=1 ; i<n ; i++)
	    pres->set(i) = 0 ;
	    }
    else {
	for (int i=0 ; i<n ; i++)
	    coloc[i] = pow(sin(M_PI*i/2/(n-1)), double(l)) ;
	
	cfrchebp(deg, deg, coloc, deg, coloc) ;
	for (int i=0 ; i<n ;i++)
	    pres->set(i) = coloc[i] ;
	}
	
    delete [] coloc ;
    delete [] deg ;
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   }
    
   return *tab[indice] ;
}
	
	
	      	//-------------------
	       //--  R_CHEBI   -----
	      //-------------------
	
Tbl _solh_r_chebi (int n, int l, double) {
       
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_solh_r_chebi : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
       }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
    assert (div(l, 2).rem == 1)  ;
   // assert (l < 2*n) ;
    
    tab[nb_dejafait] = new Tbl(n) ;
    Tbl* pres = tab[nb_dejafait] ;
    pres->set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    if (l==1) {
	pres->set(0) = 1 ;
	for (int i=1 ; i<n ; i++)
	    pres->set(i) = 0 ;
	    }
    else {
	for (int i=0 ; i<n ; i++)
	    coloc[i] = pow(sin(M_PI*i/2/(n-1)), double(l)) ;
	
	cfrchebi(deg, deg, coloc, deg, coloc) ;
	for (int i=0 ; i<n ;i++)
	    pres->set(i) = coloc[i] ;
	}
	
    delete [] coloc ;
    delete [] deg ;
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   }
    
   return *tab[indice] ;
}
	
	
	
	       	//-------------------
	       //--  R_CHEBU   -----
	      //-------------------
	
Tbl _solh_r_chebu (int n, int l, double) {
           
   const int nmax = 200 ; // Nombre de Tbl stockes
   static Tbl* tab[nmax] ;  // les Tbl calcules
   static int nb_dejafait = 0 ; // nbre de Tbl calcules
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
	   cout << "_solh_r_chebu : trop de Tbl" << endl ;
	   abort() ;
	   exit (-1) ;
	  }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    
  //  assert (l < n-1) ;
    tab[nb_dejafait] = new Tbl(n) ;
    Tbl* pres = tab[nb_dejafait] ;
    pres->set_etat_qcq() ;
    double* coloc = new double[n] ;
    
    int * deg = new int[3] ;
    deg[0] = 1 ; 
    deg[1] = 1 ;
    deg[2] = n ;
    
    for (int i=0 ; i<n ; i++)
	coloc[i] = pow(-1-cos(M_PI*i/(n-1)), double(l+1)) ;
	
    cfrcheb(deg, deg, coloc, deg, coloc) ;
    for (int i=0 ; i<n ;i++)
	pres->set(i) = coloc[i] ;
	
    delete [] coloc ;
    delete [] deg ;
    indice = nb_dejafait ;
    nb_dejafait ++ ;
   }
    
   return *tab[indice] ;
}
	
	
	
	
	      	//-------------------
	       //--  Fonction   ----
	      //-------------------
	      
	      
Tbl solh(int n, int l, double echelle, int base_r) {

		// Routines de derivation
    static Tbl (*solh[MAX_BASE])(int, int, double) ;
    static int nap = 0 ;

		// Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    solh[i] = _solh_pas_prevu ;
	}
		// Les routines existantes
	solh[R_CHEB >> TRA_R] = _solh_r_cheb ;
	solh[R_CHEBU >> TRA_R] = _solh_r_chebu ;
	solh[R_CHEBP >> TRA_R] = _solh_r_chebp ;
	solh[R_CHEBI >> TRA_R] = _solh_r_chebi ;
	solh[R_JACO02 >> TRA_R] = _solh_r_jaco02 ;
    }
    
    return solh[base_r](n, l, echelle) ;
}
}
