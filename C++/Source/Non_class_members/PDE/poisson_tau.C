/*
 *   Copyright (c) 2005 Philippe Grandclement
 
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.12  2018/11/16 14:34:36  j_novak
 * Changed minor points to avoid some compilation warnings.
 *
 * Revision 1.11  2016/12/05 16:18:10  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.10  2014/10/13 08:53:30  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/10/06 15:16:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2013/06/05 15:10:43  j_novak
 * Suppression of FINJAC sampling in r. This Jacobi(0,2) base is now
 * available by setting colloc_r to BASE_JAC02 in the Mg3d constructor.
 *
 * Revision 1.7  2008/08/27 08:51:15  jl_cornou
 * Added Jacobi(0,2) polynomials
 *
 * Revision 1.6  2007/12/14 10:19:34  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.4  2005/11/24 14:07:54  j_novak
 * Use of Matrice::annule_hard()
 *
 * Revision 1.3  2005/08/26 14:02:41  p_grandclement
 * Modification of the elliptic solver that matches with an oscillatory exterior solution
 * small correction in Poisson tau also...
 *
 * Revision 1.2  2005/08/25 12:16:01  p_grandclement
 * *** empty log message ***
 *
 * 
 * $Header$
 *
 */

// Header C : 
#include <cstdlib>
#include <cmath>

// Headers Lorene :
#include "matrice.h"
#include "map.h"
#include "proto.h"
#include "type_parite.h"



	    //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------

/*
 * 
 * Solution de l'equation de poisson with a multi-domain Tau method
 * 
 * Entree : mapping :   le mapping affine
 *	    source : les coefficients de la source qui a ete multipliee par
 *		    r^4 r^3 ou r^2 dans la ZEC.
 *		    La base de decomposition doit etre Ylm
 *	    dzpuis : exposant de r dans le factor multiplicatif dans la ZEC
 * Sortie : renvoie les coefficients de la solution dans la meme base de 
 *	    decomposition (a savoir Ylm)
 *	    
 */


namespace Lorene {
Mtbl_cf sol_poisson_tau(const Map_af& mapping, const Mtbl_cf& source, int dzpuis)
{
    
    // Verifications d'usage sur les zones
    int nz = source.get_mg()->get_nzone() ;
    assert (nz>1) ;
    assert ((source.get_mg()->get_type_r(0) == RARE) || (source.get_mg()->get_type_r(0) == FIN)) ;
    assert (source.get_mg()->get_type_r(nz-1) == UNSURR) ;
    for (int l=1 ; l<nz-1 ; l++)
	assert(source.get_mg()->get_type_r(l) == FIN) ;

     assert ((dzpuis==4) || (dzpuis==2) || (dzpuis==3)) ;
       
    // Bases spectrales
    const Base_val& base = source.base ;
    
    // Resultat
    Mtbl_cf resultat(source.get_mg(), base) ;
    resultat.annule_hard() ;
     
    // donnees sur la zone
    int nr, nt, np ;
    int base_r ;
    double alpha, beta, echelle ;
    int l_quant, m_quant;
    
    // Determination of the size of the systeme :
    int size = 0 ;
    int max_nr = 0 ;
    for (int l=0 ; l<nz ; l++) { 
    	nr = mapping.get_mg()->get_nr(l) ;
        size += nr ;
	if (nr > max_nr)
	    max_nr = nr ;
    }
	
    Matrice systeme (size, size) ;
    systeme.set_etat_qcq() ;
    Tbl sec_membre (size) ;
   
    np = mapping.get_mg()->get_np(0) ;
    nt = mapping.get_mg()->get_nt(0) ;
    Matrice* work ;
    
    // On bosse pour chaque l, m :
    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) 
	if (nullite_plm(j, nt, k, np, base) == 1) {
	
// 	for (int lig=0 ; lig<size ; lig++)
//             for (int col=0 ; col< size ; col++)
// 	    systeme.set(lig,col) = 0 ;
	systeme.annule_hard() ;
	sec_membre.annule_hard() ;
	     
	int column_courant = 0 ;
	int ligne_courant = 0 ;
	
        	//--------------------------
		//       NUCLEUS
		//--------------------------
		
	nr = mapping.get_mg()->get_nr(0) ; 
	alpha = mapping.get_alpha()[0] ;
        base.give_quant_numbers (0, k, j, m_quant, l_quant, base_r) ;  
        work = new Matrice (laplacien_mat(nr, l_quant, 0., 0, base_r)) ;

	int nbr_cl = 0 ;
	// RARE case	
	if (source.get_mg()->get_type_r(0) == RARE) {
	// regularity conditions :
	if (l_quant > 1) {
	     nbr_cl = 1 ;
	     if (l_quant%2==0) {
	        //Even case
		for (int col=0 ; col<nr ; col++)
		    if (col%2==0)
		        systeme.set(ligne_courant, col+column_courant) = 1 ;
		    else 
		        systeme.set(ligne_courant, col+column_courant) = -1 ;
		}
	     else {
	     //Odd case
	         for (int col=0 ; col<nr ; col++)
		    if (col%2==0)
		        systeme.set(ligne_courant, col+column_courant) = 2*col+1 ;
		    else 
		        systeme.set(ligne_courant, col+column_courant) = -(2*col+1) ;
		}
	  }
	}

	// JACO02 case
	else {
	  assert( base_r == R_JACO02) ;
	// regularity conditions :
	if (l_quant == 0) {
	    nbr_cl = 1 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) =  col*(col+1)*(col+2)*(col+3)/double(12)*(2*(col%2)-1);
		}
	    }
	else if (l_quant == 1) {
	    nbr_cl = 1 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) = (col+1)*(col+2)/double(2)*(1-2*(col%2)) ;
		}
	    }
	else {
	    nbr_cl = 2 ;
		for (int col=0 ; col<nr ; col++) {
		    systeme.set(ligne_courant, col+column_courant) = (col+1)*(col+2)/double(2)*(1-2*(col%2)) ;
		    systeme.set(ligne_courant+1, col+column_courant) = col*(col+1)*(col+2)*(col+3)/double(12)*(2*(col%2)-1) ;
		}
	    }
	}	
	ligne_courant += nbr_cl ;

	// L'operateur :
	for (int lig=0 ; lig<nr-1-nbr_cl ; lig++) {
	  for (int col=0 ; col<nr ; col++)
	    systeme.set(lig+ligne_courant,col+column_courant) = (*work)(lig,col) ;
	  sec_membre.set(lig+ligne_courant) = alpha*alpha*source(0, k, j, lig) ;
	}
	
	delete work ;
	ligne_courant += nr-1-nbr_cl ;
	
	// Le raccord :
	for (int col=0 ; col<nr ; col++) {
	  if (source.get_mg()->get_type_r(0) == RARE) {
	    // La fonction
	    systeme.set(ligne_courant, col+column_courant) = 1 ;
	    // Sa d�riv�e :
	    if (l_quant%2==0) {
	      systeme.set(ligne_courant+1, col+column_courant) = 4*col*col/alpha ;
	    }
	    else { 
	      systeme.set(ligne_courant+1, col+column_courant) = (2*col+1)*(2*col+1)/alpha ;
	    }
	  }
	  else {
	    // La fonction
	    systeme.set(ligne_courant, col+column_courant) = 1 ; 
	    // Sa dérivée :
	    systeme.set(ligne_courant+1, col+column_courant) = col*(col+3)/double(2)/alpha ;
	  }
	}
	
	column_courant += nr ;
	  
	  


		//--------------------------
		//       SHELLS
		//--------------------------
	for (int l=1 ; l<nz-1 ; l++) {
	
		nr = mapping.get_mg()->get_nr(l) ; 
		alpha = mapping.get_alpha()[l] ;
		beta = mapping.get_beta()[l] ;
		echelle = beta/alpha ;
		
        	base.give_quant_numbers (l, k, j, m_quant, l_quant, base_r) ;  
        	work = new Matrice (laplacien_mat(nr, l_quant, echelle, 0, base_r)) ;
	
	   // matching with previous domain :
	   for (int col=0 ; col<nr ; col++) {
	     // La fonction
	     if (col%2==0)
	          systeme.set(ligne_courant, col+column_courant) = -1 ;
	     else 
	          systeme.set(ligne_courant, col+column_courant) = 1 ;
	     // Sa d�riv�e :
	     if (col%2==0)
	         systeme.set(ligne_courant+1, col+column_courant) = col*col/alpha ;
	     else 
	         systeme.set(ligne_courant+1, col+column_courant) = -col*col/alpha ;
	  }
	  ligne_courant += 2 ;
	  
	  // L'operateur :
	  
	  // source must be multiplied by (x+echelle)^2
	  Tbl source_aux(nr) ;
	  source_aux.set_etat_qcq() ;
	  for (int i=0 ; i<nr ; i++)
	      source_aux.set(i) = source(l,k,j,i)*alpha*alpha ;
    	Tbl xso(source_aux) ;
    	Tbl xxso(source_aux) ;
    	multx_1d(nr, &xso.t, R_CHEB) ;
    	multx_1d(nr, &xxso.t, R_CHEB) ;
    	multx_1d(nr, &xxso.t, R_CHEB) ;
    	source_aux = beta*beta/alpha/alpha*source_aux+2*beta/alpha*xso+xxso ;
	  
	for (int lig=0 ; lig<nr-2 ; lig++) {
	  for (int col=0 ; col<nr ; col++)
	    systeme.set(lig+ligne_courant,col+column_courant) = (*work)(lig,col) ;
	  sec_membre.set(lig+ligne_courant) = source_aux(lig) ;
	}
	
	delete work ;
	ligne_courant += nr-2 ;
	// Matching with the next domain :
	for (int col=0 ; col<nr ; col++) {
	  // La fonction
	  systeme.set(ligne_courant, col+column_courant) = 1 ;
	  // Sa d�riv�e :
	  systeme.set(ligne_courant+1, col+column_courant) = col*col/alpha ;
	}
	
	column_courant += nr ;   
	}
	  
	  
		//--------------------------
		//       ZEC
		//--------------------------
	nr = mapping.get_mg()->get_nr(nz-1) ; 
	alpha = mapping.get_alpha()[nz-1] ;
	beta = mapping.get_beta()[nz-1] ;
		
	base.give_quant_numbers (nz-1, k, j, m_quant, l_quant, base_r) ;  
	work = new Matrice(laplacien_mat(nr, l_quant, 0., dzpuis, base_r)) ;
	
	// Matching with the previous domain :
	 for (int col=0 ; col<nr ; col++) {
	     // La fonction
	     if (col%2==0)
	          systeme.set(ligne_courant, col+column_courant) = -1 ;
	     else 
	          systeme.set(ligne_courant, col+column_courant) = 1 ;
	     // Sa d�riv�e :
	     if (col%2==0)
	         systeme.set(ligne_courant+1, col+column_courant) = -4*alpha*col*col ;
	     else 
	         systeme.set(ligne_courant+1, col+column_courant) = 4*alpha*col*col ;
	  }
	  ligne_courant += 2 ;	
	  
	  // Regularity and BC at infinity ?
	  nbr_cl =0 ;
	  switch (dzpuis) {
	       case 4 : 
	           if (l_quant==0) {
		       nbr_cl = 1 ;
		       // Only BC at infinity :
		       for (int col=0 ; col<nr ; col++)
		           systeme.set(ligne_courant, col+column_courant) = 1 ;
		       }
		   else { 
		       nbr_cl = 2 ;
		       // BC at infinity :
		       for (int col=0 ; col<nr ; col++)
		           systeme.set(ligne_courant, col+column_courant) = 1 ;
		       // Regularity :
		       for (int col=0 ; col<nr ; col++)
		           systeme.set(ligne_courant+1, col+column_courant) = -4*alpha*col*col ;   
		    }
		    break ;
	       
	        case 3 :
		    nbr_cl = 1 ;
		    // Only BC at infinity :
		    for (int col=0 ; col<nr ; col++)
		       systeme.set(ligne_courant, col+column_courant) = 1 ;
		    break ;
		
		case 2 :
		     if (l_quant==0) {
		         nbr_cl = 1 ;
		        // Only BC at infinity :
		        for (int col=0 ; col<nr ; col++)
		           systeme.set(ligne_courant, col+column_courant) = 1 ;
			}
		     break ;
		default : 
		    cout << "Unknown dzpuis in sol_poisson_tau ..." << endl ;
		    abort() ;
	}
	
	ligne_courant += nbr_cl ;
	
	// Multiplication of the source :
	double indic = 1 ;
	switch (dzpuis) {
	    case 4 : 
	        indic = alpha*alpha ;
		break ;
	    case 3 :
	        indic = alpha ;
		break ;
	default : 
	     break ;
	}
	
	// L'operateur :
	for (int lig=0 ; lig<nr-1-nbr_cl ; lig++) {
	    for (int col=0 ; col<nr ; col++)
	       systeme.set(lig+ligne_courant,col+column_courant) = (*work)(lig,col) ;
	    sec_membre.set(lig+ligne_courant) = indic*source(nz-1, k, j, lig) ;
	}
	delete work ;
	
	// Solving the system:
	systeme.set_band (max_nr, max_nr) ;
	systeme.set_lu() ;
	Tbl soluce (systeme.inverse(sec_membre)) ;
	
	// On range :
	int conte = 0 ;
	for (int l=0 ; l<nz ; l++) {
	     nr = mapping.get_mg()->get_nr(l) ;
	     for (int i=0 ; i<nr ; i++) {
	         resultat.set(l,k,j,i) = soluce(conte) ;
		 conte ++ ;
		}
	}
	
    }

    return resultat ;
}
}
