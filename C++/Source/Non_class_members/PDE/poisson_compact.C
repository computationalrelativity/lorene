/*
 *   Copyright (c) 2000-2006 Philippe Grandclement
 *   Copyright (c) 2007 Michal Bejger
 *   Copyright (c) 2007 Eric Gourgoulhon
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


char poisson_compact_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2007/10/16 21:54:23  e_gourgoulhon
 * Added new function sol_poisson_compact (multi-domain version).
 *
 * Revision 1.3  2006/09/05 13:39:46  p_grandclement
 * update of the bin_ns_bh project
 *
 * Revision 1.2  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2000/03/16  16:28:06  phil
 * Version entirement revue et corrigee
 *
 * Revision 2.1  2000/03/09  13:51:55  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/03/09  13:44:56  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Headers Lorene
#include "map.h"
#include "diff.h"
#include "matrice.h"
#include "type_parite.h"
#include "proto.h"
#include "base_val.h"
#include "utilitaires.h"

	/////////////////////////////
	//	Single domain version  //
	/////////////////////////////

/*
 * Cette fonction resout, dans le noyau :
 *	a*(1-xi^2)*lap(uu)+b*xi*duu/dxi = source
 * avec a>0 et b<0 ;
 * Pour le stokage des operateurs, il faut faire reamorce = true au
 * debut d'un nouveau calcul.
 */

Mtbl_cf sol_poisson_compact(const Mtbl_cf& source, double a, double b, 
			    bool reamorce)  {
				    
    // Verifications :
    assert (source.get_etat() != ETATNONDEF) ;
    
    assert (a>0) ;
    assert (b<0) ;
    
    // Les tableaux de stockage :
    const int nmax = 200 ;
    static Matrice* tab_op[nmax] ;
    static int nb_deja_fait = 0 ;
    static int l_deja_fait[nmax] ;
    static int n_deja_fait[nmax] ;
    
    if (reamorce) {
	for (int i=0 ; i<nb_deja_fait ; i++)
	    delete tab_op[i] ;
	nb_deja_fait = 0 ;
    }
    
    int nz = source.get_mg()->get_nzone() ;

    // Pour le confort (on ne travaille que dans le noyau) :
    int nr = source.get_mg()->get_nr(0) ;
    int nt = source.get_mg()->get_nt(0) ;
    int np = source.get_mg()->get_np(0) ;
     
    int l_quant ;
    int m_quant ;
    int base_r ;
    
    // La solution ...    
    Mtbl_cf solution(source.get_mg(), source.base) ;
    solution.set_etat_qcq() ;
    solution.t[0]->set_etat_qcq() ;
    
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++)
	    if (nullite_plm(j, nt, k, np, source.base) == 1) 
	    {
		 // calcul des nombres quantiques :
	    donne_lm(nz, 0, j, k, source.base, m_quant, l_quant, base_r) ;
	    
        
	    //On gere le cas l_quant == 0 (c'est bien simple y'en a pas !)
	    if (l_quant == 0) {
		for (int i=0 ; i<nr ; i++)
		    solution.set(0, k, j, i) = 0 ;
		}
	    
	    // cas l_quant != 0
	    else {
		// On determine si la matrice a deja ete calculee :
	    int indice = -1 ;
	    
	    // Le cas l==1 est non singulier : pas de base de Gelerkin
	    int taille = (l_quant == 1) ? nr : nr-1 ;
	    
	    Matrice operateur (taille, taille) ;
	    for (int conte=0 ; conte<nb_deja_fait ; conte++)
		if ((l_deja_fait[conte]== l_quant) && (n_deja_fait[conte] == nr))
		    indice = conte ;
	    
	    if (indice == -1) {
		if (nb_deja_fait >= nmax) {
		    cout << "sol_poisson_compact : trop de matrices ..." << endl;
		    abort() ;
		    } 
		// Calcul a faire :
		operateur = a*lap_cpt_mat (nr, l_quant, base_r) 
			    + b*xdsdx_mat(nr, l_quant, base_r) ;
		operateur = combinaison_cpt (operateur, l_quant, base_r) ;
		
		l_deja_fait[nb_deja_fait] = l_quant ;
		n_deja_fait[nb_deja_fait] = nr ;
		tab_op[nb_deja_fait] = new Matrice(operateur) ;

		nb_deja_fait++ ;
		}
	    else {
		// rien a faire :
		operateur = *tab_op[indice] ;
		}
		
	   // La source :
	    Tbl so(taille) ;
	    so.set_etat_qcq() ;
	    for (int i=0 ; i<taille ; i++)
		so.set(i) = source(0, k, j, i) ;
	    so = combinaison_cpt (so, base_r) ;
	    
	    Tbl part (operateur.inverse(so)) ;
	    
	    if (taille == nr)
		for (int i=0 ; i<nr ; i++)
		     solution.set(0, k, j, i) = part(i) ; // cas l==1
	    else {
		solution.set(0, k, j, nr-1) = 0 ;
		for (int i=nr-2 ; i>=0 ; i--)
		    if (base_r == R_CHEBP) { //Gelerkin pair
			solution.set(0, k, j, i) = part(i) ;
			solution.set(0, k, j, i+1) += part(i) ;
			}
		    else { //Gelerkin impaire
	    		solution.set(0, k, j, i) = part(i)*(2*i+3) ;
			solution.set(0, k, j, i+1) += part(i)*(2*i+1) ;
			}
		    }
		}
	    }
	    else   // cas ou nullite_plm = 0 :
		for (int i=0 ; i<nr ; i++)
		    solution.set(0, k, j, i) = 0 ; 
	    
	// Mise a zero du coefficient (inusite) k=np+1
    for (int j=0; j<nt; j++)
	for(int i=0 ; i<nr ; i++)
	    solution.set(0, np+1, j, i) = 0 ; 
	
    for (int zone=1 ; zone<nz ; zone++)
	    solution.t[zone]->set_etat_zero() ;

    return solution ;
}



	/////////////////////////////
	//	Multi domain version  //
	/////////////////////////////


Mtbl_cf sol_poisson_compact(const Map_af& mapping, const Mtbl_cf& source, 
		const Tbl& ac, const Tbl& bc, bool )  {

	// Number of domains inside the star :
	int nzet = ac.get_dim(1) ; 
	assert(nzet<=source.get_mg()->get_nzone()) ; 
	
    // Some checks
    assert (source.get_mg()->get_type_r(0) == RARE) ;
    for (int l=1 ; l<nzet ; l++)
		assert(source.get_mg()->get_type_r(l) == FIN) ;

    // Spectral bases
    const Base_val& base = source.base ;
    
    // Result
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
    for (int l=0 ; l<nzet ; l++) { 
    	nr = mapping.get_mg()->get_nr(l) ;
        size += nr ;
		if (nr > max_nr) max_nr = nr ;
    }
	
    Matrice systeme (size, size) ;
    systeme.set_etat_qcq() ;
    Tbl sec_membre (size) ;
    Tbl soluce (size) ;
    soluce.set_etat_qcq() ;

    np = mapping.get_mg()->get_np(0) ;
    nt = mapping.get_mg()->get_nt(0) ;
    Matrice* work ;
    
    // On bosse pour chaque l, m :
    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++)
      	if (nullite_plm(j, nt, k, np, base) == 1) {
	
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

    if (l_quant == 0) {
		for (int i=0 ; i<size ; i++)
 			soluce.set(i) = 0 ;
	}  
	else {

    Diff_dsdx2 		d2_n(base_r, nr) ;		// suffix _n stands for "nucleus"
    Diff_sxdsdx 	sxd_n(base_r, nr) ;
   	Diff_sx2 		sx2_n(base_r, nr) ;
    Diff_x2dsdx2	x2d2_n(base_r,nr) ;
    Diff_xdsdx		xd_n(base_r,nr) ;
    Diff_id			id_n(base_r,nr) ;  
	
	work = new Matrice( ac(0,0) * ( d2_n + 2.*sxd_n -l_quant*(l_quant+1)*sx2_n )
		+ ac(0,2) * ( x2d2_n + 2.*xd_n -l_quant*(l_quant+1)*id_n ) 
		+ alpha * bc(0,1) * xd_n ) ;

	// cout << *work << endl ; 
	// arrete() ; 

	// regularity conditions :
	int nbr_cl = 0 ;
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
	  ligne_courant ++ ;
	}
	
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
	     // La fonction
	     systeme.set(ligne_courant, col+column_courant) = 1 ;
	     // Sa dérivée :
	     if (l_quant%2==0)
	         systeme.set(ligne_courant+1, col+column_courant) = 4*col*col/alpha ;
	     else 
	         systeme.set(ligne_courant+1, col+column_courant) = (2*col+1)*(2*col+1)/alpha ;
	}
	column_courant += nr ;
	
		//--------------------------
		//       SHELLS
		//--------------------------

	for (int l=1 ; l<nzet ; l++) {
	
		nr = mapping.get_mg()->get_nr(l) ; 
		alpha = mapping.get_alpha()[l] ;
		beta = mapping.get_beta()[l] ;
		echelle = beta/alpha ;
		double bsa = echelle ; 
		double bsa2 = bsa*bsa ; 
		
        base.give_quant_numbers (l, k, j, m_quant, l_quant, base_r) ;  

		Diff_dsdx		dx(base_r, nr) ; 
	 	Diff_xdsdx		xdx(base_r, nr) ; 
	 	Diff_x2dsdx		x2dx(base_r, nr) ; 
	 	Diff_x3dsdx		x3dx(base_r, nr) ; 
	 	Diff_dsdx2		dx2(base_r, nr) ; 
	 	Diff_xdsdx2		xdx2(base_r, nr) ; 
	 	Diff_x2dsdx2	x2dx2(base_r, nr) ; 
	 	Diff_x3dsdx2	x3dx2(base_r, nr) ; 
    	Diff_id			id(base_r,nr) ;  
    	Diff_mx			mx(base_r,nr) ;  
	
        work = new Matrice ( ac(l,0) * ( x2dx2 + 2.*bsa*xdx2 + bsa2*dx2 + 2.*xdx 
			+ 2.*bsa*dx - l_quant*(l_quant+1)*id )
			+ ac(l,1) * ( x3dx2 + 2.*bsa*x2dx2 + bsa2*xdx2 + 2.*x2dx 
			+ 2.*bsa*xdx - l_quant*(l_quant+1) *mx )
			+ alpha * ( bc(l,0) * ( x2dx + 2.*bsa*xdx + bsa2*dx ) 
					  + bc(l,1) * ( x3dx + 2.*bsa*x2dx + bsa2*xdx ) ) ) ;
	
	   // matching with previous domain :
	   for (int col=0 ; col<nr ; col++) {
	     // La fonction
	     if (col%2==0)
	          systeme.set(ligne_courant, col+column_courant) = -1 ;
	     else 
	          systeme.set(ligne_courant, col+column_courant) = 1 ;
	     // Sa dérivée :
	     if (col%2==0)
	         systeme.set(ligne_courant+1, col+column_courant) = col*col/alpha ;
	     else 
	         systeme.set(ligne_courant+1, col+column_courant) = -col*col/alpha ;
		}
	  ligne_courant += 2 ;
	  
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
	  
	  // L'operateur :
	  
		for (int lig=0 ; lig<nr-1 ; lig++) {
			for (int col=0 ; col<nr ; col++)
			systeme.set(lig+ligne_courant,col+column_courant) = (*work)(lig,col) ;
			sec_membre.set(lig+ligne_courant) = source_aux(lig) ;
		}
	  
		// cout << *work << endl ; 
		// arrete() ; 

	  delete work ;
	  ligne_courant += nr-2 ;

		if (l<nzet-1) {  // if this not the last shell
	  		//  matching with the next domain 
			for (int col=0 ; col<nr ; col++) {
			// La fonction
			systeme.set(ligne_courant, col+column_courant) = 1 ;
			// Sa dérivée :
			systeme.set(ligne_courant+1, col+column_courant) = col*col/alpha ;
			}
		}
	     
	  column_courant += nr ;   

	  } // end loop on the shells (index l)
	  
	// cout << "systeme : " << systeme << endl ; 
	// arrete() ;   	
	
	// Solving the system:
 
	systeme.set_band (size, size) ;   
	systeme.set_lu() ;

	// cout << "Determinant: " << systeme.determinant() << endl ; 
	// cout << "Eigenvalues: " << systeme.val_propre() << endl ;

	soluce = systeme.inverse(sec_membre) ;
	
	} // end case l_quant != 0 

	// cout << "soluce: " << soluce << endl ; 
	// arrete() ; 

	// On range :
	int conte = 0 ;
	for (int l=0 ; l<nzet ; l++) {
	     nr = mapping.get_mg()->get_nr(l) ;
	     for (int i=0 ; i<nr ; i++) {
	        resultat.set(l,k,j,i) = soluce(conte) ;
		 	conte++ ;
		}
	}

	} // end case nullite_plm == 1

	return resultat ;

}















