/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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
#include "mtbl_cf.h"
#include "matrice.h"
#include "type_parite.h"
#include "proto.h"
#include "base_val.h"
#include "utilitaires.h"



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
    const int nmax = 50 ;
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
