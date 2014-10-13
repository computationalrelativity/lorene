
/*
 *   Copyright (c) 2004 Francois Limousin
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


char poisson_interne_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:29  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:16:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2004/11/23 12:51:42  f_limousin
 * Minor changes.
 *
 * Revision 1.1  2004/03/31 11:36:15  f_limousin
 * First version
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
#include "tbl.h"
#include "mtbl_cf.h"
#include "map.h"
#include "base_val.h"
#include "proto.h"
#include "type_parite.h"
#include "utilitaires.h"




            //----------------------------------------------
	   //		Version Mtbl_cf
	  //----------------------------------------------

namespace Lorene {
Mtbl_cf sol_poisson_interne (const Map_af& mapping, 
    const Mtbl_cf& source, const Mtbl_cf& lim_der){

    int nz = source.get_mg()->get_nzone() ;

    assert(source.get_mg()->get_type_r(0) == RARE) ;
    assert (lim_der.get_mg() == source.get_mg()->get_angu()) ;
    assert (source.get_etat() != ETATNONDEF) ;
    assert (lim_der.get_etat() != ETATNONDEF) ;
     
    // Bases spectrales
    const Base_val& base = source.base ;
    
    // donnees sur la zone
    int nr = source.get_mg()->get_nr(0) ;
    int nt = source.get_mg()->get_nt(0) ;
    int np = source.get_mg()->get_np(0) ;;
    int base_r ;
    int l_quant, m_quant;
    
    double alpha = mapping.get_alpha()[0] ;
    double beta = mapping.get_beta()[0] ;
    double facteur ;
    
    //Rangement des valeurs intermediaires 
    Tbl *so ;
    Tbl *sol_hom ;
    Tbl *sol_part ;
    Matrice *operateur ;
    Matrice *nondege ;
    
    
    Mtbl_cf resultat(source.get_mg(), base) ;
    resultat.annule_hard() ;
    
    for (int k=0 ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++) 
	    if (nullite_plm(j, nt, k, np, base) == 1)
	    {
		// calcul des nombres quantiques :
		donne_lm(nz, 0, j, k, base, m_quant, l_quant, base_r) ;
	    
		// Construction de l'operateur
		operateur = new Matrice(laplacien_mat
				    (nr, l_quant, 0., 0, base_r)) ;
		
		(*operateur) = combinaison(*operateur, l_quant, 0.,0, base_r) ;
		
		 // Operateur inversible
		nondege = new Matrice(prepa_nondege(*operateur, l_quant, 
							0., 0, base_r)) ;  
		
		// Calcul DE LA SH
		sol_hom = new Tbl(solh(nr, l_quant, 0., base_r)) ;
		
		// Calcul de la SP
		so = new Tbl(nr) ;
		so->set_etat_qcq() ;
		for (int i=0 ; i<nr ; i++)
		    so->set(i) = source(0, k, j, i) ;
		
		sol_part = new Tbl (solp(*operateur, *nondege, alpha,
					 beta, *so, 0, base_r)) ;

		//-------------------------------------------
		// On est parti pour imposer la boundary
		//-------------------------------------------    

		// Condition de raccord de type Neumann :
		double val_der_solp = 0 ;
		for (int i=0 ; i<nr ; i++)
		    if (m_quant%2 == 0)
			val_der_solp += (2*i)*(2*i)*(*sol_part)(i)/alpha ;
		    else
			val_der_solp += (2*i+1)*(2*i+1)*(*sol_part)(i)/alpha ;

		double val_der_solh = 0 ;
		for (int i=0 ; i<nr ; i++)
	 	    if (m_quant%2 == 0)
			val_der_solh += (2*i)*(2*i)*(*sol_hom)(i)/alpha ;
		    else
			val_der_solh += (2*i+1)*(2*i+1)*(*sol_hom)(i)/alpha ;

		if (l_quant != 0){
		    assert (val_der_solh != 0) ;

		    facteur = (lim_der(0, k, j, 0)-val_der_solp) /
			val_der_solh ;
		    
		    for (int i=0 ; i<nr ; i++)
			sol_part->set(i) += facteur*(*sol_hom)(i) ;
		}
		else {
		    for (int i=0 ; i<nr ; i++)
			sol_part->set(i) = 0. ;
		}
		    

		// solp contient le bon truc (normalement ...)
		for (int i=0 ; i<nr ; i++)
		    resultat.set(0, k, j, i) = (*sol_part)(i) ;
		
		delete operateur ;
		delete nondege ;
		delete so ;
		delete sol_hom ;
		delete sol_part ;
	    }
    
    return resultat ;
}
}
