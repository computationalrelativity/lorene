 /*
 *  Methods of class Iso_hor to compute sources for Psi, N y beta
 *
 *    (see file hor_isol.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004-2005  Jose Luis Jaramillo
 *                            Francois Limousin
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

char source_hor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.12  2005/03/28 19:42:39  f_limousin
 * Implement the metric and A^{ij}A_{ij} of Sergio for pertubations
 * of Kerr black holes.
 *
 * Revision 1.11  2005/03/22 13:25:36  f_limousin
 * Small changes. The angular velocity and A^{ij} are computed
 * with a differnet sign.
 *
 * Revision 1.10  2005/03/03 10:11:57  f_limousin
 * The function source_psi(), source_nn() and source_beta() are
 * now const and return an object const.
 *
 * Revision 1.9  2005/02/07 10:37:21  f_limousin
 * Minor changes.
 *
 * Revision 1.8  2004/12/31 15:35:18  f_limousin
 * Modifications to avoid warnings.
 *
 * Revision 1.7  2004/12/22 18:16:43  f_limousin
 * Many different changes.
 *
 * Revision 1.6  2004/11/05 10:59:07  f_limousin
 * Delete ener_dens, mom_dens and trace stress in functions
 * source_nn, source_psi and source_beta.
 * And some modification to avoid warnings (source_nn change to source...).
 *
 * Revision 1.5  2004/11/03 17:16:44  f_limousin
 * Delete argument trk_point for source_nn()
 *
 * Revision 1.4  2004/10/29 15:45:08  jl_jaramillo
 * Change name of functions
 *
 * Revision 1.3  2004/09/28 16:08:26  f_limousin
 * Minor modifs.
 *
 * Revision 1.2  2004/09/09 16:55:08  f_limousin
 * Add the 2 lines $Id$Log: for CVS
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "isol_hor.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "graphique.h"
#include "utilitaires.h"

const Scalar Isol_hor::source_psi() const{

    using namespace Unites ;
   
    // Initialisations
    // ---------------

    const Base_vect& triad = *(ff.get_triad()) ;
    
    Scalar tmp(mp) ;
    Scalar tmp_scal(mp) ; 
    Sym_tensor tmp_sym(mp, CON, triad) ;

    Scalar source(mp) ; 
       
    //===============================================
    //  Computations of the source for Psi 
    //===============================================
    
    const Vector& d_psi = psi().derive_cov(ff) ;       // D_i Psi
     
    // Source for Psi 
    // --------------
    tmp = 0.125* psi() * met_gamt.ricci_scal() 
      - contract(hh(), 0, 1, d_psi.derive_cov(ff), 0, 1 ) ;
    tmp.inc_dzpuis() ; // dzpuis : 3 -> 4
    
    tmp -= contract(hdirac(), 0, d_psi, 0) ;  
                
    source = tmp - psi()*psi4()* ( 0.125* aa_quad() 
				       - 8.33333333333333e-2* trK*trK ) ;
    source.annule_domain(0) ;

    return source ;

}


const Scalar Isol_hor::source_nn() const{

    using namespace Unites ;
   
    // Initialisations
    // ---------------
 
    const Base_vect& triad = *(ff.get_triad()) ;
    
    Scalar tmp(mp) ;
    Scalar tmp_scal(mp) ; 
    Sym_tensor tmp_sym(mp, CON, triad) ;

    Scalar source(mp) ; 
       
    //===============================================
    //  Computations of the source for NN 
    //===============================================
    
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& dnnn = nn().derive_cov(ff) ;         // D_i N
    
    // Source for N 
    // ------------

    source = psi4()*( nn()*( aa_quad() + 0.3333333333333333* trK*trK )
			 - trK_point ) 
	     - 2.* contract(dln_psi, 0, nn().derive_con(met_gamt), 0)  
    - contract(hdirac(), 0, dnnn, 0) ; 
        
    tmp = psi4()* contract(beta(), 0, trK.derive_cov(ff), 0) 
      - contract( hh(), 0, 1, dnnn.derive_cov(ff), 0, 1 ) ;
        
    tmp.inc_dzpuis() ; // dzpuis: 3 -> 4
        
    source += tmp ;

    source.annule_domain(0) ;

    return source ;

}



const Vector Isol_hor::source_beta() const {

    using namespace Unites ;
   
    // Initialisations
    // ---------------

    const Base_vect& triad = *(ff.get_triad()) ;
    
    Scalar tmp(mp) ;
    Scalar tmp_scal(mp) ; 
    Sym_tensor tmp_sym(mp, CON, triad) ;
    Vector tmp_vect(mp, CON, triad) ;

    Vector source(mp, CON, triad) ; 

    //===============================================
    //  Computations of the source for beta 
    //===============================================
    
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& dnnn = nn().derive_cov(ff) ;         // D_i N

    // Source for beta 
    // ---------------
    
    source = 2.* contract(aa(), 1, 
			       dnnn - 6.*nn() * dln_psi, 0) ;
                
    tmp_vect = 0.66666666666666666* trK.derive_con(met_gamt) ;
    tmp_vect.inc_dzpuis() ;

    source += 2.* nn() * ( tmp_vect
			    - contract(met_gamt.connect().get_delta(), 1, 2, 
					   aa(), 0, 1) ) ;
            
    Vector vtmp = contract(hh(), 0, 1, 
                           beta().derive_cov(ff).derive_cov(ff), 1, 2)
      + 0.3333333333333333*
      contract(hh(), 1, beta().divergence(ff).derive_cov(ff), 0) 
      - hdirac().derive_lie(beta()) 
	+ gamt_point.divergence(ff) ;      // zero in the Dirac gauge
    vtmp.inc_dzpuis() ; // dzpuis: 3 -> 4
    
    source -= vtmp ; 
        
    source += 0.66666666666666666* beta().divergence(ff) * hdirac() ;

    source.annule_domain(0) ;

    return source ;

}




  
     
        

