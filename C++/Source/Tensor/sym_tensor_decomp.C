/*
 *  Methods transverse( ) and longit_pot( ) of class Sym_tensor
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003-2004  Eric Gourgoulhon & Jerome Novak
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

char sym_tensor_decomp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2004/03/29 16:13:07  j_novak
 * New methods set_longit_trans and set_tt_trace .
 *
 * Revision 1.5  2004/02/09 12:56:27  e_gourgoulhon
 * Method longit_pot: added test of the vector Poisson equation.
 *
 * Revision 1.4  2004/02/02 09:18:11  e_gourgoulhon
 * Method longit_pot: treatment of case divergence dzpuis = 5.
 *
 * Revision 1.3  2003/12/10 10:17:54  e_gourgoulhon
 * First operational version.
 *
 * Revision 1.2  2003/11/27 16:01:47  e_gourgoulhon
 * First implmentation.
 *
 * Revision 1.1  2003/11/26 21:57:03  e_gourgoulhon
 * First version; not ready yet.
 *
 *
 * $Header$
 *
 */


// C headers
#include <stdlib.h>

// Lorene headers
#include "tensor.h"
#include "metric.h"

void Sym_tensor::set_longit_trans(const Sym_tensor_trans& ht, 
				  const Vector& v_pot) {

  assert ( v_pot.get_index_type(0) == CON ) ;

  const Metric& metre = ht.get_met_div() ;

  *this = ht + v_pot.ope_killing(metre) ;

  del_deriv() ;

  set_dependance(metre) ;
  int jp = get_place_met(metre) ;
  assert ((jp>=0) && (jp<N_MET_MAX)) ;

  p_transverse[jp] = new Sym_tensor_trans(ht) ;
  p_longit_pot[jp] = new Vector( v_pot ) ;
  
}

const Sym_tensor_trans& Sym_tensor::transverse(const Metric& metre) const {

    set_dependance(metre) ;
    int jp = get_place_met(metre) ;
    assert ((jp>=0) && (jp<N_MET_MAX)) ;

    if (p_transverse[jp] == 0x0) { // a new computation is necessary

        assert( (type_indice(0) == CON) && (type_indice(1) == CON) ) ; 

        for (int ic=0; ic<n_comp; ic++) {
            assert(cmp[ic]->check_dzpuis(4)) ;  // dzpuis=4 is assumed
        }

        const Vector& ww = longit_pot(metre) ;
        
        Tensor dww = ww.derive_con(metre) ; 
                
        for (int i=1; i<=3; i++) {
            for (int j=1; j<=3; j++) {
                dww.set(i,j).inc_dzpuis(2) ; 
            }
        }
        
        p_transverse[jp] = new Sym_tensor_trans(*mp, *triad, metre) ;
        
        for (int i=1; i<=3; i++) {
            for (int j=i; j<=3; j++) {
                p_transverse[jp]->set(i,j) = operator()(i,j) 
                    - dww(i,j) - dww(j,i) ; 
            }
        }

    }

    return *p_transverse[jp] ;
    

}




const Vector& Sym_tensor::longit_pot(const Metric& metre) const {

    set_dependance(metre) ;
    int jp = get_place_met(metre) ;
    assert ((jp>=0) && (jp<N_MET_MAX)) ;

    if (p_longit_pot[jp] == 0x0) {  // a new computation is necessary
        
        const Metric_flat* metf = dynamic_cast<const Metric_flat*>(&metre) ; 
        if (metf == 0x0) {
            cout << "Sym_tensor::longit_pot : the case of a non flat metric"
             << endl <<"  is not treated yet !" << endl ; 
            abort() ; 
        }
        
        Vector hhh = divergence(metre) ; 


        // If dpzuis is 5, it should be decreased to 4 for the Poisson equation:
        bool dzp5 = false ; 
        for (int i=1; i<=3; i++) {
            dzp5 = dzp5 || hhh(i).check_dzpuis(5) ;
        }
        if (dzp5) hhh.dec_dzpuis() ; 
                
        p_longit_pot[jp] = new Vector( hhh.poisson(double(1)) ) ; 
        
        //## Test of resolution of the vector Poisson equation:
        const Vector& vv = *(p_longit_pot[jp]) ; 

        hhh.dec_dzpuis() ; 

        Vector vtest = vv.derive_con(metre).divergence(metre) 
                        + (vv.divergence(metre)).derive_con(metre)
                        - hhh ;

        // cout << "## Sym_tensor::longit_pot : test of Poisson : \n " ; 
        // vtest.spectral_display() ;  
        cout << "## Sym_tensor::longit_pot : test of Poisson : \n " ;
        cout << 
        "  Max absolute error in each domain on the vector Poisson equation: \n" ;   
        maxabs(vtest) ; 
        
    }

    return *p_longit_pot[jp] ;
    

}

