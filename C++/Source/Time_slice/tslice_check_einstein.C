/*
 *  Methods of class Time_slice to check Einstein equation solutions.
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon, Jose Luis Jaramillo & Jerome Novak
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

char tslice_check_einstein_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.3  2004/04/07 07:58:21  e_gourgoulhon
 * Constructor as Minkowski slice: added call to std_spectral_base()
 * after setting the lapse to 1.
 *
 * Revision 1.2  2004/04/05 12:38:45  j_novak
 * Minor modifs to prevent some warnings.
 *
 * Revision 1.1  2004/04/05 11:54:20  j_novak
 * First operational (but not tested!) version of checks of Eintein equations.
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "unites.h"

Tbl Time_slice::check_hamiltonian_constraint(const Scalar* energy_density,
					     ostream& ost) const {
  using namespace Unites ;

  bool vacuum = ( energy_density == 0x0 ) ;
  
  Scalar field = trk() * trk() - contract( k_uu(), 0, 1, k_dd(), 0, 1 ) ;
  
  field.dec_dzpuis() ;  // dzpuis: 4 -> 3 
    
  field += gam().ricci_scal() ;

  const Scalar* matter ;
  if (vacuum) 
    matter = new Scalar (0*field) ;
  else
    matter = energy_density ;

  cout << endl ;

  Tbl resu = diffrelmax(field, (4*qpig) * (*matter), 
			"Check of the Hamiltonian constraint", ost ) ;
  cout << endl ;

  if (vacuum) delete matter ;

  return resu ;

}
    
Tbl Time_slice::check_momentum_constraint(const Vector* momentum_density, 
					     ostream& ost) const {
  using namespace Unites ;
  
  bool vacuum = ( momentum_density == 0x0 ) ;
  
  Vector field = - trk().derive_cov(gam()) ;

  if (k_uu_evol.is_known(jtime)) 
    field += k_uu().down(0, gam()).divergence(gam()) ;
  else
    field += k_dd().up(1, gam()).divergence(gam()) ;

  const Vector* matter ;
  if (vacuum) 
    matter = new Vector (0*field) ;
  else {
    assert (momentum_density->get_index_type(0) == COV) ;
    matter = momentum_density ;
  }
  
  cout << endl ;

  Tbl resu = diffrelmax(field, (2*qpig) * (*matter), 
			"Check of the momentum constraint", ost ) ;

  cout << endl ;

  if (vacuum) delete matter ;

  return resu ;

}

Tbl Time_slice::check_dynamical_equations(const Sym_tensor* strain_tensor,
					  const Scalar* energy_density,
					  ostream& ost) const {

  using namespace Unites ;

  bool vacuum = ( ( strain_tensor == 0x0 ) && ( energy_density == 0x0 ) ) ;

  Sym_tensor dyn_lhs = k_dd_evol.time_derive(jtime, scheme_order) 
    - k_dd().derive_lie(beta()) ;
  
  const Sym_tensor* matter ;
  if (vacuum) 
    matter = new Sym_tensor(0 * dyn_lhs) ;
  else {
    const Scalar* ener = 0x0 ;
    const Sym_tensor* sij = 0x0 ;
    bool new_e = false ;
    bool new_s = false ;
    if (energy_density == 0x0) {
      ener = new Scalar ( 0 * nn() ) ;
      new_e = true ;
    }
    else {
      ener = energy_density ;
      if (strain_tensor == 0x0) {
	sij = new Sym_tensor(0 * dyn_lhs) ;
	new_s = true ;
      }
      else
	sij = strain_tensor ;
    }
    matter = new Sym_tensor( (sij->trace(gam()) - *ener)*gam().cov() 
			      - 2*(*sij) ) ;
    if (new_e) delete ener ;
    if (new_s) delete sij ;
  }

  Sym_tensor dyn_rhs = nn()*( (gam().ricci() + trk()*k_dd() + qpig * (*matter))
			      - 2*contract(k_dd(), 1, k_dd().up(0, gam()), 0) )
    - nn().derive_cov(gam()).derive_cov(gam()) ;
    
  cout << endl ;

  Tbl resu_tmp = diffrelmax(dyn_lhs, dyn_rhs, 
			"Check of the dynamical equations", ost ) ;

  cout << endl ;

  if (vacuum) delete matter ;

  Tbl resu(4, 4, resu_tmp.get_dim(0)) ;
  resu.annule_hard() ; //## To initialize resu(0,0,...) which 
                       // does not correspond to a valid component

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      Itbl idx(2) ;
      idx.set(0) = i ;
      idx.set(1) = j ;
      int pos = dyn_lhs.position(idx) ;
      assert (dyn_rhs.position(idx) == pos) ;

      for (int lz=0; lz<resu_tmp.get_dim(0); lz++)
	resu.set(i,j,lz) = resu_tmp(pos, lz) ;
    }
  }

  return resu ;

}
