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

 

/*
 * $Id$
 * $Log$
 * Revision 1.11  2016/12/05 16:18:19  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.10  2014/10/13 08:53:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2010/10/20 07:58:09  j_novak
 * Better implementation of the explicit time-integration. Not fully-tested yet.
 *
 * Revision 1.7  2007/11/06 12:10:56  j_novak
 * Corrected some mistakes.
 *
 * Revision 1.6  2007/11/06 11:53:34  j_novak
 * Use of contravariant version of the 3+1 equations.
 *
 * Revision 1.5  2007/06/28 14:40:36  j_novak
 * Dynamical check: the fields in the last domain are set to zero to avoid dzpuis problems
 *
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
#include <cstdlib>
#include <cassert>

// Lorene headers
#include "time_slice.h"
#include "unites.h"

namespace Lorene {
Tbl Time_slice::check_hamiltonian_constraint(const Scalar* energy_density,
					     ostream& ost, bool verb) const {
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

  if (verb)  ost << endl ;
  const char* comment = 0x0 ;
  if (verb) comment = "Check of the Hamiltonian constraint" ;
  Tbl resu = maxabs(field - (4*qpig) * (*matter), comment, ost,verb ) ;
  if (verb) ost << endl ;

  if (vacuum) delete matter ;

  return resu ;

}
    
Tbl Time_slice::check_momentum_constraint(const Vector* momentum_density, 
					  ostream& ost, bool verb) const {
  using namespace Unites ;
  
  bool vacuum = ( momentum_density == 0x0 ) ;
  
  Vector field = k_uu().divergence(gam()) - trk().derive_con(gam()) ;


  const Vector* matter ;
  if (vacuum) 
    matter = new Vector (0*field) ;
  else {
    assert (momentum_density->get_index_type(0) == CON) ;
    matter = momentum_density ;
  }
  
  if (verb)  ost << endl ;
  const char* comment = 0x0 ;
  if (verb) comment = "Check of the momentum constraint" ;
  Tbl resu = maxabs(field - (2*qpig) * (*matter), comment, ost, verb ) ;

  if (verb)  ost << endl ;

  if (vacuum) delete matter ;

  return resu ;

}

Tbl Time_slice::check_dynamical_equations(const Sym_tensor* strain_tensor,
					  const Scalar* energy_density,
					  ostream& ost, bool verb) const {

  using namespace Unites ;

  bool vacuum = ( ( strain_tensor == 0x0 ) && ( energy_density == 0x0 ) ) ;

  Sym_tensor dyn_lhs = (k_dd_evol.time_derive(jtime, scheme_order)).up_down(gam()) ;
  int nz = dyn_lhs.get_mp().get_mg()->get_nzone() ;
  dyn_lhs.annule_domain(nz-1) ;
  dyn_lhs = dyn_lhs - k_dd().derive_lie(beta()).up_down(gam())  ;
  dyn_lhs.annule_domain(nz-1) ;
  
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
    matter = new Sym_tensor( (sij->trace(gam()) - *ener)*gam().con() 
			      - 2*(*sij) ) ;
    if (new_e) delete ener ;
    if (new_s) delete sij ;
  }

  Sym_tensor dyn_rhs = nn()*( (gam().ricci().up_down(gam()) + trk()*k_uu() 
						     + qpig * (*matter))) ;
  dyn_rhs.annule_domain(nz-1) ;
  dyn_rhs = dyn_rhs - 2*nn()*contract(k_uu(), 1, k_dd().up(0, gam()), 1)  ;
  dyn_rhs.annule_domain(nz-1) ;
  dyn_rhs = dyn_rhs - nn().derive_con(gam()).derive_con(gam()) ;
  dyn_rhs.annule_domain(nz-1) ;
    
  if (verb) ost << endl ;
  const char* comment = 0x0 ;
  if (verb) comment = "Check of the dynamical equations" ;
  Tbl resu_tmp = maxabs(dyn_lhs - dyn_rhs, comment, ost, verb) ;

  if (verb) ost << endl ;

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
}
