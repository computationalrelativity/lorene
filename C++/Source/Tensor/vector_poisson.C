/*
 *  Methods for solving vector Poisson equation.
 *
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

char vector_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/22 13:08:06  j_novak
 * Better handling of dzpuis flags
 *
 * Revision 1.1  2003/10/20 15:15:42  j_novak
 * New method Vector::poisson().
 *
 *
 * $Headers: $
 *
 */

//C headers
#include "math.h"

// Lorene headers
#include "metric.h"

Vector Vector::poisson(const double lambda) const {
 
  for (int i=0; i<3; i++)
    assert(cmp[i]->check_dzpuis(4)) ;

  Metric_flat met_local(*mp, *triad) ;

  Scalar poten(*mp) ;
  if (fabs(lambda+1) < 1.e-6)
    poten.set_etat_zero() ;
  else {
    Scalar tmp = potential(met_local) / (lambda + 1) ;
    tmp.inc2_dzpuis() ;
    poten = tmp.poisson() ;
  }

  Vector grad = poten.derive_con(met_local) ;
  grad.dec2_dzpuis() ;

  return ( div_free(met_local).poisson() + grad) ;
    
 
}

// Version with parameters
// -----------------------

void Vector::poisson(const double, Param&, Scalar& ) const {

  cout << "Not ready yet!" << endl ;


}
