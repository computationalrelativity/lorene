/*
 *  Definition of methods for the class Metric_flat.
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (for previous class Metrique)
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

char metric_flat_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/11 14:40:39  e_gourgoulhon
 * Suppressed declaration of unusued argument in method operator=.
 *
 * Revision 1.1  2003/10/06 15:30:33  j_novak
 * Defined methods for flat metric.
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include<stdlib.h>

// Lorene headers
#include "metric.h"
#include "utilitaires.h"

Metric_flat::Metric_flat(const Map& mpi, const Base_vect& triadi) :
  Metric(mpi), triad(&triadi) {
  
}
  

Metric_flat::Metric_flat(const Metric_flat& meti) : Metric(meti), 
triad(meti.triad) {


}

Metric_flat::Metric_flat(const Map& mpi, FILE* fd) : Metric(mpi, fd) {

  cout << "Metric_flat::Metric_flat(FILE*) : not implemented yet!" << endl ;

  abort() ;
}


Metric_flat::~Metric_flat() {

}

void Metric_flat::operator=(const Metric_flat& meti) {

  Metric::operator=(meti) ;

  triad = meti.triad ;
}

void Metric_flat::operator=(const Sym_tensor& ) {
  
  cout << "Metric_flat::operator=(const Sym_tensor& ) :" << '\n' ;
  cout << "Error: a flat metric should not be specified" << '\n' ;
  cout << "by a symmetric tensor!" << endl ;

  abort() ;
}
  

void Metric_flat::fait_cov() const {
  
  assert( p_met_cov == 0x0 ) ;

  p_met_cov = new Sym_tensor(*mp, COV, *triad) ;
  p_met_cov->set_etat_zero() ;
  for (int i=1; i<=3; i++) 
    p_met_cov->set(i,i) = 1 ;

}

void Metric_flat::fait_con() const {

  assert( p_met_con == 0x0 ) ;

  p_met_con = new Sym_tensor(*mp, CON, *triad) ;
  p_met_con->set_etat_zero() ;
  for (int i=1; i<=3; i++) 
    p_met_con->set(i,i) = 1 ;

}

void Metric_flat::fait_connection() const {

  assert( p_connect == 0x0 ) ;

  const Base_vect_spher* bvs =
    dynamic_cast<const Base_vect_spher*>(triad) ;
  
  const Base_vect_cart* bvc =
    dynamic_cast<const Base_vect_cart*>(triad) ;

  if (bvs != 0x0) {
    assert (bvc == 0x0) ;
    p_connect = new Connection_fspher(*mp, *bvs) ;
  }
  else {
    assert(bvc != 0x0) ;
    p_connect = new Connection_fcart(*mp, *bvc) ;
  }

}

void Metric_flat::fait_ricci_scal() const {

  assert( p_ricci_scal == 0x0 ) ;

  p_ricci_scal = new Scalar(*mp) ;
  p_ricci_scal->set_etat_zero() ;
}

void Metric_flat::fait_determinant() const {

  assert( p_determinant == 0x0 ) ;

  p_determinant = new Scalar(*mp) ;
  *p_determinant = 1 ;
}
 
void Metric_flat::sauve(FILE* ) const {

  cout << "Metric_flat::sauve(FILE*) : not implemented yet!" << endl ;

  abort() ; //## What to do with the connection, triad ... ?

}

ostream& Metric_flat::operator>>(ostream& ost) const {

  ost << '\n' ;

  ost << "Flat metric in an orthonormal triad" << '\n' ;
  ost << "-----------------------------------" << '\n' ;
  ost << '\n' ;

  ost << *triad ;

  ost << endl ;
  return ost ;
}

