/*
 *  Definition of methods for the class Metric.
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

char metric_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/10/02 15:45:50  j_novak
 * New class Metric
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

Metric::Metric(const Sym_tensor& symti) : mp(&symti.get_mp()), p_connect(0x0),
					  p_met_cov(0x0), p_met_con(0x0) {
  
  int type_index = symti.get_index_type(0) ;
  assert (symti.get_index_type(1) == type_index) ;

  if (type_index == COV) {
    p_met_cov = new Sym_tensor(symti) ;
  }
  else {
    assert(type_index == CON) ;
    p_met_con = new Sym_tensor(symti) ;
  }

  set_der_0x0() ;

}


Metric::Metric(const Metric& meti) : mp(meti.mp), p_connect(0x0), p_met_cov(0x0),
				     p_met_con(0x0) {

  if (meti.p_met_cov != 0x0) p_met_cov = new Sym_tensor(*meti.p_met_cov) ;

  if (meti.p_met_con != 0x0) p_met_con = new Sym_tensor(*meti.p_met_con) ;

  set_der_0x0() ;

}

Metric::Metric(const Map& mpi, FILE* ) : mp(&mpi), p_connect(0x0), 
					 p_met_cov(0x0), p_met_con(0x0) {

  cout << "Metric::Metric(FILE*) : not implemented yet!" << endl ;

  abort() ;
}

Metric::Metric(const Map& mpi) : mp(&mpi), p_connect(0x0), 
					 p_met_cov(0x0), p_met_con(0x0) {

}

Metric::~Metric() {

  if (p_connect != 0x0) delete p_connect ;

  if (p_met_cov != 0x0) delete p_met_cov ;

  if (p_met_con != 0x0) delete p_met_con ;

  del_deriv() ;

}

void Metric::del_deriv() const {
  
  if (p_ricci_scal != 0x0) delete p_ricci_scal ;

  if (p_determinant != 0x0) delete p_determinant ;

  set_der_0x0() ;

}

void Metric::set_der_0x0() const {

  p_ricci_scal = 0x0 ;

  p_determinant = 0x0 ;

}

void Metric::operator=(const Metric& meti) {

  assert( mp == meti.mp) ;

  if (p_connect != 0x0) delete p_connect ;

  if (meti.p_connect != 0x0) {
    p_connect = new Connection(*meti.p_connect) ;
  }
  else {
    p_connect = 0x0 ;
  }

  if (p_met_cov != 0x0) delete p_met_cov ;

  if (meti.p_met_cov != 0x0) {
    p_met_cov = new Sym_tensor(*meti.p_met_cov) ;
  }
  else {
    p_met_cov = 0x0 ;
  }

  if (p_met_con != 0x0) delete p_met_con ;

  if (meti.p_met_con != 0x0) {
    p_met_con = new Sym_tensor(*meti.p_met_con) ;
  }
  else {
    p_met_con = 0x0 ;
  }

  del_deriv() ;

}







