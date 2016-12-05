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

 

/*
 * $Id$
 * $Log$
 * Revision 1.8  2016/12/05 16:17:59  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:53:07  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:14  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2005/01/12 16:03:32  j_novak
 * We initialize the covariant representation of a flat metric, to avoid any
 * possible problem when building a general metric from a flat one.
 *
 * Revision 1.4  2004/02/19 10:57:18  e_gourgoulhon
 * Methods cov() and con(): set the spectral bases to the standard ones
 * after the initialization of p_met_cov and p_met_con.
 *
 * Revision 1.3  2003/12/30 23:07:29  e_gourgoulhon
 * Suppression of virtual methods fait_* : the actual computations
 *      are now performed via the virtual methods con(), cov(), connect(),
 *      ricci(), ricci_scal(), determinant().
 *
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
#include <cstdlib>

// Lorene headers
#include "metric.h"
#include "utilitaires.h"

                    //-----------------//
                    //  Constructors   //
                    //-----------------//

namespace Lorene {
Metric_flat::Metric_flat(const Map& mpi, const Base_vect& triadi) :
  Metric(mpi), triad(&triadi) {

    cov() ; //## to avoid problems when initializing general metrics
            //## to flat ones. Might be improved??
  
}
  

Metric_flat::Metric_flat(const Metric_flat& meti) : Metric(meti), 
triad(meti.triad) {


}

Metric_flat::Metric_flat(const Map& mpi, FILE* fd) : Metric(mpi, fd) {

  cout << "Metric_flat::Metric_flat(FILE*) : not implemented yet!" << endl ;

  abort() ;
}


                    //---------------//
                    //  Destructor   //
                    //---------------//

Metric_flat::~Metric_flat() {

}


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

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
  

            //----------------//
            //   Accessors    //
            //----------------//


const Sym_tensor& Metric_flat::cov() const {
  
    if (p_met_cov == 0x0) {   // a new computation is necessary

        p_met_cov = new Sym_tensor(*mp, COV, *triad) ;

        p_met_cov->set(1,1) = 1 ; 
        p_met_cov->set(1,2) = 0 ; 
        p_met_cov->set(1,3) = 0 ; 
        p_met_cov->set(2,2) = 1 ; 
        p_met_cov->set(2,3) = 0 ; 
        p_met_cov->set(3,3) = 1 ; 
        
        p_met_cov->std_spectral_base() ; 
    }

    return *p_met_cov ; 
}

const Sym_tensor& Metric_flat::con() const {
  
    if (p_met_con == 0x0) {   // a new computation is necessary
  
        p_met_con = new Sym_tensor(*mp, CON, *triad) ;

        p_met_con->set(1,1) = 1 ; 
        p_met_con->set(1,2) = 0 ; 
        p_met_con->set(1,3) = 0 ; 
        p_met_con->set(2,2) = 1 ; 
        p_met_con->set(2,3) = 0 ; 
        p_met_con->set(3,3) = 1 ; 

        p_met_con->std_spectral_base() ; 

    }


    return *p_met_con ; 
}


const Connection& Metric_flat::connect() const {

    if (p_connect == 0x0) {   // a new computation is necessary
    
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

    return *p_connect ; 

}



const Scalar& Metric_flat::ricci_scal() const {

    if (p_ricci_scal == 0x0) {   // a new computation is necessary

        p_ricci_scal = new Scalar(*mp) ;
        p_ricci_scal->set_etat_zero() ;
    }

    return *p_ricci_scal  ; 

}


const Scalar& Metric_flat::determinant() const {

    if (p_determinant == 0x0) {   // a new computation is necessary

        p_determinant = new Scalar(*mp) ;
        p_determinant->set_etat_one() ;
    }

    return *p_determinant ; 
}


                //---------//
                // Outputs //
                //---------//

 
void Metric_flat::sauve(FILE* ) const {

  cout << "Metric_flat::sauve(FILE*) : not implemented yet!" << endl ;

  abort() ; //## What to do with the  triad ... ?

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

}
