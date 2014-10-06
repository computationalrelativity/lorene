/*
 *  Methods of template class Evolution_std
 *
 *    (see file evolution.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon & Jerome Novak
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
 * Revision 1.9  2014/10/06 15:12:51  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2005/01/11 12:48:46  f_limousin
 * Implement the function operator=(const Evolution_std<TyT>& ).
 *
 * Revision 1.7  2004/05/31 09:06:12  e_gourgoulhon
 * Added protection against self-assignement in method update.
 *
 * Revision 1.6  2004/03/26 13:31:09  j_novak
 * Definition of the macro UNDEF_STEP for non-defined time-steps.
 * Changes in the way the time derivative is calculated.
 *
 * Revision 1.5  2004/03/26 08:22:13  e_gourgoulhon
 * *** Full reorganization of class Evolution ***
 * Introduction of the notion of absoluteuniversal time steps,
 * stored in the new array 'step'.
 * The new function position(int j) makes a correspondence
 * between a universal time step j and the position in the
 * arrays step, the_time and val.
 * Only method update is now virtual.
 * Methods operator[], position, is_known, downdate belong to
 * the base class.
 *
 * Revision 1.4  2004/03/23 14:50:41  e_gourgoulhon
 * Added methods is_updated, downdate, get_jlast, get_size,
 * as well as constructors without any initial value.
 * Formatted documentation for Doxygen.
 *
 * Revision 1.3  2004/02/17 22:13:34  e_gourgoulhon
 * Suppressed declaration of global char[] evolution_C = ...
 *
 * Revision 1.2  2004/02/16 12:37:34  e_gourgoulhon
 * Added an assert in method update.
 *
 * Revision 1.1  2004/02/15 21:55:33  e_gourgoulhon
 * Introduced derived classes Evolution_full and Evolution_std.
 * Evolution is now an abstract base class.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h" 

// C headers
#include <cstdlib>
#include <cassert>



                    //-------------------------//
                    //      Constructors       //
                    //-------------------------//

                    
template<typename TyT> 
Evolution_std<TyT>::Evolution_std(const TyT& initial_value, int nstored, 
                                  int initial_j, double initial_time) 
      : Evolution<TyT>(initial_value, initial_j, initial_time, nstored)
{ }                    
                                        

template<typename TyT> 
Evolution_std<TyT>::Evolution_std(int nstored) 
      : Evolution<TyT>(nstored)
{ }                    
                                        

template<typename TyT> 
Evolution_std<TyT>::Evolution_std(const Evolution_std<TyT>& evo)
      : Evolution<TyT>(evo)
{ }


                    
                    
                    //-----------------------//
                    //      Destructor       //
                    //-----------------------//

                    
template<typename TyT> 
Evolution_std<TyT>::~Evolution_std(){ }

                    
                    
                    //-----------------------//
                    //      Mutators         //
                    //-----------------------//

                    
template<typename TyT> 
void Evolution_std<TyT>::operator=(const Evolution_std<TyT>& evo) {

    size = evo.size ;
    pos_jtop = evo.pos_jtop ;
 
    for (int j=0; j<size; j++) {
        step[j] = evo.step[j] ; 
    }
    
    for (int j=0; j<size; j++) {
        the_time[j] = evo.the_time[j] ; 
    }
    

    for (int j=0; j<size; j++) {
        if (val[j] != 0x0) {
	    delete val[j] ;
	    val[j] = 0x0 ;
	}
    }

    for (int j=0; j<size; j++) {
        if (evo.val[j] != 0x0) {	    
            val[j] = new TyT( *(evo.val[j]) ) ; 
        }
        else {
            val[j] = 0x0 ; 
        }
    }
}

template<typename TyT> 
void Evolution_std<TyT>::operator=(const Evolution<TyT>& ) {

    cerr << "void Evolution_std<TyT>::operator= : not implemented yet ! \n" ; 
    abort() ; 
 
}


                    
template<typename TyT> 
void Evolution_std<TyT>::update(const TyT& new_value, int j, double time_j) {


    if (is_known(j)) {   // Case of a time step already stored
                         //-----------------------------------
        int pos = position(j) ; 
        assert( fabs(the_time[pos] - time_j) < 1.e-14 ) ;   
        assert( val[pos] != &new_value ) ; // to avoid self assignment
        delete val[pos] ; 
        val[pos] = new TyT(new_value) ; 
    }
    else {              // Storage of a new time step
                        //---------------------------

        if ( (pos_jtop != -1) && (j < step[pos_jtop]) ) {
            cerr << 
                "Evolution_std<TyT>::update : the time step j = "
                << j << " must be in the future\n" 
                << "  of the last stored time step (" <<  step[pos_jtop] << ") !"
                << endl ; 
            abort() ; 
        }
          
        pos_jtop++ ; 
            
        if (pos_jtop == size) {  // re-organization of arrays step, the_time 
                                // and val is necessary
    
            if ( val[0] != 0x0 ) delete val[0] ; 

            for (int i=0; i<size-1; i++) {
                step[i] = step[i+1] ; 
                the_time[i] = the_time[i+1] ; 
                val[i] = val[i+1] ; 
            }
        
            pos_jtop-- ;        // pos_jtop = size-1
        
        }
        else {
            assert( pos_jtop < size ) ; 
            assert( val[pos_jtop] == 0x0 ) ; 
        }
    
        step[pos_jtop] = j ;
        the_time[pos_jtop] = time_j ; 
        val[pos_jtop] = new TyT( new_value ) ; 
    }

}                   
                    

                    
