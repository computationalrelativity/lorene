/*
 *  Definition of Lorene template class Evolution
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


#ifndef __EVOLUTION_H_ 
#define __EVOLUTION_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/02/13 15:53:20  e_gourgoulhon
 * New (template) class for time evolution.
 *
 *
 *
 *
 * $Header$
 *
 */

/** Time evolution. 
 * 
 * Time evolution is managed through the template class {\tt Evolution}.
 *
 */
template<typename TyT> class Evolution {

        
    // Data:
    // -----
    
    protected: 
        /// Maximum number of stored time steps
        int size ; 
        
        /** Determines whether {\tt size} is maintained fixed
         *  during the object life
         */
        const bool fixed_size ; 
        
        /// Array of pointers onto the values (size {\tt size}) 
        TyT** val ; 
        
        /// Array of successive time steps
        double* the_time ; 
      
        /// Index of last updated time step for the object
        int jlast ; 
      
        
    // Constructors - Destructor
    // -------------------------
    public:
        /** Standard constructor for rescalable storage.
         *  The size of storage is set to {\tt global_size}
         *
         */
	    Evolution(const TyT& initial_value, double initial_time) ;			
	
        /// Standard constructor for fixed storage 
	    Evolution(const TyT& initial_value, double initial_time,
                int nstored) ;			
	
        Evolution(const Evolution<TyT>& ) ;		/// Copy constructor

	    virtual ~Evolution() ;			/// Destructor
 
    // Mutators 
    // --------
        /// Update
        void update(const TyT& new_value, double new_time) ; 
    
        /// Assignement
        void operator=(const Evolution<TyT>& ) ;

    
    // Accessors
    // ---------
        /// Returns value at local time step j
        const TyT& operator[](int j) const ;

        /// Returns value at time t
        TyT operator()(double t) const ;

    // Outputs
    // -------
    
    

};

#include "Template/evolution.C"

#endif

