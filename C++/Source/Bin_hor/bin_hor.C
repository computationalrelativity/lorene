/*
 *  Methods of class Bin_hor
 *
 *   (see file bin_hor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose Luis Jaramillo
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char bin_hor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/12/29 16:11:02  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "nbr_spx.h"
#include "tenseur.h"
#include "tensor.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

// Constucteur standard
Bin_hor::Bin_hor (Map_af& mp1, Map_af& mp2, int depth_in) :
	hole1(mp1, depth_in), hole2(mp2, depth_in), omega(0){

    holes[0] = &hole1 ;
    holes[1] = &hole2 ;
}

// Copy
Bin_hor::Bin_hor (const Bin_hor& source) :
	    hole1(source.hole1), hole2(source.hole2), omega(source.omega) {
    
    holes[0] = &hole1 ;
    holes[1] = &hole2 ;
    }
    
Bin_hor::~Bin_hor () {
}

//Affectation
void Bin_hor::operator= (const Bin_hor& source) {    
    hole1 = source.hole1 ;
    hole2 = source.hole2 ;
    
    omega = source.omega ;
}

//Initialisation : Sum of two static BH
void Bin_hor::init_bin_hor() {
    set_omega (0) ;
    hole1.init_bhole() ;
    hole2.init_bhole() ;
    
    hole1.psi_comp(hole2) ;
    hole2.psi_comp(hole1) ;
    
    hole1.n_comp(hole2) ;
    hole2.n_comp(hole1) ;
    
    decouple() ;
}


