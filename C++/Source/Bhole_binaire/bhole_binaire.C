/*
 *  Methods of class Bhole_binaire
 *
 *   (see file bhole.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char bhole_binaire_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.46  2001/05/07  09:11:43  phil
 * *** empty log message ***
 *
 * Revision 2.45  2001/04/26  12:04:23  phil
 * *** empty log message ***
 *
 * Revision 2.44  2001/04/02  12:16:32  phil
 * *** empty log message ***
 *
 * Revision 2.43  2001/03/30  15:52:15  phil
 * y
 *
 * Revision 2.42  2001/03/22  10:43:05  phil
 * *** empty log message ***
 *
 * Revision 2.41  2001/03/01  08:17:43  phil
 * mets les procedures ailleurs
 *
 * Revision 2.40  2001/02/28  15:26:44  phil
 * correction dans construction source_deux
 *
 * Revision 2.39  2001/02/28  14:36:03  phil
 * modif etat_zero dans solve_shift
 *
 * Revision 2.38  2001/02/28  13:22:51  phil
 * vire kk_auto
 *
 * Revision 2.37  2001/02/12  15:37:29  phil
 * ajout calcul de J infini
 *
 * Revision 2.36  2001/01/29  14:30:26  phil
 * ajout type_rotation
 *
 * Revision 2.35  2001/01/22  09:29:23  phil
 * vire trucs sur bare masse
 *
 * Revision 2.34  2001/01/10  11:03:27  phil
 * utilisation de homothetie_interne dans fixe_bare
 *
 * Revision 2.33  2001/01/10  09:52:24  phil
 * correction de fait_kk_auto
 *
 * Revision 2.32  2001/01/10  09:31:41  phil
 * ajout de fait_kk_auto
 *
 * Revision 2.31  2001/01/09  13:38:42  phil
 * ajout membre kk_auto
 *
 * Revision 2.30  2001/01/04  10:53:44  phil
 * modification des sources pour tenir compte sommation
 *
 * Revision 2.29  2000/12/21  10:49:34  phil
 * retour a la version 2.27
 *
 * Revision 2.27  2000/12/20  10:43:14  phil
 * correction calcul de la distance propre
 *
 * Revision 2.26  2000/12/20  09:09:29  phil
 * modification calcul de la masse ADM
 *
 * Revision 2.25  2000/12/18  17:43:21  phil
 * correction homothetie
 *
 * Revision 2.24  2000/12/18  16:38:19  phil
 * ajout de convergence vers une masse donnee
 *
 * Revision 2.23  2000/12/15  16:44:20  phil
 * *** empty log message ***
 *
 * Revision 2.22  2000/12/15  16:41:17  phil
 * ajout calcul de la separation
 *
 * Revision 2.21  2000/12/14  15:13:02  phil
 * corection equation sur Psi
 *
 * Revision 2.20  2000/12/14  14:12:14  phil
 * correction cl sur Psi
 *
 * Revision 2.19  2000/12/14  12:41:41  phil
 * on met les bases dans les sources
 *
 * Revision 2.18  2000/12/14  10:53:18  phil
 * *** empty log message ***
 *
 * Revision 2.17  2000/12/14  10:45:11  phil
 * ATTENTION : PASSAGE DE PHI A PSI
 *
 * Revision 2.16  2000/12/13  15:35:32  phil
 * ajout calcul de bare_masse
 *
 * Revision 2.15  2000/12/04  14:30:00  phil
 * ajout de grad_n_tot
 *
 * Revision 2.14  2000/12/01  16:13:04  phil
 * *** empty log message ***
 *
 * Revision 2.13  2000/12/01  14:27:24  phil
 * *** empty log message ***
 *
 * Revision 2.12  2000/12/01  14:18:19  phil
 * *** empty log message ***
 *
 * Revision 2.11  2000/12/01  14:16:29  phil
 * ahout verifie_cl
 *
 * Revision 2.10  2000/11/24  11:05:16  phil
 * calcul moment foireux
 *
 * Revision 2.9  2000/11/24  09:57:46  phil
 * ajout calcul masse et moment pour systeme
 *
 * Revision 2.8  2000/11/15  18:26:53  phil
 * vire recherches de minimum
 *
 * Revision 2.7  2000/11/15  13:03:40  phil
 * *** empty log message ***
 *
 * Revision 2.6  2000/11/15  13:00:16  phil
 * changements diverses
 *
 * Revision 2.5  2000/10/26  09:01:26  phil
 * *** empty log message ***
 *
 * Revision 2.4  2000/10/26  08:25:24  phil
 * *** empty log message ***
 *
 * Revision 2.3  2000/10/26  08:21:04  phil
 * ajout de verifie_shift
 *
 * Revision 2.2  2000/10/24  13:38:36  phil
 * on ne resout qu une fois quand omega est bloque
 *
 * Revision 2.1  2000/10/23  09:24:43  phil
 * rearangement constructeurs.
 *
 * Revision 2.0  2000/10/20  09:19:20  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>
#include <fstream.h>

// Lorene
#include "nbr_spx.h"
#include "map.h"
#include "valeur.h"
#include "cmp.h"
#include "tenseur.h"
#include "bhole.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

// Constucteur standard
Bhole_binaire::Bhole_binaire (Map_af& mp1, Map_af& mp2) :
	hole1(mp1), hole2(mp2), omega(0){

    holes[0] = &hole1 ;
    holes[1] = &hole2 ;
}

// Copy
Bhole_binaire::Bhole_binaire (const Bhole_binaire& source) :
	    hole1(source.hole1), hole2(source.hole2), omega(source.omega) {
    
    holes[0] = &hole1 ;
    holes[1] = &hole2 ;
    }
    
Bhole_binaire::~Bhole_binaire () {
}

//Affectation
void Bhole_binaire::operator= (const Bhole_binaire& source) {    
    hole1 = source.hole1 ;
    hole2 = source.hole2 ;
    
    omega = source.omega ;
}

//Initialise : somme de deux statiques
void Bhole_binaire::init_bhole_binaire() {
    set_omega (0) ;
    hole1.init_bhole() ;
    hole2.init_bhole() ;
    
    hole1.fait_psi_comp(hole2) ;
    hole2.fait_psi_comp(hole1) ;
    
    hole1.fait_n_comp(hole2) ;
    hole2.fait_n_comp(hole1) ;
    
    fait_decouple() ;
}


