/*
 *   Copyright (c) 2003 Philippe Grandclement
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
 * Revision 1.3  2016/12/05 16:18:11  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:53:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2003/12/11 14:48:50  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

#include "proto.h"
#include "ope_elementary.h"

// Standard constructor :
namespace Lorene {
Ope_elementary::Ope_elementary (int nbr, int base, double alf, double bet) : 
  nr (nbr), base_r (base), alpha(alf), beta(bet),
  ope_mat(0x0), ope_cl (0x0), non_dege(0x0) {}

// Constructor by copy:
Ope_elementary::Ope_elementary (const Ope_elementary& so) : 
  nr (so.nr), base_r(so.base_r), alpha(so.alpha), beta(so.beta), 
  ope_mat(0x0), ope_cl (0x0), non_dege(0x0) {
  
  if (so.ope_mat != 0x0)
    ope_mat = new Matrice (*so.ope_mat) ;

  if (so.ope_cl != 0x0)
    ope_cl = new Matrice (*so.ope_cl) ;

   if (so.non_dege != 0x0)
    non_dege = new Matrice (*so.non_dege) ;
}

// Destructor :
Ope_elementary::~Ope_elementary() {
  if (ope_mat != 0x0)
    delete ope_mat ;

  if (ope_cl != 0x0)
    delete ope_cl ;
  
  if (non_dege != 0x0)
    delete non_dege ;
}

}
