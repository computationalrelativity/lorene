/*
 * Method of regularization of the source of Poisson equation
 *
 * (see file cmp.h for documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.3  2016/12/05 16:17:49  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:52:48  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  2000/09/07  15:29:14  keisuke
 * Add a new argument Cmp& uu.
 *
 * Revision 2.8  2000/09/04  15:53:08  keisuke
 * Insert the polar and azimuthal parts of duu_div into Map_af::poisson_regular.
 *
 * Revision 2.7  2000/09/04  13:19:34  keisuke
 * Suppress the version without parameters.
 * Express the code by using Map_af::poisson_regular.
 *
 * Revision 2.6  2000/08/31  15:58:53  keisuke
 * Modify the polar and azimuthal derivatives of uu_div.
 *
 * Revision 2.5  2000/08/30  16:01:52  keisuke
 * Change the constant "R" into "mp_radial->dxdr".
 *
 * Revision 2.4  2000/08/29  13:52:19  keisuke
 * Add the polar and azimuthal derivatives of the diverging potential.
 * Modify the argumants.
 *
 * Revision 2.3  2000/08/29  08:39:37  keisuke
 * Minor change.
 *
 * Revision 2.2  2000/08/28  15:57:39  keisuke
 * *** empty log message ***
 *
 * Revision 2.1  2000/08/28  15:54:18  keisuke
 * Add "int nzet" in the argumant.
 *
 * Revision 2.0  2000/08/25  08:43:38  keisuke
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header Lorene
#include "cmp.h"
#include "tenseur.h"
#include "map.h"
#include "param.h"

namespace Lorene {

//******************************************************************

void Cmp::poisson_regular(int k_div, int nzet, double unsgam1, Param& par,
			  Cmp& uu, Cmp& uu_regu, Cmp& uu_div,
			  Tenseur& duu_div,
			  Cmp& source_regu, Cmp& source_div) const {

    mp->poisson_regular(*this, k_div, nzet, unsgam1, par,
			uu, uu_regu, uu_div, duu_div,
			source_regu, source_div) ;


}


}
