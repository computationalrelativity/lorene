/*
 * Main code for reading output files from stationary axisymmetric rotating stars. 
 * 
 */

/*
 *   Copyright (c) 2010 Frederic Vincent
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
 * Revision 1.3  2016/12/05 16:18:25  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:53:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2010/03/22 12:44:55  f_vincent
 * Added read_star.C for reading output files computed by nrotstar.C
 *
 *
 *
 * $Header$
 */

// C headers
#include <cstdlib>
#include <cmath>
#include <cstring>

// Lorene headers
#include "star_rot.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"	    

using namespace Lorene ;

int main(){
  FILE* res=fopen("resu.d","r");
  
  Mg3d mg(res);
  
  cout << "mg= " << mg << endl;

  Map_et mp(mg,res);

  cout << "mp= " << mp << endl;

  Eos* p_eos = Eos::eos_from_file(res);

  cout << "eos=" << *p_eos << endl;

  Star_rot star(mp,*p_eos,res);

  star.equation_of_state();
  star.update_metric();
  star.hydro_euler();
  cout << "star=" << star << endl;

  Metric gg=star.get_gamma();

  cout << "metric=" << gg << endl;

  Scalar g_11=gg.cov()(1,1);
  
  cout << "g11=" << g_11.val_point(1.,M_PI,0.) << endl;

  Scalar d1g11=g_11.dsdr();

  cout << "der " << d1g11.val_point(1.,M_PI,0.) << endl;
}
