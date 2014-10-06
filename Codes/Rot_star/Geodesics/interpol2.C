
/*
 *   Copyright (c) 2003 CHABBERT Jean-Philippe
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

char interpol2_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:12:51  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2003/02/07 17:31:52  jp_chabbert
 * First version with rotstar input data
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers

// C headers
#include <cmath>

// Lorene headers
#include "etoile.h"
#include "main.h"
  
point interpol2( double r0, double theta, const Etoile_rot& star )
{
      double r=r0/10; // passage en unite Lorene
      point p;
      double nu=star.get_logn()().val_point(r,theta,0.);
      double dzeta=star.get_dzeta()().val_point(r,theta,0.);
      double b=star.get_bbb()().val_point(r,theta,0.);
      p.omega=star.get_nphi()().val_point(r,theta,0.)/10.;
      p.alpha=dzeta-nu;
      double beta=log(b);
      p.gamma=nu+beta;
      p.rho=nu-beta;
      
      double dnudr=star.get_logn()().dsdr().val_point(r,theta,0.);
      double dnudt=r*star.get_logn()().srdsdt().val_point(r,theta,0.);
      double dbdr=star.get_bbb()().dsdr().val_point(r,theta,0.);
      double dbdt=r*star.get_bbb()().srdsdt().val_point(r,theta,0.);
      double ddzetadr=star.get_dzeta()().dsdr().val_point(r,theta,0.);
      double ddzetadt=r*star.get_dzeta()().srdsdt().val_point(r,theta,0.);
      double dbetadr=dbdr/b;
      double dbetadt=dbdt/b;

      p.drhodr=(dnudr-dbetadr)/10.;
      p.drhodtheta=(dnudt-dbetadt);
      p.dgammadr=(dnudr+dbetadr)/10.;
      p.dgammadtheta=(dnudt+dbetadt);
      p.dalphadr=(ddzetadr-dnudr)/10.;
      p.dalphadtheta=(ddzetadt-dnudt);
      p.domegadr=star.get_nphi()().dsdr().val_point(r,theta,0.)/100.;
      p.domegadtheta=r*star.get_nphi()().srdsdt().val_point(r,theta,0.)/10.;

      //     cout << p.rho << endl;
      //cout << p.gamma << endl;
      //cout << p.alpha << endl;
      //cout << p.omega << endl;


      return p;
	}
