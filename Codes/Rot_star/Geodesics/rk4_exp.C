
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

char _C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
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
// Lorene headers
#include "etoile.h"

#include "nrutil.h" 
#include "main.h"


void rk4(double y[], int n, double x, double h, double yout[]
	 , const Etoile_rot& star
	 , void (*derivs)(double, double [], double [] , const Etoile_rot&)) 
     /*Given values for the variables y[1..n] and their derivatives dydx[1..n] 
known at x, use the fourth-order Runge-Kutta method to advance the solution
 over an interval h and return the incremented variables as yout[1..n], which
 need not be a distinct array from y. The user supplies the routine 
 derivs(x,y,dydx), which returns derivatives dydx at x.*/
{
  int i; 
  double xh,hh,h6,*dym,*dyt,*yt,*dydx;
  dym=dvector(1,n);
  dyt=dvector(1,n);
  yt=dvector(1,n);
  dydx=dvector(1,n);
  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;
  #ifdef DEBUG
  printf("Entrée dans rk4 ...  ");
  #endif

  (*derivs)(x,y,dydx,star);
  /* printf("expo: %f %f %f\n",x,y[1],dydx[1]); */
  /* First step */
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];

  /* Second step */
  (*derivs)(xh,yt,dyt,star);
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];

  /* Third step */
  (*derivs)(xh,yt,dym,star);
  for (i=1;i<=n;i++)
    { 
      yt[i]=y[i]+h*dym[i];
      dym[i] += dyt[i];
    }

  /* Fourth step */
  (*derivs)(x+h,yt,dyt,star);
  for (i=1;i<=n;i++) yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

  /* Free memory */
  free_dvector(yt,1,n);
  free_dvector(dyt,1,n);
  free_dvector(dym,1,n);
  free_dvector(dydx,1,n);
  #ifdef DEBUG
  printf("Fin de rk4\n");
  #endif
}
