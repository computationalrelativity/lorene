
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

char geodesics_C[] = "$Header$" ;

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
#include <stdio.h>
#include <math.h>

// Lorene headers
#include "etoile.h"

#include "main.h"

void geodesics( double x, double *y, double *dydx , const Etoile_rot& star )
     /* Null geodesics equations in RNS metric */
     /* re = rayon de l'étoile                 */
{
  /* coordonneés et leurs dérivées */
  double t,r,theta,phi,tp,rp,thetap,phip,tpp,rpp,thetapp,phipp;
  /* Potentiels */
  double alpha,rho,gamma,omega;
  /* Dérivées des potentiels */
  double dalphadr,drhodr,dgammadr,domegadr;
  double dalphadtheta,drhodtheta,dgammadtheta,domegadtheta;
  point p;

  #ifdef DEBUG
  printf("Entrée dans geodesics ... ");
  #endif

  /* noms de variables lisibles */ 
  t     = y[1];
  r     = y[2];
  theta = y[3];
  phi   = y[4];
  tp    = y[5];
  rp    = y[6];
  thetap= y[7];
  phip  = y[8];

  /* Pour se ramener à système du premier ordre */
  dydx[1] = tp;
  dydx[2] = rp;
  dydx[3] = thetap;
  dydx[4] = phip;

  /* chargement des coefs de la metrique */
  p=interpol2(r,theta,star);

  rho=p.rho;
  gamma=p.gamma;
  alpha=p.alpha;
  omega=p.omega;

  drhodr=p.drhodr;
  dgammadr=p.dgammadr;
  dalphadr=p.dalphadr;
  domegadr=p.domegadr;

  drhodtheta=p.drhodtheta;
  dgammadtheta=p.dgammadtheta;
  dalphadtheta=p.dalphadtheta;
  domegadtheta=p.domegadtheta;


  /* Equation en t     */
  tpp=-dgammadr-drhodr-dgammadtheta-drhodtheta-exp(-2.*rho)*r*r*sin(theta)*sin(theta)*(-omega*domegadr*tp*rp-omega*domegadtheta*tp*thetap+domegadr*rp*thetap+domegadtheta*thetap*phip);

  /* Equation en r     */
  rpp=0.5*(-dgammadr*exp(gamma+rho)-drhodr*exp(gamma+rho)+2*omega*exp(gamma-rho)*r*r*sin(theta)*sin(theta)*domegadr+omega*omega*exp(gamma-rho)*r*r*sin(theta)*sin(theta)*(dgammadr-drhodr)+2*omega*omega*exp(gamma-rho)*r*sin(theta)*sin(theta))*tp*tp*exp(-2.*alpha)-exp(gamma-rho-2.*alpha)*r*sin(theta)*sin(theta)*tp*phip*(r*domegadr+omega*r*dgammadr-omega*r*drhodr+2*omega)-dalphadr*rp*rp-2.*dalphadtheta*rp*thetap+r*(1.+r*dalphadr)*thetap*thetap+0.5*exp(gamma-rho-2.*alpha)*r*sin(theta)*sin(theta)*(r*dgammadr-r*drhodr+2.)*phip*phip;

  /* Equation en theta */
  thetapp=0.5*(-dgammadtheta*exp(gamma+rho)-drhodtheta*exp(gamma+rho)+2*omega*exp(gamma-rho)*r*r*sin(theta)*sin(theta)*domegadtheta+omega*omega*exp(gamma-rho)*r*r*sin(theta)*sin(theta)*(dgammadtheta-drhodtheta)+2*omega*omega*exp(gamma-rho)*r*r*sin(theta)*cos(theta))*tp*tp*exp(-2.*alpha)*1/(r*r)-exp(gamma-rho-2.*alpha)*sin(theta)*tp*phip*(domegadtheta*sin(theta)+omega*sin(theta)*dgammadtheta-omega*sin(theta)*drhodtheta+2.*omega*cos(theta))+dalphadtheta*rp*rp/(r*r)-2./r*(1.+r*dalphadr)*rp*thetap-dalphadtheta*thetap*thetap+0.5*exp(gamma-rho-2*alpha)*sin(theta)*(sin(theta)*dgammadtheta-sin(theta)*drhodtheta+2.*cos(theta))*phip*phip;

  /* Equation en phi   */
  phipp=(-2*omega*r*exp(gamma+rho)*drhodr+omega*omega*r*r*r*exp(gamma-rho)*sin(theta)*sin(theta)*domegadr+exp(gamma+rho)*domegadr*r+2*omega*exp(gamma+rho))*tp*rp/(r*exp(gamma+rho))+(-2*omega*sin(theta)*exp(gamma+rho)*drhodtheta+omega*omega*r*r*exp(gamma-rho)*sin(theta)*sin(theta)*sin(theta)*domegadtheta+exp(gamma+rho)*domegadtheta*sin(theta)+2*omega*cos(theta)*exp(gamma+rho))*tp*thetap/(sin(theta)*exp(gamma+rho))-(omega*r*r*r*sin(theta)*sin(theta)*exp(gamma-rho)*domegadr+exp(gamma+rho)*(2+r*dgammadr-r*drhodr))*rp*phip/(r*exp(gamma+rho))-(omega*r*r*sin(theta)*sin(theta)*sin(theta)*exp(gamma-rho)*domegadtheta+exp(gamma+rho)*(2*cos(theta)+sin(theta)*dgammadtheta-sin(theta)*drhodtheta))*thetap*phip/(sin(theta)*exp(gamma+rho));

  dydx[5] = tpp;
  dydx[6] = rpp;
  dydx[7] = thetapp;
  dydx[8] = phipp;

  #ifdef DEBUG
  printf("Fin de geodesics \n");
  #endif


}
