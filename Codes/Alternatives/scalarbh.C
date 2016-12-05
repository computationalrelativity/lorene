/*
 *  Code for reading scalar-hair BH from a file
 *  Scalar-hair BH are the KBHsSH described in Herdeiro & Radu (2015)
 *  CQG, 32, 144001 
 */

/*
 *   Copyright (c) 2015 Frederic Vincent, Eric Gourgoulhon
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
 * Revision 1.6  2016/12/05 16:18:22  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2016/05/10 12:57:18  f_vincent
 * scalarBH: added the computation of quantities needed to define an accretion torus in any BS or scalar BH spacetime
 *
 * Revision 1.4  2015/12/15 06:47:11  f_vincent
 * Code cleaning in scalarBH
 *
 * Revision 1.3  2015/11/09 16:01:42  f_vincent
 * Updated scalarBH code
 *
 * Revision 1.2  2015/11/05 17:32:11  f_vincent
 * Updated code for class scalarBH.
 *
 * Revision 1.1  2015/10/22 09:20:02  f_vincent
 * New code scalarbh
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>
using namespace std ;

// Lorene headers
#include "compobj.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "proto.h"
#include "graphique.h"

using namespace Lorene ;

int main() {

  // Parameters of the computation
  // -----------------------------
    
  ifstream fpar("par_scalarbh.d") ;
  if ( !fpar.good() ) {
    cerr << "Problem in opening the file par_scalarbh.d ! " << endl ;
    abort() ;
  }
    
  char file_name[256] ; 
  fpar.getline(file_name, 256) ;
  cout << "File to be read: " << file_name << endl ;

  ifstream file(file_name) ; 
  if ( !file.good() ) {
    cerr << "Problem in opening the file " << file_name << endl ;
    abort() ;
  }
  double rHor ;
  int nrfile, nthetafile;
  file >> nrfile >> nthetafile ;
  file >> rHor ; 
  cout << "Event horizon is present at r= " << rHor << endl;

  double mass ; // M_{object,ADM} in units of m_planck^2 / m_boson
  fpar >> mass ; fpar.ignore(1000,'\n') ;
    
  int graphic_out ; // flag for graphical outputs
  fpar >> graphic_out ; fpar.ignore(1000,'\n') ; 

  int compute_qty ; // flag for computing useful doughnut quantities
  fpar >> compute_qty ; fpar.ignore(1000,'\n') ; 

  double ell ; // value of the constant angular momentum l = -u_phi/u_t
  fpar >> ell ; fpar.ignore(1000,'\n') ; 

  int nr ; // Number of collocation points in r in each domain
  fpar >> nr; fpar.ignore(1000,'\n') ;

  int nt ; // Number of collocation points in theta in each domain
  fpar >> nt; fpar.ignore(1000,'\n') ;

  int np ; // Number of collocation points in phi in each domain
  fpar >> np; fpar.ignore(1000,'\n') ;

  int nz ; // Number of domains
  fpar >> nz ; fpar.ignore(1000,'\n') ;
  int nzm1 = nz - 1 ; // Index of outermost domain

  fpar.ignore(1000,'\n') ; // skip title
  double* r_limits = new double[nz+1];  // inner boundaries of each domain in units of M      
  for (int l=0; l<nz; l++) 
    {
      fpar >> r_limits[l]; 
    }
  r_limits[nz] = __infinity ;
    
  fpar.close();


  // Setup of a multi-domain grid (Lorene class Mg3d)
  // ------------------------------------------------
  
  int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
  int symmetry_phi = SYM ; // symmetry with respect to phi --> phi + pi
  bool compact = true ; // external domain is compactified

  Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;

  //cout << mgrid << endl ; 

  
  // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
  // --------------------------------------------------------------------------
  
  Map_af map(mgrid, r_limits) ;

  // Construction of the ScalarBH object:
  // ----------------------------------

  ScalarBH object(map, file_name) ; 

  //object.update_metric() ; // TO BE ADDED WHEN WRITTEN

  cout.precision(15) ; 
  //cout << endl << "******* object ******** " << endl ;  
  //cout << object << endl ; 

  /*cout << "Lapse= " << endl;
  cout << object.get_nn() << endl;

  cout << "F0= " << endl;
  cout << object.get_ff0() << endl;*/

  cout << "betap far= " << object.get_beta()(3).val_point(1e5,M_PI/2.,0.) << endl;
  double rtest=10.;
  cout << "grr, gthth, gpp in Gyoto basis= " << object.get_gamma().cov()(1,1).val_point(rtest,M_PI/2.,0.) << " " << rtest*rtest*object.get_gamma().cov()(2,2).val_point(rtest,M_PI/2.,0.) << " " << rtest*rtest*object.get_gamma().cov()(3,3).val_point(rtest,M_PI/2.,0.) << endl;
  //cout << "lapse far= " << object.get_nn().val_point(10000,0.157,0.) << endl;
    
  // Drawings    
  if (graphic_out == 1) 
    {
      double r_min=0., r_max = 1.5;//1.5*map.val_r(nzm1,-1.,0.,0.) ; 
      //des_meridian(object.get_sfield(), 0, r_max, "Phi", 1) ; 

      
      des_meridian(object.get_nn(), r_min, r_max, "N", 1) ; 
      	
      des_meridian(object.get_gamma().cov()(1,1), r_min, r_max, "gamma_11", 2); 
      des_meridian(object.get_gamma().cov()(2,2), r_min, r_max, "gamma_22", 3) ; 
      des_meridian(object.get_gamma().cov()(3,3), r_min, r_max, "gamma_33", 4) ; 

      des_meridian(object.get_beta()(3), r_min, r_max, "Nphi", 5) ; 
      //des_meridian(object.get_beta()(3), 1e10, 1e12, "Nphi", 6) ; 
      //des_meridian(object.get_beta()(3), 1.62, 1.63, "Nphi", 6) ; 
	
      des_meridian(object.get_kk()(1,3), r_min, r_max, "K_(r)(ph)", 7) ; 
      des_meridian(object.get_kk()(2,3), r_min, r_max, "K_(th)(ph)", 8) ; 
      
      arrete() ; 
    }

  // Computing useful doughnut quantities
  /*
    This if loop computes quite a few relevant quantities that are necessary
    in order to properly define an accretion torus ("ion torus") in a given
    boson-star or scalar BH spacetime. These quantities will be saved in
    the file output.dat.
   */
  if (compute_qty == 1)
    {
      ofstream outputfile;
      outputfile.open("output.dat");
      int nsize = 50000;
      double rout=50.;
      double rstart=rHor;
      double dr=(rout-rstart)/double(nsize);
      double Madm = mass;
      //ell *= Madm; // no: ell is given in Mrond units already now
      for(int ii=1;ii<nsize;ii++){
	double rr = rstart+ii*dr;
	// equatorial plane quantities, theta=pi/2
	// Caution, switching from Lorene orthonormal to Gyoto coord basis:
	double betap_Lorene = object.get_beta()(3).val_point(rr,M_PI/2.,0.);
	double betap = 1./rr*betap_Lorene; 
	double gpp = rr*rr*object.get_gamma().cov()(3,3).val_point(rr,M_PI/2.,0.);
	double BB = sqrt(gpp/(rr*rr));
	double NN = object.get_nn().val_point(rr,M_PI/2.,0.);
	double Phi = fabs(object.get_sfield().val_point(rr,M_PI/2.,0.)); // I take fabs as the sign does not matter, I want it >0 to plot it
	double dbetapdr = 
	  1./rr*(object.get_beta()(3).dsdr().val_point(rr,M_PI/2.,0.) - betap);
	double dNdr = object.get_nn().dsdr().val_point(rr,M_PI/2.,0.);
	double gppr = 
	  rr*rr*object.get_gamma().cov()(3,3).dsdr().val_point(rr,M_PI/2.,0.)
	  +2./rr*gpp;
	double dBdr = 
	  1./(2.*BB)*object.get_gamma().cov()(3,3).dsdr().val_point(rr,M_PI/2.,0.);
	double DD = BB*BB*rr*rr/(NN*NN)*dbetapdr*dbetapdr
	  + 4.*dNdr/NN*(dBdr/BB+1./rr);
	double V_Kepler_ZAMO_plus = 
	  0.5*(-BB*rr/NN*dbetapdr+sqrt(DD))/(1./rr+dBdr/BB);
	double V_Kepler_ZAMO_minus = 
	  0.5*(-BB*rr/NN*dbetapdr-sqrt(DD))/(1./rr+dBdr/BB);
	double E_Kepler = 
	  1./sqrt(1-V_Kepler_ZAMO_plus*V_Kepler_ZAMO_plus)
	  *(NN-betap*BB*rr*V_Kepler_ZAMO_plus);
	double l_Kepler = BB*rr*V_Kepler_ZAMO_plus
	  /(NN-betap*BB*rr*V_Kepler_ZAMO_plus);

	double d2Ndr2 = object.get_nn().dsdr().dsdr().val_point(rr,M_PI/2.,0.);
	double d2lnN_over_dr2 = 1./(NN*NN)*(d2Ndr2*NN-dNdr*dNdr);
	double d2omega_over_dr2 = 
	  2./(rr*rr)*(object.get_beta()(3).dsdr().val_point(rr,M_PI/2.,0.) 
		      - betap)
	  -1./rr*object.get_beta()(3).dsdr().dsdr().val_point(rr,M_PI/2.,0.);
	double d2Bdr2 = 
	  1./BB*(0.5*object.get_gamma().cov()(3,3).dsdr().dsdr().val_point(rr,M_PI/2.,0.) - dBdr*dBdr);
	double d2lnB_over_dr2 = 1./(BB*BB)*(d2Bdr2*BB-dBdr*dBdr);
	double iscoeq = 
	  d2lnN_over_dr2 - 2.*dNdr*dNdr/(NN*NN) 
	  + V_Kepler_ZAMO_plus*BB*rr/NN*(d2omega_over_dr2 + 4.*dNdr/NN*dbetapdr)
	  + V_Kepler_ZAMO_plus*V_Kepler_ZAMO_plus*(-d2lnB_over_dr2+4./rr*dBdr/BB+2.*dBdr/BB*dBdr/BB+3./(rr*rr))
	  - V_Kepler_ZAMO_plus*V_Kepler_ZAMO_plus*BB*BB*rr*rr/(NN*NN)*dbetapdr*dbetapdr;

	double term1 = betap*BB*BB*rr*rr, term2 = -NN*NN+betap*term1;
	double delta = term1*term1-term2*BB*BB*rr*rr;
	double lmax = (-term1-sqrt(delta))/term2;
	double gtp = gpp*betap;
	double gtt = -NN*NN + gpp*betap*betap;
	double utdownstairs2_nume = gtp*gtp-gtt*gpp;
	double utdownstairs2_deno = gtt*ell*ell+2.*ell*gtp+gpp;
	double Eminus=0, Eplus=0, EminusD=0, EplusD=0, photonpot=0,
	  epsVzerominus=0, epsVzeroplus=0, particlepotential=0,
	  epsVprimeplus=0, epsVprimeminus=0, g_tt=0, g_rr=0, g_thth=0,
	  g_pp=0, g_tp=0; // just for fitting historical format

	// this should give the same output as the readKadathData.C routine for BS
	outputfile << ell << " " << Madm << " " << rr <<  " " << betap << " " << BB << " " << NN << " " << Phi << " " << DD  <<  " " << V_Kepler_ZAMO_plus << " " << E_Kepler << " " << l_Kepler  << " " << lmax << " " << utdownstairs2_nume << " " << utdownstairs2_deno << " "<< BB*BB*rr*rr/(NN*NN)*dbetapdr*dbetapdr  << " " << V_Kepler_ZAMO_minus << " " << dBdr << " " << dbetapdr << " " << dNdr << " " << particlepotential << " " << Eminus << " " << Eplus << " " << EminusD <<  " " << EplusD << " " << photonpot << " " << epsVzerominus << " " << epsVzeroplus << " " << epsVprimeminus << " " << epsVprimeplus << " " << g_tt << " " << g_rr << " " << g_thth << " " << g_pp << " " << g_tp << " " << iscoeq << endl;
      }
      
      outputfile.close();
    }
    
  //----------------------
  // Output file for GYOTO
  //----------------------
  object.gyoto_data("gyoto_scalarBH.d") ;    
    
  return EXIT_SUCCESS ; 
}













