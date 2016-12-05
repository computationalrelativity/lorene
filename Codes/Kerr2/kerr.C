/* Computes the Kerr metric in Dirac gauge.
 *
 */

/*
 *   Copyright (c) 2010 Nicolas Vasset
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


 

// headers Lorene
#include "nbr_spx.h"
#include "excised_slice.h"
#include "unites.h"

using namespace Lorene ;

int main(){ 
 
  using namespace Unites ; 
  //------------------------------------------------------------------
  //	    Parameters of the computation 
  //------------------------------------------------------------------

  int nr, nt, np, CFornot, mer_max, mer_max2;
  double relax, precis, Omega, N_hor;
  bool isCF;
  
  ifstream fpar("parkerr.d") ;
  fpar.ignore(1000, '\n') ;
  fpar.ignore(1000, '\n') ;
  fpar >> nr; fpar.ignore(1000, '\n');
  fpar >> nt; fpar.ignore(1000, '\n');
  fpar >> np; fpar.ignore(1000, '\n');
  fpar >> CFornot; fpar.ignore(1000, '\n');
  isCF = (CFornot == 1) ;
  
  int nz=6; // Number of domain; unable to change it in parameter file until now, bu can be changed inside the code.

  cout << "==========GRID PARAMETERS========" << endl;
  cout << "total number of domains :   nz = " << nz << endl ;
  cout << "number of points in r  :    nr = " << nr << endl ;
  cout << "number of points in phi :   np = " << np << endl ;
  cout << "number of points in theta : nt = " << nt << endl ;
  
  fpar >> relax; fpar.ignore(1000, '\n');
  fpar >> precis; fpar.ignore(1000, '\n');
  fpar >> mer_max; fpar.ignore(1000, '\n');
  fpar >> mer_max2; fpar.ignore(1000, '\n');
  fpar >> Omega; fpar.ignore(1000, '\n');
  fpar >> N_hor; fpar.ignore(1000, '\n');
  
  cout << "=============PHYSICAL PARAMETERS===============" << endl;
  cout << "Chosen rotation parameter :             Omega= " << Omega << endl;
  cout << "Value for the lapse at the horizon :    N=     " << N_hor << endl;
  
  fpar.close();
    
  //----------------------------------------------------------
  // Construction of a multi-grid (Mg3d) and associated mapping
  // ----------------------------------------------------------
  
  // Note that the horizon radius is set by default to 1, which is
  // the inner boundary of the first shell.
  
  int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
  int symmetry_phi = SYM ; // symmetry in phi
  int nbr[] = {nr, nr, nr, nr, nr, nr};
  int nbt[] = {nt, nt, nt, nt, nt, nt} ;
  int nbp[] = {np, np, np, np, np, np} ;
  int tipe_r[] = {RARE, FIN, FIN, FIN, FIN, UNSURR} ;
  
  Mg3d mgrid(nz, nbr, tipe_r, nbt, symmetry_theta, nbp, symmetry_phi) ;
  
  // Construct angular grid for h(theta,phi) 
  const Mg3d& g_angu = *mgrid.get_angu_1dom() ;
  
  // Construction of an affine mapping (Map_af)
  // ------------------------------------------
  
  // Boundaries of each domains
  double Rmax = (1/0.3) ;
  double r_limits[] = {0., 0.3*Rmax, 0.6*Rmax, 1.*Rmax, 2.*Rmax, 4.*Rmax, __infinity} ; 
  double r_limits2[] = {0.3*Rmax, 0.6*Rmax} ; 
 
  const Map_af map(mgrid, r_limits); 
  const Map_af map_2(g_angu, r_limits2);

  // Some helpful stuff...
  const Coord& rr = map.r;
  Scalar rrr (map) ; 
  rrr = rr ; 
  rrr.std_spectral_base(); 
  //  const Metric_flat& mets = map.flat_met_spher() ;

  // Definition of the class for black hole data containing an isolated horizon
  //---------------------------------------------------------------------------
  Scalar lapse_hor(map); lapse_hor = N_hor; lapse_hor.std_spectral_base();
  Excised_slice Kerr_hole(*(dynamic_cast<const Map*>(&map)), 1 ,1);
  
  // Compute the metric data using the Isol_hole class  
  //--------------------------------------------------
  Kerr_hole.compute_stat_metric(precis, Omega, false, lapse_hor, isCF, relax, mer_max, mer_max2);

  // Once metric fields are calculated, it is possible to perform several diagnostics on the data.

  cout << "==============================================" << endl;
  cout << "        Computation for the ADM mass" << endl;
  cout << "        (LORENE Units)              " << endl;
  cout << "        M_ADM =  " << Kerr_hole.adm_mass() << endl;
  cout << "==============================================" << endl;
  
  cout << "==============================================" << endl;
  cout << "        Total Komar angular momentum          " << endl;
  cout << "        (LORENE Units)              " << endl;
  cout << "        J_K =    " << Kerr_hole.komar_angmom() << endl;
  cout << "==============================================" << endl;

  cout << "==============================================" << endl;
  cout << "       Virial residue (stationarity marker)   " << endl;
  cout << "       Difference between Komar and ADM masses" << endl;
  cout << "       (Rescaled values)                     " << endl;
  cout << "       Vir =     " << Kerr_hole.virial_residue() << endl;
  cout << "==============================================" << endl;
  
  cout << "==============================================" << endl;
  cout << "       Violation of Einstein Equations        " << endl;

  Kerr_hole.Einstein_errors();
  
  cout << "==============================================" << endl;
  

  // The computed metric fields are saved as a whole in Isol_hole data. 
  
  FILE* Kerr_holedata = fopen("Kerr_data_Om0.05_N0.55.d", "w");
  Kerr_hole.get_mp().get_mg()->sauve(Kerr_holedata) ;		// writing of the grid
  Kerr_hole.get_mp().sauve(Kerr_holedata) ;      
  Kerr_hole.sauve(Kerr_holedata);

  fclose(Kerr_holedata);

  return EXIT_SUCCESS ; 
        
}  
  
    
