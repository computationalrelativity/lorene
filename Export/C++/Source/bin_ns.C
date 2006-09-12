/*
 * Methods of class Bin_NS (binary neutron star exportation)
 *
 * (see file bin_ns.h for documentation).
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *   Copyright (c) 2002  Keisuke Taniguchi
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


char bin_ns_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2006/09/12 08:04:06  j_novak
 * Removal of the include path Export/C++/Include, updating of the relevant
 * source files in Export/C++/Source.
 *
 * Revision 1.1  2002/01/11 17:03:02  e_gourgoulhon
 * Exportation of binary neutron stars configuration to a Cartesian grid
 *
 *
 * $Header$
 *
 */

#include "../Include/bin_ns.h"

void write_lines(ostream& fich, int dpl, const double* pdata, int np) ;


		    //-----------------------------------------//
		    //	    Constructor from a binary file     //
		    //-----------------------------------------//

Bin_NS::Bin_NS(FILE* fich) {

        fread(eos_name1, sizeof(char), 100, fich) ;
        fread(&gamma_poly1, sizeof(double), 1, fich) ;
        fread(&kappa_poly1, sizeof(double), 1, fich) ;
        fread(eos_name2, sizeof(char), 100, fich) ;
        fread(&gamma_poly2, sizeof(double), 1, fich) ;
        fread(&kappa_poly2, sizeof(double), 1, fich) ;
 	fread(&omega, sizeof(double), 1, fich) ;
 	fread(&dist, sizeof(double), 1, fich) ;
 	fread(&dist_mass, sizeof(double), 1, fich) ;
 	fread(&mass1_b, sizeof(double), 1, fich) ;
 	fread(&mass2_b, sizeof(double), 1, fich) ;
 	fread(&mass_adm, sizeof(double), 1, fich) ;
 	fread(&angu_mom, sizeof(double), 1, fich) ;
 	fread(&rad1_x_comp, sizeof(double), 1, fich) ;	
 	fread(&rad1_y, sizeof(double), 1, fich) ;	
 	fread(&rad1_z, sizeof(double), 1, fich) ;	
 	fread(&rad1_x_opp, sizeof(double), 1, fich) ;	
 	fread(&rad2_x_comp, sizeof(double), 1, fich) ;	
 	fread(&rad2_y, sizeof(double), 1, fich) ;	
 	fread(&rad2_z, sizeof(double), 1, fich) ;	
 	fread(&rad2_x_opp, sizeof(double), 1, fich) ;	
 	fread(&np, sizeof(int), 1, fich) ;
	
	alloc_memory() ;

	fread(xx, sizeof(double), np, fich) ;
	fread(yy, sizeof(double), np, fich) ;
	fread(zz, sizeof(double), np, fich) ;
	fread(nnn, sizeof(double), np, fich) ;
	fread(beta_x, sizeof(double), np, fich) ;
	fread(beta_y, sizeof(double), np, fich) ;
	fread(beta_z, sizeof(double), np, fich) ;
	fread(g_xx, sizeof(double), np, fich) ;
	fread(g_xy, sizeof(double), np, fich) ;
	fread(g_xz, sizeof(double), np, fich) ;
	fread(g_yy, sizeof(double), np, fich) ;
	fread(g_yz, sizeof(double), np, fich) ;
	fread(g_zz, sizeof(double), np, fich) ;
	fread(k_xx, sizeof(double), np, fich) ;
	fread(k_xy, sizeof(double), np, fich) ;
	fread(k_xz, sizeof(double), np, fich) ;
	fread(k_yy, sizeof(double), np, fich) ;
	fread(k_yz, sizeof(double), np, fich) ;
	fread(k_zz, sizeof(double), np, fich) ;
	fread(nbar, sizeof(double), np, fich) ;
	fread(ener_spec, sizeof(double), np, fich) ;
	fread(u_euler_x, sizeof(double), np, fich) ;
	fread(u_euler_y, sizeof(double), np, fich) ;
	fread(u_euler_z, sizeof(double), np, fich) ;


}

		    //--------------------------------------------//
		    //	    Constructor from a formatted file     //
		    //--------------------------------------------//

Bin_NS::Bin_NS(ifstream& fich) {


    char comment[120] ; 

    fich.getline(eos_name1,100) ;
    fich >> gamma_poly1 ; fich.getline(comment, 120) ;
    fich >> kappa_poly1 ; fich.getline(comment, 120) ;
    fich.getline(eos_name2,100) ;
    fich >> gamma_poly2 ; fich.getline(comment, 120) ;
    fich >> kappa_poly2 ; fich.getline(comment, 120) ;
    fich >> omega ; fich.getline(comment, 120) ;
    fich >> dist ; fich.getline(comment, 120) ;
    fich >> dist_mass ; fich.getline(comment, 120) ;
    fich >> mass1_b ; fich.getline(comment, 120) ;
    fich >> mass2_b ; fich.getline(comment, 120) ;
    fich >> mass_adm ; fich.getline(comment, 120) ;
    fich >> angu_mom ; fich.getline(comment, 120) ;
    fich >> rad1_x_comp ; fich.getline(comment, 120) ;
    fich >> rad1_y ; fich.getline(comment, 120) ;
    fich >> rad1_z ; fich.getline(comment, 120) ;
    fich >> rad1_x_opp ; fich.getline(comment, 120) ;
    fich >> rad2_x_comp ; fich.getline(comment, 120) ;
    fich >> rad2_y ; fich.getline(comment, 120) ;
    fich >> rad2_z ; fich.getline(comment, 120) ;
    fich >> rad2_x_opp ; fich.getline(comment, 120) ;
    fich >> np ; fich.getline(comment, 120) ;

    alloc_memory() ;

    double* field[] = {xx, yy, zz,
		       nnn, beta_x, beta_y, beta_z,
		       g_xx, g_xy, g_xz, g_yy, g_yz, g_zz,
		       k_xx, k_xy, k_xz, k_yy, k_yz, k_zz,
		       nbar, ener_spec, u_euler_x, u_euler_y, u_euler_z} ;

    for (int j=0; j<24; j++) {

	double* pdata = field[j] ;

	for (int i=0; i<np; i++) {
 	    fich >> *pdata ;
	    pdata++ ;
	}
    }
    
}




		    //--------------------------//
		    //	    Destructor		//
		    //--------------------------//


Bin_NS::~Bin_NS() {

    delete [] xx ; 
    delete [] yy ; 
    delete [] zz ; 
   
    delete [] nnn ;

    delete [] beta_x ;
    delete [] beta_y ;
    delete [] beta_z ;

    delete [] g_xx ;
    delete [] g_xy ;
    delete [] g_xz ;
    delete [] g_yy ;
    delete [] g_yz ;
    delete [] g_zz ;

    delete [] k_xx ;
    delete [] k_xy ;
    delete [] k_xz ;
    delete [] k_yy ;
    delete [] k_yz ;
    delete [] k_zz ;

    delete [] nbar ;
    delete [] ener_spec ;

    delete [] u_euler_x ;
    delete [] u_euler_y ;
    delete [] u_euler_z ;

}


		    //--------------------------//
		    //	    Memory allocation 	//
		    //--------------------------//

void Bin_NS::alloc_memory() {

    xx = new double[np] ;
    yy = new double[np] ;
    zz = new double[np] ;

    nnn = new double[np] ;

    beta_x = new double[np] ;
    beta_y = new double[np] ;
    beta_z = new double[np] ;

    g_xx = new double[np] ;
    g_xy = new double[np] ;
    g_xz = new double[np] ;
    g_yy = new double[np] ;
    g_yz = new double[np] ;
    g_zz = new double[np] ;

    k_xx = new double[np] ;
    k_xy = new double[np] ;
    k_xz = new double[np] ;
    k_yy = new double[np] ;
    k_yz = new double[np] ;
    k_zz = new double[np] ;

    nbar = new double[np] ;
    ener_spec = new double[np] ;

    u_euler_x = new double[np] ;
    u_euler_y = new double[np] ;
    u_euler_z = new double[np] ;

}

		    //--------------------------//
		    //	  	Outputs         //
		    //--------------------------//

ostream& operator<<(ostream& ost, const Bin_NS& bns) {

	ost << "Binary neutron star system :" << endl ;
	ost << "----------------------------" << endl ;
        ost << "  EOS of star 1 : " << bns.eos_name1 << endl ;
        ost << "      gamma = " << bns.gamma_poly1 << endl ;
        ost << "      kappa = " << bns.kappa_poly1 << " rho_nuc c^2 / n_nuc^gamma" << endl ;
        ost << "  EOS of star 2 : " << bns.eos_name2 << endl ;
        ost << "      gamma = " << bns.gamma_poly2 << endl ;
        ost << "      kappa = " << bns.kappa_poly2 << " rho_nuc c^2 / n_nuc^gamma" << endl ;
	ost << "  Omega          : " << bns.omega << " rad/s" << endl ;
	ost << "  Separation d   : " << bns.dist << " km" << endl ;
	ost << "  Separation d_G : " << bns.dist_mass << " km" << endl ;
	ost << "  Baryon mass of star 1  : " << bns.mass1_b
	    << " M_sol" << endl ;
	ost << "  Baryon mass of star 2  : " << bns.mass2_b
	    << " M_sol" << endl ;
	ost << "  ADM mass of the system : " << bns.mass_adm
	    << " M_sol" << endl ;
	ost << "  Total angular momentum : " << bns.angu_mom
	    << " G M_sol^2/c" << endl ;
    	ost << "  Radius of star 1 (x_comp) : " << bns.rad1_x_comp
	    << " km" << endl ;
    	ost << "  Radius of star 1 (y)      : " << bns.rad1_y
	    << " km" << endl ;
    	ost << "  Radius of star 1 (z)      : " << bns.rad1_z
	    << " km" << endl ;
    	ost << "  Radius of star 1 (x_opp)  : " << bns.rad1_x_opp
	    << " km" << endl ;
    	ost << "  Radius of star 2 (x_comp) : " << bns.rad2_x_comp
	    << " km" << endl ;
    	ost << "  Radius of star 2 (y)      : " << bns.rad2_y
	    << " km" << endl ;
    	ost << "  Radius of star 2 (z)      : " << bns.rad2_z
	    << " km" << endl ;
    	ost << "  Radius of star 2 (x_opp)  : " << bns.rad2_x_opp
	    << " km" << endl ;

	return ost ;
}


		    //----------------------------------//
		    //	    Save in a binary file 	//
		    //----------------------------------//

void Bin_NS::save_bin(FILE* fresu) const {

        fwrite(eos_name1, sizeof(char), 100, fresu) ;
        fwrite(&gamma_poly1, sizeof(double), 1, fresu) ;
        fwrite(&kappa_poly1, sizeof(double), 1, fresu) ;
        fwrite(eos_name2, sizeof(char), 100, fresu) ;
        fwrite(&gamma_poly2, sizeof(double), 1, fresu) ;
        fwrite(&kappa_poly2, sizeof(double), 1, fresu) ;
 	fwrite(&omega, sizeof(double), 1, fresu) ;
 	fwrite(&dist, sizeof(double), 1, fresu) ;
 	fwrite(&dist_mass, sizeof(double), 1, fresu) ;
 	fwrite(&mass1_b, sizeof(double), 1, fresu) ;
 	fwrite(&mass2_b, sizeof(double), 1, fresu) ;	
 	fwrite(&mass_adm, sizeof(double), 1, fresu) ;	
 	fwrite(&angu_mom, sizeof(double), 1, fresu) ;	
 	fwrite(&rad1_x_comp, sizeof(double), 1, fresu) ;	
 	fwrite(&rad1_y, sizeof(double), 1, fresu) ;	
 	fwrite(&rad1_z, sizeof(double), 1, fresu) ;	
 	fwrite(&rad1_x_opp, sizeof(double), 1, fresu) ;	
 	fwrite(&rad2_x_comp, sizeof(double), 1, fresu) ;	
 	fwrite(&rad2_y, sizeof(double), 1, fresu) ;	
 	fwrite(&rad2_z, sizeof(double), 1, fresu) ;	
 	fwrite(&rad2_x_opp, sizeof(double), 1, fresu) ;
 	fwrite(&np, sizeof(int), 1, fresu) ;	
	fwrite(xx, sizeof(double), np, fresu) ; 
	fwrite(yy, sizeof(double), np, fresu) ; 
	fwrite(zz, sizeof(double), np, fresu) ; 
	fwrite(nnn, sizeof(double), np, fresu) ; 
	fwrite(beta_x, sizeof(double), np, fresu) ; 
	fwrite(beta_y, sizeof(double), np, fresu) ; 
	fwrite(beta_z, sizeof(double), np, fresu) ; 
	fwrite(g_xx, sizeof(double), np, fresu) ;
	fwrite(g_xy, sizeof(double), np, fresu) ; 
	fwrite(g_xz, sizeof(double), np, fresu) ; 
	fwrite(g_yy, sizeof(double), np, fresu) ; 
	fwrite(g_yz, sizeof(double), np, fresu) ; 
	fwrite(g_zz, sizeof(double), np, fresu) ;
	fwrite(k_xx, sizeof(double), np, fresu) ; 
	fwrite(k_xy, sizeof(double), np, fresu) ; 
	fwrite(k_xz, sizeof(double), np, fresu) ; 
	fwrite(k_yy, sizeof(double), np, fresu) ; 
	fwrite(k_yz, sizeof(double), np, fresu) ; 
	fwrite(k_zz, sizeof(double), np, fresu) ; 
	fwrite(nbar, sizeof(double), np, fresu) ; 
	fwrite(ener_spec, sizeof(double), np, fresu) ; 
	fwrite(u_euler_x, sizeof(double), np, fresu) ; 
	fwrite(u_euler_y, sizeof(double), np, fresu) ; 
	fwrite(u_euler_z, sizeof(double), np, fresu) ; 
   
}

		
		    //----------------------------------//
		    //	    Save in a formatted file 	//
		    //----------------------------------//

void Bin_NS::save_form(ofstream& fich) const {

    fich.precision(13) ;
    fich.setf(ios::scientific, ios::floatfield) ; 

    fich << eos_name1 << endl ;
    fich << gamma_poly1 << " gamma_poly1" << endl ;
    fich << kappa_poly1 << " kappa_poly1 [rho_nuc c^2 / n_nuc^gamma]" << endl ;
    fich << eos_name2 << endl ;
    fich << gamma_poly2 << " gamma_poly2" << endl ;
    fich << kappa_poly2 << " kappa_poly2 [rho_nuc c^2 / n_nuc^gamma]" << endl ;
    fich << omega << " Omega [rad/s]" << endl ;
    fich << dist << " dist [km]" << endl ;
    fich << dist_mass << " dist_mass [km]" << endl ;
    fich << mass1_b << " mass1_b [M_sol]" << endl ;
    fich << mass2_b << " mass2_b [M_sol]" << endl ;
    fich << mass_adm << " mass_adm [M_sol]" << endl ;
    fich << angu_mom << " angu_mom [G M_sol^2/c]" << endl ; 
    fich << rad1_x_comp << " rad1_x_comp [km]" << endl ; 
    fich << rad1_y << " rad1_y [km]" << endl ; 
    fich << rad1_z << " rad1_z [km]" << endl ; 
    fich << rad1_x_opp << " rad1_x_opp [km]" << endl ; 
    fich << rad2_x_comp << " rad1_x_comp [km]" << endl ; 
    fich << rad2_y << " rad1_y [km]" << endl ; 
    fich << rad2_z << " rad1_z [km]" << endl ; 
    fich << rad2_x_opp << " rad1_x_opp [km]" << endl ;
    fich << np << " Number of grid points" << endl ; 
    
    int dpl = 6 ;    // number of data per line 
    
    const double* field[] = {xx, yy, zz,
			     nnn, beta_x, beta_y, beta_z,
			     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz,
			     k_xx, k_xy, k_xz, k_yy, k_yz, k_zz,
			     nbar, ener_spec, u_euler_x, u_euler_y, u_euler_z} ;
			        
    for (int j=0; j<24; j++) {

         write_lines(fich, dpl, field[j], np)  ;

    }

}

