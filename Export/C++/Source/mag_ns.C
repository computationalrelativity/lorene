/*
 * Methods of class Mag_NS (magnetized neutron star exportation)
 *
 * (see file mag_ns.h for documentation).
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *   Copyright (c) 2009 Jerome Novak
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


char mag_ns_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2009/11/22 16:05:05  j_novak
 * *** empty log message ***
 *
 * Revision 1.1  2009/11/19 16:15:21  j_novak
 * Export class for magnetized neutron stars.
 *
 *
 *
 * $Header$
 *
 */

#include "../Include/mag_ns.h"

void write_lines(ostream& fich, int dpl, const double* pdata, int np) ;


		    //-----------------------------------------//
		    //	    Constructor from a binary file     //
		    //-----------------------------------------//

Mag_NS::Mag_NS(FILE* fich) {

        fread(eos_name, sizeof(char), 100, fich) ;
        fread(&gamma_poly, sizeof(double), 1, fich) ;
        fread(&kappa_poly, sizeof(double), 1, fich) ;
 	fread(&omega, sizeof(double), 1, fich) ;
 	fread(&rho_c, sizeof(double), 1, fich) ;
 	fread(&eps_c, sizeof(double), 1, fich) ;
 	fread(&mass_b, sizeof(double), 1, fich) ;
 	fread(&mass_g, sizeof(double), 1, fich) ;
 	fread(&r_eq, sizeof(double), 1, fich) ;
 	fread(&r_p, sizeof(double), 1, fich) ;
 	fread(&angu_mom, sizeof(double), 1, fich) ;
 	fread(&T_over_W, sizeof(double), 1, fich) ;	
 	fread(&magn_mom, sizeof(double), 1, fich) ;	
 	fread(&b_z_pole, sizeof(double), 1, fich) ;	
 	fread(&b_z_eq, sizeof(double), 1, fich) ;	
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
	fread(bb_x, sizeof(double), np, fich) ;
	fread(bb_y, sizeof(double), np, fich) ;
	fread(bb_z, sizeof(double), np, fich) ;
	fread(jj_t, sizeof(double), np, fich) ;
	fread(jj_x, sizeof(double), np, fich) ;
	fread(jj_y, sizeof(double), np, fich) ;
	fread(jj_z, sizeof(double), np, fich) ;


}

		    //--------------------------------------------//
		    //	    Constructor from a formatted file     //
		    //--------------------------------------------//

Mag_NS::Mag_NS(ifstream& fich) {


    char comment[120] ; 

    fich.getline(eos_name,100) ;
    fich >> gamma_poly ; fich.getline(comment, 120) ;
    fich >> kappa_poly ; fich.getline(comment, 120) ;
    fich >> omega ; fich.getline(comment, 120) ;
    fich >> rho_c ; fich.getline(comment, 120) ;
    fich >> eps_c ; fich.getline(comment, 120) ;
    fich >> mass_b ; fich.getline(comment, 120) ;
    fich >> mass_g ; fich.getline(comment, 120) ;
    fich >> r_eq ; fich.getline(comment, 120) ;
    fich >> r_p ; fich.getline(comment, 120) ;
    fich >> angu_mom ; fich.getline(comment, 120) ;
    fich >> T_over_W ; fich.getline(comment, 120) ;
    fich >> magn_mom ; fich.getline(comment, 120) ;
    fich >> b_z_pole ; fich.getline(comment, 120) ;
    fich >> b_z_eq ; fich.getline(comment, 120) ;
    fich >> np ; fich.getline(comment, 120) ;

    alloc_memory() ;

    double* field[] = {xx, yy, zz,
		       nnn, beta_x, beta_y, beta_z,
		       g_xx, g_xy, g_xz, g_yy, g_yz, g_zz,
		       k_xx, k_xy, k_xz, k_yy, k_yz, k_zz,
		       nbar, ener_spec, u_euler_x, u_euler_y, u_euler_z,
                       bb_x, bb_y, bb_z,
                       jj_t, jj_x, jj_y, jj_z} ;

    for (int j=0; j<31; j++) {

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


Mag_NS::~Mag_NS() {

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

    delete [] bb_x ;
    delete [] bb_y ;
    delete [] bb_z ;

    delete [] jj_t ;
    delete [] jj_x ;
    delete [] jj_y ;
    delete [] jj_z ;

}


		    //--------------------------//
		    //	    Memory allocation 	//
		    //--------------------------//

void Mag_NS::alloc_memory() {

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

    bb_x = new double[np] ;
    bb_y = new double[np] ;
    bb_z = new double[np] ;

    jj_t = new double[np] ;
    jj_x = new double[np] ;
    jj_y = new double[np] ;
    jj_z = new double[np] ;

}

		    //--------------------------//
		    //	  	Outputs         //
		    //--------------------------//

ostream& operator<<(ostream& ost, const Mag_NS& mns) {

	ost << "Magnetized neutron star :" << endl ;
	ost << "-------------------------" << endl ;
        ost << "  EOS : " << mns.eos_name << endl ;
        ost << "      gamma = " << mns.gamma_poly << endl ;
        ost << "      kappa = " << mns.kappa_poly << " rho_nuc c^2 / n_nuc^gamma" << endl ;
	ost << "  Omega          : " << mns.omega << " rad/s" << endl ;
	ost << "  Central density   : " << mns.rho_c << " kg/m^3" << endl ;
	ost << "  Central specific energy : " << mns.eps_c << " c^2" << endl ;
	ost << "  Baryon mass  : " << mns.mass_b
	    << " M_sol" << endl ;
	ost << "  Gravitational mass  : " << mns.mass_g
	    << " M_sol" << endl ;
	ost << " Coordinate equatorial radius : " << mns.r_eq << " km" << endl ;
	ost << " Coordinate polar radius : " << mns.r_p << " km" << endl ;
	ost << "  Total angular momentum : " << mns.angu_mom
	    << " G M_sol^2/c" << endl ;
    	ost << "  T/W : " << mns.T_over_W << endl ;
    	ost << "  Magnetic momentum    : " << mns.magn_mom
	    << " A m^2" << endl ;
    	ost << "  Magnetic field at north pole : " << mns.b_z_pole
	    << "  x 10^9 T" << endl ;
    	ost << "  Tangent magnetic field equatorial value : " << mns.b_z_eq
	    << " x 10^9 T" << endl ;

	return ost ;
}


		    //----------------------------------//
		    //	    Save in a binary file 	//
		    //----------------------------------//

void Mag_NS::save_bin(FILE* fresu) const {

        fwrite(eos_name, sizeof(char), 100, fresu) ;
        fwrite(&gamma_poly, sizeof(double), 1, fresu) ;
        fwrite(&kappa_poly, sizeof(double), 1, fresu) ;
 	fwrite(&omega, sizeof(double), 1, fresu) ;
 	fwrite(&rho_c, sizeof(double), 1, fresu) ;
 	fwrite(&eps_c, sizeof(double), 1, fresu) ;
 	fwrite(&mass_b, sizeof(double), 1, fresu) ;
 	fwrite(&mass_g, sizeof(double), 1, fresu) ;
 	fwrite(&r_eq, sizeof(double), 1, fresu) ;
 	fwrite(&r_p, sizeof(double), 1, fresu) ;
 	fwrite(&angu_mom, sizeof(double), 1, fresu) ;
 	fwrite(&T_over_W, sizeof(double), 1, fresu) ;	
 	fwrite(&magn_mom, sizeof(double), 1, fresu) ;	
 	fwrite(&b_z_pole, sizeof(double), 1, fresu) ;	
 	fwrite(&b_z_eq, sizeof(double), 1, fresu) ;	
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
	fwrite(bb_x, sizeof(double), np, fresu) ;
	fwrite(bb_y, sizeof(double), np, fresu) ;
	fwrite(bb_z, sizeof(double), np, fresu) ;
	fwrite(jj_t, sizeof(double), np, fresu) ;
	fwrite(jj_x, sizeof(double), np, fresu) ;
	fwrite(jj_y, sizeof(double), np, fresu) ;
	fwrite(jj_z, sizeof(double), np, fresu) ;
   
}

		
		    //----------------------------------//
		    //	    Save in a formatted file 	//
		    //----------------------------------//

void Mag_NS::save_form(ofstream& fich) const {

    fich.precision(13) ;
    fich.setf(ios::scientific, ios::floatfield) ; 

    fich << eos_name << endl ;
    fich << gamma_poly << " gamma_poly" << endl ;
    fich << kappa_poly << " kappa_poly [rho_nuc c^2 / n_nuc^gamma]" << endl ;
    fich << omega << " Omega [rad/s]" << endl ;
    fich << rho_c << " central density [kg.m^{-3}]" << endl ;
    fich << eps_c << " central specific internal energy [c^2]" << endl ;
    fich << mass_b << " baryon mass [M_sol]" << endl ;
    fich << mass_g << " gravitational mass [M_sol]" << endl ;
    fich << r_eq << " equatorial coordinate radius [km]" << endl ;
    fich << r_p << " polar coordinate radius [km]" << endl ;
    fich << angu_mom << " angular momentum [G M_sol^2/c]" << endl ; 
    fich << T_over_W << " T/W" << endl ; 
    fich << magn_mom << " Magnetic momentum [A m^2]" << endl ; 
    fich << b_z_pole << " Magnetic field at the pole [10^9 T]" << endl ; 
    fich << b_z_eq << " Magnetic field at the equator [10^9 T]" << endl ; 
    fich << np << " Number of grid points" << endl ; 
    
    int dpl = 6 ;    // number of data per line 
    
    const double* field[] = {xx, yy, zz,
			     nnn, beta_x, beta_y, beta_z,
			     g_xx, g_xy, g_xz, g_yy, g_yz, g_zz,
			     k_xx, k_xy, k_xz, k_yy, k_yz, k_zz,
			     nbar, ener_spec, u_euler_x, u_euler_y, u_euler_z,
			     bb_x, bb_y, bb_z,
			     jj_t, jj_x, jj_y, jj_z} ;
			        
    for (int j=0; j<31; j++) {

         write_lines(fich, dpl, field[j], np)  ;

    }

}

