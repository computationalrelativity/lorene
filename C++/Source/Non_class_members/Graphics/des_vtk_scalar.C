/*
 *   Copyright (c) 2012-2013 Jason Penner
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
 * Written with the assistance of Ian Hawke, University of Southampton
 */

 

// Header C
#include <cmath>
#include <iostream>
#include <string>

// Header Lorene
#include "graphique_vtk.h"
#include "graphique.h"
#include "scalar.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
void des_coupe_vtk_x(const Scalar& uu, double x0, double y_min, double y_max, 
		 double z_min, double z_max, const char* title, int ny, int nz) {
		
using namespace Unites ;

    const Map& mp = uu.get_mp() ; 

    int nx = 1;
    double hx = 0;
    double x_min = x0;

    int ii=0;
    int i,j;
    char c;

    // Dump Scalar Slice to VTK Format 
    // -------------------------------
    
    double hy = (y_max - y_min) / double(ny-1) ; 
    double hza = (z_max - z_min) / double(nz-1) ; 

// output file management
    char title2[256] = {' '};
    char title3[256] = {' '};

    while(title[ii]) {
        c = char(tolower(title[ii]));
        if(c==' ') {
           c = '_';
        }
	title2[ii] = c ;
        ii++;
    }
    strcat(title3,title2); 
    strcat(title2,"_x.vtk");  
    ofstream myfile;
    myfile.open (title2);


    myfile << "# vtk DataFile Version 3.0" << endl;
    myfile << "Data produced by mf_evolve" << endl;
    myfile << "ASCII" << endl;
    myfile << "DATASET RECTILINEAR_GRID" << endl;
    myfile << "DIMENSIONS" << " " << nx << " " << ny << " " << nz << endl;

    cout.precision(12);
    cout.setf(ios::scientific);

    myfile << "X_COORDINATES " << nx << " FLOAT" << endl;
    for ( i=0; i<nx; i++){
      myfile << scientific << setprecision(12) << "  " << x_min + hx * i << " " ;
    } 
    myfile << endl;

    myfile << "Y_COORDINATES " << ny << " FLOAT" << endl;
    for ( i=0; i<ny; i++){
      myfile << "  " << y_min + hy * i << " " ;
    } 
    myfile << endl;

    myfile << "Z_COORDINATES " << nz << " FLOAT" << endl;
    for ( i=0; i<nz; i++){
      myfile << "  " << z_min + hza * i << " " ;
    } 
    myfile << endl;

    myfile << "POINT_DATA " << nx*ny*nz << endl;

    myfile << endl;

    myfile << "SCALARS " << title3 << " FLOAT" << endl;

    myfile << "LOOKUP_TABLE default" << endl;

    for ( j=0; j<nz; j++) {
	
	double z = z_min + hza * j ; 
	
	for ( i=0; i<ny; i++) {
    
	    double y = y_min + hy * i ; 

	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x0, y, z, r, theta, phi) ; 
	myfile << "  " << float(uu.val_point(r, theta, phi)) << endl;
	}
    }

    myfile.close();
} 

//******************************************************************************

void des_coupe_vtk_x(const Scalar& uu, double x0, int nzdes, const char* title, 
		 double zoom, int ny, int nz) {
		     
    const Map& mp = uu.get_mp() ; 

    double a1 = mp.val_r(nzdes-1, 1., M_PI/2., 0.) ; 		 
    double a2 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI/2.) ; 		 
    double a3 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI) ; 		 
    double ray = mp.val_r(nzdes-1, 1., 0., 0.) ; 
    
    ray = ( a1 > ray ) ? a1 : ray ; 
    ray = ( a2 > ray ) ? a2 : ray ; 
    ray = ( a3 > ray ) ? a3 : ray ; 
    		 
    ray *= zoom ; 
    
    double y_min = mp.get_ori_y() - ray ; 
    double y_max = mp.get_ori_y() + ray ; 
    double z_min = mp.get_ori_z() - ray ; 
    double z_max = mp.get_ori_z() + ray ; 

    printf("Printing %s to file\n",title);
    des_coupe_vtk_x(uu, x0, y_min, y_max, z_min, z_max, title, ny, nz) ;

}


//******************************************************************************

void des_coupe_vtk_y(const Scalar& uu, double y0, double x_min, double x_max, 
		 double z_min, double z_max, const char* title, int nx, int nz) {
		
  using namespace Unites ;
  
    const Map& mp = uu.get_mp() ; 

    int ii = 0;
    int ny = 1;
    int i,j;
    char c;

    double y_min = y0;
    double hy = 0;

    // Plot of isocontours
    // -------------------
              
    double hx = (x_max - x_min) / double(nx-1) ; 
    double hza = (z_max - z_min) / double(nz-1) ; 
    
// output file management
    char title2[256] = {' '};
    char title3[256] = {' '};

    while(title[ii]) {
        c = char(tolower(title[ii]));
        if(c==' ') {
           c = '_';
        }
	title2[ii] = c ;
        ii++;
    }
    strcat(title3,title2); 
    strcat(title2,"_y.vtk");  
    ofstream myfile;
    myfile.open (title2);


    myfile << "# vtk DataFile Version 3.0" << endl;
    myfile << "Data produced by mf_evolve" << endl;
    myfile << "ASCII" << endl;
    myfile << "DATASET RECTILINEAR_GRID" << endl;
    myfile << "DIMENSIONS" << " " << nx << " " << ny << " " << nz << endl;

    cout.precision(12);
    cout.setf(ios::scientific);

    myfile << "X_COORDINATES " << nx << " FLOAT" << endl;
    for ( i=0; i<nx; i++){
      myfile << scientific << setprecision(12) << "  " << x_min + hx * i << " " ;
    } 
    myfile << endl;

    myfile << "Y_COORDINATES " << ny << " FLOAT" << endl;
    for ( i=0; i<ny; i++){
      myfile << "  " << y_min + hy * i << " " ;
    } 
    myfile << endl;

    myfile << "Z_COORDINATES " << nz << " FLOAT" << endl;
    for ( i=0; i<nz; i++){
      myfile << "  " << z_min + hza * i << " " ;
    } 
    myfile << endl;

    myfile << "POINT_DATA " << nx*ny*nz << endl;

    myfile << endl;

    myfile << "SCALARS " << title3 << " FLOAT" << endl;

    myfile << "LOOKUP_TABLE default" << endl;


    for (j=0; j<nz; j++) {
	
	double z = z_min + hza * j ; 
	
	for (i=0; i<nx; i++) {
    
	    double x = x_min + hx * i ; 
	    
	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x, y0, z, r, theta, phi) ; 
	myfile << "  " << float(uu.val_point(r, theta, phi)) << endl;
	}
    }

   myfile.close() ;
} 

//******************************************************************************

void des_coupe_vtk_y(const Scalar& uu, double y0, int nzdes, const char* title, 
		 double zoom, int nx, int nz) {
		     
    const Map& mp = uu.get_mp() ; 

    double a1 = mp.val_r(nzdes-1, 1., M_PI/2., 0.) ; 		 
    double a2 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI/2.) ; 		 
    double a3 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI) ; 		 
    double ray = mp.val_r(nzdes-1, 1., 0., 0.) ; 
    
    ray = ( a1 > ray ) ? a1 : ray ; 
    ray = ( a2 > ray ) ? a2 : ray ; 
    ray = ( a3 > ray ) ? a3 : ray ; 
    		 
    ray *= zoom ; 
    
    double x_min = mp.get_ori_x() - ray ; 
    double x_max = mp.get_ori_x() + ray ; 
    double z_min = mp.get_ori_z() - ray ; 
    double z_max = mp.get_ori_z() + ray ; 

    printf("Printing %s to file\n",title);
    des_coupe_vtk_y(uu, y0, x_min, x_max, z_min, z_max, title, nx, nz) ;

}

//******************************************************************************

void des_coupe_vtk_z(const Scalar& uu, double z0, double x_min, double x_max, 
		 double y_min, double y_max, const char* title, int nx, int ny) {
		
  using namespace Unites ;
  
    const Map& mp = uu.get_mp() ; 

    int ii = 0;
    int nz = 1;
    int i,j;

    char c;

    double hza = 0;
    double z_min = z0;

    // Plot of isocontours
    // -------------------
    
    double hy = (y_max - y_min) / double(ny-1) ; 
    double hx = (x_max - x_min) / double(nx-1) ; 

// output file management
    char title2[256] = {' '};
    char title3[256] = {' '};

    while(title[ii]) {
        c = char(tolower(title[ii]));
        if(c==' ') {
           c = '_';
        }
	title2[ii] = c ;
        ii++;
    }
    strcat(title3,title2); 
    strcat(title2,"_z.vtk");  
    ofstream myfile;
    myfile.open (title2);


    myfile << "# vtk DataFile Version 3.0" << endl;
    myfile << "Data produced by mf_evolve" << endl;
    myfile << "ASCII" << endl;
    myfile << "DATASET RECTILINEAR_GRID" << endl;
    myfile << "DIMENSIONS" << " " << nx << " " << ny << " " << nz << endl;

    cout.precision(12);
    cout.setf(ios::scientific);

    myfile << "X_COORDINATES " << nx << " FLOAT" << endl;
    for ( i=0; i<nx; i++){
      myfile << scientific << setprecision(12) << "  " << x_min + hx * i << " " ;
    } 
    myfile << endl;

    myfile << "Y_COORDINATES " << ny << " FLOAT" << endl;
    for ( i=0; i<ny; i++){
      myfile << "  " << y_min + hy * i << " " ;
    } 
    myfile << endl;

    myfile << "Z_COORDINATES " << nz << " FLOAT" << endl;
    for ( i=0; i<nz; i++){
      myfile << "  " << z_min + hza * i << " " ;
    } 
    myfile << endl;

    myfile << "POINT_DATA " << nx*ny*nz << endl;

    myfile << endl;

    myfile << "SCALARS " << title3 << " FLOAT" << endl;

    myfile << "LOOKUP_TABLE default" << endl;
    for (j=0; j<ny; j++) {
	
	double y = y_min + hy * j ; 
	
	for (i=0; i<nx; i++) {
    
	    double x = x_min + hx * i ; 
	    
	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x, y, z0, r, theta, phi) ; 
	myfile << "  " << float(uu.val_point(r, theta, phi)) << endl;
	}
    }
   myfile.close();
} 
/***************************************************************************/

void des_coupe_vtk_z(const Scalar& uu, double z0, int nzdes, const char* title, 
		 double zoom, int nx, int ny) {
		     
    const Map& mp = uu.get_mp() ; 

    double a1 = mp.val_r(nzdes-1, 1., M_PI/2., 0.) ; 		 
    double a2 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI/2.) ; 		 
    double a3 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI) ; 		 
    double ray = mp.val_r(nzdes-1, 1., 0., 0.) ; 
    
    ray = ( a1 > ray ) ? a1 : ray ; 
    ray = ( a2 > ray ) ? a2 : ray ; 
    ray = ( a3 > ray ) ? a3 : ray ; 
    		 
    ray *= zoom ; 
    
    double x_min = mp.get_ori_x() - ray ; 
    double x_max = mp.get_ori_x() + ray ; 
    double y_min = mp.get_ori_y() - ray ; 
    double y_max = mp.get_ori_y() + ray ; 

    printf("Printing %s to file\n",title);
    des_coupe_vtk_z(uu, z0, x_min, x_max, y_min, y_max, title, nx, ny) ;

}

//******************************************************************************


void des_vtk_xyz(const Scalar& uu, double x_min, double x_max, double y_min, 
                 double y_max, double z_min, double z_max,
                 const char* title, int nx, int ny, int nz) {
		
using namespace Unites ;

    const Map& mp = uu.get_mp() ; 

    int ii=0;
    int i,j,k;
    char c;

    // Dump Scalar Slice to Ascii 
    // --------------------------
       
    double hx = (x_max - x_min) / double(nx-1) ; 
    double hy = (y_max - y_min) / double(ny-1) ; 
    double hz = (z_max - z_min) / double(nz-1) ; 

// output file management
    char title2[256] = {' '};
    char title3[256] = {' '};

    while(title[ii]) {
        c = char(tolower(title[ii]));
        if(c==' ') {
           c = '_';
        }
	title2[ii] = c ;
        ii++;
    }
    strcat(title3,title2); 
    strcat(title2,"_xyz.vtk");  
    ofstream myfile;
    myfile.open (title2);


    myfile << "# vtk DataFile Version 3.0" << endl;
    myfile << "Data produced by mf_evolve" << endl;
    myfile << "ASCII" << endl;
    myfile << "DATASET RECTILINEAR_GRID" << endl;
    myfile << "DIMENSIONS" << " " << nx << " " << ny << " " << nz << endl;

    cout.precision(12);
    cout.setf(ios::scientific);

    myfile << "X_COORDINATES " << nx << " FLOAT" << endl;
    for ( i=0; i<nx; i++){
      myfile << scientific << setprecision(12) << "  " << x_min + hx * i << " " ;
    } 
    myfile << endl;

    myfile << "Y_COORDINATES " << ny << " FLOAT" << endl;
    for ( i=0; i<ny; i++){
      myfile << "  " << y_min + hy * i << " " ;
    } 
    myfile << endl;

    myfile << "Z_COORDINATES " << nz << " FLOAT" << endl;
    for ( i=0; i<nz; i++){
      myfile << "  " << z_min + hz * i << " " ;
    } 
    myfile << endl;

    myfile << "POINT_DATA " << nx*ny*nz << endl;

    myfile << endl;

    myfile << "SCALARS " << title3 << " FLOAT" << endl;

    myfile << "LOOKUP_TABLE default" << endl;

    for ( k=0; k<nx; k++) {

        double x = x_min + hx * k ;

      for ( j=0; j<nz; j++) {
	
  	double z = z_min + hz * j ; 
	
	for ( i=0; i<ny; i++) {
    
	    double y = y_min + hy * i ; 

	    // Computation of (r,theta,phi) : 	    
	    double r, theta, phi ; 
	    mp.convert_absolute(x, y, z, r, theta, phi) ;
        // calculate scalar field uu at (r,theta,phi)
	myfile << "  " << uu.val_point(r, theta, phi) << endl;
	}
      }
   }

    myfile.close();
} 

//******************************************************************************

void des_vtk_xyz(const Scalar& uu, int nzdes, const char* title, 
	         int nx, int ny, int nz) {
		     
    const Map& mp = uu.get_mp() ; 

    double a1 = mp.val_r(nzdes-1, 1., M_PI/2., 0.) ; 		 
    double a2 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI/2.) ; 		 
    double a3 = mp.val_r(nzdes-1, 1., M_PI/2., M_PI) ; 		 
    double ray = mp.val_r(nzdes-1, 1., 0., 0.) ; 
    
    ray = ( a1 > ray ) ? a1 : ray ; 
    ray = ( a2 > ray ) ? a2 : ray ; 
    ray = ( a3 > ray ) ? a3 : ray ; 
    
    double x_min = mp.get_ori_x() - ray ; 
    double x_max = mp.get_ori_x() + ray ; 
    double y_min = mp.get_ori_y() - ray ; 
    double y_max = mp.get_ori_y() + ray ; 
    double z_min = mp.get_ori_z() - ray ; 
    double z_max = mp.get_ori_z() + ray ; 

    printf("Printing %s to file\n",title);
    des_vtk_xyz(uu, x_min, x_max, y_min, y_max, z_min, z_max, title,
		nx, ny, nz) ;

}

//******************************************************************************
}
