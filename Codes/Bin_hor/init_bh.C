/*
 * Main code for computing initial configuration
 *
 */

/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose Luis Jaramillo
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

char init_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/12/29 18:00:20  f_limousin
 * First version
 *
 * 
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>
#include <fstream.h>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "coord.h"
#include "tenseur.h"
#include "tensor.h"
#include "isol_hor.h"
#include "utilitaires.h"
#include "graphique.h"


int main(int argc, char** argv) {
    
    char blabla [120] ;
    ifstream param("par_init.d") ;
    
    double  precis, relax, radius, beta ;
    int nz, nbr, nbt, nbp ;
    
    param.getline(blabla, 120) ;
    param.getline(blabla, 120) ;
    param >> beta ; param.getline(blabla, 120) ;
    param >> nz ; param.getline(blabla, 120) ;
    param >> nbr ; param >> nbt ; param >> nbp ; param.getline(blabla, 120) ;
    
    double* bornes = new double[nz+1] ;
    for (int i=0 ; i<nz ; i++)
	param >> bornes[i] ;
    bornes[nz] = __infinity ;
    radius = bornes[1] ;
    param.getline(blabla, 120) ;
    
    param >> precis ; param.getline(blabla, 120) ;
    param >> relax ; param.getline(blabla, 120) ;
   
    double distance = radius*beta ;
    
    param.close() ;
    
    // Type of sampling in theta and phi :
    int type_t = SYM ; 
    int type_p = NONSYM ; 

    int* np = new int [nz] ;
    int* nt = new int [nz] ;
    int* nr = new int [nz] ;
    for (int l=0 ; l<nz ; l++){
	np[l] = nbp ;
    	nt[l] = nbt ;
 	nr[l] = nbr ;
    }
    
    int* type_r = new int[nz] ;
    type_r[0] = RARE ;
    for (int l=1 ; l<nz-1 ; l++)
	type_r[l] = FIN ;
    type_r[nz-1] = UNSURR ;
    
    Mg3d grid (nz, nr, type_r, nt, type_t, np, type_p) ;
    
    Map_af map_un (grid, bornes) ;
    Map_af map_deux (grid, bornes) ;
    
    map_un.set_ori (distance/2.,0, 0) ;
    map_deux.set_ori (-distance/2., 0, 0) ;
    map_deux.set_rot_phi (M_PI) ;

    int depth = 3 ;
    Bin_hor bin (map_un, map_deux, depth) ;
    bin.set_statiques(precis, relax) ;
    
    FILE* fich = fopen("static.d", "w") ;
    grid.sauve(fich) ;
    map_un.sauve(fich) ;
    map_deux.sauve(fich) ;
    bin(1).sauve(fich) ;
    bin(2).sauve(fich) ;
    fclose(fich) ;
  
    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    delete [] bornes ;

    return 1 ;
}
