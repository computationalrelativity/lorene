/*
 *  Methods of class Champ_cart
 *
 */

/*
 *   Copyright (c) 2002 Nicolas Chamel
 *   Copyright (c) 2002 Eric Gourgoulhon
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

char champ_cart_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:33  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1  2002/03/28 10:47:29  n_chamel
 * New class Champ_cart for fields on Cartesian grids.
 *
 *
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "headcpp.h"
#include "champ_cart.h"
#include "type_parite.h"
#include "utilitaires.h"

                                                               //-------------------------//
                                                                //              Constructors                      //
                                                               //-------------------------//


// Standard constructor
// ---------------

Champ_cart::Champ_cart(const Map_cart&  map_i) :
        map(map_i) ,
        etat(ETATNONDEF),
        val(0x0) ,
        cf(0x0)
{
        base[0] = NONDEF  ;
        base[1] = NONDEF  ;
        base[2] = NONDEF  ;
}

// Copy constructor
// ------------

Champ_cart::Champ_cart(const Champ_cart&  uu) :
        map(uu.map) ,
        etat(uu.etat)
{
        if (uu.etat == ETATQCQ) {

                if (uu.val != 0x0) {
                        val = new Tbl( *(uu.val) ) ;
                }
                else  {
                        val = 0x0 ;
                }

                if (uu.cf != 0x0) {
                        cf = new Tbl( *(uu.cf) ) ;
                }
                else  {
                        cf = 0x0 ;
                }

        }
        else {
                val = 0x0 ;
                cf = 0x0 ;
        }

        base[0] = uu.base[0]  ;
        base[1] = uu.base[1]  ;
        base[2] = uu.base[2]  ;
}

// Constructor from file
// --------------

Champ_cart::Champ_cart(const Map_cart&  map_i, FILE* fich) :
        map(map_i) ,
        cf(0x0)
{
        fread_be(&etat, sizeof(int), 1, fich) ;
        fread_be(base, sizeof(int), 3, fich) ;

        if (etat == ETATQCQ) {
	val = new Tbl(fich) ;
        }
        else {
                val = 0x0 ;
        }
}


// Destructor
// -------

Champ_cart::~Champ_cart(){
        del_all() ;
}


                                                                        // ---------------- //
                                                                        // Memory management   //
                                                                        // ---------------- //
void Champ_cart::del_all(){
        delete val ;
        delete cf ;
        val = 0x0 ;
        cf = 0x0 ;
        etat = ETATNONDEF ;
        base[0] = NONDEF  ;
        base[1] = NONDEF  ;
        base[2] = NONDEF  ;
}

                                                                // ---------------- //
                                                               // Mutators / assignment  //
                                                               // ---------------- //

void Champ_cart::operator=(const Champ_cart& uu)  {

        if (  &map != &(uu.map) ) {
                cout << " Champ_cart::operator= : the two Champ_cart must be defined on the same mapping !" << endl ;
                abort() ;
        }

        del_all() ;

        etat = uu.etat ;

       if (uu.etat == ETATQCQ) {

                if (uu.val != 0x0) {
                        val = new Tbl( *(uu.val) ) ;
                }
                else  {
                        val = 0x0 ;
                }

                if (uu.cf != 0x0) {
                        cf = new Tbl( *(uu.cf) ) ;
                }
                else  {
                        cf = 0x0 ;
                }

        }
        else {
                val = 0x0 ;
                cf = 0x0 ;
        }

        base[0] = uu.base[0]  ;
        base[1] = uu.base[1]  ;
        base[2] = uu.base[2]  ;

}


                                                                        // ------//
                                                                        // Outputs  //
                                                                        // ------//


// Save in file
// --------

void Champ_cart::sauve(FILE* fich) const {

        fwrite_be(&etat, sizeof(int), 1, fich) ;
        fwrite_be(base, sizeof(int), 3, fich) ;

         if (etat == ETATQCQ) {
	// Values at the grid points are required
                if (val==0x0) {
                        coef_i() ;
                }
                val->sauve(fich) ;
        }

}

/// Display
ostream& operator<<(ostream&  ost , const Champ_cart& uu) {

        ost << "Field on a Cartesian grid" << endl ;
        if (uu.etat == ETATNONDEF) {
                ost << " Undefined state" << endl ;
                return ost ;
        }

        if (uu.etat == ETATZERO) {
                ost << " Null state" << endl ;
                return ost ;
        }

        assert( uu.etat == ETATQCQ ) ;

        if (uu.val != 0x0) {
                ost << "Values at the grid points : " << endl ;
                ost << *(uu.val)  << endl ;
        }

        if (uu.cf != 0x0) {
                ost << "Values of the spectral coefficients : " << endl ;
                ost << *(uu.cf)  << endl ;
        }

        return ost ;
}

