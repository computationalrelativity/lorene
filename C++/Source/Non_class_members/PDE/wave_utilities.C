/*
 *  Miscellaneous functions for the wave equation
 *
 */

/*
 *   Copyright (c) 2008 Jerome Novak
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

char wave_utilities_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2008/07/18 12:28:41  j_novak
 * Corrected some mistakes.
 *
 * Revision 1.2  2008/07/18 09:17:35  j_novak
 * New function tilde_laplacian().
 *
 * Revision 1.1  2008/07/11 13:20:54  j_novak
 * Miscellaneous functions for the wave equation.
 *
 *
 * $Header $
 *
 */

#include"tensor.h"

void tilde_laplacian(const Scalar& B_in, Scalar& tilde_lap, int dl) {
    
    if (B_in.get_etat() == ETATZERO) {
	tilde_lap.set_etat_zero() ;
	return ;
    }
    if (dl == 0) {
	tilde_lap = B_in.laplacian(0) ;
	return ;
    }
    assert(B_in.get_etat() != ETATNONDEF) ;
    const Map_af* map =dynamic_cast<const Map_af*>(&B_in.get_mp()) ;
    assert(map != 0x0) ;
    
    tilde_lap = 2*B_in.dsdr() ;
    tilde_lap.div_r_dzpuis(3) ;
    tilde_lap += B_in.dsdr().dsdr() ;
    tilde_lap.dec_dzpuis() ;
    tilde_lap.set_spectral_va().ylm() ;
    Scalar B_over_r2 = B_in ;
    B_over_r2.div_r_dzpuis(1) ;
    B_over_r2.div_r_dzpuis(2) ;
    B_over_r2.set_spectral_va().ylm() ;

    const Base_val& base = B_in.get_spectral_base() ;
    const Mg3d& mg = *map->get_mg() ;
    int nz = mg.get_nzone() ;
    int l_q, m_q, base_r ;
    for (int lz=0; lz<nz; lz++) {
	if (B_in.domain(lz).get_etat() == ETATZERO)
	    tilde_lap.set_domain(lz).set_etat_zero() ;
	else {
	    for (int k=0; k<mg.get_np(lz)+2; k++)
		for (int j=0; j<mg.get_nt(lz); j++) {
		    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		    if (l_q > 1) {
			l_q += dl ;
			for (int i=0; i<mg.get_nr(lz); i++)
			    tilde_lap.set_spectral_va().c_cf->set(lz, k, j, i)
				-= l_q*(l_q+1)*
				(*B_over_r2.get_spectral_va().c_cf)(lz, k, j, i) ;
		    }
		}
	}
    }
    if (tilde_lap.get_spectral_va().c != 0x0) {
	delete tilde_lap.set_spectral_va().c ;
	tilde_lap.set_spectral_va().c = 0x0 ;
    }
    tilde_lap.dec_dzpuis(2) ;
    return ;
}

/* Performs one time-step integration of the wave equation, using
 * a third-order Runge-Kutta scheme.
 * phi = d fff / dt
 * \Delta fff = d phi / dt
 * Inputs are dt, fff, phi; outputs fnext, phinext.
 */
void runge_kutta3_wave_sys(double dt, const Scalar& fff, const Scalar& phi,
			   Scalar& fnext, Scalar& phinext, int dl) {

    const Map& map = fff.get_mp() ;
    Scalar k1 = phi ;
    Scalar dk1(map) ; tilde_laplacian(fff, dk1, dl) ;
    Scalar y1 = fff + 0.5*dt*k1 ;
    Scalar dy1 = phi + 0.5*dt*dk1 ;
    Scalar k2 = dy1 ; Scalar dk2(map) ; tilde_laplacian(y1, dk2, dl) ;
    Scalar y2 = fff - dt*k1 + 2*dt*k2 ;
    Scalar dy2 = phi - dt*dk1 + 2*dt*dk2 ;
    Scalar k3 = dy2 ;
    Scalar dk3(map) ; tilde_laplacian(y2, dk3, dl) ;
    fnext = fff + dt*(k1 + 4*k2 + k3)/6. ;
    phinext = phi + dt*(dk1 + 4*dk2 + dk3)/6. ;
    
    return ;
}

/* Performs one time-step integration of the quantities needed for the
 * enhanced outgoing-wave boundary condition. It DOES NOT impose the BC
 * d phi / dr + d phi / dt + phi / r = xi(theta, varphi).
 * nz_bound: index of the domain on which to impose the BC
 * phi: the field that should leave the grid
 * sphi: source of the Robin BC, without xi : a phi + b d phi / dr = sphi + xi
 * ccc: (output) total source of the Robin BC
 */ 
void evolve_outgoing_BC(double dt, int nz_bound, const Scalar& phi, Scalar& sphi, 
			Tbl& xij, Tbl& xijm1, Tbl& ccc, int dl) {
    
    const Map* map = &phi.get_mp() ;
    const Map_af* mp_aff = dynamic_cast<const Map_af*>(map) ;
    assert(mp_aff != 0x0) ;

    const Mg3d& grid = *mp_aff->get_mg() ;
#ifndef NDEBUG
    int nz = grid.get_nzone() ;
    assert(nz_bound < nz) ;
#endif
    int np2 = grid.get_np(nz_bound-1) + 2 ;
    int nt = grid.get_nt(nz_bound-1) ;
    assert(xij.get_ndim() == 2) ;
    assert(xijm1.get_ndim() == 2) ;
    assert(ccc.get_ndim() == 2) ;
    assert(xij.get_dim(0) == nt) ;
    assert(xij.get_dim(1) == np2) ;
    assert(xijm1.get_dim(0) == nt) ;
    assert(xijm1.get_dim(1) == np2) ;
    assert(ccc.get_dim(0) == nt) ;
    assert(ccc.get_dim(1) == np2) ;
    
    double Rmax = mp_aff->get_alpha()[nz_bound] + mp_aff->get_beta()[nz_bound] ;
    
    Scalar source_xi = phi ; source_xi.div_r_dzpuis(2) ;
    source_xi -= phi.dsdr() ;
    source_xi.set_spectral_va().ylm() ;
    sphi.set_spectral_va().ylm() ;
    const Base_val& base = sphi.get_spectral_base() ;
    int l_q, m_q, base_r ;
    for (int k=0; k<np2; k++) 
	for (int j=0; j<nt; j++) {
	    base.give_quant_numbers(nz_bound, k, j, m_q, l_q, base_r) ;
	    if (l_q > 1) {
		l_q += dl ;
		double fact = 8*Rmax*Rmax + dt*dt*(6+3*l_q*(l_q+1)) + 12*Rmax*dt ;
		double souphi = -4*dt*dt*l_q*(l_q+1)*
		    source_xi.get_spectral_va().c_cf->val_out_bound_jk(nz_bound, j, k) ;
		double xijp1 = ( 16*Rmax*Rmax*xij(k,j) -
				 (fact - 24*Rmax*dt)*xijm1(k,j) 
				 + souphi) / fact  ;
		ccc.set(k, j) = xijp1 
		    + sphi.get_spectral_va().c_cf->val_out_bound_jk(nz_bound, j, k) ;
		xijm1.set(k,j) = xij(k,j) ;
		xij.set(k,j) = xijp1 ;
	    }
	}

}


