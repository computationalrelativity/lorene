#include<math.h>
// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "proto.h"
#include "diff.h"


namespace Lorene {
// Inversion of the weakly degenerate eliptic operator associatd with spectral quantity tilde(B), with main characteristics
//fit and fit2 (Suitable for tensorial resolution of BH spacetime)
//See also function tilde(laplacian)

 void tensorellipticBt( Scalar source, Scalar& resu, double fit, double fit2, double fit0d2, double fit1d2, double fit0d3, double fit1d3) {


  const int nz = (*source.get_mp().get_mg()).get_nzone(); 	// Number of domains
  int nr = (*source.get_mp().get_mg()).get_nr(1); 	// Number of collocation points in r in each domain
  int nt = (*source.get_mp().get_mg()).get_nt(1); 	// Number of collocation points in theta in each domain
  int np = (*source.get_mp().get_mg()).get_np(1); 	// Number of collocation points in phi in each domain
  const Map_af* map = dynamic_cast<const Map_af*>(&source.get_mp()) ;
  const Mg3d* mgrid = (*map).get_mg();


	// Some helpful stuff...
	
	const Coord& rr = (*map).r ;
	Scalar rrr (*map) ; 
	rrr = rr ; 
	rrr.set_spectral_va().set_base(source.get_spectral_va().base); 

  	Scalar source_coq = source ;
	source_coq.mult_r() ;
	source_coq.mult_r() ;
	source.set_spectral_va().ylm() ;
	source_coq.set_spectral_va().ylm() ;
	Scalar phi(source.get_mp()) ;
	phi.annule_hard() ;
	//	phi.std_spectral_base();
       	phi.set_spectral_va().set_base(source.get_spectral_va().base) ;
	phi.set_spectral_va().ylm() ;	
	Mtbl_cf& sol_coef = (*phi.set_spectral_va().c_cf) ;

	const Base_val& base = source.get_spectral_base() ;
    Mtbl_cf sol_part(mgrid, base) ; sol_part.annule_hard() ;
    Mtbl_cf sol_hom1(mgrid, base) ; sol_hom1.annule_hard() ;
    Mtbl_cf sol_hom2(mgrid, base) ; sol_hom2.annule_hard() ;

	int l_q, m_q, base_r ;

	{ int lz = 0 ;
	    for (int k=0; k < np; k++)
		for (int j=0; j<nt; j++) {
		    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		    if (nullite_plm(j, nt, k, np, base) == 1) {
		      for (int ii=0 ; ii<nr ; ii++){
			    sol_hom1.set(lz, k, j, ii) = 0 ;
			    sol_part.set(lz, k, j, ii) = 0 ;
		      }

		    }
		}
	}


	{ int lz = 1 ;    // The first shell is a really particular case, where the operator is different, and homogeneous solutions have to be handled really carefully.
	    double alpha = (*map).get_alpha()[lz] ;
	    double beta = (*map).get_beta()[lz] ;
	    double ech = beta / alpha ; 
	    for (int k=0; k < np; k++)
	      for (int j=0; j<nt; j++) {
		base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
		if (nullite_plm(j, nt, k, np, base) == 1) {
		  
		  
		  
		  Matrice ope(nr,nr) ;
		  ope.annule_hard() ;
		  
		  Diff_dsdx dx(base_r, nr) ; const Matrice& mdx = dx.get_matrice() ;
		  Diff_dsdx2 dx2(base_r, nr) ; const Matrice& mdx2 = dx2.get_matrice() ;
		  Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;
		  Diff_xdsdx xdx(base_r, nr) ; const Matrice& mxdx = xdx.get_matrice() ;
		  Diff_xdsdx2 xdx2(base_r, nr) ; const Matrice& mxdx2 = xdx2.get_matrice() ;
		  Diff_x2dsdx2 x2dx2(base_r, nr) ; const Matrice& mx2dx2 = x2dx2.get_matrice() ;
		  Diff_x3dsdx2 x3dx2 (base_r, nr); const Matrice& mx3dx2 = x3dx2.get_matrice();
		  Diff_x4dsdx2 x4dx2 (base_r, nr); const Matrice& mx4dx2 = x4dx2.get_matrice();
		  
		  
	
// 		  	ope = mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2 + 2*(mxdx + ech*mdx)
// 		  - l_q*(l_q+1)*mid  - (mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) - (fit)*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) + (fit)*beta*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) + (fit)*alpha*(mx3dx2 + 2*ech*mx2dx2 + ech*ech*mxdx2)  ;
		  

// 		  ope = mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2 + 2*(mxdx + ech*mdx)
// 		    - l_q*(l_q+1)*mid  - ((mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) -(fit)*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) + (fit)*beta*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) 
// 					  + (fit)*alpha*(mx3dx2 + 2*ech*mx2dx2 + ech*ech*mxdx2) + fit2*alpha*alpha*(mx4dx2 + 2*ech*mx3dx2 + ech*ech*mx2dx2) 
// 					  + 2*fit2*alpha*(beta-1.)*(mx3dx2 + 2*ech*mx2dx2 + ech*ech*mxdx2) + fit2*(beta-1.)*(beta-1.)*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2));
		  


		  ope = mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2 + 2*(mxdx + ech*mdx)
		    - l_q*(l_q+1)*mid + 2*l_q*mid - ((mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) -(fit)*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) + (fit)*beta*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) 
					  + (fit)*alpha*(mx3dx2 + 2*ech*mx2dx2 + ech*ech*mxdx2) + fit2*alpha*alpha*(mx4dx2 + 2*ech*mx3dx2 + ech*ech*mx2dx2) 
					  + fit2*alpha*(beta-1.)*(mx3dx2 + 2*ech*mx2dx2 + ech*ech*mxdx2) + fit2*alpha*(beta- rrr.val_grid_point(1,0,0, nr-1))*(mx3dx2 + 2*ech*mx2dx2 + ech*ech*mxdx2)+ fit2*(beta-1.)*(beta- rrr.val_grid_point(1,0,0,nr -1))*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2));
		  



		  for (int col=0; col<nr; col++)
		    ope.set(nr-1, col) = 0 ;
		  ope.set(nr-1, 0) = 1 ;

		  
		  Tbl rhs(nr); 
		  rhs.annule_hard() ;
		  for (int i=0; i<nr; i++)
		    rhs.set(i) = (*source_coq.get_spectral_va().c_cf)(1, k, j, i) ;
		  rhs.set(nr-1) = 0 ;
		  ope.set_lu() ;
		  Tbl sol = ope.inverse(rhs) ;

		  
		  for (int i=0; i<nr; i++)
		    sol_part.set(1, k, j, i) = sol(i) ;
		  
		  rhs.annule_hard();
		  rhs.set(nr-1) = 1 ;
		  sol = ope.inverse(rhs) ;


		  for (int i=0; i<nr; i++)
		    sol_hom1.set(1, k, j, i) = sol(i) ;              
		  
		  
		}
	      }
	}
	
	// Attention! zones 2 et 3traitee separement egalement!!!

	
	// Current implementations only allow grids with more than 3 shells.
	
          { int lz = 2 ;

	    double alpha = (*map).get_alpha()[lz] ;
	    double beta = (*map).get_beta()[lz] ;
	    double ech = beta / alpha ; 
	
	  for (int k=0; k < np; k++)
	    for (int j=0; j<nt; j++) {
	      base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	      if (nullite_plm(j, nt, k, np, base) == 1) {
		
		Matrice ope(nr,nr) ;
		ope.annule_hard() ;
		
		Diff_dsdx dx(base_r, nr) ; const Matrice& mdx = dx.get_matrice() ;
		Diff_dsdx2 dx2(base_r, nr) ; const Matrice& mdx2 = dx2.get_matrice() ;
		Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;
		Diff_xdsdx xdx(base_r, nr) ; const Matrice& mxdx = xdx.get_matrice() ;
		Diff_xdsdx2 xdx2(base_r, nr) ; const Matrice& mxdx2 = xdx2.get_matrice() ;
		Diff_x2dsdx2 x2dx2(base_r, nr) ; const Matrice& mx2dx2 = x2dx2.get_matrice() ;
		Diff_x3dsdx2 x3dx2 (base_r, nr); const Matrice& mx3dx2 = x3dx2.get_matrice();
	
			ope = mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2 + 2*(mxdx + ech*mdx)
		- l_q*(l_q+1)*mid + 2*l_q*mid - fit0d2*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) + (fit1d2)*(rrr.val_grid_point(lz, 0, 0, 0))*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) - (fit1d2)*beta*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) - (fit1d2)*alpha*(mx3dx2 + 2*ech*mx2dx2 + ech*ech*mxdx2)  ;



// 				ope = mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2 + 2*(mxdx + ech*mdx)
// 				  - l_q*(l_q+1)*mid ;
		
		for (int col=0; col<nr; col++)
		  ope.set(nr-1, col) = 0 ;
		ope.set(nr-1, 0) = 1 ;
		for (int col=0; col<nr; col++) {
		  ope.set(nr-2, col) = 0 ;
		}
		ope.set(nr-2, 1) = 1 ;
	
		Tbl rhs(nr) ;
		rhs.annule_hard() ;
		for (int i=0; i<nr; i++)
		  rhs.set(i) = (*source_coq.get_spectral_va().c_cf)(lz, k, j, i) ;
		rhs.set(nr-2) = 0 ;
		rhs.set(nr-1) = 0 ;
		ope.set_lu() ;
		Tbl sol = ope.inverse(rhs) ;

		for (int i=0; i<nr; i++)
			    sol_part.set(lz, k, j, i) = sol(i) ;
		
		rhs.annule_hard() ;
		rhs.set(nr-2) = 1 ;
		sol = ope.inverse(rhs) ;
		for (int i=0; i<nr; i++)
		  sol_hom1.set(lz, k, j, i) = sol(i) ;
		
		rhs.set(nr-2) = 0 ;
		rhs.set(nr-1) = 1 ;
		sol = ope.inverse(rhs) ;
		for (int i=0; i<nr; i++)
		  sol_hom2.set(lz, k, j, i) = sol(i) ;
		
	      }
	    }
	  
	  }	


	
          { int lz = 3 ;

	    double alpha = (*map).get_alpha()[lz] ;
	    double beta = (*map).get_beta()[lz] ;
	    double ech = beta / alpha ; 
	
	  for (int k=0; k < np; k++)
	    for (int j=0; j<nt; j++) {
	      base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	      if (nullite_plm(j, nt, k, np, base) == 1) {
		
		Matrice ope(nr,nr) ;
		ope.annule_hard() ;
		
		Diff_dsdx dx(base_r, nr) ; const Matrice& mdx = dx.get_matrice() ;
		Diff_dsdx2 dx2(base_r, nr) ; const Matrice& mdx2 = dx2.get_matrice() ;
		Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;
		Diff_xdsdx xdx(base_r, nr) ; const Matrice& mxdx = xdx.get_matrice() ;
		Diff_xdsdx2 xdx2(base_r, nr) ; const Matrice& mxdx2 = xdx2.get_matrice() ;
		Diff_x2dsdx2 x2dx2(base_r, nr) ; const Matrice& mx2dx2 = x2dx2.get_matrice() ;
		Diff_x3dsdx2 x3dx2 (base_r, nr); const Matrice& mx3dx2 = x3dx2.get_matrice();
	
			ope = mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2 + 2*(mxdx + ech*mdx)
		- l_q*(l_q+1)*mid  +2*l_q*mid - fit0d3*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) + (fit1d3)*(rrr.val_grid_point(lz, 0, 0, 0))*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) - (fit1d3)*beta*(mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2) - (fit1d3)*alpha*(mx3dx2 + 2*ech*mx2dx2 + ech*ech*mxdx2)  ;



// 				ope = mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2 + 2*(mxdx + ech*mdx)
// 				  - l_q*(l_q+1)*mid ;
		
		for (int col=0; col<nr; col++)
		  ope.set(nr-1, col) = 0 ;
		ope.set(nr-1, 0) = 1 ;
		for (int col=0; col<nr; col++) {
		  ope.set(nr-2, col) = 0 ;
		}
		ope.set(nr-2, 1) = 1 ;
	
		Tbl rhs(nr) ;
		rhs.annule_hard() ;
		for (int i=0; i<nr; i++)
		  rhs.set(i) = (*source_coq.get_spectral_va().c_cf)(lz, k, j, i) ;
		rhs.set(nr-2) = 0 ;
		rhs.set(nr-1) = 0 ;
		ope.set_lu() ;
		Tbl sol = ope.inverse(rhs) ;

		for (int i=0; i<nr; i++)
			    sol_part.set(lz, k, j, i) = sol(i) ;
		
		rhs.annule_hard() ;
		rhs.set(nr-2) = 1 ;
		sol = ope.inverse(rhs) ;
		for (int i=0; i<nr; i++)
		  sol_hom1.set(lz, k, j, i) = sol(i) ;
		
		rhs.set(nr-2) = 0 ;
		rhs.set(nr-1) = 1 ;
		sol = ope.inverse(rhs) ;
		for (int i=0; i<nr; i++)
		  sol_hom2.set(lz, k, j, i) = sol(i) ;
		
	      }
	    }
	  
	  }	





	
	
	// Current implementations only allow grids with more than 2 shells.
	
	for (int lz=4; lz<nz-1; lz++) {
	  double ech = (*map).get_beta()[lz] / (*map).get_alpha()[lz] ;
	  for (int k=0; k < np; k++)
	    for (int j=0; j<nt; j++) {
	      base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	      if (nullite_plm(j, nt, k, np, base) == 1) {
		
		Matrice ope(nr,nr) ;
		ope.annule_hard() ;
		
		Diff_dsdx dx(base_r, nr) ; const Matrice& mdx = dx.get_matrice() ;
		Diff_dsdx2 dx2(base_r, nr) ; const Matrice& mdx2 = dx2.get_matrice() ;
		Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;
		Diff_xdsdx xdx(base_r, nr) ; const Matrice& mxdx = xdx.get_matrice() ;
		Diff_xdsdx2 xdx2(base_r, nr) ; const Matrice& mxdx2 = xdx2.get_matrice() ;
		Diff_x2dsdx2 x2dx2(base_r, nr) ; const Matrice& mx2dx2 = x2dx2.get_matrice() ;
		ope = mx2dx2 + 2*ech*mxdx2 + ech*ech*mdx2 + 2*(mxdx + ech*mdx)
		  - l_q*(l_q+1)*mid +2*l_q*mid ;
		
		for (int col=0; col<nr; col++)
		  ope.set(nr-1, col) = 0 ;
		ope.set(nr-1, 0) = 1 ;
		for (int col=0; col<nr; col++) {
		  ope.set(nr-2, col) = 0 ;
		}
		ope.set(nr-2, 1) = 1 ;
	
		Tbl rhs(nr) ;
		rhs.annule_hard() ;
		for (int i=0; i<nr; i++)
		  rhs.set(i) = (*source_coq.get_spectral_va().c_cf)(lz, k, j, i) ;
		rhs.set(nr-2) = 0 ;
		rhs.set(nr-1) = 0 ;
		ope.set_lu() ;
		Tbl sol = ope.inverse(rhs) ;

		for (int i=0; i<nr; i++)
			    sol_part.set(lz, k, j, i) = sol(i) ;
		
		rhs.annule_hard() ;
		rhs.set(nr-2) = 1 ;
		sol = ope.inverse(rhs) ;
		for (int i=0; i<nr; i++)
		  sol_hom1.set(lz, k, j, i) = sol(i) ;
		
		rhs.set(nr-2) = 0 ;
		rhs.set(nr-1) = 1 ;
		sol = ope.inverse(rhs) ;
		for (int i=0; i<nr; i++)
		  sol_hom2.set(lz, k, j, i) = sol(i) ;
		
	      }
	    }
	  
	}
	{ int lz = nz-1 ;
	double alpha = (*map).get_alpha()[lz] ;
	for (int k=0; k < np; k++)
	  for (int j=0; j<nt; j++) {
	    base.give_quant_numbers(lz, k, j, m_q, l_q, base_r) ;
	    if (nullite_plm(j, nt, k, np, base) == 1) {
	      
	      Matrice ope(nr,nr) ;
	      ope.annule_hard() ;
	      Diff_dsdx2 dx2(base_r, nr) ; const Matrice& mdx2 = dx2.get_matrice() ;
	      Diff_sx2 sx2(base_r, nr) ; const Matrice& ms2 = sx2.get_matrice() ;
	      
			ope = (mdx2 - l_q*(l_q+1)*ms2 + 2*l_q*ms2)/(alpha*alpha) ;
			
			for (int i=0; i<nr; i++)
			  ope.set(nr-1, i) = 0 ;
			ope.set(nr-1, 0) = 1 ; //for the true homogeneous solution
			
			for (int i=0; i<nr; i++) {
			  ope.set(nr-2, i) = 1 ; //for the limit at inifinity
			}

			if (l_q > 0) {
			    for (int i=0; i<nr; i++) {
				ope.set(nr-3, i) = i*i ; //for the finite part (derivative = 0 at infty)
			    }
			}
			//	cout << "l: " << l_q << endl ;

			Tbl rhs(nr) ;
			rhs.annule_hard() ;
			for (int i=0; i<nr; i++)
			    rhs.set(i) = (*source.get_spectral_va().c_cf)(lz, k, j, i) ;
			if (l_q>0) rhs.set(nr-3) = 0 ;
			rhs.set(nr-2) = 0 ;
			rhs.set(nr-1) = 0 ;
			ope.set_lu() ;
			Tbl sol = ope.inverse(rhs) ;
		       
			for (int i=0; i<nr; i++)
			    sol_part.set(lz, k, j, i) = sol(i) ;

			rhs.annule_hard() ;
			rhs.set(nr-1) = 1 ;
			sol = ope.inverse(rhs) ;
			for (int i=0; i<nr; i++)
			    sol_hom2.set(lz, k, j, i) = sol(i) ;

		    }
		}
	}

	Mtbl_cf dpart = sol_part ; dpart.dsdx() ;
	Mtbl_cf dhom1 = sol_hom1 ; dhom1.dsdx() ;
	Mtbl_cf dhom2 = sol_hom2 ; dhom2.dsdx() ;


	// Now matching the homogeneous solutions between the different domains...
	
	for (int k=0; k < np; k++)
	    for (int j=0; j<nt; j++) {
		base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
		if (nullite_plm(j, nt, k, np, base) == 1) {
		    Matrice systeme(2*nz-4, 2*nz-4) ;
		    systeme.annule_hard() ;
		    Tbl rhs(2*nz-4) ;
		    rhs.annule_hard() ;

		    //First shell 
		    int lin = 0 ;
		    int col = 0 ;
		    
		    double alpha = (*map).get_alpha()[1] ;
		    
		    systeme.set(lin, col) += sol_hom1.val_out_bound_jk(1, j, k) ;
		    rhs.set(lin) -= sol_part.val_out_bound_jk(1, j, k) ;
		    
		    lin++ ;
		    systeme.set(lin, col) += dhom1.val_out_bound_jk(1, j, k) / alpha ;
		    rhs.set(lin) -=  dpart.val_out_bound_jk(1, j, k) / alpha ;
		    col += 1 ;
		    
		    //Shells
		    for (int lz=2; lz<nz-1; lz++) {
		      alpha = (*map).get_alpha()[lz] ;
		      lin-- ;
		      systeme.set(lin,col) -= sol_hom1.val_in_bound_jk(lz, j, k) ;
		      systeme.set(lin,col+1) -= sol_hom2.val_in_bound_jk(lz, j, k) ;
		      rhs.set(lin) += sol_part.val_in_bound_jk(lz, j, k) ;
		      
		      lin++ ;
		      systeme.set(lin,col) -= dhom1.val_in_bound_jk(lz, j, k) / alpha ;
		      systeme.set(lin,col+1) -= dhom2.val_in_bound_jk(lz, j, k) / alpha ;
		      rhs.set(lin) += dpart.val_in_bound_jk(lz, j, k) / alpha;
		      
		      lin++ ;
		      systeme.set(lin, col) += sol_hom1.val_out_bound_jk(lz, j, k) ;
		      systeme.set(lin, col+1) += sol_hom2.val_out_bound_jk(lz, j, k) ;
		      rhs.set(lin) -= sol_part.val_out_bound_jk(lz, j, k) ;
		      
		      lin++ ;
		      systeme.set(lin, col) += dhom1.val_out_bound_jk(lz, j, k) / alpha ;
		      systeme.set(lin, col+1) += dhom2.val_out_bound_jk(lz, j, k) / alpha ;
		      rhs.set(lin) -= dpart.val_out_bound_jk(lz, j, k) / alpha ;
		      col += 2 ;
		    }
		    
		    //CED
		    alpha = (*map).get_alpha()[nz-1] ;
		    lin-- ;
		    systeme.set(lin,col) -= sol_hom2.val_in_bound_jk(nz-1, j, k) ;
		    rhs.set(lin) += sol_part.val_in_bound_jk(nz-1, j, k) ;

		    lin++ ;
		    systeme.set(lin,col) -= (-4*alpha)*dhom2.val_in_bound_jk(nz-1, j, k) ;
		    rhs.set(lin) += (-4*alpha)*dpart.val_in_bound_jk(nz-1, j, k) ;

		    systeme.set_lu() ;

		    // cout << systeme << endl;

		    Tbl coef = systeme.inverse(rhs);
		    int indice = 0 ;
		    		    
		    //	    int tryr; cin >> tryr;
		    // 		    for (int i=0; i<mgrid.get_nr(0); i++)
		    // 			sol_coef.set(0, k, j, i) = 0 ;
		    // 		    sol_coef.set(0, k, j, 0) =  (*bound.get_spectral_va().c_cf)(0, k, j, 0) ;
		    
		    for (int i=0; i<(*mgrid).get_nr(1); i++)
		      sol_coef.set(1, k, j, i) = sol_part(1, k, j, i)
			+coef(indice)*sol_hom1(1, k, j, i) ;
		    indice +=1;
		    
		    
		    for (int lz=2; lz<nz-1; lz++) {
		      for (int i=0; i<(*mgrid).get_nr(lz); i++)
			sol_coef.set(lz, k, j, i) = sol_part(lz, k, j, i)
			  +coef(indice)*sol_hom1(lz, k, j, i) 
			  +coef(indice+1)*sol_hom2(lz, k, j, i) ;
		      indice += 2 ;
		    }
		    for (int i=0; i<(*mgrid).get_nr(nz-1); i++)
		      sol_coef.set(nz-1, k, j, i) = sol_part(nz-1, k, j, i)
			+coef(indice)*sol_hom2(nz-1, k, j, i) ;

		}
	    }
	
	
	delete phi.set_spectral_va().c ;
	phi.set_spectral_va().c = 0x0 ;
	//	phi.set_spectral_va().ylm_i() ;
	
	phi.annule_domain(nz-1);
	
	resu = phi;

}
}
