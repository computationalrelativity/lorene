/*
 *  Definition of Lorene classes Ope_elementary
 *
 */

/*
 *   Copyright (c) 2003 Philippe Grandclement
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

#ifndef __OPE_ELEMENTARY_H_ 
#define __OPE_ELEMENTARY_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/12/11 14:57:00  p_grandclement
 * I had forgotten the .h (sorry folks...)
 *
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 *
 * $Header$
 *
 */

#include "matrice.h"

/**
 * Basic class for elementary elliptic operators.
 *
 * Such objects describe a type of elliptic operator, in a given domain and 
 * for a given spherical harmonics.
 * They are called by the general elliptic solver. {\tt Ope_elementary} 
 * objects 
 * know how to compute the approriate particular solution for a given source 
 * and the associated homogeneous ones.
 * 
 * The class {\tt Ope_elementary} is an abstract one: 
 * it cannot be instanciated. 
 * Specific implementation of coordinate mappings will be performed by derived
 * classes of {\tt Ope_elementary}. 
 * 
 *  @version #$Id:#
 **/

class Ope_elementary {

 protected:

  int nr ; /// Number of radial points
  int base_r ; /// Radial basis of decomposition
  double alpha ; /// Parameter $\alpha$ of the associated mapping.
  double beta ; /// Parameter $\beta$ of the associated mapping.
  
  /**
   * Pointer on the matrix representation of the operator.
   **/
  mutable Matrice* ope_mat ;
  /**
   * Pointer on the banded-matrix of the operator.
   **/
  mutable Matrice* ope_cl ;
  /**
   * Pointer on the non-degenerated matrix of the operator.
   **/
  mutable Matrice* non_dege ;

  /**
   * Value of the first homogeneous solution at the outer boundary.
   **/
  mutable double s_one_plus ;
  /**
   * Value of the first homogeneous solution at the inner boundary.
   **/
  mutable double s_one_minus ;
  /**
   * Value of the derivative  of the 
   * first homogeneous solution at the outer boundary.
   **/
  mutable double ds_one_plus ;
  /**
   * Value of the derivative  of the 
   * first homogeneous solution at the inner boundary.
   **/
  mutable double ds_one_minus ;
  
  /**
   * Value of the second homogeneous solution at the outer boundary.
   **/
  mutable double s_two_plus ;
  /**
   * Value of the second homogeneous solution at the inner boundary.
   **/
  mutable double s_two_minus ;
  /**
   * Value of the derivative  of the 
   * second homogeneous solution at the outer boundary.
   **/
  mutable double ds_two_plus ;
  /**
   * Value of the derivative  of the 
   * second homogeneous solution at the inner boundary.
   **/
  mutable double ds_two_minus ;
 
  /**
   * Value of the particular solution at the inner boundary.
   **/
  mutable double sp_minus ;
  /**
   * Value of the particular solution at the outer boundary.
   **/
  mutable double sp_plus ;
  /**
   * Value of the derivative of the particular solution at the inner boundary.
   **/
  mutable double dsp_minus ;
  /**
   * Value of the derivative of the particular solution at the outer boundary.
   **/
  mutable double dsp_plus ;

  // Constructors destructor
 protected:
  /**
   * Standard constructor, protected because the class is an abstract one.
   *
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter $\alpha$ of the mapping.
   * @param bet [input] parameter $\beta$ of the mapping.
   **/
  explicit Ope_elementary (int nbr , int baser , double alf, double eta) ;
  Ope_elementary (const Ope_elementary&) ; /// Constructor by copy

 public:
  virtual ~Ope_elementary() ; /// Destructor

 public: 
  /**
   * Returns the value of the first homogeneous solution at the inner
   * boundary.
   **/
  double val_sh_one_minus() const {return s_one_minus ;} ; 
  /**
   * Returns the value of the first homogeneous solution at the outer
   * boundary.
   **/
  double val_sh_one_plus() const {return s_one_plus ;} ; 
  /**
   * Returns the value of the derivative of the 
   * first homogeneous solution at the inner
   * boundary.
   **/
  double der_sh_one_minus() const {return ds_one_minus ;} ; 
  /**
   * Returns the value of the derivative of the 
   * first homogeneous solution at the outer
   * boundary.
   **/
  double der_sh_one_plus() const {return ds_one_plus ;} ;
 
  /**
   * Returns the value of the second homogeneous solution at the inner
   * boundary.
   **/
  double val_sh_two_minus() const {return s_two_minus ;} ;
  /**
   * Returns the value of the second homogeneous solution at the outer
   * boundary.
   **/
  double val_sh_two_plus() const {return s_two_plus ;} ;
  /**
   * Returns the value of the derivative of the 
   * second homogeneous solution at the inner
   * boundary.
   **/
  double der_sh_two_minus() const {return ds_two_minus ;} ;
  /**
   * Returns the value of the derivative of the 
   * second homogeneous solution at the outer
   * boundary.
   **/
  double der_sh_two_plus() const {return ds_two_plus ;} ;

  /**
   * Returns the value of the particular solution at the inner boundary.
   **/
  double val_sp_minus() const {return sp_minus ;} ;
  /**
   * Returns the value of the particular solution at the outer boundary.
   **/
  double val_sp_plus() const {return sp_plus ;} ; 
  /**
   * Returns the value of the derivative 
   * particular solution at the inner boundary.
   **/
  double der_sp_minus() const {return dsp_minus ;} ;
  /**
   * Returns the value of the derivative 
   * particular solution at the outer boundary.
   **/
  double der_sp_plus() const {return dsp_plus ;} ;

  double get_alpha() const {return alpha ;} ; /// Returns {\tt alpha}.
  double get_beta() const {return beta ;} ; /// Returns {\tt beta}.
  int get_base_r() const {return base_r ;} ; /// Returns {\tt base_r}.

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const = 0 ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const = 0 ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const = 0 ;  
 
 public:
  /**
   * Computes the particular solution, given the source {\tt so}.
   **/
  virtual Tbl get_solp(const Tbl& so) const = 0 ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const = 0 ;
  /**
   * Increases the quatum number $l$ by one unit.
   **/
  virtual void inc_l_quant() = 0 ;

} ;

/**
 * Class for the operator of the Poisson equation (i.e. the Laplacian !).
 *
 * It is implemented in every type of domains.
 **/

class Ope_poisson : public Ope_elementary {

 protected:
  int l_quant ; /// quantum number
  int dzpuis ; /// the associated dzpuis, if in the compactified domain. 
  
 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter $\alpha$ of the mapping.
   * @param bet [input] parameter $\beta$ of the mapping.
   * @param lq [input] quantum number $l$.
   * @param dz [input] dzpuis of the source.
   **/
  Ope_poisson (int nbr, int baser, double alf, double bet, int lq, int dz) ;
  Ope_poisson (const Ope_poisson&) ; /// Constructor by copy
  virtual ~Ope_poisson() ; /// Destructor

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source {\tt so}.
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number $l$ by one unit.
   **/
  virtual void inc_l_quant() ;
} ;

/**
 * Class for the Helmholtz operator $\Delta - m^2$ ($m > 0$).
 * 
 * It is implemented only in the shells and in the compactified domain.
 **/
class Ope_helmholtz_minus : public Ope_elementary {

 protected:
  double masse ; /// The mass parameter $m$.

 public:
   /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter $\alpha$ of the mapping.
   * @param bet [input] parameter $\beta$ of the mapping.
   * @param mas [input] mass parameter $m$.
   **/
  Ope_helmholtz_minus (int nbr, int baser, double alf, double bet, 
		       double mas) ;
  Ope_helmholtz_minus (const Ope_helmholtz_minus&) ; /// Constructor by copy
  virtual ~Ope_helmholtz_minus() ; /// Destructor
  
 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source {\tt so}.
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number $l$ by one unit (CURRENTLY NOT IMPLEMENTED)
   **/
  virtual void inc_l_quant() ;
} ;

/**
 * Class for the Helmholtz operator $\Delta + m^2$ ($m > 0$).
 * 
 * It is implemented only in the shells.
 **/
class Ope_helmholtz_plus : public Ope_elementary {

 protected:
  double masse ; /// The mass parameter $m$.

 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter $\alpha$ of the mapping.
   * @param bet [input] parameter $\beta$ of the mapping.
   * @param mas [input] mass parameter $m$.
   **/
  Ope_helmholtz_plus (int nbr, int baser, double alf, double bet, 
		      double mas) ;
  Ope_helmholtz_plus (const Ope_helmholtz_plus&) ; /// Constructor by copy
  virtual ~Ope_helmholtz_plus() ; /// Destructor
  
 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source {\tt so}.
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number $l$ by one unit (CURRENTLY NOT IMPLEMENTED)
   **/
  virtual void inc_l_quant() ;
} ;

# endif
