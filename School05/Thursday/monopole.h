#ifndef __MONOPOLE_H_ 
#define __MONOPOLE_H_ 

#include "map.h"
#include "scalar.h"

using namespace Lorene ;

class Monopole {

 protected:
  const Map_af& mp ; /// Mapping 

  // radius
  Scalar radius ; /// Scalar containing the radius
  int nz ; ///Number of domains
  
  // DATA : 
  double beta ; /// beta parameter
  
  Scalar big_W ; /// Function W
  Scalar small_w ; /// Representation of W in terms of w
  
  Scalar big_H ; /// Function H
  Scalar small_h ; /// Representation of h in terms of h;
  
  // Constructors destructor
 public:
  explicit Monopole (const Map_af&, double) ; /// Standard constructor 
  Monopole (const Monopole&) ; /// Constructor by copy
  
  virtual ~Monopole() ; /// Destructor

  // Assigment
 public:
  void operator= (const Monopole&) ; /// Assignment to another monopole

  // Accessor :
 public:
  const Map_af& get_mp() const {return mp ;} ; /// Return the mapping
  
  Scalar get_big_W() const {return big_W ;} ; /// Return the field W
  
  Scalar get_small_w() const {return small_w ;} ; /// Return the field w
  
  Scalar get_big_H() const {return big_H ;} ; /// Return the field H
 
  Scalar get_small_h() const {return small_h ;} ; /// Return the field h
    
  // The numerical methods (at last) :
 public:
  void init_big_W() ; /// Initialize big W to an initial guess
  
  void init_big_H() ; /// Initialize big H to an initial guess
  
  /*
   * Computes the auxiliary variable w.
   */
  void do_small_w() ;

  /*
   * Computes big_W from small_w
   */
  void do_big_W() ;

  /*
   * Computes the auxiliary variable h.
   */
  void do_small_h() ;

  /*
   * Computes big_H from small_h
   */
  void do_big_H() ;

  /*
   * Compute the source of the equation for W , using either the small_w version 
   * or the W one, depending on the region
   */

  Scalar compute_source_W() const ;
 
  /*
   * Compute the source of the equation for H , using either the small_h version 
   * or the H one, depending on the region
   */

  Scalar compute_source_H() const ;
 
 public:
  /**
   * Solve the system by iteration : 
   * prec : required precision.
   * itemax : number of iterations.
   * relax : relaxation parameter
   * returns the actual number of iterations performed
   **/
  int solve_config (double prec = 1e-10, int itemax= 200, double relax=0.5) ;

  
  double give_a () const ; /// Computes a
  double give_b() const ; /// Computes b
} ;

# endif
