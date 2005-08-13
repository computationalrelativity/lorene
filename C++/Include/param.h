/*
 *  Definition of Lorene class Param
 *
 */

/*
 *   Copyright (c) 1999-2005 Eric Gourgoulhon
 *   Copyright (c) 2000-2003 Jerome Novak
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


#ifndef __PARAM_H_ 
#define __PARAM_H_ 


/*
 * $Id$
 * $Log$
 * Revision 1.7  2005/08/13 16:08:20  m_saijo
 * Corrected the documents related to the Star
 *
 * Revision 1.6  2005/08/13 16:03:36  m_saijo
 * Added storage of a Star
 *
 * Revision 1.5  2005/03/24 21:55:58  e_gourgoulhon
 * Added storage of a Scalar.
 *
 * Revision 1.4  2004/03/22 13:12:43  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.3  2003/09/25 12:08:02  j_novak
 * Tensors can be stored in Param objects
 *
 * Revision 1.2  2002/09/19 09:52:42  j_novak
 * Added objects Qtenseur and Qmetrique for 4D tensor and metric handling.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.10  2001/10/27  09:26:24  novak
 * *** empty log message ***
 *
 * Revision 1.9  2001/10/11 07:44:12  eric
 * Ajout du stokage des Etoile's
 *
 * Revision 1.8  2000/10/24  14:54:49  novak
 * Added the function clean_all()
 *
 * Revision 1.7  2000/05/25 12:39:19  eric
 * MODIFICATION MAJEURE: pour les int et les double, ce sont desormais les
 * adresses qui sont stokees, et non plus les nombres eux-memes
 * (le traitement des int et des double est donc desormais completement
 * aligne sur celui des Tbl, Cmp, etc...)
 *
 * Revision 1.6  1999/12/29  13:10:39  eric
 * Ajout du stokage des Mtbl_cf.
 *
 * Revision 1.5  1999/12/27  12:16:43  eric
 * Ajout du stokage des mappings (class Map).
 *
 * Revision 1.4  1999/12/16  10:27:40  eric
 * Ajout des membres modifiables.
 * Par defaut, les objets listes sont const.
 *
 * Revision 1.3  1999/12/15  16:49:52  eric
 * *** empty log message ***
 *
 * Revision 1.2  1999/12/15  16:22:36  eric
 * Changement de l'ordre des arguments dans add_*
 * Argument par defaut: position = 0
 * Ajout du stokage des int et des double.
 *
 * Revision 1.1  1999/12/13  14:35:56  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

class Tbl ; 
class Itbl ; 
class Mtbl_cf ; 
class Map ; 
class Cmp ; 
class Tenseur ;
class Qtenseur ;
class Tensor ;
class Scalar ; 
class Etoile ;
class Star ;

/** Parameter storage.
 * 
 *  This class is intended to store addresses of various Lorene objects to
 *  pass them as parameters in some subroutines. 
 * \ingroup (util)
 * 
 */
class Param {

    // Data : 
    // -----
    private:
	int n_int ;	///< Number of \c int 's (integers). 
	/// Array (size \c n_int ) of the \c int 's addresses.
	const int** p_int ;	

	int n_int_mod ;	///< Number of modifiable \c int 's (integers). 
	/// Array (size \c n_int_mod ) of the modifiable \c int 's addresses.
	int** p_int_mod ;	

	int n_double ; ///< Number of \c double 's (double precis. numbers). 
	/// Array (size \c n_double ) of the \c double 's addresses.
	const double** p_double ; 

	int n_double_mod ; ///< Number of modifiable \c double 's (double precis. numbers). 
	/// Array (size \c n_double_mod ) of the \c double 's addresses
	double** p_double_mod ; 


	int n_tbl ;	///< Number of \c Tbl 's 
	/// Array (size \c n_tbl ) of the \c Tbl 's addresses
	const Tbl** p_tbl ;	
    
	int n_tbl_mod ;	///< Number of modifiable \c Tbl 's 
	/// Array (size \c n_tbl_mod ) of the modifiable \c Tbl 's addresses
	Tbl** p_tbl_mod ;	
    
	int n_itbl ;	///< Number of \c Itbl 's 
	/// Array (size \c n_itbl ) of the \c Itbl 's addresses
	const Itbl** p_itbl ;	
    
	int n_itbl_mod ;    ///< Number of modifiable \c Itbl 's 
	/// Array (size \c n_itbl_mod ) of the modifiable \c Itbl 's addresses
	Itbl** p_itbl_mod ;	
    
	int n_cmp ;	///< Number of \c Cmp 's 
	/// Array (size \c n_cmp ) of the \c Cmp 's addresses
	const Cmp** p_cmp ;	
    
	int n_cmp_mod ;	///< Number of modifiable \c Cmp 's 
	/// Array (size \c n_cmp_mod ) of the modifiable \c Cmp 's addresses
	Cmp** p_cmp_mod ;	
    
	int n_tenseur ;	///< Number of \c Tenseur 's 
	/// Array (size \c n_tenseur ) of the \c Tenseur 's addresses
	const Tenseur** p_tenseur ;	
    
	int n_tenseur_mod ;	///< Number of modifiable \c Tenseur 's 
	/// Array (size \c n_tenseur_mod ) of the modifiable \c Tenseur 's addresses
	Tenseur** p_tenseur_mod ;	

	int n_qtenseur ;	///< Number of \c Qtenseur 's 
	/// Array (size \c n_qtenseur ) of the \c Qtenseur 's addresses
	const Qtenseur** p_qtenseur ;	
    
	int n_qtenseur_mod ;	///< Number of modifiable \c Qtenseur 's 
	/// Array (size \c n_qtenseur_mod ) of the modifiable \c Qtenseur 's addresses
	Qtenseur** p_qtenseur_mod ;	

	int n_map ;	///< Number of \c Map 's 
	/// Array (size \c n_map ) of the \c Map 's addresses
	const Map** p_map ;	
    
	int n_mtbl_cf ;	///< Number of \c Mtbl_cf 's 
	/// Array (size \c n_mtbl_cf ) of the \c Mtbl_cf 's addresses
	const Mtbl_cf** p_mtbl_cf ;	
    
	int n_scalar ;	///< Number of \c Scalar 's 
	/// Array (size \c n_scalar ) of the \c Scalar 's addresses
	const Scalar** p_scalar ;	
    
	int n_scalar_mod ;	///< Number of modifiable \c Scalar 's 
	/// Array (size \c n_scalar_mod ) of the modifiable \c Scalar 's addresses
	Scalar** p_scalar_mod ;
	
	int n_tensor ;	///< Number of \c Tensor 's 
	/// Array (size \c n_tensor ) of the \c Tensor 's addresses
	const Tensor** p_tensor ;	
    
	int n_tensor_mod ;	///< Number of modifiable \c Tensor 's 
	/// Array (size \c n_tensor_mod ) of the modifiable \c Tensor 's addresses
	Tensor** p_tensor_mod ;
	
	int n_etoile ;	///< Number of \c Etoile 's
	/// Array (size \c n_etoile ) of the \c Etoile 's addresses
	const Etoile** p_etoile ;	

	int n_star ;	///< Number of \c Star 's
	/// Array (size \c n_star ) of the \c Star 's addresses
	const Star** p_star ;	

    // Constructors - Destructor
    // -------------------------
	
    public:
	Param() ;	///< Default constructor is the only constructor

    private:
	/** Copy constructor (private and not implemented to make \c Param 
	 * a non-copyable class)
	 */ 
	Param(const Param& ) ;
	
    public: 
	~Param() ;	///< Destructor

	/** Deletes all the objects stored as modifiables,
	 *  i.e. all quantities with the suffix \c mod .
	 */
	void clean_all() ; 	



    // Assignment
    // -----------
    private: 
	/** Assignment operator (private and not implemented to make 
	 *   \c Param  a non-copyable class)
	 */
	void operator=(const Param& ) ;
	 	

    // Addition/Extraction of one element
    // ----------------------------------
    public:
    
	///Returns the number of stored \c int 's addresses.
	int get_n_int() const ; 
    
	/** Adds the address of a new \c int  to the list.
	 * 
	 *  @param n [input] \c int  the address of which is 
	 *                             to be stored
	 *  @param position [input] position of the \c int  in the list
	 *			    of stored \c int  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_int(const int& n, int position = 0) ;
	
	/** Returns the reference of a \c int  stored in the list.
	 * 
	 *  @param position [input] position of the \c int  in the list
	 *			    of stored \c int  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c int  the address of 
	 *           which is stored at the location  \c position  in the 
	 *           list
	 */
	const int& get_int(int position = 0) const; 

	///Returns the number of modifiable \c int 's addresses in the list.
	int get_n_int_mod() const ; 
    
	/** Adds the address of a new modifiable \c int  to the list.
	 * 
	 *  @param n [input] modifiable \c int  the address of which is 
	 *                   to be stored
	 *  @param position [input] position of the modifiable \c int  
	 *                in the list of stored modifiable \c int  addresses
	 *                (default value = 0)
	 * 
	 */
	void add_int_mod(int& n, int position = 0) ;
	
	/** Returns the reference of a modifiable \c int  stored in the list.
	 * 
	 *  @param position [input] position of the modifiable \c int  
	 *              in the list of stored modifiable \c int  addresses 
	 *                    (default value = 0)
	 *  @return Reference to the modifiable \c int  the address of 
	 *           which is stored at the location  \c position  in the 
	 *           list
	 */
	int& get_int_mod(int position = 0) const; 
	

	///Returns the number of stored \c double 's addresses.
	int get_n_double() const ; 
    
	/** Adds the the address of a new \c double  to the list.
	 * 
	 *  @param x [input] \c double  the address of which is 
	 *                             to be stored
	 *  @param position [input] position of the \c double  in the list
	 *			    of stored \c double  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_double(const double& x, int position = 0) ;
	
	/** Returns the reference of a \c double  stored in the list.
	 * 
	 *  @param position [input] position of the \c double  in the list
	 *			    of stored \c double  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c double  the address of 
	 *           which is stored at the location  \c position  in the 
	 *           list
	 */
	const double& get_double(int position = 0) const; 
	

	///Returns the number of stored modifiable \c double 's addresses.
	int get_n_double_mod() const ; 
    
	/** Adds the address of a new modifiable \c double  to the list.
	 * 
	 *  @param x [input] modifiable \c double  the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the modifiable \c double  
	 *                          in the list of stored modifiable 
	 *                      \c double  addresses (default value = 0)
	 * 
	 */
	void add_double_mod(double& x, int position = 0) ;
	
	/** Returns the reference of a stored modifiable \c double .
	 * 
	 *  @param position [input] position of the modifiable \c double  
	 *                          in the list of stored modifiable 
	 *                      \c double  addresses (default value = 0)
	 * 
	 *  @return Reference to the modifiable \c double  the address of 
	 *                  which is stored at 
	 *		    the location  \c position  in the list
	 */
	double& get_double_mod(int position = 0) const; 
	

	///Returns the number of \c Tbl 's addresses in the list.
	int get_n_tbl() const ; 
    
	/** Adds the address of a new \c Tbl  to the list.
	 * 
	 *  @param ti [input] \c Tbl  the address of which is to be stored
	 *  @param position [input] position of the \c Tbl  in the list
	 *			    of stored \c Tbl  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_tbl(const Tbl& ti, int position = 0) ;
	
	/** Returns the reference of a \c Tbl  stored in the list.
	 * 
	 *  @param position [input] position of the \c Tbl  in the list
	 *			    of stored \c Tbl  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Tbl  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Tbl& get_tbl(int position = 0) const; 
	

	///Returns the number of modifiable \c Tbl 's addresses in the list.
	int get_n_tbl_mod() const ; 
    
	/** Adds the address of a new modifiable \c Tbl  to the list.
	 * 
	 *  @param ti [input] modifiable \c Tbl  the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the \c Tbl  in the list
	 *			    of stored modifiable \c Tbl  addresses 
	 *                          (default value = 0)
	 */
	void add_tbl_mod(Tbl& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable \c Tbl  stored in the list.
	 * 
	 *  @param position [input] position of the \c Tbl  in the list
	 *			    of stored modifiable \c Tbl  addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable \c Tbl  the address of 
	 *           which is stored at  the location  \c position  in the 
	 *           list
	 */
	 Tbl& get_tbl_mod(int position = 0) const; 
	

	///Returns the number of \c Itbl 's addresses in the list.
	int get_n_itbl() const ; 
    
	/** Adds the address of a new \c Itbl  to the list.
	 * 
	 *  @param ti [input] \c Itbl  the address of which is to be stored
	 *  @param position [input] position of the \c Itbl  in the list
	 *			    of stored \c Itbl  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_itbl(const Itbl& ti, int position = 0) ;
	
	/** Returns the reference of a \c Itbl  stored in the list.
	 * 
	 *  @param position [input] position of the \c Itbl  in the list
	 *			    of stored \c Itbl  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Itbl  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Itbl& get_itbl(int position = 0) const; 
	

	///Returns the number of modifiable \c Itbl 's addresses in the list.
	int get_n_itbl_mod() const ; 
    
	/** Adds the address of a new modifiable \c Itbl  to the list.
	 * 
	 *  @param ti [input] modifiable \c Itbl  the address of which is 
	 *         to be stored
	 *  @param position [input] position of the \c Itbl  in the list
	 *			    of stored modifiable \c Itbl  addresses 
	 *                          (default value = 0)
	 * 
	 */
	void add_itbl_mod(Itbl& ti, int position = 0) ;
	
	/** Returns the reference of a stored modifiable \c Itbl .
	 * 
	 *  @param position [input] position of the \c Itbl  in the list
	 *			    of stored modifiable \c Itbl  addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable \c Itbl  the address of 
	 *            which is stored at the location  \c position  
	 *            in the list
	 */
	Itbl& get_itbl_mod(int position = 0) const; 
	

	///Returns the number of \c Cmp 's addresses in the list.
	int get_n_cmp() const ; 
    
	/** Adds the address of a new \c Cmp  to the list.
	 * 
	 *  @param ti [input] \c Cmp  the address of which is to be stored
	 *  @param position [input] position of the \c Cmp  in the list
	 *			    of stored \c Cmp  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_cmp(const Cmp& ti, int position = 0) ;
	
	/** Returns the reference of a \c Cmp  stored in the list.
	 * 
	 *  @param position [input] position of the \c Cmp  in the list
	 *			    of stored \c Cmp  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Cmp  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Cmp& get_cmp(int position = 0) const; 
	

	///Returns the number of modifiable \c Cmp 's addresses in the list.
	int get_n_cmp_mod() const ; 
    
	/** Adds the address of a new modifiable \c Cmp  to the list.
	 * 
	 *  @param ti [input] modifiable \c Cmp  the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the \c Cmp  in the list
	 *			    of stored modifiable \c Cmp  addresses 
	 *                          (default value = 0)
	 */
	void add_cmp_mod(Cmp& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable \c Cmp  stored in the list.
	 * 
	 *  @param position [input] position of the \c Cmp  in the list
	 *			    of stored modifiable \c Cmp  addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable \c Cmp  the address of 
	 *           which is stored at  the location  \c position  in the 
	 *           list
	 */
	 Cmp& get_cmp_mod(int position = 0) const; 
	

	///Returns the number of \c Tenseur 's addresses in the list.
	int get_n_tenseur() const ; 
    
	/** Adds the address of a new \c Tenseur  to the list.
	 * 
	 *  @param ti [input] \c Tenseur  the address of which is to be stored
	 *  @param position [input] position of the \c Tenseur  in the list
	 *			    of stored \c Tenseur  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_tenseur(const Tenseur& ti, int position = 0) ;
	
	/** Returns the reference of a \c Tenseur  stored in the list.
	 * 
	 *  @param position [input] position of the \c Tenseur  in the list
	 *			    of stored \c Tenseur  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Tenseur  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Tenseur& get_tenseur(int position = 0) const; 
	

	///Returns the number of modifiable \c Tenseur 's addresses in the list.
	int get_n_tenseur_mod() const ; 
    
	/** Adds the address of a new modifiable \c Tenseur  to the list.
	 * 
	 *  @param ti [input] modifiable \c Tenseur  the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the \c Tenseur  in the list
	 *			    of stored modifiable \c Tenseur  addresses 
	 *                          (default value = 0)
	 */
	void add_tenseur_mod(Tenseur& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable \c Tenseur  stored in the list.
	 * 
	 *  @param position [input] position of the \c Tenseur  in the list
	 *			    of stored modifiable \c Tenseur  addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable \c Tenseur  the address of 
	 *           which is stored at  the location  \c position  in the 
	 *           list
	 */
	 Tenseur& get_tenseur_mod(int position = 0) const; 
	
	///Returns the number of \c Qtenseur 's addresses in the list.
	int get_n_qtenseur() const ; 
    
	/** Adds the address of a new \c Qtenseur  to the list.
	 * 
	 *  @param ti [input] \c Qtenseur  the address of which is to be stored
	 *  @param position [input] position of the \c Qtenseur  in the list
	 *			    of stored \c Qtenseur  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_qtenseur(const Qtenseur& ti, int position = 0) ;
	
	/** Returns the reference of a \c Qtenseur  stored in the list.
	 * 
	 *  @param position [input] position of the \c Qtenseur  in the list
	 *			    of stored \c Qtenseur  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Qtenseur  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Qtenseur& get_qtenseur(int position = 0) const; 
	

	///Returns the number of modifiable \c Qtenseur 's addresses in the list.
	int get_n_qtenseur_mod() const ; 
    
	/** Adds the address of a new modifiable \c Qtenseur  to the list.
	 * 
	 *  @param ti [input] modifiable \c Qtenseur  the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the \c Qtenseur  in the list
	 *			    of stored modifiable \c Qtenseur  addresses 
	 *                          (default value = 0)
	 */
	void add_qtenseur_mod(Qtenseur& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable \c Qtenseur  stored in the list.
	 * 
	 *  @param position [input] position of the \c Qtenseur  in the list
	 *			    of stored modifiable \c Qtenseur  addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable \c Qtenseur  the address of 
	 *           which is stored at  the location  \c position  in the 
	 *           list
	 */
	 Qtenseur& get_qtenseur_mod(int position = 0) const; 
	
	///Returns the number of \c Map 's addresses in the list.
	int get_n_map() const ; 
    
	/** Adds the address of a new \c Map  to the list.
	 * 
	 *  @param mi [input] \c Map  the address of which is to be stored
	 *  @param position [input] position of the \c Map  in the list
	 *			    of stored \c Map  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_map(const Map& mi, int position = 0) ;
	
	/** Returns the reference of a \c Map  stored in the list.
	 * 
	 *  @param position [input] position of the \c Map  in the list
	 *			    of stored \c Map  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Map  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Map& get_map(int position = 0) const; 
	
	///Returns the number of \c Mtbl_cf 's addresses in the list.
	int get_n_mtbl_cf() const ; 
    
	/** Adds the address of a new \c Mtbl_cf  to the list.
	 * 
	 *  @param mi [input] \c Mtbl_cf  the address of which is to be stored
	 *  @param position [input] position of the \c Mtbl_cf  in the list
	 *			    of stored \c Mtbl_cf  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_mtbl_cf(const Mtbl_cf& mi, int position = 0) ;
	
	/** Returns the reference of a \c Mtbl_cf  stored in the list.
	 * 
	 *  @param position [input] position of the \c Mtbl_cf  in the list
	 *			    of stored \c Mtbl_cf  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Mtbl_cf  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Mtbl_cf& get_mtbl_cf(int position = 0) const; 
	
	///Returns the number of \c Scalar 's addresses in the list.
	int get_n_scalar() const ; 
    
	/** Adds the address of a new \c Scalar  to the list.
	 * 
	 *  @param ti [input] \c Scalar  the address of which is to be stored
	 *  @param position [input] position of the \c Scalar  in the list
	 *			    of stored \c Scalar  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_scalar(const Scalar& ti, int position = 0) ;
	
	/** Returns the reference of a \c Scalar  stored in the list.
	 * 
	 *  @param position [input] position of the \c Scalar  in the list
	 *			    of stored \c Scalar  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Scalar  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Scalar& get_scalar(int position = 0) const; 
	

	///Returns the number of modifiable \c Scalar 's addresses in the list.
	int get_n_scalar_mod() const ; 
    
	/** Adds the address of a new modifiable \c Scalar  to the list.
	 * 
	 *  @param ti [input] modifiable \c Scalar  the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the \c Scalar  in the list
	 *			    of stored modifiable \c Scalar  addresses 
	 *                          (default value = 0)
	 */
	void add_scalar_mod(Scalar& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable \c Scalar  stored in the list.
	 * 
	 *  @param position [input] position of the \c Scalar  in the list
	 *			    of stored modifiable \c Scalar  addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable \c Scalar  the address of 
	 *           which is stored at  the location  \c position  in the 
	 *           list
	 */
	 Scalar& get_scalar_mod(int position = 0) const; 


	///Returns the number of \c Tensor 's addresses in the list.
	int get_n_tensor() const ; 
    
	/** Adds the address of a new \c Tensor  to the list.
	 * 
	 *  @param ti [input] \c Tensor  the address of which is to be stored
	 *  @param position [input] position of the \c Tensor  in the list
	 *			    of stored \c Tensor  addresses (default
	 *			    value = 0)
	 * 
	 */
	void add_tensor(const Tensor& ti, int position = 0) ;
	
	/** Returns the reference of a \c Tensor  stored in the list.
	 * 
	 *  @param position [input] position of the \c Tensor  in the list
	 *			    of stored \c Tensor  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Tensor  the address of which is stored at 
	 *		    the location  \c position  in the list
	 */
	const Tensor& get_tensor(int position = 0) const; 
	

	///Returns the number of modifiable \c Tensor 's addresses in the list.
	int get_n_tensor_mod() const ; 
    
	/** Adds the address of a new modifiable \c Tensor  to the list.
	 * 
	 *  @param ti [input] modifiable \c Tensor  the address of which is 
	 *                    to be stored
	 *  @param position [input] position of the \c Tensor  in the list
	 *			    of stored modifiable \c Tensor  addresses 
	 *                          (default value = 0)
	 */
	void add_tensor_mod(Tensor& ti, int position = 0) ;
	
	/** Returns the reference of a modifiable \c Tensor  stored in the list.
	 * 
	 *  @param position [input] position of the \c Tensor  in the list
	 *			    of stored modifiable \c Tensor  addresses 
	 *                          (default value = 0)
	 *  @return Reference to the modifiable \c Tensor  the address of 
	 *           which is stored at  the location  \c position  in the 
	 *           list
	 */
	 Tensor& get_tensor_mod(int position = 0) const; 

	///Returns the number of \c Etoile 's addresses in the list.
	int get_n_etoile() const ;

	/** Adds the address of a new \c Etoile  to the list.
	 *
	 *  @param mi [input] \c Etoile  the address of which is to be stored
	 *  @param position [input] position of the \c Etoile  in the list
	 *			    of stored \c Etoile  addresses (default
	 *			    value = 0)
	 *
	 */
	void add_etoile(const Etoile& eti, int position = 0) ;
	
	/** Returns the reference of a \c Etoile  stored in the list.
	 *
	 *  @param position [input] position of the \c Etoile  in the list
	 *			    of stored \c Etoile  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Etoile  the address of which is stored at
	 *		    the location  \c position  in the list
	 */
	const Etoile& get_etoile(int position = 0) const;
	
	///Returns the number of \c Star 's addresses in the list.
	int get_n_star() const ;

	/** Adds the address of a new \c Star  to the list.
	 *
	 *  @param mi [input] \c Star  the address of which is to be stored
	 *  @param position [input] position of the \c Star  in the list
	 *			    of stored \c Etoile  addresses (default
	 *			    value = 0)
	 *
	 */
	void add_star(const Star& eti, int position = 0) ;
	
	/** Returns the reference of a \c Star  stored in the list.
	 *
	 *  @param position [input] position of the \c Star  in the list
	 *			    of stored \c Star  addresses (default
	 *			    value = 0)
	 *  @return Reference to the \c Star  the address of which is stored at
	 *		    the location  \c position  in the list
	 */
	const Star& get_star(int position = 0) const;
	
 };

#endif
