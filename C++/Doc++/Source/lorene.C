/**@name LORENE --- Reference manual
 *
 * \URL[Lorene]{http://www.lorene.obspm.fr/}
 * (Langage Objet pour la RElativit\'e Num\'eriquE)
 * is a C++ based language for numerical relativity.
 *
 * Lorene home page is 
 * \URL[http://www.lorene.obspm.fr/]{http://www.lorene.obspm.fr/}
 *
 */

//@{
    /**@name Utilities
     *
     * The classes {\tt Tbl} and {\tt Itbl} implement the 1-D, 2-D and 3-D
     * array representation in Lorene, whereas the class {\tt Matrice}
     * is devoted to matrix representation and {\tt Param} to generic
     * parameter storage.  
     * All these classes are independent of the numerical method 
     * (spectral method).
     *
     */
    //@{
    	//@Include:dim_tbl.h_r
	//@Include:itbl.h_r
    	//@Include:tbl.h_r
    	//@Include:matrice.h_r
    	//@Include:param.h_r
	//@Include:utilitaires.h_r
    //@}

    /**@name Spectral representation
     *
     * These classes correspond to the implementation of spectral
     * methods in Lorene. They describe the collocation points
     * (classes {\tt Grille3d} and {\tt Mg3d}), the values
     * of a function at these points (class {\tt Mtbl}), and the coefficients
     * of the spectral expansions of a function (class {\tt Mtbl\_cf}). 
     *
     */
    //@{
    	//@Include:grilles.h_r
    	//@Include:mtbl.h_r
    	//@Include:base_val.h_r
    	//@Include:mtbl_cf.h_r
    	//@Include:valeur.h_r
    //@}

   /**@name Mapping grid -> physical space (spherical coordinates)
     *
     * These classes implement the mapping between the
     * grid coordinates $(\xi, \theta', \phi')$ (described by the {\tt Mg3d}
     * class) and the physical coordinates
     * $(r, \theta, \phi)$ [cf. Bonazzola, Gourgoulhon \& Marck, {\it Phys. Rev. D}
     * {\bf 58}, 104020 (1998)].
     *
     * The class {\tt Map} and its derived classes determine the methods for
     * partial derivatives with respect to the physical coordinates, as well
     * as resolution of basic partial differential equations (e.g. Poisson
     * equations).
     *
     */
    //@{
	//@Include:map.h_r
	//@Include:coord.h_r
    //@}

    /**@name General elliptic solver (***under development***)
     *
     * These classes are needed for using the general elliptic solver, for
     * which the variables and the operators can be different from one domain
     * to the other.
     */
    //@{
	//@Include:change_var.h_r
        //@Include:ope_elementary.h_r
        //@Include:param_elliptic.h_r
    //@}  

   /**@name Tensorial fields
     *
     * These classes implement the tensorial calculus in Lorene. 
     * They are high level classes and therefore are independent of the
     * actual numerical method (spectral method). 
     *
     */
    //@{
	//@Include:base_vect.h_r
	//@Include:tensor.h_r
	//@Include:scalar.h_r
	//@Include:vector.h_r
	//@Include:sym_tensor.h_r
	//@Include:connection.h_r
	//@Include:metric.h_r
    //@}  
    
   /**@name Old tensorial fields (*** Deprecated ***)
     *
     * These classes have been used up to 2003 to treat scalar 
     * (class {\tt Cmp}) and tensorial fields (class {\tt Tenseur}). 
     * They are now deprecated and
     * have been replaced by the class {\tt Tensor} and its various derived
     * classes, among which {\tt Scalar}.
     */
    //@{
	//@Include:cmp.h_r
	//@Include:tenseur.h_r
    //@}  

    /**@name Time evolution (***under development***)
     *
     * The storage and manipulation (e.g. time derivation) of an
     * evolving quantity is performed through the template class
     * {\tt Evolution}. 
     */
    //@{
	//@Include:evolution.h_r
    //@}  

    /**@name Grid wedding
     *
     * These classes are used to make the interface with Godunov-type
     * methods, used by the Valencia group.
     *
     */
    //@{
    	//@Include:grille_val.h_r
	//@Include:tbl_val.h_r
    //@}

    /**@name Physical units
     *
     */
    //@{
	//@Include:unites.h_r
	//@Include:unites_mag.h_r
    //@}

    /**@name Equations of state
     *
     */
    //@{
	//@Include:eos.h_r
        //@Include:eos_tabul.h_r
        //@Include:eos_bifluid.h_r
    //@}

    /**@name Stars and black holes
     *
     */
    //@{
	//@Include:etoile.h_r
        //@Include:et_rot_bifluid.h_r
        //@Include:et_rot_diff.h_r
        //@Include:et_rot_mag.h_r
        //@Include:et_bin_ncp.h_r
        //@Include:star.h_r
	//@Include:binaire.h_r
	//@Include:binary.h_r
	//@Include:bhole.h_r
	//@Include:bin_ns_bh.h_r
        //@Include:bin_ns_ncp.h_r
    //@}
    
    /**@name Graphical outputs
     *
     * The 2-D graphical outputs of various Lorene objects 
     * are performed via the PGPLOT library:
     *
     * \URL[http://astro.caltech.edu/~tjp/pgplot/]{http://astro.caltech.edu/~tjp/pgplot/}
     *
     *
     * To open an X11 display area, use the command 'pgdisp'.
     *
     * 3-D visualization of various Lorene objects is performed by
     * means of OpenDX:
    *
    * \URL[http://www.opendx.org/]{http://www.opendx.org/}
    *
    * OpenDX is called from the methods {\tt Scalar::visu\_section}, 
    *  {\tt Scalar::visu\_box}, {\tt Vector::visu\_arrows}, etc...
    * The corresponding OpenDX scripts are provided in the 
    * directory {\tt Lorene/Visu/OpenDX}. They must be copied to the
    * working directory in order to use the above methods.  
    *
    */
    //@{
	//@Include:graphique.h_r
    //@}


//@}
