/*
 *  Methods of class Eos_GlendNH3
 *
 *  (see file eos_tabul.h for documentation).
 *
 */

// Headers Lorene
#include "headcpp.h"
#include "eos.h"

			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
namespace Lorene {
Eos_GlendNH3::Eos_GlendNH3(const char* path)
		: Eos_tabul("EOS GlendNH3",
		            "eos_glendnh3.d", path)
{}


// Constructor from binary file
// ----------------------------
Eos_GlendNH3::Eos_GlendNH3(FILE* fich) : Eos_tabul(fich) {}



// Constructor from a formatted file
// ---------------------------------
Eos_GlendNH3::Eos_GlendNH3(ifstream& fich) :
			Eos_tabul(fich, "eos_glendnh3.d") {}



			//--------------//
			//  Destructor  //
			//--------------//

Eos_GlendNH3::~Eos_GlendNH3(){

    // does nothing

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_GlendNH3::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_GlendNH3 !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_GlendNH3::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_GlendNH3::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_GlendNH3 : Glendenning, N.K. 1985 ApJ293, 470, case 3 "
    	<< endl ;
    	
    ost << "  composition :  n,p,H" << endl ;
    ost << "  model : Lagrangian field theory; MFT"
        << endl ;
    ost << "  BPS EOS below neutron drip point" << endl ;
    ost << "  Sly4 EOS up to the liquid core" << endl ;
    ost << "  Crust bottom at n = 0.1 fm^{-3}" << endl ;

    return ost ;

}
}
