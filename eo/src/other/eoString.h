// -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

//-----------------------------------------------------------------------------
// eoString.h
// (c) GeNeura Team, 1998
/* 
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Contact: todos@geneura.ugr.es, http://geneura.ugr.es
 */
//-----------------------------------------------------------------------------

#ifndef _eoString_H
#define _eoString_H

// STL libraries
#include <string>		
#include <stdexcept>

using namespace std;

//-----------------------------------------------------------------------------
// eoString
//-----------------------------------------------------------------------------

/** Adaptor that turns an STL string into an EO */
template <class fitnessT >
class eoString: public EO<fitnessT>, public string 
{
public:

    typedef char Type;

	/// Canonical part of the objects: several ctors, copy ctor, dtor and assignment operator
	//@{
	/// ctor
	eoString( const string& _str ="" )
		: string( _str ) {};

	/** @name Methods from eoObject
	readFrom and printOn are directly inherited from eo1d
	*/
	//@{
	/** Inherited from eoObject 
		  @see eoObject
	*/
	virtual string className() const {return "eoString";};
    //@}
	

};

#endif

