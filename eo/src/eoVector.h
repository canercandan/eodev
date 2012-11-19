/* (C) GeNeura Team, 2000 - EEAAX 1999 - Maarten Keijzer 2000

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

Contact:     eodev-main@lists.sourceforge.net
Old contact: todos@geneura.ugr.es, http://geneura.ugr.es
             Marc.Schoenauer@polytechnique.fr
             mak@dhi.dk
*/

#ifndef _eoVector_h
#define _eoVector_h

#include <vector>
#include <iterator>
#include <EO.h>
#include <utils/eoLogger.h>

#ifdef WITH_MPI
#include <serial/eoSerial.h>

// #ifdef WITH_BOOST
#include <boost/mpi.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/assume_abstract.hpp>
// #endif // !WITH_BOOST

#endif // !WITH_MPI

/**
 @defgroup Representations Representations

 Solution to a given optimization problem are using a given representation, and are called individuals.

 Some of the most classical representations are proposed, but you can create your own one, providing
 that it inherits from the EO base class. It will be used as a template parameter in almost all operators.

 @{
*/

/** Base class for fixed length chromosomes

It just derives from EO and std::vector and redirects the smaller than
operator to EO (fitness based comparison).

GeneType must have the following methods: void ctor (needed for the
std::vector<>), copy ctor,

*/
template <class FitT, class GeneType>
class eoVector : public EO<FitT>, public std::vector<GeneType>
#ifdef WITH_MPI
	       , public eoserial::Persistent
#endif // !WITH_MPI
{
public:

#if defined(WITH_MPI) and defined(WITH_BOOST)
    friend class boost::serialization::access;
#endif

    using EO<FitT>::invalidate;
    using std::vector<GeneType>::operator[];
    using std::vector<GeneType>::begin;
    using std::vector<GeneType>::end;
    using std::vector<GeneType>::resize;
    using std::vector<GeneType>::size;

    typedef GeneType                AtomType;
    typedef std::vector<GeneType>   ContainerType;

    /** default constructor

    @param _size Length of vector (default is 0)
    @param _value Initial value of all elements (default is default value of type GeneType)
    */
    eoVector(unsigned _size = 0, GeneType _value = GeneType())
        : EO<FitT>(), std::vector<GeneType>(_size, _value)
        {}

    /// copy ctor abstracting from the FitT
    template <class OtherFitnessType>
    eoVector(const eoVector<OtherFitnessType, GeneType>& _vec) : std::vector<GeneType>(_vec)
        {}

    // we can't have a Ctor from a std::vector, it would create ambiguity
    //  with the copy Ctor
    void value(const std::vector<GeneType>& _v)
        {
            if (_v.size() != size())       // safety check
            {
                if (size())                // NOT an initial empty std::vector
                    eo::log << eo::warnings << "Warning: Changing size in eoVector assignation" << std::endl;
                resize(_v.size());
            }

            std::copy(_v.begin(), _v.end(), begin());
            invalidate();
        }

    /// to avoid conflicts between EO::operator< and std::vector<GeneType>::operator<
    bool operator<(const eoVector<FitT, GeneType>& _eo) const
        {
            return EO<FitT>::operator<(_eo);
        }

    /// printing...
    virtual void printOn(std::ostream& os) const
        {
            EO<FitT>::printOn(os);
            os << ' ';

            os << size() << ' ';

            std::copy(begin(), end(), std::ostream_iterator<AtomType>(os, " "));
        }

    /// reading...
    virtual void readFrom(std::istream& is)
        {
            EO<FitT>::readFrom(is);

            unsigned sz;
            is >> sz;

            resize(sz);
            unsigned i;

            for (i = 0; i < sz; ++i)
            {
                AtomType atom;
                is >> atom;
                operator[](i) = atom;
            }
        }

#ifdef WITH_MPI
    void unpack( const eoserial::Object* obj )
	{
	    this->clear();
	    eoserial::unpackArray< std::vector<GeneType>, eoserial::Array::UnpackAlgorithm >( *obj, "vector", *this );

	    bool invalidFitness;
            eoserial::unpack( *obj, "invalid_fitness", invalidFitness );
	    if( invalidFitness )
		{
		    this->invalidate();
		}
	    else
		{
		    FitT f;
		    eoserial::unpack( *obj, "fitness", f );
		    this->fitness( f );
		}

	    // DIM
	    lastIslands.clear();
	    eoserial::unpackArray< std::vector< size_t >, eoserial::Array::UnpackAlgorithm >( *obj, "lastIslands", lastIslands );
	    lastIslandsCount.clear();
	    eoserial::unpackArray< std::vector< size_t >, eoserial::Array::UnpackAlgorithm >( *obj, "lastIslandsCount", lastIslandsCount );
	    lastFitnesses.clear();
	    eoserial::unpackArray< std::vector< FitT >, eoserial::Array::UnpackAlgorithm >( *obj, "lastFitnesses", lastFitnesses );
	    // !DIM
	}

    eoserial::Object* pack( eoserial::Object* json = NULL ) const
	{
	    eoserial::Object* obj = json;
	    if (NULL == obj)
		{
		    obj = new eoserial::Object;
		}
	    obj->add( "vector", eoserial::makeArray< std::vector<GeneType>, eoserial::MakeAlgorithm >( *this ) );

	    bool invalidFitness = this->invalid();
	    obj->add( "invalid_fitness", eoserial::make( invalidFitness ) );
	    if( !invalidFitness )
		{
		    obj->add( "fitness", eoserial::make( this->fitness() ) );
		}

	    // DIM
	    obj->add( "lastIslands", eoserial::makeArray< std::vector< size_t >, eoserial::MakeAlgorithm >( lastIslands ) );
	    obj->add( "lastIslandsCount", eoserial::makeArray< std::vector< size_t >, eoserial::MakeAlgorithm >( lastIslandsCount ) );
	    obj->add( "lastFitnesses", eoserial::makeArray< std::vector< FitT >, eoserial::MakeAlgorithm >( lastFitnesses ) );
	    // !DIM

	    return obj;
	}

/* DIM */
public:
    void addIsland( size_t isl )
    {
	if ( getLastIsland() == int(isl) )
	{
	    ++lastIslandsCount.back();
	}
	else
	{
	    lastIslands.push_back( isl );
	    lastIslandsCount.push_back( 0 );
	}
    }

    void addFitness()
    {
	lastFitnesses.push_back( this->invalid() ? -1 : this->fitness() );
    }

    inline int getLastIsland() const { return lastIslands.empty() ? -1 : lastIslands.back(); }
    inline int getLastIslandCount() const { return lastIslandsCount.empty() ? -1 : lastIslandsCount.back(); }
    inline FitT getLastFitness() const { return lastFitnesses.empty() ? -1 : lastFitnesses.back(); }

    inline const std::vector<size_t>& getLastIslands() const { return lastIslands; }
    inline const std::vector<size_t>& getLastIslandsCount() const { return lastIslandsCount; }
    inline const std::vector<FitT>& getLastFitnesses() const { return lastFitnesses; }

    void printLastIslands() const
    {
	if (!lastIslands.size()) { return; }

	if ( lastIslandsCount.back() > 0 )
	{
	    std::cout << "(" << lastIslands.back() << "," << lastIslandsCount.back() << ")";
	}
	else
	{
	    std::cout << lastIslands.back();
	}
	std::cout.flush();

	for ( size_t i = lastIslands.size()-1; i > 0; --i )
	    {
		if ( lastIslandsCount[i-1] > 0 )
		{
		    std::cout << " < (" << lastIslands[i-1] << "," << lastIslandsCount[i-1] << ")";
		}
		else
		{
		    std::cout << " < " << lastIslands[i-1];
		}
		std::cout.flush();
	    }
	std::cout << std::endl;
	std::cout.flush();
    }

    size_t getLastIslandsSize() const
    {
	return lastIslands.size();
    }

private:
    std::vector< size_t > lastIslands;
    std::vector< size_t > lastIslandsCount;
    std::vector< FitT > lastFitnesses;
/* !DIM */

#endif // !WITH_MPI

public:
#if defined(WITH_MPI)// and defined(WITH_BOOST)
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/)
	{
	    ar & boost::serialization::base_object< EO<FitT> >(*this);
	    ar & boost::serialization::base_object< std::vector<GeneType> >(*this);
	    ar & lastIslands & lastIslandsCount & lastFitnesses;
	}
#endif

};
/** @example t-eoVector.cpp
 */

/** Less than

This is impemented to avoid conflicts between EO::operator< and
std::vector<GeneType>::operator<
*/
template <class FitT, class GeneType>
bool operator<(const eoVector<FitT, GeneType>& _eo1, const eoVector<FitT, GeneType>& _eo2)
{
    return _eo1.operator<(_eo2);
}


/** Greater than

This is impemented to avoid conflicts between EO::operator> and
std::vector<GeneType>::operator>
*/
template <class FitT, class GeneType>
bool operator>(const eoVector<FitT, GeneType>& _eo1, const eoVector<FitT, GeneType>& _eo2)
{
    return _eo1.operator>(_eo2);
}

#if defined(WITH_MPI)// and defined(WITH_BOOST)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(eoVector)
#endif

#endif

/** @} */

// Local Variables:
// coding: iso-8859-1
// mode: C++
// c-file-offsets: ((c . 0))
// c-file-style: "Stroustrup"
// fill-column: 80
// End:
