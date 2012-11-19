#ifndef __LOCALSEARCH_H__
#define __LOCALSEARCH_H__

#include <eo>

#include <ga/make_ga.h>
#include <ga/eoBitOp.h>

namespace my
{

    template <class EOT>
    class LocalSearch: public eoSequentialOp<EOT>
    {
    public:
	LocalSearch(unsigned nbIter = 100) : _nbIter(nbIter) {}

	/// The class name.
	virtual std::string className() const { return "LocalSearch"; }

	void apply(eoPopulator<EOT>& pop)
	{
	    for ( unsigned i = 0; i < _nbIter; ++i )
		{
		    eoSequentialOp<EOT>::apply( pop );
		}
	}

    private:
	unsigned _nbIter;
    };

}

#endif // !__LOCALSEARCH_H__
