#ifndef __MAKE_LS_OP_H__
#define __MAKE_LS_OP_H__

#include <eo>

//#include <ga/make_ga.h>
#include <ga/eoBitOp.h>

#include "LocalSearch.h"

namespace my
{
    /* Create local search operators and there parameters */
    template <class EOT>
    eoGenOp<EOT>& make_ls_op(eoParser& parser, eoState& state, eoInit<EOT>& init)
    {
	eoValueParam<unsigned>& kflipParam = parser.createParam(unsigned(0), "kflip", "K-Flip, 0: disabled", 0, "Problem");
	eoValueParam<bool>& bitflipParam = parser.createParam(false, "bitflip", "Bit-Flip", 0, "Problem");
	eoValueParam<bool>& uniformParam = parser.createParam(true, "uniform", "Uniform", 0, "Problem");
	eoValueParam<unsigned>& uniformNbOpsParam = parser.createParam(unsigned(5), "uniform-nb-ops", "Uniform: number of operators (n + bitflip)", 0, "Problem");

	eoValueParam<unsigned>& nbMoveLSParam = parser.createParam(unsigned(100), "nbMoveLS", "Numero of moves by individual for Local Search", 0, "Problem");

	eoSequentialOp<EOT>* ptSeqOp = new LocalSearch<EOT>(nbMoveLSParam.value());
	state.storeFunctor(ptSeqOp);

	eoMonOp<EOT>* ptMon = NULL;

	if ( kflipParam.value() )
	    {
		std::cout << "kflip(" << kflipParam.value() << ")" << std::endl;
		ptMon = new eoDetBitFlip<EOT>(kflipParam.value());
		state.storeFunctor(ptMon);
		ptSeqOp->add(*ptMon, 1);
	    }
	else if ( bitflipParam.value() )
	    {
		std::cout << "bitflip" << std::endl;
		ptMon = new eoBitMutation<EOT>(1, true);
		state.storeFunctor(ptMon);
		ptSeqOp->add(*ptMon, 1);
	    }
	else if ( uniformParam.value() )
	    {
		std::cout << "uniform" << std::endl;
		eoOpContainer<EOT>* ptPropOp = new eoProportionalOp<EOT>;
		state.storeFunctor(ptPropOp);

		ptMon = new eoBitMutation<EOT>(1, true);
		state.storeFunctor(ptMon);
		double rate = 1/(double)uniformNbOpsParam.value();
		ptPropOp->add(*ptMon, rate);

		for (unsigned i = 1; i <= uniformNbOpsParam.value(); ++i)
		    {
			ptMon = new eoDetBitFlip<EOT>( i );
			state.storeFunctor(ptMon);
			ptPropOp->add(*ptMon, rate);
		    }

		ptSeqOp->add(*ptPropOp, 1);
	    }
	else
	    {
		// exception to throw
	    }

	return *ptSeqOp;
    }

}

#endif // !__MAKE_LS_OP_H__
