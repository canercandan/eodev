#include <eompi.h>
#include <eo>
#include <ga.h>

#include "t-mpi-common.h"

using namespace eo::mpi;

template <typename EOT>
double binary_value(const EOT& _chrom)
{
    double sum = 0.0;
    for (unsigned i=0; i<_chrom.size(); i++)
	sum += _chrom[i];
    return sum;
}

typedef eoBit<double> EOT;

int main(int ac, char** av)
{
    // Node::init( ac, av, MPI_THREAD_MULTIPLE );
    // Node::init( ac, av, MPI_THREAD_SERIALIZED );
    Node::init( ac, av );
    eoParser parser(ac, av);

    // S
    unsigned seed = parser.createParam(unsigned(time(0)), "seed", "Random number seed", 'S').value(); // will be in default section General
    // P
    size_t popSize = parser.createParam(size_t(20), "popSize", "Population size",'P', "Evolution engine").value();
    // D
    size_t dimSize = parser.createParam(size_t(10), "dimSize", "Dimension size",'D', "Evolution engine").value();
    // G
    unsigned maxGen = parser.createParam(unsigned(50), "maxGen", "Maximum number of iterations",'G', "Stopping criterion").value(); // also know as the number of migrations for IM

    make_parallel( parser );
    make_verbose(parser);
    make_help(parser);

    rng.reseed(seed);

    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    mpi_debug();

    if ( 1 == RANK )
	{
	    {
		eoUniformGenerator<typename EOT::AtomType> gen(0, 2);
		eoInitFixedLength<EOT> init(dimSize, gen);

		eoPop<EOT> pop;
		pop.append( popSize, init );

		Node::comm().isend(2, pop);
	    }

	    {
		eoPop<EOT> pop;
		Node::comm().irecv(2, pop);

		eo::log << pop;
	    }
	}
    else if ( 2 == RANK )
	{
	    eoPop<EOT> pop;
	    Node::comm().irecv(1, pop);

	    eoEvalFuncPtr<EOT> eval(binary_value);
	    apply<EOT>(eval, pop);

	    Node::comm().isend(1, pop);
	}

    return 0;
}
