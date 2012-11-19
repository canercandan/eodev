#include <thread>
#include <mutex>
#include <atomic>
#include <queue>

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

size_t popSize;
size_t dimSize;

void send( eoGenContinue<EOT>& cont, eoPop<EOT>& pop )
{
    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    if ( 1 == RANK )
	{

	    while ( cont.thisGeneration < cont.totalGenerations() )
		{
		    {
			eoUniformGenerator<typename EOT::AtomType> gen(0, 2);
			eoInitFixedLength<EOT> init(dimSize, gen);

			eoPop<EOT> pop;
			pop.append( popSize, init );

			Node::comm().send(2, pop);
		    }

		    {
			eoPop<EOT> pop;
			Node::comm().recv(2, pop);

			eo::log << pop;
		    }

		    std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}

	}
 }

void recv( eoGenContinue<EOT>& cont, eoPop<EOT>& pop )
{
    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    if ( 2 == RANK )
	{

	    while ( cont.thisGeneration < cont.totalGenerations() )
		{
		    eoPop<EOT> pop;
		    Node::comm().recv(1, pop);

		    eoEvalFuncPtr<EOT> eval(binary_value);
		    apply<EOT>(eval, pop);

		    Node::comm().send(1, pop);

		    std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}

	}
}

void mainThread( eoGenContinue<EOT>& cont, eoPop<EOT>& pop )
{
    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    while ( cont(pop) )
	{
	}
}

int main(int ac, char** av)
{
    Node::init( ac, av );
    eoParser parser(ac, av);

    // S
    unsigned seed = parser.createParam(unsigned(time(0)), "seed", "Random number seed", 'S').value(); // will be in default section General
    // P
    /*size_t */popSize = parser.createParam(size_t(20), "popSize", "Population size",'P', "Evolution engine").value();
    // D
    /*size_t */dimSize = parser.createParam(size_t(10), "dimSize", "Dimension size",'D', "Evolution engine").value();
    // G
    unsigned maxGen = parser.createParam(unsigned(50), "maxGen", "Maximum number of iterations",'G', "Stopping criterion").value(); // also know as the number of migrations for IM

    make_parallel( parser );
    make_verbose(parser);
    make_help(parser);

    rng.reseed(seed);

    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    if( ALL < 3 )
	{
	    throw std::runtime_error("Needs at least 3 processes to be launched!");
	}

    mpi_debug();

    eoGenContinue<EOT> cont( maxGen );

    eoUniformGenerator<typename EOT::AtomType> gen(0, 2);
    eoInitFixedLength<EOT> init(dimSize, gen);

    eoPop<EOT> pop;
    pop.append(popSize, init);

    eoEvalFuncPtr<EOT> eval(binary_value);
    apply<EOT>(eval, pop);

    std::thread t1(send, std::ref(cont), std::ref(pop));
    std::thread t2(recv, std::ref(cont), std::ref(pop));
    std::thread t3(mainThread, std::ref(cont), std::ref(pop));

    t1.join();
    t2.join();
    t3.join();

    return 0;
}
