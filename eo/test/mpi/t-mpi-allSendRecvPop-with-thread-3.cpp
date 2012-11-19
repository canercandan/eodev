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

std::queue< std::pair< size_t, eoPop<EOT> > > imm;
std::queue< std::pair< size_t, eoPop<EOT> > > em;

std::mutex imm_mutex;
std::mutex em_mutex;
std::mutex pop_mutex;

std::mutex mutex;

void emigrate( eoPop<EOT>& pop )
{
    const size_t ALL = Node::comm().size();
    // const size_t RANK = Node::comm().rank();

    eo::log << eo::quiet << "e";
    eo::log.flush();

    for ( size_t i = 0; i < ALL; ++i )
	{
	    // if ( i == RANK ) { continue; }

	    em.push( std::make_pair(i, pop) );

	    // ++count;
	}
}

void immigrate( eoPop<EOT>& pop )
{
    // const size_t ALL = Node::comm().size();
    // const size_t RANK = Node::comm().rank();

    eo::log << eo::quiet << "i";
    eo::log.flush();

    while (!imm.empty())
	{
	    eoPop<EOT> newpop = imm.front().second;
	    // pop.resize( pop.size() + newpop.size() );
	    // std::copy( newpop.begin(), newpop.end(), pop.begin() + newpop.size() );
	    imm.pop();
	}
}

void mainThread( eoGenContinue<EOT>& cont, eoPop<EOT>& pop )
{
    // const size_t ALL = Node::comm().size();
    const size_t RANK = Node::comm().rank();

    while ( cont(pop) )
	{
	    emigrate( pop );
	    immigrate( pop );

	    eo::log << eo::debug << "[" << RANK << ":" << pop.size() << "] ";
	    eo::log.flush();
	}
}

void sendThread( eoGenContinue<EOT>& cont, std::atomic< size_t >& count )
{
    // const int ALL = Node::comm().size();
    // const int RANK = Node::comm().rank();

    while ( cont.thisGeneration < cont.totalGenerations() )
	{

	    while ( !em.empty() )
		{
		    eo::log << eo::quiet << "s" << count;
		    eo::log.flush();

		    size_t dst = em.front().first;
		    eoPop<EOT> pop = em.front().second;
		    // Node::comm().send(dst, true);
		    Node::comm().send(dst, pop);
		    em.pop();

		    ++count;

		    // std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}

	}
}

void recvThread( eoGenContinue<EOT>& cont, std::atomic< size_t >& count )
{
    // const size_t ALL = Node::comm().size();
    // const size_t RANK = Node::comm().rank();

    while ( cont.thisGeneration < cont.totalGenerations() )
	{

	    if ( count > 0 )
		{
		    eo::log << eo::quiet << "r" << count;
		    eo::log.flush();

		    mpi::status stat = Node::comm().probe(mpi::any_source, 0);

		    eoPop<EOT> pop;
		    Node::comm().recv(stat.source(), pop);
		    imm.push( std::make_pair(stat.source(), pop) );

		    --count;

		    // std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}

	}
}

int main(int ac, char** av)
{
    Node::init( ac, av, MPI_THREAD_MULTIPLE );
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

    // const int ALL = Node::comm().size();
    // const int RANK = Node::comm().rank();

    mpi_debug();

    eoGenContinue<EOT> cont( maxGen );

    eoUniformGenerator<typename EOT::AtomType> gen(0, 2);
    eoInitFixedLength<EOT> init(dimSize, gen);

    eoPop<EOT> pop;
    pop.append(popSize, init);

    eoEvalFuncPtr<EOT> eval(binary_value);
    apply<EOT>(eval, pop);

    std::atomic< size_t > count( 0 );

    std::thread t1(mainThread, std::ref(cont), std::ref(pop));
    std::thread t2(sendThread, std::ref(cont), std::ref(count));
    std::thread t3(recvThread, std::ref(cont), std::ref(count));

    t1.join();
    t2.join();
    t3.join();

    return 0;
}
