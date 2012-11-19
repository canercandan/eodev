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

void emigrate( eoGenContinue<EOT>& cont, eoPop<EOT>& pop )
{
    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    if ( RANK != 1 ) { return; }

    while ( cont.thisGeneration < cont.totalGenerations() )
	{
	    // eo::log << eo::progress << "pack ";
	    // eo::log.flush();

	    // if (!pop.empty())
	    // 	{
	    em.push( std::make_pair(2, pop) );
	    // }
	}
}

void immigrate( eoGenContinue<EOT>& cont, eoPop<EOT>& pop )
{
    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    if ( RANK != 2 ) { return; }

    while ( cont.thisGeneration < cont.totalGenerations() )
	{
	    // eo::log << eo::progress << "unpack ";
	    // eo::log.flush();

	    while (!imm.empty())
	    	{
	    	    eoPop<EOT> newpop = imm.front().second;
	    	    imm.pop();
	    	}
	}
}

void sendRecv( eoGenContinue<EOT>& cont, eoPop<EOT>& pop )
{
    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    while ( cont(pop) )
	{
	    if ( 1 == RANK )
		{
		    if ( !em.empty() )
			{
			    size_t dst = em.front().first;
			    eoPop<EOT> pop = em.front().second;
			    Node::comm().send(dst, true);
			    Node::comm().send(dst, pop);
			    em.pop();
			}
		    else
			{
			    Node::comm().send(2, false);
			}
		}
	    else if ( 2 == RANK )
		{
		    mpi::status stat = Node::comm().probe(1, 0);
		    // eo::log << "stat: " << stat.source() << " ";
		    // eo::log.flush();

		    bool isValid = false;
		    Node::comm().recv(stat.source(), isValid);

		    // eo::log << "v(" << isValid << ") ";
		    // eo::log.flush();

		    if (isValid)
			{
			    eoPop<EOT> pop;
			    Node::comm().recv(stat.source(), pop);
			    imm.push( std::make_pair(stat.source(), pop) );
			}
		}

	    if ( RANK == 1 || RANK == 2 )
		{
		    // eo::log << "[" << RANK << ":" << pop.size() << "] ";
		    // eo::log.flush();
		}
	}
}

int main(int ac, char** av)
{
    Node::init( ac, av );
    eoParser parser(ac, av);

    // S
    unsigned seed = parser.createParam(unsigned(time(0)), "seed", "Random number seed", 'S').value(); // will be in default section General
    // s
    std::string schema = parser.createParam(std::string(""), "schema","an xml file mapping process",'s', "Persistence").value();

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

    std::thread t1(emigrate, std::ref(cont), std::ref(pop));
    std::thread t2(immigrate, std::ref(cont), std::ref(pop));
    std::thread t3(sendRecv, std::ref(cont), std::ref(pop));

    t1.join();
    t2.join();
    t3.join();

    return 0;
}
