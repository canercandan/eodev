#include <mpi/eoMpi.h>
#include <mpi/eoParallelApply.h>
#include <mpi/eoTerminateJob.h>

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

class Cooperative
{
public:
    virtual void pack() = 0;
    virtual void unpack() = 0;
};

class Topology
{
public:
    virtual ~Topology()
    {
	/* Nothing ! */
    }

    void add(Cooperative & __mig)
    {
	mig.push_back (& __mig);
    }

    virtual void setNeighbors (Cooperative * __mig,
			       std::vector <Cooperative *> & __from,
			       std::vector <Cooperative *> & __to) = 0;

    operator std::vector<Cooperative *>& ()
    {
	return mig;
    }

protected:
    std::vector<Cooperative *> mig;
};

class RingTopology : public Topology
{
public:
    void setNeighbors (Cooperative * __mig,
                       std::vector <Cooperative *> & __from,
                       std::vector <Cooperative *> & __to)
    {
	__from.clear() ;
	__to.clear() ;

	int len = mig.size() ;

	for (int i = 0; i < len; i ++)
	    {
		if (mig[i] == __mig)
		    {
			__from.push_back(mig [(i - 1 + len) % len]) ;
			__to.push_back(mig [(i + 1) % len]) ;
			break;
		    }
	    }
    }
};

class CompleteTopology : public Topology
{
public:
    void setNeighbors (Cooperative * __mig,
		       std::vector <Cooperative *> & __from,
		       std::vector <Cooperative *> & __to)
    {
	__from.clear() ;
	__to.clear() ;

	for (unsigned i = 0; i < mig.size(); i ++)
	    {
		if (mig [i] != __mig)
		    {
			__from.push_back(mig [i]);
			__to.push_back(mig [i]);
		    }
	    }
    }
};

template < typename EOT >
class Island : public Cooperative, public eoUF< eoPop<EOT>&, void >
{
public:
    Island( Topology& __topology, int __size, int __rank )
	: topology( __topology ), size( __size ), rank( __rank )
    {
	__topology.add( *this );
    }

    void pack()
    {
	// empty
    }


    void unpack()
    {
	// empty

    }

    void operator()( eoPop<EOT>& pop )
    {
	for ( size_t i = 0; i < size; ++i )
	    {
		// if ( i == rank ) { continue; }
		if ( pop.size() - 1 )
		    {
			pop.resize( pop.size() - 1 );
			Node::comm().send(i, 1); // it works
			Node::comm().send(i, pop);
		    }
		else
		    {
			Node::comm().send(i, 0); // doesnot works
		    }
	    }

	for ( size_t i = 0; i < size; ++i )
	    {
		// if ( i == rank ) { continue; }
		int stat;
		Node::comm().recv(i, stat);
		if ( stat )
		    {
			Node::comm().recv(i, pop);
		    }
	    }

	eo::log << "[" << rank << ":" << pop.size() << "] ";
	eo::log.flush();
    }

private:
    Topology& topology;
    int size;
    int rank;
};

// an eoQuadOp that does nothing
template <typename EOT>
class DummyQuadOp : public eoQuadOp<EOT>
{
public :
    std::string className() const {return "DummyQuadOp";}
    bool operator()(EOT&, EOT&) {return false;}
};

typedef eoBit<double> EOT;

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

    // C
    double pCross = parser.createParam(double(0.6), "pCross", "Probability of Crossover", 'C', "Genetic Operators").value();
    // M
    double pMut = parser.createParam(double(0.1), "pMut", "Probability of Mutation", 'M', "Genetic Operators").value();

    make_parallel( parser );
    make_verbose(parser);
    make_help(parser);

    rng.reseed(seed);

    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    mpi_debug();

    CompleteTopology topology;

    Island<EOT> isl( topology, ALL, RANK );

    eoUniformGenerator<typename EOT::AtomType> gen(0, 2);
    eoInitFixedLength<EOT> init(dimSize, gen);

    eoPop<EOT> pop;
    pop.append( popSize, init );

    eoDetSelect<EOT> select;
    eoEvalFuncPtr<EOT> eval(binary_value);
    DummyQuadOp<EOT> xover;
    eoOneBitFlip<EOT> bitflip;
    eoSGATransform<EOT> transform(xover, pCross, bitflip, pMut);
    eoPlusReplacement<EOT> replace;

    eoPeriodicContinue<EOT> periodicCont( 1 );
    eoCheckPoint<EOT> checkpoint( periodicCont );

    eoEasyEA<EOT> algo( checkpoint, eval, select, transform, replace );

    apply(eval, pop);

    for ( size_t i = 0; i < maxGen; ++i )
	{
	    algo( pop );
	    isl( pop );
 	    Node::comm().barrier();
	}

    return 0;
}
