#include <thread>
#include <mutex>
#include <queue>

#include <mpi/eoMpi.h>

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

// class Communicable
// {
// public:
//     Communicable() {}
//     // Communicable() : t1(&Communicable::pack, this), t2(&Communicable::unpack, this) {}
//     virtual ~Communicable() {}

//     virtual void pack() = 0;
//     virtual void unpack() = 0;

//     void run()
//     {
// 	t1 = std::thread( &Communicable::pack, this );
// 	t2 = std::thread( &Communicable::unpack, this );
// 	// t1.detach();
// 	// t2.detach();
//     }

// protected:
//     std::thread t1;
//     std::thread t2;
// };

class Cooperative/* : public Communicable*/ {};

// class Topology
// {
// public:
//     virtual ~Topology()
//     {
// 	/* Nothing ! */
//     }

//     void add(Cooperative & __mig)
//     {
// 	mig.push_back (& __mig);
//     }

//     virtual void setNeighbors (Cooperative * __mig,
// 			       std::vector <Cooperative *> & __from,
// 			       std::vector <Cooperative *> & __to) = 0;

//     operator std::vector<Cooperative *>& ()
//     {
// 	return mig;
//     }

// protected:
//     std::vector<Cooperative *> mig;
// };

// class RingTopology : public Topology
// {
// public:
//     void setNeighbors (Cooperative * __mig,
//                        std::vector <Cooperative *> & __from,
//                        std::vector <Cooperative *> & __to)
//     {
// 	__from.clear() ;
// 	__to.clear() ;

// 	int len = mig.size() ;

// 	for (int i = 0; i < len; i ++)
// 	    {
// 		if (mig[i] == __mig)
// 		    {
// 			__from.push_back(mig [(i - 1 + len) % len]) ;
// 			__to.push_back(mig [(i + 1) % len]) ;
// 			break;
// 		    }
// 	    }
//     }
// };

// class CompleteTopology : public Topology
// {
// public:
//     void setNeighbors (Cooperative * __mig,
// 		       std::vector <Cooperative *> & __from,
// 		       std::vector <Cooperative *> & __to)
//     {
// 	__from.clear() ;
// 	__to.clear() ;

// 	for (unsigned i = 0; i < mig.size(); i ++)
// 	    {
// 		if (mig [i] != __mig)
// 		    {
// 			__from.push_back(mig [i]);
// 			__to.push_back(mig [i]);
// 		    }
// 	    }
//     }
// };

template < typename EOT >
class Island : public Cooperative, public eoF< void >
{
public:
    Island( /*Topology& __topology, */int __size, int __rank )
	: /*topology( __topology ),*/ size( __size ), rank( __rank )//, launched( false ), toStop( false )
	  //, pack(em, em_mutex), unpack(imm, imm_mutex, __size)
    {
	//__topology.add( *this );
    }

    virtual ~Island()
    {
	// t1.join();
	// t2.join();
    }

    // void run()
    // {
    // 	// t1 = std::thread( &Island<EOT>::pack, this );
    // 	// t2 = std::thread( &Island<EOT>::unpack, this );

    // 	// t1 = std::thread( pack );
    // 	// t2 = std::thread( unpack );

    // 	// t1.detach();
    // 	// t2.detach();
    // }

    void pack( eoPop<EOT>& pop )
    {
	std::cout << "pack ";
	std::cout.flush();
	// std::this_thread::sleep_for(std::chrono::milliseconds(10));

	// em_mutex.lock();
	for ( size_t i = 0; i < size; ++i )
	    {
		// coop_em.push( i );
		// if ( pop.size() - 1 )
		//     {
		// 	pop.resize( pop.size() - 1 );
		//     }
		// em.push( pop );

		// em_mutex.lock();
		coop_em.push( i );
		em.push( pop );
		// em_mutex.unlock();
	    }
	// em_mutex.unlock();
    }

    void unpack( eoPop<EOT>& pop )
    {
	std::cout << "unpack ";
	std::cout.flush();
	// std::this_thread::sleep_for(std::chrono::milliseconds(10));

	// imm_mutex.lock();
	while ( !imm.empty() )
	    {
		eoPop<EOT>& source = imm.front();
		// pop = source;
		imm.pop();
	    }
	// imm_mutex.unlock();

	eo::log << "[" << rank << ":" << pop.size() << "] ";
	eo::log.flush();
    }

    void operator()()
    {
	emigrate();
	immigrate();
    }

    void emigrate()
    {
	// em_mutex.lock();
	while ( !em.empty() )
	    {
		size_t dst = coop_em.front();
		eoPop<EOT> pop = em.front();

		if (!pop.empty())
		    {
			Node::comm().send( dst, true ); // it works
			Node::comm().send( dst, pop );
		    }
		else
		    {
			Node::comm().send( dst, false ); // doesnot works
		    }

		em.pop();
		coop_em.pop();
	    }
	// em_mutex.unlock();
    }

    void immigrate()
    {
	// imm_mutex.lock();
	// for ( size_t i = 0; i < size; ++i )
	//     {
	bool isValid = false;
	Node::comm().recv( mpi::any_source, isValid );
	std::cout << "v(" << isValid << ") ";
	std::cout.flush();
	if ( isValid )
	    {
		eoPop<EOT> pop;
		Node::comm().recv( mpi::any_source, pop );
		imm.push(pop);
	    }
	// }
	// imm_mutex.unlock();
    }

private:
    // Topology& topology;
    int size;
    int rank;

    std::queue< eoPop<EOT> > imm;
    // std::queue< std::pair< size_t, eoPop<EOT> > > em;
    std::queue< eoPop<EOT> > em;
    std::queue< size_t > coop_em;
    // std::queue< Cooperative* > coop_em;

    std::mutex imm_mutex;
    std::mutex em_mutex;
    std::mutex pop_mutex;

    std::thread t1;
    std::thread t2;

    // bool launched;
    // bool toStop;

    // std::mutex toStop_mutex;
};

// an eoQuadOp that does nothing
template <typename EOT>
class DummyQuadOp : public eoQuadOp<EOT>
{
public :
    std::string className() const {return "DummyQuadOp";}
    bool operator()(EOT&, EOT&) {return false;}
};

template <typename EOT>
class Algo : public eoUF< eoPop<EOT>&, void >
{
public:
    Algo( size_t maxGen, eoEasyEA<EOT>& algo, Island<EOT>& isl ) : _maxGen( maxGen ), _algo( algo ), _isl( isl ) {}

    void operator()( eoPop<EOT>& pop )
    {
	for ( size_t i = 0; i < _maxGen; ++i )
	    {
		_algo( pop );
		// _isl( pop );
		_isl.pack( pop );
		_isl.unpack( pop );
		// Node::comm().barrier();
	    }
    }

private:
    size_t _maxGen;
    eoEasyEA<EOT>& _algo;
    Island<EOT>& _isl;
};

typedef eoBit<double> EOT;

int main(int ac, char** av)
{
    Node::init( ac, av, MPI_THREAD_FUNNELED );
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

    // CompleteTopology topology;

    Island<EOT> isl( /*topology, */ALL, RANK );

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

    Algo<EOT> algom( maxGen, algo, isl );

    std::thread th = std::thread( algom, std::ref(pop) );

    isl();

    th.join();

    return 0;
}
