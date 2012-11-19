#include <stdexcept>
#include <iostream>
#include <iterator>
#include <vector>
#include <queue>
#include <algorithm>
#include <assert.h>
#include <mpi.h>

#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <condition_variable>

#include <eompi.h>
#include <eo>
#include <ga.h>

#include "t-mpi-common.h"

using namespace std;
using namespace eo::mpi;
using namespace MPI;

template < typename Atom >
class SquareMatrix : public std::vector< Atom >, public eoObject, public eoPersistent
{
public:
    // Ctor : sets size
    SquareMatrix(unsigned _s = 0) : std::vector< Atom >(_s*_s), rSize(_s) {}

    /** simple accessor */
    const Atom& operator()(unsigned _i, unsigned _j) const
    {
	return this->operator[](_i*rSize + _j);
    }

    /** reference - to set values */
    Atom& operator()(unsigned _i, unsigned _j)
    {
	return this->operator[](_i*rSize + _j);
    }

    /** returns a vector value in the row _i */
    std::vector< Atom > operator()(unsigned _i)
    {
	typename std::vector< Atom >::iterator begin = this->begin();
	std::advance(begin, _i * rSize);
	typename std::vector< Atom >::iterator end = begin;
	std::advance(end, rSize);
	std::vector< Atom > vec(begin, end);
	// std::copy(vec.begin(), vec.end(), std::ostream_iterator<Atom>(std::cout, " "));
	// std::cout << std::endl;
	return vec;
    }

    /** just in case */
    virtual void printOn(std::ostream & _os) const
    {
	unsigned index=0;
	for (unsigned i=0; i<rSize; i++)
	    {
		for (unsigned j=0; j<rSize; j++)
		    {
			_os << this->operator[](index++) << " " ;
		    }
		_os << std::endl;
	    }
	_os << std::endl;
    }

    virtual void readFrom(std::istream& /*_is*/)
    {
    }

    virtual std::string className() const {return "SquareMatrix";};

    size_t size() const {return rSize;}

private:
    unsigned rSize;              // row size (== number of columns!)
};

class MigrationMatrix : public SquareMatrix< double >
{
public:
    MigrationMatrix(unsigned s = 0) : SquareMatrix< double >(s) {}

    virtual void printOn(std::ostream & os) const
    {
	os << "Migration probability (in %) among islands" << std::endl;
	os << "\t" << "1/l";

	for (size_t i = 0; i < this->size() - 1; ++i)
	    {
		os << "\t" << i * 2 + 1 << "f";
	    }

	os << std::endl;

	for (size_t i = 0; i < this->size(); ++i)
	    {
		if (i == 0)
		    os << "1/l";
		else
		    os << i*2-1 << "f";

		for (size_t j = 0; j < this->size(); ++j)
		    {
			os << "\t" << floor((*this)(i,j) * 100) / 1000;
		    }

		os << std::endl;
	    }

	os << std::endl;
	os << "sum";

	double sum;

	for (size_t i = 0; i < this->size(); ++i)
	    {
		sum = 0;
		for (size_t j = 0; j < this->size(); ++j)
		    {
			sum = sum + (*this)(j,i);
		    }
		os << "\t" << floor(sum * 100) / 1000;
	    }
	os << std::endl;
    }
};

class InitMatrix : public eoUF< SquareMatrix<double>&, void >
{
public:
    InitMatrix(int initG = 0, double same = 90) : _initG(initG), _same(same) {}

    void operator()(SquareMatrix<double>& matrix)
    {
	double sum;

	for (size_t i = 0; i < matrix.size(); ++i)
	    {
		sum = 0;

		for (size_t j = 0; j < matrix.size(); ++j)
		    {
			if (i == j)
			    matrix(i,j) = _same * 10;
			else
			    {
				if (_initG == 0)
				    matrix(i,j) = rand();
				else
				    matrix(i,j) = (1000 - _same * 10) / (matrix.size() - 1);
				sum = sum + matrix(i,j);
			    }
		    }

		for (size_t j = 0; j < matrix.size(); ++j)
		    {
			if (i != j)
			    {
				if (sum == 0)
				    matrix(i,j) = 0;
				else
				    matrix(i,j) = matrix(i,j) / sum * (1000 - _same * 10);
			    }
		    }
	    }
    }

private:
    int _initG;
    double _same;
};

template <typename EOT>
double binary_value(const EOT& _chrom)
{
    double sum = 0.0;
    for (unsigned i=0; i<_chrom.size(); i++)
	sum += _chrom[i];
    return sum;
}

typedef eoBit<double> EOT;

queue< pair< size_t, eoPop<EOT> > > imm;
queue< pair< size_t, eoPop<EOT> > > em;
condition_variable cv;
mutex cv_m;

bool ready = false;

int get_sum( eoPop<EOT>& pop, mpi::communicator& topology )
{
    int size = pop.size();
    int sum = -1;
    for (int i = 0; i < topology.size(); ++i)
	{
	    MPI_Reduce( &size, &sum, 1, MPI_INT, MPI_SUM, i, topology.mpi_comm() );
	}
    return sum;
}

void print_sum( eoPop<EOT>& pop, mpi::communicator& topology )
{
    cout << "F" << pop.size() << " "; cout.flush();
    int sum = get_sum(pop, topology);
    if ( 0 == topology.rank() )
    	{
	    cout << "sum: " << sum << endl; cout.flush();
	}
}

void mainThread( eoGenContinue<EOT>& cont, eoPop<EOT>& pop, vector< double >& vecProba, mpi::communicator& topology, vector<int>& to, vector<int>& from )
{
    while ( !ready )
	{

	    // emigrate
	    {
		// for (size_t i = 0; i < to.size(); ++i)
		//     {
		// 	em.push( make_pair(to[i], pop) );
		//     }
		// pop.clear();

		vector< eoPop<EOT> > pops( to.size(), eoPop<EOT>() );

		for (size_t i = 0; i < pop.size(); ++i)
		    {
			double s = 0;
			int r = rng.rand() % 1000 + 1;

			size_t j;
			for ( j = 0; j < to.size() && r > s; ++j )
			    {
				s += vecProba[j];
			    }
			--j;

			pops[j].push_back(pop[i]);
		    }

		for ( size_t i = 0; i < to.size(); ++i )
		    {
			em.push( make_pair( to[i], pops[i] ) );
			// cout << "e" << pops[i].size() << " "; cout.flush();
		    }

		pop.clear();
	    }

	    // immigrate
	    {
		while ( !imm.empty() )
		    {
			// eoPop<EOT> newpop = imm.front().second;
			// pop.resize( pop.size() + (newpop.size()/2) );
			// copy( newpop.begin(), newpop.end(), pop.begin() + (newpop.size()/2) );
			// pop.assign( imm.front().second.begin(), imm.front().second.end() );
			// copy( imm.front().second.begin(), imm.front().second.end(), pop.begin() );
			eoPop<EOT>& newpop = imm.front().second;
			// cout << "i" << newpop.size() << " "; cout.flush();
			for (size_t i = 0; i < newpop.size(); ++i)
			    {
				pop.push_back( newpop[i] );
			    }
			imm.pop();
		    }
	    }

	    // // cout << "f" << pop.size() << " "; cout.flush();

	    this_thread::sleep_for(chrono::milliseconds(10));
	    // this_thread::sleep_for(chrono::seconds(1));
	    // cerr << "Notifying...\n";
	    cv.notify_all();

	}

    // immigrate
    {
	while ( !imm.empty() )
	    {
		// eoPop<EOT> newpop = imm.front().second;
		// pop.resize( pop.size() + (newpop.size()/2) );
		// copy( newpop.begin(), newpop.end(), pop.begin() + (newpop.size()/2) );
		// pop.assign( imm.front().second.begin(), imm.front().second.end() );
		// copy( imm.front().second.begin(), imm.front().second.end(), pop.begin() );
		eoPop<EOT>& newpop = imm.front().second;
		// cout << "i" << newpop.size() << " "; cout.flush();
		for (size_t i = 0; i < newpop.size(); ++i)
		    {
			pop.push_back( newpop[i] );
		    }
		imm.pop();
	    }
    }
}

void sendRecvThread( eoGenContinue<EOT>& cont, eoPop<EOT>& pop, vector< double >& vecProba, mpi::communicator& topology, vector<int>& to, vector<int>& from )
{
    while ( cont( pop ) )
	{
	    unique_lock<mutex> lk(cv_m);

	    // cerr << "Waiting... \n";
	    cv.wait(lk);
	    // cv.wait_for(lk, chrono::milliseconds(10));
	    // cerr << "...finished waiting.\n";

	    mpi::requests reqs;
	    size_t count = 0;

	    while ( !em.empty() )
		{
		    // cout << "s" << em.front().second.size() << " "; cout.flush();
		    int dst = em.front().first;
		    eoserial::Object* obj = em.front().second.pack();
		    stringstream ss;
		    obj->print( ss );
		    delete obj;

		    const string& str = ss.str();

		    int size = str.size() + 1;
		    int pos = 0;
		    vector< char > buf( size, 0 );

		    MPI_Pack( &size, 1, MPI_INT, buf.data(), buf.size(), &pos, topology.mpi_comm() );
		    buf.resize(size + pos);

		    MPI_Pack( (char*)str.c_str(), size, MPI_CHAR, buf.data(), buf.size(), &pos, topology.mpi_comm() );

		    MPI_Request req;
		    MPI_Isend( buf.data(), pos, MPI_PACKED, dst, 0, topology.mpi_comm(), &req );
		    reqs.push_back(req);

		    em.pop();
		}

	    vector< vector< char > > from_buf( from.size(), vector< char >( 10000, 0 ) );
	    vector< int > from_pos( from.size(), 0 );

	    for (size_t i = 0; i < from.size(); ++i)
		{
		    MPI_Request req;
		    MPI_Irecv( from_buf[i].data(), from_buf[i].size(), MPI_PACKED, from[i], 0, topology.mpi_comm(), &req);
		    reqs.push_back(req);
		}

	    vector<MPI_Status> status(reqs.size());
	    MPI_Waitall(reqs.size(), reqs.data(), status.data());

	    for (size_t i = 0; i < from.size(); ++i)
		{
		    int size = -1;
		    MPI_Unpack( from_buf[i].data(), from_buf[i].size(), &from_pos[i], &size, 1, MPI_INT, topology.mpi_comm() );
		    string str(size, 0);
		    MPI_Unpack( from_buf[i].data(), from_buf[i].size(), &from_pos[i], (char*)str.data(), size, MPI_CHAR, topology.mpi_comm() );

		    eoserial::Object* obj = eoserial::Parser::parse( str );
		    eoPop<EOT> newpop;
		    newpop.unpack( obj );
		    // cout << "u" << newpop.size() << " "; cout.flush();
		    imm.push( make_pair(from[i], newpop) );
		    delete obj;
		}
	}
    ready = true;
}

int main(int argc, char *argv[])
{
    // Node::init( ac, av, MPI_THREAD_MULTIPLE );
    Node::init( argc, argv/*, MPI_THREAD_SERIALIZED*/ );
    eoParser parser(argc, argv);

    // S
    unsigned seed = parser.createParam(unsigned(time(0)), "seed", "Random number seed", 'S').value(); // will be in default section General
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
    // unsigned nbMigrants = parser.createParam(unsigned(popSize), "nbMigrants", "Number of migrants",'n', "Migration Policy").value();
    // N
    unsigned manGeneration = parser.createParam(unsigned(1), "manGeneration", "Migration at N generations",'N', "Migration Policy").value();
    // i
    // size_t nbIslands = parser.createParam(size_t(4), "nbIslands", "Number of islands",'i', "Island Model").value();

    make_parallel( parser );
    make_verbose(parser);
    make_help(parser);

    rng.reseed(seed);

    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    // mpi_debug();

    if ( ALL < 4 )
    	{
    	    if ( 0 == RANK )
    		{
    		    cerr << "Needs at least 4 processes to be launched!" << endl;
    		}
    	    return 0;
    	}

    eoGenContinue<EOT> cont( maxGen );

    eoUniformGenerator<typename EOT::AtomType> gen(0, 2);
    eoInitFixedLength<EOT> init(dimSize, gen);

    eoPop<EOT> pop;
    pop.append(popSize, init);

    eoEvalFuncPtr<EOT> eval(binary_value);
    apply<EOT>(eval, pop);

    Intracomm comm = COMM_WORLD;

    // mpi::ring topology;
    mpi::complete topology;

    const size_t GALL = topology.size();
    const size_t GRANK = topology.rank();

    if ( 0 == GRANK )
	{
	    // topology.print();
	    topology.test();
	}

    size_t size = topology.neighbors_count();

    vector<int> to = topology.to();
    vector<int> from = topology.from();

    cout << to.size() << " " << from.size() << endl; cout.flush();

    topology.barrier();

    MigrationMatrix probabilities( /*nbIslands*/ GALL );
    InitMatrix initmatrix;

    initmatrix( probabilities );

    if ( 0 == GRANK )
    	{
    	    cout << probabilities; cout.flush();
    	}

    topology.barrier();

    print_sum(pop, topology);

    topology.barrier();

    vector< double > vecProba = probabilities(GRANK);

    int sum = get_sum(pop, topology);

    do
    	{
	    thread t1(mainThread,	ref(cont), ref(pop), ref(vecProba), ref(topology), ref(to), ref(from));
	    thread t2(sendRecvThread,	ref(cont), ref(pop), ref(vecProba), ref(topology), ref(to), ref(from));

	    t1.join();
	    t2.join();

	    topology.barrier();

	    cont.reset();
	    ready = false;

	    sum = get_sum(pop, topology);
	    if (0 == GRANK)
		{
		    cout << "S: " << sum << " "; cout.flush();
		}
	}
    while ( sum != GALL*popSize );

    print_sum(pop, topology);

    return 0;
}
