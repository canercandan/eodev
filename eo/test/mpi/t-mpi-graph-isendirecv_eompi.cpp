#include <stdexcept>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <mpi.h>

#include <thread>
#include <mutex>
#include <atomic>
#include <queue>

#include <eompi.h>
#include <eo>
#include <ga.h>

#include "t-mpi-common.h"

using namespace std;
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

    // eoPop<EOT> pop;
    // pop.append(popSize, init);

    eoEvalFuncPtr<EOT> eval(binary_value);
    // apply<EOT>(eval, pop);

    MPI::Intracomm comm = MPI::COMM_WORLD;


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

    cout << to.size() << " " << from.size() << endl;

    MPI_Barrier(topology.mpi_comm());

    vector< eoPop<EOT> > to_data( to.size(), eoPop<EOT>(popSize, init) );
    vector< eoPop<EOT> > from_data( from.size(), eoPop<EOT>() );

    // vector< SerializableBase< int > > to_data( to.size(), 42 );
    // vector< SerializableBase< int > > from_data( from.size(), 0 );

    // vector< vector< int > > to_data( to.size(), vector<int>( popSize, 42 ) );
    // vector< vector< int > > from_data( from.size(), vector<int>( popSize, 0 ) );

    for (int j = 0; j < maxGen; ++j)
	{

	    mpi::requests reqs;

	    for (size_t i = 0; i < to.size(); ++i)
		{
		    eoserial::Object* obj = to_data[i].pack();
		    std::stringstream ss;
		    obj->print( ss );
		    delete obj;

		    const std::string& str = ss.str();

		    int size = str.size() + 1;
		    int pos = 0;
		    std::vector< char > buf( size, 0 );

		    MPI_Pack( &size, 1, MPI_INT, buf.data(), buf.size(), &pos, topology.mpi_comm() );
		    buf.resize(size + pos);

		    MPI_Pack( (char*)str.c_str(), size, MPI_CHAR, buf.data(), buf.size(), &pos, topology.mpi_comm() );

		    MPI_Request req;
		    MPI_Isend( buf.data(), pos, MPI_PACKED, to[i], 0, topology.mpi_comm(), &req );

		    // MPI_Send( buf.data(), pos, MPI_PACKED, to[i], 0, topology.mpi_comm() );

		    reqs.push_back(req);

		    // mpi::requests r = topology.isend(to[i], to_data[i]);
		    // for (int i = 0; i < r.size(); ++i)
		    // 	{
		    // 	    reqs.push_back(r[i]);
		    // 	}
		}

	    vector< vector< char > > from_buf( from.size(), vector< char >( 10000, 0 ) );
	    vector< int > from_pos( from.size(), 0 );

	    for (size_t i = 0; i < from.size(); ++i)
		{
		    MPI_Request req;
		    MPI_Irecv( from_buf[i].data(), from_buf[i].size(), MPI_PACKED, from[i], 0, topology.mpi_comm(), &req);
		    reqs.push_back(req);

		    // MPI_Status stat_prob;
		    // MPI_Probe( MPI_ANY_SOURCE, 0, topology.mpi_comm(), &stat_prob );

		    // MPI_Status stat;
		    // MPI_Recv( from_buf[i].data(), from_buf[i].size(), MPI_PACKED, stat_prob.MPI_SOURCE, 0, topology.mpi_comm(), &stat );

		    // mpi::requests r = topology.irecv(from[i], from_data[i]);
		    // for (int i = 0; i < r.size(); ++i)
		    // 	{
		    // 	    reqs.push_back(r[i]);
		    // 	}
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
		    from_data[i].unpack( obj );
		    delete obj;
		}

	    for (size_t i = 0; i < from.size(); ++i)
		{
		    // copy(from_data[i].begin(), from_data[i].end(), ostream_iterator< int >(cout, " "));
		    cout << from_data[i].size() << " ";
		    cout.flush();
		}

	}

    return 0;
}
