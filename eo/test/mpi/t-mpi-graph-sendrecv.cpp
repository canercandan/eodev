#include <stdexcept>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <mpi.h>

using namespace MPI;
using namespace std;

class Topology
{
public:
    Topology( Intracomm& comm = COMM_WORLD,
	      int nnode = COMM_WORLD.Get_size(),
	      int nedge = COMM_WORLD.Get_size(),
	      bool reorder = true )
	: _index(nnode), _edges(nedge)
    {}

    Topology( Graphcomm& g )
    {
	_graph = g;
    }

    virtual ~Topology(){}

    operator Graphcomm() { return _graph; }
    Graphcomm& comm() { return _graph; }

    inline size_t size() { return _graph.Get_size(); }
    inline size_t rank() { return _graph.Get_rank(); }

    void print()
    {
	copy(_index.begin(), _index.end(), ostream_iterator< int >(cout, " "));
	copy(_edges.begin(), _edges.end(), ostream_iterator< int >(cout, " "));
	cout << endl;
    }

    void test()
    {
	assert( _graph.Get_topology() == MPI_GRAPH );
    }

    inline int neighbors_count() { return _graph.Get_neighbors_count(rank()); }
    inline int neighbors_count(int rank) { return _graph.Get_neighbors_count(rank); }

    vector<int> to(int rank)
    {
	const size_t size = neighbors_count(rank);
	vector<int> neighbors(size);
	_graph.Get_neighbors(rank, size, neighbors.data());
	return neighbors;
    }

    inline vector<int> to() { return to(rank()); }

    vector<int> from(int rank)
    {
	vector<int> vfrom;

	for (int i = 0; i < size(); ++i)
	    {
		vector<int> v = to(i);
		if ( find(v.begin(), v.end(), rank) != v.end() )
		    {
			vfrom.push_back(i);
		    }
	    }

	copy(vfrom.begin(), vfrom.end(), ostream_iterator< int >(cout, " "));
	cout.flush();
	cout << endl;
	cout.flush();

	return vfrom;
    }

    inline vector<int> from() { return from(rank()); }

protected:
    vector<int> _index;
    vector<int> _edges;
    Graphcomm _graph;
};

class Complete : public Topology
{
public:
    Complete( Intracomm& comm = COMM_WORLD, int nnode = COMM_WORLD.Get_size(), bool reorder = true )
	: Topology(comm, nnode, nnode*nnode, reorder)
    {
	for (int i = 0; i < nnode; ++i)
	    {
		_index[i] = (i+1) * nnode;

		for (int j = 0; j < nnode; ++j)
		    {
			_edges[ j + (i*nnode) ] = j;
		    }
	    }

	_graph = comm.Create_graph(nnode, _index.data(), _edges.data(), reorder);
    }
};

class Ring : public Topology
{
public:
    Ring( Intracomm& comm = COMM_WORLD, int nnode = COMM_WORLD.Get_size(), bool reorder = true )
	: Topology(comm, nnode, nnode, reorder)
    {
	for (int i = 0; i < nnode; ++i)
	    {
		_index[i] = i+1;
		_edges[i] = (i+1) % nnode;
	    }

	_graph = comm.Create_graph(nnode, _index.data(), _edges.data(), reorder);
    }

    virtual ~Ring() {}
};

int main(int argc, char *argv[])
{
    Init(argc, argv);

    Intracomm comm = COMM_WORLD;

    const size_t ALL = comm.Get_size();
    const size_t RANK = comm.Get_rank();

    if ( ALL < 4 )
	{
	    if ( 0 == RANK )
		{
		    cerr << "Needs at least 4 processes to be launched!" << endl;
		}
	    Finalize();
	    return 0;
	}

    // Complete topology;
    Ring topology;

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

    topology.comm().Barrier();

    vector< vector<int> > to_data(to.size(), vector<int>(10, 42) );
    vector< vector<int> > from_data(from.size(), vector<int>(10, 0) );

    vector<Request> reqs;

    for (size_t i = 0; i < to.size(); ++i)
    	{
    	    topology.comm().Send(to_data[i].data(), to_data[i].size(), INT, to[i], 0);
    	}

    for (size_t i = 0; i < from.size(); ++i)
    	{
    	    topology.comm().Recv(from_data[i].data(), from_data[i].size(), INT, from[i], 0);
    	}

    for (size_t i = 0; i < from_data.size(); ++i)
    	{
    	    copy(from_data[i].begin(), from_data[i].end(), ostream_iterator< int >(cout, " "));
    	    cout << endl;
    	    cout.flush();
    	}

    Finalize();
    return 0;
}
