/*
(c) Thales group, 2012

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation;
    version 2 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
Contact: http://eodev.sourceforge.net

Authors:
    Benjamin Bouvier <benjamin.bouvier@gmail.com>
*/

#include <iterator>
#include <assert.h>
#include <algorithm>

# include "implMpi.h"

namespace mpi
{
    const int any_source = MPI_ANY_SOURCE;
    const int any_tag = MPI_ANY_TAG;

    environment::environment()
    {
	MPI::Init();
    }

    environment::environment(int argc, char**argv)
    {
        MPI_Init(&argc, &argv);
    }

    environment::environment(int argc, char**argv, int required)
    {
	int provided;
	MPI_Init_thread(&argc, &argv, required, &provided);
    }

    environment::~environment()
    {
        MPI_Finalize();
    }

    status::status( const MPI_Status & s, bool flag /*= true*/ )
    {
        _source = s.MPI_SOURCE;
        _tag = s.MPI_TAG;
        _error = s.MPI_ERROR;
	_flag = flag;
    }

    communicator::communicator(MPI_Comm comm /*= MPI_COMM_WORLD*/)
    {
	_comm = comm;

        _rank = -1;
        _size = -1;

        _buf = 0;
        _bufsize = -1;
    }

    communicator::~communicator()
    {
        if( _buf )
        {
            delete _buf;
            _buf = 0;
        }
    }

    int communicator::rank()
    {
        if ( _rank == -1 )
        {
            MPI_Comm_rank( _comm, &_rank );
        }
        return _rank;
    }

    int communicator::size()
    {
        if ( _size == -1 )
        {
            MPI_Comm_size( _comm, &_size );
        }
        return _size;
    }

    MPI_Comm communicator::mpi_comm()
    {
        return _comm;
    }

    /*
     * SEND / RECV INT
     */
    void communicator::send( int dest, int tag, int n )
    {
        MPI_Send( &n, 1, MPI_INT, dest, tag, _comm );
    }

    void communicator::send( int dest, int n )
    {
	send(dest, 0, n);
    }

    void communicator::recv( int src, int tag, int& n )
    {
        MPI_Status stat;
        MPI_Recv( &n, 1, MPI_INT, src, tag, _comm , &stat );
    }

    void communicator::recv( int src, int& n )
    {
	recv(src, 0, n);
    }

    /*
     * ISEND / IRECV INT
     */
    requests communicator::isend( int dest, int tag, int n )
    {
	MPI_Request req;
        MPI_Isend( &n, 1, MPI_INT, dest, tag, _comm, &req );
	int flag;
	MPI_Status stat;
	MPI_Request_get_status( req, &flag, &stat );
	requests reqs;
	reqs.push_back(req);
	return reqs;
    }

    requests communicator::isend( int dest, int n )
    {
	return isend(dest, 0, n);
    }

    requests communicator::irecv( int src, int tag, int& n )
    {
	MPI_Request req;
        MPI_Irecv( &n, 1, MPI_INT, src, tag, _comm , &req );
	int flag;
	MPI_Status stat;
	MPI_Request_get_status( req, &flag, &stat );
	requests reqs;
	reqs.push_back(req);
	return reqs;
    }

    requests communicator::irecv( int src, int& n )
    {
	return irecv(src, 0, n);
    }

    /*
     * SEND / RECV BOOL
     */
    void communicator::send( int dest, int tag, bool b )
    {
        MPI_Send( &b, 1, MPI::BOOL, dest, tag, _comm );
    }

    void communicator::send( int dest, bool b )
    {
	send(dest, 0, b);
    }

    void communicator::recv( int src, int tag, bool& b )
    {
        MPI_Status stat;
        MPI_Recv( &b, 1, MPI::BOOL, src, tag, _comm , &stat );
    }

    void communicator::recv( int src, bool& b )
    {
	recv(src, 0, b);
    }

    /*
     * ISEND / IRECV BOOL
     */
    requests communicator::isend( int dest, int tag, bool b )
    {
	MPI_Request req;
        MPI_Isend( &b, 1, MPI::BOOL, dest, tag, _comm, &req );
	int flag;
	MPI_Status stat;
	MPI_Request_get_status( req, &flag, &stat );
	requests reqs;
	reqs.push_back(req);
	return reqs;
    }

    requests communicator::isend( int dest, bool b )
    {
	return isend(dest, 0, b);
    }

    requests communicator::irecv( int src, int tag, bool& b )
    {
        MPI_Request req;
        MPI_Irecv( &b, 1, MPI::BOOL, src, tag, _comm, &req );
	int flag;
	MPI_Status stat;
	MPI_Request_get_status( req, &flag, &stat );
	requests reqs;
	reqs.push_back(req);
	return reqs;
    }

    requests communicator::irecv( int src, bool& b )
    {
	return irecv(src, 0, b);
    }

    /*
     * SEND / RECV STRING
     */
    void communicator::send( int dest, int tag, const std::string& str )
    {
        int size = str.size() + 1;
        send( dest, tag, size );
        MPI_Send( (char*)str.c_str(), size, MPI_CHAR, dest, tag, _comm);
    }

    void communicator::send( int dest, const std::string& str )
    {
	send(dest, 0, str);
    }

    void communicator::recv( int src, int tag, std::string& str )
    {
        int size = -1;
        MPI_Status stat;
        recv( src, tag, size );

        if( _buf == 0 )
        {
            _buf = new char[ size ];
            _bufsize = size;
        } else if( _bufsize < size )
        {
            delete [] _buf;
            _buf = new char[ size ];
            _bufsize = size;
        }
        MPI_Recv( _buf, size, MPI_CHAR, src, tag, _comm, &stat );
        str.assign( _buf );
    }

    void communicator::recv( int src, std::string& str )
    {
	recv(src, 0, str);
    }

    /*
     * ISEND / IRECV STRING
     */
    requests communicator::isend( int dest, int tag, const std::string& str )
    {
        int size = str.size() + 1;
	int pos = 0;
	std::vector< char > buf( size, 0 );

	MPI_Pack( &size, 1, MPI_INT, buf.data(), buf.size(), &pos, _comm );
	buf.resize(size + pos);

	MPI_Pack( (char*)str.c_str(), size, MPI_CHAR, buf.data(), buf.size(), &pos, _comm );

	MPI_Request req;
        MPI_Isend( (char*)buf.data(), pos, MPI_PACKED, dest, tag, _comm, &req);

	int flag;
	MPI_Status stat;
	MPI_Request_get_status( req, &flag, &stat );

	requests reqs;
	reqs.push_back(req);
	return reqs;
    }

    requests communicator::isend( int dest, const std::string& str )
    {
	return isend(dest, 0, str);
    }

    requests communicator::irecv( int src, int tag, std::string& str )
    {
        int size = -1;
        requests reqs = irecv( src, tag, size );

        if( _buf == 0 )
        {
            _buf = new char[ size ];
            _bufsize = size;
        } else if( _bufsize < size )
        {
            delete [] _buf;
            _buf = new char[ size ];
            _bufsize = size;
        }

	MPI_Request req;
        MPI_Irecv( _buf, size, MPI_CHAR, src, tag, _comm, &req );

        str.assign( _buf );
	int flag;
	MPI_Status stat;
	MPI_Request_get_status( req, &flag, &stat );
	reqs.push_back(req);
	return reqs;
    }

    requests communicator::irecv( int src, std::string& str )
    {
	return irecv(src, 0, str);
    }

    /*
     * SEND / RECV Objects
     */
    void communicator::send( int dest, int tag, const eoserial::Persistent & persistent )
    {
        eoserial::Object* obj = persistent.pack();
        std::stringstream ss;
        obj->print( ss );
        delete obj;
        send( dest, tag, ss.str() );
    }

    void communicator::send( int dest, const eoserial::Persistent & persistent )
    {
	send(dest, 0, persistent);
    }

    void communicator::recv( int src, int tag, eoserial::Persistent & persistent )
    {
        std::string asText;
        recv( src, tag, asText );
        eoserial::Object* obj = eoserial::Parser::parse( asText );
        persistent.unpack( obj );
        delete obj;
    }

    void communicator::recv( int src, eoserial::Persistent & persistent )
    {
	recv(src, 0, persistent);
    }

    /*
     * ISEND / IRECV Objects
     */
    requests communicator::isend( int dest, int tag, const eoserial::Persistent & persistent )
    {
        eoserial::Object* obj = persistent.pack();
        std::stringstream ss;
        obj->print( ss );
        delete obj;
        return isend( dest, tag, ss.str() );
    }

    requests communicator::isend( int dest, const eoserial::Persistent & persistent )
    {
	return isend(dest, 0, persistent);
    }

    requests communicator::irecv( int src, int tag, eoserial::Persistent & persistent )
    {
        std::string asText;
        requests req = irecv( src, tag, asText );
        eoserial::Object* obj = eoserial::Parser::parse( asText );
        persistent.unpack( obj );
        delete obj;
	return req;
    }

    requests communicator::irecv( int src, eoserial::Persistent & persistent )
    {
	return irecv(src, 0, persistent);
    }

    /*
     * Other methods
     */
    status communicator::probe( int src, int tag )
    {
        MPI_Status stat;
        MPI_Probe( src, tag, _comm , &stat );
        return status( stat );
    }

    status communicator::iprobe( int src, int tag )
    {
	int flag;
        MPI_Status stat;
        MPI_Iprobe( src, tag, _comm , &flag, &stat );
        return status( stat, flag );
    }

    void communicator::barrier()
    {
        MPI_Barrier( _comm );
    }

    void broadcast( communicator & comm, int value, int root )
    {
        MPI_Bcast( &value, 1, MPI_INT, root, comm.mpi_comm() );
    }

    topology::topology() {}

    graph_topology::graph_topology( MPI_Comm comm /*= MPI_COMM_WORLD*/,
				    int nnode /*= COMM_WORLD.Get_size()*/,
				    int nedge /*= COMM_WORLD.Get_size()*/,
				    bool reorder /*= true*/ )
	: _index(nnode), _edges(nedge)
    {}

    void graph_topology::print()
    {
	std::copy(_index.begin(), _index.end(), std::ostream_iterator< int >(std::cout, " "));
	std::copy(_edges.begin(), _edges.end(), std::ostream_iterator< int >(std::cout, " "));
	std::cout << std::endl;
    }

    void graph_topology::test()
    {
	int top_type;
	MPI_Topo_test( _comm, &top_type );
	assert( top_type == MPI_GRAPH );
    }

    int graph_topology::neighbors_count(int rank)
    {
	int nneighbors;
	MPI_Graph_neighbors_count(_comm, rank, &nneighbors);
	return nneighbors;
    }

    int graph_topology::neighbors_count()
    {
	return neighbors_count(rank());
    }

    std::vector<int> graph_topology::to(int rank)
    {
	const size_t size = neighbors_count(rank);
	std::vector<int> neighbors(size);
	MPI_Graph_neighbors( _comm, rank, size, neighbors.data() );
	return neighbors;
    }

    std::vector<int> graph_topology::to()
    {
	return to(rank());
    }

    std::vector<int> graph_topology::from(int rank)
    {
	std::vector<int> vfrom;

	for (int i = 0; i < size(); ++i)
	    {
		std::vector<int> v = to(i);
		if ( find(v.begin(), v.end(), rank) != v.end() )
		    {
			vfrom.push_back(i);
		    }
	    }

	return vfrom;
    }

    std::vector<int> graph_topology::from()
    {
	return from(rank());
    }

    ring::ring(MPI_Comm comm /*= MPI_COMM_WORLD*/, int nnode /*= MPI::COMM_WORLD.Get_size()*/, bool reorder /*= true*/)
	: graph_topology(comm, nnode, nnode, reorder)
    {
	for (int i = 0; i < nnode; ++i)
	    {
		_index[i] = i+1;
		_edges[i] = (i+1) % nnode;
	    }

	MPI_Graph_create(comm, nnode, _index.data(), _edges.data(), reorder, &_comm);
    }

    complete::complete(MPI_Comm comm /*= MPI_COMM_WORLD*/, int nnode /*= MPI::COMM_WORLD.Get_size()*/, bool reorder /*= true*/)
	: graph_topology(comm, nnode, nnode*nnode, reorder)
    {
	for (int i = 0; i < nnode; ++i)
	    {
		_index[i] = (i+1) * nnode;

		for (int j = 0; j < nnode; ++j)
		    {
			_edges[ j + (i*nnode) ] = j;
		    }
	    }

	MPI_Graph_create(comm, nnode, _index.data(), _edges.data(), reorder, &_comm);
    }
}
