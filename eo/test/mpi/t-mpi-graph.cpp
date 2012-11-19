#include <stdexcept>
#include <iostream>
#include <iterator>
#include <vector>
#include <assert.h>
#include <mpi.h>

using namespace MPI;
using namespace std;

vector<int> getNeighbors(Graphcomm& graph_comm)
{
    const size_t GRANK = graph_comm.Get_rank();
    const size_t size = graph_comm.Get_neighbors_count(GRANK);
    vector<int> neighbors(size);
    graph_comm.Get_neighbors(GRANK, size, neighbors.data());
    return neighbors;
}

void printNeighbors(vector<int>& neighbors)
{
    copy(neighbors.begin(), neighbors.end(), ostream_iterator< int >(cout, " "));
    cout.flush();
    cout << endl;
    cout.flush();
}

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

    /*
      Process    Neighbors
      0          1, 3
      1          0
      2          3
      3          0, 2
    */

    // int nnode = 4;
    // int index[] = {2, 3, 4, 6};
    // int edges[] = {1, 3, 0, 3, 0, 2};
    // bool reorder = true;
    // Graphcomm graph_comm = comm.Create_graph(nnode, index, edges, reorder);

    /*
      Process    Neighbors
      0          0, 1, 3
      1          0, 1
      2          2, 3
      3          0, 2, 3
    */

    int nnode = 4;
    int index[] = {3, 5, 7, 10};
    int edges[] = {0, 1, 3, 0, 1, 2, 3, 0, 2, 3};
    bool reorder = true;
    Graphcomm graph_comm = comm.Create_graph(nnode, index, edges, reorder);

    const size_t GALL = graph_comm.Get_size();
    const size_t GRANK = graph_comm.Get_rank();

    size_t size = graph_comm.Get_neighbors_count(GRANK);

    vector<int> v = getNeighbors(graph_comm);
    printNeighbors(v);

    assert( graph_comm.Get_topology() == MPI_GRAPH );

    Finalize();
    return 0;
}
