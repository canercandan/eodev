#include <eo>
#include <eompi.h>

#include "t-mpi-common.h"

using namespace eo::mpi;

int main(void)
{
    Node::init();

    const int ALL = Node::comm().size();
    const int RANK = Node::comm().rank();

    mpi_debug();

    eoPop<EOT> pop;

    Node::comm().isend(0, pop);
    Node::comm().irecv(1, pop);

    return 0;
}
