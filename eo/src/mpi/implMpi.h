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
# ifndef __EO_IMPL_MPI_HPP__
# define __EO_IMPL_MPI_HPP__

# include <mpi.h>
# include <serial/eoSerial.h>

/**
 * This namespace contains reimplementations of some parts of the Boost::MPI API in C++, so as to be used in
 * EO without any dependance to Boost. Historically, EO's parallelization module used the
 * boost library to add a layer over MPI. After having noticed that just some functions
 * were really used, we decided to reimplement our own C++-like implementation of MPI.
 *
 * Because the Boost::MPI API is really clean, we reused it in this module. However, all
 * the functions of Boost::MPI were not used, hence a subset of the API is reused. For
 * instance, users can just send integer, std::string or eoserial::Persistent objects;
 * furthermore, only eoserial::Persistent objects can sent in a table.
 *
 * The documentation of the functions is exactly the same as the official Boost::MPI
 * documentation. You can find it on www.boost.org/doc/libs/1_49_0/doc/html/mpi/
 * The entities are here shortly described, if you need further details, don't hesitate
 * to visit the boost URL.
 */

namespace mpi
{
    /**
     * @ingroup Parallel
     * @{
     */

    /**
     * @brief Constant indicating that a message can come from any process.
     */
    extern const int any_source;

    /**
     * @brief Constant indicating that a message can come from any tag (channel).
     */
    extern const int any_tag;

    /**
     * @brief Wrapper class to have a MPI environment.
     *
     * Instead of calling MPI_Init and MPI_Finalize, it is only necessary to instantiate
     * this class once, in the global context.
     */
    class environment
    {
        public:

        /**
         * @brief Inits MPI context without arguments.
         *
         */
        environment();

        /**
         * @brief Inits MPI context.
         *
         * @param argc Number of params in command line (same as one in main)
         * @param argv Strings containing params (same as one in main)
         */
        environment(int argc, char**argv);

        /**
         * @brief Inits MPI context with thread support.
         *
         * @param argc Number of params in command line (same as one in main)
         * @param argv Strings containing params (same as one in main)
	 * @param required Desired level of thread support (integer): : MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE
         */
        environment(int argc, char**argv, int required);

        /**
         * @brief Closes MPI context.
         */
        ~environment();
    };

    /**
     * @brief Wrapper class for MPI_Status
     *
     * Consists only in a C++ wrapper class, giving getters on status attributes.
     */
    class status
    {
        public:

        /**
         * @brief Converts a MPI_Status into a status.
         */
        status( const MPI_Status & s, bool flag = true );

        /**
         * @brief Returns the tag of the associated communication.
         */
        int tag() { return _tag; }

        /**
         * @brief Indicates which error number we obtained in the associated communication.
         */
        int error() { return _error; }

        /**
         * @brief Returns the MPI rank of the source of the associated communication.
         */
        int source() { return _source; }

        private:
            int _source;
            int _tag;
            int _error;
	    bool _flag;
    };

    /**
     * @brief Wrapper type for MPI_Request
     */
    typedef MPI_Request request;

    /**
     * @brief Type for a set of request
     */
    typedef std::vector< request > requests;

    /**
     * @brief Main object, used to send / receive messages, get informations about the rank and the size of the world,
     * etc.
     */
    class communicator
    {
        public:

        /**
         * Creates the communicator, using the whole world as a MPI_Comm.
         *
         * @todo Allow the user to precise which MPI_Comm to use
         */
        communicator(MPI_Comm _comm = MPI_COMM_WORLD);

        ~communicator();

        /**
         * @brief Returns the MPI rank of the current process.
         */
        int rank();

        /**
         * @brief Returns the size of the MPI cluster.
         */
        int size();

        /**
         * @brief Returns the handled MPI communicator.
         */
        MPI_Comm mpi_comm();

        /*
         * SEND / RECV INT
         */

        /**
         * @brief Sends an integer to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param n The integer to send
         */
        void send( int dest, int tag, int n );

        /**
         * @brief Sends an integer to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param n The integer to send
         */
        void send( int dest, int n );

        /*
         * @brief Receives an integer from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param n Where to save the received integer
         */
        void recv( int src, int tag, int& n );

        /*
         * @brief Receives an integer from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param n Where to save the received integer
         */
        void recv( int src, int& n );

        /*
         * ISEND / IRECV INT
         */

        /**
         * @brief Non-blocking method to send an integer to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param n The integer to send
         */
        requests isend( int dest, int tag, int n );

        /**
         * @brief Non-blocking method to send an integer to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param n The integer to send
         */
        requests isend( int dest, int n );

        /*
         * @brief Non-blocking method to receive an integer from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param n Where to save the received integer
         */
        requests irecv( int src, int tag, int& n );

        /*
         * @brief Non-blocking method to receive an integer from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param n Where to save the received integer
         */
        requests irecv( int src, int& n );

        /*
         * SEND / RECV BOOL
         */

        /**
         * @brief Sends a boolean to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param b The boolean to send
         */
        void send( int dest, int tag, bool b );

        /**
         * @brief Sends a boolean to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param b The boolean to send
         */
        void send( int dest, bool b );

        /*
         * @brief Receives a boolean from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param b Where to save the received boolean
         */
        void recv( int src, int tag, bool& b );

        /*
         * @brief Receives a boolean from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param b Where to save the received boolean
         */
        void recv( int src, bool& b );

        /*
         * ISEND / IRECV BOOL
         */

        /**
         * @brief Non-blocking method to send a boolean to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param b The boolean to send
         */
        requests isend( int dest, int tag, bool b );

        /**
         * @brief Non-blocking method to send a boolean to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param b The boolean to send
         */
        requests isend( int dest, bool b );

        /*
         * @brief Non-blocking method to receive a boolean from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param b Where to save the received boolean
         */
        requests irecv( int src, int tag, bool& b );

        /*
         * @brief Non-blocking method to receive a boolean from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param b Where to save the received boolean
         */
        requests irecv( int src, bool& b );

        /*
         * SEND / RECV STRING
         */

        /**
         * @brief Sends a string to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param str The std::string to send
         */
        void send( int dest, int tag, const std::string& str );

        /**
         * @brief Sends a string to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param str The std::string to send
         */
        void send( int dest, const std::string& str );

        /*
         * @brief Receives a string from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param std::string Where to save the received string
         */
        void recv( int src, int tag, std::string& str );

        /*
         * @brief Receives a string from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param std::string Where to save the received string
         */
        void recv( int src, std::string& str );

        /*
         * ISEND / IRECV STRING
         */

        /**
         * @brief Non-blocking method to send a string to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param str The std::string to send
         */
        requests isend( int dest, int tag, const std::string& str );

        /**
         * @brief Non-blocking method to send a string to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param str The std::string to send
         */
        requests isend( int dest, const std::string& str );

        /*
         * @brief Non-blocking method to receive a string from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param std::string Where to save the received string
         */
        requests irecv( int src, int tag, std::string& str );

        /*
         * @brief Non-blocking method to receive a string from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param std::string Where to save the received string
         */
        requests irecv( int src, std::string& str );

        /*
         * SEND / RECV Objects
         */

        /**
         * @brief Sends an eoserial::Persistent to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param persistent The object to send (it must absolutely implement eoserial::Persistent)
         */
        void send( int dest, int tag, const eoserial::Persistent & persistent );

        /**
         * @brief Sends an eoserial::Persistent to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param persistent The object to send (it must absolutely implement eoserial::Persistent)
         */
        void send( int dest, const eoserial::Persistent & persistent );

        /**
         * @brief Sends an array of eoserial::Persistent to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param table The array of eoserial::Persistent objects
         * @param size The number of elements to send (no check is done, the user has to be sure that the size won't
         * overflow!)
         */
        template< class T >
        void send( int dest, int tag, T* table, int size )
        {
            // Puts all the values into an array
            eoserial::Array* array = new eoserial::Array;

            for( int i = 0; i < size; ++i )
            {
                array->push_back( table[i].pack() );
            }

            // Encapsulates the array into an object
            eoserial::Object* obj = new eoserial::Object;
            obj->add( "array", array );
            std::stringstream ss;
            obj->print( ss );
            delete obj;

            // Sends the object as a string
            send( dest, tag, ss.str() );
        }

        /**
         * @brief Sends an array of eoserial::Persistent to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param table The array of eoserial::Persistent objects
         * @param size The number of elements to send (no check is done, the user has to be sure that the size won't
         * overflow!)
         */
        template< class T >
        void send( int dest, T* table, int size )
        {
	    send(dest, 0, table, size);
	}

        /*
         * @brief Receives an eoserial::Persistent object from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param persistent Where to unpack the serialized object?
         */
        void recv( int src, int tag, eoserial::Persistent & persistent );

        /*
         * @brief Receives an eoserial::Persistent object from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param persistent Where to unpack the serialized object?
         */
        void recv( int src, eoserial::Persistent & persistent );

        /*
         * @brief Receives an array of eoserial::Persistent from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param table The table in which we're saving the received objects. It must have been allocated by the user,
         * as no allocation is performed here.
         * @param size The number of elements to receive (no check is done, the user has to be sure that the size won't
         * overflow!)
         */
        template< class T >
        void recv( int src, int tag, T* table, int size )
        {
            // Receives the string which contains the object
            std::string asText;
            recv( src, tag, asText );

            // Parses the object and retrieves the table
            eoserial::Object* obj = eoserial::Parser::parse( asText );
            eoserial::Array* array = static_cast<eoserial::Array*>( (*obj)["array"] );

            // Retrieves all the values from the array
            for( int i = 0; i < size; ++i )
            {
                eoserial::unpackObject( *array, i, table[i] );
            }
            delete obj;
        }

        /*
         * @brief Receives an array of eoserial::Persistent from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param table The table in which we're saving the received objects. It must have been allocated by the user,
         * as no allocation is performed here.
         * @param size The number of elements to receive (no check is done, the user has to be sure that the size won't
         * overflow!)
         */
        template< class T >
        void recv( int src, T* table, int size )
        {
	    recv(src, 0, table, size);
	}

        /*
         * ISEND / IRECV Objects
         */

        /**
         * @brief Non-blocking method to send an eoserial::Persistent to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param persistent The object to send (it must absolutely implement eoserial::Persistent)
         */
        requests isend( int dest, int tag, const eoserial::Persistent & persistent );

        /**
         * @brief Non-blocking method to send an eoserial::Persistent to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param persistent The object to send (it must absolutely implement eoserial::Persistent)
         */
        requests isend( int dest, const eoserial::Persistent & persistent );

        /**
         * @brief Non-blocking method to send an array of eoserial::Persistent to dest on channel "tag".
         *
         * @param dest MPI rank of the receiver
         * @param tag MPI tag of message
         * @param table The array of eoserial::Persistent objects
         * @param size The number of elements to send (no check is done, the user has to be sure that the size won't
         * overflow!)
         */
        template< class T >
        requests isend( int dest, int tag, T* table, int size )
        {
            // Puts all the values into an array
            eoserial::Array* array = new eoserial::Array;

            for( int i = 0; i < size; ++i )
            {
                array->push_back( table[i].pack() );
            }

            // Encapsulates the array into an object
            eoserial::Object* obj = new eoserial::Object;
            obj->add( "array", array );
            std::stringstream ss;
            obj->print( ss );
            delete obj;

            // Sends the object as a string
            return isend( dest, tag, ss.str() );
        }

        /**
         * @brief Non-blocking method to send an array of eoserial::Persistent to dest on default channel "commands".
         *
         * @param dest MPI rank of the receiver
         * @param table The array of eoserial::Persistent objects
         * @param size The number of elements to send (no check is done, the user has to be sure that the size won't
         * overflow!)
         */
        template< class T >
        requests isend( int dest, T* table, int size )
        {
	    return isend(dest, 0, table, size);
	}

        /*
         * @brief Non-blocking method to receive an eoserial::Persistent object from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param persistent Where to unpack the serialized object?
         */
        requests irecv( int src, int tag, eoserial::Persistent & persistent );

        /*
         * @brief Non-blocking method to receive an eoserial::Persistent object from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param persistent Where to unpack the serialized object?
         */
        requests irecv( int src, eoserial::Persistent & persistent );

        /*
         * @brief Non-blocking method to receive an array of eoserial::Persistent from src on channel "tag".
         *
         * @param src MPI rank of the sender
         * @param tag MPI tag of message
         * @param table The table in which we're saving the received objects. It must have been allocated by the user,
         * as no allocation is performed here.
         * @param size The number of elements to receive (no check is done, the user has to be sure that the size won't
         * overflow!)
         */
        template< class T >
        requests irecv( int src, int tag, T* table, int size )
        {
            // Receives the string which contains the object
            std::string asText;
            requests req = irecv( src, tag, asText );

            // Parses the object and retrieves the table
            eoserial::Object* obj = eoserial::Parser::parse( asText );
            eoserial::Array* array = static_cast<eoserial::Array*>( (*obj)["array"] );

            // Retrieves all the values from the array
            for( int i = 0; i < size; ++i )
            {
                eoserial::unpackObject( *array, i, table[i] );
            }
            delete obj;

	    return req;
        }

        /*
         * @brief Non-blocking method to receive an array of eoserial::Persistent from src on default channel "commands".
         *
         * @param src MPI rank of the sender
         * @param table The table in which we're saving the received objects. It must have been allocated by the user,
         * as no allocation is performed here.
         * @param size The number of elements to receive (no check is done, the user has to be sure that the size won't
         * overflow!)
         */
        template< class T >
        requests irecv( int src, T* table, int size )
        {
	    return irecv(src, 0, table, size);
	}

        /*
         * Other methods
         */

        /**
         * @brief Wrapper for MPI_Probe
         *
         * Waits for a message to come from process having rank src, on the channel
         * tag.
         *
         * @param src MPI rank of the sender (any_source if it can be any sender)
         * @param tag MPI tag of the expected message (any_tag if it can be any tag)
         */
        status probe( int src = any_source, int tag = any_tag );


        /**
         * @brief Wrapper for MPI_Iprobe
         *
         * Waits for a message to come from process having rank src, on the channel
         * tag.
         *
         * @param src MPI rank of the sender (any_source if it can be any sender)
         * @param tag MPI tag of the expected message (any_tag if it can be any tag)
         */
        status iprobe( int src = any_source, int tag = any_tag );

        /**
         * @brief Wrapper for MPI_Barrier
         *
         *
         */
        void barrier();

        protected:
	    MPI_Comm _comm;

        private:
            int _rank;
            int _size;

            char* _buf; // temporary buffer for sending and receiving strings. Avoids reallocations
            int _bufsize; // size of the above temporary buffer
    };

    /**
     * @brief Wrapper for MPI_Bcast
     *
     * Broadcasts an integer value on the communicator comm, from the process having the MPI rank root.
     *
     * @param comm The communicator on which to broadcast
     * @param value The integer value to send
     * @param root The MPI rank of the broadcaster
     *
     * @todo Actually comm isn't used and broadcast is performed on the whole MPI_COMM_WORLD. TODO: Use comm instead
     */
    void broadcast( communicator & comm, int value, int root );


    class topology : public communicator
    {
    public:
	topology();
	// virtual ~topology();
    };


    class graph_topology : public topology
    {
    public:
	graph_topology( MPI_Comm comm = MPI_COMM_WORLD,
			int nnode = /*COMM_WORLD.Get_size()*/ 0,
			int nedge = /*COMM_WORLD.Get_size()*/ 0,
			bool reorder = true );
	// virtual ~graph_topology();

	void print();
	void test();

	int neighbors_count(int rank);
	int neighbors_count();

	std::vector<int> to(int rank);
	std::vector<int> to();
	std::vector<int> from(int rank);
	std::vector<int> from();

    protected:
	std::vector<int> _index;
	std::vector<int> _edges;
    };


    class ring : public graph_topology
    {
    public:
	ring(MPI_Comm comm = MPI_COMM_WORLD, int nnode = MPI::COMM_WORLD.Get_size(), bool reorder = true);
	// virtual ~ring();
    };


    class complete : public graph_topology
    {
    public:
	complete(MPI_Comm comm = MPI_COMM_WORLD, int nnode = MPI::COMM_WORLD.Get_size(), bool reorder = true);
	// virtual ~complete();
    };


    /**
     * @}
     */
} // namespace mpi

# endif //__EO_IMPL_MPI_HPP__
