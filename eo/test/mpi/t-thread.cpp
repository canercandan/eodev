#include <thread>
#include <mutex>
#include <queue>
#include <eo>

class Communicable
{
public:
    Communicable()
    {
	// t1 = std::thread( &Communicable::pack, this );
	// t2 = std::thread( &Communicable::unpack, this );
    }

    // Communicable() : t1(&Communicable::pack, this), t2(&Communicable::unpack, this) {}

    virtual ~Communicable()
    {
	t1.join();
	t2.join();
    }

    virtual void pack() = 0;
    virtual void unpack() = 0;

    void run()
    {
    	t1 = std::thread( &Communicable::pack, this );
    	t2 = std::thread( &Communicable::unpack, this );

    	// t1.detach();
    	// t2.detach();
    }

protected:
    std::thread t1;
    std::thread t2;
};

class Test : public Communicable, public eoF<void>
{
public:
    void pack()
    {
    	// while (1)
    	//     {
    		// std::cout << "pack ";
    		// std::cout.flush();
    		// std::this_thread::sleep_for(std::chrono::milliseconds(10));

		em_mutex.lock();
		if ( !em.empty() )
		    {
			std::vector< int > v = em.front();
			em.pop();
		    }
		em_mutex.unlock();
    	    // }
    }

    void unpack()
    {
    	// while (1)
    	//     {
    		// std::cout << "unpack ";
    		// std::cout.flush();
    		// std::this_thread::sleep_for(std::chrono::milliseconds(10));

		em_mutex.lock();
		std::vector< int > v( 10, 10 );
		em.push(v);
		em_mutex.unlock();
    	    // }
    }

    void operator()() {}

private:
    std::queue< std::vector< int > > em;
    std::mutex em_mutex;
};

int main(void)
{
    Test comm;

    comm.run();

    for ( int i = 0; i < 100; ++i )
	{
	    // std::cout << "main ";
	    // std::cout.flush();
	    std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}

    return 0;
}
