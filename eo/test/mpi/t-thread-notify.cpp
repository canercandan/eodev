#include <iostream>
#include <condition_variable>
#include <thread>
#include <chrono>

using namespace std;

condition_variable cv;
mutex cv_m;
mutex m;
bool ready = false;

void waits()
{
    while (1)
	{
	    unique_lock<mutex> lk(cv_m);
	    cerr << "Waiting... \n";
	    cv.wait(lk);
	    cerr << "...finished waiting.";
	}
}

void signals()
{
    unique_lock<mutex> lk(m);
    while (1)
	{
	    this_thread::sleep_for(chrono::seconds(1));
	    cerr << "Notifying...\n";
	    cv.notify_all();
	}
    ready = true;
    std::notify_all_at_thread_exit(cv, std::move(lk));
}

int main()
{
    thread t1(waits), t2(waits), t3(waits), t4(signals);
    t1.join(); t2.join(), t3.join(), t4.join();
}
