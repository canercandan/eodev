#include <mutex>
#include <iostream>
#include <condition_variable>
#include <thread>
#include <chrono>

using namespace std;

mutex m;
condition_variable cv;

bool ready = false;

void thread_func()
{
    unique_lock<mutex> lk(m);
    ready = true;
    //notify_all_at_thread_exit(cv, move(lk));
    lk.unlock();
    cv.notify_all();
} // destroy thread_locals, notify cv, unlock mutex

int main()
{
    thread t(thread_func);
    t.detach();
    // do other work
    unique_lock<mutex> lk(m);
    while(!ready)
	{
	    cv.wait(lk); // wait for the detached thread
	}
}
