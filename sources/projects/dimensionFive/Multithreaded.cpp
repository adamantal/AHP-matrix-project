#include <vector>
#include <thread>

typedef unsigned long long int Ulli;

std::vector<bool> check;

void excludePermutations(Ulli index) {

}

class ThreadController {
private:
    unsigned int maxThreads;

    std::mutex threadLock;
    std::list<std::thread> threads;

public:
    ThreadController(unsigned int numOfThreads):maxThreads(numOfThreads) {
    }
    void start() {
        if () {

        }
    }
    void threadExited(Ulli index) {
        if (threads.size() < numOfThreads) {
            threadLock.lock();
            std::thread thread(excludePermutations, index);
            threads.
            threadLock.unlock();
        }
    }
}

int main() {
    unsigned int numOfThreads = std::thread::hardware_concurrency();
    std::cout << "The number of threads that can run concurrently is " << numOfThreads << std::endl;


}
