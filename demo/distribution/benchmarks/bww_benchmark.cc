#include "bww_benchmark.h"

using namespace bdm;

zmqpp::context Info::ctx {};
std::string Info::worker1 ("W1");
std::string Info::worker2 ("W2");
size_t Info::n_messages (1000);
bool Info::verbose (false);

int main(int argc, char *argv[]) {
    Info::n_messages = (argc >= 2 ? atoi(argv[1]) : Info::n_messages);
    Info::verbose = (argc == 2 && strcmp(argv[1], "-v") == 0) ||
                    (argc == 3 && strcmp(argv[2], "-v") == 0);

    std::thread broker_t (BrokerTask);
    std::thread worker1_t (Worker1Task);
    std::thread worker2_t (Worker2Task);

    std::this_thread::sleep_for( std::chrono::seconds(1) );
    std::thread client_t (ClientTask);

    client_t.join();
    worker1_t.join();
    worker2_t.join();


    return 0;
}
