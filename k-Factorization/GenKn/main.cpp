#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif
#include "generateKn.h"

// Helper to get elapsed time in [HH:MM:SS] format
static std::string get_timestamp() {
    static const auto start_time = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    auto elapsed_secs = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
    long long hours = elapsed_secs / 3600;
    long long minutes = (elapsed_secs % 3600) / 60;
    long long seconds = elapsed_secs % 60;
    std::ostringstream oss;
    oss << "[" << std::setfill('0')
        << std::setw(2) << hours << ":"
        << std::setw(2) << minutes << ":"
        << std::setw(2) << seconds << "] ";
    return oss.str();
}

int main(int argc, char* argv[]) {
#if defined(_WIN32) || defined(_WIN64)
    SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
    // Initialize start time for elapsed tracking
    get_timestamp();

    auto print_usage = []() {
        std::cerr << "Usage: GenKn <N> [p] [k] [m] or keyword arguments (e.g. n=16 nthreads=10)" << std::endl;
        std::cerr << "Available Keyword Parameters:" << std::endl;
        std::cerr << "  n        : Graph size N (must be a positive even integer)" << std::endl;
        std::cerr << "  p        : Symmetry parameter p" << std::endl;
        std::cerr << "  k        : Symmetry parameter k" << std::endl;
        std::cerr << "  m        : Symmetry parameter m" << std::endl;
        std::cerr << "  nthreads : Number of worker threads (default: 8)" << std::endl;
        std::cerr << "  verbose  : Verbosity level (0 = monitor only, 1 = thread details, default: 0)" << std::endl;
        std::cerr << "Examples:" << std::endl;
        std::cerr << "  GenKn n=16 nthreads=12" << std::endl;
        std::cerr << "  GenKn 16 7 2 2 verbose=1" << std::endl;
        std::cerr << "  GenKn help" << std::endl;
    };

    // Check for help flags
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string lower_arg = arg;
        for (char &c : lower_arg) {
            c = std::tolower(static_cast<unsigned char>(c));
        }
        if (lower_arg == "?" || lower_arg == "/?" || lower_arg == "help" || lower_arg == "--help" || lower_arg == "-h") {
            print_usage();
            return 0;
        }
    }

    if (argc < 2) {
        print_usage();
        return 1;
    }

    int N = -1;
    int p = 0;
    int k = 0;
    int m = 0;
    int num_threads = 8;
    int verbose = 0;
    int pos_count = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        size_t eq_pos = arg.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = arg.substr(0, eq_pos);
            std::string val_str = arg.substr(eq_pos + 1);
            for (char &c : key) {
                c = std::tolower(static_cast<unsigned char>(c));
            }
            int val = std::atoi(val_str.c_str());
            if (key == "n") {
                N = val;
            } else if (key == "nthreads") {
                num_threads = val;
            } else if (key == "p") {
                p = val;
            } else if (key == "k") {
                k = val;
            } else if (key == "m") {
                m = val;
            } else if (key == "verbose") {
                verbose = val;
            } else {
                std::cerr << "Error: Unknown parameter '" << key << "'" << std::endl << std::endl;
                print_usage();
                return 1;
            }
        } else {
            bool is_numeric = true;
            for (char c : arg) {
                if (!std::isdigit(static_cast<unsigned char>(c))) {
                    is_numeric = false;
                    break;
                }
            }
            if (!is_numeric) {
                std::cerr << "Error: Unknown parameter '" << arg << "'" << std::endl << std::endl;
                print_usage();
                return 1;
            }
            int val = std::atoi(arg.c_str());
            if (pos_count == 0) {
                N = val;
            } else if (pos_count == 1) {
                p = val;
            } else if (pos_count == 2) {
                k = val;
            } else if (pos_count == 3) {
                m = val;
            }
            pos_count++;
        }
    }

    if (N <= 0 || N % 2 != 0) {
        std::cerr << "Error: N must be a positive even integer." << std::endl;
        return 1;
    }

    if (num_threads <= 0) {
        std::cerr << "Error: nthreads must be a positive integer." << std::endl;
        return 1;
    }

    if (p != 0) {
        std::cout << get_timestamp() << "Starting generateKn for N = " << N 
                  << " (Target combination: p=" << p << ", k=" << k << ", m=" << m 
                  << ") using " << num_threads << " threads" << std::endl;
    } else {
        std::cout << get_timestamp() << "Starting generateKn for N = " << N 
                  << " using " << num_threads << " threads" << std::endl;
    }
    
    int result = generateKn(N, p, k, m, num_threads, verbose);
    return result;
}
