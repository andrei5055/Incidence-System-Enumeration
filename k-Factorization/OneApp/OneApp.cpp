#include <sycl/sycl.hpp>
#include "OneApp.h"

#include <string>

#if defined(FPGA_HARDWARE) || defined(FPGA_EMULATOR) || defined(FPGA_SIMULATOR)
#include <sycl/ext/intel/fpga_extensions.hpp>
#endif

using namespace sycl;

//************************************
// Define a bitwise matrix multiplication kernel using an nd_range parallel_for
void binary_matmul(sycl::queue& q,
    const unsigned int* A,
    const unsigned int* B,
    unsigned int* C,
    size_t M, size_t N, size_t K) {
    q.submit([&](sycl::handler& h) {
        // Define a 2D nd_range for matrix C, with one work-item per output element
        sycl::range<2> global_range(M, N);
        sycl::range<2> local_range(16, 16); // Tune local work group size for performance

        h.parallel_for(sycl::nd_range<2>(global_range, local_range), [=](sycl::nd_item<2> item) {
            size_t row = item.get_global_id(0);
            size_t col = item.get_global_id(1);

            unsigned int sum = 0;
            for (size_t i = 0; i < K; ++i) {
                // Perform the binary AND operation and accumulate
                // Use a standard binary multiplication logic
                sum += A[row * K + i] & B[i * N + col];
            }
            // A final binary operation (e.g., population count, check for non-zero) might be needed
            // depending on the exact semantic. Here we sum the bits for simplicity.
            C[row * N + col] = sum;
            });
        });
    q.wait();
}

void MultArrays(const long long* a, const long long* b,
    long long* res, int num, int nc) {

    // Create the range object for the vectors managed by the buffer.
    range<1> num_items{ num };
    // Create buffers that hold the data shared between the host and the devices.
    // The buffer destructor is responsible to copy the data back to host when it
    // goes out of scope.
    static queue* Queue = NULL;
    if (!Queue) {
        // SYCL_UR_TRACE=2
        // Create a SYCL queue targeting a GPU device
        Queue = new queue(gpu_selector_v);
        std::cout << "Running on device: " << Queue->get_device().get_info<sycl::info::device::name>() << std::endl;
    }

    buffer<long long> bufA(a, num_items);
    buffer<long long> bufB(b, num_items);
    buffer<long long> bufRes(res, num_items);

    // Submit a command group to the queue by a lambda function that contains the
    // data access permission and device computation (kernel).
    auto tm = clock();
    range<1> padded_range{ ((num_items[0] + 31) / 32) * 32 };

    Queue->submit([&](handler& h) {
        // Create an accessor for each buffer with access permission: read, write or
        // read/write. The accessor is a mean to access the memory in the buffer.
        accessor A(bufA, h, read_only);
        accessor B(bufB, h, read_only);

        // The sum_accessor is used to store (with write permission) the sum data.
        accessor AandB(bufRes, h, write_only, no_init);

        // Use parallel_for to run vector addition in parallel on device. This
        // executes the kernel.
        //    1st parameter is the number of work items.
        //    2nd parameter is the kernel, a lambda that specifies what to do per
        //    work item. The parameter of the lambda is the work item id.
        // SYCL supports unnamed lambda kernel by default.
#define USE_ChatGPT_code 0
#if !USE_ChatGPT_code
        // That one is faster 561 vs 753 ms for nc = 1
        //          and much more faster for nc = 5000000
        h.parallel_for<class VectorAdd>(padded_range, [=](id<1> i) {
            for (int k = 0; k < nc; k++)
                AandB[i] = A[i] & B[i];
        });
#else
        constexpr int VL = 8;  // vector length (matches AVX2 256-bit for int32)
        h.parallel_for<class VectorAdd>(
            sycl::range<1>(padded_range / VL), [=](sycl::id<1> idx) {
            for (int k = 0; k < nc; k++) {
                size_t base = idx[0] * VL;

                // Load 8 elements at once into a vector
                sycl::vec<int, VL> va;
                sycl::vec<int, VL> vb;
                for (int j = 0; j < VL; ++j) {
                    va[j] = A[base + j];
                    vb[j] = B[base + j];
                }

                // Portable vectorized AND
                sycl::vec<int, VL> vr = va & vb;

                // Store back results
                for (int j = 0; j < VL; ++j) {
                    AandB[base + j] = vr[j];
                }
            }
        });
#endif
    });
    // Wait until compute tasks on GPU done
    Queue->wait();
    printf("gpu x %d iterations  lenArray = %d: %ld ms\n", nc, num, clock() - tm);
}

//************************************
// Initialize the vector from 0 to vector_size - 1
//************************************
void InitializeVector(IntVector &a) {
    for (size_t i = 0; i < a.size(); i++) a.at(i) = i;
}
