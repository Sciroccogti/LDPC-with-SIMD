#include <iostream>

#include "LDPC/LDPC.hpp"
#include "alist/Alist.hpp"

int main(int, char**) {
    // LDPC ldpc("73.alist");
    LDPC ldpc("example/CCSDS_ldpc_n128_k64.alist");
    std::cout << ldpc.getK() << std::endl;
}
