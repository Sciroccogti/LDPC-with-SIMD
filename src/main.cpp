#include <iostream>

#include "alist/Alist.hpp"
// #include "alist/alist_matrix.h"

int main(int, char**) {
    Alist<alist_matrix> alist("example/CCSDS_ldpc_n128_k64.alist");
    std::cout << "ok\n";
    
    alist.save("new");
}
