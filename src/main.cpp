#include <iostream>

#include "LDPC/LDPC.hpp"
// #include "alist/Alist.hpp"
// #include <Eigen/Eigen>

int main(int, char**) {
    LDPC ldpc("example/simon.alist");

    // LDPC ldpc("example/CCSDS_ldpc_n128_k64.alist");

    Eigen::SparseMatrix<int> H = ldpc.getH(), G = ldpc.getG();
    std::cout << G << H << std::endl;
    std::cout << G * H.transpose() << std::endl;
    // std::cout << mat << std::endl;
}

// int main() {
//     double a;
//     Eigen::Vector3i index1(11, 21, 31);

//     a = index1.norm();

//     std::cout << "a is " << a << std::endl;

//     return 0;
// }