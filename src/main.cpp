#include <iostream>

#include "LDPC/LDPC.hpp"
// #include "alist/Alist.hpp"
// #include <Eigen/Eigen>

int main(int, char**) {
    // LDPC ldpc("73.alist");
    // LDPC ldpc("example/CCSDS_ldpc_n128_k64.alist");
    // std::cout << ldpc.getK() << std::endl;
    Alist<alist_matrix> a("example/73.alist");
    Eigen::SparseMatrix<int> mat = a.getMat();
    std::cout << mat << std::endl;
}

// int main() {
//     double a;
//     Eigen::Vector3i index1(11, 21, 31);

//     a = index1.norm();

//     std::cout << "a is " << a << std::endl;

//     return 0;
// }