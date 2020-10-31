#include <iostream>

#include "LDPC/LDPC.hpp"
// #include "alist/Alist.hpp"
// #include <Eigen/Eigen>

int main(int, char**) {
    LDPC ldpc("example/73.alist");

    // LDPC ldpc("example/CCSDS_ldpc_n128_k64.alist");
    std::cout << ldpc.getH() << std::endl;
    Eigen::SparseMatrix<int> mat = ldpc.getG();
    // Eigen::MatrixXf m(3, 4);
    // m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
    // std::cout << m << std::endl;
    // // m.resize(4, 3);
    // std::cout << m.reshaped(4, 3) << std::endl;
    std::cout << mat << std::endl;
}

// int main() {
//     double a;
//     Eigen::Vector3i index1(11, 21, 31);

//     a = index1.norm();

//     std::cout << "a is " << a << std::endl;

//     return 0;
// }