# LDPC-with-SIMD

LDPC with SIMD.

## Environment

- OS: Ubuntu 20.04.1 or 18.04.4
- Instruction: AVX2 or AVX512
- Dependency:
  - [CMake](https://cmake.org/)
  - [Conon](https://conan.io/)
  - [Eigen3](http://eigen.tuxfamily.org)
  - [Python3-dev](https://www.python.org/)

## Functions

### Binary LDPC

- [x] Encode
- [x] Modulate(BPSK)
- [x] Channel(AWGN)
- [x] Demodulate(BPSK)
- [x] Decode(NMS/SPA)

### NonBinary LDPC

- [x] Encode
- [x] Modulate(BPSK)
- [x] Channel(AWGN)
- [x] Demodulate(BPSK)
- [x] Decode(EMS)

### Tools

- [x] Draw Tanner
- [x] Type in Alist

### Install Conan

```Bash
sudo apt install python3-pip
pip3 install conan
```

### Install Python3-dev

```Bash
sudo apt install python3-dev
```

## How to use

### 1 Simulation

```Bash
mkdir build
cd build
conan install ..
conan profile update settings.compiler.libcxx=libstdc++11 default
cmake ..
cmake --build .
../bin/LDPC-with-SIMD --dec-h-path ../example/H.alist
```

### 2 Tools

#### 2-1 Draw Tanner

```Bash
python3 script/tanner.py --dec-h-path example/H.alist
```

![](assets/tanner.png)

#### 2-1 Type in Alist

```Bash
python3 script/makeAlist.py
# nRow is: 4
# nCol is: 6
# 1 1 0 1 0 0 
# 0 1 1 0 1 0
# 1 0 0 0 1 1
# 0 0 1 1 0 1
# filename is: aaa.alist
```

## References

### Binary LDPC and Basic

- [pyldpc](https://github.com/hichamjanati/pyldpc.git)
- [5G-SIMD-LDPC](https://github.com/SherlockHsu/5G-SIMD-LDPC)
- [Python处理alist文件](https://www.cnblogs.com/lingr7/p/13038410.html)
- [aff3ct](https://github.com/aff3ct/aff3ct)
- [paper: NMS](https://www.researchgate.net/publication/3160637_Near_optimum_universal_belief_propagation_based_decoding_of_low-density_parity_check_codes)

### Nonbinary LDPC

- [NB_LDPC_FB](https://github.com/cedricomarchando/NB_LDPC_FB)
- [Kaiserslautern database](https://www.uni-kl.de/channel-codes/channel-codes-database/non-binary-ldpc/)

### C++ tutorial

- [c++模板类(一)理解编译器的编译模板过程](http://blog.csdn.net/onafioo/article/details/29857281)
- [Conan-CMake transparent integration](https://blog.conan.io/2018/06/11/Transparent-CMake-Integration.html)
