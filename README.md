# LPDC-with-SIMD

LPDC with SIMD aimed at AMD CPUs.

## Environment

- OS: Ubuntu 20.04.1
- CPU: AMD 3700X

## How to use

```Bash
mkdir build
cd build
cmake ..
make
./LDPC-with-SIMD --dec-h-path ../example/H.alist
```

## References

### Binary LDPC and Basic

- [pyldpc](https://github.com/hichamjanati/pyldpc.git)
- [5G-SIMD-LDPC](https://github.com/SherlockHsu/5G-SIMD-LDPC)
- [Python处理alist文件](https://www.cnblogs.com/lingr7/p/13038410.html) 及其配套代码 [LDPC_alist](https://github.com/lingr7/LDPC_alist)

### Nonbinary LDPC

- [NB_LDPC_FB](https://github.com/cedricomarchando/NB_LDPC_FB)
