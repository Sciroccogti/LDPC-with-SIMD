/**
 * @file pyplot.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-12-13 22:16:26
 * @modified: 2020-12-13 22:16:36
 */

#ifndef PYPLOT_HPP
#define PYPLOT_HPP

#define PY_SSIZE_T_CLEAN
#include <Python.h>

int pyplot(double* x, double* y, size_t length);

#endif
