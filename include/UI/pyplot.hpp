/**
 * @file pyplot.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-12-13 22:16:26
 * @modified: 2021-04-27 13:10:44
 */

#ifndef PYPLOT_HPP
#define PYPLOT_HPP

#define PY_SSIZE_T_CLEAN
#include <Python.h>

int pyplot(float* x, float* y, size_t length, const char* label);

#endif
