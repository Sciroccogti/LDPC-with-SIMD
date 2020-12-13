#include "UI/pyplot.hpp"

/**
 * @brief use python to plot
 * https://docs.python.org/3.9/extending/embedding.html
 * @param x
 * @param y
 * @param length
 * @return int
 */
int pyplot(double* x, double* y, size_t length) {
    Py_Initialize();
    if (!Py_IsInitialized()) {
        return -1;
    }

    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\"./scripts\")");
    PyRun_SimpleString("import os");
    PyRun_SimpleString("print(os.getcwd())");

    PyObject* pModule = PyImport_ImportModule("waveshow");
    if (pModule == NULL) {
        return -1;
    }

    PyObject* pFunc = PyObject_GetAttrString(pModule, "plot");
    if (!pFunc || !PyCallable_Check(pFunc)) {
        return -1;
    }

    PyObject *x_list = PyList_New(length), *y_list = PyList_New(length);
    for (size_t i = 0; i < length; i++) {
        PyList_SetItem(x_list, i, PyFloat_FromDouble(x[i]));
        PyList_SetItem(y_list, i, PyFloat_FromDouble(y[i]));
    }
    PyObject* label = Py_BuildValue("s", "test");

    PyObject* pArgs = PyTuple_New(3);
    PyTuple_SetItem(pArgs, 0, x_list);
    PyTuple_SetItem(pArgs, 1, y_list);
    PyTuple_SetItem(pArgs, 2, label);

    PyObject_CallObject(pFunc, pArgs);

    Py_Finalize();
    return 0;
}