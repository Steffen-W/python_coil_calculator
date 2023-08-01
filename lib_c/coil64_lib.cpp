// g++ -shared -o coil64_lib.so -fPIC -I /usr/include/python3.10 -lpython3.10 *.cpp

#include <Python.h>

#include "bessel.h"
#include "resolv.h"
#include "resolves.h"
#include "resolve_q.h"
#include "resolve_srf_cs.h"

#include <map>
#include <cstdio>
using namespace std;

static PyObject *convert2Dict(std::map<std::string, double> result);

// Python-Wrapper für die C++ Funktion
static PyObject *py_getOneLayerN_byWindingLength(PyObject *self, PyObject *args)
{
    double D, L, I;
    unsigned int accuracy = 1;
    if (!PyArg_ParseTuple(args, "ddd", &D, &L, &I))
        return NULL;

    _CoilResult result;
    double res = getOneLayerN_byWindingLength(D, L, I, &result, accuracy);

    std::map<std::string, double>
        data = {
            {"lw", result.sec},
            {"dw", result.five},
            {"N", res}};

    return convert2Dict(data);
}

static PyObject *py_getOneLayerN_Poligonal(PyObject *self, PyObject *args)
{
    double I, Dk, dw, p, n;
    unsigned int accuracy = 1;
    if (!PyArg_ParseTuple(args, "ddddd", &I, &Dk, &dw, &p, &n))
        return NULL;

    _CoilResult result;
    double res = getOneLayerN_Poligonal(I, Dk, dw, p, n, &result, accuracy);

    std::map<std::string, double> data = {
        {"p * N", result.sec},
        {"lw", result.thd},
        {"iDk", result.seven},
        {"N", res}};

    return convert2Dict(data);
}

static PyObject *py_odCalc(PyObject *self, PyObject *args)
{
    double id;
    if (!PyArg_ParseTuple(args, "d", &id))
        return NULL;

    double result = odCalc(id);

    return PyFloat_FromDouble(result);
}

static PyObject *py_find_Cs(PyObject *self, PyObject *args)
{
    double p, Dk, lk;
    if (!PyArg_ParseTuple(args, "ddd", &p, &Dk, &lk))
        return NULL;

    double result = find_Cs(p, Dk, lk);

    return PyFloat_FromDouble(result);
}

static PyObject *py_solve_Qr(PyObject *self, PyObject *args)
{
    double I, Df, pm, dw, fa, N, Cs;
    int mt;
    if (!PyArg_ParseTuple(args, "dddddddi", &I, &Df, &pm, &dw, &fa, &N, &Cs, &mt))
        return NULL;

    _CoilResult result;
    double res = solve_Qr(I, Df, pm, dw, fa, N, Cs, (Material)mt, &result);

    std::map<std::string, double> data = {
        {"Rac", result.seven},
        {"R_ind / Rac", res}};

    return convert2Dict(data);
}

static PyObject *py_findSRF(PyObject *self, PyObject *args)
{
    double lk, Dk, lw;
    if (!PyArg_ParseTuple(args, "ddd", &lk, &Dk, &lw))
        return NULL;

    double result = findSRF(lk, Dk, lw);

    return PyFloat_FromDouble(result);
}

static PyObject *py_getMultiLayerN(PyObject *self, PyObject *args)
{
    double I, D, dw, k, lk, g, Ng;
    if (!PyArg_ParseTuple(args, "ddddddd", &I, &D, &dw, &k, &lk, &g, &Ng))
        return NULL;

    _CoilResult result;
    getMultiLayerN(I, D, dw, k, lk, g, Ng, &result);

    std::map<std::string, double> data = {
        {"Number turns", result.six},
        {"Thickness", result.fourth},
        {"Length", result.sec},
        {"R_DC", result.N},
        {"Ng", result.five},
        {"Number layers", result.thd}};

    // print("Number of interlayers Ng = {}".format(result.five))

    return convert2Dict(data);
}

static PyObject *py_getMultiLayerN_rectFormer(PyObject *self, PyObject *args)
{
    double Ind, a, b, l, dw, k;
    if (!PyArg_ParseTuple(args, "dddddd", &Ind, &a, &b, &l, &dw, &k))
        return NULL;

    _CoilResult result;
    getMultiLayerN_rectFormer(Ind, a, b, l, dw, k, &result);

    std::map<std::string, double> data = {
        {"Number turns", result.N},
        {"Number layers", result.sec},
        {"Length wire", result.thd},
        {"Rdc", result.fourth},
        {"thickness", result.five}};

    return convert2Dict(data);
}

static PyObject *py_getMultilayerN_Foil(PyObject *self, PyObject *args)
{
    double D, w, t, ins, I;
    if (!PyArg_ParseTuple(args, "ddddd", &D, &w, &t, &ins, &I))
        return NULL;

    _CoilResult result;
    getMultilayerN_Foil(D, w, t, ins, I, &result);

    std::map<std::string, double> data = {
        {"Number turns", result.N},
        {"Length spiral", result.sec},
        {"Rdcc", result.fourth},
        {"Rdca", result.five},
        {"Do", result.thd}};

    return convert2Dict(data);
}

static PyObject *py_getFerriteN(PyObject *self, PyObject *args)
{
    double L, Do, Di, h, dw, mu, Ch;
    if (!PyArg_ParseTuple(args, "ddddddd", &L, &Do, &Di, &h, &dw, &mu, &Ch))
        return NULL;

    _CoilResult result;
    getFerriteN(L, Do, Di, h, dw, mu, Ch, &result);

    std::map<std::string, double> data = {
        {"Number turns", result.N},
        {"Length wire", result.sec},
        {"Al", result.thd}};

    return convert2Dict(data);
}

static PyObject *py_calc_getPCB_N(PyObject *self, PyObject *args)
{
    double I, D, d, ratio;
    int layout;
    if (!PyArg_ParseTuple(args, "ddddi", &I, &D, &d, &ratio, &layout))
        return NULL;

    _CoilResult result;
    getPCB_N(I, D, d, ratio, layout, &result);

    std::map<std::string, double> data = {
        {"Number turns", result.N},
        {"Winding pitch", result.sec},
        {"Width", result.thd}};

    return convert2Dict(data);
}

static PyObject *py_getPCB_RectN(PyObject *self, PyObject *args)
{
    double I, A, B, _a, th, ratio;
    if (!PyArg_ParseTuple(args, "dddddd", &I, &A, &B, &_a, &th, &ratio))
        return NULL;

    _CoilResult result;
    getPCB_RectN(I, A, B, _a, th, ratio, &result);

    std::map<std::string, double> data = {
        {"Number turns", result.N},
        {"Winding pitch", result.sec},
        {"Width", result.thd}};

    return convert2Dict(data);
}

static PyObject *py_solve_Qpcb(PyObject *self, PyObject *args)
{
    long N;
    double _I, _D, _d, _W, _t, _s, _f;
    int layout;
    if (!PyArg_ParseTuple(args, "ldddddddi", &N, &_I, &_D, &_d, &_W, &_t, &_s, &_f, &layout))
        return NULL;

    double result = solve_Qpcb(N, _I, _D, _d, _W, _t, _s, _f, layout);
    return PyFloat_FromDouble(result);
}

static PyObject *py_getSpiralN(PyObject *self, PyObject *args)
{
    double I, Di, dw, s;
    if (!PyArg_ParseTuple(args, "dddd", &I, &Di, &dw, &s))
        return NULL;

    _CoilResult result;
    getSpiralN(I, Di, dw, s, &result);

    std::map<std::string, double> data = {
        {"Number turns", result.N},
        {"Length spiral", result.sec},
        {"Do", result.thd}};
    return convert2Dict(data);
}

// Methode-Definitionen für das Python-Modul
static PyMethodDef methods[] = {
    {"calc_getPCB_N", py_calc_getPCB_N, METH_VARARGS, "calc_getPCB_N"},
    {"solve_Qpcb", py_solve_Qpcb, METH_VARARGS, "solve_Qpcb"},
    {"getSpiralN", py_getSpiralN, METH_VARARGS, "getSpiralN"},
    {"getPCB_RectN", py_getPCB_RectN, METH_VARARGS, "getPCB_RectN"},
    {"getFerriteN", py_getFerriteN, METH_VARARGS, "getFerriteN"},
    {"getMultiLayerN", py_getMultiLayerN, METH_VARARGS, "getMultiLayerN"},
    {"getMultilayerN_Foil", py_getMultilayerN_Foil, METH_VARARGS, "getMultilayerN_Foil"},
    {"getMultiLayerN_rectFormer", py_getMultiLayerN_rectFormer, METH_VARARGS, "getMultiLayerN_rectFormer"},
    {"getOneLayerN_byWindingLength", py_getOneLayerN_byWindingLength, METH_VARARGS, "getOneLayerN_byWindingLength"},
    {"getOneLayerN_Poligonal", py_getOneLayerN_Poligonal, METH_VARARGS, "getOneLayerN_Poligonal"},
    {"odCalc", py_odCalc, METH_VARARGS, "odCalc"},
    {"find_Cs", py_find_Cs, METH_VARARGS, "find_Cs"},
    {"solve_Qr", py_solve_Qr, METH_VARARGS, "solve_Qr"},
    {"findSRF", py_findSRF, METH_VARARGS, "findSRF"},
    {NULL, NULL, 0, NULL} // Sentinel-Wert am Ende der Methodenliste
};

// Modul-Definition
static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "coil64_lib",
    NULL,
    -1,
    methods};

// Modul-Initialisierungsfunktion
PyMODINIT_FUNC PyInit_coil64_lib()
{
    return PyModule_Create(&module);
}

static PyObject *convert2Dict(std::map<std::string, double> result)
{
    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    // Füge die Schlüssel-Wert-Paare aus dem C++-Dict in das Python-Dict ein
    for (const auto &entry : result)
    {
        PyObject *key = PyUnicode_FromString(entry.first.c_str());
        PyObject *value = PyFloat_FromDouble(entry.second);
        PyDict_SetItem(py_result, key, value);
        Py_DECREF(key);
        Py_DECREF(value);
    }

    return py_result;
}