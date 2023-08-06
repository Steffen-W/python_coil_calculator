// g++ -shared -o coil64_lib.so -fPIC -I /usr/include/python3.10 -lpython3.10 *.cpp

#include <Python.h>
#include <stdbool.h>

// #include "bessel.h"
// // #include "resolv.h"
// #include "resolves.h"
// #include "resolve_q.h"
// #include "resolve_srf_cs.h"

struct _CoilResult
{
    double N;
    double sec;
    double thd;
    double fourth;
    double five;
    unsigned long int six;
    double seven;
};
typedef struct _CoilResult _CoilResult;

extern double getOneLayerI_withRoundWire(double Dk, double dw, double p, double N, double *lw, unsigned int accuracy);
extern void getMultilayerI_Foil(double D, double w, double t, double ins, int _N, _CoilResult *result);
extern double getFerriteI(double N, double Do, double Di, double h, double mu, double Ch, _CoilResult *result);
extern double getPCB_I(double N, double _d, double _s, int layout, _CoilResult *result);
extern double getPCB_RectI(int N, double A, double B, double s, double w, double th, _CoilResult *result);
extern void getSpiralI(double Do, double Di, double dw, int _N, _CoilResult *result);
extern double getOneLayerN_byWindingLength(double D, double L, double I, _CoilResult *result, unsigned int accuracy);
extern double getOneLayerN_Poligonal(double I, double Dk, double dw, double p, double n, _CoilResult *result, unsigned int accuracy);
extern double getOneLayerN_withRectWire(double Dk, double w, double t, double p, double I, double *lw, unsigned int accuracy);
extern unsigned long int solve_Qc(double I, double Df, double pm, double _w, double _t, double fa, double N, double Cs, int mt, _CoilResult *result);
extern double odCalc(double id);
extern double find_Cs(double p, double Dk, double lk);
extern unsigned long int solve_Qr(double I, double Df, double pm, double dw, double fa, double N, double Cs, int mt, _CoilResult *result);
extern double findSRF(double lk, double Dk, double lw);
extern void getMultiLayerN(double I, double D, double dw, double k, double lk, double gap, long Ng, _CoilResult *result);
extern void getMultiLayerN_rectFormer(double Ind, double a, double b, double l, double dw, double k, _CoilResult *result);
extern void getMultilayerN_Foil(double D, double w, double t, double ins, double I, _CoilResult *result);
extern void getFerriteN(double L, double Do, double Di, double h, double dw, double mu, double Ch, _CoilResult *result);
extern void getPCB_RectN(double I, double A, double B, double _a, double th, double ratio, _CoilResult *result);
extern void getPCB_N(double I, double D, double d, double ratio, int layout, _CoilResult *result);
extern double solve_Qpcb(long N, double _I, double _D, double _d, double _W, double _t, double _s, double _f, int layout);
extern void getSpiralN(double I, double Di, double dw, double s, _CoilResult *result);

void add_dict_item(PyObject *dict, const char *key, double value)
{
    PyObject *py_key = PyUnicode_FromString(key);
    PyObject *py_value = PyFloat_FromDouble(value);
    PyDict_SetItem(dict, py_key, py_value);
    Py_DECREF(py_key);
    Py_DECREF(py_value);
}

// Python-Wrapper für die C++ Funktion

static PyObject *py_getOneLayerI_withRoundWire(PyObject *self, PyObject *args)
{
    // getOneLayerI_withRoundWire(D, d, p, N, accuracy)
    double Dk, dw, p, N;
    unsigned int accuracy = 1;
    if (!PyArg_ParseTuple(args, "dddd", &Dk, &dw, &p, &N))
        return NULL;

    double lw;
    double res = getOneLayerI_withRoundWire(Dk, dw, p, N, &lw, accuracy);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "lw", lw);
    add_dict_item(py_result, "L", res);

    return py_result;
}

static PyObject *py_getMultilayerI_Foil(PyObject *self, PyObject *args)
{
    // getMultilayerI_Foil(D, w, t, ins, N, result)
    double D, w, t, ins;
    int _N;
    if (!PyArg_ParseTuple(args, "ddddi", &D, &w, &t, &ins, &_N))
        return NULL;

    _CoilResult result;
    getMultilayerI_Foil(D, w, t, ins, _N, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "L", result.N);
    add_dict_item(py_result, "Length", result.sec);
    add_dict_item(py_result, "Do", result.thd);
    add_dict_item(py_result, "R_DC", result.fourth);
    add_dict_item(py_result, "R_AC", result.five);

    return py_result;
}

static PyObject *py_getFerriteI(PyObject *self, PyObject *args)
{
    // getFerriteI(N, OD, ID, h, mu, C, result)
    double N, Do, Di, h, mu, Ch;
    if (!PyArg_ParseTuple(args, "dddddd", &N, &Do, &Di, &h, &mu, &Ch))
        return NULL;

    _CoilResult result;
    double res = getFerriteI(N, Do, Di, h, mu, Ch, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "L", res);
    add_dict_item(py_result, "Al", result.thd);

    return py_result;
}

static PyObject *py_getPCB_I(PyObject *self, PyObject *args)
{
    // getPCB_I(N, d, s, layoutPCB.value, result)
    double N, _d, _s;
    int layout;
    if (!PyArg_ParseTuple(args, "dddi", &N, &_d, &_s, &layout))
        return NULL;

    _CoilResult result;
    double res = getPCB_I(N, _d, _s, layout, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Do", result.five);
    add_dict_item(py_result, "L", res);

    return py_result;
}

static PyObject *py_getPCB_RectI(PyObject *self, PyObject *args)
{
    // getPCB_RectI(N, A, B, s, W, t, result)
    int N;
    double A, B, s, w, th;
    if (!PyArg_ParseTuple(args, "iddddd", &N, &A, &B, &s, &w, &th))
        return NULL;

    _CoilResult result;
    double res = getPCB_RectI(N, A, B, s, w, th, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Do", result.five);
    add_dict_item(py_result, "L", res);

    return py_result;
}

static PyObject *py_getSpiralI(PyObject *self, PyObject *args)
{
    // getSpiralI(OD, ID, d, N, result)
    double Do, Di, dw;
    int _N;
    if (!PyArg_ParseTuple(args, "dddi", &Do, &Di, &dw, &_N))
        return NULL;

    _CoilResult result;
    getSpiralI(Do, Di, dw, _N, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Number turns", result.N);
    add_dict_item(py_result, "Length spiral", result.sec);

    return py_result;
}

static PyObject *py_getOneLayerN_byWindingLength(PyObject *self, PyObject *args)
{
    double D, L, I;
    unsigned int accuracy = 1;
    if (!PyArg_ParseTuple(args, "ddd", &D, &L, &I))
        return NULL;

    _CoilResult result;
    double res = getOneLayerN_byWindingLength(D, L, I, &result, accuracy);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "lw", result.sec);
    add_dict_item(py_result, "dw", result.five);
    add_dict_item(py_result, "N", res);

    return py_result;
}

static PyObject *py_getOneLayerN_Poligonal(PyObject *self, PyObject *args)
{
    double I, Dk, dw, p, n;
    unsigned int accuracy = 1;
    if (!PyArg_ParseTuple(args, "ddddd", &I, &Dk, &dw, &p, &n))
        return NULL;

    _CoilResult result;
    double res = getOneLayerN_Poligonal(I, Dk, dw, p, n, &result, accuracy);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "p * N", result.sec);
    add_dict_item(py_result, "lw", result.thd);
    add_dict_item(py_result, "iDk", result.seven);
    add_dict_item(py_result, "N", res);

    return py_result;
}

static PyObject *py_getOneLayerN_withRectWire(PyObject *self, PyObject *args)
{
    double Dk, w, t, p, I;
    if (!PyArg_ParseTuple(args, "ddddd", &Dk, &w, &t, &p, &I))
        return NULL;

    unsigned int accuracy = 1;
    double lw;
    double res = getOneLayerN_withRectWire(Dk, w, t, p, I, &lw, accuracy);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "lw", lw);
    add_dict_item(py_result, "N", res);

    return py_result;
}

static PyObject *py_solve_Qc(PyObject *self, PyObject *args)
{
    double I, Df, pm, _w, _t, fa, N, Cs;
    int mt;
    if (!PyArg_ParseTuple(args, "ddddddddi", &I, &Df, &pm, &_w, &_t, &fa, &N, &Cs, &mt))
        return NULL;

    _CoilResult result;
    unsigned long int res = solve_Qc(I, Df, pm, _w, _t, fa, N, Cs, mt, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Rac", result.seven);
    add_dict_item(py_result, "Q", res);

    return py_result;
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
    double res = solve_Qr(I, Df, pm, dw, fa, N, Cs, mt, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Rac", result.seven);
    add_dict_item(py_result, "R_ind / Rac", res);

    return py_result;
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

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Number turns", result.six);
    add_dict_item(py_result, "Thickness", result.fourth);
    add_dict_item(py_result, "Length", result.sec);
    add_dict_item(py_result, "R_DC", result.N);
    add_dict_item(py_result, "Ng", result.five);
    add_dict_item(py_result, "Number layers", result.thd);

    return py_result;
}

static PyObject *py_getMultiLayerN_rectFormer(PyObject *self, PyObject *args)
{
    double Ind, a, b, l, dw, k;
    if (!PyArg_ParseTuple(args, "dddddd", &Ind, &a, &b, &l, &dw, &k))
        return NULL;

    _CoilResult result;
    getMultiLayerN_rectFormer(Ind, a, b, l, dw, k, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Number turns", result.N);
    add_dict_item(py_result, "Number layers", result.sec);
    add_dict_item(py_result, "Length wire", result.thd);
    add_dict_item(py_result, "Rdc", result.fourth);
    add_dict_item(py_result, "thickness", result.five);

    return py_result;
}

static PyObject *py_getMultilayerN_Foil(PyObject *self, PyObject *args)
{
    double D, w, t, ins, I;
    if (!PyArg_ParseTuple(args, "ddddd", &D, &w, &t, &ins, &I))
        return NULL;

    _CoilResult result;
    getMultilayerN_Foil(D, w, t, ins, I, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Number turns", result.N);
    add_dict_item(py_result, "Length spiral", result.sec);
    add_dict_item(py_result, "Rdcc", result.fourth);
    add_dict_item(py_result, "Rdca", result.five);
    add_dict_item(py_result, "Do", result.thd);

    return py_result;
}

static PyObject *py_getFerriteN(PyObject *self, PyObject *args)
{
    double L, Do, Di, h, dw, mu, Ch;
    if (!PyArg_ParseTuple(args, "ddddddd", &L, &Do, &Di, &h, &dw, &mu, &Ch))
        return NULL;

    _CoilResult result;
    getFerriteN(L, Do, Di, h, dw, mu, Ch, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Number turns", result.N);
    add_dict_item(py_result, "Length wire", result.sec);
    add_dict_item(py_result, "Al", result.thd);

    return py_result;
}

static PyObject *py_calc_getPCB_N(PyObject *self, PyObject *args)
{
    double I, D, d, ratio;
    int layout;
    if (!PyArg_ParseTuple(args, "ddddi", &I, &D, &d, &ratio, &layout))
        return NULL;

    _CoilResult result;
    getPCB_N(I, D, d, ratio, layout, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Number turns", result.N);
    add_dict_item(py_result, "Winding pitch", result.sec);
    add_dict_item(py_result, "Width", result.thd);

    return py_result;
}

static PyObject *py_getPCB_RectN(PyObject *self, PyObject *args)
{
    double I, A, B, _a, th, ratio;
    if (!PyArg_ParseTuple(args, "dddddd", &I, &A, &B, &_a, &th, &ratio))
        return NULL;

    _CoilResult result;
    getPCB_RectN(I, A, B, _a, th, ratio, &result);

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Number turns", result.N);
    add_dict_item(py_result, "Winding pitch", result.sec);
    add_dict_item(py_result, "Width", result.thd);

    return py_result;
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

    // Erstelle das Python-Dict
    PyObject *py_result = PyDict_New();
    if (py_result == NULL)
    {
        return NULL;
    }

    add_dict_item(py_result, "Number turns", result.N);
    add_dict_item(py_result, "Length spiral", result.sec);
    add_dict_item(py_result, "Do", result.thd);

    return py_result;
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
    {"solve_Qc", py_solve_Qc, METH_VARARGS, "solve_Qc"},
    {"getOneLayerN_withRectWire", py_getOneLayerN_withRectWire, METH_VARARGS, "getOneLayerN_withRectWire"},
    {"find_Cs", py_find_Cs, METH_VARARGS, "find_Cs"},
    {"solve_Qr", py_solve_Qr, METH_VARARGS, "solve_Qr"},
    {"findSRF", py_findSRF, METH_VARARGS, "findSRF"},
    {"getOneLayerI_withRoundWire", py_getOneLayerI_withRoundWire, METH_VARARGS, "getOneLayerI_withRoundWire"},
    {"getMultilayerI_Foil", py_getMultilayerI_Foil, METH_VARARGS, "getMultilayerI_Foil"},
    {"getFerriteI", py_getFerriteI, METH_VARARGS, "getFerriteI"},
    {"getPCB_I", py_getPCB_I, METH_VARARGS, "getPCB_I"},
    {"getPCB_RectI", py_getPCB_RectI, METH_VARARGS, "getPCB_RectI"},
    {"getSpiralI", py_getSpiralI, METH_VARARGS, "getSpiralI"},
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