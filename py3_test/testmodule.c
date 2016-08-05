#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

static PyObject* foobar(PyObject *self, PyObject * args);


struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyMethodDef test23_methods[] = {
    {"foobar",  foobar, METH_VARARGS, "Test function"},
    {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int test23_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int test23_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "test23",
        NULL,
        sizeof(struct module_state),
        test23_methods,
        NULL,
        test23_traverse,
        test23_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_test23(void)

#else
#define INITERROR return

void
inittest23(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("test23", test23_methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("test23.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    import_array();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

static PyObject*
foobar(PyObject *self, PyObject * args)
{
  return Py_BuildValue("");
}
