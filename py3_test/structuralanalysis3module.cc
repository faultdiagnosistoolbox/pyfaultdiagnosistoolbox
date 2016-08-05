extern "C" {
  #include "Python.h"
  #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
  #include <numpy/ndarrayobject.h>
}

#include <iostream>
#include <cstring>
#include "StructuralAnalysisModel.h"
#include "MSOAlg.h"
//#include "timer.h"



int DictToCS( PyObject *csDict, cs **sm );
PyObject* CreateDMpermOutput( csd* dm, cs* sm );

class MSOResultPython : public MSOResult {
public:
  PyObject * CreateOutput( void );
};

PyObject *
MSOResultPython::CreateOutput()
{
  npy_intp n = (npy_intp)this->Size();
  PyArrayObject *out = (PyArrayObject *)PyArray_SimpleNew(1, &n, NPY_OBJECT);
  PyObject *pymso;
  npy_intp idx = 0; 
  for(MSOList::iterator e=msos.begin(); e!= msos.end(); e++ ) {
    npy_intp dims = (npy_intp)e->size();
    pymso = PyArray_SimpleNew(1, &dims, NPY_INT);
    std::copy( e->begin(), e->end(), (int *)PyArray_DATA((PyArrayObject *)pymso) );
    
    PyArray_SETITEM(out, (char *)PyArray_GetPtr(out, &idx), pymso);
    Py_INCREF(pymso);
    idx++;
  }
  
  return (PyObject *)out;
}

static PyObject*
structuralanalysis3_findmsointernal(PyObject *self, PyObject * args)
{
  PyObject *x;

  if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &x)) {
    return NULL;
  }

  Py_INCREF( x ); // Protect sparse matrix dictionary from Python

  // Create CSparse matrix representation
  cs *sma;
  DictToCS( x, &sma );

  // Declare MSO analysis objects
  StructuralAnalysisModel sm(sma);
  MSOAlg msoalg = sm;
  MSOResultPython msos;

  // Run MSO algorithm
  msoalg.MSO( msos );

  // Collect the results
  PyObject *res = msos.CreateOutput();

  Py_DECREF(x); // Release lock on sparse matrix dictionary from Python
  
  return res;
}

static PyObject*
structuralanalysis3_dmperminternal(PyObject *self, PyObject * args)
{
  PyObject *x;

  if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &x)) {
    return NULL;
  }

  Py_INCREF( x );
  cs* a;
  DictToCS( x, &a );
  csd *dm = cs_dmperm(a, 0);
  PyObject* res=CreateDMpermOutput(dm,a);
  
  Py_DECREF(x);  

  return res;
}

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

static PyMethodDef structuralanalysis3_methods[] = {
    {"dmperm_internal",  structuralanalysis3_dmperminternal, METH_VARARGS, "Dulmage-Mendelsohn decomposition"},
    {"findmso_internal",  structuralanalysis3_findmsointernal, METH_VARARGS, "Compute MSO sets"},
    {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int structuralanalysis3_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int structuralanalysis3_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "structuralanalysis3",
        NULL,
        sizeof(struct module_state),
        structuralanalysis3_methods,
        NULL,
        structuralanalysis3_traverse,
        structuralanalysis3_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_structuralanalysis3(void)

#else
#define INITERROR return

void
initstructuralanalysis3(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("structuralanalysis3", structuralanalysis3_methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    char errName[] = "structuralanalysis3.Error";
    st->error = PyErr_NewException(errName, NULL, NULL);
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

int
DictToCS( PyObject *csDict, cs **sm )
{
  // Alloc space for sparse matrix representation
  *sm = (cs *)cs_malloc (1, sizeof(cs));

  // Extract data

  // nzmax - number of non-zero elements
  PyObject *pnzmax = PyDict_GetItemString(csDict, "nzmax");
  //  (*sm)->nzmax = PyInt_AS_LONG(pnzmax);
  (*sm)->nzmax = PyLong_AsLong(pnzmax);
  
  // m - number of rows
  PyObject *pm = PyDict_GetItemString(csDict, "m");
  // (*sm)->m = PyInt_AS_LONG(pm);
  (*sm)->m = PyLong_AsLong(pm);
  
  // n - number of columns
  PyObject *pn = PyDict_GetItemString(csDict, "n");
  // (*sm)->n = PyInt_AS_LONG(pn);
  (*sm)->n = PyLong_AsLong(pn);
  
  // nz - -1 for compressed formate
  (*sm)->nz = -1;
  
  // x - data
  PyArrayObject *px = (PyArrayObject *)PyDict_GetItemString(csDict, "x");
  (*sm)->x = (double *)PyArray_DATA( px );

  // p - column pointers
  PyArrayObject *pp = (PyArrayObject *)PyDict_GetItemString(csDict, "p");
  (*sm)->p = (long int *)PyArray_DATA( pp );

  // i - row positions
  PyArrayObject *pi = (PyArrayObject *)PyDict_GetItemString(csDict, "i");
  (*sm)->i = (long int *)PyArray_DATA( pi );
  
  return 1;
}

PyObject*
CreateDMpermOutput( csd* dm, cs* sm )
{
  PyObject* p=PyArray_SimpleNewFromData(1, &(sm->m), NPY_INT64, (void *)dm->p); 
  PyObject* q=PyArray_SimpleNewFromData(1, &(sm->n), NPY_INT64, (void *)dm->q);

  npy_intp numBlocks = dm->nb+1;
  PyObject* r=PyArray_SimpleNewFromData(1, &numBlocks, NPY_INT64, (void *)dm->r);  
  PyObject* s=PyArray_SimpleNewFromData(1, &numBlocks, NPY_INT64, (void *)dm->s);

  numBlocks = 5;
  PyObject* rr=PyArray_SimpleNewFromData(1, &numBlocks, NPY_INT64, (void *)dm->rr);
  PyObject* cc=PyArray_SimpleNewFromData(1, &numBlocks, NPY_INT64, (void *)dm->cc);

  return Py_BuildValue("(O,O,O,O,O,O)", p, q, r, s, rr, cc);  
  //  return Py_BuildValue("{s:O,s:O,s:O,s:O,s:O,s:O}", "p", p, "q", q, "r", r, "s", s, "rr", rr, "cc", cc);
}
