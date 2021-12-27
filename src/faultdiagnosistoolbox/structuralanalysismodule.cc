extern "C" {
  #include "Python.h"
  #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
  #include <numpy/ndarrayobject.h>
}

#include <iostream>
#include <cstring>
#include "StructuralAnalysisModel.h"
#include "MSOAlg.h"


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
structuralanalysis_findmsointernal(PyObject *self, PyObject * args)
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
structuralanalysis_dmperminternal(PyObject *self, PyObject * args)
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

static PyMethodDef structuralanalysis_methods[] = {
    {"dmperm_internal",  structuralanalysis_dmperminternal, METH_VARARGS, "Dulmage-Mendelsohn decomposition"},
    {"findmso_internal",  structuralanalysis_findmsointernal, METH_VARARGS, "Compute MSO sets"},
    {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int structuralanalysis_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int structuralanalysis_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "structuralanalysis",
        NULL,
        sizeof(struct module_state),
        structuralanalysis_methods,
        NULL,
        structuralanalysis_traverse,
        structuralanalysis_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_structuralanalysis(void)

#else
#define INITERROR return

PyMODINIT_FUNC
initstructuralanalysis(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("structuralanalysis", structuralanalysis_methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    char errName[] = "structuralanalysis.Error";
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

int
DictToCS( PyObject *csDict, cs **sm )
{
  // Alloc space for sparse matrix representation
  *sm = (cs *)cs_malloc (1, sizeof(cs));

  // Extract data

  // nzmax - number of non-zero elements
  PyObject *pnzmax = PyDict_GetItemString(csDict, "nzmax");
  (*sm)->nzmax = PyLong_AsLong(pnzmax);
  
  // m - number of rows
  PyObject *pm = PyDict_GetItemString(csDict, "m");
  (*sm)->m = PyLong_AsLong(pm);
  
  // n - number of columns
  PyObject *pn = PyDict_GetItemString(csDict, "n");
  (*sm)->n = PyLong_AsLong(pn);
  
  // nz - -1 for compressed format
  (*sm)->nz = -1;
  
  // x - data
  PyArrayObject *px = (PyArrayObject *)PyDict_GetItemString(csDict, "x");
  (*sm)->x = (double *)PyArray_DATA( px );

  // p - column pointers
  PyArrayObject *pp = (PyArrayObject *)PyDict_GetItemString(csDict, "p");
  (*sm)->p = (ptrdiff_t *)PyArray_DATA( pp );

  // i - row positions
  PyArrayObject *pi = (PyArrayObject *)PyDict_GetItemString(csDict, "i");
  (*sm)->i = (ptrdiff_t *)PyArray_DATA( pi );
  
  return 1;
}

PyObject*
CreateDMpermOutput( csd* dm, cs* sm )
{
  PyObject* p=PyArray_SimpleNewFromData(1, (npy_intp *)&(sm->m), NPY_INT64, (void *)dm->p); 
  PyObject* q=PyArray_SimpleNewFromData(1, (npy_intp *)&(sm->n), NPY_INT64, (void *)dm->q);

  npy_intp numBlocks = dm->nb+1;
  PyObject* r=PyArray_SimpleNewFromData(1, &numBlocks, NPY_INT64, (void *)dm->r);  
  PyObject* s=PyArray_SimpleNewFromData(1, &numBlocks, NPY_INT64, (void *)dm->s);

  numBlocks = 5;
  PyObject* rr=PyArray_SimpleNewFromData(1, &numBlocks, NPY_INT64, (void *)dm->rr);
  PyObject* cc=PyArray_SimpleNewFromData(1, &numBlocks, NPY_INT64, (void *)dm->cc);

  // Create matching vector m[j] = i if variable j is matched in equation i, -1 otherwise
  npy_intp *data_m = (npy_intp *)cs_malloc (sm->m, sizeof(npy_intp));
  for(auto it=data_m; it < data_m + sm->m; it++) {
    *it = -1;
  }

  // m[p[rr[0]:rr[1]]] = q[cc[1]:cc[2]]
  for(auto i=0; i < dm->rr[1] - dm->rr[0]; i++) {
    data_m[dm->p[i + dm->rr[0]]] = dm->q[dm->cc[1] + i];
  }
  // m[p[rr[1]:rr[2]]] = q[cc[2]:cc[3]]
  for(auto i=0; i < dm->rr[2] - dm->rr[1]; i++) {
    data_m[dm->p[i + dm->rr[1]]] = dm->q[dm->cc[2] + i];
  }
  // m[p[rr[2]:rr[3]]] = q[cc[3]:cc[4]]
  for(auto i=0; i < dm->rr[3] - dm->rr[2]; i++) {
    data_m[dm->p[i + dm->rr[2]]] = dm->q[dm->cc[3] + i];
  }

  PyObject* m=PyArray_SimpleNewFromData(1, (npy_intp *)&(sm->m), NPY_INT64, (void *)data_m);

  return Py_BuildValue("(O,O,O,O,O,O,O)", p, q, r, s, cc, rr, m);
}
