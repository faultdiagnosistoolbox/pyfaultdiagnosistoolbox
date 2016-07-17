
extern "C" {
  #include <Python.h>
  #include "cs.h"
  #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
  #include <numpy/ndarrayobject.h>
}
#include <iostream>
#include <cstring>
#include "StructuralAnalysisModel.h"
#include "MSOAlg.h"
#include "timer.h"

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
  //  return Py_BuildValue("[i]",10);
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

PyMODINIT_FUNC
initstructuralanalysis(void)
{
  static PyMethodDef StructuralAnalysisMethods[] = {
    {"dmperm_internal",  structuralanalysis_dmperminternal, METH_VARARGS, "Dulmage-Mendelsohn decomposition"},
    {"findmso_internal",  structuralanalysis_findmsointernal, METH_VARARGS, "Compute MSO sets"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
  };
  (void) Py_InitModule("structuralanalysis", StructuralAnalysisMethods);
  import_array();
}

int
main(int argc, char *argv[])
{
  /* Pass argv[0] to the Python interpreter */
  Py_SetProgramName(argv[0]);
  
  /* Initialize the Python interpreter.  Required. */
  Py_Initialize();
  
  /* Add a static module */
  initstructuralanalysis();

  return 0;
}


int
DictToCS( PyObject *csDict, cs **sm )
{
  // Alloc space for sparse matrix representation
  *sm = (cs *)cs_malloc (1, sizeof(cs));

  // Extract data

  // nzmax - number of non-zero elements
  PyObject *pnzmax = PyDict_GetItemString(csDict, "nzmax");
  (*sm)->nzmax = PyInt_AS_LONG(pnzmax);

  // m - number of rows
  PyObject *pm = PyDict_GetItemString(csDict, "m");
  (*sm)->m = PyInt_AS_LONG(pm);

  // n - number of columns
  PyObject *pn = PyDict_GetItemString(csDict, "n");
  (*sm)->n = PyInt_AS_LONG(pn);

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

  return Py_BuildValue("{s:O,s:O,s:O,s:O,s:O,s:O}", "p", p, "q", q, "r", r, "s", s, "rr", rr, "cc", cc);
}
