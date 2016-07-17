
extern "C" {
  #include <Python.h>
  #include "cs.h"
  #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
  #include <numpy/ndarrayobject.h>
}
#include <iostream>


int DictToCS( PyObject *csDict, cs **sm );
PyObject* CreateOutput( csd* dm, cs* sm );

static PyObject*
structuralanalysis_dmperminternal(PyObject *self, PyObject * args)
{
  PyObject *x;

  if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &x)) {
    return NULL;
  }

  cs* a;
  DictToCS( x, &a );
  csd *dm = cs_dmperm(a, 0);
  PyObject* res=CreateOutput(dm,a);
  
  //  free( a );

  Py_DECREF((PyArrayObject *)PyDict_GetItemString(x, "x"));  
  Py_DECREF((PyArrayObject *)PyDict_GetItemString(x, "i"));  
  Py_DECREF((PyArrayObject *)PyDict_GetItemString(x, "p"));  

  // free( dm );

  return res;
}

PyMODINIT_FUNC
initstructuralanalysis(void)
{
  static PyMethodDef StructuralAnalysisMethods[] = {
    {"dmperm_internal",  structuralanalysis_dmperminternal, METH_VARARGS, "Dulmage-Mendelsohn decomposition"},
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
  //  Py_DECREF(pnzmax);

  // m - number of rows
  PyObject *pm = PyDict_GetItemString(csDict, "m");
  (*sm)->m = PyInt_AS_LONG(pm);
  // Py_DECREF(pm);

  // n - number of columns
  PyObject *pn = PyDict_GetItemString(csDict, "n");
  (*sm)->n = PyInt_AS_LONG(pn);
  // Py_DECREF(pn);

  // nz - -1 for compressed formate
  (*sm)->nz = -1;
  
  // x - data
  PyArrayObject *px = (PyArrayObject *)PyDict_GetItemString(csDict, "x");
  (*sm)->x = (double *)PyArray_DATA( px );
  Py_INCREF(px);

  // p - column pointers
  PyArrayObject *pp = (PyArrayObject *)PyDict_GetItemString(csDict, "p");
  (*sm)->p = (long int *)PyArray_DATA( pp );
  Py_INCREF(pp);

  // i - row positions
  PyArrayObject *pi = (PyArrayObject *)PyDict_GetItemString(csDict, "i");
  (*sm)->i = (long int *)PyArray_DATA( pi );
  Py_INCREF(pi);
  
  return 1;
}

PyObject*
CreateOutput( csd* dm, cs* sm )
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

