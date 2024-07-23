// Compile with: python setup.py build_ext --inplace
extern "C" {
  #include "Python.h"
  #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
  #include <numpy/ndarrayobject.h>
}
#include <math.h>
#include <vector>
#include <string>

typedef PyArrayObject* PyArrayObjectp;

static PyObject* BatchFaultIsolation(PyObject *self, PyObject * args);

// BatchFaultIsolation(r, J, FSM, dc, N)
static PyObject*
BatchFaultIsolation(PyObject *self, PyObject * args)
{
  // BatchFaultIsolation(r, J, FSM, dc, N)
  PyObject *py_r;
  PyArrayObject *py_J;
  PyArrayObject *py_FSM;
  PyObject *py_dc;
  long int N_streak;

  if (!PyArg_ParseTuple(args, "O!O!O!O!l", &PyList_Type, &py_r, &PyArray_Type, &py_J,
			&PyArray_Type, &py_FSM, &PyList_Type, &py_dc, &N_streak)) {
    return NULL;
  }

  npy_intp *FSMshape = PyArray_SHAPE( py_FSM );
  long int n_r = FSMshape[0];
  long int n_f = FSMshape[1];
  long int n_dc = (long int)PyList_Size(py_dc);

  // J
  std::vector<double> J;
  for( npy_intp i=0; i < n_r; i++ ) {
    double *Jip = (double *)PyArray_GETPTR1(py_J, i);
    J.push_back(*Jip);
  }

  // dc
  std::vector<std::string> dc;
  //  PySys_WriteStdout("|dc| = %ld\n", n_dc);
  for( Py_ssize_t i=0; i < n_dc; i++ ) {
    const char *fm = PyBytes_AsString(PyUnicode_AsASCIIString(PyList_GetItem(py_dc, i )));
    dc.push_back(std::string(fm));
  }

  PyObject *dx = PyList_New(n_dc);
  PyArrayObjectp* res = new PyArrayObjectp[n_r];
  long* alarm_streak = new long[n_r];  
  for( int tc = 0; tc < n_dc; tc++ ) { // Iterate over all test cases, tc - test case index
    std::string testCase = dc[tc];
    //    PySys_WriteStdout( "Testcase: %s\n", testCase.c_str());

    // Get object pointers to all residuals
    //    PyArrayObjectp res[n_r];
    //    long alarm_streak[n_r];

    for( int j=0; j < n_r; j++ ) {
      res[j] = (PyArrayObject *)PyDict_GetItemString(PyList_GetItem(py_r, j ), dc[tc].c_str());
      alarm_streak[j] = 0;
    }
    npy_intp* Np = PyArray_SHAPE(res[0]);
    long N = (long)(*Np);

    // Create output space for this test case
    npy_intp tc_shape[2] = {n_f, N};
    PyArrayObject* dec_k = (PyArrayObject *)PyArray_SimpleNew(2, tc_shape, NPY_INT64);
    
    // Iterate over all time points, k - time
    for( long k=0; k < N; k++ ) {
      for( int f=0; f < n_f; f++ ) {
	long int *dec_kf = (long int *)PyArray_GETPTR2(dec_k, f, k);
	*dec_kf = 1; // Initialize all to diagnosis to present
      }

      bool any_alarm = false;
      for( int r=0; r < n_r; r++ ) { // Iterate over all residuals r at time k
	double *rkp = (double *)PyArray_GETPTR1(res[r], k);
	double rk = *rkp;

	if( fabs(rk)>J[r] ) {
	  alarm_streak[r] += 1;

	  // We have an alarm in residual r, remove all faults in dec_k corresponding to 0:s in FSM(r,:)
	  if( alarm_streak[r] >= N_streak ) { 
	    any_alarm = true;
	    for( int f=0; f < n_f; f++ ) {
	      long int *FSMfrp = (long int *)PyArray_GETPTR2(py_FSM, r, f);
	      if( *FSMfrp == 0 ) {
		long int *dec_kf = (long int *)PyArray_GETPTR2(dec_k, f, k);
		*dec_kf = 0;
	      }
	    }
	  }
	} else {
	  alarm_streak[r] = 0;
	}	  
      }
      if( !any_alarm ) {
	// If no alarm, set all diagoses to 0
	for( int f=0; f < n_f; f++ ) {
	  long int *dec_kf = (long int *)PyArray_GETPTR2(dec_k, f, k);
	  *dec_kf = 0;
	}
      }
    }
    PyList_SetItem(dx, (Py_ssize_t)tc, (PyObject *)dec_k);
  }
  delete[] res;
  delete[] alarm_streak;
  return dx;
}

struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyMethodDef batchfaultisolation_methods[] = {
    {"BatchFaultIsolation",  BatchFaultIsolation, METH_VARARGS, "Batch Fault Isolation"},
    {NULL, NULL}
};

static int batchfaultisolation_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int batchfaultisolation_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "batchfaultisolation",
        NULL,
        sizeof(struct module_state),
        batchfaultisolation_methods,
        NULL,
        batchfaultisolation_traverse,
        batchfaultisolation_clear,
        NULL
};

#define INITERROR return NULL
PyMODINIT_FUNC PyInit_batchfaultisolation(void)
{
  PyObject *module = PyModule_Create(&moduledef);
  
  if (module == NULL)
    INITERROR;
  struct module_state *st = GETSTATE(module);
  
  char errName[] = "batchfaultisolation.Error";
  st->error = PyErr_NewException(errName, NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    INITERROR;
  }
  
  import_array();
  
  return module;
}
