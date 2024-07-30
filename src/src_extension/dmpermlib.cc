extern "C" {
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <cs.h>
#include <numpy/ndarrayobject.h>
}

#include "MSOAlg.h"
#include <array>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/list.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;
using nd1vec = nb::ndarray<nb::numpy, const npy_int64, nb::ndim<1>>; // convenience type alias for 1D numpy array

/// \brief Initializes the NumPy library.
///
/// This function imports the NumPy library using the import_array() function.
/// \return A null pointer.
void *
init_numpy()
{
    import_array();
    return nullptr;
}

/// Converts a numpy array object to a compressed sparse column (CSC) matrix.
/// The function extracts the necessary information from the numpy array object
/// and creates a CSC matrix with the corresponding data.
///
/// @param o The numpy array object to convert to a CSC matrix.
/// @return The converted CSC matrix.
cs
get_cs_matrix(const nb::object &o)
{
    PyObject *obj = o.ptr();
    auto shape = nb::tuple(nb::getattr(obj, "shape"));
    auto nnz = nb::cast<csi>(nb::getattr(obj, "nnz"));
    auto m = nb::cast<csi>(shape[0]);
    auto n = nb::cast<csi>(shape[1]);
    auto indptr = nb::getattr(obj, "indptr").ptr();
    auto indices = nb::getattr(obj, "indices").ptr();

    if (!PyArray_Check(indptr) || !PyArray_Check(indices)) {
        PyErr_SetString(PyExc_TypeError, "Expected numpy array");
        throw nb::python_error();
    }

    auto *indptr_pyarray = reinterpret_cast<PyArrayObject *>(indptr);   // Now numpy array object
    auto *indices_pyarray = reinterpret_cast<PyArrayObject *>(indices); // Now numpy array object

    if (PyArray_TYPE(indptr_pyarray) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected numpy array of type NPY_INT32");
        throw nb::python_error();
    }
    auto indptr_data = static_cast<int *>(PyArray_DATA(indptr_pyarray));   // Get pointer to data
    auto indices_data = static_cast<int *>(PyArray_DATA(indices_pyarray)); // Get pointer to data

    cs sparse_matrix{.nzmax = nnz, .m = m, .n = n, .p = nullptr, .i = nullptr, .x = nullptr, .nz = -1};

    // Copy data
    sparse_matrix.p = new csi[n + 1];
    sparse_matrix.i = new csi[nnz];
    sparse_matrix.x = new double[nnz];
    std::copy(indptr_data, indptr_data + n + 1, sparse_matrix.p);
    std::copy(indices_data, indices_data + nnz, sparse_matrix.i);
    std::fill(sparse_matrix.x, sparse_matrix.x + nnz, 1.0);

    return sparse_matrix;
}

std::tuple<nd1vec, nd1vec, nd1vec, nd1vec, nd1vec, nd1vec, nd1vec>
create_dmperm_output(const cs &sm, const csd *dm)
{
    // Create matching vector m[j] = i if variable j is matched in equation i, -1 otherwise
    auto data_m = new npy_intp[sm.m];
    std::fill(data_m, data_m + sm.m, -1);

    // m[p[rr[0]:rr[1]]] = q[cc[1]:cc[2]]
    for (auto i = 0; i < dm->rr[1] - dm->rr[0]; i++) {
        data_m[dm->p[i + dm->rr[0]]] = dm->q[dm->cc[1] + i];
    }
    // m[p[rr[1]:rr[2]]] = q[cc[2]:cc[3]]
    for (auto i = 0; i < dm->rr[2] - dm->rr[1]; i++) {
        data_m[dm->p[i + dm->rr[1]]] = dm->q[dm->cc[2] + i];
    }
    // m[p[rr[2]:rr[3]]] = q[cc[3]:cc[4]]
    for (auto i = 0; i < dm->rr[3] - dm->rr[2]; i++) {
        data_m[dm->p[i + dm->rr[2]]] = dm->q[dm->cc[3] + i];
    }

    // Create output vectors p, q, r, s, rr, cc, m
    // size_t shape[1] = {static_cast<size_t>(sm.m)};
    std::array<size_t, 1> shape = {static_cast<size_t>(sm.m)};

    nd1vec p(dm->p, 1, shape.data(), nb::handle());

    shape[0] = static_cast<size_t>(sm.n);
    nd1vec q(dm->q, 1, shape.data(), nb::handle());

    shape[0] = static_cast<size_t>(dm->nb + 1);
    nd1vec r(dm->r, 1, shape.data(), nb::handle());
    nd1vec s(dm->s, 1, shape.data(), nb::handle());

    shape[0] = 5;
    nd1vec rr(dm->rr, 1, shape.data(), nb::handle());
    nd1vec cc(dm->cc, 1, shape.data(), nb::handle());

    shape[0] = static_cast<size_t>(sm.m);
    nd1vec m(data_m, 1, shape.data(), nb::handle());

    // Do not cleanup data, ndarray is responsible for it
    // delete[] data_m;
    return {p, q, r, s, cc, rr, m};
}

/// @brief Compute the Dulmage-Mendelsohn permutation for a given sparse matrix.
/// @param o  - nanonbind python object representing the sparse matrix X.
/// @return p, q, r, s, cc, rr, m where
/// p and q are row and column permutation vectors, respectively, such that X(p, q) has a block upper triangular form. r
/// and s are index vectors indicating the block boundaries for the fine decomposition. cc and rr are vectors of length
/// five indicating the block boundaries of the coarse decomposition.
nb::tuple
dmperm(const nb::object &o)
{
    const cs sm = get_cs_matrix(o);
    csd *dm = cs_dmperm(&sm, 0);

    // Do not cleanup data, ndarray is responsible for it
    // cs_dfree(dm);
    const auto [p, q, r, s, cc, rr, m] = create_dmperm_output(sm, dm);
    return nb::make_tuple(p, q, r, s, cc, rr, m);
}

/// @brief Compute the MSO (Minimally Structural Overdetermined) sets for a given sparse matrix.
/// @param o - nanonbind python object representing the sparse matrix.
MSOList
MSO(const nb::object &o)
{
    auto sparse_matrix = get_cs_matrix(o);

    const StructuralAnalysisModel model(&sparse_matrix);
    MSOResult result;
    MSOAlg alg(model);
    alg.MSO(result);

    // cleanup
    delete[] sparse_matrix.p;
    delete[] sparse_matrix.i;
    delete[] sparse_matrix.x;
    return result.msos;
}

NB_MODULE(dmpermlib, m)
{
    m.def("MSO", &MSO);
    m.def("dmperm", &dmperm);
    init_numpy(); // Initialize numpy before using the C-API:
                  // https://numpy.org/doc/stable/reference/c-api/array.html#importing-the-api
}