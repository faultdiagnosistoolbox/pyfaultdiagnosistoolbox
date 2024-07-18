#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>
#include <list>
#include <string>
#include <vector>
#ifdef MATIO
#include <matio.h>
#endif

extern "C" {
#include "cs.h"
//  #include <stdlib.h>
}

struct DMPermResult {
    std::vector<int> p{};
    std::vector<int> q{};
    std::vector<int> r{};
    std::vector<int> s{};
    int nb{};
    std::vector<int> rr{};
    std::vector<int> cc{};
};

class SparseMatrix {
protected:
    cs *sm;

    // Private membership functions
    cs *CSCopy(cs *a);
    cs *RowSelectionMatrix(std::list<int>::iterator startRow, std::list<int>::iterator stopRow, int nrows);
    cs *ColSelectionMatrix(std::list<int>::iterator startCol, std::list<int>::iterator stopCol, int ncols);
    cs *RowDropMatrix(std::list<int>::iterator startRow, std::list<int>::iterator stopRow, int nrows);
    cs *ColDropMatrix(std::list<int>::iterator startCol, std::list<int>::iterator stopCol, int ncols);
    friend std::ostream &operator<<(std::ostream &s, SparseMatrix m);

    void RemoveRow(int e);
    void RemoveRow(cs *lsm, int e);

public:
    SparseMatrix() : sm(nullptr){};
    explicit SparseMatrix(const std::string &s);
    explicit SparseMatrix(cs *a);
    SparseMatrix(const SparseMatrix &y);
    SparseMatrix &operator=(const SparseMatrix &y);
    virtual ~SparseMatrix() { cs_spfree(sm); };

    bool IsEmpty();
    void RawPrint();
    void Print() const { Print(std::cout); };
    void Print(std::ostream &s) const;

    void GetRows(std::list<int>::iterator startRow, std::list<int>::iterator stopRow);
    void GetCols(std::list<int>::iterator startCol, std::list<int>::iterator stopCol);
    void Get(std::list<int>::iterator startRow, std::list<int>::iterator stopRow, std::list<int>::iterator startCol,
             std::list<int>::iterator stopCol);
    //! Get rows and columns
    inline void Get(std::list<int> &rows, std::list<int> &cols)
    {
        Get(rows.begin(), rows.end(), cols.begin(), cols.end());
    };
    inline void GetRows(std::list<int> &rows) { GetRows(rows.begin(), rows.end()); };
    inline void GetCols(std::list<int> &cols) { GetCols(cols.begin(), cols.end()); };

    void DropRows(std::list<int>::iterator startRow, std::list<int>::iterator stopRow);
    void DropCols(std::list<int>::iterator startCol, std::list<int>::iterator stopCol);
    void Drop(std::list<int>::iterator startRow, std::list<int>::iterator stopRow, std::list<int>::iterator startCol,
              std::list<int>::iterator stopCol);
    inline void Drop(std::list<int> &rows, std::list<int> &cols)
    {
        Drop(rows.begin(), rows.end(), cols.begin(), cols.end());
    };
    inline void DropRows(std::list<int> &rows) { DropRows(rows.begin(), rows.end()); };
    inline void DropCols(std::list<int> &cols) { DropCols(cols.begin(), cols.end()); };

    inline void Permute(const DMPermResult &res) { Permute(res.p, res.q); };

    void Permute(const std::vector<int> &rowp, const std::vector<int> &colp);

    void DMPerm(DMPermResult &res);

    int SRank() const;

    void FullIncidenceMatrix(int *im);

    int ExportToMatlab(const std::string &fileName, const std::string &varName) const;

    int NRows() const { return ((int)sm->m); };
    int NCols() const { return ((int)sm->n); };
    int nz() const { return ((int)sm->nzmax); };
    void GetInfo(int &nRows, int &nCols, int &nz) const;

    int Plus(std::list<int> &rows, std::list<int> &cols);
    void Plus();
    void DropNonCausal();
};

#endif // SPARSEMATRIX_H
