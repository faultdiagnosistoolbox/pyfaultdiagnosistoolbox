//
//  StructuralAnalysisModel.h
//  MSO
//
//  Created by Erik Frisk on 24/08/15.
//  Copyright (c) 2015 Linköping University. All rights reserved.
//

#ifndef StructuralAnalysisModel_H
#define StructuralAnalysisModel_H

#include "SparseMatrix.h"
#include <algorithm>
#ifdef MATIO
#include "matio.h"
#endif

using EqList = std::list<int>;

struct EqClass {
    std::list<int> eq{};
    std::list<int> var{};
};

class StructuralAnalysisModel : public SparseMatrix {
protected:
    //! Which equations correspond to each row in the matrix.
    std::list<EqList> eqList{};

    void InitEqList();

public:
    // Constructors
    StructuralAnalysisModel() : SparseMatrix() { eqList.clear(); };
    explicit StructuralAnalysisModel(const std::string &s) : SparseMatrix(s) { InitEqList(); };
    explicit StructuralAnalysisModel(cs *a) : SparseMatrix(a) { InitEqList(); };
    explicit StructuralAnalysisModel(const SparseMatrix a) : SparseMatrix(a) { InitEqList(); };

    StructuralAnalysisModel &operator=(const SparseMatrix &x);
    virtual ~StructuralAnalysisModel() = default; // Needed?

    void RawPrint() { SparseMatrix::RawPrint(); };
    void Print();

    std::list<EqList>::iterator EqBegin() { return eqList.begin(); };
    std::list<EqList>::iterator EqEnd() { return eqList.end(); };

    // Member functions
    void GetEqClass(int e, EqClass &res);
    int Redundancy();

    void RemoveRow(int e);
    void RemoveRow(cs *lsm, int e);
    void LumpRows(std::list<int> &rows);
    void LumpEqClass(EqClass &res);

    void GetRows(std::list<int>::iterator startRow, std::list<int>::iterator stopRow);
    void GetCols(std::list<int>::iterator startCol, std::list<int>::iterator stopCol);
    void Get(std::list<int>::iterator startRow, std::list<int>::iterator stopRow, std::list<int>::iterator startCol,
             std::list<int>::iterator stopCol);
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

    void Plus();
    int Plus(std::list<int> &rows, std::list<int> &cols);

#ifdef MATIO
    int ExportToMatlab(std::string fileName, std::string XvarName, std::string EvarName);
#endif
};

#endif // StructuralAnalysisModel_H
