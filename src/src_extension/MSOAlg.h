//
//  msoalg.h
//  MSO
//
//  Created by Erik Frisk on 10/09/15.
//  Copyright (c) 2015 Linköping University. All rights reserved.
//

#ifndef MSOAlg_H
#define MSOAlg_H

#include "StructuralAnalysisModel.h"
#include <algorithm>
#include <iostream>

#ifndef NDEBUG
#include <cassert>
#endif

using MSOList = std::list<EqList>;

class MSOResult {
protected:
    int mode{0};
    unsigned long numMSOs{0};
    int verbN{0};

public:
    MSOList msos{};
    MSOResult() : mode(0), numMSOs(0), verbN(-1){};
    void Clear()
    {
        msos.clear();
        numMSOs = 0;
    };
    void AddMSO(std::list<EqList>::iterator start, std::list<EqList>::iterator stop);
    void AddMSO(StructuralAnalysisModel &m) { AddMSO(m.EqBegin(), m.EqEnd()); };
    size_t Size() const
    {
        if (mode == 0) {
            return (msos.size());
        } else {
            return numMSOs;
        }
    };
#ifdef MATIO
    int ExportToMatlab(std::string s, std::string varname) const;
#endif
    //  void Print( );
    void CountMode() { mode = 1; };
    void MSOMode() { mode = 0; };
    void VerboseN(int n) { verbN = n; };

    void RemoveNonCausal(const SparseMatrix &m);

    MSOList::iterator begin() { return (msos.begin()); };
    MSOList::iterator end() { return (msos.end()); };
};

class MSOAlg {
protected:
    StructuralAnalysisModel SM{};

    EqList R{};

    // private member functions
    bool SubsetQ(const EqList &R1, int e);
    bool SubsetQ(const EqList &R1, const EqList &R2);
    void SetDiff(EqList &R1, EqList R2);
    void UpdateIndexListAfterLump(EqList &R, EqList &lEq);
    void InitR();
    // void RemoveEquation(int e);
    void RemoveNextEquation();
    void FindMSO(MSOResult &msos);
    void LumpModel();
    bool CausalPSO();

public:
    MSOAlg() : SM(){};
    explicit MSOAlg(const std::string &s) : SM(s) { InitR(); };
    explicit MSOAlg(const SparseMatrix &a) : SM(a) { InitR(); };
    explicit MSOAlg(const StructuralAnalysisModel &a) : SM(a) { InitR(); };

    MSOAlg &operator=(const StructuralAnalysisModel &x);
    MSOAlg &operator=(const SparseMatrix &x);
    virtual ~MSOAlg() = default;

    void MSO(MSOResult &msos);
};

#endif // MSOAlg_H
