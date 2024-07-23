"""Code generation."""

import numpy as np
import sympy as sym
import sys
import faultdiagnosistoolbox as fdt
from sympy.printing import ccode, octave_code
import time
import os
import platform


def UsedVars(rels, variables):
    """Determine used variables in expression."""
    uv = set([])
    for expr in rels:
        if not fdt.IsDifferentialConstraint(expr):
            uv |= set([str(v) for v in expr.atoms() if str(v) in variables])
    return [v for v in variables if v in uv]


def ExprToCode(expr, language, user_functions=None):
    """Convert symbolic expression to code."""
    genCode = ""
    usr_fun = {} if user_functions is None else user_functions

    if language == "C":
        genCode = ccode(expr, user_functions=usr_fun)
    elif language == "Matlab":
        genCode = octave_code(expr, user_functions=usr_fun)
    elif language == "Python":
        genCode = str(expr)
    else:
        print("Unknown language")
    return genCode


def CodeEOL(language):
    """Return EOL character for different languages."""
    if language == "C" or language == "Matlab":
        return ";"
    elif language == "Python":
        return ""
    return ""


def CodeComment(language):
    """Return comment characters for different languages."""
    if language == "C":
        return "//"
    elif language == "Matlab":
        return "%"
    elif language == "Python":
        return "#"
    return "#"


def DVar(diffConstraint):
    """Return differentiated variable from differential constraint."""
    return diffConstraint[0]


def IVar(diffConstraint):
    """Return integrated variable from differential constraint."""
    return diffConstraint[1]


def SeqResGenCausality(Gamma, e, diffres):
    """Return causality of matching."""
    if not fdt.IsDifferentialConstraint(e):
        return Gamma.matchType
    elif diffres == "int":
        # Treat differential residual constraint in integral causality
        if Gamma.matchType == "der":
            return "mixed"
        elif Gamma.matchType == "int" or Gamma.matchType == "mixed":
            return Gamma.matchType
        elif Gamma.matchType == "algebraic":
            return "int"
    else:
        # Treat differential residual constraint in derivative causality
        if Gamma.matchType == "int":
            return "mixed"
        elif Gamma.matchType == "der" or Gamma.matchType == "mixed":
            return Gamma.matchType
        elif Gamma.matchType == "algebraic":
            return "der"


def AlgebraicHallComponent(model, g, language, user_functions=None):
    """Generate code for Algebraic Hall component."""
    fzSub = list(zip(model.f, np.zeros(len(model.f), dtype=np.int64)))
    eqs = list(map(lambda eq: eq.subs(fzSub), model.syme[g.row]))
    #    hcVars = list(np.array(model.x)[g.col])
    hcVars = [sym.var(vi) for vi in np.array(model.x)[g.col]]
    sol = sym.solve(eqs, hcVars, dict=True)
    sol = sol[0]  # Take the first solution
    resGen = []
    for varSol in sol:
        if len(hcVars) == 1:
            genCode = "{} = {}{} {} {}".format(
                varSol,
                ExprToCode(sol[varSol], language, user_functions),
                CodeEOL(language),
                CodeComment(language),
                model.e[g.row[0]],
            )
        else:
            genCode = "{} = {}{}".format(varSol, ExprToCode(sol[varSol], language, user_functions), CodeEOL(language))
        resGen.append(genCode)
    return resGen


def CodeApproxInt(v, dv, enum, language):
    """Generate code for simple approximate integrator."""
    if language == "Matlab":
        return "ApproxInt({}, state.{}, Ts){} {} {}".format(dv, v, CodeEOL(language), CodeComment(language), enum)
    elif language == "C":
        return "ApproxInt({}, state->{}, Ts){} {} {}".format(dv, v, CodeEOL(language), CodeComment(language), enum)
    else:
        return "ApproxInt({}, state['{}'], Ts){} {} {}".format(dv, v, CodeEOL(language), CodeComment(language), enum)


def CodeApproxDer(v, enum, language):
    """Generate code for simple approximate differentiator."""
    if language == "Matlab":
        return "ApproxDiff({}, state.{}, Ts){} {} {}".format(v, v, CodeEOL(language), CodeComment(language), enum)
    elif language == "C":
        return "ApproxDiff({}, state->{}, Ts){} {} {}".format(v, v, CodeEOL(language), CodeComment(language), enum)
    else:
        return "ApproxDiff({}, state['{}'], Ts){} {} {}".format(v, v, CodeEOL(language), CodeComment(language), enum)


def IntegralHallComponent(model, g, language, user_functions=None):
    """Generate code for Integral Hall component."""
    fzSub = list(zip(model.f, np.zeros(len(model.f), dtype=np.int64)))
    resGen = []
    integ = []
    iState = []
    for e, v, enum in zip(model.syme[g.row], np.array(model.x)[g.col], np.array(model.e)[g.row]):
        if not fdt.IsDifferentialConstraint(e):
            sol = sym.solve(e.subs(fzSub), sym.var(v))
            genCode = "{} = {}{} {} {}".format(
                v, ExprToCode(sol[0], language, user_functions), CodeEOL(language), CodeComment(language), enum
            )
            resGen.append(genCode)
        else:
            dv = DVar(e)
            #            iv = IVar(e)
            genCode = v + " = " + CodeApproxInt(v, dv, enum, language)
            integ.append(genCode)
            iState.append(v)
    return resGen, integ, iState


def MixedHallComponent(model, g, language, user_functions=None):
    """Generate code for Mixed Hall component."""
    fzSub = list(zip(model.f, np.zeros(len(model.f), dtype=np.int64)))
    resGen = []
    iState = []
    dState = []
    integ = []
    for e, v, enum in zip(model.syme[g.row], np.array(model.x)[g.col], np.array(model.e)[g.row]):
        if not fdt.IsDifferentialConstraint(e):
            sol = sym.solve(e.subs(fzSub), sym.var(v))
            genCode = "{} = {}{} {} {}".format(
                v, ExprToCode(sol[0], language, user_functions), CodeEOL(language), CodeComment(language), enum
            )
            resGen.append(genCode)
        elif v == IVar(e):
            dv = DVar(e)
            #            iv = IVar(e)
            genCode = v + " = " + CodeApproxInt(v, dv, enum, language)
            integ.append(genCode)
            iState.append(v)
        else:  # v == DVar(e)
            genCode = v + " = " + CodeApproxDer(IVar(e), enum, language)
            resGen.append(genCode)
            dState.append(IVar(e))

    return resGen, integ, iState, dState


def GenerateResidualEquations(model, resEq, diffres, language, user_functions=None):
    """Generate code for the residual equations."""
    if language == "C":
        resvar = "r[0]"
    else:
        resvar = "r"

    e = model.syme[resEq]
    if not fdt.IsDifferentialConstraint(e):
        fzSub = list(zip(model.f, np.zeros(len(model.f), dtype=np.int64)))
        e = e.subs(fzSub)
        resExpr = e.lhs - e.rhs
        genCode = [
            "{} = {}{} {} {}".format(
                resvar,
                ExprToCode(resExpr, language, user_functions=user_functions),
                CodeEOL(language),
                CodeComment(language),
                np.array(model.e)[resEq],
            )
        ]
        iState = []
        dState = []
        integ = []
    else:  # IsDifferentialConstraint(e)
        if diffres == "der":
            iv = IVar(e)
            dv = DVar(e)
            genCode = ["{} = {}-{}".format(resvar, dv, CodeApproxDer(iv, np.array(model.e)[resEq], language))]
            iState = []
            dState = [dv]
            integ = []
        else:  # diffres == 'int'
            iv = IVar(e)
            dv = DVar(e)
            if language == "Python":
                genCode = [
                    "{} = {}-state['{}']{} {} {}".format(
                        resvar, iv, iv, CodeEOL(language), CodeComment(language), np.array(model.e)[resEq]
                    )
                ]
            elif language == "C":
                genCode = [
                    "{} = {}-state->{}{} {} {}".format(
                        resvar, iv, iv, CodeEOL(language), CodeComment(language), np.array(model.e)[resEq]
                    )
                ]
            else:
                genCode = [
                    "{} = {}-state.{}{} {} {}".format(
                        resvar, iv, iv, CodeEOL(language), CodeComment(language), np.array(model.e)[resEq]
                    )
                ]
            iState = [iv]
            dState = []
            integ = [iv + " = " + CodeApproxInt(iv, dv, np.array(model.e)[resEq], language)]
    return np.array(genCode), iState, dState, integ


def GenerateExactlyDetermined(model, Gamma, language, user_functions=None):
    """Generate code for exactly determined part."""
    resGenM0 = np.array([])
    integ = np.array([])
    iState = np.array([])
    dState = np.array([])

    for g in Gamma.matching:
        sys.stdout.write(".")
        sys.stdout.flush()
        if g.matchType == "algebraic":
            codeGen = AlgebraicHallComponent(model, g, language, user_functions)
            resGenM0 = np.concatenate((resGenM0, codeGen))
        elif g.matchType == "int":
            codeGen, gInteg, giState = IntegralHallComponent(model, g, language, user_functions)
            iState = np.concatenate((iState, giState))
            resGenM0 = np.concatenate((resGenM0, codeGen))
            integ = np.concatenate((integ, gInteg))
        elif g.matchType == "der":
            dc = model.syme[g.row[0]]
            codeGen = DVar(dc) + " = " + CodeApproxDer(IVar(dc), model.e[g.row[0]], language)
            resGenM0 = np.concatenate((resGenM0, [codeGen]))
            dState = np.concatenate((dState, [IVar(dc)]))
        elif g.matchType == "mixed":
            codeGen, gInteg, giState, gdState = MixedHallComponent(model, g, language, user_functions)
            resGenM0 = np.concatenate((resGenM0, codeGen))
            iState = np.concatenate((iState, giState))
            dState = np.concatenate((dState, gdState))
            integ = np.concatenate((integ, gInteg))
    sys.stdout.write("\n")
    return resGenM0, iState, dState, integ


def WriteApproxIntFunction(f, language):
    """Write code for approximate integration function."""
    if language == "C":
        f.write("double\n")
        f.write("ApproxInt(double dx, double x0, double Ts)\n")
        f.write("{\n")
        f.write("  return x0 + Ts*dx;\n")
        f.write("}\n")
    elif language == "Python":
        f.write("    def ApproxInt(dx, x0, Ts):\n")
        f.write("        return x0 + Ts*dx\n")
    elif language == "Matlab":
        pass


def WriteApproxDerFunction(f, language):
    """Write code for approximate differentiation function."""
    if language == "C":
        f.write("double\nApproxDiff(double x, double xold, double Ts)\n")
        f.write("{\n")
        f.write("  return (x-xold)/Ts;\n")
        f.write("}\n")
    elif language == "Python":
        f.write("    def ApproxDiff(x, xold, Ts):\n")
        f.write("        return (x - xold) / Ts\n")
    elif language == "Matlab":
        f.write("function dx=ApproxDiff(x, xold, Ts)\n")
        f.write("  if length(xold) == 1\n")
        f.write("    dx = (x - xold) / Ts;\n")
        f.write("  end\n")
        f.write("end\n")


def WriteResGenCorePython(f, name, model, resGen, state, integ, resGenEqs):
    """Write python skeleton for residual generator."""
    # Function definition
    f.write("    def " + name + "_core(z, state, params, Ts):\n")

    tab = "        "
    # Parameters
    usedParams = UsedVars(model.syme[resGenEqs], model.parameters)
    if len(usedParams) > 0:
        f.write(tab + "# Parameters\n")
        for p in usedParams:
            f.write(tab + p + " = params['" + p + "']\n")
        f.write("\n")

    # Known variables
    usedKnownSignals = UsedVars(model.syme[resGenEqs], model.z)
    f.write(tab + "# Known signals\n")
    for zv in usedKnownSignals:
        f.write(tab + zv + " = z[" + str(model.z.index(zv)) + "]\n")
    f.write("\n")

    # Initialize state variables
    if len(state) > 0:
        f.write(tab + "# Initialize state variables\n")
        for sv in state:
            f.write(tab + sv + " = state['" + sv + "']\n")
        f.write("\n")

    # Residual generator body
    f.write(tab + "# Residual generator body\n")
    for pyExpr in resGen:
        f.write(tab + pyExpr + "\n")
    f.write("\n")

    # Update integrator variables
    if len(integ) > 0:
        f.write(tab + "# Update integrator variables\n")
        for iExpr in integ:
            f.write(tab + iExpr + "\n")
        f.write("\n")

    # Update state variables
    if len(state) > 0:
        f.write(tab + "# Update state variables\n")
        for sv in state:
            f.write(tab + "state['" + sv + "'] = " + sv + "\n")
        f.write("\n")

    # Return value
    f.write(tab + "return (r, state)\n")


def WriteResGenPython(model, resGen, state, integ, name, batch, resGenCausality, resGenEqs, external_src):
    """Write residual generator code in Python."""
    f = open(name + ".py", "w")
    tab = "    "
    f.write("import numpy as np\n")
    f.write("from numpy import * # For access to all fundamental functions, constants etc.\n")
    if batch:
        f.write("from copy import deepcopy\n\n")

    if len(external_src) > 0:
        f.write("# External modules\n")
        for ext_mod in external_src:
            f.write("from " + os.path.splitext(ext_mod)[0] + " import *\n")
        f.write("\n")

    if not batch:
        f.write("def " + name + "(z,state,params,Ts):\n")
        f.write(tab + '""" ' + name.upper() + " Sequential residual generator for model '" + model.name + "'\n")
        f.write(tab + "Causality: " + resGenCausality + "\n")
        f.write("\n")

        fSens = np.array(model.f)[np.any(np.array(model.F[resGenEqs, :]), axis=0)]
        f.write(tab + "Structurally sensitive to faults: ")
        for fault in fSens[:-1]:
            f.write(fault + ", ")
        f.write(fSens[-1] + "\n")
        f.write("\n")
        f.write(tab + "Example of basic usage:\n")
        f.write(
            tab
            + "Let z be the observations matrix, each column corresponding to "
            + "a known signal and Ts the sampling time,\n"
        )
        f.write(tab + "then the residual generator can be simulated by:\n")
        f.write("\n")
        f.write(tab + "r = np.zeros(N) # N number of data points\n")
        f.write(tab + "state = {")
        for s in state[:-1]:
            f.write("'" + s + "': " + s + "_0, ")
        if len(state) > 0:
            f.write("'" + state[-1] + "': " + state[-1] + "_0}\n")
        f.write(tab + "for k,zk in enumerate(z):\n")
        f.write(tab + tab + "r[k], state = " + name + "( zk, state, params, Ts )\n")
        f.write("\n")
        f.write(tab + "State is a dictionary with the keys: ")
        for s in state[:-1]:
            f.write(s + ", ")
        if len(state) > 0:
            f.write(state[-1] + "\n")
        f.write("\n")
        f.write(tab + "File generated " + time.strftime("%c") + "\n")
        f.write(tab + '"""')
        f.write("\n")

        if resGenCausality == "int" or resGenCausality == "mixed":
            WriteApproxIntFunction(f, "Python")
            f.write("\n")

        if resGenCausality == "der" or resGenCausality == "mixed":
            WriteApproxDerFunction(f, "Python")
            f.write("\n")

        WriteResGenCorePython(f, name, model, resGen, state, integ, resGenEqs)

        f.write("\n")
        f.write("    return " + name + "_core(z, state, params, Ts)\n")
    else:  # batch
        f.write("def " + name + "(z,state,params,Ts):\n")
        f.write(tab + '""" ' + name.upper() + " Sequential residual generator for model '" + model.name + "'\n")
        f.write(tab + "Causality: " + resGenCausality + "\n")
        f.write("\n")

        fSens = np.array(model.f)[np.any(np.array(model.F[resGenEqs, :]), axis=0)]
        f.write(tab + "Structurally sensitive to faults: ")
        for fault in fSens[:-1]:
            f.write(fault + ", ")
        f.write(fSens[-1] + "\n")
        f.write("\n")
        f.write(tab + "Example of basic usage:\n")
        f.write(
            tab
            + "Let z be the observations matrix, each column corresponding to a known signal "
            + "and Ts the sampling time,\n"
        )
        f.write(tab + "then the residual generator can be simulated by:\n")
        f.write("\n")
        f.write(tab + "state = {")
        for s in state[:-1]:
            f.write("'" + s + "': " + s + "_0, ")
        if len(state) > 0:
            f.write("'" + state[-1] + "': " + state[-1] + "_0}\n")
        f.write(tab + "r = " + name + "( zk, state, params, Ts )\n")
        f.write("\n")
        f.write(tab + "State is a dictionary with the keys: ")
        for s in state[:-1]:
            f.write(s + ", ")
        if len(state) > 0:
            f.write(state[-1] + "\n")
        f.write("\n")
        f.write(tab + "File generated " + time.strftime("%c") + "\n")
        f.write(tab + '"""')
        f.write("\n")

        if resGenCausality == "int" or resGenCausality == "mixed":
            WriteApproxIntFunction(f, "Python")
            f.write("\n")

        if resGenCausality == "der" or resGenCausality == "mixed":
            WriteApproxDerFunction(f, "Python")
            f.write("\n")

        WriteResGenCorePython(f, name, model, resGen, state, integ, resGenEqs)

        f.write("\n")
        f.write(tab + "N = z.shape[0] # Number of data points\n")
        f.write(tab + "r = np.zeros(N)\n")
        f.write(tab + "dynState = deepcopy(state)\n")
        f.write(tab + "for k,zk in enumerate(z):\n")
        f.write(tab + tab + "r[k], dynState = " + name + "_core(zk, dynState, params, Ts)\n")
        f.write("\n")
        f.write(tab + "return r\n")
    f.close()


def SeqResGen(
    model,
    Gamma,
    resEq,
    name,
    diffres="int",
    language="Python",
    batch=False,
    api="Python",
    user_functions=None,
    external_src=None,
    external_headers=None,
):
    """Generate code for sequential residual generator."""
    if not (model.modelType == "Symbolic"):
        print("Code generation only possible for symbolic models")
        return []

    usr_fun = user_functions if user_functions is not None else {}
    ext_src = external_src if external_src is not None else []
    ext_head = external_headers if external_headers is not None else []

    sys.stdout.write("Generating residual generator " + name + " (" + language + ", ")
    if batch:
        sys.stdout.write("batch")
    else:
        sys.stdout.write("no batch")
    sys.stdout.write(")\n")

    sys.stdout.write("  Generating code for the exactly determined part: ")
    m0Code, m0iState, m0dState, m0integ = GenerateExactlyDetermined(model, Gamma, language, usr_fun)

    sys.stdout.flush()
    print("  Generating code for the residual equations")
    resCode, resiState, resdState, resinteg = GenerateResidualEquations(model, resEq, diffres, language, usr_fun)

    # Collect code, integrators, and states
    resGenCode = np.concatenate((m0Code, [" "], resCode))
    resGeniState = np.concatenate((m0iState, resiState))
    resGendState = np.concatenate((m0dState, resdState))
    resGenState = np.concatenate((resGeniState, resGendState))
    resGenInteg = np.concatenate((m0integ, resinteg))

    print("  Writing residual generator file")
    resGenCausality = SeqResGenCausality(Gamma, model.syme[resEq], diffres)
    resGenEqs = np.array([], dtype=np.int64)
    for g in Gamma.matching:
        resGenEqs = np.concatenate((resGenEqs, g.row))
    resGenEqs = np.concatenate((resGenEqs, [resEq]))
    if language == "Python":
        WriteResGenPython(model, resGenCode, resGenState, resGenInteg, name, batch, resGenCausality, resGenEqs, ext_src)
        print("File " + name + ".py generated.")
    elif language == "Matlab":
        pass
    elif language == "C":
        if api == "Python":
            WriteResGenCPython(
                model, resGenCode, resGenState, resGenInteg, name, batch, resGenCausality, resGenEqs, ext_src, ext_head
            )
            print("Files " + name + ".cc and " + name + "_setup.py generated")
            print("Compile by running: python " + name + "_setup.py build_ext --inplace")


def WriteResGenCoreC(f, name, model, resGen, state, integ, resGenEqs, batch):
    """Write code for residual generator core C function."""
    usedParams = UsedVars(model.syme[resGenEqs], model.parameters)
    tab = "  "

    # Function definition
    f.write("void\n")
    if not batch:
        f.write(name + "_core(double* r, PyArrayObject *pyZ, ")
    else:
        f.write(name + "_core(double* r, PyArrayObject *pyZ, npy_intp k, ")
    if len(state) > 0:
        f.write("ResState* state, ")
    if len(usedParams) > 0:
        f.write("const Parameters *params, ")
    f.write("double Ts)\n")
    f.write("{\n")
    # Parameters
    if len(usedParams) > 0:
        f.write(tab + "// Parameters\n")
        for p in usedParams:
            f.write(tab + "double " + p + " = params->" + p + ";\n")
        f.write("\n")

    usedVars = UsedVars(model.syme[resGenEqs], model.x)
    if len(usedVars) > 0:
        f.write(tab + "// Residual generator variables\n")
        for p in usedVars:
            f.write(tab + "double " + p + ";\n")
        f.write("\n")

    # Initialize state variables
    if len(state) > 0:
        f.write(tab + "// Initialize state variables\n")
        for sv in state:
            f.write(tab + sv + " = state->" + sv + ";\n")
        f.write("\n")

    # Known variables
    usedKnownSignals = UsedVars(model.syme[resGenEqs], model.z)
    f.write(tab + "// Known signals\n")
    for zv in usedKnownSignals:
        if not batch:
            f.write(tab + "double " + zv + " = *((double *)PyArray_GETPTR1(pyZ, " + str(model.z.index(zv)) + "));\n")
        else:
            f.write(tab + "double " + zv + " = *((double *)PyArray_GETPTR2(pyZ, k, " + str(model.z.index(zv)) + "));\n")
    f.write("\n")

    # Residual generator body
    f.write(tab + "// Residual generator body\n")
    for pyExpr in resGen:
        f.write(tab + pyExpr + "\n")
    f.write("\n")

    # Update integrator variables
    if len(integ) > 0:
        f.write(tab + "// Update integrator variables\n")
        for iExpr in integ:
            f.write(tab + iExpr + "\n")
        f.write("\n")

    # Update state variables
    if len(state) > 0:
        f.write(tab + "// Update state variables\n")
        for sv in state:
            f.write(tab + "state->" + sv + " = " + sv + ";\n")
    f.write("}\n")


def WriteCPythonAPI(f, name):
    """Write code fro C-Python API."""
    f.write("struct module_state {\n")
    f.write("    PyObject *error;\n")
    f.write("};\n")
    f.write("\n")
    f.write("#if PY_MAJOR_VERSION >= 3\n")
    f.write("#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))\n")
    f.write("#else\n")
    f.write("#define GETSTATE(m) (&_state)\n")
    f.write("static struct module_state _state;\n")
    f.write("#endif\n")
    f.write("\n")
    f.write("static PyObject *\n")
    f.write("error_out(PyObject *m) {\n")
    f.write("    struct module_state *st = GETSTATE(m);\n")
    f.write('    PyErr_SetString(st->error, "something bad happened");\n')
    f.write("    return NULL;\n")
    f.write("}\n")
    f.write("\n")
    f.write("static PyMethodDef " + name + "_methods[] = {\n")
    f.write('    {"' + name + '",  ' + name + ', METH_VARARGS, "Residual generator ' + name + '"},\n')
    f.write("    {NULL, NULL}\n")
    f.write("};\n")
    f.write("\n")
    f.write("#if PY_MAJOR_VERSION >= 3\n")
    f.write("\n")
    f.write("static int " + name + "_traverse(PyObject *m, visitproc visit, void *arg) {\n")
    f.write("    Py_VISIT(GETSTATE(m)->error);\n")
    f.write("    return 0;\n")
    f.write("}\n")
    f.write("\n")
    f.write("static int " + name + "_clear(PyObject *m) {\n")
    f.write("    Py_CLEAR(GETSTATE(m)->error);\n")
    f.write("    return 0;\n")
    f.write("}\n")
    f.write("\n")
    f.write("static struct PyModuleDef moduledef = {\n")
    f.write("        PyModuleDef_HEAD_INIT,\n")
    f.write('        "' + name + '",\n')
    f.write("        NULL,\n")
    f.write("        sizeof(struct module_state),\n")
    f.write("        " + name + "_methods,\n")
    f.write("        NULL,\n")
    f.write("        " + name + "_traverse,\n")
    f.write("        " + name + "_clear,\n")
    f.write("        NULL\n")
    f.write("};\n")
    f.write("\n")
    f.write("#define INITERROR return NULL\n")
    f.write("\n")
    f.write("PyMODINIT_FUNC\n")
    f.write("PyInit_" + name + "(void)\n")
    f.write("\n")
    f.write("#else\n")
    f.write("#define INITERROR return\n")
    f.write("\n")
    f.write("PyMODINIT_FUNC\n")
    f.write("init" + name + "(void)\n")
    f.write("#endif\n")
    f.write("{\n")
    f.write("#if PY_MAJOR_VERSION >= 3\n")
    f.write("    PyObject *module = PyModule_Create(&moduledef);\n")
    f.write("#else\n")
    f.write('    PyObject *module = Py_InitModule("' + name + '", ' + name + "_methods);\n")
    f.write("#endif\n")
    f.write("\n")
    f.write("    if (module == NULL)\n")
    f.write("        INITERROR;\n")
    f.write("    struct module_state *st = GETSTATE(module);\n")
    f.write("\n")
    f.write('    char errName[] = "' + name + '.Error";\n')
    f.write("    st->error = PyErr_NewException(errName, NULL, NULL);\n")
    f.write("    if (st->error == NULL) {\n")
    f.write("        Py_DECREF(module);\n")
    f.write("        INITERROR;\n")
    f.write("    }\n")
    f.write("\n")
    f.write("    import_array();\n")
    f.write("\n")
    f.write("#if PY_MAJOR_VERSION >= 3\n")
    f.write("    return module;\n")
    f.write("#endif\n")
    f.write("}\n")


def WriteResGenCPython(
    model, resGenCode, resGenState, resGenInteg, name, batch, resGenCausality, resGenEqs, external_src, external_headers
):
    """Write residual generator code for C."""
    f = open(name + ".cc", "w")
    tab = "  "
    if not batch:
        f.write("#include <Python.h>\n")
        f.write("#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION\n")
        f.write("#include <numpy/ndarrayobject.h>\n")

        f.write("\n#ifdef _WIN32\n")
        f.write("#define _USE_MATH_DEFINES\n")
        f.write("#include <cmath>\n")
        f.write("#endif\n\n")

        for header in external_headers:
            f.write('#include "' + header + '"\n')
        f.write("\n")

        f.write("// " + name.upper() + " Sequential residual generator for model '" + model.name + "'\n")
        f.write("// Causality: " + resGenCausality + "\n")
        f.write("//\n")

        fSens = np.array(model.f)[np.any(np.array(model.F[resGenEqs, :]), axis=0)]
        f.write("// Structurally sensitive to faults: ")
        for fault in fSens[:-1]:
            f.write(fault + ", ")
        f.write(fSens[-1] + "\n")
        f.write("//\n")
        f.write("// Example of basic usage:\n")
        f.write(
            "// Let z be the observations matrix, each column corresponding to a known signal "
            + "and Ts the sampling time,\n"
        )
        f.write("// then the residual generator can be simulated by:\n")
        f.write("//\n")
        f.write("// r = np.zeros(N) # N number of data points\n")
        f.write("// state = {")
        for s in resGenState[:-1]:
            f.write("'" + s + "': " + s + "_0, ")
        if len(resGenState) > 0:
            f.write("'" + resGenState[-1] + "': " + resGenState[-1] + "_0}\n")
        f.write("// for k,zk in enumerate(z):\n")
        f.write("//     r[k], state = " + name + "( zk, state, params, Ts )\n")
        f.write("//\n")
        f.write("// State is a dictionary with the keys: ")
        for s in resGenState[:-1]:
            f.write(s + ", ")
        if len(resGenState) > 0:
            f.write(resGenState[-1] + "\n")
        f.write("\n")
        f.write("// File generated " + time.strftime("%c") + "\n")
        f.write("\n")

        # Declare Parameters struct
        usedParams = UsedVars(model.syme[resGenEqs], model.parameters)
        if len(usedParams) > 0:
            f.write("typedef struct {\n")
            for p in usedParams:
                f.write(tab + "double " + p + ";\n")
            f.write("} Parameters;\n")
            f.write("\n")

        # Declare state struct
        if len(resGenState) > 0:
            f.write("typedef struct {\n")
            for sv in resGenState:
                f.write(tab + "double " + sv + ";\n")
            f.write("} ResState;\n")
            f.write("\n")

        # Declare function prototypes
        if resGenCausality == "int" or resGenCausality == "mixed":
            f.write("double ApproxInt(double dx, double x0, double Ts);\n")
        if resGenCausality == "der" or resGenCausality == "mixed":
            f.write("double ApproxDiff(double x, double xold, double Ts);\n")
        if len(usedParams) > 0:
            f.write("void GetParameters( PyObject *pyParams, Parameters* params );\n")
        if len(resGenState) > 0:
            f.write("void GetState( PyObject *pyStates, ResState* state );\n")

        f.write("void " + name + "_core( double* r, PyArrayObject *pyZ, ")
        if len(resGenState) > 0:
            f.write("ResState *state, ")
        if len(usedParams) > 0:
            f.write("const Parameters *params, ")
        f.write("double Ts );\n")

        f.write("\n")
        f.write("static PyObject*\n")
        f.write(name + "(PyObject *self, PyObject * args)\n")
        f.write("{\n")
        f.write("  PyArrayObject *pyZ;\n")
        f.write("  PyObject *pyState;\n")
        f.write("  PyObject *pyParams;\n")
        f.write("  double Ts;\n")
        f.write("\n")
        f.write(
            '  if (!PyArg_ParseTuple(args, "O!O!O!d", &PyArray_Type, &pyZ, &PyDict_Type, &pyState, '
            + "&PyDict_Type, &pyParams, &Ts)) {\n"
        )
        f.write("    return NULL;\n")
        f.write("  }\n")
        f.write("\n")
        f.write("  double r;\n")
        f.write("\n")
        if len(resGenState) > 0:
            f.write("  ResState state;\n")
            f.write("  GetState(pyState, &state);\n")
            f.write("\n")
        if len(usedParams) > 0:
            f.write("  Parameters params;\n")
            f.write("  GetParameters( pyParams, &params );\n")
            f.write("\n")
        f.write("  " + name + "_core( &r, pyZ, ")
        if len(resGenState) > 0:
            f.write("&state, ")
        if len(usedParams) > 0:
            f.write("&params, ")
        f.write("Ts );\n")
        f.write("\n")
        f.write("  PyObject *newState = PyDict_New();\n\n")
        for sv in resGenState:
            f.write('  PyDict_SetItemString(newState, "' + sv + '", PyFloat_FromDouble(state.' + sv + "));\n")
        f.write("\n")
        f.write('  return Py_BuildValue("(f,O)", r, newState);\n')
        f.write("}\n")

        f.write("\n")
        WriteResGenCoreC(f, name, model, resGenCode, resGenState, resGenInteg, resGenEqs, batch)
        f.write("\n")

        if resGenCausality == "int" or resGenCausality == "mixed":
            WriteApproxIntFunction(f, "C")
            f.write("\n")

        if resGenCausality == "der" or resGenCausality == "mixed":
            WriteApproxDerFunction(f, "C")
            f.write("\n")

        if len(usedParams) > 0:
            f.write("void\n")
            f.write("GetParameters( PyObject *pyParams, Parameters* params )\n")
            f.write("{\n")
            for p in usedParams:
                f.write(tab + "PyObject *" + p + ' = PyDict_GetItemString(pyParams, "' + p + '");\n')
                f.write(tab + "params->" + p + " = PyFloat_AsDouble(" + p + ");\n")
                f.write("\n")
            f.write("}\n\n")

        if len(resGenState) > 0:
            f.write("void\n")
            f.write("GetState( PyObject *pyState, ResState* state )\n")
            f.write("{\n")
            for sv in resGenState:
                f.write(tab + "PyObject *" + sv + ' = PyDict_GetItemString(pyState, "' + sv + '");\n')
                f.write(tab + "state->" + sv + " = PyFloat_AsDouble(" + sv + ");\n")
                f.write("\n")
            f.write("}\n\n")

        WriteCPythonAPI(f, name)
        WriteSetupBuild(name, external_src)
    else:  # batch
        f.write("#include <Python.h>\n")
        f.write("#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION\n")
        f.write("#include <numpy/ndarrayobject.h>\n")

        f.write("\n#ifdef _WIN32\n")
        f.write("#define _USE_MATH_DEFINES\n")
        f.write("#include <cmath>\n")
        f.write("#endif\n\n")

        for header in external_headers:
            f.write('#include "' + header + '"\n')
        f.write("\n")

        f.write("// " + name.upper() + " Sequential residual generator for model '" + model.name + "'\n")
        f.write("// Causality: " + resGenCausality + "\n")
        f.write("//\n")

        fSens = np.array(model.f)[np.any(np.array(model.F[resGenEqs, :]), axis=0)]
        f.write("// Structurally sensitive to faults: ")
        for fault in fSens[:-1]:
            f.write(fault + ", ")
        f.write(fSens[-1] + "\n")
        f.write("//\n")
        f.write("// Example of basic usage:\n")
        f.write(
            "// Let z be the observations matrix, each column corresponding to a known signal and "
            + "Ts the sampling time,\n"
        )
        f.write("// then the residual generator can be simulated by:\n")
        f.write("//\n")

        if len(resGenState) > 0:
            f.write("// state = {")
            for s in resGenState[:-1]:
                f.write("'" + s + "': " + s + "_0, ")
            f.write("'" + resGenState[-1] + "': " + resGenState[-1] + "_0}\n")
            f.write("// r = " + name + "( z, state, params, Ts )\n")
            f.write("//\n")
            f.write("// State is a dictionary with the keys: ")
            for s in resGenState[:-1]:
                f.write(s + ", ")
            f.write(resGenState[-1] + "\n")
        f.write("\n")
        f.write("// File generated " + time.strftime("%c") + "\n")
        f.write("\n")

        # Declare Parameters struct
        usedParams = UsedVars(model.syme[resGenEqs], model.parameters)
        if len(usedParams) > 0:
            f.write("typedef struct {\n")
            for p in usedParams:
                f.write(tab + "double " + p + ";\n")
            f.write("} Parameters;\n")
            f.write("\n")

        # Declare state struct
        if len(resGenState) > 0:
            f.write("typedef struct {\n")
            for sv in resGenState:
                f.write(tab + "double " + sv + ";\n")
            f.write("} ResState;\n")
            f.write("\n")

        # Declare function prototypes
        if resGenCausality == "int" or resGenCausality == "mixed":
            f.write("double ApproxInt(double dx, double x0, double Ts);\n")
        if resGenCausality == "der" or resGenCausality == "mixed":
            f.write("double ApproxDiff(double x, double xold, double Ts);\n")
        if len(usedParams) > 0:
            f.write("void GetParameters( PyObject *pyParams, Parameters* params );\n")
        if len(resGenState) > 0:
            f.write("void GetState( PyObject *pyStates, ResState* state );\n")

        f.write("void " + name + "_core( double* r, PyArrayObject *pyZ, npy_intp k, ")
        if len(resGenState) > 0:
            f.write("ResState *state, ")
        if len(usedParams) > 0:
            f.write("const Parameters *params, ")
        f.write("double Ts );\n")

        f.write("\n")
        f.write("static PyObject*\n")
        f.write(name + "(PyObject *self, PyObject * args)\n")
        f.write("{\n")
        f.write("  PyArrayObject *pyZ;\n")
        f.write("  PyObject *pyState;\n")
        f.write("  PyObject *pyParams;\n")
        f.write("  double Ts;\n")
        f.write("\n")
        f.write(
            '  if (!PyArg_ParseTuple(args, "O!O!O!d", &PyArray_Type, &pyZ, &PyDict_Type, &pyState, '
            + "&PyDict_Type, &pyParams, &Ts)) {\n"
        )
        f.write("    return NULL;\n")
        f.write("  }\n")
        f.write("\n")

        f.write("  if( PyArray_NDIM(pyZ) != 2 ) {\n")
        f.write('    return Py_BuildValue("");\n')
        f.write("  }\n")
        f.write("  npy_intp *shape = PyArray_SHAPE( pyZ );\n")
        f.write("  npy_intp N = shape[0];\n\n")

        f.write("  PyObject* pyR = PyArray_SimpleNew(1, &N, NPY_FLOAT64);\n")
        f.write("  double *r = (double *)PyArray_DATA((PyArrayObject *)pyR);\n\n")
        f.write("\n")
        if len(resGenState) > 0:
            f.write("  ResState state;\n")
            f.write("  GetState(pyState, &state);\n")
            f.write("\n")
        if len(usedParams) > 0:
            f.write("  Parameters params;\n")
            f.write("  GetParameters( pyParams, &params );\n")
            f.write("\n")

        f.write("  // Main computational loop\n")
        f.write("  for( npy_intp k=0; k < N; k++ ) {\n")
        f.write("    " + name + "_core( &(r[k]), pyZ, k, ")
        if len(resGenState) > 0:
            f.write("&state, ")
        if len(usedParams) > 0:
            f.write("&params, ")
        f.write("Ts );\n")
        f.write("}\n")
        f.write("\n")
        f.write('  return Py_BuildValue("O", pyR);\n')
        f.write("}\n")

        f.write("\n")
        WriteResGenCoreC(f, name, model, resGenCode, resGenState, resGenInteg, resGenEqs, batch)
        f.write("\n")

        if resGenCausality == "int" or resGenCausality == "mixed":
            WriteApproxIntFunction(f, "C")
            f.write("\n")

        if resGenCausality == "der" or resGenCausality == "mixed":
            WriteApproxDerFunction(f, "C")
            f.write("\n")

        if len(usedParams) > 0:
            f.write("void\n")
            f.write("GetParameters( PyObject *pyParams, Parameters* params )\n")
            f.write("{\n")
            for p in usedParams:
                f.write(tab + "PyObject *" + p + ' = PyDict_GetItemString(pyParams, "' + p + '");\n')
                f.write(tab + "params->" + p + " = PyFloat_AsDouble(" + p + ");\n")
                f.write("\n")
            f.write("}\n\n")

        if len(resGenState) > 0:
            f.write("void\n")
            f.write("GetState( PyObject *pyState, ResState* state )\n")
            f.write("{\n")
            for sv in resGenState:
                f.write(tab + "PyObject *" + sv + ' = PyDict_GetItemString(pyState, "' + sv + '");\n')
                f.write(tab + "state->" + sv + " = PyFloat_AsDouble(" + sv + ");\n")
                f.write("\n")
            f.write("}\n\n")
        WriteCPythonAPI(f, name)
        WriteSetupBuild(name, external_src)

    f.close()


def WriteSetupBuild(name, external_src):
    """Write setup script for building residual generator."""
    f = open(name + "_setup.py", "w")
    f.write("# Build file for compiling " + name + ".c\n")
    f.write("#\n")
    f.write("# Compile by executing:\n")
    f.write("# > python " + name + "_setup.py build\n")
    f.write("# at a shell prompt\n\n")
    f.write("# File generated " + time.strftime("%c") + "\n\n")
    f.write("from setuptools import setup, Extension\n")
    f.write("import numpy as np\n")
    f.write("\n")
    f.write("incdir = np.get_include()\n")
    f.write("\n")
    f.write("module1 = Extension('" + name + "',\n")
    f.write("                    sources=['" + name + ".cc'")
    for src in external_src:
        f.write(", '" + src + "'")
    f.write("],\n")
    f.write("                    include_dirs=[incdir]")
    if platform.system() != "Windows":
        f.write(",\n                    extra_compile_args=['-Wno-unused-function'])\n")
    else:
        f.write(")\n")
    f.write("\n")
    f.write("setup(name='" + name + "',\n")
    f.write("      version='0.1',\n")
    f.write("      description='Module for residual generator " + name + "',\n")
    f.write("      ext_modules=[module1])\n")
