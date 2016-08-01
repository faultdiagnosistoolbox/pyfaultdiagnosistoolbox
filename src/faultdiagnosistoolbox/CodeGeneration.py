import numpy as np
import dmperm as dmperm
import sympy as sym
import sys
import faultdiagnosistoolbox as fdt
from sympy.printing import ccode, octave_code
import time

def UsedVars(rels, variables):
    uv = set([])
    for expr in rels:
        if not fdt.IsDifferentialConstraint(expr):
            uv |= set([str(v) for v in expr.atoms() if str(v) in variables])
    return [v for v in variables if v in uv]

def ExprToCode( expr, language):
    if language is 'C':
        genCode = ccode(expr)
    elif language is 'Matlab':
        genCode = octave_code(expr)
    elif language is 'Python':
        genCode = str(expr)
    else:
        print "Unknown language"
    return genCode

def CodeEOL(language):
    if language is 'C' or language is 'Matlab':
        return ';'
    elif language is 'Python':
        return ''
    return ''

def CodeComment(language):
    if language is 'C':
        return '//'
    elif language is 'Matlab':
        return '%'
    elif language is 'Python':
        return '#'
    return '#'

def DVar(diffConstraint):
    return diffConstraint[0]
def IVar(diffConstraint):
    return diffConstraint[1]

def SeqResGenCausality( Gamma, e, diffres ):
    if not fdt.IsDifferentialConstraint(e):
        return Gamma.matchType
    elif diffres is 'int':
        # Treat differential residual constraint in integral causality
        if Gamma.matchType is 'der':
            return 'mixed'
        elif Gamma.matchType is 'int' or Gamma.matchType is 'mixed':
            return Gamma.type
        elif Gamma.matchType is 'algebraic':
            return 'int'
    else:
        # Treat differential residual constraint in derivative causality
        if Gamma.matchType is 'int':
            return 'mixed'
        elif Gamma.matchType is 'der' or Gamma.matchType is 'mixed':
            return Gamma.matchType
        elif Gamma.matchType is 'algebraic':
            return 'der'

def AlgebraicHallComponent(model, g, language):
    fzSub = zip(model.f, np.zeros(len(model.f),dtype=np.int64))
    eqs = map(lambda eq: eq.subs(fzSub), model.syme[g.row])    
    hcVars = list(np.array(model.x)[g.col])
    sol = sym.solve(eqs,hcVars)
    resGen = []
    for varSol in sol:
        if len(hcVars)==1:
            genCode = "%s = %s%s %s %s" %(varSol,ExprToCode(sol[varSol], language), CodeEOL(language),
                                         CodeComment(language),model.e[g.row[0]])
        else:
            genCode = "%s = %s%s" %(varSol,ExprToCode(sol[varSol], language), CodeEOL(language))
        resGen.append(genCode)
    return resGen

def CodeApproxInt(v,dv,enum,language):
    if language is 'Matlab':
         return 'ApproxInt(%s,state.%s,Ts)%s %s %s' % (dv,v,CodeEOL(language),CodeComment(language),enum)
    elif language is 'C':
         return 'ApproxInt(%s,state->%s,Ts)%s %s %s' % (dv,v,CodeEOL(language),CodeComment(language),enum)

    else:
        return "ApproxInt(%s,state['%s'],Ts)%s %s %s" % (dv,v,CodeEOL(language),CodeComment(language),enum)    
    
def CodeApproxDer(v,enum,language):
    if language is 'Matlab':
        return 'ApproxDiff(%s,state.%s,Ts)%s %s %s' % (v,v,CodeEOL(language),CodeComment(language),enum)
    elif language is 'C':
        return 'ApproxDiff(%s,state->%s,Ts)%s %s %s' % (v,v,CodeEOL(language),CodeComment(language),enum)

    else:
        return "ApproxDiff(%s,state['%s'],Ts)%s %s %s" % (v,v,CodeEOL(language),CodeComment(language),enum)

def IntegralHallComponent(model, g, language):
    fzSub = zip(model.f, np.zeros(len(model.f),dtype=np.int64))
    resGen = []
    integ = []
    iState = []
    for e,v,enum in zip(model.syme[g.row],np.array(model.x)[g.col],np.array(model.e)[g.row]):
        if not fdt.IsDifferentialConstraint(e):
            sol=sym.solve(e.subs(fzSub),v)
            genCode = "%s = %s%s %s %s" %(v,ExprToCode(sol[0], language), CodeEOL(language),CodeComment(language),enum)
            resGen.append(genCode)
        else:
            dv = DVar(e)
            iv = IVar(e)
            genCode = v + " = " + CodeApproxInt(v,dv,enum,language)
            integ.append(genCode)
            iState.append(v)
    return (resGen,integ,iState)

def MixedHallComponent(model, g, language):
    fzSub = zip(model.f, np.zeros(len(model.f),dtype=np.int64))
    resGen = []
    iState = []
    dState = []
    integ = []
    for e,v,enum in zip(model.syme[g.row],np.array(model.x)[g.col],np.array(model.e)[g.row]):
        if not fdt.IsDifferentialConstraint(e):
            sol=sym.solve(e.subs(fzSub),v)
            genCode = "%s = %s%s %s %s" %(v,ExprToCode(sol[0], language), CodeEOL(language),CodeComment(language),enum)
            resGen.append(genCode)
        elif v is IVar(e):
            dv = DVar(e)
            iv = IVar(e)
            genCode = v + " = " + CodeApproxInt(v,dv,enum,language)
            integ.append(genCode)
            iState.append(v)
        else: # v is DVar(e)
            genCode = v + " = " + CodeApproxDer(IVar(e),enum,language)
            resGenM0.append(genCode)
            dState.append(Ivar(e))

    return (resGen,integ,iState,dState)

def GenerateResidualEquations( model, resEq, diffres, language):
    if language is 'C':
        resvar = 'r[0]'
    else:
        resvar = 'r'
    fzSub = zip(model.f, np.zeros(len(model.f),dtype=np.int64))
    e = model.syme[resEq].subs(fzSub)
    
    if not fdt.IsDifferentialConstraint(e):
        resExpr = e.lhs - e.rhs
        genCode = ["%s = %s%s %s %s" % (resvar,ExprToCode(resExpr,language), CodeEOL(language), 
                                       CodeComment(language), np.array(model.e)[resEq])]
        iState = []
        dState = []
        integ = []
    else: # IsDifferentialConstraint(e)
        if diffres is 'der':
            iv = IVar(e)
            dv = DVar(e)
            genCode = ['%s = %s-%s' % (resvar,dv,CodeApproxDer(iv,np.array(model.e)[resEq],language))]
            iState = []
            dState = [dv]
            integ = []
        else: # diffres is 'int'
            iv = IVar(e)
            dv = DVar(e)
            if language is 'Python':
                genCode = ["%s = %s-state['%s']%s %s %s" % (resvar,v,v,CodeEOL(language),
                                                         CodeComment(language), np.array(model.e)[resEq])]
            else:
                genCode = ['%s = %s-state.%s%s %s %s' % (resvar,v,v,CodeEOL(language),
                                                         CodeComment(language), np.array(model.e)[resEq])]
            iState = [iv];
            dState = []
            integ = [iv + " = " + CodeApproxInt(v,dv,np.array(model.e)[resEq],language)]
    return (np.array(genCode), iState, dState, integ)

def GenerateExactlyDetermined( model, Gamma, language):
    resGenM0 = np.array([])
    integ = np.array([])
    iState = np.array([])
    dState = np.array([])
    
    for g in Gamma.matching:
        sys.stdout.write('.')
        sys.stdout.flush()
        if g.matchType is 'algebraic':
            codeGen = AlgebraicHallComponent(model, g, language)
            resGenM0 = np.concatenate((resGenM0, codeGen))
        elif g.matchType is 'int':
            codeGen,gInteg,giState = IntegralHallComponent(model, g, language)
            iState = np.concatenate((iState, giState))
            resGenM0 = np.concatenate((resGenM0, codeGen))
            integ = np.concatenate((integ, gInteg))
        elif g.matchType is 'der':
            dc = model.syme[g.row[0]];
            codeGen = DVar(dc) + " = " + CodeApproxDer(IVar(dc),model.e[g.row[0]],language)
            resGenM0 = np.concatenate((resGenM0, [codeGen]))
            dState = np.concatenate((dState, [IVar(dc)]))
        elif g.matchType is 'mixed':
            codeGen,gInteg,giState,gdState = IntegralHallComponent(model, g, language)
            resGenM0 = np.concatenate((resGenM0, codeGen))    
            iState = np.concatenate((iState, giState))
            dState = np.concatenate((dState, gdState))
            integ = np.concatenate((integ, gInteg))
    sys.stdout.write('\n')
    return (resGenM0, iState, dState, integ)

def WriteApproxIntFunction(f,language):
    if language is 'C':
        f.write('double\n')
        f.write('ApproxInt(double dx, double x0, double Ts)\n')
        f.write('{\n')
        f.write('  return x0 + Ts*dx;\n')
        f.write('}\n')
    elif language is 'Python':
        f.write('    def ApproxInt(dx,x0,Ts):\n')
        f.write('        return x0 + Ts*dx\n')
    elif language is 'Matlab':
        pass

def WriteApproxDerFunction(f,language):
    if language is 'C':
        f.write('double\nApproxDiff(double x, double xold, double Ts)\n')
        f.write('{\n')
        f.write('  return (x-xold)/Ts;\n')
        f.write('}\n')
    elif language is 'Python':
        f.write('    def ApproxDiff(x,xold,Ts):\n')
        f.write('        return (x-xold)/Ts\n')
    elif language is 'Matlab':
        f.write('function dx=ApproxDiff(x,xold,Ts)\n')
        f.write('  if length(xold)==1\n')
        f.write('    dx = (x-xold)/Ts;\n')
        f.write('  end\n')
        f.write('end\n')    
    
def WriteResGenCorePython(f, name, model, resGen, state, integ, resGenEqs):
    # Function definition
    f.write('    def ' + name + '_core(z, state, params, Ts):\n')

    tab = '        '
    # Parameters
    usedParams = UsedVars(model.syme[resGenEqs], model.parameters)
    if len(usedParams)>0:
        f.write(tab + '# Parameters\n')
        for p in usedParams:
            f.write(tab + p + " = params['" + p + "']\n")
        f.write('\n')
            
    # Known variables
    usedKnownSignals = UsedVars(model.syme[resGenEqs], model.z)
    f.write(tab + '# Known signals\n')
    for zIdx,zv in enumerate(usedKnownSignals):
        f.write(tab + zv + ' = z[' + str(zIdx) + ']\n')
    f.write('\n') 
            
    # Initialize state variables
    if len(state)>0:
        f.write(tab + '# Initialize state variables\n')
        for sv in state:
            f.write(tab + sv + " = state['" + sv + "']\n")
        f.write('\n')            
            
    # Residual generator body
    f.write(tab + '# Residual generator body\n')
    for pyExpr in resGen:
        f.write(tab + pyExpr + '\n')
    f.write('\n')
        
    # Update integrator variables
    if len(integ)>0:
        f.write(tab + '# Update integrator variables\n')
        for iExpr in integ:
            f.write(tab + iExpr + '\n')
        f.write('\n')            

    # Update state variables
    if len(state)>0:
        f.write(tab + '# Update state variables\n')
        for sv in state:
            f.write(tab + "state['" + sv + "'] = " + sv + '\n')
        f.write('\n')

    # Return value
    f.write(tab + 'return (r, state)\n')
    
def WriteResGenPython( model, resGen, state, integ, name, batch, resGenCausality, resGenEqs ):
    f=open(name + ".py", 'w')
    tab = '    '
    if not batch:
        f.write('def ' + name + "(z,state,params,Ts):\n")
        f.write(tab + '""" ' + name.upper() + " Sequential residual generator for model '" + model.name + "'\n")
        f.write(tab + 'Causality: ' + resGenCausality + '\n')
        f.write('\n')
        
        fSens = np.array(model.f)[np.any(np.array(model.F[resGenEqs,:].todense()),axis=0)]
        f.write(tab + 'Structurally sensitive to faults: ')
        for fault in fSens[:-1]:
            f.write(fault + ', ')
        f.write(fSens[-1] + '\n')
        f.write('\n')
        f.write(tab + 'Example of basic usage:\n')
        f.write(tab + 'Let z be the observations matrix, each column corresponding to a known signal and Ts the sampling time,\n')
        f.write(tab + 'then the residual generator can be simulated by:\n')
        f.write('\n')
        f.write(tab + "r = np.zeros(N) # N number of data points\n")
        f.write(tab + "state = {")
        for s in state[:-1]:
            f.write("'" + s + "': " + s + '_0, ')
        f.write("'" + state[-1] + "': " + state[-1] + '_0}\n')
        f.write(tab + 'for k,zk in enumerate(z):\n')
        f.write(tab + tab + 'r[k], state = ' + name + '( zk, state, params, Ts )\n')
        f.write('\n')
        f.write(tab + 'State is a dictionary with the keys: ')
        for s in state[:-1]:
            f.write(s + ', ')
        f.write(state[-1] + '\n')
        f.write('\n')
        f.write(tab+'File generated ' + time.strftime('%c') + '\n')
        f.write(tab + '"""') 
        f.write('\n')
        
        if resGenCausality is 'int' or resGenCausality is 'mixed':
            WriteApproxIntFunction(f,'Python')
            f.write('\n')

        if resGenCausality is 'der' or resGenCausality is 'mixed':
            WriteApproxDerFunction(f,'Python')
            f.write('\n')
    
        WriteResGenCorePython(f,name, model,resGen,state, integ,resGenEqs)
        
        f.write('\n')
        f.write('    return ' + name + '_core(z, state, params, Ts)\n')
    else: # batch
        f.write('import numpy as np\n');
        f.write('from copy import copy\n');
        f.write('def ' + name + "(z,state,params,Ts):\n")
        f.write(tab + '""" ' + name.upper() + " Sequential residual generator for model '" + model.name + "'\n")
        f.write(tab + 'Causality: ' + resGenCausality + '\n')
        f.write('\n')
        
        fSens = np.array(model.f)[np.any(np.array(model.F[resGenEqs,:].todense()),axis=0)]
        f.write(tab + 'Structurally sensitive to faults: ')
        for fault in fSens[:-1]:
            f.write(fault + ', ')
        f.write(fSens[-1] + '\n')
        f.write('\n')
        f.write(tab + 'Example of basic usage:\n')
        f.write(tab + 'Let z be the observations matrix, each column corresponding to a known signal and Ts the sampling time,\n')
        f.write(tab + 'then the residual generator can be simulated by:\n')
        f.write('\n')
        f.write(tab + "state = {")
        for s in state[:-1]:
            f.write("'" + s + "': " + s + '_0, ')
        f.write("'" + state[-1] + "': " + state[-1] + '_0}\n')
        f.write(tab + 'r = ' + name + '( zk, state, params, Ts )\n')
        f.write('\n')
        f.write(tab + 'State is a dictionary with the keys: ')
        for s in state[:-1]:
            f.write(s + ', ')
        f.write(state[-1] + '\n')
        f.write('\n')
        f.write(tab+'File generated ' + time.strftime('%c') + '\n')
        f.write(tab + '"""') 
        f.write('\n')
        
        if resGenCausality is 'int' or resGenCausality is 'mixed':
            WriteApproxIntFunction(f,'Python')
            f.write('\n')

        if resGenCausality is 'der' or resGenCausality is 'mixed':
            WriteApproxDerFunction(f,'Python')
            f.write('\n')
    
        WriteResGenCorePython(f,name, model,resGen,state, integ,resGenEqs)
        
        f.write('\n')
        f.write(tab + 'N = z.shape[0] # Number of datapoints\n')
        f.write(tab + 'r = np.zeros(N)\n');
        f.write(tab + 'state = copy(state)\n')
        f.write(tab + 'for k,zk in enumerate(z):\n');
        f.write(tab + tab + 'r[k], state = ' + name + '_core(zk, state, params, Ts)\n')
        f.write('\n')
        f.write(tab + 'return r\n');
    f.close()

def SeqResGen(model, Gamma, resEq, name, diffres='int', language='Python', batch=False, api='Python', external=[]):
    if not model.type is 'Symbolic':
        print "Code generation only possible for symbolic models"
        return []
    
    sys.stdout.write("Generating residual generator " + name + " (" + language + ', ')
    if batch:
        sys.stdout.write("batch")
    else:
        sys.stdout.write("no batch")
    sys.stdout.write(")\n")
    
    sys.stdout.write("  Generating code for the exactly determined part: ")
    m0Code, m0iState, m0dState, m0integ = GenerateExactlyDetermined( model, Gamma, language)

    print "  Generating code for the residual equations"
    resCode, resiState, resdState, resinteg = GenerateResidualEquations( model, resEq, diffres, language)

    # Collect code, integrators, and states
    resGenCode = np.concatenate((m0Code,[' '],resCode))
    resGeniState = np.concatenate((m0iState,resiState))
    resGendState = np.concatenate((m0dState,resdState))
    resGenState  = np.concatenate((resGeniState,resGendState))
    resGenInteg = np.concatenate((m0integ,resinteg))
        
    print "  Writing residual generator file"
    resGenCausality = SeqResGenCausality( Gamma, model.syme[resEq], diffres )
    resGenEqs = np.array([],dtype=np.int64)
    for g in Gamma.matching:
        resGenEqs = np.concatenate((resGenEqs,g.row))
    resGenEqs = np.concatenate((resGenEqs,[resEq]))
    if language is 'Python':
        WriteResGenPython( model, resGenCode, resGenState, resGenInteg, name, batch, resGenCausality, resGenEqs )
        print 'File ' + name + '.py generated.'
    elif language is 'Matlab':
        pass
    elif language is 'C':
        if api is 'Python':
            WriteResGenCPython( model, resGenCode, resGenState, resGenInteg, name, batch, resGenCausality, resGenEqs )
            print 'Files ' + name + '.c and ' + name + '_setup.py generated'
            print 'Compile by running: python ' + name + '_setup.py build'

def WriteResGenCoreC(f,name, model, resGen, state, integ, resGenEqs):
    usedParams = UsedVars(model.syme[resGenEqs], model.parameters)
    tab = '  '

    # Function definition
    f.write('void\n')
    f.write(name + '_core(double* r, const double *z, ')
    if len(resGen)>0:
        f.write('ResState* state, ')
    if len(usedParams)>0:
        f.write('const Parameters *params, ')
    f.write('double Ts)\n')
    f.write('{\n')
    # Parameters
    if len(usedParams)>0:
        f.write(tab + '// Parameters\n')
        for p in usedParams:
            f.write(tab + 'double ' + p + " = params->" + p + ";\n")
        f.write('\n')

    usedVars = UsedVars(model.syme[resGenEqs], model.x)
    if len(usedVars)>0:
        f.write(tab + '// Residual generator variables\n')
        for p in usedVars:
            f.write(tab + 'double ' + p + ';\n')
        f.write('\n')
        
    # Initialize state variables
    if len(state)>0:
        f.write(tab + '// Initialize state variables\n')
        for sv in state:
            f.write(tab + sv + " = state->" + sv + ";\n")
        f.write('\n')            
            
    # Known variables
    usedKnownSignals = UsedVars(model.syme[resGenEqs], model.z)
    f.write(tab + '// Known signals\n')
    for zIdx,zv in enumerate(usedKnownSignals):
        f.write(tab + 'double ' + zv + ' = z[' + str(zIdx) + '];\n')
    f.write('\n') 

    # Residual generator body
    f.write(tab + '// Residual generator body\n')
    for pyExpr in resGen:
        f.write(tab + pyExpr + '\n')
    f.write('\n')
        
    # Update integrator variables
    if len(integ)>0:
        f.write(tab + '// Update integrator variables\n')
        for iExpr in integ:
            f.write(tab + iExpr + '\n')
        f.write('\n')            

    # Update state variables
    if len(state)>0:
        f.write(tab + '// Update state variables\n')
        for sv in state:
            f.write(tab + "state->" + sv + " = " + sv + ';\n')
    f.write('}\n')

def WriteResGenCPython( model, resGenCode, resGenState, resGenInteg, name, batch, resGenCausality, resGenEqs ):
    f=open(name + ".c", 'w')
    tab = '  '
    if not batch:
        f.write('#include <Python.h>\n')
        f.write('#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION\n')
        f.write('#include <numpy/ndarrayobject.h>\n')
        f.write('\n')

        f.write('// ' + name.upper() + " Sequential residual generator for model '" + model.name + "'\n")
        f.write('// Causality: ' + resGenCausality + '\n')
        f.write('//\n')
        
        fSens = np.array(model.f)[np.any(np.array(model.F[resGenEqs,:].todense()),axis=0)]
        f.write('// Structurally sensitive to faults: ')
        for fault in fSens[:-1]:
            f.write(fault + ', ')
        f.write(fSens[-1] + '\n')
        f.write('//\n')
        f.write('// Example of basic usage:\n')
        f.write('// Let z be the observations matrix, each column corresponding to a known signal and Ts the sampling time,\n')
        f.write('// then the residual generator can be simulated by:\n')
        f.write('//\n')
        f.write("// r = np.zeros(N) # N number of data points\n")
        f.write("// state = {")
        for s in resGenState[:-1]:
            f.write("'" + s + "': " + s + '_0, ')
        f.write("'" + resGenState[-1] + "': " + resGenState[-1] + '_0}\n')
        f.write('// for k,zk in enumerate(z):\n')
        f.write('//     r[k], state = ' + name + '( zk, state, params, Ts )\n')
        f.write('//\n')
        f.write('// State is a dictionary with the keys: ')
        for s in resGenState[:-1]:
            f.write(s + ', ')
        f.write(resGenState[-1] + '\n')
        f.write('\n')
        f.write('// File generated ' + time.strftime('%c') + '\n')
        f.write('\n')
        
        # Declare Parameters struct
        usedParams = UsedVars(model.syme[resGenEqs], model.parameters)
        if len(usedParams)>0:
            f.write( 'typedef struct {\n')
            for p in usedParams:
                f.write( tab + 'double ' + p + ';\n')
            f.write( '} Parameters;\n')
            f.write('\n')
            
        # Declare state struct
        if len(resGenState)>0:
            f.write( 'typedef struct {\n')
            for sv in resGenState:
                f.write( tab + 'double ' + sv + ';\n')
            f.write( '} ResState;\n')
            f.write('\n')

        # Declare function prototypes
        if resGenCausality is 'int' or resGenCausality is 'mixed':
            f.write('double ApproxInt(double dx, double x0, double Ts);\n')
        if resGenCausality is 'der' or resGenCausality is 'mixed':
            f.write('double ApproxDiff(double x, double xold, double Ts);\n')
        if len(usedParams)>0:
            f.write('void GetParameters( PyObject *pyParams, Parameters* params );\n')
        if len(resGenState)>0:
            f.write('void GetState( PyObject *pyStates, ResState* state );\n')
        f.write('void ' + name + '_core( double* r, const double *z, ')
        if len(resGenState)>0:
            f.write('ResState *state, ')
        if len(usedParams)>0:
            f.write('const Parameters *params, ')
        f.write('double Ts );\n')

        f.write('\n')
        f.write('static PyObject*\n')
        f.write( name + '(PyObject *self, PyObject * args)\n')
        f.write('{\n')
        f.write('  PyArrayObject *pyZ;\n')
        f.write('  PyObject *pyState;\n')
        f.write('  PyObject *pyParams;\n')
        f.write('  double Ts;\n')
        f.write('\n')
        f.write('  if (!PyArg_ParseTuple(args, "O!O!O!d", &PyArray_Type, &pyZ, &PyDict_Type, &pyState, &PyDict_Type, &pyParams, &Ts)) {\n')
        f.write('    return NULL;\n')
        f.write('  }\n')
        f.write('\n')
        f.write('  double r;\n')
        f.write('  double *z = (double *)PyArray_DATA( pyZ ); // column first\n')
        f.write('\n')
        if len(resGenState)>0:
            f.write('  ResState state;\n')
            f.write('  GetState(pyState, &state);\n')
            f.write('\n')
        if len(usedParams)>0:
            f.write('  Parameters params;\n')
            f.write('  GetParameters( pyParams, &params );\n')
            f.write('\n')
        f.write('  ' + name + '_core( &r, z, ')
        if len(resGenState)>0:
            f.write('&state, ')
        if len(usedParams)>0:
            f.write('&params, ')
        f.write('Ts );\n')
        f.write('\n')
        f.write('  PyObject *newState = PyDict_New();\n\n')
        for sv in resGenState:
            f.write('  PyDict_SetItemString(newState, "' + sv + '", PyFloat_FromDouble(state.' + sv + '));\n')
        f.write('\n')
        f.write('  return Py_BuildValue("(f,O)", r, newState);\n')
        f.write('}\n')
        
        f.write('\n')
        WriteResGenCoreC(f, name, model, resGenCode, resGenState, resGenInteg, resGenEqs)
        f.write('\n')
        
        if resGenCausality is 'int' or resGenCausality is 'mixed':
            WriteApproxIntFunction(f,'C')
            f.write('\n')

        if resGenCausality is 'der' or resGenCausality is 'mixed':
            WriteApproxDerFunction(f,'C')
            f.write('\n')
                        
        if len(usedParams)>0:
            f.write('void\n')
            f.write('GetParameters( PyObject *pyParams, Parameters* params )\n')
            f.write('{\n')
            for p in usedParams:
                f.write( tab + 'PyObject *' + p + ' = PyDict_GetItemString(pyParams, "' + p + '");\n');
                f.write( tab + 'params->' + p + ' = PyFloat_AsDouble(' + p + ');\n')
                f.write('\n')                
            f.write('}\n\n')
            
        if len(resGenState)>0:
            f.write('void\n')
            f.write('GetState( PyObject *pyState, ResState* state )\n')
            f.write('{\n')
            for sv in resGenState:
                f.write( tab + 'PyObject *' + sv + ' = PyDict_GetItemString(pyState, "' + sv + '");\n');
                f.write( tab + 'state->' + sv + ' = PyFloat_AsDouble(' + sv + ');\n')
                f.write('\n')
            f.write('}\n\n')
        
            f.write('PyMODINIT_FUNC\n')
            f.write('init' + name + '(void)\n')
            f.write('{\n')
            f.write('  static PyMethodDef ' + name + 'Methods[] = {\n')
            f.write('    {"' + name + '",  ' + name + ', METH_VARARGS, "Residual generator"},\n')
            f.write('    {NULL, NULL, 0, NULL}\n')
            f.write('  };\n')
            f.write('  (void) Py_InitModule("' + name + '", ' + name + 'Methods);\n')
            f.write('  import_array();\n')
            f.write('}\n')
            f.write('\n')
            f.write('int\n')
            f.write('main(int argc, char *argv[])\n')
            f.write('{\n')
            f.write('  Py_SetProgramName(argv[0]);\n')
            f.write('  \n')
            f.write('  Py_Initialize();\n')
            f.write('  \n')
            f.write('  init' + name + '();\n')
            f.write('\n')
            f.write('  return 0;\n')
            f.write('}\n')
            
            WriteSetupBuild(name)
    else: # batch
        f.write('#include <Python.h>\n')
        f.write('#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION\n')
        f.write('#include <numpy/ndarrayobject.h>\n')
        f.write('\n')

        f.write('// ' + name.upper() + " Sequential residual generator for model '" + model.name + "'\n")
        f.write('// Causality: ' + resGenCausality + '\n')
        f.write('//\n')
        
        fSens = np.array(model.f)[np.any(np.array(model.F[resGenEqs,:].todense()),axis=0)]
        f.write('// Structurally sensitive to faults: ')
        for fault in fSens[:-1]:
            f.write(fault + ', ')
        f.write(fSens[-1] + '\n')
        f.write('//\n')
        f.write('// Example of basic usage:\n')
        f.write('// Let z be the observations matrix, each column corresponding to a known signal and Ts the sampling time,\n')
        f.write('// then the residual generator can be simulated by:\n')
        f.write('//\n')
        f.write("// state = {")
        for s in resGenState[:-1]:
            f.write("'" + s + "': " + s + '_0, ')
        f.write("'" + resGenState[-1] + "': " + resGenState[-1] + '_0}\n')
        f.write("// r = " + name + '( z, state, params, Ts )\n')
        f.write('//\n')
        f.write('// State is a dictionary with the keys: ')
        for s in resGenState[:-1]:
            f.write(s + ', ')
        f.write(resGenState[-1] + '\n')
        f.write('\n')
        f.write('// File generated ' + time.strftime('%c') + '\n')
        f.write('\n')
        
        # Declare Parameters struct
        usedParams = UsedVars(model.syme[resGenEqs], model.parameters)
        if len(usedParams)>0:
            f.write( 'typedef struct {\n')
            for p in usedParams:
                f.write( tab + 'double ' + p + ';\n')
            f.write( '} Parameters;\n')
            f.write('\n')
            
        # Declare state struct
        if len(resGenState)>0:
            f.write( 'typedef struct {\n')
            for sv in resGenState:
                f.write( tab + 'double ' + sv + ';\n')
            f.write( '} ResState;\n')
            f.write('\n')

        # Declare function prototypes
        if resGenCausality is 'int' or resGenCausality is 'mixed':
            f.write('double ApproxInt(double dx, double x0, double Ts);\n')
        if resGenCausality is 'der' or resGenCausality is 'mixed':
            f.write('double ApproxDiff(double x, double xold, double Ts);\n')
        if len(usedParams)>0:
            f.write('void GetParameters( PyObject *pyParams, Parameters* params );\n')
        if len(resGenState)>0:
            f.write('void GetState( PyObject *pyStates, ResState* state );\n')
        f.write('void ' + name + '_core( double* r, const double *z, ')
        if len(resGenState)>0:
            f.write('ResState *state, ')
        if len(usedParams)>0:
            f.write('const Parameters *params, ')
        f.write('double Ts );\n')

        f.write('\n')
        f.write('static PyObject*\n')
        f.write( name + '(PyObject *self, PyObject * args)\n')
        f.write('{\n')
        f.write('  PyArrayObject *pyZ;\n')
        f.write('  PyObject *pyState;\n')
        f.write('  PyObject *pyParams;\n')
        f.write('  double Ts;\n')
        f.write('\n')
        f.write('  if (!PyArg_ParseTuple(args, "O!O!O!d", &PyArray_Type, &pyZ, &PyDict_Type, &pyState, &PyDict_Type, &pyParams, &Ts)) {\n')
        f.write('    return NULL;\n')
        f.write('  }\n')
        f.write('\n')
        
        f.write('  if( PyArray_NDIM(pyZ) != 2 ) {\n')
        f.write('    return Py_BuildValue("");\n')
        f.write('  }\n')
        f.write('  npy_intp *shape = PyArray_SHAPE( pyZ );\n')
        f.write('  npy_intp N = shape[0];\n')
        f.write('  npy_intp nz = shape[1];\n\n')

        f.write('  PyObject* pyR = PyArray_SimpleNew(1, &N, NPY_FLOAT64);\n')
        f.write('  double *r = (double *)PyArray_DATA((PyArrayObject *)pyR);\n\n')
        f.write('  double *z = (double *)PyArray_DATA( pyZ ); // column first\n')
        f.write('\n')
        if len(resGenState)>0:
            f.write('  ResState state;\n')
            f.write('  GetState(pyState, &state);\n')
            f.write('\n')
        if len(usedParams)>0:
            f.write('  Parameters params;\n')
            f.write('  GetParameters( pyParams, &params );\n')
            f.write('\n')
            
        f.write('  // Main computational loop\n')
        f.write('  for( npy_intp k=0; k < N; k++ ) {\n')
        f.write('    ' + name + '_core( &(r[k]), &(z[k*nz]), ')
        if len(resGenState)>0:
            f.write('&state, ')
        if len(usedParams)>0:
            f.write('&params, ')
        f.write('Ts );\n')
        f.write('}\n')
        f.write('\n')
        f.write('  return Py_BuildValue("O", pyR);\n')
        f.write('}\n')
        
        f.write('\n')
        WriteResGenCoreC(f, name, model, resGenCode, resGenState, resGenInteg, resGenEqs)
        f.write('\n')
        
        if resGenCausality is 'int' or resGenCausality is 'mixed':
            WriteApproxIntFunction(f,'C')
            f.write('\n')

        if resGenCausality is 'der' or resGenCausality is 'mixed':
            WriteApproxDerFunction(f,'C')
            f.write('\n')
                        
        if len(usedParams)>0:
            f.write('void\n')
            f.write('GetParameters( PyObject *pyParams, Parameters* params )\n')
            f.write('{\n')
            for p in usedParams:
                f.write( tab + 'PyObject *' + p + ' = PyDict_GetItemString(pyParams, "' + p + '");\n');
                f.write( tab + 'params->' + p + ' = PyFloat_AsDouble(' + p + ');\n')
                f.write('\n')                
            f.write('}\n\n')
            
        if len(resGenState)>0:
            f.write('void\n')
            f.write('GetState( PyObject *pyState, ResState* state )\n')
            f.write('{\n')
            for sv in resGenState:
                f.write( tab + 'PyObject *' + sv + ' = PyDict_GetItemString(pyState, "' + sv + '");\n');
                f.write( tab + 'state->' + sv + ' = PyFloat_AsDouble(' + sv + ');\n')
                f.write('\n')
            f.write('}\n\n')
        
            f.write('PyMODINIT_FUNC\n')
            f.write('init' + name + '(void)\n')
            f.write('{\n')
            f.write('  static PyMethodDef ' + name + 'Methods[] = {\n')
            f.write('    {"' + name + '",  ' + name + ', METH_VARARGS, "Residual generator"},\n')
            f.write('    {NULL, NULL, 0, NULL}\n')
            f.write('  };\n')
            f.write('  (void) Py_InitModule("' + name + '", ' + name + 'Methods);\n')
            f.write('  import_array();\n')
            f.write('}\n')
            f.write('\n')
            f.write('int\n')
            f.write('main(int argc, char *argv[])\n')
            f.write('{\n')
            f.write('  Py_SetProgramName(argv[0]);\n')
            f.write('  \n')
            f.write('  Py_Initialize();\n')
            f.write('  \n')
            f.write('  init' + name + '();\n')
            f.write('\n')
            f.write('  return 0;\n')
            f.write('}\n')
            
            WriteSetupBuild(name)

    f.close()
    
def WriteSetupBuild(name):
    f = open(name + '_setup.py', 'w')
    f.write('# Build file for compiling ' + name + '.c\n')
    f.write('#\n')
    f.write('# Compile by executing:\n')
    f.write('# > python ' + name + '_setup.py build\n')
    f.write('# at a shell prompt\n\n')
    f.write('# File generated ' + time.strftime('%c') + '\n\n')
    f.write('from distutils.core import setup, Extension\n')
    f.write('import numpy as np\n')
    f.write('\n')
    f.write('incdir = np.get_include()\n')
    f.write('\n')
    f.write("module1 = Extension('" + name + "',\n")
    f.write("                    sources = ['" + name + ".c'],\n")
    f.write("                    include_dirs=[incdir],\n")
    f.write("                    extra_compile_args=['-Wno-unused-function'])\n")
    f.write('\n')
    f.write("setup (name = '" + name + "',\n")
    f.write("       version = '0.1',\n")
    f.write("       description = 'Module for residual generator " + name + "',\n")
    f.write('       ext_modules = [module1]\n')
    f.write(')\n')

