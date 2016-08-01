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
    if language is 'C' or language is 'Matlab':
         return 'ApproxInt(%s,state.%s,Ts)%s %s %s' % (dv,v,CodeEOL(language),CodeComment(language),enum)
    else:
        return "ApproxInt(%s,state['%s'],Ts)%s %s %s" % (dv,v,CodeEOL(language),CodeComment(language),enum)    
    
def CodeApproxDer(v,enum,language):
    if language is 'C' or language is 'Matlab':
        return 'ApproxDiff(%s,state.%s,Ts)%s %s %s' % (v,v,CodeEOL(language),CodeComment(language),enum)
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
        pass
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

def SeqResGen(model, Gamma, resEq, name, diffres='int', language='Python', batch=False, external=[]):
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
        pass

