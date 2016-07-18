import numpy as np
import dmperm as dmperm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.sparse as sp

class DiagnosisModel(object):
    X = np.zeros((0,0),dtype='float64')
    F = np.zeros((0,0),dtype='float64')
    Z = np.zeros((0,0),dtype='float64')
    x = []
    f = []
    z = []
    name = ""

    def __init__(self, modeldef) : 
        self.X = self._ModelStructure( modeldef['rels'], modeldef['x'])
        self.x = modeldef['x']
        self.F = self._ModelStructure( modeldef['rels'], modeldef['f'])
        self.f = modeldef['f']
        self.Z = self._ModelStructure( modeldef['rels'], modeldef['z'])
        self.z = modeldef['z']
        self.e = map(lambda x:"e"+np.str(x+1),np.arange(0,self.ne()))
    
    def ne(self):
        return self.X.shape[0]

    def nx(self):
        return self.X.shape[1]

    def nf(self):
        return self.F.shape[1]

    def nz(self):
        return self.Z.shape[1]
    
    def GetDMParts(self):
        return dmperm.GetDMParts(self.X)
    
    def _ModelStructure(self,rels,x) : 
        ne=len(rels)
        nx=len(x)
    
        X = np.zeros((ne,nx),dtype='int64')
        for k,rel in enumerate(rels):
            if self._IsDifferentialConstraint(rel):
                if (rel[0] in x) and (rel[1] in x):
                    dvIdx = x.index(rel[0])
                    ivIdx = x.index(rel[1])
                    X[k,dvIdx] = 3;
                    X[k,ivIdx] = 2;
            else:        
                X[k, self._RelModelStructure(rel, x)]=1
        return sp.csc_matrix(X)

    def _RelModelStructure(self,rel,x):
        return [xi for xi in range(0,len(x)) if x[xi] in rel]

    def _IsDifferentialConstraint( self, rel ):
        return (len(rel)==3 and rel[2]=="diff");

    def MSO(self):
        return dmperm.MSO(self.X)

    def sprank(self):
        return dmperm.sprank(self.X)
  
    def IsPSO( self, *args ):
        if len(args)>0:
            eq=args[0]
        else:
            eq=np.arange(0,self.X.shape[0])
        
        dm = dmperm.GetDMParts(self.X[eq,:])    
        return (len(dm.Mm.row)==0) and (len(dm.M0)==0)
  
    def PlotDM(self, **options) :

        X = self.X
        labelVars = False
        if options.has_key('verbose'):
            labelVars = options['verbose']
        elif X.shape[0]<40:
            labelVars = True
        if options.has_key('eqclass'):
            eqclass=options['eqclass']
        else:
            eqclass=False
        if options.has_key('fault'):
            fault=options['fault']
        else:
            fault=False

        dm = dmperm.GetDMParts(X)
        rowp = dm.rowp
    
        if eqclass and len(dm.Mp.row)>0:
            # Perform PSO decomposition of M+
            Xp = X[dm.Mp.row,:][:,dm.Mp.col]
            P = dmperm.PSODecomposition(Xp)
     
            # Update PSO decomposition description to correspond to global equation
            # indices
            rowp = dm.Mp.row[P['p']]
            colp = dm.Mp.col[P['q']]

            for idx, ec in enumerate(P['eqclass']):
                P['eqclass'][idx].row = dm.Mp.row[P['eqclass'][idx].row]
                P['eqclass'][idx].col = dm.Mp.col[P['eqclass'][idx].col]

            P['trivclass'] = dm.Mp.row[P['trivclass']]
            P['X0']        = dm.Mp.col[P['X0']]
            P['p']         = dm.Mp.row[P['p']]
            P['q']         = dm.Mp.row[P['q']]
    
            # Update dm.rowp and dm.colp according to PSO decomposition
            prowstart = len(dm.rowp)-len(P['p'])
            dm.rowp[prowstart:] = rowp
    
            pcolstart = len(dm.colp)-len(P['q']);    
            dm.colp[pcolstart:] = colp;    

            plt.spy(X[dm.rowp,:][:,dm.colp]==1,markersize=4, marker="o")
            for idx,val in enumerate(np.argwhere(X[dm.rowp,:][:,dm.colp]==3)):
                plt.text(val[1]-0.06,val[0]+0.15, 'I',color="b")
      
            for idx,val in enumerate(np.argwhere(X[dm.rowp,:][:,dm.colp]==2)):
                plt.text(val[1]-0.06,val[0]+0.15, 'D',color="b")
        
            if labelVars:
                plt.xticks(np.arange(0,X.shape[1]),dm.colp)
                plt.yticks(np.arange(0,X.shape[0]),dm.rowp)

            if len(dm.Mm.row)>0:
                r = len(dm.Mm.row);
                c = len(dm.Mm.col);
                x1 = -0.5;
                x2 = x1+c;
                y1 = -0.5;
                y2 = y1+r;
                plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'b')    

            # Plot exactly determined part
            r = len(dm.Mm.row);
            c = len(dm.Mm.col);
            for hc in dm.M0:
                n = len(hc.row);
                x1 = c-0.5;
                x2 = x1+n;
                y1 = r-0.5;
                y2 = y1+n;
                plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'b')
                r = r+n;
                c = c+n;
    
            # Plot over determined part  
            if len(dm.Mp.row)>0:
                nr = len(dm.Mp.row);
                nc = len(dm.Mp.col);
                x1 = c-0.5;
                x2 = x1+nc;
                y1 = r-0.5;
                y2 = y1+nr;
                plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'b')    

            # Plot equivalence classes in over determined part  
            if eqclass and len(dm.Mp.row)>0:
                r1 = r;
                c1 = c;
                for ec in P['eqclass']:
                    nr = len(ec.row);
                    nc = len(ec.col);
                    x1 = c-0.5;
                    x2 = x1+nc;
                    y1 = r-0.5;
                    y2 = y1+nr;
                    plt.gca().add_patch(mpatches.Rectangle((x1,y1), x2-x1, y2-y1, facecolor='0.7'))
                    r = r+nr;
                    c = c+nc;

            plt.plot([c1-0.5, len(dm.colp)+0.5], [r-0.5, r-0.5], 'k--')
            plt.plot([c-0.5, c-0.5], [r1-0.5, len(dm.rowp)+0.5], 'k--')   

            if fault:
                fPlotRowIdx = map(lambda f: np.argwhere(dm.rowp==f[0])[0][0],np.argwhere(self.F))
                nVars = len(dm.colp)

                for ff in np.unique(fPlotRowIdx):
                    fstr=''
                    faultlist = [self.f[x[1]] for x in np.argwhere(self.F[rowp[ff],:])]
                    for fvidx,fv in enumerate(faultlist):
                        if fvidx==0:
                            fstr = fv
                        else:
                            fstr = "%s, %s" % (fstr,fv)
                        plt.plot([-1,nVars],[ff,ff],'r--')
                        plt.text(nVars-0.1,ff+0.17,fstr, color='r')
    
            # Plot axis ticks
            if labelVars:
                plt.xticks(np.arange(0,X.shape[1]),[self.x[xidx] for xidx in dm.colp])
                plt.yticks(np.arange(0,X.shape[0]),[self.e[eidx] for eidx in dm.rowp])

            # Change plot range
            plt.axis([-0.7,X.shape[1]-0.3,X.shape[0]-0.3,-0.7])

            plt.xlabel('Variables')
            plt.ylabel('Equations')
  
    def PlotModel(self, **options):
        labelVars = False;
        if options.has_key('verbose'):
            labelVars = options['verbose'];
        elif self.nx()+self.nf()+self.nz()<40:
            labelVars = True;
    
        X = self.X.todense()
        F = self.F.todense()
        Z = self.Z.todense()
 
        plt.spy(np.concatenate((X==1,np.zeros(F.shape),np.zeros(Z.shape)),axis=1),markersize=4,marker="o", color="b")
        plt.spy(np.concatenate((np.zeros(X.shape),F,np.zeros(Z.shape)),axis=1),markersize=4,marker="o",color="r")
        plt.spy(np.concatenate((np.zeros(X.shape),np.zeros(F.shape),Z),axis=1),markersize=4,marker="o",color="k")

        for idx,val in enumerate(np.argwhere(X==3)):
            plt.text(val[1]-0.06,val[0]+0.15, 'I',color="b")
      
        for idx,val in enumerate(np.argwhere(X==2)):
            plt.text(val[1]-0.06,val[0]+0.15, 'D',color="b")

        # Plot axis ticks
        if labelVars:
            plt.xticks(np.arange(0,self.nx()+self.nf()+self.nz()),self.x + self.f + self.z)
            plt.yticks(np.arange(0,self.X.shape[0]),self.e)

        # Plot variable set divisors
        plt.plot( [self.nx()-1+0.5, self.nx()-1+.5], [0, self.ne()-1 ], color="k", linestyle="dashed")
        plt.plot( [self.nx()+self.nf()-1+0.5, self.nx()+self.nf()-1+.5], [0, self.ne()-1 ], color="k", linestyle="dashed")
    
        # Change plot range
        #plt.axis([-0.7,self.X.shape[1]-0.3,self.X.shape[0]-0.3, -0.7])

        if len(self.name)>0:
            plt.title(self.name)

        plt.xlabel('Variables')
        plt.ylabel('Equations')


    def IsolabilityAnalysis( self, **options ):
        plot = False
        permute = True
        caus='mixed'
        
        if options.has_key('plot'):
            plot = options['plot'];
        if options.has_key('permute'):
            permute = options['permute'];
        if options.has_key('causality'):
            caus = options['causality'];

        ne = self.X.shape[0]
        nf = len(self.f)

        dm = dmperm.Mplus( self.X )

        # Determine equations for each fault
        feq = map(lambda fi: np.argwhere(self.F[:,fi])[0][0],np.arange(0,nf))

        # Determine non-detectable faults
        ndrows = [x for x in np.arange(0,ne) if x not in dm['row']]
        ndf = [self.f[fi] for fi in np.arange(0,len(self.f)) if feq[fi] in ndrows]        
        df = [self.f[fi] for fi in np.arange(0,len(self.f)) if not feq[fi] in ndrows]

        im = np.ones((nf,nf),dtype=np.int64)
        for fi in np.arange(0,nf):
            # Decouple fi
            fieqs = [x[0] for x in np.argwhere(self.F[:,fi].todense()==0)]
            X = self.X[fieqs,:]
            G = dmperm.Mplus(X, caus)
            fisolrows = [fieqs[ei] for ei in G['row']]
            idx = [fj for fj in np.arange(0,nf) if np.any(self.F[fisolrows,:].todense(),axis=0)[(0,fj)]]
            im[fi,idx]=0;

        if plot:
            if permute:
                dm = dmperm.GetDMParts(im)
                p = dm.rowp
                q = dm.colp
            else:
                p = np.arange(0,nf)
                q = p
            plt.spy(im[p,:][:,q],markersize=8,marker="o", color="b")
            plt.xticks(np.arange(0,nf),[self.f[fi] for fi in p])
            plt.yticks(np.arange(0,nf),[self.f[fi] for fi in p])
            if len(self.name)>0:
                titleString = "Isolability matrix for '%s'" % self.name
            else:
                titleString = "Isolability matrix"
            if caus is 'der':
                titleString = "%s (derivative causality)" % titleString
            elif caus is 'int':
                titleString = "%s (integral causality)" % titleString
            plt.title( titleString )
        return im        

def DiffConstraint(dvar,ivar):
    return [dvar, ivar, "diff"];
            


            
