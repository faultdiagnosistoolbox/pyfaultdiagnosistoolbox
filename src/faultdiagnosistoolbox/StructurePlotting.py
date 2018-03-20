import numpy as np
import faultdiagnosistoolbox.dmperm as dmperm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def PlotModel(model, **options):
    labelVars = False;
    if 'verbose' in options:
        labelVars = options['verbose'];
    elif model.nx()+model.nf()+model.nz()<40:
        labelVars = True;

    X = model.X
    F = model.F
    Z = model.Z

    plt.spy(np.concatenate((X==1,np.zeros(F.shape),np.zeros(Z.shape)),axis=1),markersize=4,marker="o", color="b")
    plt.spy(np.concatenate((np.zeros(X.shape),F,np.zeros(Z.shape)),axis=1),markersize=4,marker="o",color="r")
    plt.spy(np.concatenate((np.zeros(X.shape),np.zeros(F.shape),Z),axis=1),markersize=4,marker="o",color="k")

    if model.nx()<50:
        fSize = 12
    elif model.nx() > 50 and model.nx() < 75:
        fSize = 10
    else:
        fSize = 8

    for idx,val in enumerate(np.argwhere(X==2)):
        plt.text(val[1],val[0], 'I', color="b", fontsize=fSize, horizontalalignment="center", verticalalignment="center")
    for idx,val in enumerate(np.argwhere(X==3)):
        plt.text(val[1],val[0], 'D', color="b", fontsize=fSize, horizontalalignment="center", verticalalignment="center")

    # Plot axis ticks
    if labelVars:
        plt.xticks(np.arange(0,model.nx()+model.nf()+model.nz()),model.x + model.f + model.z, rotation='vertical')
        bottomMargin = 0.1 + (np.max(list(map(lambda v:len(v),model.x)))-3)*0.1/6 # ad-hoc expression
        plt.subplots_adjust(bottom=np.max([0.1,bottomMargin]))
        plt.yticks(np.arange(0,model.X.shape[0]),model.e)

    # Plot variable set divisors
    plt.plot( [model.nx()-1+0.5, model.nx()-1+.5], [0, model.ne()-1 ], color="k", linestyle="dashed")
    plt.plot( [model.nx()+model.nf()-1+0.5, model.nx()+model.nf()-1+.5], [0, model.ne()-1 ], color="k", linestyle="dashed")

    # Change plot range
    #plt.axis([-0.7,model.X.shape[1]-0.3,model.X.shape[0]-0.3, -0.7])

    if len(model.name)>0:
        plt.title(model.name)

    plt.gca().xaxis.tick_bottom()
    plt.xlabel('Variables')
    plt.ylabel('Equations')

def PlotMatching( model, Gamma, **options):
    p = np.array([],dtype=np.int64)
    q = np.array([],dtype=np.int64)
    for g in Gamma.matching:
        p = np.concatenate((p,g.row))
        q = np.concatenate((q,g.col))
    p = np.flipud(p)
    q = np.flipud(q)
    Xm = model.X[p,:][:,q]

    # Determine if axis should be labeled
    labelVars = False;
    if 'verbose' in options:
        labelVars = options['verbose'];
    elif len(q)<40:
        labelVars = True;

    # Plot structure
    if model.nx()<50:
        fSize = 12
    elif model.nx() > 50 and model.nx() < 75:
        fSize = 10
    else:
        fSize = 8

    plt.spy(Xm==1,markersize=4,marker="o", color="b")
    for idx,val in enumerate(np.argwhere(Xm==3)):
        plt.text(val[1],val[0], 'D', color="b", fontsize=fSize, horizontalalignment="center", verticalalignment="center")
        for idx,val in enumerate(np.argwhere(Xm==2)):
            plt.text(val[1],val[0], 'I', color="b", fontsize=fSize, horizontalalignment="center", verticalalignment="center")

    # Plot axis ticks
    if labelVars:
        bottomMargin = 0.1 + (np.max(list(map(lambda v:len(v),[model.x[xi] for xi in q])))-3)*0.1/6 # ad-hoc expression
        plt.xticks(np.arange(0,len(q)), [model.x[xi] for xi in q],rotation='vertical')
        plt.subplots_adjust(bottom=np.max([0.1,bottomMargin]))
        plt.yticks(np.arange(0,len(p)), [model.e[ei] for ei in p])
    else:
        plt.xticks(np.arange(0,len(q)))
        plt.yticks(np.arange(0,len(p)))

    # Draw lines for Hall components
    pos = len(q)-1 # Lower right corner of spy-plot
    for gi in Gamma.matching:
        n = len(gi.row)
        x1 = pos+0.5
        x2 = pos - n + 0.5
        y1 = pos + 0.5
        y2 = pos - n + 0.5
        plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'k')
        pos = pos - n

    plt.axis([-1, len(q), len(q), -1])
    plt.gca().xaxis.tick_bottom()
    plt.xlabel('Variables')
    plt.ylabel('Equations')


def PlotDM(model, **options) :

    X = model.X
    labelVars = False
    if 'verbose' in options:
        labelVars = options['verbose']
    elif X.shape[0]<40:
        labelVars = True
    if 'eqclass' in options:
        eqclass=options['eqclass']
    else:
        eqclass=False
    if 'fault' in options:
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

    if model.nx()<50:
        fSize = 12
    elif model.nx() > 50 and model.nx() < 75:
        fSize = 10
    else:
        fSize = 8
    plt.spy(X[dm.rowp,:][:,dm.colp]==1,markersize=4, marker="o")
    for idx,val in enumerate(np.argwhere(X[dm.rowp,:][:,dm.colp]==3)):
        plt.text(val[1],val[0]+0.15, 'D', color="b", fontsize=fSize, horizontalalignment="center", verticalalignment="center")

    for idx,val in enumerate(np.argwhere(X[dm.rowp,:][:,dm.colp]==2)):
        plt.text(val[1],val[0]+0.15, 'I', color="b", fontsize=fSize, horizontalalignment="center", verticalalignment="center")

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
        fPlotRowIdx = list(map(lambda f: np.argwhere(dm.rowp==f[0])[0][0],np.argwhere(model.F)))
        nVars = len(dm.colp)
        for ff in np.unique(fPlotRowIdx):
            fstr=''
            faultlist = [model.f[x] for x in np.argwhere(model.F[dm.rowp[ff],:])[0]]
            for fvidx,fv in enumerate(faultlist):
                if fvidx==0:
                    fstr = fv
                else:
                    fstr = "%s, %s" % (fstr,fv)
                plt.plot([-1,nVars],[ff,ff],'r--')
                plt.text(nVars-0.1,ff+0.17,fstr, color='r')

    # Plot axis ticks
    if labelVars:
        plt.xticks(np.arange(0,X.shape[1]),[model.x[xidx] for xidx in dm.colp], rotation='vertical')

        bottomMargin = 0.1 + (np.max(list(map(lambda v:len(v),[model.x[xidx] for xidx in dm.colp])))-3)*0.1/6 # ad-hoc expression
        plt.subplots_adjust(bottom=np.max([0.1,bottomMargin]))
        plt.yticks(np.arange(0,model.X.shape[0]),model.e)


        plt.yticks(np.arange(0,X.shape[0]),[model.e[eidx] for eidx in dm.rowp])

    # Change plot range
    plt.axis([-0.7,X.shape[1]-0.3,X.shape[0]-0.3,-0.7])

    plt.gca().xaxis.tick_bottom()
    plt.xlabel('Variables')
    plt.ylabel('Equations')
