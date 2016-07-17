import numpy as np
import dmperm as dmperm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def PlotDM(X, **options) :
  labelVars = False
  if options.has_key('verbose'):
    labelVars = options['verbose']
  elif X.shape[0]<30:
    labelVars = True
  if options.has_key('eqclass'):
    eqclass=options['eqclass']
  else:
    eqclass=False

  dm = dmperm.GetDMParts(X)

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

  if eqclass and len(dm.Mp.row)>0:
    # Plot equivalence classes in over determined part  
    r1 = r;
    c1 = c;
    for ec in P['eqclass']:
      nr = len(ec.row);
      nc = len(ec.col);
      x1 = c-0.5;
      x2 = x1+nc;
      y1 = r-0.5;
      y2 = y1+nr;
      #plot( [x1 x1 x2 x2 x1],[y1 y2 y2 y1 y1],'g')
      plt.gca().add_patch(mpatches.Rectangle((x1,y1), x2-x1, y2-y1, facecolor='0.7'))
      r = r+nr;
      c = c+nc;

    plt.plot([c1-0.5, len(dm.colp)+0.5], [r-0.5, r-0.5], 'k--')
    plt.plot([c-0.5, c-0.5], [r1-0.5, len(dm.rowp)+0.5], 'k--')   
    
  # Plot axis ticks
  if labelVars:
    plt.xticks(np.arange(0,X.shape[1]),dm.colp+1)
    plt.yticks(np.arange(0,X.shape[0]),dm.rowp+1)

  # Change plot range
  plt.axis([-0.7,X.shape[1]-0.3,X.shape[0]-0.3,-0.7])

  plt.xlabel('Variables')
  plt.ylabel('Equations')

  
