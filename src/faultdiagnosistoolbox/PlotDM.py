import dmperm as dmperm
import matplotlib.pyplot as plt

def PlotDM(X) :
  dm = dmperm.GetDMParts(X)

  plt.spy(X[dm.rowp,:][:,dm.colp],markersize=2, marker="o")

  if len(dm.Mm.row)>0:
    r = len(dm.Mm.row);
    c = len(dm.Mm.col);
    x1 = 0.5;
    x2 = x1+c;
    y1 = 0.5;
    y2 = y1+r;
    plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'b')    

  # Plot exactly determined part
  r = 1+len(dm.Mm.row);
  c = 1+len(dm.Mm.col);
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

  plt.xticks(np.arange(0,X.shape[0]),[self.x[i] for i in dm.colp])
  plt.yticks(np.arange(0,X.shape[1]),[self.e[i] for i in dm.rowp])
    
