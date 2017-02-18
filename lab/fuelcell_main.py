# %matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import gaussian_kde
import sys

new_paths = ['../src', '../Misc']
[sys.path.append(d) for d in new_paths if not d in sys.path];
 
from misc import BoxOff
