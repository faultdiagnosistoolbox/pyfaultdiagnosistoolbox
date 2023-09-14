# __init__.py
""" Fault Diagnosis Toolbox is a toolbox for analysis and design of
fault diagnosis systems for dynamic systems, primarily described by
differential equations. In particular, the toolbox is focused on
techniques that utilize structural analysis, i.e., methods that
analyze and utilize the model structure. The model structure is the
interconnections of model variables and is often described as a
bi-partite graph or an incidence matrix. Key features of the toolbox
are:

* Defining diagnosis models, using only model structure or full
  symbolic expressions.
* Diagnosability analysis - analyze a given model to
  determine which faults that can be detected and which faults that
  can be isolated
* Model exploration and analysis, e.g., plotting model properties,
  Dulmage-Mendelsohn decomposition, DAE index analysis, ...
* Finding overdetermined sets of equations (MSO sets), which are
  minimal submodels that can be used to design fault detectors
* Sensor placement - determine minimal sets of sensors needed to
  be able to detect and isolate faults
* Code generation (C and Python) for residual generators. 

"""
from . DiagnosisModel import DiagnosisModel, DiffConstraint, IsDifferentialConstraint
from . testselection_rf import RandomForestTestSelection
from . diag_util import IsolabilityMatrix, DiagnosesAndConfusionMatrix, PlotConfusionMatrix
from . MHS import MHS
from . import models
from . dmperm import srank
from . _version import __version__
