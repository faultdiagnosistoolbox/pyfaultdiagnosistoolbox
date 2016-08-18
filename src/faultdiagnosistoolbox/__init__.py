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

The toolbox available under a MIT license. The latest
version can always be downloaded from our website at
http://www.fs.isy.liu.se/Software/PyFaultDiagnosisToolbox/ and
links to relevant publications can be found also at our list of
publications http://www.fs.isy.liu.se/Publications.
"""
from . DiagnosisModel import DiagnosisModel, DiffConstraint
from . dmperm import dmperm, srank, MSO, GetDMParts, Mplus

#__all__ = ["DiagnosisModel", "dmperm", "PlotDM", "Matching"]
#__all__ = ["DiagnosisModel", "dmperm", "structureanalysis"]
#__all_  = ["DiagnosisModel"]
