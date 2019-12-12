# EWPTstudy
Code repository used to study the phenomenology of the EWPT

Thermal fluctuations in the field content of the standard model are expected to lead to electroweak symmetry restoration at high temperatures. The code in this repository lets you play with the phase transition in example models by tuning different parameters, such as coupling constants and the zero-temperature or vacuum expectation value of the Higgs' field.

This project was performed for a final paper in my Advanced Statistical Mechanics class in Fall 2019, see ./paper/ for the text.

The main code can be found in MakePaperPlots.C. This must be compiled against ROOT (available at root.cern.ch) or can be run interactively within a ROOT shell. Parameters are as defined in paper/Paper.pdf. The default values for standard model parameters are their measured values, but one can alter these to see, e.g., how the phase transition would be more strongly first order if the Higgs' boson and top quark masses were lower.
