"""
Formulation described in:
    Jabr, R. A. (2012). Polyhedral formulations and loop elimination constraints 
    for distribution network expansion planning.     
    IEEE Transactions on Power Systems, 28(2), 1888-1897.
"""
module Jabr2012

import ..PowerFlowFormulation 
import ..TimeFormulation

struct BIM <: PowerFlowFormulation end
struct TimeIndependent <: TimeFormulation end

end