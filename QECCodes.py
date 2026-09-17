from Z2 import Z2
from Z2 import Z2_vector
from Z2 import Z2_matrix
from enum import IntEnum

class Phase(IntEnum):
    # Note that the fourth roots of unity under multiplication is isomorphic to Z4 under addition (mod 4). 
    # We use the aformentioned isomorphism to decide which element maps where below. 
    PLUS_ONE=0
    PLUS_I=1
    MINUS_ONE=2
    MINUS_I=3

class StabilizerCode:
    def __init__(self, stab_gens: Z2_matrix, phases: list[Phase]):
        if (stab_gens.M%2 != 0):
            raise ValueError(f"StabilizerCode Error: Matrix must have an even number of columns.")
        self.phases = phases
        self.K = stab_gens.M
        self.N = stab_gens.N/2
        self.D_BOUND = None # compute this later via heuristic (maybe add support to exact computing of it) 
        self.is_css = None
        self.stabilizers = None 
    
