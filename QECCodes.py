from Z2 import Z2
from Z2 import Z2_vector
from Z2 import Z2_matrix
from enum import IntEnum
from enum import Enum

class Phase(IntEnum):
    # Note that the fourth roots of unity under multiplication is isomorphic to Z4 under addition (mod 4). 
    # We use the aformentioned isomorphism to decide which element maps where below. 
    PLUS_ONE=0
    PLUS_I=1
    MINUS_ONE=2
    MINUS_I=3

class Pauli(Enum):
    # We use a similar isomorphism here 
    I=Z2_vector([0,0])
    X=Z2_vector([1,0])
    Z=Z2_vector([0,1])
    Y=Z2_vector([1,1])

    def __getitem__(self, index: int):
        return self.value[index]

class StabGens():
    def __init__(self, _gens: list[str]):
        if len(_gens) == 0:
            raise ValueError(f"StabGen Error. _gens must be non-empty.")
        valid_symbols = {"I", "X", "Y", "Z"}
        self.N = len(_gens[0])
        self.gens = []
        for _g in _gens:
            culprit = next((char for char in _g if char not in valid_symbols), None)
            if culprit:
                raise ValueError(f"StabGen Error. Invalid character detected: '{culprit}'")
            if self.N != len(_g):
                raise ValueError(f"StabGen Error. Inhomogeneous data, all stab gens must be the same length.")
            self.gens.append([])
            for _p in _g:
                if _p == 'I':
                    self.gens[-1].append(Pauli.I)
                if _p == 'X':
                    self.gens[-1].append(Pauli.X)
                if _p == 'Z':
                    self.gens[-1].append(Pauli.Z)
                if _p == 'Y':
                    self.gens[-1].append(Pauli.Y)

class StabilizerCode:
    def __init__(self, stab_gens: StabGens, phases: list[Phase]):
        # Make the stabilizer matrix. 
        m = len(stab_gens.gens)
        gens_mat_rows = []
        for gen in stab_gens.gens:
            row_x = []
            row_z = []
            for p in gen:
                row_x.append(p[0])
                row_z.append(p[1])
            gens_mat_rows.append(row_x + row_z)

        self.generator_matrix = Z2_matrix.create(gens_mat_rows)
        self.phases = phases
        self.K = self.generator_matrix.NUM_COLS
        self.N = self.generator_matrix.NUM_ROWS/2
        self.D_BOUND = None # compute this later via heuristic (maybe add support to exact computing of it) 
        self.is_css = None
        self.stabilizers = None 
        
    
