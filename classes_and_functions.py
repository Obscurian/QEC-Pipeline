import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import warnings
import pytest
import qiskit
import qiskit_aer
import array


class Z2(int):
    def __sanity__(self) -> bool:
        return 
    
    def __new__(cls, val: int):
        if val != 0 and val  != 1:
            raise ValueError("Z2 Construction Error: val must be 0 or 1.")
        return super().__new__(cls, int(val) % 2)

    # Arithmetic Methods 
    def __add__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        return Z2(self ^ other)

    def __sub__(self, other):
        self.add(other)

    def __mul__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        return Z2(self & other)

    def __truediv__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        if other.val == 0: 
            raise ZeroDivisionError(f"Z2 Divide by Zero Error: you cannot divide by zero.")
        else: 
            return Z2(self.val)

    __floordiv__ = __truediv__

    def __repr__(self):
        return f"Z2({super().__repr__()})"


class Z2_vector:
    def __init__(self, _list: list[int]):
        if not isinstance(_list, list):
            raise ValueError(f"Z2_vector Construction Error: expected a list[int], got {_list}")

        self.DIM = len(_list)
        self._vec = array.array('i', [0] * self.DIM)

        for i, x in enumerate(_list):
            try:
                self._vec[i] = Z2(x)
            except ValueError as original_error:
                raise ValueError(f"Z2_vector Construction Error: invalid input data.") from original_error

        # if we get here, then the input is a valid list of Z2 objects   
    
    def __add__(self, other):
        if not isinstance(other, Z2_vector): return NotImplemented
        if self.DIM != other.DIM:
            raise ValueError(f"Z2_vector Addition Error: you cannot add or subtract vectors with different dimensions.")
        _new_vec = list()
        for x, y in zip(self._vec, other._vec):
            _new_vec.append(x + y)
        return Z2_vector(_new_vec)

    def __sub__(self, other):
        self.add(other)

    def __mul__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        _new_vec = list()
        for x in self._vec:
            _new_vec.append(x*other)
        return Z2_vector(_new_vec)

    def __truediv__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        return self.__mul__(1/other)
        
class Z2_matrix:
    def __init__(self, _mat: list[list[Z2]]):
        pass
            

class StabilizerCode:
    pass
