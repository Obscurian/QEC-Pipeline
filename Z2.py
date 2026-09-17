from __future__ import annotations
from enum import IntEnum as Enum


class Z2(int):
    def __new__(cls, val: int | Z2):
        return super().__new__(cls, int(val) % 2)

    # Arithmetic Methods 
    def __add__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        return Z2(self ^ other)

    def __sub__(self, other):
        self.__add__(other)

    def __mul__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        return Z2(self & other)

    def __truediv__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        if other.val == 0: 
            raise ZeroDivisionError(f"Z2 Divide by Zero Error: you cannot divide by zero.")
        else: 
            return Z2(self)

    __floordiv__ = __truediv__

class VEC_TYPE(Enum):
    COLUMN = 0
    ROW = 1

class Z2_vector:
    def __init__(self, _list: list[int]):
        if not isinstance(_list, list):
            raise ValueError(f"Z2_vector Construction Error: expected a list[int], got {_list}")

        self.type = VEC_TYPE.COLUMN # assumed to be column by default
        self.DIM = len(_list)
        self._vec = [Z2(x) for x in _list]

    def T(self):
        transposed = Z2_vector(self._vec)
        transposed.type = 1 - transposed.type
        return transposed

    def __getitem__(self, index: int) -> Z2:
        return self._vec[index]

    def __setitem__(self, index: int, value: int | Z2):
        try:
            self._vec[index] = Z2(value)
        except ValueError as original_error:
            raise ValueError("Z2_vector Assignment Error: Invalid input data.") from original_error

    def __add__(self, other):
        if not isinstance(other, Z2_vector): return NotImplemented
        if self.DIM != other.DIM or self.type != other.type:
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

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        return self.__mul__(1/other)

    def __matmul__(self, other) -> Z2 | Z2_matrix:
        if not isinstance(other, Z2_vector): return NotImplemented
        if self.type == other.type: 
            raise ValueError(f"Z2_vector Product Error: Inner/Outer product dimensions dont match.")
        if self.type == VEC_TYPE.COLUMN:
            # Outer product 
            M = [[None for _ in range(self.DIM)] for _ in range(other.DIM)]
            for i, ai in enumerate(self._vec):
                for j, bj in enumerate(other._vec):
                    M[i][j] = ai*bj

            return Z2_matrix.create(M)

        # If we get here, its def an inner product
        s = 0
        for i in range(len(self._vec)):
            s += self._vec[i]*other._vec[i]

        return Z2(s%2)

    def shape(self):
        if self.type == VEC_TYPE.ROW:
            return (1,self.DIM)
        else:
            return (self.DIM,1)

    def __repr__(self):
        if self.type == VEC_TYPE.ROW:
            return f"[ {' '.join(map(str, self._vec))} ]"
        lines = [f"[ {val} ]" for val in self._vec]
        return "\n".join(lines)

        
class Z2_matrix:
    __PRIVATE_KEY = object()
    def __init__(self, _table: list[list[Z2]], _auth_key=None):
        if _auth_key is not self.__PRIVATE_KEY:
            raise RuntimeError(
                "Z2_matrix Privacy Error: Do not use Z2_matrix() directly. "
                "Please use the factory method: Matrix.create(data)"
            )
        self.N = len(_table)
        self.M = len(_table[0])

        for row in _table:
            pass
        self._mat = _table
        
    @classmethod
    def create(cls: Z2_matrix, _table: list[list[Z2]]) -> Z2_matrix | Z2_vector:
        if not _table:
            raise ValueError(f"Z2_matrix Error: Matrix cannot be empty")

        t0 = len(_table[0])
        for t in _table:
            if len(t) != t0:
                raise ValueError(f"Z2_matrix Error: Inhomogenous data (every row must be the same length!.")
        
        if len(_table) == 1:
            data = _table[0]
            return Z2_vector(data).T()

        if len(_table[0]) == 1:
            data = []
            for t in _table:
                data.append(t[0])
            return Z2_vector(data)

        if len(_table) == 1 and len(_table[0]) == 1:
            return Z2(_table[0][0])

        return cls(_table, _auth_key=cls.__PRIVATE_KEY)

    def __getitem__(self, index: int) -> list[Z2]:
        return self._mat[index]
    
    def __add__(self, other):
            if not isinstance(other, Z2_matrix): return NotImplemented
            if self.N != other.N or self.M != other.M:
                raise ValueError(f"Z2_matrix Addition Error: you cannot add or subtract matrices with different dimensions.")
            _new_data = list()
            for row_self, row_other in zip(self._mat, other._mat):
                _new_data.append(list())
                for x, y in zip(row_self,row_other):
                    _new_data[-1].append(x + y)
            return Z2_matrix.create(_new_data)
    
    def __sub__(self, other):
        self.__add__(other)

    def __mul__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        _new_data = list()
        for row in self._mat:
            _new_data.append(list())
            for x in row:
                _new_data[-1].append(x*other)
        return Z2_matrix.create(_new_data)
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, Z2): return NotImplemented
        return self.__mul__(1/other)

    def T(self) -> Z2_matrix | Z2_vector:
        transposed_data = [list(col) for col in zip(*self._mat)]
        return Z2_matrix.create(transposed_data)

    def __matmul__(self, other) -> Z2_vector | Z2_matrix:
        if (not isinstance(other, Z2_vector)) and (not isinstance(other, Z2_matrix)): 
            return NotImplemented
        
        # Case A: Matrix x Vector
        if isinstance(other, Z2_vector):
            if other.type != VEC_TYPE.COLUMN:
                raise ValueError("Z2_matrix Product Error: Cannot right multiply Matrix by a Row Vector.")
            if self.M != other.DIM:
                raise ValueError(f"Z2_matrix Product Error: Dimension mismatch. Matrix columns ({self.M}) must match Vector dimension ({other.DIM}).")
            
            result_vec = []
            for row in self._mat:
                dot_product = Z2(0)
                for x, y in zip(row, other._vec):
                    dot_product += x * y
                result_vec.append(dot_product)

            return Z2_vector(result_vec)
            
        # Case B: Matrix x Matrix
        if self.M != other.N:
            raise ValueError(f"Z2_matrix Product Error: Dimension mismatch. Columns of A ({self.M}) must match Rows of B ({other.N}).")
            
        new_data = [[Z2(0) for _ in range(other.M)] for _ in range(self.N)]
        for i in range(self.N):
            for j in range(other.M):
                dot_product = Z2(0)
                for k in range(self.M):
                    dot_product += self._mat[i][k] * other._mat[k][j]
                new_data[i][j] = dot_product
                
        return Z2_matrix.create(new_data)

    def __rmatmul__(self, other) -> Z2_vector:
        if not isinstance(other, Z2_vector): 
            return NotImplemented
        
        if other.type != VEC_TYPE.ROW:
            raise ValueError("Z2_matrix Product Error: Left operand must be a Row Vector for pre-multiplication.")
        if other.DIM != self.N:
            raise ValueError(f"Z2_matrix Product Error: Dimension mismatch. Vector dimension ({other.DIM}) must match Matrix rows ({self.N}).")
        
        result_vec = []
        # Compute dot product of the row vector with each column of the matrix
        for j in range(self.M):
            dot_product = Z2(0)
            for i in range(self.N):
                dot_product += other._vec[i] * self._mat[i][j]
            result_vec.append(dot_product)
            
        # Return a Z2_vector explicitly forced to be a ROW vector
        return Z2_vector(result_vec).T()

    def shape(self):
        return (self.N, self.M)

    def __repr__(self) -> str:
        lines = []
        for row in self._mat:
            # Join elements with a space, then frame them inside brackets
            row_str = " ".join(map(str, row))
            lines.append(f"[ {row_str} ]")
        return "\n".join(lines)