from __future__ import annotations
from enum import IntEnum as Enum


class Z2(int):
    def __new__(cls, val: int | Z2):
        return super().__new__(cls, int(val) % 2)

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
        '''Creates a column vector with entries from the _list argument'''
        if not isinstance(_list, list):
            raise ValueError(f"Z2_vector Construction Error: expected a list[int], got {_list}")

        self.type = VEC_TYPE.COLUMN # assumed to be column by default
        self.DIM = len(_list)
        self._vec = [Z2(x) for x in _list]

    def T(self):
        '''Returns the transpose of the vector.'''
        transposed = Z2_vector(self._vec)
        transposed.type = 1 - transposed.type
        return transposed

    @staticmethod
    def ones(n: int):
        '''Returns an n by 1 column vector of 1s.'''
        _list = [1 for _ in range(n)]
        return Z2_vector(_list)

    @staticmethod
    def zeroes(n: int):
        '''Returns an n by 1 column vector of 0s.'''
        _list = [0 for _ in range(n)]
        return Z2_vector(_list)

    def __getitem__(self, index: int) -> Z2:
        return self._vec[index]

    def __setitem__(self, index: int, value: int | Z2):
        try:
            self._vec[index] = Z2(value)
        except ValueError as original_error:
            raise ValueError("Z2_vector Assignment Error: Invalid input data.") from original_error

    def __eq__(self, value):
        if isinstance(value, Z2_vector):
            if self._vec == value._vec:
                return self.type == value.type
        return False

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
    def __init__(self, _table: list[list[int|Z2]], _auth_key=None):
        '''Do not use this. This is for private use only. call Z2_matrix.create() instead.'''
        if _auth_key is not self.__PRIVATE_KEY:
            raise RuntimeError(
                "Z2_matrix Privacy Error: Do not use Z2_matrix() directly. "
                "Please use the factory method: Matrix.create(data)"
            )
        self.NUM_ROWS = len(_table)
        self.NUM_COLS = len(_table[0])

        self._mat =  [[Z2(y) for y in x] for x in _table]

    def __eq__(self, value):
        if isinstance(value, Z2_matrix):
            return self._mat == value._mat

    def copy(self):
        return Z2_matrix.create(self._mat)

    @property
    def rref(self):
        '''The reduced row echelon form.'''
        # Only compute the RREF if we haven't done it yet
        if not hasattr(self, '_rref_cache'):
            reduced_matrix = self.__row_reduce__(self.__PRIVATE_KEY) 
            self._rref_cache = reduced_matrix
            
        return self._rref_cache

    @property
    def rank(self):
        '''The rank of the matrix (dim(rowspace) or dim(colspace)).'''
        if not hasattr(self, '_rank_cache'):
            self._rank_cache = len(self.row_space)

        return self._rank_cache
        
    @property
    def nullity(self):
        '''The dimension of the kernel.'''
        if not hasattr(self, '_nul_cache'):
            self._nul_cache = self.NUM_COLS - self.rank # via rank-nullity theorem 

        return self._nul_cache

    @property
    def row_space(self):
        '''List of basis vectors for the row space.'''
        if not hasattr(self, '_rowsp_cache'):
            rows = []
            for row in self.rref:
                if Z2(1) in row:
                    rows.append(Z2_vector(row))

            self._rowsp_cache = rows

        return self._rowsp_cache

    @property
    def column_space(self):
        '''List of basis vectors for the column space.'''
        if not hasattr(self, '_colsp_cache'):
            pivots = []
            for i, row in enumerate(self.rref):
                j_max = 0
                while j_max < self.NUM_COLS and row[j_max] == Z2(0):
                    j_max += 1
                if j_max != self.NUM_COLS:
                    pivots.append(i)

            self._colsp_cache = [self.get_col(i) for i in pivots]

        return self._colsp_cache

    @property
    def kernel(self):
        '''List of basis vectors for the null space.'''
        if not hasattr(self, '_ker_cache'):
            pivot_row_map = {}
            for i, row in enumerate(self.rref):
                j_max = 0
                while j_max < self.NUM_COLS and row[j_max] == Z2(0):
                    j_max += 1
                if j_max != self.NUM_COLS:
                    pivot_row_map[i] = j_max

            frees = [j for j in range(self.NUM_COLS) if j not in pivot_row_map.values()]

            basis = []

            for f in frees:
                _vec = [Z2(0) for _ in range(self.NUM_COLS)]
                _vec[f] = Z2(1)
                for i, j in pivot_row_map.items():
                    _vec[j] = self.rref[i][f]

                basis.append(Z2_vector(_vec))
            
            self._ker_cache = basis

        return self._ker_cache
    
    def __row_reduce__(self, _auth_key=None):
        '''Do not use this. It is a private method that is called implicitly. Algorithm from: https://en.wikipedia.org/wiki/Gaussian_elimination'''
        if _auth_key is not self.__PRIVATE_KEY:
                raise RuntimeError(
                    "Z2_matrix Privacy Error: Do not call __row_reduce__ directly."
                    "It is implicitly called by the construction procedure."
                )
        A = self
        h = 0 # Initialization of the pivot row
        k = 0 # Initialization of the pivot column

        while h < self.NUM_ROWS and k < self.NUM_COLS:
            # Find the k-th pivot 
            i_max = h
            while i_max < self.NUM_ROWS and A[i_max][h] == Z2(0):
                i_max += 1

            if i_max == self.NUM_ROWS:
                # No pivot in this column, pass to next column
                k += 1
            else:
                if h != i_max:
                    A = Z2_matrix.get_ero_swap(self.NUM_ROWS, h, i_max) @ A
                # Do for all rows below pivot:
                for i in range(self.NUM_ROWS):
                    if i != h and A[i][k] == Z2(1):
                        A = Z2_matrix.get_ero_rowadd(self.NUM_ROWS, i, h) @ A
                h += 1
                k += 1

        return A

    @classmethod
    def create(cls: Z2_matrix, _table: list[list[int|Z2]]) -> Z2_matrix | Z2_vector:
        '''Creates a matrix with values specified in the _table argument.'''
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

    @staticmethod
    def ones(n: int, m: int):
        '''Returns an m by n matrix of 1s.'''
        _data = [[1 for _ in range(n)] for _ in range(m)]
        return Z2_matrix.create(_data)
    
    @staticmethod
    def zeroes(m: int, n: int):
        '''Returns an m by n matrix of 0s.'''
        _data = [[0 for _ in range(n)] for _ in range(m)]
        return Z2_matrix.create(_data)

    @staticmethod
    def id(m: int):
        '''Returns an m by m identity matrix.'''
        zeroes = Z2_matrix.zeroes(m,m)
        for i in range(m):
            zeroes[i][i] = Z2(1)
        return zeroes
    
    @staticmethod
    def get_ero_swap(m: int, i: int, j: int):
        '''Provides an m x m elementary matrix that swaps rows i and j.'''
        if i >= m or j >= m:
            raise ValueError(f"Z2_matrix Error: neither row ({i}) nor column ({j}) index can be greater than m ({m}).")
        E = Z2_matrix.id(m)
        E._mat[i], E._mat[j] = E._mat[j], E._mat[i]
        return E

    @staticmethod
    def get_ero_rowadd(m: int, target: int, source: int):
        '''Provides an m x m elementary matrix that replaces R_target with R_target + R_source.'''
        if target >= m or source >= m:
            raise ValueError(f"Z2_matrix Error: neither target ({target}) nor source ({source}) can be greater than m ({m}).")
        if target == source:
            raise ValueError(f"Z2_matrix Error: target ({target}) cannot equal source ({source}).")
        E = Z2_matrix.id(m)
        # Add a 1 at the intersection of the target row and the source column
        E._mat[target][source] = Z2(1)
        return E

    def get_row(self, i: int):
        '''Returns row i'''
        return Z2_vector(self._mat[i]).T()

    def get_col(self, i: int):
        '''Returns column i'''
        _list = []
        for row in self._mat:
            _list.append(row[i])
        return Z2_vector(_list)
    
    def __getitem__(self, key) -> list[Z2]:
        if isinstance(key, tuple):
            row, col = key
            return self.data[row][col]
        if isinstance(key, int):
            return self._mat[key]
        else: raise ValueError(f"Z2_matrix Error. Invalid index {key}. Expected int or tuple.")
    
    def __add__(self, other):
            if not isinstance(other, Z2_matrix): return NotImplemented
            if self.NUM_ROWS != other.NUM_ROWS or self.NUM_COLS != other.NUM_COLS:
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
            if self.NUM_COLS != other.DIM:
                raise ValueError(f"Z2_matrix Product Error: Dimension mismatch. Matrix columns ({self.NUM_COLS}) must match Vector dimension ({other.DIM}).")
            
            result_vec = []
            for row in self._mat:
                dot_product = Z2(0)
                for x, y in zip(row, other._vec):
                    dot_product += x * y
                result_vec.append(dot_product)

            return Z2_vector(result_vec)
            
        # Case B: Matrix x Matrix
        if self.NUM_COLS != other.NUM_ROWS:
            raise ValueError(f"Z2_matrix Product Error: Dimension mismatch. Columns of A ({self.NUM_COLS}) must match Rows of B ({other.NUM_ROWS}).")
            
        new_data = [[Z2(0) for _ in range(other.NUM_COLS)] for _ in range(self.NUM_ROWS)]
        for i in range(self.NUM_ROWS):
            for j in range(other.NUM_COLS):
                dot_product = Z2(0)
                for k in range(self.NUM_COLS):
                    dot_product += self._mat[i][k] * other._mat[k][j]
                new_data[i][j] = dot_product
                
        return Z2_matrix.create(new_data)

    def __rmatmul__(self, other) -> Z2_vector:
        if not isinstance(other, Z2_vector): 
            return NotImplemented
        
        if other.type != VEC_TYPE.ROW:
            raise ValueError("Z2_matrix Product Error: Left operand must be a Row Vector for pre-multiplication.")
        if other.DIM != self.NUM_ROWS:
            raise ValueError(f"Z2_matrix Product Error: Dimension mismatch. Vector dimension ({other.DIM}) must match Matrix rows ({self.NUM_ROWS}).")
        
        result_vec = []
        # Compute dot product of the row vector with each column of the matrix
        for j in range(self.NUM_COLS):
            dot_product = Z2(0)
            for i in range(self.NUM_ROWS):
                dot_product += other._vec[i] * self._mat[i][j]
            result_vec.append(dot_product)
            
        # Return a Z2_vector explicitly forced to be a ROW vector
        return Z2_vector(result_vec).T()

    def shape(self):
        return (self.NUM_ROWS, self.NUM_COLS)

    def __repr__(self) -> str:
        lines = []
        for row in self._mat:
            # Join elements with a space, then frame them inside brackets
            row_str = " ".join(map(str, row))
            lines.append(f"[ {row_str} ]")
        return "\n".join(lines)