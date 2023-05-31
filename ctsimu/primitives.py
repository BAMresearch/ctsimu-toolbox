# -*- coding: UTF-8 -*-
"""Basic structures for linear algebra and geometry: vectors, matrices, lines and polygons."""

import math
import numpy
import numbers

from .helpers import *

class Matrix:
    """Simple matrix class.

    Attributes
    ----------
    cols : int
        Number of matrix columns.

    rows : int
        Number of matrix rows.

    n_entries : int
        Number of matrix elements. (Computed internally whenever matrix size is set or changed.)

    value : numpy.ndarray
        NumPy array that contains the matrix values.
    """

    def __init__(self, cols:int = None, rows:int = None, numpy_data = None):
        """Initialize matrix by given size or numpy data array.

        If `cols` and `rows` are not `None`, a matrix of the given size
        will be created with all elements being zero. Otherwise,
        the matrix will be set up from the `numpy_data` array if provided.

        Parameters
        ----------
        cols : int, optional
            Number of matrix columns.

        rows : int, optional
            Number of matrix rows.

        numpy_data : numpy.ndarray, optional
            2-dimensional NumPy array. The number of columns and rows is determined from the array.
        """

        if numpy_data is not None:
            self.set_numpy_data_array(numpy_data)
        elif (cols is not None) and (rows is not None):
            self.reset(cols=cols, rows=rows)
        else:
            raise Exception("Matrix initialization failed. Number of rows and columns must be provided, or an array of values.")

    def __str__(self):
        return f"{self.value}"

    def __add__(self, x):
        result = self.get_copy()
        result.add(x)
        return result

    def __sub__(self, x):
        result = self.get_copy()
        result.subtract(x)
        return result

    def __mul__(self, x):
        result = self.get_copy()
        result.multiply(x)
        if isinstance(x, Vector):
            # We need to return a vector
            return Vector(numpy_data=result.value)
        
        return result

    def __truediv__(self, x):
        result = self.get_copy()
        result.divide(x)
        return result

    def __floordiv__(self, x):
        result = self.get_copy()
        result.floor_divide(x)
        return result

    def __radd__(self, x):
        return self.__add__(x)

    def __rsub__(self, x):
        if isinstance(x, numbers.Number):
            result = self.get_copy()
            result.multiply(-1)
            result.add(x)
            return result

    def __rmul__(self, x):
        if isinstance(x, numbers.Number):
            return self.__mul__(x)

    def size(self) -> int:
        """Get the number of matrix elements.

        Returns
        -------
        n_entries : int
            Number of matrix elements.
        """
        return self.n_entries

    def same_size(self, M:'Matrix') -> bool:
        """Check if this matrix has the same size as the given matrix `M`.

        Parameters
        ----------
        M : Matrix

        Returns
        -------
        same_size : bool
            `True` if `M` has the same number of rows and columns as this matrix, otherwise `False`.
        """
        if (self.cols == M.cols) and (self.rows == M.rows):
            return True

        return False

    def reset(self, cols:int, rows:int):
        """Set matrix size to given number of columns and rows, and set all matrix elements to zero.

        Parameters
        ----------
        cols : int
            Number of matrix columns.

        rows : int
            Number of matrix rows.
        """
        self.cols = cols
        self.rows = rows
        self.n_entries = cols*rows
        self.value = numpy.zeros((rows, cols), dtype=numpy.float64)

    def make_identity(self):
        """Make this an identity matrix."""
        self.value = numpy.zeros((self.rows, self.cols), dtype=numpy.float64)
        for i in range(min(self.cols, self.rows)):
            self.value[i][i] = float(1)

    def set(self, col:int, row:int, value:float):
        """Set a matrix element's value.

        Parameters
        ----------
        col : int
            The element's column position.

        row : int
            The element's row position.

        value : float
            The new value of the matrix element.
        """
        self.value[row][col] = value

    def set_numpy_data_array(self, numpy_data:'numpy.ndarray'):
        """Set up the matrix from a given NumPy array.

        Gets the number of columns and rows from the array and
        uses the array for its values.

        Warning: the array content is not copied, only referenced.
        If you change the array afterwards, the matrix content will
        change as well. Use `numpy.copy()` if you want to pass
        a copy of your array.

        Parameters
        ----------
        numpy_data : numpy.ndarray
            NumPy array of the new matrix values.
        """
        self.rows = len(numpy_data)
        self.cols = len(numpy_data[0])
        self.n_entries = self.cols * self.rows
        self.value = numpy_data.astype(float)

    def copy(self, x:'Matrix'):
        """Make this matrix a copy of the given matrix `x`.

        Parameters
        ----------
        x : Matrix
            The matrix whose contents shall be copied.
        """
        self.set_numpy_data_array(numpy_data=numpy.copy(x.value))

    def get(self, col:int, row:int) -> float:
        """Get a matrix element's value.

        Parameters
        ----------
        col : int
            Element's column position.

        row : int
            Element's row position.

        Returns
        -------
        value : float
            Value of the requested matrix element.
        """
        if len(self.value) > row:
            if len(self.value[row]) > col:
                return self.value[row][col]

        raise Exception(f"Matrix.get(col={col}, row={row}): requested index does not exist.")

    def get_copy(self) -> 'Matrix':
        """Get a copy of this matrix object.

        Returns
        -------
        copy : Matrix
            Copy of this matrix object.
        """
        new_values = numpy.copy(self.value)
        return Matrix(numpy_data=new_values)

    def add(self, x):
        """Add a scalar to all matrix elements, or add another compatible matrix of the same size.

        Parameters
        ----------
        x : Matrix or float
            Matrix to be added to this matrix, or scalar to be added to all matrix elements.
        """
        if isinstance(x, Matrix):
            # 'x' is another matrix:
            if self.same_size(x):
                self.value = numpy.add(self.value, x.value)
            else:
                raise Exception("Incompatible matrix sizes. Can only add matrices of same size.")
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.add(self.value, x)
        else:
            raise Exception(f"Matrix addition failed: M+A is only supported when A is of type 'Matrix' or a scalar. The given type is not supported: {type(x)}.")

    def subtract(self, x):
        """Subtract a scalar from all matrix elements, or subtract another compatible matrix.

        Parameters
        ----------
        x : Matrix or float
            Matrix to be subtracted from this matrix, or scalar to be subtracted from all matrix elements.
        """
        if isinstance(x, Matrix):
            # 'x' is another matrix:
            if self.same_size(x):
                self.value = numpy.subtract(self.value, x.value)
            else:
                raise Exception("Incompatible matrix sizes. Can only subtract matrices of same size.")
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.subtract(self.value, x)
        else:
            raise Exception(f"Matrix subtraction failed: M-A is only supported when A is of type 'Matrix' or a scalar. The given type is not supported: {type(x)}.")

    def multiply(self, x):
        """Multiply a scalar to all matrix elements, or perform matrix multiplication with a compatible matrix.

        Parameters
        ----------
        x : Matrix or float
            Compatible matrix to be multiplied with this matrix, or scalar to be multiplied to all matrix elements.
        """
        if isinstance(x, Matrix) or isinstance(x, Vector):
            # 'x' is another matrix or a vector:
            self.value = numpy.matmul(self.value, x.value)
            self.n_entries = self.value.size
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.multiply(self.value, x)
        else:
            raise Exception(f"Matrix multiplication failed: M*A is only supported when A is of type 'Matrix' or 'Vector' or a scalar. The given type is not supported: {type(x)}.")

    def divide(self, x):
        """Divide all matrix elements by a scalar, or perform element-wise division with a matrix of the same size.

        Parameters
        ----------
        x : Matrix or float
            Matrix of same size for element-wise division with this matrix, or scalar to divide all matrix elements.
        """
        if isinstance(x, Matrix):
            # 'x' is another matrix... element-wise division:
            if self.same_size(x):
                self.value = numpy.divide(self.value, x.value)
                self.n_entries = self.value.size
            else:
                raise Exception("Incompatible matrix sizes. Can only run an element-wise division for matrices of same size.")
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.divide(self.value, x)
        else:
            raise Exception(f"Matrix division failed: M/A is only supported when A is of type 'Matrix' or a scalar. The given type is not supported: {type(x)}.")

    def floor_divide(self, x):
        """Floor-divide all matrix elements by a scalar, or perform element-wise division with a matrix of the same size.

        Parameters
        ----------
        x : Matrix or float
            Matrix of same size for element-wise floor-division with this matrix, or scalar to floor-divide all matrix elements.
        """
        if isinstance(x, Matrix):
            # 'x' is another matrix... element-wise division:
            if self.same_size(x):
                self.value = numpy.floor_divide(self.value, x.value)
            else:
                raise Exception("Incompatible matrix sizes. Can only run an element-wise floor-division for matrices of same size.")
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.floor_divide(self.value, x)
        else:
            raise Exception(f"Matrix floor-division failed: M//A is only supported when A is of type 'Matrix' or a scalar. The given type is not supported: {type(x)}.")

    def scale(self, factor:float):
        """Scale matrix by a scalar factor.

        Parameters
        ----------
        factor : float
            Factor that scales all matrix elements.
        """
        self.value = numpy.multiply(self.value, factor)

class Vector:
    """A vector in space, arbitrary number of dimensions.

    Attributes
    ----------
    n_entries : int
        Number of vector elements.

    value : numpy.ndarray
        NumPy array that contains the vector elements.
    """

    def __init__(self, x:float=None, y:float=None, z:float=None, w:float=None, n:int=None, numpy_data=None):
        """Initialize vector by providing elements `x`, (`y`, `z`, `w`) and possibly the number of elements `n`, or a NumPy data array.

        Parameters
        ----------
        x : float, optional
            Value for first vector element.

        y : float, optional
            Value for second vector element.

        z : float, optional
            Value for third vector element.

        w : float, optional
            Value for fourth vector element.

        n : int, optional
            Number of vector elements. Must be provided if no NumPy data array
            is given and if the number of vector elements cannot be automatically
            determined from the given parameters `x`, `y`, `z` and `w`
            (those which are set to `None` are not considered in the
            determination of the dimension).

        numpy_data : numpy.ndarray, optional
            One-dimensional NumPy data array. Must be provided if number of vector elements `n` is not given.
        """

        # The following member variables are private and should not
        # be accessed from outside. They are invalidated whenever the
        # vector changes.
        self._unit_vector = None  # the unit vector that corresponds to this vector
        self._length = None  # vector's length

        if numpy_data is not None:
            # Initialize from given numpy data array:
            self.set_numpy_data_array(numpy_data)
        else:
            # Or from regular value arguments:
            if n is not None:
                # If the number of vector components is specified:
                self.n_entries = n
            else:
                # Determine number of vector components from given values:
                if x is not None:
                    self.n_entries = 1
                    if y is not None:
                        self.n_entries = 2
                        if z is not None:
                            self.n_entries = 3
                            if w is not None:
                                self.n_entries = 4

            self.reset(n=self.n_entries)
            self.set(x=x, y=y, z=z, w=w)

        self.update()

    def __str__(self):
        return f"{self.value}"

    def __add__(self, x:'Vector'):
        result = self.get_copy()
        result.add(x)
        return result

    def __sub__(self, x:'Vector'):
        result = self.get_copy()
        result.subtract(x)
        return result

    def __mul__(self, x:'Vector'):
        result = self.get_copy()
        result.multiply(x)
        return result

    def __truediv__(self, x:'Vector'):
        result = self.get_copy()
        result.divide(x)
        return result

    def __floordiv__(self, x:'Vector'):
        result = self.get_copy()
        result.floor_divide(x)
        return result

    def __radd__(self, x:'Vector'):
        if isinstance(x, numbers.Number):
            return self.__add__(x)

    def __rsub__(self, x:'Vector'):
        if isinstance(x, numbers.Number):
            result = self.get_copy()
            result.multiply(-1)
            result.add(x)
            return result

    def __rmul__(self, x:'Vector'):
        if isinstance(x, numbers.Number):
            return self.__mul__(x)

    def size(self) -> int:
        """Get the number of vector elements.

        Returns
        -------
        n_entries : int
            Number of vector elements.
        """
        return self.n_entries

    def reset(self, n:int):
        """Set vector to given size and initialize all values to zero.

        Parameters
        ----------
        n : int
            Number of vector elements.
        """
        self.n_entries = n
        self.value = numpy.zeros(n, dtype=numpy.float64)

    def same_size(self, x:'Vector') -> bool:
        """Check if this vector has the same size (i.e., number of vector elements, not length!) as the given vector `x`.

        Returns
        -------
        same_dim : bool
            `True` if number of vector elements matches, `False` otherwise.
        """

        if self.n_entries == x.n_entries:
            return True

        return False

    def update(self):
        """Called when vector is changed.

        Invalidates private member variables for length and unit vector,
        so they will be re-calculated the next time their public getter
        functions are called.
        """
        self._unit_vector = None
        self._length = None

    def x(self) -> float:
        """Get first vector element.

        Returns
        -------
        x : float
            First vector element.
        """
        return self.value[0]

    def y(self) -> float:
        """Get second vector element.

        Returns
        -------
        y : float
            Second vector element.
        """
        return self.value[1]

    def z(self) -> float:
        """Get third vector element.

        Returns
        -------
        z : float

            Third vector element.
        """
        return self.value[2]

    def w(self) -> float:
        """Get fourth vector element.

        Returns
        -------
        w : float

            Fourth vector element.
        """
        return self.value[3]

    def get(self, i:int) -> float:
        """Get value at vector index `i`. Indices start at `0`.

        Parameters
        ----------
        i : int
            Element index of value.

        Returns
        -------
        value : float
            Element value at index position `i`.
        """
        return self.value[i]

    def get_copy(self) -> 'Vector':
        """Get copy of this `Vector` object.

        Returns
        -------
        v : Vector
            Copy of this vector object.
        """
        new_values = numpy.copy(self.value)
        return Vector(numpy_data=new_values)

    def set_x(self, value:float):
        """Set vector's x value (i.e., value of first vector element).

        Parameters
        ----------
        value : float
            Value for the first vector element.
        """
        self.value[0] = float(value)
        self.update()

    def set_y(self, value:float):
        """Set vector's y value (i.e., value of second vector element).

        Parameters
        ----------
        value : float
            Value for the second vector element.
        """
        self.value[1] = float(value)
        self.update()

    def set_z(self, value:float):
        """Set vector's z value (i.e., value of third vector element).

        Parameters
        ----------
        value : float
            Value for the third vector element.
        """
        self.value[2] = float(value)
        self.update()

    def set_w(self, value:float):
        """Set vector's w value (i.e., value of fourth vector element).

        Parameters
        ----------
        value : float
            Value for the fourth vector element.
        """
        self.value[3] = float(value)
        self.update()

    def set_x_y(self, x:float=0, y:float=0):
        """Set x and y component (relevant for 2D computations).

        Parameters
        ----------
        x : float
            Value for the first vector element.

        y : float
            Value for the second vector element.
        """
        self.value[0] = float(x)
        self.value[1] = float(y)
        self.update()

    def set(self, x:float=0, y:float=0, z:float=0, w:float=None):
        """ Set all vector components.

        Non-existing vector components (such as `w` for a 3D vector) should be set to `None`.

        Parameters
        ----------
        x : float
            Value for the first vector element.

        y : float, optional
            Value for the second vector element.

        z : float, optional
            Value for the third vector element.

        w : float, optional
            Value for the fourth vector element.
        """
        if self.n_entries > 0 and x is not None:
            self.value[0] = float(x)

        if self.n_entries > 1 and y is not None:
            self.value[1] = float(y)

        if self.n_entries > 2 and z is not None:
            self.value[2] = float(z)

        if self.n_entries > 3 and w is not None:
            self.value[3] = float(w)

        self.update()

    def set_value_for_index(self, i:int, value:float):
        """Set value for element at vector index `i`.

        Parameters
        ----------
        i : int
            Index of the vector element.

        value : float
            New value for the vector element.
        """
        self.value[i] = float(value)
        self.update()

    def set_numpy_data_array(self, numpy_data):
        """Provide a one-dimensional NumPy array for the vector data.
        The new vector size will be the same as the size of the NumPy array.

        Parameters
        ----------
        numpy_data : numpy.ndarray
            NumPy array for the new vector values.
        """
        self.n_entries = len(numpy_data)
        self.value = numpy_data.astype(float)
        self.update()

    def min(self) -> float:
        """Minimum value of all vector components."""
        return self.value.min()

    def max(self) -> float:
        """Maximum value of all vector components."""
        return self.value.min()

    def absmin(self) -> float:
        """Minimum value of absolute of all vector components."""
        return numpy.absolute(self.value).min()

    def absmax(self) -> float:
        """Maximum value of absolute of all vector components."""
        return numpy.absolute(self.value).max()

    def absmin_nonzero(self) -> float:
        """Minimum non-zero value of absolute of all vector components."""
        nonzeros = self.value[numpy.nonzero(self.value)]
        if len(nonzeros) > 0:
            return numpy.absolute(nonzeros).min()
        else:
            return None

    def absmax_nonzero(self) -> float:
        """Maximum non-zero value of absolute of all vector components."""
        nonzeros = self.value[numpy.nonzero(self.value)]
        if len(nonzeros) > 0:
            return numpy.absolute(nonzeros).max()
        else:
            return None

    def make_crystal_vector(self):
        """Scale this vector into a crystal direction vector,
        using Miller indices."""
        absmin_number = self.absmin_nonzero()
        if absmin_number is not None:
            if absmin_number > 0:
                self.scale(1.0 / absmin_number)

    def crystal_direction(self) -> 'Vector':
        """Create a vector using Miller indices."""
        miller = self.get_copy()
        miller.make_crystal_vector()

        return miller

    def length(self) -> float:
        """Get the length of the vector.

        Returns
        -------
        length : float
            Length of the vector.
        """
        if self._length is None:
            self._length = numpy.linalg.norm(self.value)

        return self._length

    def angle(self, x:'Vector') -> float:
        """Calculate angle between this vector and the given vector `x`.

        Parameters
        ----------
        x : Vector
            Second vector for angle calculation.

        Returns
        -------
        angle : float
            Angle (in radians) between this vector and given vector `x`.
        """
        dotProd = self.dot(x)
        l1 = self.length()
        l2 = x.length()
        lp = l1 * l2

        if lp > 0:
            cs = dotProd / lp
            alpha = 0

            # Avoid out-of-domain due to rounding errors:
            if cs >= 1.0:
                alpha = 0
            elif cs <= -1.0:
                alpha = math.pi
            else:
                alpha = math.acos(cs)

            return alpha
        else:
            return 0

    def make_unit_vector(self):
        """ Normalize vector length to 1. """
        vector_length = self.length()
        if vector_length != 0:
            if vector_length != 1.0:
                self.value /= vector_length
                self.update()
        else:
            raise Exception("Unit vector: a zero length vector cannot be converted into a unit vector.")

    def unit_vector(self) -> 'Vector':
        """ Get a unit vector that points in the same direction as this vector.

        Returns
        -------
        unit_vector : Vector
            Unit vector for this vector.

        Note
        ----
        If this vector changes, the returned unit vector object will be invalidated. To store the unit vector permanently, get a copy of the unit vector (using `get_copy()`) so it won't change:

        ```python
        unitv = v.unit_vector().get_copy()
        ```
        """
        if self._unit_vector is None:
            self._unit_vector = self.get_copy()
            self._unit_vector.make_unit_vector()

        return self._unit_vector

    def add(self, x):
        """Add a scalar to all vector elements, or add another compatible vector.

        Parameters
        ----------
        x : Vector or float
            Vector of same size to be added (element-wise) or scalar to be added to all vector elements.
        """
        if isinstance(x, Vector):
            # 'x' is another vector:
            if self.same_size(x):
                self.value = numpy.add(self.value, x.value)
            else:
                raise Exception(f"Incompatible vector sizes. Can only add vectors of same size. Vectors: {self} ({self.n_entries}) and {x} ({x.n_entries})")
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.add(self.value, x)
        else:
            raise Exception(f"Vector addition failed: M+A is only supported when A is of type 'Vector' or a scalar. The given type is not supported: {type(x)}.")

    def subtract(self, x):
        """Subtract a scalar from all vector elements, or subtract another compatible vector.

        Parameters
        ----------
        x : Vector or float
            Vector of same size to be subtracted (element-wise) or scalar to be subtracted from all vector elements.
        """
        if isinstance(x, Vector):
            # 'x' is another vector:
            if self.same_size(x):
                self.value = numpy.subtract(self.value, x.value)
            else:
                raise Exception("Incompatible vector sizes. Can only subtract vectors of same size.")
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.subtract(self.value, x)
        else:
            raise Exception(f"Vector subtraction failed: M-A is only supported when A is of type 'Vector' or a scalar. The given type is not supported: {type(x)}.")

    def multiply(self, x):
        """Multiply a scalar to all vector elements, or perform vector multiplication with a compatible vector.

        Parameters
        ----------
        x : Vector or float
            Vector of same size to be multiplied (element-wise) or scalar to be multiplied to all vector elements.
        """
        if isinstance(x, Vector):
            # 'x' is another vector:
            self.value = numpy.matmul(self.value, x.value)
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.multiply(self.value, x)
        else:
            raise Exception(f"Vector multiplication failed: M*A is only supported when A is of type 'Vector' or 'Vector' or a scalar. The given type is not supported: {type(x)}.")

    def divide(self, x):
        """Divide all vector elements by a scalar, or perform element-wise division with a vector of the same size.

        Parameters
        ----------
        x : Vector or float
            Vector of same size to divide this vector (element-wise) or scalar to divide all vector elements.
        """
        if isinstance(x, Vector):
            # 'x' is another vector... element-wise division:
            if self.same_size(x):
                self.value = numpy.divide(self.value, x.value)
            else:
                raise Exception("Incompatible vector sizes. Can only run an element-wise division for vectors of same size.")
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.divide(self.value, x)
        else:
            raise Exception(f"Vector division failed: M/A is only supported when A is of type 'Vector' or a scalar. The given type is not supported: {type(x)}.")

    def floor_divide(self, x):
        """Divide all vector elements by a scalar, or perform element-wise division with a vector of the same size.

        Parameters
        ----------
        x : Vector or float
            Vector of same size to floor-divide this vector (element-wise) or scalar to floor-divide all vector elements.
        """
        if isinstance(x, Vector):
            # 'x' is another vector... element-wise division:
            if self.same_size(x):
                self.value = numpy.floor_divide(self.value, x.value)
            else:
                raise Exception("Incompatible vector sizes. Can only run an element-wise division for vectors of same size.")
        elif isinstance(x, numbers.Number):
            # 'x' is a scalar:
            self.value = numpy.floor_divide(self.value, x)
        else:
            raise Exception(f"Vector division failed: M/A is only supported when A is of type 'Vector' or a scalar. The given type is not supported: {type(x)}.")

    def scale(self, factor:float):
        """Scale vector by a scalar factor.

        Parameters
        ----------
        factor : float
            Factor that scales all vector elements.
        """
        self.value = numpy.multiply(self.value, factor)
        self.update()

    def scaled(self, factor:float) -> 'Vector':
        """ Get a copy of this vector, scaled by the given `factor`.

        Parameters
        ----------
        factor : float
            Factor that scales all vector elements.

        Returns
        -------
        result : Vector
            Scaled copy of this vector.
        """
        result = self.get_copy()
        result.scale(factor)
        return result

    def square(self):
        """ Square all elements of this vector. """
        self.value = numpy.square(self.value)
        self.update()

    def squared(self) -> 'Vector':
        """ Get a squared copy of this vector.

        Returns
        -------
        result : Vector
            Squared copy of this vector.
        """
        result = self.get_copy()
        result.square()
        return result

    def sqrt(self):
        """ Set all elements of this vector to their square roots. """
        self.value = numpy.sqrt(self.value)
        self.update()

    def distance(self, p:'Vector') -> float:
        """ Distance between target points of this and another vector.

        Parameters
        ----------
        p : Vector
            Second vector.

        Returns
        -------
        distance : float
            Distance between the target point of this vector and the target point of the second vector `p`.
        """
        return (self-p).length()

    def dot(self, x:'Vector') -> float:
        """ Calculate vector dot product.

        Parameters
        ----------
        x : Vector
            Another vector to calculate the dot product with this vector.

        Returns
        -------
        dotp : float
            Dot product of this vector with the vector `x`.
        """
        return numpy.dot(self.value, x.value)

    def cross_z(self, x:'Vector') -> float:
        """ Calculate the z component of the 3D cross product.

        Parameters
        ----------
        x : Vector
            The second vector to calculate the cross product. Must contain three elements.

        Returns
        -------
        crossz : float
            z component of 3D cross product of this vector with the second vector `x`.
        """
        return self.x()*x.y() - self.y()*x.x()

    def cross(self, x:'Vector') -> 'Vector':
        """ Calculate 3D vector cross product.

        Parameters
        ----------
        x : Vector
            The second vector to calculate the cross product. Must contain three elements.

        Returns
        -------
        crossp : Vector
            3D cross product of this vector with the second vector `x`.
        """
        cp = numpy.cross(self.value, x.value)
        return Vector(numpy_data=cp)

    def sum(self) -> float:
        """ Get the sum of all vector elements.

        Returns
        -------
        sum : float
            Sum of all vector elements.
        """
        return numpy.sum(self.value)

    def invert(self):
        """ Invert this vector: v to -v. """
        self.value = -self.value
        self.update()

    def inverse(self) -> 'Vector':
        """ Get the inverse of this vector: -v.

        Returns
        -------
        inv : Vector
            Inverse of this vector, i.e. vector of same length pointing in the opposite direction.
        """
        result = self.get_copy()
        result.invert()
        return result

    def rotate(self, axis:'Vector', angle:float):
        """ Rotate vector around given axis by given angle (in rad).

        Parameters
        ----------
        axis : Vector
            Rotation axis.

        angle : float
            Rotation angle (in rad).
        """

        # Implementing a general rotation matrix.
        cs = math.cos(angle)
        sn = math.sin(angle)

        vx = self.x()
        vy = self.y()
        vz = self.z()

        nx = axis.unit_vector().x()
        ny = axis.unit_vector().y()
        nz = axis.unit_vector().z()

        rx = vx*(nx*nx*(1.0-cs) + cs)    + vy*(nx*ny*(1.0-cs) - nz*sn) + vz*(nx*nz*(1.0-cs) + ny*sn)
        ry = vx*(ny*nx*(1.0-cs) + nz*sn) + vy*(ny*ny*(1.0-cs) + cs)    + vz*(ny*nz*(1.0-cs) - nx*sn)
        rz = vx*(nz*nx*(1.0-cs) - ny*sn) + vy*(nz*ny*(1.0-cs) + nx*sn) + vz*(nz*nz*(1.0-cs) + cs)

        self.set(x=rx, y=ry, z=rz)

    def transform(self, M:'Matrix'):
        """Apply the transformation given by matrix `M` to this vector.

        Parameters
        ----------
        M : Matrix
            Transformation matrix.
        """
        result = M*self
        self.n_entries = result.n_entries
        self.value     = result.value
        self.update()

    def to(self, x:'Vector') -> 'Vector':
        """ Get a vector that points from this location to the given location `x`.

        Parameters
        ----------
        x : Vector
            Second location vector.

        Returns
        -------
        pointer : Vector
            Vector that points from this vector's target point to the target point of the given vector `x`.
        """
        return self.connection(self, x)

    @staticmethod
    def connection(p0:'Vector', p1:'Vector') -> 'Vector':
        """ Connection vector between two points (represented by vectors).

        Parameters
        ----------
        p0 : Vector
            Origin point of the connection.

        p1 : Vector
            Target point of the connection.

        Returns
        -------
        connecting_vector : Vector
            Vector pointing from the origin `p0` to the target point `p1`.
        """
        return p1 - p0

class Line2D:
    """A mathematical line in 2D space, with a slope and an offset.

    `y = m*x + n`

    Attributes
    ----------
    m : float
        Slope

    n : float
        Vertical offset (i.e., intersection with y axis)
    """

    def __init__(self, m:float=None, n:float=None):
        """Initialize with specified slope and offset (both optional).

        Parameters
        ----------
        m : float, optional
            Slope

        n : float, optional
            Vertical offset (i.e., intersection with y axis)
        """

        self.m = m
        self.n = n

    def __str__(self):
        return "m={}, n={}".format(self.m, self.n)

    def set(self, m:float, n:float):
        """Set line parameters.

        Parameters
        ----------
        m : float
            Slope

        n : float
            Vertical offset (i.e., intersection with y axis)
        """
        self.m = m
        self.n = n

    def set_from_points(self, p0:'Vector', p1:'Vector'):
        """Set line slope `m` and offset `n` from two given points `p0` and `p1`.
        Both points are assumed to be on the line.

        Parameters
        ----------
        p0 : Vector
            First point, given by 2D vector. For higher-dimensional vectors, only the first two coordinates are considered.

        p1 : Vector
            Second point, given by 2D vector. For higher-dimensional vectors, only the first two coordinates are considered.
        """

        # points are defined by 2D vectors
        x0 = p0.x()
        y0 = p0.y()
        x1 = p1.x()
        y1 = p1.y()

        if x0 != x1:
            self.m = (y1-y0) / (x1-x0)
            self.n = y0 - self.m*x0
        else:
            # Vertical line
            self.m = math.inf
            self.n = x0  # Store x intersection in n if line is vertical

    def intersection(self, v:'Line2D') -> 'Vector':
        """Intersection point with another line.

        Parameters
        ----------
        v : Line2D
            Another line that intersects this line.

        Returns
        -------
        intersection : Vector
            2D vector that contains the coordinates of the intersection point.

        Raises
        ------
        Exception : "Lines are parallel."
            If the lines don't intersect.
        """
        m0 = self.m
        n0 = self.n

        m1 = v.m
        n1 = v.n

        if m0 == m1:
            # Lines are parallel
            raise Exception("Lines are parallel.")

        if m0 != math.inf and m1 != math.inf:
            xs = (n1-n0)/(m0-m1)
            ys = m0*xs + n0
            return Vector(x=xs, y=ys, n=2)
        elif m0 == math.inf and m1 != math.inf:
            xs = n0
            ys = m1*xs + n1
            return Vector(x=xs, y=ys, n=2)
        elif m0 != math.inf and m1 == math.inf:
            xs = n1
            ys = m0*xs + n0
            return Vector(x=xs, y=ys, n=2)


class Polygon:
    """ A general 2D polygon with N points in space.

    Attributes
    ----------
    points : list
        List of points that make up the polygon, each point represented by a `Vector` (at least 2D).

    vertex_order_CCW : bool
        `True` if the vertex order is counter-clockwise (standard) or `False` if clockwise.
    """

    def __init__(self, *points:'Vector'):
        """ Initialize with a series of points.
        They should be defined in counter-clockwise direction. They are objects of class `Vector`.

        If the points are specified in clockwise direction, the parameter `vertex_order_CCW` should be set to `False` manually.

        Parameters
        ----------
        *points : Vector
            An arbitrary number of points given.
        """
        self.points = []
        self._area = None

        # Vertices defined in counter-clockwise order (True) or clockwise order (False).
        self.vertex_order_CCW = True

        self.set(*points)

    def __str__(self):
        s = ""
        for i, p in enumerate(self.points):
            s += "P{i}: ({x}, {y})\n".format(i=i, x=p.x, y=p.y)

        return s

    def make_3D(self, z_component:float):
        """Convert all points in xy-plane from 2D vectors to 3D vectors,
        using the provided z_component.

        Parameters
        ----------
        z_component : float
            The new z component for each point's 3D vector.
        """

        for i, p in enumerate(self.points):
            newPoint = Vector(x=p.x, y=p.y, z=z_component, n=3)
            self.points[i] = newPoint

    def set(self, *points:'Vector'):
        """ Set polygon from an arbitrary number of points.
        The points should be defined in counter-clockwise direction. They are objects of class `Vector`.

        If the points are specified in clockwise direction, the parameter `vertex_order_CCW` should be set to `False` manually.

        Parameters
        ----------
        *points : Vector
            An arbitrary number of points.
        """

        self.points = []
        self.points.extend(points)

        self._area = None

    def area(self) -> float:
        """Get the area enclosed by the polygon.

        Returns
        -------
        area : float
            Area enclosed by the polygon.
        """

        if self._area is None:
            self._calculate_area()

        return self._area

    def _calculate_area(self):
        """Calculate the area enclosed by the polygon.

        The area will be stored as an internal parameter. The function `area()` can be used to get this value.
        """

        self._area = 0

        # Split polygon into triangles and calculate area of each
        # triangle using the trapezoid method.
        if len(self.points) >= 3:
            # Start at first point
            p1 = self.points[0]
            x1 = p1.x()
            y1 = p1.y()

            for i in range(1, len(self.points)-1):
                p2 = self.points[i]
                p3 = self.points[i+1]
                x2 = p2.x()
                y2 = p2.y()
                x3 = p3.x()
                y3 = p3.y()
                self._area += 0.5 * ( (y1+y3)*(x3-x1) + (y2+y3)*(x2-x3) - (y1+y2)*(x2-x1) )

    def get_bounding_box(self) -> tuple[int, int, int, int]:
        """Get the polygon's bounding box values.

        Returns
        -------
        leftmost : float
            Leftmost coordinate.

        upmost : float
            Upmost coordinate.

        rightmost : float
            Rightmost coordinate.

        downmost : float
            Downmost coordinate.
        """

        leftmost  = self.points[0].x()
        rightmost = -1
        upmost    = self.points[0].y()
        downmost  = -1

        for p in self.points:
            if p.x() < leftmost:
                leftmost = math.floor(p.x())
            if p.x() > rightmost:
                rightmost = math.ceil(p.x())

            if p.y() < upmost:
                upmost = math.floor(p.y())
            if p.y() > downmost:
                downmost = math.ceil(p.y())

        return int(leftmost), int(upmost), int(rightmost), int(downmost)

    def is_inside_2D(self, point:'Vector') -> bool:
        """ Check if the given point is inside the polygon or on an edge.
        Only the xy plane is considered (2D projection).

        Parameters
        ----------
        point : Vector
            Point coordinates to check if they are inside the polygon.

        Returns
        -------
        inside : bool
            `True` if inside, `False` otherwise.
        """
        x = point.x()
        y = point.y()

        if len(self.points) >= 3:
            p1 = self.points[0]
            x1 = p1.x()
            y1 = p1.y()

            # Set up sub-triangles and check if point is in any of those:
            for i in range(1, len(self.points)-1):
                p2 = self.points[i]
                p3 = self.points[i+1]
                x2 = p2.x()
                y2 = p2.y()
                x3 = p3.x()
                y3 = p3.y()

                # Calculate the barycentric coordinates of the point with respect to the triangle:
                D = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)

                lambda1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / D
                lambda2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / D
                lambda3 = 1 - lambda1 - lambda2

                #print("Is {} inside triangle?   D: {}, l1: {}, l2: {}, l3: {}".format(point, D, lambda1, lambda2, lambda3))

                if (lambda1>=0 and lambda2>=0 and lambda3>=0):
                    return True

        return False

    def _inside_edge(self, edgePoint0:'Vector', edgePoint1:'Vector', vertexToTest:'Vector') -> bool:
        """ Helper function for clip():
            decide if vertex point is on the "inside" of the clipping edge.
            Inside means "to the left" if vertices are in counter-clockwise direction,
            otherwise "to the right". """

        edge  = edgePoint1 - edgePoint0
        point = vertexToTest - edgePoint0

        cpz = edge.cross_z(point)

        if cpz >= 0:  # on the edge
            return True

        return False

    def clip(self, clipping_polygon:'Polygon') -> 'Polygon':
        """ Clips the polygon using the given clipping polygon.

        Implementation of the Sutherland-Hodgman clipping algorithm.

        The given `clipping_polygon` must be convex.

        Parameters
        ----------
        clipping_polygon : Polygon
            Clipping polygon. The result will be a polygon that only exists within the boundaries of this clipping polygon.

        Returns
        -------
        clipped : Polygon
            Clipped polygon.
        """

        # Make a list of edges (lines) of the clipping polygon:
        outputVertices = self.points # copy.deepcopy(self.points)

        nPoints = len(clipping_polygon.points)
        for i in range(nPoints):
            edgePoint0 = clipping_polygon.points[i]
            edgePoint1 = clipping_polygon.points[int((i+1)%nPoints)]

            edgeLine = Line2D()
            edgeLine.set_from_points(p0=edgePoint0, p1=edgePoint1)

            inputVertices = outputVertices
            outputVertices = []

            #print("## Clip line: {}".format(edgeLine))

            for i in range(len(inputVertices)):
                currentPoint  = inputVertices[i]
                previousPoint = inputVertices[(i+len(inputVertices)-1)%len(inputVertices)]

                #print("  Current Point:  {}".format(currentPoint))
                #print("  Previous Point: {}".format(previousPoint))

                currentLine = Line2D()
                currentLine.set_from_points(p0=currentPoint, p1=previousPoint)

                #print("  -> current Line: {}".format(currentLine))

                if self._inside_edge(edgePoint0, edgePoint1, currentPoint):
                    #print("  Current point is inside clipping polygon.")
                    if not self._inside_edge(edgePoint0, edgePoint1, previousPoint):
                        try:
                            intersectionPoint = currentLine.intersection(edgeLine)
                            #print("  Added intersection to output list.")
                            #print("  -> intersection: {}".format(intersectionPoint))
                            outputVertices.append(intersectionPoint)
                        except:  # Parallel lines -> no intersection
                            pass

                    #print("  Added currentPoint to output list.")
                    outputVertices.append(currentPoint)
                elif self._inside_edge(edgePoint0, edgePoint1, previousPoint):
                    #print("  Current point is not inside clipping polygon, but previous point is.")
                    try:
                        intersectionPoint = currentLine.intersection(edgeLine)
                        #print("  Only added intersection point to output list.")
                        #print("  -> intersection: {}".format(intersectionPoint))
                        outputVertices.append(intersectionPoint)
                    except:
                        pass
                #else:
                    #print("Neither current nor previous point is inside clipping polygon.")


                #print("\n")

        result = Polygon(*outputVertices)
        return result

def rotation_matrix(axis:'Vector', angle:float) -> 'Matrix':
    """A matrix that performs a 3D vector rotation around the
    given `axis` vector by the given `angle` (in rad).

    Note that this is only a rotation matrix; translations
    are not taken into account. This means that the pivot point
    will always be the origin of the object that you rotate, i.e.,
    the rotation axis vector is attached to the object's origin.

    Parameters
    ----------
    axis : Vector
        Rotation axis. Must not be a unit vector.

    angle : float
        Rotation angle.

    Returns
    -------
    R : Matrix
        Rotation matrix.
    """
    R = Matrix(3, 3)
    cs = math.cos(angle)
    sn = math.sin(angle)

    nx = axis.unit_vector().x()
    ny = axis.unit_vector().y()
    nz = axis.unit_vector().z()

    # Row 1:
    R.value[0][0] = nx*nx*(1.0-cs) + cs
    R.value[0][1] = nx*ny*(1.0-cs) - nz*sn
    R.value[0][2] = nx*nz*(1.0-cs) + ny*sn

    # Row 2:
    R.value[1][0] = ny*nx*(1.0-cs) + nz*sn
    R.value[1][1] = ny*ny*(1.0-cs) + cs
    R.value[1][2] = ny*nz*(1.0-cs) - nx*sn

    # Row 3:
    R.value[2][0] = nz*nx*(1.0-cs) - ny*sn
    R.value[2][1] = nz*ny*(1.0-cs) + nx*sn
    R.value[2][2] = nz*nz*(1.0-cs) + cs

    return R
