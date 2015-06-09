import math
from random import random


class ValueNotRepresentable(Exception):
    def __init__(self, value, qubits, *args, **kwargs):
        super(Exception, self).__init__('Can\'t represent value %s'
        ' with %s qubits' % (value, qubits))


def tensor(a, b):
    cols_a, rows_a, cols_b, rows_b = len(a), len(a[0]), len(b), len(b[0])
    return [[a[i//cols_b][j//rows_b]*b[i%cols_b][j%rows_b]
            for j in range(rows_a*rows_b)]
            for i in range(cols_a*cols_b)]


def weighted_choice(choices):
   r = random()
   for i, p in enumerate(choices):
      r -= p
      if r < 0:
        return i


def rnd(v, digits=2):
    '''rounds floats and complex numbers and
    removes imaginary part if it's 0'''
    try:
        if round(v.imag, digits):
            return complex(round(v.real, digits), round(v.imag, digits))
        else:
            return round(v.real, digits)
    except:
        return round(v, digits)


def col(v):
    '''Print Qubit amplitudes (or any vector) in a column form'''
    try:
        print("\n".join(map(lambda x:str(rnd(x)), v._values)))
    except:
        print("\n".join(map(lambda x:str(rnd(x)), v)))


class Q:
    def __init__(self, n=1, value=0, values=None):
        if values:
            self._values = values
            self.n = math.log(len(values), 2)
        else:
            self.n = n
            self._values = [0j]*2**n
            try:
                self._values[value] = 1+0j
            except:
                raise ValueNotRepresentable(value, n)

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, v):
        self._values = v
        self.n = int(math.log(len(v), 2))

    def prob(self, bit=None, digits=2):
        if bit != None:
            pos = 1 << bit
            p = sum((v*v).real for i,v in enumerate(self.values) if i & pos)
            if digits != None:
                return (round(1 - p, digits), round(p, digits))
            else:
                return (1 - p, p)
        else:
            try:
                return [round((v*v).real, digits) for v in self._values]
            except:
                return [(v*v).real for v in self._values]

    def measure(self, bit=None):
        if bit != None:
            result = int(random() > self.prob(bit=bit, digits=None)[0])
            #Collapsing state
            pos = 1 << bit
            d = 0
            for i,v in enumerate(self.values):
                if bool(i & pos) != bool(result): #XOR
                    self.values[i] = 0
                else:
                    d += (v*v).real
            #Normalizing values
            self.values = [v/math.sqrt(d) for v in self.values]
            return result
        else:
            result = weighted_choice(self.prob(digits=None))
            #collapsing to state i
            self._values = [0j]*2**self.n
            self._values[result] = 1+0j
            return result    

    def __getitem__(self, i):
        try:
            return QLabel(self, i, i + 1)
        except:
            return QLabel(self, i[0], i[1])

    def __pow__(self, q):
        return Q(values=tensor([self._values], [q._values])[0])

    def __str__(self):
        return str([rnd(v) for v in self._values])

    def __repr__(self):
        return str([rnd(v) for v in self._values])


class QLabel:
    def __init__(self, q, begin, end):
        self._q = q
        self._begin = begin or 0
        self._end = end or q.n
        self.n = self._end - self._begin

    def reverse(self):
        self._begin, self._end = self._end, self._begin

    @property
    def values(self):
        return self._q.values

    @values.setter
    def values(self, v):
        self._q.values = v

    def __str__(self):
        return str(self._q)

    def __repr__(self):
        return str(self._q)


class Operator:
    def __init__(self, matrix):
        '''matrix is a list of columns'''
        self._matrix = matrix

    @property
    def cols(self):
        return len(self._matrix)

    @property
    def rows(self):
        return len(self._matrix[0])

    def inverse(self):
        n = len(self._matrix)
        m = [[self[i][j].conjugate() for j in range(n)] for i in range(n)]
        return Operator(m)

    def __getitem__(self, i):
        return self._matrix[i]

    def __pow__(self, o):
        return Operator(tensor(self._matrix, o._matrix))

    def __mul__(self, q):
        if isinstance(q, Q):
            v = [sum(self[j][i]*q.values[j]
                    for j in range(self.rows))
                    for i in range(self.cols)]
            return Q(values=v)
        else:
            m = [[sum(a*b for a,b in zip(X_row,Y_col))
                    for Y_col in zip(*q._matrix)]
                    for X_row in self._matrix]
            return Operator(m)

    def __call__(self, q):
        return self*q

    def __repr__(self):
        return repr(self._matrix)


class Identity(Operator):
    def __init__(self, size=1):
        m = [[int(i==j) for j in range(2**size)] for i in range(2**size)]
        super(Identity, self).__init__(m)


class Hadamard(Operator):
    def __init__(self, size=1):
        a = 1/math.sqrt(2) + 0j
        base = Operator([[a,a], [a,-a]])
        m = base
        for i in range(size - 1):
            m = m**base
        super(Hadamard, self).__init__(m._matrix)


class Cnot(Operator):
    def __init__(self):
        m = [[1, 0, 0, 0],
             [0, 1, 0, 0],
             [0, 0, 0, 1],
             [0, 0, 1, 0]]
        super(Cnot, self).__init__(m)


def operation(o, q):
    if isinstance(q, QLabel):
        if q._begin > 0:
            op = Identity(q._begin)
            op = op**o
        else:
            op = o
        if q._end < q._q.n:
            op = op**Identity(q._q.n - q._end)
        q.values = op(q._q).values
        return q._q
    else:
        q.values = o(q).values
        return q


def H(q):
    h = Hadamard(q.n)
    return operation(h, q)


def CNOT(q, qt=None):
    '''qt is the target qubit label, if set, assuming q is the
    control qubit label'''
    if qt:
        #Performing a Cnot without building the operator
        val = list(q.values)
        ctrl_mask, target_mask = 1 << q._begin, 1 << qt._begin
        for i,v in enumerate(q.values):
            if i & ctrl_mask:
                val[i ^ target_mask] = q.values[i]
        q.values = val
        return q._q
    else:
        c = Cnot()
        return operation(c, q)


def measure(q):
    if isinstance(q, QLabel):
        result = 0
        for i in range(q._begin, q._end, 1 if q._begin < q._end else -1):
            result += q._q.measure(i) << (i - q._begin)
        return result
    else:
        return q.measure()
