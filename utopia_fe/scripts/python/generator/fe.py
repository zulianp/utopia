from symengine import *
# import sympy as sympy
import sympy

# from sympy import *
from time import perf_counter
from sympy.codegen.ast import AddAugmentedAssignment
from sympy.codegen.ast import Assignment
from sympy.printing.c import C99CodePrinter
from sympy.utilities.codegen import CCodeGen
from sympy.utilities.codegen import codegen

# import sys
# from sage.all import *

import os
import rich
import sys

from rich.syntax import Syntax
console = rich.get_console()

# class SymbolicEngine:

class SymPyEngine:
    def symbols(self, args):
        return sympy.symbols(args)

    def simplify(self, expr):
        return sympy.simplify(expr)

    def diff(self, f, x):
        return sympy.diff(f, x)

    def point_symbols(self, num_points, dim):
        prefix = [ 'x', 'y', 'z', 't', 'a', 'b', 'c', 'd']
        ret = []
        for i in range(0, num_points):
            p = []
            for j in range(0, dim):
                p.append(sympy.symbols(f'p{prefix[j]}[{i}]'))
                
            ret.append(sympy.Matrix(dim, 1, p))
        return ret

    def point2(self, x, y):
        return sympy.Matrix(2, 1, [x, y])

    def point3(self, x, y, z):
        return sympy.Matrix(3, 1, [x, y, z])

    def point4(self, x, y, z, t):
        return sympy.Matrix(4, 1, [x, y, z, t])

    def array(self, args):
        return sympy.Matrix(len(args), 1, args)

    def zeros(self, rows, cols):
        return sympy.zeros(rows, cols)

    def matrix(self, rows, cols, data):
        return sympy.Matrix(rows, cols, data)

    def det3(self, mat):
          return mat[0, 0] * mat[1, 1] * mat[2, 2] + mat[0, 1] * mat[1, 2] * mat[2, 0] + mat[0, 2] * mat[1, 0] * mat[2, 1] - mat[0, 0] * mat[1, 2] * mat[2, 1] - mat[0, 1] * mat[1, 0] * mat[2, 2] - mat[0, 2] * mat[1, 1] * mat[2, 0]

    def inv3(self, mat):
        mat_inv = self.zeros(3, 3)

        det = self.det3(mat)
        # det = self.simplify(det)
        mat_inv[0, 0] = (mat[1, 1] * mat[2, 2] - mat[1, 2] * mat[2, 1]) / det
        mat_inv[0, 1] = (mat[0, 2] * mat[2, 1] - mat[0, 1] * mat[2, 2]) / det
        mat_inv[0, 2] = (mat[0, 1] * mat[1, 2] - mat[0, 2] * mat[1, 1]) / det
        mat_inv[1, 0] = (mat[1, 2] * mat[2, 0] - mat[1, 0] * mat[2, 2]) / det
        mat_inv[1, 1] = (mat[0, 0] * mat[2, 2] - mat[0, 2] * mat[2, 0]) / det
        mat_inv[1, 2] = (mat[0, 2] * mat[1, 0] - mat[0, 0] * mat[1, 2]) / det
        mat_inv[2, 0] = (mat[1, 0] * mat[2, 1] - mat[1, 1] * mat[2, 0]) / det
        mat_inv[2, 1] = (mat[0, 1] * mat[2, 0] - mat[0, 0] * mat[2, 1]) / det
        mat_inv[2, 2] = (mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]) / det

        return mat_inv
        # return self.simplify(mat_inv)

    def inv4(self, m):
        mat_inv = self.zeros(4, 4)

        mat_inv[0, 0] = m[1, 1] * m[2, 2] * m[3, 3] - m[1, 1] * m[2, 3] * m[3, 2] - m[2, 1] * m[1, 2] * m[3, 3] + m[2, 1] * m[1, 3] * m[3, 2] + m[3, 1] * m[1, 2] * m[2, 3] - m[3, 1] * m[1, 3] * m[2, 2]

        mat_inv[1, 0] = -m[1, 0] * m[2, 2] * m[3, 3] + m[1, 0] * m[2, 3] * m[3, 2] + m[2, 0] * m[1, 2] * m[3, 3] - m[2, 0] * m[1, 3] * m[3, 2] - m[3, 0] * m[1, 2] * m[2, 3] + m[3, 0] * m[1, 3] * m[2, 2]

        mat_inv[2, 0] = m[1, 0] * m[2, 1] * m[3, 3] - m[1, 0] * m[2, 3] * m[3, 1] - m[2, 0] * m[1, 1] * m[3, 3] + m[2, 0] * m[1, 3] * m[3, 1] + m[3, 0] * m[1, 1] * m[2, 3] - m[3, 0] * m[1, 3] * m[2, 1]

        mat_inv[3, 0] = -m[1, 0] * m[2, 1] * m[3, 2] + m[1, 0] * m[2, 2] * m[3, 1] + m[2, 0] * m[1, 1] * m[3, 2] - m[2, 0] * m[1, 2] * m[3, 1] - m[3, 0] * m[1, 1] * m[2, 2] + m[3, 0] * m[1, 2] * m[2, 1]

        mat_inv[0, 1] = -m[0, 1] * m[2, 2] * m[3, 3] + m[0, 1] * m[2, 3] * m[3, 2] + m[2, 1] * m[0, 2] * m[3, 3] - m[2, 1] * m[0, 3] * m[3, 2] - m[3, 1] * m[0, 2] * m[2, 3] + m[3, 1] * m[0, 3] * m[2, 2]

        mat_inv[1, 1] = m[0, 0] * m[2, 2] * m[3, 3] - m[0, 0] * m[2, 3] * m[3, 2] - m[2, 0] * m[0, 2] * m[3, 3] + m[2, 0] * m[0, 3] * m[3, 2] + m[3, 0] * m[0, 2] * m[2, 3] - m[3, 0] * m[0, 3] * m[2, 2]

        mat_inv[2, 1] = -m[0, 0] * m[2, 1] * m[3, 3] + m[0, 0] * m[2, 3] * m[3, 1] + m[2, 0] * m[0, 1] * m[3, 3] - m[2, 0] * m[0, 3] * m[3, 1] - m[3, 0] * m[0, 1] * m[2, 3] + m[3, 0] * m[0, 3] * m[2, 1]

        mat_inv[3, 1] = m[0, 0] * m[2, 1] * m[3, 2] - m[0, 0] * m[2, 2] * m[3, 1] - m[2, 0] * m[0, 1] * m[3, 2] + m[2, 0] * m[0, 2] * m[3, 1] + m[3, 0] * m[0, 1] * m[2, 2] - m[3, 0] * m[0, 2] * m[2, 1]

        mat_inv[2] = m[0, 1] * m[1, 2] * m[3, 3] - m[0, 1] * m[1, 3] * m[3, 2] - m[1, 1] * m[0, 2] * m[3, 3] + m[1, 1] * m[0, 3] * m[3, 2] + m[3, 1] * m[0, 2] * m[1, 3] - m[3, 1] * m[0, 3] * m[1, 2]

        mat_inv[1, 2] = -m[0, 0] * m[1, 2] * m[3, 3] + m[0, 0] * m[1, 3] * m[3, 2] + m[1, 0] * m[0, 2] * m[3, 3] - m[1, 0] * m[0, 3] * m[3, 2] - m[3, 0] * m[0, 2] * m[1, 3] + m[3, 0] * m[0, 3] * m[1, 2]

        mat_inv[2, 2] = m[0, 0] * m[1, 1] * m[3, 3] - m[0, 0] * m[1, 3] * m[3, 1] - m[1, 0] * m[0, 1] * m[3, 3] + m[1, 0] * m[0, 3] * m[3, 1] + m[3, 0] * m[0, 1] * m[1, 3] - m[3, 0] * m[0, 3] * m[1, 1]

        mat_inv[3, 2] = -m[0, 0] * m[1, 1] * m[3, 2] + m[0, 0] * m[1, 2] * m[3, 1] + m[1, 0] * m[0, 1] * m[3, 2] - m[1, 0] * m[0, 2] * m[3, 1] - m[3, 0] * m[0, 1] * m[1, 2] + m[3, 0] * m[0, 2] * m[1, 1]

        mat_inv[0, 3] = -m[0, 1] * m[1, 2] * m[2, 3] + m[0, 1] * m[1, 3] * m[2, 2] + m[1, 1] * m[0, 2] * m[2, 3] - m[1, 1] * m[0, 3] * m[2, 2] - m[2, 1] * m[0, 2] * m[1, 3] + m[2, 1] * m[0, 3] * m[1, 2]

        mat_inv[1, 3] = m[0, 0] * m[1, 2] * m[2, 3] - m[0, 0] * m[1, 3] * m[2, 2] - m[1, 0] * m[0, 2] * m[2, 3] + m[1, 0] * m[0, 3] * m[2, 2] + m[2, 0] * m[0, 2] * m[1, 3] - m[2, 0] * m[0, 3] * m[1, 2]

        mat_inv[2, 3] = -m[0, 0] * m[1, 1] * m[2, 3] + m[0, 0] * m[1, 3] * m[2, 1] + m[1, 0] * m[0, 1] * m[2, 3] - m[1, 0] * m[0, 3] * m[2, 1] - m[2, 0] * m[0, 1] * m[1, 3] + m[2, 0] * m[0, 3] * m[1, 1]

        mat_inv[3, 3] = m[0, 0] * m[1, 1] * m[2, 2] - m[0, 0] * m[1, 2] * m[2, 1] - m[1, 0] * m[0, 1] * m[2, 2] + m[1, 0] * m[0, 2] * m[2, 1] + m[2, 0] * m[0, 1] * m[1, 2] - m[2, 0] * m[0, 2] * m[1, 1]
        
        det = m[0, 0] * mat_inv[0, 0] + m[0, 1] * mat_inv[1, 0] + m[0, 2] * mat_inv[2, 0] + m[0, 3] * mat_inv[3, 0]

        return mat_inv / det


    def inverse(self, mat):
        # return sympy.Inverse(mat)

        if(mat.shape[0] == 3):
            return self.inv3(mat)
        elif(mat.shape[0] == 4):
            return self.inv4(mat)
        else:
            return mat.inv()

    def c_gen(self, expr):
        console.print("--------------------------")
        console.print(f'Running cse')
        start = perf_counter()

        sub_expr, simpl_expr = sympy.cse(expr)
        
        printer = sympy.printing.c.C99CodePrinter()
        lines = []

        for var,expr in sub_expr:
            lines.append(f'T {var} = {printer.doprint(expr)};')

        for v in simpl_expr:
                lines.append(printer.doprint(v))

        code_string='\n'.join(lines)

        stop = perf_counter()
        console.print(f'Elapsed  {stop - start} seconds')
        console.print("--------------------------")
        console.print(f'generated code')

        console.print(code_string)
        return code_string

se = SymPyEngine()

class FE:
    def __init__(self, name, dim, symplify_expr = False):
        self.name = name
        self.dim = dim
        self.symplify_expr = symplify_expr

    def transform(self, x_ref):
        start = perf_counter()

        f = self.fun(x_ref)
        n = self.nodes

        n_fun = f.shape[0]
        x = se.zeros(self.dim, 1)

        for i in range(0, n_fun):
            x += f[i] * n[i]

        stop = perf_counter()
        console.print(f'\t- FE.transform: {stop - start} seconds')
        return x

    def jacobian(self, x_ref):
        start = perf_counter()

        x = self.transform(x_ref)
        rows = self.dim
        cols = x_ref.shape[0]

        G = se.zeros(rows, cols) 

        for i in range(0, rows):
            for j in range(0, cols):
                G[i, j] = se.diff(x[i], x_ref[j])

        if self.symplify_expr:
            ret = se.simplify(G)
        else:
            ret = G

        stop = perf_counter()
        console.print(f'\t- FE.jacobian: {stop - start} seconds')
        return ret

    def grad(self, x_ref):
        start = perf_counter()

        f = self.fun(x_ref)
        G_inv = self.jacobian_inverse(x_ref)

        n_fun = f.shape[0]
        grads = []

        for i in range(0, n_fun):
            
            gs = []

            for d in range(0, self.dim):
                gd = se.diff(f[i], x_ref[d]) 
                gs.append(gd)

            g_ref = se.array(gs)   

            # g = G_inv.T * g_ref
            g = G_inv * g_ref

            if self.symplify_expr:
                g = se.simplify(g)

            grads.append(g)
        
        stop = perf_counter()
        console.print(f'\t- FE.grad: {stop - start} seconds')
        return grads

    def jacobian_inverse(self, x_ref):
        start = perf_counter()

        G = self.jacobian(x_ref)
        G_inv = se.inverse(G)
       
        if self.symplify_expr:
            ret = se.simplify(G_inv)
        else:
            ret = G_inv

        stop = perf_counter()
        console.print(f'\t- FE.jacobian_inverse {stop - start} seconds')
        return ret

    def generate_code(self, x):
        console.print("--------------------------")
        console.print(f'{self.name}: creating expressions')
        start = perf_counter()

        f = self.fun(x)
        n_fun = f.shape[0]
    
        grads = self.grad(x)
        expr = []

        stop = perf_counter()
        console.print(f'Elapsed {stop - start} seconds')
        console.print("--------------------------")
        console.print(f'{self.name}: assignements')

        start = perf_counter()

        for i in range(0, n_fun):
            for d in range(0, self.dim):
                expr.append(Assignment(sympy.symbols(f"g{d}[{i}]"), grads[i][d]))

        stop = perf_counter()
        console.print(f'Elapsed {stop - start} seconds')
        
        return se.c_gen(expr)

class Simplex(FE):
    def __init__(self, name, dim, symplify_expr = False):
        super().__init__(name, dim, symplify_expr)

    def jacobian(self, x_ref):        
        start = perf_counter()

        n = self.nodes
        n_fun = len(n)
        G = se.zeros(n_fun-1, self.dim) 

        for i in range(1, n_fun):
            for j in range(0, self.dim):
                G[i-1, j] = n[i][j] - n[0][j]

        if self.symplify_expr:
            ret = se.simplify(G)
        else:
            ret = G

        stop = perf_counter()
        console.print(f'\t- Simplex.jacobian: {stop - start} seconds')
        return ret

class Tri3(Simplex):
    def __init__(self, symplify_expr = False):
        super().__init__('Tri3', 2, symplify_expr)
        self.nodes = se.point_symbols(3, 2)

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([ 1 - x[0] - x[1], x[0], x[1] ])
        return ret

class Tet4(Simplex):
    def __init__(self, symplify_expr = False):
        super().__init__('Tet4', 3, symplify_expr)
        self.nodes = se.point_symbols(4, 3)

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([ 1 - x[0] - x[1] - x[2], x[0], x[1], x[2] ])
        return ret

class Pentatope5(Simplex):
    def __init__(self, symplify_expr = False):
        super().__init__('Pentatope5', 4, symplify_expr)
        self.nodes = se.point_symbols(5, 4)

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([ 1 - x[0] - x[1] - x[2] - x[3], x[0], x[1], x[2], x[3]])
        return ret

class Quad4(FE):
    def __init__(self, symplify_expr = False):
        super().__init__('Quad4', 2, symplify_expr)
        self.nodes = se.point_symbols(4, 2)

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([ (1 - x[0]) * (1 - x[1]), x[0] * (1 - x[1]), x[0] * x[1], (1 - x[0] * x[1]) ])
        return ret

class Hex8(FE):
    def __init__(self, symplify_expr = False):
        super().__init__('Hex8', 3, symplify_expr)
        self.nodes = se.point_symbols(8, 3)

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([
            (1.0 - x[0]) * (1.0 - x[1]) * (1.0 - x[2]), 
            x[0] * (1.0 - x[1]) * (1.0 - x[2]), 
            x[0] * x[1] * (1.0 - x[2]), 
            (1.0 - x[0]) * x[1] * (1.0 - x[2]), 
            (1.0 - x[0]) * (1.0 - x[1]) * x[2], 
            x[0] * (1.0 - x[1]) * x[2], 
            x[0] * x[1] * x[2], 
            (1.0 - x[0]) * x[1] * x[2]
            ])
        return ret

def main(args):
    x, y, z, t = se.symbols('x y z t')
    p2 = se.point2(x, y);
    p3 = se.point3(x, y, z);
    p4 = se.point4(x, y, z, t);

    use_simplify = False

    tri3 = Tri3(use_simplify)
    quad4 = Quad4(use_simplify)
    tet4 = Tet4(use_simplify)
    hex8 = Hex8(use_simplify)
    pentatope5 = Pentatope5(use_simplify)

    # tri3.generate_code(p2)
    # quad4.generate_code(p2)
    # tet4.generate_code(p3)
    hex8.generate_code(p3)
    pentatope5.generate_code(p4)

if __name__ == '__main__':
    main(sys.argv[1:])

