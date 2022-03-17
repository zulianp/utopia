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

import os
import rich
import sys
import pdb

from rich.syntax import Syntax
console = rich.get_console()

#  Work in progress output
output_dir = './workspace'
material_output_dir = './workspace'

# Where we save the final result
# output_dir = '../../../backend/kokkos/fe/generated'
material_output_dir = '../../../backend/kokkos/assembly/generated'

class SymPyEngine:
    def __init__(self):
        self.prefix = [ 'x', 'y', 'z', 't', 'a', 'b', 'c', 'd']

    def symbols(self, args, finite=True, real=True):
        return sympy.symbols(args, finite=finite, real=real)

    def pointd(self, d, prefix='', postfix=''):
        vars = []

        for k in range(0, d):
            vars.append(f'{prefix}{self.prefix[k]}{postfix}')

        return sympy.Matrix(d, 1, vars)

    def simplify(self, expr):
        return sympy.simplify(expr)

    def diff(self, f, x):
        return sympy.diff(f, x)

    def rational(self, num, denum):
        return sympy.Rational(num, denum)

    def point_symbols(self, num_points, dim):
        ret = []
        for i in range(0, num_points):
            p = []
            for j in range(0, dim):
                p.append(sympy.symbols(f'p{self.prefix[j]}[{i}]', finite=True, real=True))
                
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
    
    def det(self, mat):
        if(mat.shape[0] == 3):
            return self.det3(mat)
        elif(mat.shape[0] == 4):
            return self.det4(mat)
        else:
            # Slow
            return sympy.det(mat)

    def inv3(self, mat):
        mat_inv = self.zeros(3, 3)

        det = self.det3(mat)
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
        
    def det4(self, m):
        mat_inv00 = m[1, 1] * m[2, 2] * m[3, 3] - m[1, 1] * m[2, 3] * m[3, 2] - m[2, 1] * m[1, 2] * m[3, 3] + m[2, 1] * m[1, 3] * m[3, 2] + m[3, 1] * m[1, 2] * m[2, 3] - m[3, 1] * m[1, 3] * m[2, 2]
        mat_inv10 = -m[1, 0] * m[2, 2] * m[3, 3] + m[1, 0] * m[2, 3] * m[3, 2] + m[2, 0] * m[1, 2] * m[3, 3] - m[2, 0] * m[1, 3] * m[3, 2] - m[3, 0] * m[1, 2] * m[2, 3] + m[3, 0] * m[1, 3] * m[2, 2]
        mat_inv20 = m[1, 0] * m[2, 1] * m[3, 3] - m[1, 0] * m[2, 3] * m[3, 1] - m[2, 0] * m[1, 1] * m[3, 3] + m[2, 0] * m[1, 3] * m[3, 1] + m[3, 0] * m[1, 1] * m[2, 3] - m[3, 0] * m[1, 3] * m[2, 1]
        mat_inv30 = -m[1, 0] * m[2, 1] * m[3, 2] + m[1, 0] * m[2, 2] * m[3, 1] + m[2, 0] * m[1, 1] * m[3, 2] - m[2, 0] * m[1, 2] * m[3, 1] - m[3, 0] * m[1, 1] * m[2, 2] + m[3, 0] * m[1, 2] * m[2, 1]

        det = m[0, 0] * mat_inv00 + m[0, 1] * mat_inv10 + m[0, 2] * mat_inv20 + m[0, 3] * mat_inv30
        return det

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
        mat_inv[3, 2] = m[0, 1] * m[1, 2] * m[3, 3] - m[0, 1] * m[1, 3] * m[3, 2] - m[1, 1] * m[0, 2] * m[3, 3] + m[1, 1] * m[0, 3] * m[3, 2] + m[3, 1] * m[0, 2] * m[1, 3] - m[3, 1] * m[0, 3] * m[1, 2]
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

    def c_gen(self, expr, dump=False):
        console.print("--------------------------")
        console.print(f'Running cse')
        start = perf_counter()

        sub_expr, simpl_expr = sympy.cse(expr)

        cost = f'FLOATING POINT OPS!\n//\t- Result: {sympy.count_ops(simpl_expr, visual=True)}\n//\t- Subexpressions: {sympy.count_ops(sub_expr, visual=True)}'
        
        printer = sympy.printing.c.C99CodePrinter()
        lines = []

        for var,expr in sub_expr:
            lines.append(f'T {var} = {printer.doprint(expr)};')

        for v in simpl_expr:
                lines.append(printer.doprint(v))

        code_string=f'\n'.join(lines)

        stop = perf_counter()
        console.print(f'Elapsed  {stop - start} seconds')
        console.print("--------------------------")
        console.print(f'generated code')

        code_string = f'//{cost}\n' + code_string

        if dump:
            console.print(code_string)

        return code_string

se = SymPyEngine()

class Template:
    def __init__(self, header_tpl_path, impl_tpl_path, header_output_path, impl_output_path):
    
        self.header_tpl_path = header_tpl_path
        self.impl_tpl_path = impl_tpl_path
        self.header_output_path = header_output_path
        self.impl_output_path = impl_output_path

        with open(header_tpl_path, 'r') as f:
            self.header_tpl = f.read()

        with open(impl_tpl_path, 'r') as f:
            self.impl_tpl = f.read()

class FE:
    def __init__(self, name, dim, symplify_expr = False, symbolic_map_inversion = True):
        self.name = name
        self.dim = dim
        self.symplify_expr = symplify_expr
        self.symbolic_map_inversion = symbolic_map_inversion
        self.order = 1  # FIXME
        
        self.tpl = Template(
            "templates/fe/utopia_tpl_fe.hpp",
            f"templates/fe/utopia_tpl_fe_{dim}_impl.hpp",
            f"{output_dir}/../utopia_fe_{name}.hpp",
            f"{output_dir}/utopia_fe_{name}_{dim}.hpp")

        with open("templates/fe/utopia_tpl_fe_map_inversion_snippet.hpp", 'r') as f:
            self.map_inversion_tpl = f.read()

    def coeff(self, name):
        coeffs = []

        n_fun = self.n_shape_functions()
        
        for i in range(0, n_fun):
            coeffs.append(se.symbols(f'{name}[{i}]'))

        return se.array(coeffs)
    
    def interpolate(self, coeff, x_ref):
        f = self.fun(x_ref)

        ret = 0
        n_fun = self.n_shape_functions()

        for i in range(0, n_fun):
            ret += coeff[i] * f[i]

        return ret

    def grad_interpolate(self, coeff, x_ref):
        g = self.grad(x_ref)

        ret = se.zeros(self.dim, 1)

        n_fun = self.n_shape_functions()

        for i in range(0, n_fun):
            for d in range(0, self.dim):
                ret[d] += coeff[i] * g[i][d]

        return ret

    def reference_measure(self):
        return 1

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

    def measure(self, x_ref):
        return se.det(self.jacobian(x_ref))

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

    def generate_code(self, x_ref):
        console.print("--------------------------")
        console.print(f'{self.name}: creating expressions', style="magenta")
        start = perf_counter()

        f = self.fun(x_ref)
        n_fun = f.shape[0]
    
        grads = self.grad(x_ref)
        
        stop = perf_counter()
        console.print(f'Elapsed {stop - start} seconds')
        console.print("--------------------------")
        console.print(f'{self.name}: assignements')

        start = perf_counter()

        grad_expr = []
        for i in range(0, n_fun):
            for d in range(0, self.dim):
                grad_expr.append(Assignment(se.symbols(f"g{se.prefix[d]}[{i}]"), grads[i][d]))

        value_expr = []
        for i in range(0, n_fun):
            value_expr.append(Assignment(se.symbols(f"f[{i}]"), f[i]))

        stop = perf_counter()
        console.print(f'Elapsed {stop - start} seconds')

        measure_expr = [Assignment(se.symbols(f"measure_value"), self.measure(x_ref))]

        combined_expression = []
        combined_expression.extend(value_expr)
        combined_expression.extend(measure_expr)
        combined_expression.extend(grad_expr)
            
        #############################################################
        # Standard assembly values
        #############################################################

        grad_code = se.c_gen(grad_expr)
        value_code = se.c_gen(value_expr)
        measure_code = se.c_gen(measure_expr)
        combined_code = se.c_gen(combined_expression)

        #############################################################
        # Jacobian
        #############################################################

        J = self.jacobian(x_ref)
        J_inv = self.jacobian_inverse(x_ref)

        jacobian_expr = []
        jacobian_inverse_expr = []

        for d1 in range(0, self.dim):
            for d2 in range(0, self.dim):
                jacobian_expr.append(Assignment(se.symbols(f"J[{d1*self.dim + d2}]"), J[d1, d2]))
                jacobian_inverse_expr.append(Assignment(se.symbols(f"J_inv[{d1*self.dim + d2}]"), J_inv[d1, d2]))

        jacobian_code = se.c_gen(jacobian_expr)
        jacobian_inverse_code = se.c_gen(jacobian_inverse_expr)

        #############################################################
        # Geometric transformation
        #############################################################

        transform_expr = []
        trafo = self.transform(x_ref)

        for d in range(0, self.dim):
            transform_expr.append(Assignment(se.symbols(f"t{se.prefix[d]}"), trafo[d]))

        inverse_transform_expr = []
        transform_code = se.c_gen(transform_expr)

        # Point in physical space
        x = se.pointd(self.dim, 't')

        if self.symbolic_map_inversion:
            
            inv_trafo = sympy.solve(self.transform(x_ref) - x, x_ref, positive=True, simplify=self.symplify_expr, dict=True)

            for d in range(0, self.dim):
                inverse_transform_expr.append(Assignment(se.symbols(f"{se.prefix[d]}"), inv_trafo[0][x_ref[d]]))
                inverse_transform_code = se.c_gen(inverse_transform_expr)
        else:
            g = x - self.transform(x_ref)
            c = J_inv * g
            x = x + c 

            norm_g = 0

            for d in range(0, self.dim):
                norm_g += g[d] * g[d]

            for d in range(0, self.dim):
                inverse_transform_expr.append(Assignment(se.symbols(f"{se.prefix[d]}"), x[d]))

            inverse_transform_expr.append(Assignment(se.symbols("residual_norm"), norm_g))
            update = se.c_gen(inverse_transform_expr)


            inverse_transform_code = self.map_inversion_tpl.format(
                update=update
                )

        kernel = self.tpl.impl_tpl.format(
            name=self.name,
            measure=measure_code,
            value=value_code,
            gradient=grad_code,
            combined=combined_code,
            transform=transform_code,
            inverse_transform=inverse_transform_code,
            jacobian=jacobian_code,
            jacobian_inverse=jacobian_inverse_code,
            hessian='',
            dim=self.dim,
            nnodes=self.n_shape_functions(),
            order=self.order)

        with open(self.tpl.impl_output_path, 'w') as f:
            f.write(kernel)

class Simplex(FE):
    def __init__(self, name, dim, symplify_expr = False, symbolic_map_inversion = True):
        super().__init__(name, dim, symplify_expr, symbolic_map_inversion)

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
    def __init__(self, symplify_expr = True):
        super().__init__('Tri3', 2, symplify_expr)
        self.nodes = se.point_symbols(3, 2)

    def reference_measure(self):
        return 0.5

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([ 1 - x[0] - x[1], x[0], x[1] ])
        return ret

class Tet4(Simplex):
    def reference_measure(self):
        return se.rational(1, 6)

    def __init__(self, symplify_expr = False):
        super().__init__('Tet4', 3, symplify_expr)
        self.nodes = se.point_symbols(4, 3)

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([ 1 - x[0] - x[1] - x[2], x[0], x[1], x[2] ])
        return ret

class Pentatope5(Simplex):
    def __init__(self, symplify_expr = False, symbolic_map_inversion = True):
        super().__init__('Pentatope5', 4, symplify_expr, symbolic_map_inversion)
        self.nodes = se.point_symbols(5, 4)

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([ 1 - x[0] - x[1] - x[2] - x[3], x[0], x[1], x[2], x[3]])
        return ret

class Quad4(FE):
    def __init__(self, symplify_expr = True):
        super().__init__('Quad4', 2, symplify_expr)
        self.nodes = se.point_symbols(4, 2)

    def n_shape_functions(self):
        return len(self.nodes)

    def fun(self, x):
        ret = se.array([ (1 - x[0]) * (1 - x[1]), x[0] * (1 - x[1]), x[0] * x[1], (1 - x[0] * x[1]) ])
        return ret

class AxisAlignedQuad4(FE):
    def __init__(self, symplify_expr = True):
        super().__init__('AxisAlignedQuad4', 2, symplify_expr)
        self.nodes = se.point_symbols(4, 2)
        self.quad4 = Quad4(symplify_expr)

    def n_shape_functions(self):
        return len(self.nodes)

    def transform(self, x_ref):
        # pdb.set_trace()
        start = perf_counter()

        scaling = self.nodes[2] - self.nodes[0]

        ret = se.point2(0, 0) + self.nodes[0]
        for d in range(0, self.dim):
            ret[d] += x_ref[d] * scaling[d]

        ret = se.simplify(ret)

        stop = perf_counter()
        return ret

    def fun(self, x):
        return self.quad4.fun(x)

class Hex8(FE):
    def __init__(self, symplify_expr = False, symbolic_map_inversion = False):
        super().__init__('Hex8', 3, symplify_expr, symbolic_map_inversion)
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

class AxisAlignedHex8(FE):
    def __init__(self, symplify_expr = True, symbolic_map_inversion = True):
        super().__init__('AxisAlignedHex8', 3, symplify_expr, symbolic_map_inversion)
        self.nodes = se.point_symbols(8, 3)
        self.hex8 = Hex8(symplify_expr, symbolic_map_inversion)

    def n_shape_functions(self):
        return len(self.nodes)

    def transform(self, x_ref):
        scaling = self.nodes[6] - self.nodes[0]

        ret = se.point3(0, 0, 0) + self.nodes[0]
        for d in range(0, self.dim):
            ret[d] += x_ref[d] * scaling[d]

        return se.simplify(ret)

    def fun(self, x):
        return self.hex8.fun(x)

class Material:
    def __init__(self, name):
        self.name = name
        self.tpl = ''

    def is_linear(self):
        return False

    def load_tpl(self):
        with open(f"templates/material/utopia_tpl_material_{self.test.dim}_impl.hpp", 'r') as f:
            self.tpl = f.read()

        with open(f"templates/material/utopia_tpl_material.hpp", 'r') as f:
            self.header_tpl = f.read()

    def generate_code(self, x_ref, w):
        H = self.hessian(x_ref, w)

        H_expr = []
        for i in range(0, H.shape[0]):
            for j in range(0, H.shape[1]):
                H_expr.append(AddAugmentedAssignment(se.symbols(f"H[{i*H.shape[1] + j}]"), H[i, j]))

        H_code = se.c_gen(H_expr)

        trial = self.trial
        u = trial.coeff("u")
        Hx = self.apply(u, x_ref, w)

        Hx_expr = []

        for i in range(0, Hx.shape[0]):
                Hx_expr.append(AddAugmentedAssignment(se.symbols(f"Hx[{i}]"), Hx[i]))

        Hx_code = se.c_gen(Hx_expr)


        # FIXME
        g_expr = []
        for i in range(0, Hx.shape[0]):
                g_expr.append(AddAugmentedAssignment(se.symbols(f"g[{i}]"), Hx[i]))


        g_code = se.c_gen(g_expr)

        combined_expr = []
        combined_expr.extend(H_expr)
        combined_expr.extend(g_expr)

        combined_code = se.c_gen(combined_expr)

        self.load_tpl()

        header_code = self.header_tpl.format(name=self.name)

        code = self.tpl.format(
            name=self.name,
            dim=self.test.dim,
            trial=self.trial.name,
            test=self.test.name,
            hessian = H_code,
            apply_hessian = Hx_code,
            gradient = g_code,
            value = '//TODO',
            combined = combined_code,
            get_params = '//TODO',
            set_params = '//TODO',
            fields = '//TODO'
        )

        with open(material_output_dir + f'/utopia_material_{self.name}.hpp', 'w') as f:
            f.write(header_code)

        with open(material_output_dir + f'/utopia_material_{self.name}_{self.test.name}_{self.test.dim}_impl.hpp', 'w') as f:
            f.write(code)

        return code

class LaplaceOperator(Material):
    def __init__(self, trial, test):
        super().__init__('LaplaceOperator')

        self.trial = trial
        self.test = test

    def is_linear(self):
        return True

    def hessian(self, x_ref, w):
        trial = self.trial
        test = self.test

        g_trial = trial.grad(x_ref)
        g_test = test.grad(x_ref)
        dx = test.measure(x_ref) * w

        rows = test.n_shape_functions()
        cols = trial.n_shape_functions()
        
        M = se.zeros(rows, cols)
        dim = test.dim

        for i in range(0, rows):
            for j in range(0, cols):
                scalar_prod = 0

                for d in range(0, dim):
                    scalar_prod += g_test[i][d] * g_trial[j][d]

                scalar_prod *=  dx
                M[i, j] = scalar_prod

        return M

    def apply(self, coeff, x_ref, w):
        trial = self.trial
        test = self.test

        g_trial = trial.grad_interpolate(coeff, x_ref)
        g_test = test.grad(x_ref)
        dx = test.measure(x_ref) * w

        rows = test.n_shape_functions()
        cols = trial.n_shape_functions()
        
        Mx = se.zeros(rows, cols)
        dim = test.dim

        for i in range(0, rows):
            for j in range(0, cols):
                scalar_prod = 0

                for d in range(0, dim):
                    scalar_prod += g_test[i][d] * g_trial[d]

                scalar_prod *=  dx
                Mx[i] += scalar_prod

        return Mx

class Mass(Material):
    def __init__(self, trial, test):
        super().__init__('Mass')

        self.trial = trial
        self.test = test

    def is_linear(self):
        return True

    def hessian(self, x_ref, w):
        trial = self.trial
        test = self.test

        f_trial = trial.fun(x_ref)
        f_test = test.fun(x_ref)
        dx = test.measure(x_ref) * w

        rows = test.n_shape_functions()
        cols = trial.n_shape_functions()
        
        M = se.zeros(rows, cols)

        for i in range(0, rows):
            for j in range(0, cols):
                M[i, j] = f_test[i] * f_trial[j] * dx

        return M

    def apply(self, coeff, x_ref, w):
        trial = self.trial
        test = self.test

        f_trial = trial.interpolate(coeff, x_ref)
        f_test = test.fun(x_ref)
        dx = test.measure(x_ref) * w

        rows = test.n_shape_functions()
        cols = trial.n_shape_functions()
        
        Mx = se.zeros(rows, cols)
        dim = test.dim

        for i in range(0, rows):
            Mx[i] += f_test[i] * f_trial * dx

        return Mx

def main(args):
    x, y, z, t = se.symbols('x y z t')
    p2 = se.point2(x, y)
    p3 = se.point3(x, y, z)
    p4 = se.point4(x, y, z, t)

    # console.print(x.assumptions0)
    use_simplify = False
    symbolic_map_inversion = False
    generate_fe = False

    tri3 = Tri3()
    quad4 = Quad4()
    aaquad4 = AxisAlignedQuad4()

    tet4 = Tet4()
    hex8 = Hex8(use_simplify, symbolic_map_inversion)
    aahex8 = AxisAlignedHex8()
    pentatope5 = Pentatope5()

    elements_2D = [tri3, quad4, aaquad4]
    elements_3D = [tet4, hex8, aahex8]
    elements_4D = [pentatope5]

    elements = []
    elements.extend(elements_2D)
    elements.extend(elements_3D)
    elements.extend(elements_4D)

    if generate_fe:
        for e in elements_2D:
            e.generate_code(p2)

        for e in elements_3D:
            e.generate_code(p3)

        for e in elements_4D:
            e.generate_code(p4)


    w = se.symbols('weight')

    for e in elements_2D:
        mass = Mass(e, e)
        mass.generate_code(p2, w)

    for e in elements_3D:
        mass = Mass(e, e)
        mass.generate_code(p3, w)

    for e in elements_4D:
        mass = Mass(e, e)
        mass.generate_code(p4, w)

if __name__ == '__main__':
    main(sys.argv[1:])

