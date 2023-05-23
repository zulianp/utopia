import numpy as np
import sympy as sp
from sympy.utilities.codegen import codegen
import sympy.codegen.ast as ast
import rich
from time import perf_counter

from rich.syntax import Syntax
console = rich.get_console()

def c_gen(expr, dump=False):
    console.print("--------------------------")
    console.print(f'Running cse')
    start = perf_counter()

    sub_expr, simpl_expr = sp.cse(expr)

    sub_ops = sp.count_ops(sub_expr, visual=True)
    result_ops = sp.count_ops(simpl_expr, visual=True)
    # total_ops = sub_ops
    cost = f'FLOATING POINT OPS!\n//\t- Result: {result_ops}\n//\t- Subexpressions: {sub_ops}'
    
    printer = sp.printing.c.C99CodePrinter()
    lines = []

    for var,expr in sub_expr:
        lines.append(f'const Scalar {var} = {printer.doprint(expr)};')

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

def c_code(expr):
	code_string = c_gen(expr)
	console.print(code_string)


d, d_hat = sp.symbols('d d_hat')

param = (d - d_hat) / d_hat;

# Quadratic barrier
# f = 0.5 * (param * param)

# High-order barrier
# f = sp.Rational(1,4) * (param * param * param * param)


# Arbitary order 2*N
N=3

p1, p2, p3, p4 = sp.symbols('p1 p2 p3 p4')
p = [p1, p2, p3, p4]

f = 0

mul = (param * param)
scale = 0.5
f = scale * mul
for i in range(1, N):
	pi = p[i]
	scale = scale * scale
	mul = mul * mul
	f += pi * scale * mul


grad_f = sp.simplify(sp.diff(f, d))
H_f = sp.simplify(sp.diff(grad_f, d))


console.print(f)
console.print(grad_f)
console.print(H_f)

console.print("-------------------------")
printer = sp.printing.c.C99CodePrinter()


c_code


console.print('f')
c_code(f)
console.print('grad_f')
c_code(grad_f)
console.print('H_f')
c_code(H_f)
