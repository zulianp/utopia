import sympy as sp
from sympy.printing.c import C99CodePrinter
import rich 


console = rich.get_console()

d, d_hat = sp.symbols('d d_hat')

param = (d - d_hat) / d_hat;

# Quadratic barrier
# f = 0.5 * (param * param)

# High-order barrier
f = sp.Rational(1,4) * (param * param * param * param)

grad_f = sp.simplify(sp.diff(f, d))
H_f = sp.simplify(sp.diff(grad_f, d))


# console.print(f)
# console.print(grad_f)
# console.print(H_f)

printer = sp.printing.c.C99CodePrinter()


console.print('f')
console.print(printer.doprint(f))
console.print('grad_f')
console.print(printer.doprint(grad_f))
console.print('H_f')
console.print(printer.doprint(H_f))