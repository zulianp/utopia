from sympy import symbols
from sympy import Matrix
from sympy import shape
from sympy import diff
from sympy.printing.c import C99CodePrinter
from sympy.utilities.codegen import codegen
from sympy import cse


# def deformation_gradient(d):
#     if d == 1:
#         f00 = symbols('f_00')
#         F = Matrix(1,1,[f00])
#     elif d == 2:
#         f00, f01, f10, f11 = symbols('f_00 f_01 f_10 f_11')
#         F = Matrix(2,2,[f00,f01,f10,f11])
#     else:
#         f00, f01, f02, f10, f11, f12, f20, f21, f22  = symbols('f_00 f_01 f_02 f_10 f_11 f_12 f_20 f_21 f_22')
#         F = Matrix(3,3,[f00, f01, f02, f10, f11, f12, f20, f21, f22])
#     return F

def deformation_gradient(d):
    if d == 1:
        f00 = symbols('f[0]')
        F = Matrix(1,1,[f00])
    elif d == 2:
        f00, f01, f10, f11 = symbols('f[0] f[1] f[2] f[3]')
        F = Matrix(2,2,[f00,f01,f10,f11])
    else:
        f00, f01, f02, f10, f11, f12, f20, f21, f22  = symbols('f[0] f[1] f[2] f[3] f[4] f[5] f[6] f[7] f[8]')
        F = Matrix(3,3,[f00, f01, f02, f10, f11, f12, f20, f21, f22])
    return F

def trial_function(d):
    if d == 1:
        return symbols('trial_0')
    elif d == 2:
       trial_00, trial_01 = symbols('trial_00 trial_01')
       trial = Matrix(2,1,[trial_00, trial_01])
    elif d == 3:
       trial_00, trial_01, trial_02 = symbols('trial_00 trial_01 trial_02')
       trial = Matrix(3,1,[trial_00, trial_01, trial_02])
    return trial

def trial_gradient(d):
    if d == 1:
        return symbols('trial_grad_0')
    elif d == 2:
       trial_grad_00, trial_grad_01, trial_grad_10, trial_grad_11 = symbols('trial_grad_00 trial_grad_01 trial_grad_10 trial_grad_11')
       trial = Matrix(2,2,[trial_grad_00, trial_grad_01, trial_grad_10, trial_grad_11])
    elif d == 3:
       trial_grad_00, trial_grad_01, trial_grad_02, trial_grad_10, trial_grad_11, trial_grad_12, trial_grad_20, trial_grad_21, trial_grad_22 = symbols('trial_grad_00 trial_grad_01 trial_grad_02 trial_grad_10 trial_grad_11 trial_grad_12 trial_grad_20 trial_grad_21 trial_grad_22')
       trial = Matrix(3,3,[
        trial_grad_00, trial_grad_01, trial_grad_02,
        trial_grad_10, trial_grad_11, trial_grad_12,
        trial_grad_20, trial_grad_21, trial_grad_22])
    return trial

def test_function(d):
    if d == 1:
        return symbols('test_0')
    elif d == 2:
       test_00, test_01 = symbols('test_00 test_01')
       test = Matrix(2,1,[test_00, test_01])
    elif d == 3:
       test_00, test_01, test_02 = symbols('test_00 test_01 test_02')
       test = Matrix(3,1,[test_00, test_01, test_02])
    return trial

def test_gradient(d):
    if d == 1:
        return symbols('test_grad_0')
    elif d == 2:
       test_grad_00, test_grad_01, test_grad_10, test_grad_11 = symbols('test_grad_00 test_grad_01 test_grad_10 test_grad_11')
       test = Matrix(2,2,[test_grad_00, test_grad_01, test_grad_10, test_grad_11])
    elif d == 3:
       test_grad_00, test_grad_01, test_grad_02,test_grad_10, test_grad_11, test_grad_12, test_grad_20, test_grad_21, test_grad_22 = symbols('test_grad_00 test_grad_01 test_grad_02 test_grad_10 test_grad_11 test_grad_12 test_grad_20 test_grad_21 test_grad_22')
       test = Matrix(3,3,[
        test_grad_00, test_grad_01, test_grad_02,
        test_grad_10, test_grad_11, test_grad_12,
        test_grad_20, test_grad_21, test_grad_22])
    return test


class TensorProductBasis:
    def __init__(self, dim):
      
        self.dim = dim

    def bilinear_subs(self, i_trial, d_trial, i_test, d_test, expr):
        return self.subs_gradient_test(i_test, d_test, self.subs_gradient_trial(i_trial, d_trial, expr))

    def linear_subs(self, i_test, d_test, expr):
        return self.subs_gradient_test(i_test, d_test, expr);

    def subs_gradient_trial(self, i, d, f):
        dim = self.dim

        trial = f"grad_{i}";
        
        ret = f;
        for k in range(0, dim):
            ret = ret.subs(f'trial_grad_{d}{k}', symbols(f'{trial}[{k}]'))

        for k in range(0, dim):
            if k != d:
                for j in range(0, dim):
                    ret = ret.subs(f'trial_grad_{k}{j}', 0)

        return ret

    def subs_gradient_test(self, i, d,f):
        dim = self.dim

        test = f"grad_{i}";

        ret = f;
        for k in range(0, dim):
            ret = ret.subs(f'test_grad_{d}{k}', symbols(f'{test}[{k}]'))

        for k in range(0, dim):
            if k != d:
                for j in range(0, dim):
                    ret = ret.subs(f'test_grad_{k}{j}', 0)

        return ret


def first_piola(strain_energy, F):
    shape = F.shape
    # print(shape)

    P = Matrix.zeros(shape[0], shape[1])

    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            P[i, j] = diff(strain_energy, F[i, j])

    return P


class KernelGenerator:
    def __init__(self, dim):
        self.dim = dim

    def generate(self, expression_list, tpl_path, output_path):
        print("Generating code")

        sub_expr, simpl_expr = cse(expression_list)

        lines = []
        printer = C99CodePrinter()

        for var,expr in sub_expr:
            lines.append(f'T {var} = {printer.doprint(expr)};')

        for v in simpl_expr:
                lines.append(printer.doprint(v))

        code_string='\n'.join(lines)

        with open(tpl_path, 'r') as f:
            tpl = f.read()
            kernel = tpl.format(code=code_string,dim=self.dim)

            with open(output_path, 'w') as f:
                f.write(kernel)

