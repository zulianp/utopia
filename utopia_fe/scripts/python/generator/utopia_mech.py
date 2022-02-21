from sympy import *

# from sympy import Matrix
# from sympy import cse
# from sympy import diff
# from sympy import shape
# from sympy import symbols

from sympy.codegen.ast import AddAugmentedAssignment
from sympy.codegen.ast import Assignment
from sympy.printing.c import C99CodePrinter
from sympy.utilities.codegen import CCodeGen
from sympy.utilities.codegen import codegen


import os
import rich
import sys

from rich.syntax import Syntax

console = rich.get_console()



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

    def bilinear_subs(self, i_test, d_test, i_trial, d_trial, expr):
        return self.subs_test(i_test, d_test, self.subs_trial(i_trial, d_trial, expr))

    def linear_subs(self, i_test, d_test, expr):
        return self.subs_test(i_test, d_test, expr);

    def subs_test(self, i_test, d_test, expr):
        ret = expr

        for i in range(0, self.dim):
            if d_test == i:
                ret = ret.subs(f'test_{i}', f'fun_{i_test}')
            else:
                ret = ret.subs(f'test_{i}', 0)

        return ret
        

    def subs_trial(self, i_trial, d_trial, expr):
        ret = expr

        for i in range(0, self.dim):
            if d_trial == i:
                ret = ret.subs(f'trial_{i}', f'fun_{i_trial}')
            else:
                ret = ret.subs(f'trial_{i}', 0)
        return ret


    def bilinear_subs_gradients(self, i_trial, d_trial, i_test, d_test, expr):
        return self.subs_gradient_test(i_test, d_test, self.subs_gradient_trial(i_trial, d_trial, expr))

    def linear_subs_gradients(self, i_test, d_test, expr):
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


class KernelGenerator:
    def __init__(self, dim):
        self.dim = dim

    def generate_fields(self, params):
        class_fields = []
        for p in params:
            class_fields.append(f'T {p[0]}{{{p[1]}}};')

        return "\n".join(class_fields)

    def generate_read_parameters(self, params):
        class_fields = []
        for p in params:
            class_fields.append(f'in.get(\"{p[0]}\", {p[0]});')

        return "\n".join(class_fields)

    def generate_set_parameters(self, params):
        class_fields = []
        for p in params:
            class_fields.append(f'{p[0]} = params.{p[0]};')

        return "\n".join(class_fields)


    def generate_class(
        self, tpl, model,
        value_expression,
        gradient_expression,
        hessian_expression,
        output_dir):

        combined_expression = []
        combined_expression.extend(value_expression)
        combined_expression.extend(gradient_expression)
        combined_expression.extend(hessian_expression)

        value_string = self.generate_string(value_expression)
        gradient_string = self.generate_string(gradient_expression)
        hessian_string = self.generate_string(hessian_expression)
        combined_string = self.generate_string(combined_expression)

        header = tpl.header_tpl.format(
            name=model.name,
            dim=self.dim)

        fields = self.generate_fields(model.params)
        set_params = self.generate_set_parameters(model.params)

        if(model.use_default_parameter_reader):
            get_params = """StressStrainParameters<T, T> ssp;
                            ssp.read(in);

                            lambda = ssp.first_lame_parameter.get();
                            mu = ssp.shear_modulus.get();"""
        else:
            get_params = self.generate_read_parameters(model.params)

        kernel = tpl.impl_tpl.format(
            name=model.name,
            value=value_string,
            gradient=gradient_string,
            hessian=hessian_string,
            combined=combined_string,
            dim=self.dim,
            fields=fields,
            get_params=get_params,
            set_params=set_params)

        with open(tpl.impl_output_path, 'w') as f:
            f.write(kernel)

        with open(tpl.header_output_path, 'w') as f:
            f.write(header)

    def generate_string(self, expression):
        sub_expr, simpl_expr = cse(expression)

        lines = []
        printer = C99CodePrinter()

        for var,expr in sub_expr:
            lines.append(f'T {var} = {printer.doprint(expr)};')

        for v in simpl_expr:
                lines.append(printer.doprint(v))

        code_string='\n'.join(lines)
        return code_string



    def generate_files(self, expression_list, tpl_path, output_path):
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





class HyperElasticModel:
    def __init__(self, d):
        self.d = d
        self.block_size = d
        self.F = deformation_gradient(d)
        self.J = det(self.F)
        self.F_inv = Inverse(self.F)
        self.F_inv_t = self.F_inv.T
        self.C = self.F.T*self.F
        self.I_C = trace(self.C)
        self.fun = 0

        self.form = zeros(3)
        self.use_default_parameter_reader = False
        self.is_block_system = False

    def compute_forms(self):
        params = self.params
        W = self.fun
        name = self.name
        F = self.F
        d = self.d

        #############################################
        # FE
        #############################################

        d = F.shape[0]
        grad_trial = trial_gradient(d)
        grad_test = test_gradient(d)

        dx = symbols('dx')
        self.dx = dx

        #############################################
        # Energy
        #############################################

        energy = W * dx

        self.form[0] = energy
        #############################################
        # First-Piola-Kirchoff stress tensor
        #############################################

        P = first_piola(W, F)
        self.P = P

        #############################################
        # Gradient
        #############################################
        console.print("Creating Gradient")

        linear_form = 0
        for i in range(0, d):
            for j in range(0, d):
                linear_form += P[i,j] * grad_test[i,j]

        linear_form *= dx

        self.form[1] = linear_form

        #############################################
        # Hessian
        #############################################
        console.print("Creating Hessian")

        contraction = 0
        for i in range(0, d):
            for j in range(0, d):
                contraction += P[i,j] * grad_trial[i,j]

        bilinear_form = 0

        for i in range(0, d):
            for j in range(0, d):
                Hij = diff(contraction, F[i, j]);
                bilinear_form += Hij * grad_test[i,j]

        bilinear_form *= dx

        self.form[2] = bilinear_form

    def generate_files(self, output_dir, simplify_expressions):
        model = self
        #############################################
        # UI
        #############################################

        if not os.path.exists(output_dir):
            console.print(f"Creating directory {output_dir}")
            os.mkdir(output_dir)

        #############################################
        # FE
        #############################################

        model.compute_forms()
        tp = TensorProductBasis(model.d)

        hessian_expression_list = []
        gradient_expression_list = []
        energy_expression_list = []

        for d1 in range(0, model.d):
            subsituted = tp.linear_subs_gradients("test", d1, model.form[1])

            if simplify_expressions:
                subsituted = simplify(subsituted)

            gradient_expression_list.append(AddAugmentedAssignment(symbols(f"lf[{d1}]"), subsituted))

            for d2 in range(0, model.d):
                subsituted = tp.bilinear_subs_gradients("test", d1, "trial", d2, model.form[2])

                if simplify_expressions:
                    subsituted = simplify(subsituted)

                hessian_expression_list.append(AddAugmentedAssignment(symbols(f"bf[{d1*model.block_size + d2}]"), subsituted))


        energy_expression_list.append(AddAugmentedAssignment(symbols("e"), model.form[0]))

        full_expression_list = []
        full_expression_list.extend(hessian_expression_list)
        full_expression_list.extend(gradient_expression_list)
        full_expression_list.extend(energy_expression_list)

        #############################################
        # Generate code
        #############################################

        console.print("Generating code")

        generator = KernelGenerator(model.d)
        generator.generate_class(
            Template(
                "templates/utopia_tpl_hyperelasticity.hpp",
                "templates/utopia_tpl_hyperelasticity_impl.hpp",
                f"{output_dir}/../utopia_hyperelasticity_{model.name}.hpp",
                f"{output_dir}/utopia_hyperelasticity_{model.name}_{model.d}.hpp"),
            model,
            energy_expression_list,
            gradient_expression_list,
            hessian_expression_list,
            output_dir)

class IncompressibleHyperElasticModel(HyperElasticModel):
    def __init__(self, d):
        super().__init__(d)
        self.block_size = d + 1 # Add one element for the pressure
        self.p = symbols('p')
        # self.J = symbols('J')

        self.form = [0,0,0]
        self.form[1] = zeros(2) # Two blocks
        self.form[2] = zeros(2, 2) # Four blocks

        self.is_block_system = True
        self.independent_variables = [self.F, self.p]

    def compute_forms(self):
        params = self.params
        W = self.fun
        name = self.name
        F = self.F
        J = self.J
        d = self.d
        p = self.p
        
        dx = symbols('dx')

        trial = trial_function(1)
        test = test_function(1)

        self.trial = trial
        self.test = test

        grad_trial = trial_gradient(d)
        grad_test = test_gradient(d)

        P = first_piola(self.fun, F)

        #############################################
        # Energy
        #############################################
        
        console.print("Function(u,p)")
        energy = self.fun * dx
        self.form[0] = energy

        #############################################
        # Gradient
        #############################################
        console.print("Gradient(u,p)")

        linear_form_0 = 0
        for i in range(0, d):
            for j in range(0, d):
                linear_form_0 += P[i,j] * grad_test[i,j]

        linear_form_0 *= dx

        self.form[1][0] = linear_form_0

        # console.print(linear_form_0)

        dWdp = simplify(diff(self.fun, p))
        linear_form_1 = dWdp * test * dx

        self.form[1][1] = linear_form_1
        
        #############################################
        # Hessian
        #############################################
        console.print("Hessian(u,p)")

        contraction = 0
        for i in range(0, d):
            for j in range(0, d):
                contraction += P[i,j] * grad_trial[i,j]

        bilinear_form_00 = 0

        for i in range(0, d):
            for j in range(0, d):
                Hij = diff(contraction, F[i, j]);
                bilinear_form_00 += Hij * grad_test[i,j]

        bilinear_form_00 *= dx

        dWdpdF = first_piola(dWdp * trial, F)

        bilinear_form_10 = simplify(diff(contraction, p) * test) * dx
        bilinear_form_11 = simplify(diff(dWdp * trial, p)) * test * dx

        contraction = 0
        for i in range(0, d):
            for j in range(0, d):
                contraction += dWdpdF[i,j] * grad_test[i,j]

        bilinear_form_01 = simplify(contraction) * dx

        self.form[2][0,0] = bilinear_form_00
        self.form[2][0,1] = bilinear_form_01
        self.form[2][1,0] = bilinear_form_10
        self.form[2][1,1] = bilinear_form_11

        # console.print('bilinear_form_00')
        # console.print(simplify(bilinear_form_00))

        # console.print('bilinear_form_01')
        # console.print(simplify(bilinear_form_01))

        # console.print('bilinear_form_10')
        # console.print(simplify(bilinear_form_10))

        # console.print('bilinear_form_11')
        # console.print(bilinear_form_11)

        # console.print('linear_form_0')
        # console.print(linear_form_0)

        # console.print('linear_form_1')
        # console.print(linear_form_1)

    def generate_files(self, output_dir, simplify_expressions):
        model = self

        if not os.path.exists(output_dir):
            console.print(f"Creating directory {output_dir}")
            os.mkdir(output_dir)

        #############################################
        # FE
        #############################################

        model.compute_forms()
        tp = TensorProductBasis(model.d)

        hessian_expression_list = []
        gradient_expression_list = []
        energy_expression_list = []


        lf0 = model.form[1][0]
        lf1 = model.form[1][1]

        bf00 = model.form[2][0,0]
        bf01 = model.form[2][0,1]
        bf10 = model.form[2][1,0]
        bf11 = model.form[2][1,1]


        for i in range(0, 2):
            for j in range(0, 2):
                console.print(f'{i},{j})\n', style='bold blue')
                # console.print( Syntax(f'{model.form[2][i,j]}\n', 'C'))
                console.print(f'{model.form[2][i,j]}\n')

        # Displacement

        for d1 in range(0, model.d):
            subsituted_0 = tp.linear_subs_gradients("test", d1, lf0)

            if simplify_expressions:
                subsituted_0 = simplify(subsituted_0)

            gradient_expression_list.append(AddAugmentedAssignment(symbols(f"lf[{d1}]"), subsituted_0))

            for d2 in range(0, model.d):
                subsituted = tp.bilinear_subs_gradients("test", d1, "trial", d2, bf00)

                if simplify_expressions:
                    subsituted = simplify(subsituted)

                hessian_expression_list.append(AddAugmentedAssignment(symbols(f"bf[{d1*model.block_size + d2}]"), subsituted))

        # Pressure 
        subsituted_1 = tp.linear_subs('test', 0, lf1)

        if simplify_expressions:
            subsituted_1 = simplify(subsituted_1)

        gradient_expression_list.append(AddAugmentedAssignment(symbols(f"lf[{model.d}]"), subsituted_1))


        subsituted_11 = tp.bilinear_subs("test", 0, "trial", 0, bf11)

        if simplify_expressions:
            subsituted_11 = simplify(subsituted_11)

        hessian_expression_list.append(AddAugmentedAssignment(symbols(f"bf[{model.d*model.block_size +model.d}]"), subsituted_11))

        # Mixed

        # bf(p, delta_u)
        mixed_10 = tp.subs_test('test', 0, bf10)

        for d1 in range(0, model.d):
            subsituted_01 = tp.subs_gradient_trial("trial", d1, mixed_10)

            if simplify_expressions:
                subsituted_01 = simplify(subsituted_01)

            hessian_expression_list.append(AddAugmentedAssignment(symbols(f"bf[{d1+ model.d*model.block_size}]"), subsituted_01))


        # bf(u,delta_p)
        mixed_01 = tp.subs_trial('trial', 0, bf01)

        for d1 in range(0, model.d):
            subsituted_10 = tp.subs_gradient_test("test", d1, mixed_01)

            if simplify_expressions:
                subsituted_10 = simplify(subsituted_10)

            hessian_expression_list.append(AddAugmentedAssignment(symbols(f"bf[{d1*model.block_size + model.d}]"), subsituted_10))


        # Model Energy
        energy_expression_list.append(AddAugmentedAssignment(symbols("e"), model.form[0]))

        full_expression_list = []
        full_expression_list.extend(hessian_expression_list)
        full_expression_list.extend(gradient_expression_list)
        full_expression_list.extend(energy_expression_list)

        #############################################
        # Generate code
        #############################################

        console.print("Generating code")

        generator = KernelGenerator(model.d)
        generator.generate_class(
            Template(
                "templates/utopia_tpl_incompressible_hyperelasticity.hpp",
                "templates/utopia_tpl_incompressible_hyperelasticity_impl.hpp",
                f"{output_dir}/../utopia_hyperelasticity_{model.name}.hpp",
                f"{output_dir}/utopia_hyperelasticity_{model.name}_{model.d}.hpp"),
            model,
            energy_expression_list,
            gradient_expression_list,
            hessian_expression_list,
            output_dir)