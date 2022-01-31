from sympy import symbols
from sympy import Matrix
from sympy import shape
from sympy import diff

def deformation_gradient(d):
    if d == 1:
        f00 = symbols('f00')
        F = Matrix(1,1,[f00])
    elif d == 2:
        f00, f01, f10, f11 = symbols('f_00 f_01 f_10 f_11')
        F = Matrix(2,2,[f00,f01,f10,f11])
    else:
        f00, f01, f02, f10, f11, f12, f20, f21, f22  = symbols('f_00 f_01 f_02 f_10 f_11 f_12 f_20 f_21 f_22')
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
       trial_grad_00, trial_grad_01, trial_grad_02,
       trial_grad_10, trial_grad_11, trial_grad_12,
       trial_grad_20, trial_grad_21, trial_grad_22 = symbols('trial_grad_00 trial_grad_01 trial_grad_02 trial_grad_10 trial_grad_11 trial_grad_12 trial_grad_20 trial_grad_21 trial_grad_22')
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
       test_grad_00, test_grad_01, test_grad_02,
       test_grad_10, test_grad_11, test_grad_12,
       test_grad_20, test_grad_21, test_grad_22 = symbols('test_grad_00 test_grad_01 test_grad_02 test_grad_10 test_grad_11 test_grad_12 test_grad_20 test_grad_21 test_grad_22')
       test = Matrix(3,3,[
        test_grad_00, test_grad_01, test_grad_02,
        test_grad_10, test_grad_11, test_grad_12,
        test_grad_20, test_grad_21, test_grad_22])
    return test

def subs_gradient_trial(d,f):
    dim = 2
    ret = f;
    for k in range(0, dim):
        ret = ret.subs(f'trial_grad_{d}{k}', symbols(f'trial_g_{k}'))

    for k in range(0, dim):
        if k != d:
            for j in range(0, dim):
                ret = ret.subs(f'trial_grad_{k}{j}', 0)

    return ret

def subs_gradient_test(d,f):
    dim = 2
    ret = f;
    for k in range(0, dim):
        ret = ret.subs(f'test_grad_{d}{k}', symbols(f'test_g_{k}'))

    for k in range(0, dim):
        if k != d:
            for j in range(0, dim):
                ret = ret.subs(f'test_grad_{k}{j}', 0)

    return ret


def first_piola(strain_energy, F):
    shape = F.shape
    print(shape)

    P = Matrix.zeros(shape[0], shape[1])

    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            P[i, j] = diff(strain_energy, F[i, j])

    return P

