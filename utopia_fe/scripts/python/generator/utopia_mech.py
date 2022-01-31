from sympy import symbols
from sympy import Matrix

def deformation_gradient(d):
    if d == 1:
        f00 = symbols('f00')
        F = Matrix(1,1,[f00])
    elif d == 2:
        f00, f01, f10, f11 = symbols('f00 f01 f10 f11')
        F = Matrix(2,2,[f00,f01,f10,f11])
    else:
        f00, f01, f02, f10, f11, f12, f20, f21, f22  = symbols('f00 f01 f02 f10 f11 f12 f20 f21 f22')
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

