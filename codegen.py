from sympy import *
from constant import constants

def eigen_type(M):
    if M.rows == 1 and M.cols == 1:
        return 'double'
    else:
        return 'Eigen::Matrix<double,%u,%u>' % (M.rows, M.cols)

def generate_eigen_matrix_func_proto(func_name, subx, sympy_matrix, params):
    paramlist = ['%s %s' % (eigen_type(Matrix(val)), name) for name,val in params]

    return '%s %s(%s);' % (eigen_type(sympy_matrix), func_name, ', '.join(paramlist))

def generate_eigen_matrix_func(func_name, subx, sympy_matrix, params):
    paramlist = ', '.join(['%s %s' % (eigen_type(Matrix(val)), name) for name,val in params])

    subs = []
    for name, val in params:
        val = Matrix(val)
        for r in range(val.rows):
            for c in range(val.cols):
                subs.append((val[r,c], Symbol('%s(%u,%u)' % (name, r, c))))
    subs = dict(subs)

    ret = '%s %s(%s) {\n' % (eigen_type(sympy_matrix), func_name, paramlist)


    for subx_name, subx_expr in subx:
        expr = subx_expr.xreplace(subs)
        ret += '    double %s = %s;\n' % (subx_name, ccode(expr, user_functions={'get_blade_thrust':'get_blade_thrust', 'get_blade_torque':'get_blade_torque'}).replace('\n',''))
    ret += '\n'
    ret += '    %s ret;\n' % (eigen_type(sympy_matrix),)
    for r in range(sympy_matrix.rows):
        for c in range(sympy_matrix.cols):
            expr = sympy_matrix[r,c].xreplace(subs)
            ret += '    ret(%u,%u) = %s;\n' % (r,c,ccode(expr, user_functions={'get_blade_thrust':'get_blade_thrust', 'get_blade_torque':'get_blade_torque'}).replace('\n',''))
    ret += '    return ret;\n'
    ret += '}'
    return ret


with open('output/derivation.srepr', 'rb') as f:
    derivation = sympify(f.read())

with open('output/quadcoptereqns.h', 'wb') as f:
    f.truncate()
    f.write('#pragma once\n')
    f.write('#include <Eigen/Dense>\n\n')
    for name, value in constants.items():
        f.write('static const double %s = %f;\n' % (name,value))
    f.write('\n')
    f.write(generate_eigen_matrix_func_proto('get_mm', derivation['subx'], derivation['mm'], [('states',derivation['states']), ('inputs',derivation['inputs'])]))
    f.write('\n')
    f.write(generate_eigen_matrix_func_proto('get_fo', derivation['subx'], derivation['fo'], [('states',derivation['states']), ('inputs',derivation['inputs'])]))
    f.write('\n')
    f.write(generate_eigen_matrix_func_proto('get_blade_v_z', [], derivation['blade_v_z'], [('states',derivation['states']), ('inputs',derivation['inputs'])]))
    f.write(generate_eigen_matrix_func_proto('get_blade_v_2', [], derivation['blade_v_2'], [('states',derivation['states']), ('inputs',derivation['inputs'])]))

with open('output/quadcoptereqns.cpp', 'wb') as f:
    f.truncate()
    f.write('#include \"quadcoptereqns.h\"\n')
    f.write('#include <specifieds.h>\n\n')
    f.write(generate_eigen_matrix_func('get_mm', derivation['subx'], derivation['mm'], [('states',derivation['states']), ('inputs',derivation['inputs'])]))
    f.write('\n\n')
    f.write(generate_eigen_matrix_func('get_fo', derivation['subx'], derivation['fo'], [('states',derivation['states']), ('inputs',derivation['inputs'])]))
    f.write('\n\n')        
    f.write(generate_eigen_matrix_func('get_blade_v_z', [], derivation['blade_v_z'], [('states',derivation['states']), ('inputs',derivation['inputs'])]))
    f.write(generate_eigen_matrix_func('get_blade_v_2', [], derivation['blade_v_2'], [('states',derivation['states']), ('inputs',derivation['inputs'])]))
