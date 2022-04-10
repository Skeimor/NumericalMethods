import math
import numpy as np
import scipy.integrate as integrate


def f(x):
    return 1/math.sqrt(1+x**2)


def simpson(a, b):
    return ((b-a)/6)*(f(a)+4*f((a+b)/2)+f(b))


def complex_simpson(points):
    simpson_sum = 0
    for i in range(len(points)-1):
        simpson_sum += simpson(points[i], points[i + 1])
    return simpson_sum


def three_eighths(a, b):
    return ((b-a)/8)*(f(a)+3*f((2*a+b)/3)+3*f((a+2*b)/3)+f(b))


def complex_three_eighths(points):
    simpson_sum = 0
    for i in range(len(points)-1):
        simpson_sum += three_eighths(points[i], points[i + 1])
    return simpson_sum


def split_segment(a, b, step_size):
    points_count = round((b-a) / step_size)
    points = []
    current_point = a
    for i in range(points_count):
        points.append(current_point)
        current_point += step_size
    if points[len(points) - 1] != b:
        points.append(b)
    return points


def calculate_A_coeff(roots, k):
    roots = list(map(float, roots))
    multiplier_functions = []
    for i in range(len(roots)):
        if i == k:
            continue
        multiplier_functions.append(lambda x, root_i=roots[i], root_k=roots[k]: (x-root_i)/(root_k-root_i))

    function_to_integrate = lambda x: product_of_functions(x, multiplier_functions)
    result = integrate.quad(function_to_integrate, 0.5, 2)
    return result[0]


def product_of_functions(x, functions):
    product = 1
    for function in functions:
        product *= function(x)
    return product


def gauss_method():
    # Решаем систему относительно a1, a2, a3, указанную в отчёте (readme.md)
    A = np.array([[2.625, 1.875, 1.5], [3.98438, 2.625, 1.875], [6.39375, 3.98438, 2.625]])
    B = np.array([-3.98438, -6.39375, -10.6641])
    a1_a2_a3 = np.linalg.solve(A, B)

    omega_polynomial = list(a1_a2_a3)
    omega_polynomial.insert(0, 1)

    # Нашли точки, в которых будем считать значение функции
    roots = np.roots(omega_polynomial)

    A_coeffs = [calculate_A_coeff(roots, k) for k in range(len(roots))]

    integral_approximation = 0
    for i in range(len(roots)):
        integral_approximation += float(f(roots[i]))*float(A_coeffs[i])

    return integral_approximation


def main():
    points = split_segment(0.5, 2, 0.1)
    simpson_result1 = complex_simpson(points)
    print('Симпсон для шага 0.1:', simpson_result1)
    three_eighths_result1 = complex_three_eighths(points)
    print('3/8 для шага 0.1:', three_eighths_result1)
    print()

    points = split_segment(0.5, 2, 0.05)
    simpson_result2 = complex_simpson(points)
    print('Симпсон для шага 0.05:', simpson_result2)
    three_eighths_result2 = complex_three_eighths(points)
    print('3/8 для шага 0.05:', three_eighths_result2)
    print()

    points = split_segment(0.5, 2, 0.025)
    simpson_result3 = complex_simpson(points)
    print('Симпсон для шага 0.025:', simpson_result3)
    three_eighths_result3 = complex_three_eighths(points)
    print('3/8 для шага 0.025:', three_eighths_result3)
    print()

    print('Погрешность для Симпсона 0.1:', (simpson_result1-simpson_result2)/(2**4-1))
    print('Погрешность для Симпсона 0.05:', (simpson_result2-simpson_result3)/(2**4-1))
    print('Погрешность для Симпсона 0.025:', (simpson_result3-complex_simpson(split_segment(0.5, 2, 0.0125)))/(2**4-1))
    print()

    print('Погрешность для 3/8 0.1:', (three_eighths_result1-three_eighths_result2)/(2**5-1))
    print('Погрешность для 3/8 0.05:', (three_eighths_result2-three_eighths_result3)/(2**5-1))
    print('Погрешность для 3/8 0.025:', (three_eighths_result3-complex_three_eighths(split_segment(0.5, 2, 0.0125)))/(2**5-1))
    print()

    print('Квадратура Гаусса', gauss_method())


if __name__ == '__main__':
    main()
