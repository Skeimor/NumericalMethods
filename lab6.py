import math
import matplotlib.pyplot as plt

"""
Граничные условия:
y(0) = 0
y(1) = e + 1/e -2
"""
y_0 = 0
y_1 = math.e + 1/math.e - 2


def second_derivative_y(x, y):
    alpha = 3.2  # 2 + 0.1 * 12
    return y + 2 * alpha + 2 + alpha * x * (1 - x)


def exact_y(x):
    alpha = 3.2
    return alpha*x*x - alpha*x + math.exp(-x) + math.exp(x) - 2


def shooting_method(points, step_size):
    """
    Метод стрельбы
    """
    mu_0 = -5
    mu_1 = 1
    phi_0, _ = solve_euler_with_recalculation(points, step_size, mu_0)
    phi_1, _ = solve_euler_with_recalculation(points, step_size, mu_1)
    closest_mu, y_values = solve_half_division_method((mu_0, mu_1), 10e-9, phi_0, phi_1, points, step_size)
    return y_values


def solve_euler_with_recalculation(points, step_size, mu):
    z_values = [mu]
    y_values = [y_0]

    for x in points[1:]:
        y_with_overline = y_values[-1] + step_size * z_values[-1]

        special_z = z_values[-1] + (step_size / 2) * (
                    second_derivative_y(points[-1], y_values[-1]) + second_derivative_y(x, y_with_overline))

        recalculated_y = y_values[-1] + (step_size / 2) * (z_values[-1] + special_z)
        recalculated_z = z_values[-1] + (step_size / 2) * (
                    second_derivative_y(points[-1], y_values[-1]) + second_derivative_y(x, recalculated_y))

        z_values.append(recalculated_z)
        y_values.append(recalculated_y)

    # phi(mu) = y(b, mu) - y(1) = 0

    return y_values[-1] - y_1, y_values


def solve_half_division_method(segment, epsilon, phi_0, phi_1, points, step_size):
    iteration_count = math.ceil(math.log2((segment[1] - segment[0])/epsilon))
    (left, right) = segment
    f_left = phi_0
    f_right = phi_1

    for _ in range(iteration_count):
        middle = (left + right) / 2
        f_middle, y_values = solve_euler_with_recalculation(points, step_size, middle)
        if f_left*f_middle < 0:
            right = middle
            f_right = f_middle
        elif f_middle*f_right < 0:
            left = middle
            f_left = f_middle
        elif f_left == 0:
            return left
        elif f_right == 0:
            return right
        elif f_middle == 0:
            return middle
        else:
            raise Exception('Impossible situation')

    return middle, y_values


def sweep_method(points, step_size):
    """
    Метод прогонки
    """
    check_diagonal_dominance(points, step_size)
    lambdas, mu_values = find_sweep_coefficients(points, step_size)

    # Обратная прогонка
    y_values = [float('nan') for _ in range(len(points))]
    y_values[-1] = y_1
    y_values[0] = y_0
    for i in range(1, len(points) - 1):
        y = lambdas[-i]*y_values[-i] + mu_values[-i]
        y_values[-i-1] = y

    return y_values


def check_diagonal_dominance(points, step_size):
    """
    Проверяет наличие диагонального преобладания в матрице системы для метода прогонки
    """
    for x in points:
        if abs(2+step_size**2*p(x)) <= 2:
            raise Exception("Диагональное преобладание не выполняется")


def find_sweep_coefficients(points, step_size):
    lambdas = [0]
    mu_values = [y_0]

    # Прямая прогонка
    for x_i_plus_1 in points[1:]:
        A = 2 + (step_size ** 2) * p(x_i_plus_1 - step_size)
        B = (step_size ** 2) * q(x_i_plus_1 - step_size)
        lambda_i_plus_1 = 1 / (A - lambdas[-1])
        mu_i_plus_1 = (mu_values[-1] - B) / (A - lambdas[-1])

        lambdas.append(lambda_i_plus_1)
        mu_values.append(mu_i_plus_1)

    return lambdas, mu_values


def p(x):
    """
    Функция p для метода прогонки
    """
    return 1


def q(x):
    """
    Функция q для метода прогонки
    """
    alpha = 3.2
    return 2 * alpha + 2 + alpha * x * (1 - x)


def main():
    show_results(10)
    # show_results(20)


def show_results(N):
    points, step_size = split_segment(0, 1, N)
    sweep_solution = sweep_method(points, step_size)
    shooting_solution = shooting_method(points, step_size)
    exact_solution = list(map(exact_y, points))

    plt.plot(points, shooting_solution, label="Метод стрельбы")
    plt.plot(points, sweep_solution, label="Метод прогонки")
    plt.plot(points, exact_solution, label="Точное решение")
    plt.legend()
    plt.show()


def split_segment(a, b, points_count):
    step_size = (b - a) / points_count
    points = []
    current_point = a
    for i in range(points_count):
        points.append(current_point)
        current_point += step_size
    points.append(current_point)

    if abs(points[len(points) - 1] - b) > 10e-9:
        points.append(b)
    return points, step_size


if __name__ == '__main__':
    main()
