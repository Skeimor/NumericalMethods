import math
import matplotlib.pyplot as plt


def derivative_y(x, y):
    """
    y(0) = 1/2
    """
    return y ** 3 / 2


def exact_y(x):
    """
    Точное решение задачи Коши
    """
    return 1 / math.sqrt(4 - x)


def solve_explicit_euler(points, step_size):
    """
    Явный метод Эйлера
    """
    # Из начального условия y(0)=2
    y_values = [0.5]
    current_y = 0.5
    for x in points[1:]:
        current_y = current_y + derivative_y(x, current_y) * step_size
        y_values.append(current_y)

    return y_values


def solve_runge_kutta_fourth_order(points, step_size):
    """
    Метод Рунге-Кутта четвертого порядка
    """
    y_values = [0.5]
    current_y = 0.5
    for x in points[1:]:
        K1 = step_size * derivative_y(x, current_y)
        K2 = step_size * derivative_y(x + step_size / 2, current_y + K1 / 2)
        K3 = step_size * derivative_y(x + step_size / 2, current_y + K2 / 2)
        K4 = step_size * derivative_y(x + step_size, current_y + K3)
        current_y = current_y + (1 / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
        y_values.append(current_y)

    return y_values


def solve_explicit_adams_three_step(points, step_size):
    """
    Явный метод Адамса (трехшаговый)
    """
    # Для разгона используем метод Рунге-Кутта
    # Вычисляем первые 3 значения, т.к. в дальшейшем метод использует 3 предыдущих значения на каждом шаге
    x_values = points[:3]
    y_values = solve_runge_kutta_fourth_order(x_values, step_size)

    # Метод Адамса
    for x in points[3:]:
        current_y = (y_values[-1] + (step_size / 12 *
                                     (23 * derivative_y(x_values[-1], y_values[-1])
                                      - 16 * derivative_y(x_values[-2], y_values[-2])
                                      + 5 * derivative_y(x_values[-3], y_values[-3]))))

        x_values.append(x)
        y_values.append(current_y)

    return y_values


def get_exact_solution(points):
    return list(map(lambda x: exact_y(x), points))


def main():
    show_results(10)
    show_results(20)
    show_results(30)


def show_results(N):
    """
    :param N: Количество шагов
    """
    points, step_size = split_segment(0, 1, N)
    euler_solution = solve_explicit_euler(points, step_size)
    print("N = " + str(N))
    print("X", points)
    print("Y Эйлера", euler_solution)
    plt.plot(points, euler_solution)
    plt.show()

    runge_kutta_solution = solve_runge_kutta_fourth_order(points, step_size)
    print("Y Рунге-Кутта", runge_kutta_solution)
    plt.plot(points, runge_kutta_solution)
    plt.show()

    adams_solution = solve_explicit_adams_three_step(points, step_size)
    print("Y Адамса", adams_solution)
    plt.plot(points, adams_solution)
    plt.show()

    print("Y Точного решения", get_exact_solution(points))
    plt.plot(points, get_exact_solution(points))
    plt.show()

    print()


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
