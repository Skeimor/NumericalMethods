import math


def solution_fits_required_accuracy(solution, old_solution, accuracy):
    """
    Сравнивает различие решений по евклидовой метрике
    """
    difference_sum = sum(((solution[i] - old_solution[i])**2 for i in range(len(solution))))
    return math.sqrt(difference_sum) < accuracy


def gauss_method(matrix):
    """
    Метод Гаусса
    """
    # Приводим к ступенчатому виду
    for main_line_number in range(len(matrix) - 1):
        for i in range(main_line_number + 1, len(matrix)):
            line = matrix[i]
            # Определяем множитель для строки, делим элемент текущей строки так,
            # чтобы занулить его при вычитании из нее другой строки для построения
            # ступенчатой матрицы.
            multiplier = line[main_line_number] / matrix[main_line_number][main_line_number]
            line[main_line_number] = 0
            for j in range(main_line_number + 1, len(line)):
                line[j] = line[j] - matrix[main_line_number][j] * multiplier

    # Осталось решить ступенчатую систему
    roots = []
    for i in range(len(matrix)-1, -1, -1):
        line = matrix[i]
        sum_by_main_matrix_line = 0
        for j in range(len(line) - 2, i, -1):
            sum_by_main_matrix_line += line[j]
        # Выразили переменную, соответствующую номеру строки
        x = (line[len(line) - 1] - sum_by_main_matrix_line) / line[i]
        roots.append(x)

        # Корретируем значение в предыдущей строке
        # Вместо коэффициента пишем коэффициент, домноженный на только что найденный корень
        if i == 0:
            continue
        matrix[i - 1][i] *= roots[len(roots) - 1]

    roots.reverse()
    return roots


def deep_clone_matrix(matrix):
    matrix_copy = []
    for i in range(len(matrix)):
        line = matrix[i].copy()
        matrix_copy.append(line)
    return matrix_copy


def gauss_with_main_element_method(matrix):
    """
    Метод Гаусса с выбором главного элемента
    """

    # Всю сортировку можно реализовать в более общем виде и намного эффективнее,
    # но это будет нецелесообразно для данной задачи, поэтому воспользуемся знанием,
    # что у нас матрица 3x3

    matrix_copy = deep_clone_matrix(matrix)

    first = max(matrix_copy, key=lambda line: abs(line[0]))
    matrix_copy.remove(first)
    second = max(matrix_copy, key=lambda line: abs(line[1]))
    matrix_copy.remove(second)
    third = max(matrix_copy, key=lambda line: abs(line[3]))

    new_matrix = [first, second, third]
    roots = gauss_method(new_matrix)
    return roots


def jacobi_method(matrix, epsilon=0.5*10e-4):
    """
    Метод Якоби
    """
    solutions = [tuple((1 for _ in range(len(matrix))))]

    n = 1
    while True:
        new_solution = []
        for i in range(len(matrix)):
            line = matrix[i]
            line_sum = line[len(line) - 1]
            for j in range(len(line) - 1):
                if i == j:
                    continue
                line_sum -= line[j] * solutions[len(solutions)-1][j]
            new_solution.append(line_sum / line[i])

        if solution_fits_required_accuracy(new_solution, solutions[len(solutions) - 1], epsilon):
            print('Число итераций: ', n)
            return new_solution

        solutions.append(tuple(new_solution))
        n += 1


def gauss_seidel_method(matrix, epsilon=0.5*10e-4):
    """
    Метод Гаусса-Зейделя
    """
    solutions = [tuple((1 for _ in range(len(matrix))))]

    n = 1
    while True:
        new_solution = [None for _ in range(len(matrix))]
        for i in range(len(matrix)):
            line = matrix[i]
            line_sum = line[len(line) - 1]
            for j in range(len(line) - 1):
                if i == j:
                    continue
                x_value = solutions[len(solutions) - 1][j] if new_solution[j] is None else new_solution[j]
                line_sum -= line[j] * x_value
            new_solution[i] = line_sum / line[i]

        if solution_fits_required_accuracy(new_solution, solutions[len(solutions) - 1], epsilon):
            print('Число итераций: ', n)
            return new_solution

        solutions.append(tuple(new_solution))
        n += 1


def main():
    # Полная матрица системы
    matrix = [
        [0.2, -0.08, 1.62, 2.06],
        [0.24, 0.88, -0.10, 2.38],
        [-0.50, -0.32, -0.16, -2.30]
    ]
    print('Метод Гаусса')
    roots = gauss_method(deep_clone_matrix(matrix))
    print(roots)
    print()

    print('Метод Гаусса с выбором главного элемента')
    roots = gauss_with_main_element_method((deep_clone_matrix(matrix)))
    print(roots)
    print()

    print('Метод Якоби')
    jacobi_matrix = deep_clone_matrix(matrix)
    # Меняем строчки местами, чтобы метод сходился
    jacobi_matrix[0], jacobi_matrix[2] = jacobi_matrix[2], jacobi_matrix[0]
    roots = jacobi_method(deep_clone_matrix(jacobi_matrix))
    print(roots)
    print()

    print('Метод Гаусса-Зейделя')
    roots = gauss_seidel_method(deep_clone_matrix(jacobi_matrix))
    print(roots)
    print()


if __name__ == '__main__':
    main()
