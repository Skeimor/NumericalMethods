import math


def f(x):
    return 1.2*math.cos(x) - math.e**x


def derivative_f(x):
    return -1.2*math.sin(x) - math.e**x


def second_derivative_f(x):
    return -1.2*math.cos(x) - math.e**x


def phi(x):
    return math.log(1.2 * math.cos(x))


class Solver:
    def __init__(self):
        self.epsilon = 0.5*10e-5
        self.segment = (0, 0.3)
        self.x_0 = 0.3

    def solve_half_division_method(self):
        iteration_count = math.ceil(math.log2((self.segment[1] - self.segment[0])/self.epsilon))
        print('Количество итераций:', iteration_count)
        (left, right) = self.segment
        for _ in range(iteration_count):
            middle = (left + right) / 2
            if f(left)*f(middle) < 0:
                right = middle
            elif f(middle)*f(right) < 0:
                left = middle
            elif f(left) == 0:
                return f(left)
            elif f(right) == 0:
                return f(right)
            elif f(middle) == 0:
                return f(middle)
            else:
                raise Exception('Impossible situation')

        return left, right

    def solve_newton_method(self):
        x_n = self.x_0
        counter = 0
        while True:
            counter += 1
            new_x_n = x_n - f(x_n)/derivative_f(x_n)
            if math.fabs(new_x_n - x_n) <= self.epsilon:
                print('Количество итераций:', counter)
                return new_x_n
            x_n = new_x_n

    def solve_modified_newton_method(self):
        x_n = self.x_0
        counter = 0
        while True:
            counter += 1
            new_x_n = x_n - f(x_n) / derivative_f(self.x_0)
            if math.fabs(new_x_n - x_n) <= self.epsilon:
                print('Количество итераций:', counter)
                return new_x_n
            x_n = new_x_n

    def solve_chord_method(self):
        x_n = self.segment[1] if self.x_0 == self.segment[0] else self.segment[0]
        counter = 0
        while True:
            counter += 1
            new_x_n = x_n - (f(x_n) * (x_n - self.x_0)) / (f(x_n) - f(self.x_0))
            if math.fabs(new_x_n - x_n) <= self.epsilon:
                print('Количество итераций:', counter)
                return new_x_n
            x_n = new_x_n

    def solve_movable_chord_method(self):
        x_n_previous = self.x_0
        x_n = self.segment[1] if self.x_0 == self.segment[0] else self.segment[0]
        counter = 0
        while True:
            counter += 1
            new_x_n = x_n - (f(x_n) * (x_n - x_n_previous)) / (f(x_n) - f(x_n_previous))
            if math.fabs(new_x_n - x_n) <= self.epsilon:
                print('Количество итераций:', counter)
                return new_x_n
            x_n_previous = x_n
            x_n = new_x_n

    def solve_simple_iteration_method(self):
        x_n = self.x_0
        counter = 0
        while True:
            counter += 1
            new_x_n = phi(x_n)
            if math.fabs(new_x_n - x_n) <= self.epsilon:
                print('Количество итераций:', counter)
                return new_x_n
            x_n = new_x_n


def main():
    solver = Solver()

    print('Обоснование выбора отрезка: f(left)*f(right) < 0')
    print(f(solver.segment[0])*f(solver.segment[1]))
    print()

    print('Обоснование выбора x_0: f(x_0)*f\'\'(x_0) > 0')
    print(f(solver.x_0)*second_derivative_f(solver.x_0))
    print()

    print('Метод половинного деления')
    print('Интервал', solver.solve_half_division_method())
    print()

    print('Метод Ньютона')
    print(solver.solve_newton_method())
    print()

    print('Модифицированный метод Ньютона')
    print(solver.solve_modified_newton_method())
    print()

    print('Метод неподвижных хорд')
    print(solver.solve_chord_method())
    print()

    print('Метод подвижных хорд')
    print(solver.solve_movable_chord_method())
    print()

    print('Метод простой итерации')
    print(solver.solve_simple_iteration_method())
    print()


if __name__ == '__main__':
    main()
