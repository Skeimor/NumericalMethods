# Лабораторные по численным методам

## lab2 (Вариант 12)

### Обоснование выбора φ(x)
![Обоснование выбора φ](lab2.png)

Все остальные обоснования представлены в выводе программы в коде.
```
Обоснование выбора отрезка: f(left)*f(right) < 0
-0.040691004125055164

Обоснование выбора x_0: f(x_0)*f''(x_0) > 0
0.5078771576555405

Метод половинного деления
Количество итераций: 13
Интервал (0.1680908203125, 0.16812744140625002)

Метод Ньютона
Количество итераций: 4
0.1681219754789549

Модифицированный метод Ньютона
Количество итераций: 6
0.16812475684379047

Метод неподвижных хорд
Количество итераций: 5
0.16811976205572976

Метод подвижных хорд
Количество итераций: 5
0.16812197548582009

Метод простой итерации
Количество итераций: 6
0.16812603886100677
```

Для каждого метода выводится количество итераций (к пункту о скорости сходимости методов).

## lab3 (Вариант 12)

### Условие
![Условие задачи](lab3.png)

### Сравнение результатов метода Гаусса и метода Гаусса с выбором главного элемента
Получились следующие корни:
```
Метод Гаусса
[3.0000000000000004, 2.0000000000000004, 1.0000000000000002]

Метод Гаусса с выбором главного элемента
[2.9999999999999996, 2.0, 1.0]
```
У второго метода 2 корня получились более точно, чем у первого. 
Оба достаточно близко к реальному решению: (3, 2, 1)

### Cходимость метода Якоби (пункт 3)
Воспользуемся достаточным условием сходимости метода Якоби (теорема 4.7)

Из него видно, что в исходной матрице условие строгого диагонального преобладания
не выполняется для 3-й строчки: `0.5 + 0.32 = 0.82 > 0.16`. Поэтому нельзя быть
уверенными в том, что метод сойдётся. 

Поменяем 1 и 3 строчки местами. Теперь легко заметить, что условие выполняется
для всех строчек. Следовательно, метод будет сходиться для любого начального x.
```
0.32 + 0.16 = 0.48 < 0.5
0.24 + 0.1 = 0.34 < 0.88
0.2 + 0.08 = 0.28 < 1.62
```
### Cходимость метода Гаусса-Зейделя (пункт 4)

Рассуждения абсолютно аналогичны рассуждениям для метода Якоби.

### Сравнение количества итераций метода Якоби и метода Гаусса-Зейделя

```
Метод Якоби
Число итераций:  12
[2.999853346812849, 1.999992421725982, 1.0000016683831427]

Метод Гаусса-Зейделя
Число итераций:  8
[3.000042789918389, 1.9999856991515688, 0.9999940110792886]
```
Как видно из результатов, метод Гаусса-Зейделя сходится быстрее.
