import math
import numpy as np  # Используем numpy для работы с массивами


def remainder(z, M):
    return z - M * math.trunc(z / M)


def multiplicative_congruential(a0, b, n, M):
    array1 = [float(a0)]
    array2 = []

    for i in range(n):
        array1.append(remainder(float(b) * array1[i], M))
        array2.append(array1[i + 1] / M)

    return array2


def linear_congruential(a0, a, c, n, M):
    array = [float(a0) / M]

    for i in range(1, n + 1):
        z = float(a * int(array[i - 1] * M) + c)
        z = z % M
        array.append(z / M)

    return array


def maclaren_marsaglia_method(a0, a, b, c, n, K, M):
    seqC = multiplicative_congruential(a0, b, n, M)
    seqB = linear_congruential(a0, a, c, n + K, M)

    seqA = []
    tableV = list(seqB)

    for i in range(n):
        s = math.trunc(seqC[i] * K)
        seqA.append(tableV[int(s)])
        newIndex = i + K
        tableV[int(s)] = seqB[newIndex]

    return seqA


def hi_criteria(seq, n):
    K = 1000
    intervals = [i / K for i in range(K + 1)]

    match_interval = [0] * K

    for j in range(n):
        for i in range(K):
            if intervals[i] <= seq[j] < intervals[i + 1]:
                match_interval[i] += 1

    expected = n / K
    Xtheor = sum((i - expected) ** 2 / expected for i in match_interval)

    Xtable = 1044.55  # Критическое значение для уровня значимости 0.05

    return Xtheor, Xtable, "H0" if Xtheor < Xtable else "H1"


def colmogorov_criteria(seq, n):
    K = 1000
    seq.sort()

    Feps = [i / len(seq) for i in range(K)]
    modules = []

    for i in range(2):
        buf = Feps[i] - seq[i]
        modules.append(abs(buf))

    sup = max(modules)
    sqrtN = math.sqrt(n)

    nDn = sup * sqrtN
    delta = 1.36 / sqrtN  # Критическое значение для уровня значимости 0.05

    return nDn, delta, "H0" if nDn < delta else "H1"


def main():
    M = 2 ** 31
    a0 = 24389  
    b = 79507
    K = 32  
    n = 1000
    a = 22695477
    c = 1

    # Генерация последовательностей
    marsaglia_method = maclaren_marsaglia_method(a0, a, b, c, n, K, M)  # Обратите внимание на добавление 'c'
    linear_gen1 = linear_congruential(a0, a, c, n, M)
    linear_gen2 = linear_congruential(a0 + 1, a, c, n, M)  # Пример второго генератора с изменением a0

    # Тест Хи-квадрат
    hi_statistic_marsaglia, hi_critical_marsaglia, hi_result_marsaglia = hi_criteria(marsaglia_method, n)
    hi_statistic_gen1, hi_critical_gen1, hi_result_gen1 = hi_criteria(linear_gen1, n)
    hi_statistic_gen2, hi_critical_gen2, hi_result_gen2 = hi_criteria(linear_gen2, n)

    # Критерий Колмогорова
    kolmogorov_statistic_marsaglia, kolmogorov_critical_marsaglia, kolmogorov_result_marsaglia = colmogorov_criteria(
        marsaglia_method, n)
    kolmogorov_statistic_gen1, kolmogorov_critical_gen1, kolmogorov_result_gen1 = colmogorov_criteria(linear_gen1, n)
    kolmogorov_statistic_gen2, kolmogorov_critical_gen2, kolmogorov_result_gen2 = colmogorov_criteria(linear_gen2, n)

    # Вывод значений для таблицы "Генераторы"
    print("Таблица Генераторы:")
    print(
        "| Генератор                         | M       | Генератор индексов массива | Генератор чисел из [0, 1] | K (длина массива) |")
    print(
        "|-----------------------------------|---------|----------------------------|----------------------------|--------------------|")
    print(
        f"| Генератор 1                       | {M}    | {a0}                      | {linear_gen1[:5]}         | {K}                 |")
    print(
        f"| Генератор 2                       | {M}    | {a0 + 1}                  | {linear_gen2[:5]}         | {K}                 |")
    print(
        f"| Генератор Маклорена-Марсальи     | {M}    | {a0}                      | {marsaglia_method[:5]}     | {K}                 |")

    # Вывод значений для таблицы "Тест Хи-квадрат"
    print("\nТаблица Тест Хи-квадрат:")
    print(
        "| Генератор                          | Значение статистики | Критическое значение | Уровень значимости | Принятая гипотеза |")
    print(
        "|------------------------------------|---------------------|----------------------|---------------------|--------------------|")
    print(
        f"| Генератор 1                       | {hi_statistic_gen1:.4f}           | {hi_critical_gen1:.4f}              | 0.05                | {hi_result_gen1}                 |")
    print(
        f"| Генератор 2                       | {hi_statistic_gen2:.4f}           | {hi_critical_gen2:.4f}              | 0.05                | {hi_result_gen2}                 |")
    print(
        f"| Генератор Маклорена-Марсальи     | {hi_statistic_marsaglia:.4f}           | {hi_critical_marsaglia:.4f}              | 0.05                | {hi_result_marsaglia}                 |")

    # Вывод значений для таблицы "Критерий согласия Колмогорова"
    print("\nТаблица Критерий согласия Колмогорова:")
    print(
        "| Генератор                          | Значение статистики | Критическое значение | Уровень значимости | Принятая гипотеза |")
    print(
        "|------------------------------------|---------------------|----------------------|---------------------|--------------------|")
    print(
        f"| Генератор 1                       | {kolmogorov_statistic_gen1:.4f} | {kolmogorov_critical_gen1:.4f} | 0.05                | {kolmogorov_result_gen1} |")
    print(
        f"| Генератор 2                       | {kolmogorov_statistic_gen2:.4f} | {kolmogorov_critical_gen2:.4f} | 0.05                | {kolmogorov_result_gen2} |")
    print(
        f"| Генератор Маклорена-Марсальи     | {kolmogorov_statistic_marsaglia:.4f} | {kolmogorov_critical_marsaglia:.4f} | 0.05                | {kolmogorov_result_marsaglia} |")


if __name__ == "__main__":
    main()
