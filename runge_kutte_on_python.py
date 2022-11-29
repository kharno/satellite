import math
import matplotlib.pyplot as plt
import numpy as np


def calculate_function_mmd(coord):
    # Earth's gravitational constant, metres/sec^2
    mu = 398600.44 * 1000 * 1000 * 1000
    r = math.sqrt(coord[3] ** 2 + coord[4] ** 2 + coord[5] ** 2)

    # movement model of satellite in central gravity field, metres
    func = [-mu / (r ** 3) * coord[3],
            -mu / (r ** 3) * coord[4],
            -mu / (r ** 3) * coord[5],
            coord[0],
            coord[1],
            coord[2]]

    return func


def calculate_move_runge(coord, h):
    result_coord = []
    # k1 = h * f(x)
    k1 = [h * k for k in calculate_function_mmd(coord)]
    # k2 = h * f(x + 0.5*k1)
    k2 = [h * k for k in calculate_function_mmd([a + b for a, b in zip(coord, [0.5 * k for k in k1])])]
    # k3 = h * f(x + 0.5*k2)
    k3 = [h * k for k in calculate_function_mmd([a + b for a, b in zip(coord, [0.5 * k for k in k2])])]
    # k4 = h * f(x + k3)
    k4 = [h * k for k in calculate_function_mmd([a + b for a, b in zip(coord, k3)])]
    for i in range(6):
        result_coord.append(coord[i] + (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]) / 6)

    return result_coord


def calculate_agesk_to_keo(agesk):
    # Earth gravitational constant 398600.44 km^2/sec^5
    mu = 398600.44 * 1000 * 1000 * 1000
    result_keo = [0 for item in range(11)]

    # проекции и модуль вектора площадей
    cx = agesk[4] * agesk[2] - agesk[5] * agesk[1]
    cy = agesk[5] * agesk[0] - agesk[3] * agesk[2]
    cz = agesk[3] * agesk[1] - agesk[4] * agesk[0]
    c = math.sqrt(cx * cx + cy * cy + cz * cz)
    r = math.sqrt(agesk[3] * agesk[3] + agesk[4] * agesk[4] + agesk[5] * agesk[5])

    # проекции и модуль интеграла энергии
    lx = -mu * agesk[3] / r + cz * agesk[1] - cy * agesk[2]
    ly = -mu * agesk[4] / r + cx * agesk[2] - cz * agesk[0]
    lz = -mu * agesk[5] / r + cy * agesk[0] - cx * agesk[1]
    l_ = math.sqrt(lx * lx + ly * ly + lz * lz)

    # наклонение, радианы
    result_keo[0] = math.acos(cz / c)

    # долгота восходящего узла, радианы
    if -cy / c / math.sin(result_keo[0]) > 0:
        result_keo[1] = math.fmod(math.asin(cx / c / math.sin(result_keo[0])) + 2 * math.pi, 2 * math.pi)
    else:
        result_keo[1] = math.pi - math.asin(cx / c / math.sin(result_keo[0]))

    # аргумент широты перигея, радианы
    if (cx * ly - cy * lx) / l_ / c / math.sin(result_keo[0]) > 0:
        result_keo[2] = math.fmod(math.asin(lz / l_ / math.sin(result_keo[0])) + 2 * math.pi, 2 * math.pi)
    else:
        result_keo[2] = math.pi - math.asin(lz / l_ / math.sin(result_keo[0]))

    # эксцентриситет
    result_keo[3] = l_ / mu

    # большая полуось, метры
    result_keo[4] = (c * c / mu) / (1 - result_keo[3] * result_keo[3])

    if (c * c - mu * r) / l_ / r > 0:
        # истинная аномалия КА
        result_keo[6] = math.fmod(
            math.asin((c * (agesk[3] * agesk[0] + agesk[4] * agesk[1] + agesk[5] * agesk[2])) / (l_ * r)) + 2 * math.pi,
            2 * math.pi)
        # аргумент широты КА
        result_keo[5] = math.fmod(result_keo[6] + result_keo[2], 2 * math.pi)
    else:
        # истинная аномалия КА
        result_keo[6] = math.pi - math.asin(
            (c * (agesk[3] * agesk[0] + agesk[4] * agesk[1] + agesk[5] * agesk[2])) / (l_ * r))
        # аргумент широты КА
        result_keo[5] = math.fmod(result_keo[6] + result_keo[2], 2 * math.pi)

    # фокальный параметр, м
    result_keo[7] = result_keo[4] * (1 - result_keo[3] * result_keo[3])

    # апогейное расстояние, м
    result_keo[8] = result_keo[7] / (1 - result_keo[3])

    # перигейное расстояние, м
    result_keo[9] = result_keo[7] / (1 + result_keo[3])

    # радиус-вектор КА, м
    result_keo[10] = result_keo[7] / (1 + result_keo[3] * math.cos(result_keo[6]))

    return result_keo


def main():
    # example
    h = 10
    # coordinates in AGSK for round orbit with altitude H = 500 km and zero inclination
    prev_coord = [0, 3454.842, 6787.935, 6871000, 0, 0]
    for iteration in range(10000):
        next_coord = calculate_move_runge(prev_coord, h)
        prev_coord = next_coord




if __name__ == "__main__":
    main()
