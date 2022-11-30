import math
import matplotlib.pyplot as plt
from datetime import datetime


class Satellite:
    # Earth's gravitational constant, metres^3/sec^2
    GM = 398600.44e+9
    # Earth's rotation velocity, rad/sec
    wz = 2*math.pi/86164

    def __init__(self, agesk=None, gsk=None, keo=None, dt=datetime(2022, 8, 26)):
        if agesk is not None:
            self.agesk = agesk
            self.keo = Satellite.calculate_agesk_to_keo(self.agesk)
            self.gsk = Satellite.calculate_agesk_to_gsk(self.agesk, dt)
        elif gsk is not None:
            self.gsk = gsk
            self.agesk_coord = Satellite.calculate_gsk_to_agesk(self.gsk, dt)
            self.keo = Satellite.calculate_gsk_to_keo(self.gsk, dt)
        elif keo is not None:
            self.keo = keo
            self.agesk = Satellite.calculate_keo_to_agesk(keo)
            self.gsk = Satellite.calculate_gsk_to_keo(keo, dt)

    def calculate_move(self, iteration_count=8700):
        prev_coord = self.agesk_coord
        x = []
        for iteration in range(iteration_count):
            next_coord = self.calculate_move_runge(prev_coord)
            prev_coord = next_coord
            x.append(prev_coord[3])
        plt.plot(x)
        plt.show()
        self.agesk_coord = prev_coord

    @staticmethod
    def __calculate_function_mmd(coord, inert=True):

        r = math.sqrt(coord[3] ** 2 + coord[4] ** 2 + coord[5] ** 2)

        if inert:
            func = [
                -Satellite.GM / (r ** 3) * coord[3],
                -Satellite.GM / (r ** 3) * coord[4],
                -Satellite.GM / (r ** 3) * coord[5],
                coord[0],
                coord[1],
                coord[2]
            ]
        else:
            func = [
                -Satellite.GM / (r ** 3) * coord[3] + coord[3] * Satellite.wz**2 * 2*Satellite.wz*coord[1],
                -Satellite.GM / (r ** 3) * coord[4] + coord[4] * Satellite.wz**2 - 2*Satellite.wz*coord[0],
                -Satellite.GM / (r ** 3) * coord[5],
                coord[0],
                coord[1],
                coord[2]
            ]

        return func

    @staticmethod
    def calculate_move_runge(coord, h=10):
        result_coord = []
        # k1 = h * f(x)
        k1 = [h * k for k in Satellite.__calculate_function_mmd(coord)]
        # k2 = h * f(x + 0.5*k1)
        k2 = [h * k for k in Satellite.__calculate_function_mmd([a + b for a, b in zip(coord, [0.5 * k for k in k1])])]
        # k3 = h * f(x + 0.5*k2)
        k3 = [h * k for k in Satellite.__calculate_function_mmd([a + b for a, b in zip(coord, [0.5 * k for k in k2])])]
        # k4 = h * f(x + k3)
        k4 = [h * k for k in Satellite.__calculate_function_mmd([a + b for a, b in zip(coord, k3)])]
        for i in range(6):
            result_coord.append(coord[i] + (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]) / 6)

        return result_coord

    @staticmethod
    def calculate_agesk_to_keo(agesk):
        # Earth gravitational constant 398600.44 km^2/sec^5
        result_keo = [0 for coord in range(11)]

        # проекции и модуль вектора площадей
        cx = agesk[4] * agesk[2] - agesk[5] * agesk[1]
        cy = agesk[5] * agesk[0] - agesk[3] * agesk[2]
        cz = agesk[3] * agesk[1] - agesk[4] * agesk[0]
        c = math.sqrt(cx * cx + cy * cy + cz * cz)
        r = math.sqrt(agesk[3] * agesk[3] + agesk[4] * agesk[4] + agesk[5] * agesk[5])

        # проекции и модуль интеграла энергии
        lx = -Satellite.GM * agesk[3] / r + cz * agesk[1] - cy * agesk[2]
        ly = -Satellite.GM * agesk[4] / r + cx * agesk[2] - cz * agesk[0]
        lz = -Satellite.GM * agesk[5] / r + cy * agesk[0] - cx * agesk[1]
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
        result_keo[3] = l_ / Satellite.GM

        # большая полуось, метры
        result_keo[4] = (c * c / Satellite.GM) / (1 - result_keo[3] * result_keo[3])

        if (c * c - Satellite.GM * r) / l_ / r > 0:
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

        # радиус-вектор КА, м
        result_keo[8] = result_keo[7] / (1 + result_keo[3] * math.cos(result_keo[6]))

        # апогейное расстояние, м
        result_keo[9] = result_keo[7] / (1 - result_keo[3])

        # перигейное расстояние, м
        result_keo[10] = result_keo[7] / (1 + result_keo[3])

        return result_keo

    @staticmethod
    def calculate_keo_to_agesk(keo):

        agesk = [None for coord in range(6)]

        if len(keo) == 6:
            keo.append(keo[5] - keo[2])
            keo.append(keo[4] * (1 - keo[3] * keo[3]))
            keo.append(keo[7] / (1 + keo[3] * math.cos(keo[6])))

        # нахождение радиальной и трансверсальной составляющей скорости
        vr = math.sqrt(Satellite.GM/keo[7])*keo[3]*math.sin(keo[6])
        vt = math.sqrt(Satellite.GM/keo[7])*(1+keo[3]*math.cos(keo[6]))

        # нахождение скоростей в АГЭСК
        agesk[0] = vr*(math.cos(keo[1])*math.cos(keo[5]) - math.sin(keo[1])*math.sin(keo[5])*math.cos(keo[0])) - vt*(math.cos(keo[1])*math.sin(keo[5])+math.sin(keo[1])*math.cos(keo[5])*math.cos(keo[0]))
        agesk[1] = vr*(math.sin(keo[1])*math.cos(keo[5]) + math.cos(keo[1])*math.sin(keo[5])*math.cos(keo[0])) - vt*(math.sin(keo[1])*math.sin(keo[5])-math.cos(keo[1])*math.cos(keo[5])*math.cos(keo[0]))
        agesk[2] = vr*math.sin(keo[5])*math.sin(keo[0]) + vt*math.cos(keo[5])*math.sin(keo[0])

        # нахождение координат в АГЭСК
        agesk[3] = keo[8] * (math.cos(keo[1])*math.cos(keo[5])-math.sin(keo[1])*math.sin(keo[5])*math.cos(keo[0]))
        agesk[4] = keo[8] * (math.sin(keo[1])*math.cos(keo[5])+math.cos(keo[1])*math.sin(keo[5])*math.cos(keo[0]))
        agesk[5] = keo[8] * math.sin(keo[5])*math.sin(keo[0])

        return agesk

    @staticmethod
    def calculate_agesk_to_gsk(agesk, dt):

        gsk = [0 for coord in range(6)]

        s = Satellite.calculate_star_time(dt)

        gsk[3] = agesk[3]*math.cos(s)+agesk[4]*math.sin(s)
        gsk[4] = agesk[4]*math.cos(s)-agesk[3]*math.sin(s)
        gsk[5] = agesk[5]
        gsk[0] = agesk[0]*math.cos(s)+agesk[1]*math.sin(s)+Satellite.wz*gsk[4]
        gsk[1] = agesk[1]*math.cos(s)-agesk[0]*math.sin(s)-Satellite.wz*gsk[3]
        gsk[2] = agesk[2]

        return gsk

    @staticmethod
    def calculate_gsk_to_agesk(gsk, dt):
        agesk = [0 for coord in range(6)]
        s = Satellite.calculate_star_time(dt)

        agesk[3] = gsk[3] * math.cos(s) - gsk[4] * math.sin(s)
        agesk[4] = gsk[4] * math.cos(s) + gsk[3] * math.sin(s)
        agesk[5] = gsk[5]
        agesk[0] = gsk[0] * math.cos(s) - gsk[1] * math.sin(s) - Satellite.wz * agesk[4]
        agesk[1] = gsk[1] * math.cos(s) + gsk[0] * math.sin(s) + Satellite.wz * agesk[3]
        agesk[2] = gsk[2]

        return agesk

    @staticmethod
    def calculate_keo_to_gsk(keo, dt):
        return Satellite.calculate_agesk_to_gsk(Satellite.calculate_keo_to_agesk(keo), dt)

    @staticmethod
    def calculate_gsk_to_keo(gsk, dt):
        return Satellite.calculate_agesk_to_keo(Satellite.calculate_gsk_to_agesk(gsk, dt))

    @staticmethod
    def calculate_star_time(dt):

        d = dt.day
        m = dt.month
        y = dt.year
        h = dt.hour
        mi = dt.minute + (dt.second - 1) / 60 + (dt.microsecond + 12)/1000/60

        a = int((14 - m) / 12)
        yy = y + 4800 - a
        mm = m + 12 * a - 3

        jd = (d + int((153 * mm + 2) / 5) + 365 * yy + int(yy / 4) - int(yy / 100) + int(yy / 400) - 32045) + (h - 12) / 24 + mi / 1440
        t = (jd - 2451545) / 36525.
        s = ((6 * 3600 + 41 * 60 + 50.54841 + 8640184.812866 * t + 0.093104 * pow(t, 2) - 6.2 / 1000000 * pow(t, 3)) + h * 3600 + mi * 60)

        return math.fmod(s * 2 * math.pi / 86400, math.pi * 2)


def main():
    # example
    h = 10
    # coordinates in AGESK for round orbit with altitude H = 500 km and zero inclination
    # sat = Satellite(keo=[62.7*math.pi/180, 287*math.pi/180, 271*math.pi/180, 0.72, 26000000, 0])
    sat = Satellite(agesk=[0, 3454.842, 6787.935, 6871000, 0, 0])
    print(sat.keo)
    print(sat.gsk)



if __name__ == "__main__":
    main()
