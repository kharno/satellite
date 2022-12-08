from math import pi, sqrt, sin, cos, asin, acos, atan, fmod
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from src.constants import GM, wz


class Coordinates:
    def __init__(self, coords, epoch, coord_system):
        if isinstance(coords, list):
            if len(coords) == 6:
                self.coords = coords
            else:
                raise ValueError('the number of coordinates must be six')
        else:
            raise ValueError('input list of six coordinates')
        if isinstance(epoch, datetime):
            self.epoch = epoch
        else:
            raise ValueError('the epoch must be datetime')
        if coord_system in ('keo', 'agesc', 'gsc'):
            self.coord_system = coord_system
        else:
            raise ValueError('input keo, agesc or gsc')

    def calculate_agesc_to_keo(self):

        if self.coord_system != 'agesc':
            raise ValueError('the coordinate system must be AGESC')

        keo = []

        # проекции и модуль вектора площадей
        cx = self.coords[4] * self.coords[2] - self.coords[5] * self.coords[1]
        cy = self.coords[5] * self.coords[0] - self.coords[3] * self.coords[2]
        cz = self.coords[3] * self.coords[1] - self.coords[4] * self.coords[0]
        c = sqrt(cx * cx + cy * cy + cz * cz)
        r = sqrt(self.coords[3] * self.coords[3] + self.coords[4] * self.coords[4] + self.coords[5] * self.coords[5])

        # проекции и модуль интеграла энергии
        lx = -GM * self.coords[3] / r + cz * self.coords[1] - cy * self.coords[2]
        ly = -GM * self.coords[4] / r + cx * self.coords[2] - cz * self.coords[0]
        lz = -GM * self.coords[5] / r + cy * self.coords[0] - cx * self.coords[1]
        l = sqrt(lx * lx + ly * ly + lz * lz)

        # keo[0]: inclination, rad
        keo.append(acos(cz / c))

        # keo[1]: longitude of the ascending node, rad
        if -cy / c / sin(keo[0]) > 0:
            keo.append(fmod(asin(cx / c / sin(keo[0])) + 2 * pi, 2 * pi))
        else:
            keo.append(pi - asin(cx / c / sin(keo[0])))

        # keo[2]: argument of periсenter, rad
        if (cx * ly - cy * lx) / l / c / sin(keo[0]) > 0:
            keo.append(fmod(asin(lz / l / sin(keo[0])) + 2 * pi, 2 * pi))
        else:
            keo.append(pi - asin(lz / l / sin(keo[0])))

        # keo[3]: eccentricity
        keo.append(l / GM)

        # koe[4]: semimajor axis, metres
        keo.append((c * c / GM) / (1 - keo[3] * keo[3]))

        if (c * c - GM * r) / l / r > 0:
            true_anomaly = fmod(
                asin((c * (self.coords[3] * self.coords[0] + self.coords[4] * self.coords[1] + self.coords[5] *
                           self.coords[2])) / (l * r)) + 2 * pi,
                2 * pi)
            # keo[5]: argument of satellite, rad
            keo.append(fmod(true_anomaly + keo[2], 2 * pi))
        else:
            true_anomaly = pi - asin(
                (c * (self.coords[3] * self.coords[0] + self.coords[4] * self.coords[1] + self.coords[5] * self.coords[
                    2])) / (l * r))
            # keo[5]: argument of satellite, rad
            keo.append(fmod(true_anomaly + keo[2], 2 * pi))

        self.coords = keo
        self.coord_system = 'keo'


class Satellite:

    def __init__(self, agesk=None, gsk=None, keo=None, dt=datetime(2022, 8, 26)):
        if agesk is not None:
            self.agesk = agesk
            self.keo = Satellite.calculate_agesk_to_keo(self.agesk)
            self.gsk = Satellite.calculate_agesk_to_gsk(self.agesk, dt)
        elif gsk is not None:
            self.gsk = gsk
            self.agesk = Satellite.calculate_gsk_to_agesk(self.gsk, dt)
            self.keo = Satellite.calculate_gsk_to_keo(self.gsk, dt)
        elif keo is not None:
            self.keo = keo
            self.agesk = Satellite.calculate_keo_to_agesk(keo)
            self.gsk = Satellite.calculate_keo_to_gsk(keo, dt)

    @staticmethod
    def calculate_move(coord, iteration_count=8700, inert=False):
        prev_coord = coord
        next_coord = []
        for iteration in range(iteration_count):
            next_coord = Satellite.calculate_move_runge(prev_coord, inert=inert)
            prev_coord = next_coord
        return next_coord

    @staticmethod
    def __calculate_function_mmd(coord, inert=False):

        r = sqrt(coord[3] ** 2 + coord[4] ** 2 + coord[5] ** 2)

        if inert:
            func = [
                -GM / (r ** 3) * coord[3],
                -GM / (r ** 3) * coord[4],
                -GM / (r ** 3) * coord[5],
                coord[0],
                coord[1],
                coord[2]
            ]
        else:
            func = [
                -GM / (r ** 3) * coord[3] + coord[3] * wz ** 2 + 2 * wz * coord[1],
                -GM / (r ** 3) * coord[4] + coord[4] * wz ** 2 - 2 * wz * coord[0],
                -GM / (r ** 3) * coord[5],
                coord[0],
                coord[1],
                coord[2]
            ]

        return func

    @staticmethod
    def calculate_move_runge(coord, h=10, inert=False):
        result_coord = []
        # k1 = h * f(x)
        k1 = [h * k for k in Satellite.__calculate_function_mmd(coord, inert=inert)]
        # k2 = h * f(x + 0.5*k1)
        k2 = [h * k for k in
              Satellite.__calculate_function_mmd([a + b for a, b in zip(coord, [0.5 * k for k in k1])], inert=inert)]
        # k3 = h * f(x + 0.5*k2)
        k3 = [h * k for k in
              Satellite.__calculate_function_mmd([a + b for a, b in zip(coord, [0.5 * k for k in k2])], inert=inert)]
        # k4 = h * f(x + k3)
        k4 = [h * k for k in Satellite.__calculate_function_mmd([a + b for a, b in zip(coord, k3)], inert=inert)]
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
        c = sqrt(cx * cx + cy * cy + cz * cz)
        r = sqrt(agesk[3] * agesk[3] + agesk[4] * agesk[4] + agesk[5] * agesk[5])

        # проекции и модуль интеграла энергии
        lx = -GM * agesk[3] / r + cz * agesk[1] - cy * agesk[2]
        ly = -GM * agesk[4] / r + cx * agesk[2] - cz * agesk[0]
        lz = -GM * agesk[5] / r + cy * agesk[0] - cx * agesk[1]
        l_ = sqrt(lx * lx + ly * ly + lz * lz)

        # наклонение, радианы
        result_keo[0] = acos(cz / c)

        # долгота восходящего узла, радианы
        if -cy / c / sin(result_keo[0]) > 0:
            result_keo[1] = fmod(asin(cx / c / sin(result_keo[0])) + 2 * pi, 2 * pi)
        else:
            result_keo[1] = pi - asin(cx / c / sin(result_keo[0]))

        # аргумент широты перигея, радианы
        if (cx * ly - cy * lx) / l_ / c / sin(result_keo[0]) > 0:
            result_keo[2] = fmod(asin(lz / l_ / sin(result_keo[0])) + 2 * pi, 2 * pi)
        else:
            result_keo[2] = pi - asin(lz / l_ / sin(result_keo[0]))

        # эксцентриситет
        result_keo[3] = l_ / GM

        # большая полуось, метры
        result_keo[4] = (c * c / GM) / (1 - result_keo[3] * result_keo[3])

        if (c * c - GM * r) / l_ / r > 0:
            # истинная аномалия КА
            result_keo[6] = fmod(
                asin((c * (agesk[3] * agesk[0] + agesk[4] * agesk[1] + agesk[5] * agesk[2])) / (l_ * r)) + 2 * pi,
                2 * pi)
            # аргумент широты КА
            result_keo[5] = fmod(result_keo[6] + result_keo[2], 2 * pi)
        else:
            # истинная аномалия КА
            result_keo[6] = pi - asin(
                (c * (agesk[3] * agesk[0] + agesk[4] * agesk[1] + agesk[5] * agesk[2])) / (l_ * r))
            # аргумент широты КА
            result_keo[5] = fmod(result_keo[6] + result_keo[2], 2 * pi)

        # фокальный параметр, м
        result_keo[7] = result_keo[4] * (1 - result_keo[3] * result_keo[3])

        # радиус-вектор КА, м
        result_keo[8] = result_keo[7] / (1 + result_keo[3] * cos(result_keo[6]))

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
            keo.append(keo[7] / (1 + keo[3] * cos(keo[6])))

        # нахождение радиальной и трансверсальной составляющей скорости
        vr = sqrt(GM / keo[7]) * keo[3] * sin(keo[6])
        vt = sqrt(GM / keo[7]) * (1 + keo[3] * cos(keo[6]))

        # нахождение скоростей в АГЭСК
        agesk[0] = vr * (cos(keo[1]) * cos(keo[5]) - sin(keo[1]) * sin(keo[5]) * cos(keo[0])) - vt * (
                cos(keo[1]) * sin(keo[5]) + sin(keo[1]) * cos(keo[5]) * cos(keo[0]))
        agesk[1] = vr * (sin(keo[1]) * cos(keo[5]) + cos(keo[1]) * sin(keo[5]) * cos(keo[0])) - vt * (
                sin(keo[1]) * sin(keo[5]) - cos(keo[1]) * cos(keo[5]) * cos(keo[0]))
        agesk[2] = vr * sin(keo[5]) * sin(keo[0]) + vt * cos(keo[5]) * sin(keo[0])

        # нахождение координат в АГЭСК
        agesk[3] = keo[8] * (cos(keo[1]) * cos(keo[5]) - sin(keo[1]) * sin(keo[5]) * cos(keo[0]))
        agesk[4] = keo[8] * (sin(keo[1]) * cos(keo[5]) + cos(keo[1]) * sin(keo[5]) * cos(keo[0]))
        agesk[5] = keo[8] * sin(keo[5]) * sin(keo[0])

        return agesk

    @staticmethod
    def calculate_agesk_to_gsk(agesk, dt):

        gsk = [0 for coord in range(6)]

        s = Satellite.calculate_star_time(dt)

        gsk[3] = agesk[3] * cos(s) + agesk[4] * sin(s)
        gsk[4] = agesk[4] * cos(s) - agesk[3] * sin(s)
        gsk[5] = agesk[5]
        gsk[0] = agesk[0] * cos(s) + agesk[1] * sin(s) + wz * gsk[4]
        gsk[1] = agesk[1] * cos(s) - agesk[0] * sin(s) - wz * gsk[3]
        gsk[2] = agesk[2]

        return gsk

    @staticmethod
    def calculate_gsk_to_agesk(gsk, dt):
        agesk = [0 for coord in range(6)]
        s = Satellite.calculate_star_time(dt)

        agesk[3] = gsk[3] * cos(s) - gsk[4] * sin(s)
        agesk[4] = gsk[4] * cos(s) + gsk[3] * sin(s)
        agesk[5] = gsk[5]
        agesk[0] = gsk[0] * cos(s) - gsk[1] * sin(s) - wz * agesk[4]
        agesk[1] = gsk[1] * cos(s) + gsk[0] * sin(s) + wz * agesk[3]
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
        mi = dt.minute + (dt.second - 1) / 60 + (dt.microsecond + 12) / 1000 / 60

        a = int((14 - m) / 12)
        yy = y + 4800 - a
        mm = m + 12 * a - 3

        jd = (d + int((153 * mm + 2) / 5) + 365 * yy + int(yy / 4) - int(yy / 100) + int(yy / 400) - 32045) + (
                h - 12) / 24 + mi / 1440
        t = (jd - 2451545) / 36525.
        s = ((6 * 3600 + 41 * 60 + 50.54841 + 8640184.812866 * t + 0.093104 * pow(t, 2) - 6.2 / 1000000 * pow(t,
                                                                                                              3)) + h * 3600 + mi * 60)

        return fmod(s * 2 * pi / 86400, pi * 2)

    @staticmethod
    def get_subsat_point(coord, dt):
        from src.constants import a, b
        s = Satellite.calculate_star_time(dt)
        longitude = asin(coord[4] / sqrt(coord[3] ** 2 + coord[4] ** 2))
        if coord[3] > 0 and coord[4] > 0:
            longitude = longitude
        elif coord[3] < 0 and coord[4] > 0:
            longitude = pi - longitude
        elif coord[3] < 0 and coord[4] < 0:
            longitude = pi - longitude
        elif coord[3] > 0 and coord[4] < 0:
            longitude = 2 * pi + longitude

        # coord = Satellite.calculate_gsk_to_agesk(coord, dt)
        # ksi = coord[3] * cos(s + longitude) + coord[4] * sin(s + longitude)
        # eta = coord[5]
        # latitude = atan(b**2/a**2 * eta/ksi)
        latitude = asin(coord[5] / sqrt(coord[3] ** 2 + coord[4] ** 2 + coord[5] ** 2))

        return latitude * 180 / pi, longitude * 180 / pi

    @staticmethod
    def get_satroute(start_coord, start_dt, end_dt):
        prev_coord = start_coord
        satroute_latitude = [Satellite.get_subsat_point(coord=start_coord, dt=start_dt)[0], ]
        satroute_longitude = [Satellite.get_subsat_point(coord=start_coord, dt=start_dt)[1], ]
        dt = start_dt
        while dt < end_dt:
            next_coord = Satellite.calculate_move_runge(prev_coord, inert=False)
            prev_coord = next_coord
            dt = dt + timedelta(seconds=10)
            satroute_latitude.append(Satellite.get_subsat_point(coord=prev_coord, dt=dt)[0])
            satroute_longitude.append(Satellite.get_subsat_point(coord=prev_coord, dt=dt)[1])

        return satroute_longitude, satroute_latitude


def main():
    # example
    h = 10
    # coordinates in AGESK for round orbit with altitude H = 500 km and zero inclination
    # sat = Satellite(keo=[63.7*pi/180, 287*pi/180, 271*pi/180, 0.72, 26000000, 0])
    # print(sat.agesk)
    # print(sat.gsk)
    # sat = Satellite(agesk=[0, 3454.842, 6787.935, 6871000, 0, 0])
    # sat = Satellite(gsk=[-1284.4554889715778, 2659.907832592989, 6787.935, 6187360.526329692, 2987843.824100776, 0])
    # # Satellite.calculate_move(sat.gsk, inert=False)
    # satroute = Satellite.get_satroute(start_coord=sat.gsk,
    #                              start_dt=datetime(2022, 8, 26),
    #                              end_dt=datetime(2022, 8, 26, 2))
    # # print(min(satroute[0]))
    # plt.plot(*satroute)
    # plt.show()
    coords = Coordinates(coords=[0, 3454.842, 6787.935, 6871000, 0, 0],
                         epoch=datetime(2022, 12, 1),
                         coord_system='agesc')
    coords.calculate_agesc_to_keo()
    print(coords.coords)


if __name__ == "__main__":
    main()
