# -*- coding: utf-8 -*-
"""Классы для расчета коэффициента отражения от слоистой среды"""

from abc import ABCMeta, abstractmethod
import math
import cmath


class LayerParams:
    """Параметры одного слоя"""

    def __init__(self, eps=1.0, mu=1.0, sigma=0.0, thickness=0.0):
        self.eps = eps
        self.mu = mu
        self.sigma = sigma
        self.thickness = thickness


class LayersR(metaclass=ABCMeta):
    """Общие формулы для вычисления коэффициента отражения"""

    @abstractmethod
    def _calcRBody(toplayer, w):
        pass

    @abstractmethod
    def _calcZ(self, z, ksi):
        pass

    @abstractmethod
    def _calcR2(self, z1, z2):
        """Коэффициент отражения между двумя слоями"""

    def getR(self, w, angle=0.0):
        """ Посчитаем все параметры для слоев"""
        if w == 0.0:
            return 0.0

        self._calc_ksi_gamma(w, angle)
        return self._calcRBody(0, w)

    def _calc_ksi_gamma(self, w, angle):
        self.ksi = []
        self.gamma = []
        self.z = []			# Волновые сопротивления
        self.r = []			# Коэффициенты отражения между слоями

        c = 299792458.0
        wavelength = c / (w / (2.0 * math.pi))

        self.ksi.append(angle)

        for n in range(len(self.layers)):
            layer = self.layers[n]

            currgamma = self._getGamma(w, layer.eps, layer.mu, layer.sigma)
            self.gamma.append(currgamma)

            if n != 0:
                self.ksi.append(cmath.asin(
                    self.gamma[n - 1] / currgamma * cmath.sin(self.ksi[n - 1])))

            self.z.append(self._calcZ(
                self._zwl(layer.eps, layer.sigma, wavelength), self.ksi[n]))

            if n != 0:
                self.r.append(self._calcR2(self.z[n - 1], self.z[n]))

    def _getGamma(self, w, eps, mu, sigma):
        """Посчитать gamma = alpha + i * beta"""
        eps0 = 1.0 / (36.0 * math.pi) * 1.0e-9
        mu0 = (4.0 * math.pi) * 1.0e-7

        tmp = sigma / (eps * eps0 * w)

        alpha = math.sqrt((w * w) * eps * eps0 * mu * mu0 *
                          (-0.5 + 0.5 * math.sqrt(1 + tmp * tmp)))
        beta = math.sqrt(w * w * eps * eps0 * mu * mu0 *
                         (0.5 + 0.5 * math.sqrt(1 + tmp * tmp)))

        return complex(alpha, beta)

    def _zwl(self, eps, sigma, wavelength):
        tmp = 60.0 * wavelength * sigma / eps
        y = math.sqrt(1.0 + tmp * tmp)
        b = math.sqrt(2.0 * eps * (1.0 + y))

        if y >= -1.0:
            a = math.sqrt(2.0 * eps * (-1.0 + y))
        else:
            a = 0.0

        return(240.0 * math.pi) / (complex(b, -a))


class LayersRnNormal(LayersR):
    def __init__(self, layers):
        self.layers = layers

    def _calcZ(self, z, ksi):
        return z / cmath.cos(ksi)

    def _calcR2(self, z1, z2):
        """Коэффициент отражения между двумя слоями"""
        return (z2 - z1) / (z1 + z2)

    def _calcRBody(self, topLayer, w):
        """Расчет общего коэффициента отражения"""

        if len(self.layers) - topLayer == 2:
            # Если два слоя
            res = self._calcR2(self.z[topLayer], self.z[topLayer + 1])
        else:
            # Коэффициент отражения от всех более нижних слоев
            Rbottom = self._calcRBody(topLayer + 1, w)

            tmp = (Rbottom * cmath.exp(-2.0 * self.gamma[topLayer + 1] *
                                       self.layers[topLayer + 1].thickness *
                                       cmath.cos(self.ksi[topLayer + 1])))

            res = (self.r[topLayer] + tmp) / (1.0 + self.r[topLayer] * tmp)
        return res


class LayersRnParallel(LayersR):
    def __init__(self, layers):
        self.layers = layers

    def _calcZ(self, z, ksi):
        return z * cmath.cos(ksi)

    def _calcR2(self, z1, z2):
        """Коэффициент отражения между двумя слоями"""
        return -(z2 - z1) / (z1 + z2)

    def _calcRBody(self, topLayer, w):
        """Расчет общего коэффициента отражения"""

        if len(self.layers) - topLayer == 2:
            # Если два слоя
            res = self._calcR2(self.z[topLayer], self.z[topLayer + 1])
        else:
            # Коэффициент отражения от всех более нижних слоев
            Rbottom = self._calcRBody(topLayer + 1, w)

            tmp = Rbottom * \
                cmath.exp(-2.0 * self.gamma[topLayer + 1] *
                          self.layers[topLayer + 1].thickness *
                          cmath.cos(self.ksi[topLayer + 1]))

            res = (self.r[topLayer] + tmp) / (1.0 + self.r[topLayer] * tmp)
        return res
