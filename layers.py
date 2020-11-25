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


class LayersCalculator(metaclass=ABCMeta):
    """Общие формулы для вычисления коэффициентов отражения"""

    def getRE(self, w, angle=0.0):
        """Рассчитать коэффициент отражения по полю E для многослойной среды"""
        if w == 0.0:
            return 0.0

        self.ksi = self._calc_ksi(w, angle)
        self.z = self._calc_z(w, self.ksi)
        self.RE = self._calc_RE(self.z)
        return self._calcREBody(0, w)

    def getRH(self, w, angle=0.0):
        """Рассчитать коэффициент отражения по полю H для многослойной среды"""
        return self.getRE(w, angle)

    @abstractmethod
    def _calcREBody(toplayer, w):
        """Расчет коэффициента отражения по полю E"""

    @abstractmethod
    def _calcZ(self, z, ksi):
        pass

    @abstractmethod
    def _calcRE2(self, z1, z2):
        """Коэффициент отражения по полю E между двумя слоями"""

    def _calc_ksi(self, w, angle):
        ksi = []
        ksi.append(angle)

        prev_gamma = None
        curr_gamma = None

        for n in range(len(self.layers)):
            layer = self.layers[n]
            # Коэффициент распространения для текущего слоя
            curr_gamma = self._getGamma(w, layer.eps, layer.mu, layer.sigma)

            if n != 0:
                ksi.append(cmath.asin(
                    prev_gamma / curr_gamma * cmath.sin(ksi[n - 1])))

            # Коэффициент распространения для предыдущего слоя
            prev_gamma = curr_gamma

        return ksi

    def _calc_z(self, w, ksi):
        z = []		# Волновые сопротивления

        c = 299792458.0
        wavelength = c / (w / (2.0 * math.pi))

        for n in range(len(self.layers)):
            layer = self.layers[n]
            z.append(self._calcZ(
                self._zwl(layer.eps, layer.sigma, wavelength), ksi[n]))

        return z

    def _calc_RE(self, z):
        return [self._calcRE2(z[n - 1], z[n])
                for n in range(1, len(self.layers))]

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

        return (240.0 * math.pi) / (complex(b, -a))


class LayersCalculatorNormal(LayersCalculator):
    def __init__(self, layers):
        self.layers = layers

    def _calcZ(self, z, ksi):
        return z / cmath.cos(ksi)

    def _calcRE2(self, z1, z2):
        """Коэффициент отражения по полю E между двумя слоями"""
        return (z2 - z1) / (z1 + z2)

    def _calcREBody(self, topLayerNumber, w):
        """Расчет общего коэффициента отражения по полю E"""

        if len(self.layers) - topLayerNumber == 2:
            # Если два слоя
            res = self._calcRE2(
                self.z[topLayerNumber], self.z[topLayerNumber + 1])
        else:
            # Коэффициент отражения от всех более нижних слоев
            Rbottom = self._calcREBody(topLayerNumber + 1, w)
            bottomLayer = self.layers[topLayerNumber + 1]
            gamma = self._getGamma(w, bottomLayer.eps,
                                   bottomLayer.mu, bottomLayer.sigma)

            tmp = (Rbottom * cmath.exp(-2.0 * gamma *
                                       self.layers[topLayerNumber + 1].thickness *
                                       cmath.cos(self.ksi[topLayerNumber + 1])))

            res = ((self.RE[topLayerNumber] + tmp) /
                   (1.0 + self.RE[topLayerNumber] * tmp))
        return res


class LayersCalculatorParallel(LayersCalculator):
    def __init__(self, layers):
        self.layers = layers

    def _calcZ(self, z, ksi):
        return z * cmath.cos(ksi)

    def _calcRE2(self, z1, z2):
        """Коэффициент отражения по полю E между двумя слоями"""
        return -(z2 - z1) / (z1 + z2)

    def _calcREBody(self, topLayerNumber, w):
        """Расчет общего коэффициента отражения по полю E"""

        if len(self.layers) - topLayerNumber == 2:
            # Если два слоя
            res = self._calcRE2(
                self.z[topLayerNumber], self.z[topLayerNumber + 1])
        else:
            # Коэффициент отражения от всех более нижних слоев
            Rbottom = self._calcREBody(topLayerNumber + 1, w)
            bottomLayer = self.layers[topLayerNumber + 1]
            gamma = self._getGamma(w, bottomLayer.eps,
                                   bottomLayer.mu, bottomLayer.sigma)

            tmp = (Rbottom *
                   cmath.exp(-2.0 * gamma *
                             self.layers[topLayerNumber + 1].thickness *
                             cmath.cos(self.ksi[topLayerNumber + 1])))

            res = ((self.RE[topLayerNumber] + tmp) /
                   (1.0 + self.RE[topLayerNumber] * tmp))
        return res
