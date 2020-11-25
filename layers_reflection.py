# -*- coding: utf-8 -*-
"""Аналитический расчет коэффициента отражения от плоскослоистой среды"""

import numpy
import pylab

from layers import LayersRnNormal, LayerParams


if __name__ == '__main__':
    # Минимальная частота, Гц
    fmin = 0e9

    # Максимальная частота, Гц
    fmax = 1e9

    # Шаг расчета по частоте, Гц
    df = 0.005e9

    # Параметры плоскослоистой среды
    layers = []

    # Параметры среды, из которой падает плоская волна
    layers.append(LayerParams(eps=1.0))

    # Параметры слоев плоскослоистой среды
    layers.append(LayerParams(eps=6.3, thickness=0.31))
    layers.append(LayerParams(eps=2.7, thickness=0.34))
    # layers.append(layerParams(eps = 2.5, thickness = 0.20))
    layers.append(LayerParams(eps=9.7))

    freq = numpy.arange(0, fmax, df)
    w_list = freq * 2.0 * numpy.pi

    layers_list = [layers]

    for layers in layers_list:
        rn_normal = LayersRnNormal(layers)
        r = [rn_normal.getR(w, 0.0) for w in w_list]
        r_abs = numpy.abs(r)

    pylab.plot(freq * 1e-9, r_abs)
    pylab.grid(True)
    pylab.xlim(fmin * 1e-9, fmax * 1e-9)
    pylab.xlabel('f, ГГц')
    pylab.ylabel('|Г|')

    pylab.show()
