import math

from layers import LayerParams, LayersRnNormal, LayersRnParallel


def test_layers2_normal():
    layers = []
    layers.append(LayerParams(eps=1.0))
    layers.append(LayerParams(eps=4.0))

    rn_normal = LayersRnNormal(layers)

    # Для двух слоев на всех частотах коэффициент отражения одинаковый
    # (по модулю)
    assert(
        abs(rn_normal.getR(2 * math.pi * 1.0e9)) - 1.0 / 3.0 < 1.0e-10)
    assert(
        abs(rn_normal.getR(2 * math.pi * 3.0e9)) - 1.0 / 3.0 < 1.0e-10)
    assert(
        abs(rn_normal.getR(2 * math.pi * 5.0e9)) - 1.0 / 3.0 < 1.0e-10)


def test_layers2_parallel():
    layers = []
    layers.append(LayerParams(eps=1.0))
    layers.append(LayerParams(eps=4.0))

    rn_parallel = LayersRnParallel(layers)

    # Для двух слоев на всех частотах коэффициент отражения одинаковый
    # (по модулю)
    assert(abs(rn_parallel.getR(
        2 * math.pi * 1.0e9)) - 1.0 / 3.0 < 1.0e-10)
    assert(abs(rn_parallel.getR(
        2 * math.pi * 3.0e9)) - 1.0 / 3.0 < 1.0e-10)
    assert(abs(rn_parallel.getR(
        2 * math.pi * 5.0e9)) - 1.0 / 3.0 < 1.0e-10)


def test_layers4_parallel():
    layers = []
    layers.append(LayerParams(eps=1.0))
    layers.append(LayerParams(eps=4.0, thickness=0.1))
    layers.append(LayerParams(eps=1.0, thickness=0.105))
    layers.append(LayerParams(eps=4.0))

    rn_parallel = LayersRnParallel(layers)

    assert(abs(rn_parallel.getR(
        2 * math.pi * 0.4688e9)) - 0.7362 < 1.0e-4)

    assert(abs(rn_parallel.getR(
        2 * math.pi * 0.7422e9)) - 0.3329 < 1.0e-4)
    assert(abs(rn_parallel.getR(
        2 * math.pi * 2.93e9)) - 0.3878 < 1.0e-4)


def test_layers4_normal():
    layers = []
    layers.append(LayerParams(eps=1.0))
    layers.append(LayerParams(eps=4.0, thickness=0.1))
    layers.append(LayerParams(eps=1.0, thickness=0.105))
    layers.append(LayerParams(eps=4.0))

    rn_normal = LayersRnNormal(layers)

    assert(
        abs(rn_normal.getR(2 * math.pi * 0.4688e9)) - 0.7362 < 1.0e-4)
    assert(
        abs(rn_normal.getR(2 * math.pi * 0.7422e9)) - 0.3329 < 1.0e-4)
    assert(
        abs(rn_normal.getR(2 * math.pi * 2.93e9)) - 0.3878 < 1.0e-4)
