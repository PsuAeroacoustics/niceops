"""Independent closed-form oracle for the wind_profile_is_its_law CTest: writes
tests/data/wind_profile_expected.csv (python3 tests/data/wind_profile_expected.py > that file).
The run is --wind 12kt,250,6.5 --roughness_length 0.1 --surface_layer_depth 3000 --profile_levels 10
--temperature 288.15 --humidity 40 --pressure 100 in feet (the default --units).  The surface layer
(3000 ft) reaches above 1.05 x the demo track's highest point, so the top is the surface layer's."""
import math
FT = 0.3048
z0, zr, depth, per = 0.1 * FT, 6.5 * FT, 3000 * FT, 10
speed = 12 * 1852 / 3600
d = math.radians(250)
east, north = -math.sin(d), -math.cos(d)  # blowing FROM 250: toward 70
T0, RH0, p0 = 288.15, 40.0, 100.0
lapse = 0.0098
g, R = 9.80665, 287.05


def psat(T):
    return 101.325 * 10 ** (-6.8346 * (273.16 / T) ** 1.261 + 4.6151)


h = RH0 * psat(T0) / p0
anchors = [10 * z0, zr, depth]
z = [0.0, anchors[0]]
for a, b in zip(anchors, anchors[1:]):
    n = max(1, math.ceil(per * math.log10(b / a) - 1e-9))
    z += [a * (b / a) ** (j / n) for j in range(1, n)] + [b]
print('z_ft,T_K,RH,p_kPa,u_mps,v_mps')
for zz in z:
    T = T0 - lapse * zz
    p = p0 * (T / T0) ** (g / (R * lapse))
    rh = RH0 if zz == 0 else min(100.0, h * p / psat(T))
    u = 0.0 if zz == 0 else speed * math.log(zz / z0) / math.log(zr / z0)
    print(','.join(repr(x) for x in (zz / FT, T, rh, p, u * east, u * north)))
