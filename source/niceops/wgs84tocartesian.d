module wgs84tocartesian;

import std.math;

/++
 +	Ported with love from https://github.com/chrberger/WGS84toCartesian
 +/

double[2] toCartesian(immutable(double[2]) WGS84Reference, immutable(double[2]) WGS84Position) {
    enum double DEG_TO_RAD = PI / 180.0;
    enum double PI_2 = PI / 2.0;
    enum double EPSILON10 = 1.0e-10;
    enum double EPSILON12 = 1.0e-12;

    enum double EQUATOR_RADIUS = 6378137.0;
    enum double FLATTENING = 1.0 / 298.257223563;
    enum double SQUARED_ECCENTRICITY = 2.0 * FLATTENING - FLATTENING * FLATTENING;
    enum double SQUARE_ROOT_ONE_MINUS_ECCENTRICITY = 0.996647189335;
    enum double POLE_RADIUS = EQUATOR_RADIUS * SQUARE_ROOT_ONE_MINUS_ECCENTRICITY;

    enum double C00 = 1.0;
    enum double C02 = 0.25;
    enum double C04 = 0.046875;
    enum double C06 = 0.01953125;
    enum double C08 = 0.01068115234375;
    enum double C22 = 0.75;
    enum double C44 = 0.46875;
    enum double C46 = 0.01302083333333333333;
    enum double C48 = 0.00712076822916666666;
    enum double C66 = 0.36458333333333333333;
    enum double C68 = 0.00569661458333333333;
    enum double C88 = 0.3076171875;

    enum double R0 = C00 - SQUARED_ECCENTRICITY * (C02 + SQUARED_ECCENTRICITY * (C04 + SQUARED_ECCENTRICITY * (C06 + SQUARED_ECCENTRICITY * C08)));
    enum double R1 = SQUARED_ECCENTRICITY * (C22 - SQUARED_ECCENTRICITY * (C04 + SQUARED_ECCENTRICITY * (C06 + SQUARED_ECCENTRICITY * C08)));
    enum double R2T = SQUARED_ECCENTRICITY * SQUARED_ECCENTRICITY;
    enum double R2 = R2T * (C44 - SQUARED_ECCENTRICITY * (C46 + SQUARED_ECCENTRICITY * C48));
    enum double R3T = R2T * SQUARED_ECCENTRICITY;
    enum double R3 = R3T * (C66 - SQUARED_ECCENTRICITY * C68);
    enum double R4 = R3T * SQUARED_ECCENTRICITY * C88;

    double mlfn(immutable double lat) {
        immutable double sin_phi = sin(lat);
        immutable double cos_phi = cos(lat) * sin_phi;
        immutable double squared_sin_phi = sin_phi * sin_phi;
        return (R0 * lat - cos_phi * (R1 + squared_sin_phi * (R2 + squared_sin_phi * (R3 + squared_sin_phi * R4))));
    }

    immutable double ML0 = mlfn(WGS84Reference[0] * DEG_TO_RAD);

    auto msfn = (immutable double sinPhi, immutable double cosPhi, immutable double es) { return (cosPhi / sqrt(1.0 - es * sinPhi * sinPhi)); };

    double[2] project(double lat, double lon) {
        double[2] retVal = [lon, -1.0 * ML0];
        if (!(abs(lat) < EPSILON10)) {
            immutable double ms = (abs(sin(lat)) > EPSILON10) ? msfn(sin(lat), cos(lat), SQUARED_ECCENTRICITY) / sin(lat) : 0.0;
            retVal[0] = ms * sin(lon *= sin(lat));
            retVal[1] = (mlfn(lat) - ML0) + ms * (1.0 - cos(lon));
        }
        return retVal;
    }

    double[2] fwd(double lat, double lon) {
        immutable double D = abs(lat) - PI_2;
        if ((D > EPSILON12) || (abs(lon) > 10.0)) {
            //return std::array<double, 2>{0.0, 0.0};
			return [0.0, 0.0];
        }
        if (abs(D) < EPSILON12) {
            lat = (lat < 0.0) ? -1.0 * PI_2 : PI_2;
        }
        lon -= WGS84Reference[1] * DEG_TO_RAD;
        auto projectedRetVal = project(lat, lon);
        //return std::array<double, 2>{EQUATOR_RADIUS * projectedRetVal[0], EQUATOR_RADIUS * projectedRetVal[1]};
		return [EQUATOR_RADIUS * projectedRetVal[0], EQUATOR_RADIUS * projectedRetVal[1]];
    }

    return fwd(WGS84Position[0] * DEG_TO_RAD, WGS84Position[1] * DEG_TO_RAD);
}


/**
 * @param WGS84Reference WGS84 position to be used as reference.
 * @param CartesianPosition Cartesian position to be transformed.
 * @return std::array<double, 2> Approximating a WGS84 position from a given CartesianPosition based on a given WGS84Reference using Mercator projection.
 */
double[2] fromCartesian(immutable(double[2]) WGS84Reference, immutable(double[2]) CartesianPosition) {
    enum double EPSILON10 = 1.0e-2;
    enum double incLon = 1e-5;
    const int signLon = (CartesianPosition[0] < 0) ? -1 : 1;
    enum double incLat = incLon;
    const int signLat = (CartesianPosition[1] < 0) ? -1 : 1;

    double[2] approximateWGS84Position = WGS84Reference;
    double[2] cartesianResult = toCartesian(WGS84Reference, approximateWGS84Position);

	double dPrev = double.max;
    double d = abs(CartesianPosition[1] - cartesianResult[1]);

    while ((d < dPrev) && (d > EPSILON10)) {
        approximateWGS84Position[0] = approximateWGS84Position[0] + signLat * incLat;
        cartesianResult             = toCartesian(WGS84Reference, approximateWGS84Position);
        dPrev                       = d;
        d                           = abs(CartesianPosition[1] - cartesianResult[1]);
    }

	dPrev = double.max;
    d     = abs(CartesianPosition[0] - cartesianResult[0]);
    while ((d < dPrev) && (d > EPSILON10)) {
        approximateWGS84Position[1] = approximateWGS84Position[1] + signLon * incLon;
        cartesianResult             = toCartesian(WGS84Reference, approximateWGS84Position);
        dPrev                       = d;
        d                           = abs(CartesianPosition[0] - cartesianResult[0]);
    }

    return approximateWGS84Position;
}
