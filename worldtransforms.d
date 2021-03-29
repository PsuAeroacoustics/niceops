module worldtransforms;

import numd.linearalgebra.matrix;

import std.math;

class WorldTransform {
	
	enum double EarthMajorAxis = 6_378_137.0;
	enum double EarthMinorAxis = 6_356_752.3;
	enum double MinorMajorRatioSquared = (EarthMinorAxis*EarthMinorAxis)/(EarthMajorAxis*EarthMajorAxis);
	enum double EarthEcentricity = 1.0 - MinorMajorRatioSquared;

	immutable Matrix!(3, 3, double) rot_mat;
	immutable Vector!(3, double) reference_point;

	this() {
		reference_point = Vector!(3, double)(0, 0, 0);
		rot_mat = Matrix!(3, 3, double).identity;
	}

	this(Vector!(3, double) lat_long_alt_reference) {
		//lat_long_alt_reference[2] /= 1000.0;
		reference_point = lat_lon_alt_to_ecef(lat_long_alt_reference);
		rot_mat = generate_rot_mat(lat_long_alt_reference);
	}

	static double surface_height(double latitude) {
		immutable lat_sin = sin(latitude);
		return EarthMajorAxis/sqrt(1.0 - EarthEcentricity*lat_sin*lat_sin);
	}

	static double surface_height_derivative(double latitude) {
		immutable lat_sin = sin(latitude);
		return EarthMajorAxis*EarthEcentricity*sin(latitude)*cos(latitude)/pow(1.0 - EarthEcentricity*lat_sin*lat_sin, 1.5);
	}

	static Vector!(3, double) lat_lon_alt_to_ecef(Vector!(3, double) lat_long_alt) {
		immutable cos_lat = cos(lat_long_alt[0]);
		immutable surf_height = surface_height(lat_long_alt[0]);
		immutable r = (surf_height + lat_long_alt[2]);
		return Vector!(3, double)(
			r*cos_lat*cos(lat_long_alt[1]),
			r*cos_lat*sin(lat_long_alt[1]),
			(MinorMajorRatioSquared*surf_height + lat_long_alt[2])*sin(lat_long_alt[1])
		);
	}

	static Matrix!(3, 3, double) generate_rot_mat(Vector!(3, double) lat_long_alt_reference) {
		immutable sin_lat = sin(lat_long_alt_reference[0]);
		immutable cos_lat = cos(lat_long_alt_reference[0]);
		immutable sin_long = sin(lat_long_alt_reference[1]);
		immutable cos_long = cos(lat_long_alt_reference[1]);

		return Matrix!(3, 3, double)(
			-sin_long, cos_long, 0,
			-sin_lat*cos_long, -sin_lat*sin_long, cos_lat,
			cos_lat*cos_long, cos_lat*sin_long, sin_lat
		);
	}

	auto lat_long_alt_to_enu(Vector!(3, double) lat_long_alt_point) {
		//lat_long_alt_point[2] /= 1000.0; // To preserve some precision we go to km
		auto point = lat_lon_alt_to_ecef(lat_long_alt_point);

		auto xyz = rot_mat*(point - reference_point);
		//xyz *= 1000.0; // Back to m
		return xyz;
		//return rot_mat*(point - reference_point);
	}

	auto enu_to_ecef(Vector!(3, double) xyz) {
		return rot_mat.transpose*xyz + reference_point;
	}

	auto enu_to_lat_lon(Vector!(3, double) xyz) {
		xyz /= 1000.0;
		auto XYZ = enu_to_ecef(xyz);

		// Ferrari's solution (Heikkinen)

		immutable double r = sqrt(XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1]);
		immutable double e_p_2 = (EarthMajorAxis*EarthMajorAxis - EarthMinorAxis*EarthMinorAxis)/(EarthMinorAxis*EarthMinorAxis);
		immutable double F = 54.0*EarthMinorAxis*EarthMinorAxis*XYZ[2]*XYZ[2];
		immutable double G = r*r + (1.0 - EarthEcentricity)*XYZ[2]*XYZ[2] - EarthEcentricity*(EarthMajorAxis*EarthMajorAxis - EarthMinorAxis*EarthMinorAxis);
		immutable double c = EarthEcentricity*EarthEcentricity*F*r*r/(G*G*G);
		immutable double s = pow(1.0 + c + sqrt(c*c + 2.0*c), 1.0/3.0);
		immutable double s_p = (s + 1.0 + 1.0/s);
		immutable double P = F/(3.0*s_p*s_p*G*G);
		immutable double Q = sqrt(1.0 + 2.0*EarthEcentricity*EarthEcentricity*P);
		immutable double r_0 = -P*EarthEcentricity*r/(1.0 + Q) + sqrt(0.5*EarthMajorAxis*EarthMajorAxis*(1.0 + 1.0/Q) - P*(1.0 - EarthEcentricity)*XYZ[2]*XYZ[2]/(Q*(1.0 + Q)) - 0.5*P*r*r);		
		immutable double r_p = (r - EarthEcentricity*r_0);
		immutable double U = sqrt(r_p*r_p + XYZ[2]*XYZ[2]);
		immutable double V = sqrt(r_p*r_p + (1.0 - EarthEcentricity)*XYZ[2]*XYZ[2]);
		immutable double z_0 = EarthMinorAxis*EarthMinorAxis*XYZ[2]/(EarthMajorAxis*V);

		Vector!(3, double) lat_lon_alt;
		lat_lon_alt[2] = U*(1.0 - EarthMinorAxis*EarthMinorAxis/(EarthMajorAxis*V))*1000.0;
		lat_lon_alt[0] = atan((XYZ[2] + e_p_2*z_0)/r);
		lat_lon_alt[1] = atan2(XYZ[1], XYZ[0]);
		return lat_lon_alt;
	}

	auto lat_long_alt_rate_to_enu_rate(Vector!(3, double) lat_long_alt, Vector!(3, double) lat_long_alt_rate) {
		//lat_long_alt_rate[2] /= 1000.0;
		//lat_long_alt[2] /= 1000.0;
		immutable sin_lat = sin(lat_long_alt[0]);
		immutable cos_lat = cos(lat_long_alt[0]);
		immutable sin_long = sin(lat_long_alt[1]);
		immutable cos_long = cos(lat_long_alt[1]);
		immutable N_hat = surface_height_derivative(lat_long_alt[0]);
		immutable N = surface_height(lat_long_alt[0]);

		auto V = Vector!(3, double)(
			(N_hat*lat_long_alt_rate[0] + lat_long_alt_rate[2])*cos_lat*cos_long - (N + lat_long_alt[2])*(lat_long_alt_rate[0]*sin_lat*cos_long + lat_long_alt_rate[1]*cos_lat*sin_long),
			(N_hat*lat_long_alt_rate[0] + lat_long_alt_rate[2])*cos_lat*sin_long + (N + lat_long_alt[2])*(lat_long_alt_rate[1]*cos_lat*cos_long - lat_long_alt_rate[0]*sin_lat*sin_long),
			(MinorMajorRatioSquared*N_hat*lat_long_alt_rate[0] + lat_long_alt_rate[2])*sin_lat + (MinorMajorRatioSquared*N + lat_long_alt[2])*lat_long_alt_rate[0]*cos_lat
		);

		//V *= 1000.0;

		return rot_mat*V;
	}
}
