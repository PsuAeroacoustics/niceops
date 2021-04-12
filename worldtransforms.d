module worldtransforms;

import numd.linearalgebra.matrix;

import wgs84tocartesian;

import std.math;

class WorldTransform {
	
	immutable Vector!(3, double) reference_point_lat_lon_alt;

	this() {
		reference_point_lat_lon_alt = Vector!(3, double)(0, 0, 0);
	}

	this(Vector!(3, double) lat_long_alt_reference) {
		reference_point_lat_lon_alt = lat_long_alt_reference;
	}

	auto lat_long_alt_to_enu(Vector!(3, double) lat_long_alt_point) {
		double[2] lat_lon_ref = [reference_point_lat_lon_alt[0], reference_point_lat_lon_alt[1]];
		double[2] lat_lon = [lat_long_alt_point[0], lat_long_alt_point[1]];
		auto xy_point = toCartesian(lat_lon_ref, lat_lon);

		return Vector!(3, double)(
			xy_point[0],
			xy_point[1],
			lat_long_alt_point[2]
		);
	}
}
