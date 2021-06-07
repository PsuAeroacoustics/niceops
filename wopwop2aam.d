module wopwop2aam;

import plt = matplotlibd.pyplot;

import netcdf;

import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.getopt;
import std.file;
import std.math;
import std.path;
import std.range;
import std.stdio;
import std.string;
import std.typecons;

double safe_log10(double x) {
    if (x == 0) {
        return -double.infinity;
    } 
	if (x < 0) {
        return double.nan;
	} else {
        return log10(x);
	}
}

struct Dataset {
	float mach_number;
	float flight_path_angle;
	float tilt_angle;
	float sphere_radius;
	Channel[] channels;

	@nogc uint hash() {
		import std.digest.murmurhash : MurmurHash3;
		import std.bitmanip : read, write, Endian;
		ubyte[12] buff;

		buff[].write!float(mach_number, 0);
		buff[].write!float(flight_path_angle, 4);
		buff[].write!float(tilt_angle, 8);

		auto hashfunc = MurmurHash3!(uint.sizeof*8, 64)();

		hashfunc.put(buff);
		
		auto _hash = hashfunc.finish;
		// well aint this the dumbest shit. Need for read to work
		ubyte[] hash = _hash[];

		return hash.read!uint;
	}
}

struct Function {
	float[] data;
}

struct Channel {
	
	float[] frequency;
	Function[] functions;
}

auto read_dataset(string dir) {
	int[] header_buffer = new int[4];

	auto spectrum_file = File(dir.buildPath("spl_spectrum.fn"), "rb");
	scope(exit) spectrum_file.close;
	
	float[] flight_condition = readText(dir.buildPath("flight_condition.txt")).split(",").map!(a => a.strip.to!float).array;

	Dataset dataset;

	dataset.mach_number = flight_condition[0];
	dataset.flight_path_angle = flight_condition[1];
	dataset.tilt_angle = flight_condition[2];
	dataset.sphere_radius = flight_condition[3];

	//writeln("mach number: ", mach_number);
	//writeln("flight path angle: ", flight_path_angle);
	//writeln("tilt angle: ", tilt_angle);
	//writeln("sphere radius: ", sphere_radius);

	header_buffer = spectrum_file.rawRead(header_buffer);

	int i_max = header_buffer[0];
	int j_max = header_buffer[1];
	int samples = header_buffer[2];
	int functions = header_buffer[3];

	//writeln("i_max: ", i_max);
	//writeln("j_max: ", j_max);
	//writeln("samples: ", samples);
	//writeln("functions: ", functions);

	size_t total_obs = i_max*j_max;

	dataset.channels = new Channel[total_obs];

	foreach(ref channel; dataset.channels) {
		channel.frequency = new float[samples];
		channel.functions = new Function[functions - 1];
		foreach(ref func; channel.functions) {
			func.data = new float[samples];
		}
	}

	float[] sample_buff = new float[1];

	foreach(s_idx; 0..samples) {
		foreach(ref channel; dataset.channels) {
			sample_buff = spectrum_file.rawRead(sample_buff);
			channel.frequency[s_idx] = sample_buff[0];
		}
	}

	foreach(f_idx; 0..(functions - 1)) {
		foreach(s_idx; 0..samples) {
			foreach(ref channel; dataset.channels) {
				sample_buff = spectrum_file.rawRead(sample_buff);
				channel.functions[f_idx].data[s_idx] = sample_buff[0];
			}
		}
	}

	return dataset;
}

void to_aam_sphere(ref Dataset dataset, string dir, size_t idx, bool hemisphere = false) {
	
	float[] thetas = iota(0., 185., 5.).map!(a => a.to!float).array;
	float[] phis = iota(-180., 180., 10.).map!(a => a.to!float).array;

	float[] amplitude;

	size_t channel_idx = 0;

	if(!hemisphere) {
		amplitude = new float[thetas.length*phis.length*dataset.channels[0].frequency.length];
		foreach(p_idx, phi; phis) {
			foreach(t_idx, theta; thetas) {
				foreach(f_idx; 0..dataset.channels[0].frequency.length) {
					size_t amp_idx = f_idx + t_idx*dataset.channels[0].frequency.length + p_idx*dataset.channels[0].frequency.length*thetas.length;
					amplitude[amp_idx] = dataset.channels[channel_idx].functions[2].data[f_idx];
				}
				channel_idx++;
			}
		}
	} else {
		float[] real_phis;
		foreach(phi; phis) {
			if((phi >= -90) && (phi <= 90)) {
				real_phis ~= phi;
			}
		}
		amplitude = new float[thetas.length*real_phis.length*dataset.channels[0].frequency.length];

		size_t rp_idx = 0;
		foreach(p_idx, phi; phis) {
			foreach(t_idx, theta; thetas) {
				if((phi >= -90) && (phi <= 90)) {
					foreach(f_idx; 0..dataset.channels[0].frequency.length) {
						size_t amp_idx = f_idx + t_idx*dataset.channels[0].frequency.length + rp_idx*dataset.channels[0].frequency.length*thetas.length;
						amplitude[amp_idx] = dataset.channels[channel_idx].functions[2].data[f_idx];
					}
				}
				channel_idx++;
			}

			if((phi >= -90) && (phi <= 90)) {
				rp_idx++;
			}
		}

		phis = real_phis;
	}

	import std.format : format;

	auto filename = format!("%s%06d.nc")("sphere_", idx);

	int file_handle;
	auto nc_ret = nc_create(dir.buildPath(filename).toStringz, NC_WRITE, &file_handle);
	enforce(nc_ret == NC_NOERR, "Failed to create netCDF file");
	scope(exit) nc_close(file_handle);

	int phi_dimid;
	nc_ret = nc_def_dim(file_handle, "PHI_DIM", phis.length, &phi_dimid);
	enforce(nc_ret == NC_NOERR, "Failed to define phi dimension");

	int theta_dimid;
	nc_ret = nc_def_dim(file_handle, "THETA_DIM", thetas.length, &theta_dimid);
	enforce(nc_ret == NC_NOERR, "Failed to define theta dimension");

	int frequency_dimid;
	nc_ret = nc_def_dim(file_handle, "FREQUENCY_DIM", dataset.channels[0].frequency.length, &frequency_dimid);
	enforce(nc_ret == NC_NOERR, "Failed to define frequency dimension");

	int mach_dimid;
	nc_ret = nc_def_dim(file_handle, "MACH_DIM", 1, &mach_dimid);
	enforce(nc_ret == NC_NOERR, "Failed to define Mach number dimension");

	int fpa_dimid;
	nc_ret = nc_def_dim(file_handle, "FLIGHT_PATH_ANGLE_DIM", 1, &fpa_dimid);
	enforce(nc_ret == NC_NOERR, "Failed to define flight path angle dimension");

	int tilt_dimid;
	nc_ret = nc_def_dim(file_handle, "TILT_ANGLE_DIM", 1, &tilt_dimid);
	enforce(nc_ret == NC_NOERR, "Failed to define tilt angle dimension");

	int radius_dimid;
	nc_ret = nc_def_dim(file_handle, "RADIUS_DIM", 1, &radius_dimid);
	enforce(nc_ret == NC_NOERR, "Failed to define radius dimension");

	int phi_varid;
	nc_ret = nc_def_var(file_handle, "PHI", NC_FLOAT, 1, &phi_dimid, &phi_varid);
	enforce(nc_ret == NC_NOERR, "Failed to define phi variable");

	int theta_varid;
	nc_ret = nc_def_var(file_handle, "THETA", NC_FLOAT, 1, &theta_dimid, &theta_varid);
	enforce(nc_ret == NC_NOERR, "Failed to define theta variable");

	int mach_varid;
	nc_ret = nc_def_var(file_handle, "MACH", NC_FLOAT, 1, &mach_dimid, &mach_varid);
	enforce(nc_ret == NC_NOERR, "Failed to define Mach number variable");

	int fpa_varid;
	nc_ret = nc_def_var(file_handle, "FLIGHT_PATH_ANGLE", NC_FLOAT, 1, &fpa_dimid, &fpa_varid);
	enforce(nc_ret == NC_NOERR, "Failed to define flight path angle variable");

	int tilt_varid;
	nc_ret = nc_def_var(file_handle, "TILT_ANGLE", NC_FLOAT, 1, &tilt_dimid, &tilt_varid);
	enforce(nc_ret == NC_NOERR, "Failed to define tilt angle variable");

	int radius_varid;
	nc_ret = nc_def_var(file_handle, "RADIUS", NC_FLOAT, 1, &radius_dimid, &radius_varid);
	enforce(nc_ret == NC_NOERR, "Failed to define radius variable");

	int amplitude_varid;
	int[3] amplitude_dimids = [phi_dimid, theta_dimid, frequency_dimid];
	nc_ret = nc_def_var(file_handle, "AMPLITUDE", NC_FLOAT, amplitude_dimids.length, amplitude_dimids.ptr, &amplitude_varid);
	enforce(nc_ret == NC_NOERR, "Failed to define amplitude variable");

	int frequency_varid;
	nc_ret = nc_def_var(file_handle, "FREQUENCY", NC_FLOAT, 1, &frequency_dimid, &frequency_varid);
	enforce(nc_ret == NC_NOERR, "Failed to define frequency variable");

	nc_ret = nc_enddef(file_handle);
	enforce(nc_ret == NC_NOERR, "Failed to exit define mode");

	nc_ret = nc_put_var_float(file_handle, radius_varid, &dataset.sphere_radius);
	enforce(nc_ret == NC_NOERR, "Failed to put sphere radius with error "~nc_ret.to!string);
	nc_ret = nc_put_var_float(file_handle, mach_varid, &dataset.mach_number);
	enforce(nc_ret == NC_NOERR, "Failed to put mach number with error "~nc_ret.to!string);
	nc_ret = nc_put_var_float(file_handle, fpa_varid, &dataset.flight_path_angle);
	enforce(nc_ret == NC_NOERR, "Failed to put flight path angle with error "~nc_ret.to!string);
	nc_ret = nc_put_var_float(file_handle, tilt_varid, &dataset.tilt_angle);
	enforce(nc_ret == NC_NOERR, "Failed to put tilt angle with error "~nc_ret.to!string);

	nc_ret = nc_put_var_float(file_handle, frequency_varid, dataset.channels[0].frequency.ptr);
	enforce(nc_ret == NC_NOERR, "Failed to put frequency with error "~nc_ret.to!string);
	nc_ret = nc_put_var_float(file_handle, phi_varid, phis.ptr);
	enforce(nc_ret == NC_NOERR, "Failed to put phis with error "~nc_ret.to!string);
	nc_ret = nc_put_var_float(file_handle, theta_varid, thetas.ptr);
	enforce(nc_ret == NC_NOERR, "Failed to put thetas with error "~nc_ret.to!string);


	size_t[3] start = 0;
	size_t[3] count = [phis.length, thetas.length, dataset.channels[0].frequency.length];

	nc_ret = nc_put_vara_float(file_handle, amplitude_varid, start.ptr, count.ptr, amplitude.ptr);
	enforce(nc_ret == NC_NOERR, "Failed to put amplitudes with error "~nc_ret.to!string);
}

int main(string[] args) {

	string directory = "./";
	string output_directory = "./";
	bool hemispheres = false;

	arraySep = ",";
	auto help_information = getopt(
		args,
		std.getopt.config.bundling,
		"directory|d", "Directory of wopwop results to process", &directory,
		"output|o", "Output directory for aam spheres", &output_directory,
		"hemisphere|e", "Only convert and output bottom hemisphere", &hemispheres
	);

	if(help_information.helpWanted) {
		defaultGetoptPrinter("Program options:", help_information.options);
		import core.stdc.stdlib : exit;
		exit(-1);
	}

	if(!output_directory.exists) {
		mkdir(output_directory);
	}

	size_t[] dataset_hashes;

	auto directories = dirEntries(directory, SpanMode.shallow).array;
	
	alias Condition = Tuple!(string, "directory", float, "m", float, "fpa", float, "tilt");

	Condition[] conditions = new Condition[directories.length];

	foreach(d_idx, dir; directories) {
		if(!dir.isDir) {
			continue;
		}
		float[] flight_condition = readText(dir.buildPath("flight_condition.txt")).split(",").map!(a => a.strip.to!float).array;

		conditions[d_idx] = Condition(dir, flight_condition[0], flight_condition[1], flight_condition[2]);
	}

	// multiSort!("a.x < b.x", "a.y < b.y", SwapStrategy.unstable)(pts1);

	conditions.multiSort!("a.m < b.m", "a.fpa < b.fpa", "a.tilt < b.tilt");

	size_t d_idx = 0;
	foreach(condition; conditions) {
		/+if(!condition.directory.isDir) {
			continue;
		}+/

		if(condition.directory.buildPath("spl_spectrum.fn").exists) {

			auto dataset = read_dataset(condition.directory);

			size_t hash = dataset.hash;

			enforce(!dataset_hashes.canFind(hash), "Found duplicate database hash");

			writeln("Converting ", condition.directory, ". Flight condition hash: ", hash);
			dataset_hashes ~= hash;

			dataset.to_aam_sphere(output_directory, d_idx, hemispheres);

			d_idx++;
		}
	}

	return 0;
	
}
