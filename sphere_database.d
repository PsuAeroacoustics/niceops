module sphere_database;

import std.stdio;
import std.math;
import std.array;
import std.algorithm.sorting;
import std.algorithm.iteration;
import std.algorithm.searching;
import std.string;
import std.experimental.logger;
import netcdf;

// Declare "C" functions from interpolator
extern (C)
{
    int set_interpolator3(size_t n, const double* x1, const double* x2, const double* x3);
    int barycentric(double x1, double x2, double x3, double* weights, size_t* indicies);
    double bilinear_interp(double q11, double q12, double q21, double q22,
    double x1, double x2, double y1, double y2, double x, double y);
}

private int get_scalar(T)(int ncid, string variable_name, ref T value)
{
    int varid;
    auto status = nc_inq_varid( ncid, variable_name.toStringz, &varid);
    if (NC_NOERR != status)
    {
        log( "Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        status = nc_get_var( ncid, varid, &value);
        if (NC_NOERR != status)
        {
            fatal( "Found but could not read variable ", variable_name, "in ncid ", ncid, ".");
        }
    }
    return status;
}

private void require_scalar(T)(int ncid, string variable_name, ref T value)
{
    if (NC_NOERR != get_scalar( ncid, variable_name, value))
    {
        fatal( variable_name, " required but not found.");
    }
}

private int get_array(T)(int ncid, string variable_name, size_t countp, ref T[] values)
{
    int varid;
    auto status = nc_inq_varid( ncid, variable_name.toStringz, &varid);
    if (NC_NOERR != status)
    {
        log( "Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        size_t startp = 0;
        values.length = countp;
        status = nc_get_vara( ncid, varid, &startp, &countp, values.ptr);
        if (NC_NOERR != status)
        {
            fatal( "Found but could not read variable ", variable_name, "in ncid ", ncid, ".");
        }
    }
    return status;
}

private void require_array(T)(int ncid, string variable_name, size_t countp, ref T[] values)
{
    if (NC_NOERR != get_array( ncid, variable_name, countp, values))
    {
        fatal( variable_name, " required but not found.");
    }
}

private size_t[2] bracket(T)(T[] a, T val)
{
    // Grab the range up to target value
    auto b = a.until!((a, b) => a > b)( val);
    // Index is one less the length of this range
    auto index = b.array.length - 1;

    size_t[2] bracket;
    //Wrap around if range of a exceeded
    if (index < 0)
    {
        bracket[0] = a.length - 1;
        bracket[1] = 0;
    }
    else if (index >= (a.length - 1))
    {
        bracket[0] = a.length - 1;
        bracket[1] = 0;
    }
    else
    {
        bracket[0] = index;
        bracket[1] = index + 1;
    }
    return bracket;
}

private T clamp(T)(T x)
{
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

private double barycentric_interp(in double[][][] metric, size_t i_phi,
size_t i_theta, in double[4] weights, in size_t[4] indicies)
{
    return weights[0] * metric[i_phi][i_theta][indicies[0]] + weights[1]
    * metric[i_phi][i_theta][indicies[1]] + weights[2]
    * metric[i_phi][i_theta][indicies[2]] + weights[3] * metric[i_phi][i_theta][indicies[3]];
}

class Database
{
    double[] phis;
    double[] thetas;
    double[][][] dBA;
    double[][][] EAA;
    double sphere_radius, advance_ratio_min, advance_ratio_max, flight_path_angle_min,
    flight_path_angle_max, thrust_coefficient_min, thrust_coefficient_max;
    this(string datafile_path)
    {
        //Activate logging
        sharedLog = new FileLogger( "sphere_database.log");

        // Open netCDF file
        int file_ncid;
        log( "Opening netCDF file: " ~ datafile_path);
        if (NC_NOERR != nc_open( datafile_path.toStringz, NC_NOWRITE, &file_ncid))
        {
            fatal( "Can't open netCDF file at " ~ datafile_path);
        }
        scope (exit)
        {
            log("Closing netCDF file: " ~ datafile_path);
            nc_close(file_ncid);
        }

        // Check for "same_grid" flag to detirmine if interpolation is required
        auto same_grid = false;
        int same_grid_flag;
        get_scalar( file_ncid, "same_grid", same_grid_flag);
        same_grid = same_grid_flag ? true : false;
        log( "same_grid is set to: ", same_grid);

        // TODO: Enable interpolation when same_grid == false
        if (!same_grid)
            fatal( "same_grid == false not supported");

        // Get number of groups (spheres) in database
        int numgroups;
        if (NC_NOERR != nc_inq_grps( file_ncid, &numgroups, null))
        {
            fatal( "Can't inquire netCDF groups.");
        }
        // Get the group ids
        auto group_ids = new int[numgroups];
        if (NC_NOERR != nc_inq_grps( file_ncid, &numgroups, group_ids.ptr))
        {
            fatal("Cannot get number of netCDF groups.");
        }
        log("There are ", numgroups, " groups.");

        // Allocate storage for operating conditions
        double[] advance_ratios, flight_path_angles, thrust_coefficients;
        advance_ratios.reserve(numgroups);
        flight_path_angles.reserve(numgroups);
        thrust_coefficients.reserve(numgroups);

        // Local operating parameters
        double advance_ratio, flight_path_angle, thrust_coefficient, radius, scale;
        double[] phis_scattered, thetas_scattered, dBAs_scattered, EAAs_scattered;
        size_t nchannels0;
        // Cycle through sphere groups and load data
        foreach (i, group_id; group_ids)
        {
            require_scalar( group_id, "advance_ratio", advance_ratio);
            require_scalar( group_id, "flight_path_angle", flight_path_angle);
            require_scalar( group_id, "thrust_coefficient", thrust_coefficient);
            require_scalar( group_id, "radius", radius);
            require_scalar( group_id, "rotor_scale", scale);
            advance_ratios ~= advance_ratio;
            flight_path_angles ~= flight_path_angle;
            thrust_coefficients ~= thrust_coefficient;
            log( "Sphere ", group_id, " loaded with \t μ = ", advance_ratios[i],
            ",\t γ = ", flight_path_angles[i], ",\t CT = ",
            thrust_coefficients[i], ",\t scale = ", scale, ",\t radius = ", radius);
            if (i == 0)
            {
                sphere_radius = radius * scale;
            }
            else if (!approxEqual( radius * scale, sphere_radius))
            {
                fatal( "Database error: sphere radius*scale not constant.");
            }
            // Get number of channels
            int channel_id;
            if (NC_NOERR != nc_inq_dimid( group_id, "channels", &channel_id))
            {
                fatal( "Unable to find channel dimension.");
            }
            size_t nchannels;
            if (NC_NOERR != nc_inq_dimlen( group_id, channel_id, &nchannels))
            {
                fatal( "Unable to find channel dimension.");
            }
            if (i == 0)
            {
                // Record number of channels to check for changes
                nchannels0 = nchannels;
                log( "First sphere contains ", nchannels, " channels.");
                // Allocate arrays to store scattered data
                phis_scattered = new double[]( nchannels);
                thetas_scattered = new double[]( nchannels);
                dBAs_scattered = new double[]( nchannels);
                EAAs_scattered = new double[]( nchannels);
                // Get scattered phi and theta
                require_array( group_id, "phi", phis_scattered.length, phis_scattered);
                log( "The original phi values are:\n", phis_scattered);
                require_array( group_id, "theta", thetas_scattered.length, thetas_scattered);
                log( "The original theta values are:\n", thetas_scattered);
                //Sort and identify unique elements
                phis = phis_scattered.dup.sort().uniq().array();
                thetas = thetas_scattered.dup.sort().uniq().array();
                log( "The ", phis.length, " unique phi values are:\n", phis);
                log( "The ", thetas.length, " unique theta values are:\n", thetas);
                // Allocate noise level matricies
                dBA = new double[][][]( phis.length, thetas.length, group_ids.length);
                EAA = new double[][][]( phis.length, thetas.length, group_ids.length);
            }
            else if (nchannels != nchannels0)
            {
                fatal( "Database error: number of channels is not constant.");
            }
            // Fetch scattered acoustic data for each sphere
            require_array( group_id, "dBA", dBAs_scattered.length, dBAs_scattered);
            require_array( group_id, "EAA", EAAs_scattered.length, EAAs_scattered);
            if (same_grid)
            {
                //Search through scattered data for grided values
                foreach (m, phi_target; phis)
                {
                    foreach (n, theta_target; thetas)
                    {
                        foreach (c, phi_val; phis_scattered)
                        {
                            if (approxEqual( phi_target, phi_val)
                            && approxEqual( theta_target, thetas_scattered[c]))
                            {
                                //Load sphere data into grid
                                dBA[m][n][i] = dBAs_scattered[c];
                                EAA[m][n][i] = EAAs_scattered[c];
                            }
                        }
                    }
                }
            }
            else
            {
                // TODO: interpolation. Should be unreachable
                fatal( "Not implemented");
            }
        }
        // Find range of conditions
        advance_ratio_max = advance_ratios.maxElement;
        advance_ratio_min = advance_ratios.minElement;
        flight_path_angle_max = flight_path_angles.maxElement;
        flight_path_angle_min = flight_path_angles.minElement;
        thrust_coefficient_max = thrust_coefficients.maxElement;
        thrust_coefficient_min = thrust_coefficients.minElement;
        log( "Operating condition ranges are:\n", "μ  = [", advance_ratio_min,
        " , ", advance_ratio_max, "]\n", "γ  = [", flight_path_angle_min,
        " , ", flight_path_angle_max, "]\n", "CT = [",
        thrust_coefficient_min, " , ", thrust_coefficient_max, "]\n");
        // Normalize operating conditions
        advance_ratios.each!((ref x) => x = (x - advance_ratio_min) / (
        advance_ratio_max - advance_ratio_min));
        log( "Normalized μ   = ", advance_ratios);
        flight_path_angles.each!((ref x) => x = (x - flight_path_angle_min) / (
        flight_path_angle_max - flight_path_angle_min));
        log( "Normalized γ   = ", flight_path_angles);
        thrust_coefficients.each!((ref x) => x = (x - thrust_coefficient_min) / (
        thrust_coefficient_max - thrust_coefficient_min));
        log( "Normalized CT  = ", thrust_coefficients);

        if (set_interpolator3( advance_ratios.length, advance_ratios.ptr,
        flight_path_angles.ptr, thrust_coefficients.ptr))
        {
            fatal( "Could not set interpolator.");
        }
    }

    double eval(double range, double phi, double theta, double advance_ratio,
    double flight_path_angle, double thrust_coefficient)
    {
        // Find the sphere grid points bounding the desired phi, theta
        auto phi_range = bracket( phis, phi);
        auto theta_range = bracket( thetas, theta);
        log( "Found phi (", phi, ") in range ", phis[phi_range[0]], " -- ", phis[phi_range[1]]);
        log( "Found theta (", theta, ") in range ", thetas[theta_range[0]],
        " -- ", thetas[theta_range[1]]);

        // Normalize flight condition and limit range to [0.0,1.0]
        advance_ratio = clamp( (advance_ratio - advance_ratio_min) / (
        advance_ratio_max - advance_ratio_min));
        flight_path_angle = clamp( (flight_path_angle - flight_path_angle_min) / (
        flight_path_angle_max - flight_path_angle_min));
        thrust_coefficient = clamp( (thrust_coefficient - thrust_coefficient_min) / (
        thrust_coefficient_max - thrust_coefficient_min));
        log( "Normalized operating condition is: ", "μ  = ", advance_ratio,
        ", γ  = ", flight_path_angle, ", CT = ", thrust_coefficient);

        // Get barycentric weights and indicies
        size_t[4] indicies;
        double[4] weights;
        if (barycentric( advance_ratio, flight_path_angle, thrust_coefficient,
        weights.ptr, indicies.ptr))
        {
            fatal( "Invalid barycentric interpolation.");
        }
        else
        {
            log( "Barycentric weights: ", weights, " and indicies ", indicies);
        }
        // Barycentric interpolation between conditions at each gridpoint
        auto d11 = barycentric_interp( dBA, phi_range[0], theta_range[0], weights, indicies);
        auto d12 = barycentric_interp( dBA, phi_range[0], theta_range[1], weights, indicies);
        auto d21 = barycentric_interp( dBA, phi_range[1], theta_range[0], weights, indicies);
        auto d22 = barycentric_interp( dBA, phi_range[1], theta_range[1], weights, indicies);
        auto e11 = barycentric_interp( EAA, phi_range[0], theta_range[0], weights, indicies);
        auto e12 = barycentric_interp( EAA, phi_range[0], theta_range[1], weights, indicies);
        auto e21 = barycentric_interp( EAA, phi_range[1], theta_range[0], weights, indicies);
        auto e22 = barycentric_interp( EAA, phi_range[1], theta_range[1], weights, indicies);
        // Bilinear interolation between gridpoints
        auto p1 = phis[phi_range[0]];
        auto p2 = phis[phi_range[1]];
        auto t1 = thetas[theta_range[0]];
        auto t2 = thetas[theta_range[1]];
        auto sphere_level = bilinear_interp( d11, d12, d21, d22, p1, p2, t1, t2, phi, theta);
        auto excess_factor = bilinear_interp( e11, e12, e21, e22, p1, p2, t1, t2, phi, theta);
        log( "Sphere level (dB): ", sphere_level, "\tExcess Factor (dB/km): ", excess_factor);
        //Apply spherical spreading and excess atmospheric absorption to surface level
        auto excess_attenuation = excess_factor * (range - sphere_radius) / 1000.0;
        auto spreading_loss = 20.0 * log10( range / sphere_radius);
        log( "Spreading Loss (dB): ", spreading_loss,
        "\tExcess Attenuation (dB): ", excess_attenuation);
        assert(isFinite( sphere_level));
        assert(isFinite( excess_factor));
        assert(isFinite( excess_attenuation));
        assert(isFinite( spreading_loss));
        auto ground_level = sphere_level - spreading_loss - excess_attenuation;
        log( "Ground level: ", ground_level);
        return ground_level;
    }
}
