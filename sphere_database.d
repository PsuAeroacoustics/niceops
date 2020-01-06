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

private int get_scalar(T)(int ncid, string variable_name, ref T value)
{
    int varid;
    auto status = nc_inq_varid(ncid, variable_name.toStringz, &varid);
    if(NC_NOERR != status)
    {
        log("Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        status = nc_get_var(ncid, varid, &value);
        if (NC_NOERR != status)
        {
            fatal( "Found but could not read variable ", variable_name, "in ncid ", ncid, ".");
        }
    }
    return status;
}

private void require_scalar(T)(int ncid, string variable_name, ref T value)
{
    if(NC_NOERR != get_scalar(ncid, variable_name, value))
    {
        fatal(variable_name, " required but not found.");
    }
}

private int get_array(T)(int ncid, string variable_name, size_t countp, ref T[] values)
{
    int varid;
    auto status = nc_inq_varid(ncid, variable_name.toStringz, &varid);
    if(NC_NOERR != status)
    {
        log("Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        size_t startp = 0;
        values.length = countp;
        status = nc_get_vara(ncid, varid, &startp, &countp, values.ptr);
        if (NC_NOERR != status)
        {
            fatal( "Found but could not read variable ", variable_name, "in ncid ", ncid, ".");
        }
    }
    return status;
}

private void require_array(T)(int ncid, string variable_name, size_t countp, ref T[] values)
{
    if(NC_NOERR != get_array(ncid, variable_name,countp, values))
    {
        fatal(variable_name, " required but not found.");
    }
}

class Database
{
    double[] phis;
    double[] thetas;
    double[][][] dBA;
    double[][][] EAA;
    double advance_ratio_min, advance_ratio_max, flight_path_angle_min, flight_path_angle_max, thrust_coefficient_min, thrust_coefficient_max;
    this(string datafile_path)
    {
        //Activate logging
        sharedLog = new FileLogger("sphere_database.log");

        // Open netCDF file
        int file_ncid;
        log("Opening netCDF file: " ~ datafile_path);
        if (NC_NOERR != nc_open(datafile_path.toStringz, NC_NOWRITE, &file_ncid))
        {
			fatal("Can't open netCDF file at " ~ datafile_path);
        }
        scope (exit)
        {
            log("Closing netCDF file: " ~ datafile_path);
            nc_close(file_ncid);
        }

        // Check for "same_grid" flag to detirmine if interpolation is required
        auto same_grid = false;
        int same_grid_flag;
        get_scalar(file_ncid,"same_grid",same_grid_flag);
        same_grid = same_grid_flag ? true : false;
        log("same_grid is set to: ",  same_grid);

        // TODO: Enable interpolation when same_grid == false
        if(!same_grid) fatal("same_grid == false not supported");

        // Get number of groups (spheres) in database
        int numgroups;
        if (NC_NOERR != nc_inq_grps(file_ncid, &numgroups, null))
        {
            fatal("Can't inquire netCDF groups.");
        }
        // Get the group ids
        auto group_ids = new int[numgroups];
        if (NC_NOERR != nc_inq_grps(file_ncid, &numgroups, group_ids.ptr))
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
        double advance_ratio, flight_path_angle, thrust_coefficient, radius, scale, sphere_radius;
        double[] phis_scattered, thetas_scattered, dBAs_scattered, EAAs_scattered;
        size_t nchannels0;
        // Cycle through sphere groups and load data
        foreach(i, group_id; group_ids)
        {
            require_scalar(group_id,"advance_ratio", advance_ratio);
            require_scalar(group_id,"flight_path_angle", flight_path_angle);
            require_scalar(group_id,"thrust_coefficient", thrust_coefficient);
            require_scalar(group_id,"radius", radius);
            require_scalar(group_id,"rotor_scale", scale);
            advance_ratios ~= advance_ratio;
            flight_path_angles ~= flight_path_angle;
            thrust_coefficients ~= thrust_coefficient;
            log("Sphere ", group_id, " loaded with \t μ = ", advance_ratios[i], ",\t γ = ", flight_path_angles[i],
            ",\t CT = ", thrust_coefficients[i], ",\t scale = ", scale, ",\t radius = ", radius);
            if(i == 0)
            {
                sphere_radius = radius*scale;
            }
            else if(!approxEqual(radius*scale,sphere_radius))
            {
                fatal("Database error: sphere radius*scale not constant.");
            }
            // Get number of channels
            int channel_id;
            if(NC_NOERR != nc_inq_dimid(group_id,"channels",&channel_id))
            {
                fatal("Unable to find channel dimension.");
            }
            size_t nchannels;
            if(NC_NOERR != nc_inq_dimlen(group_id,channel_id,&nchannels))
            {
                fatal("Unable to find channel dimension.");
            }
            if(i == 0)
            {
                // Record number of channels to check for changes
                nchannels0 = nchannels;
                log("First sphere contains ", nchannels, " channels.");
                // Allocate arrays to store scattered data
                phis_scattered = new double[](nchannels);
                thetas_scattered = new double[](nchannels);
                dBAs_scattered = new double[](nchannels);
                EAAs_scattered = new double[](nchannels);
                // Get scattered phi and theta
                require_array(group_id,"phi",phis_scattered.length,phis_scattered);
                log("The original phi values are:\n", phis_scattered);
                require_array(group_id,"theta",thetas_scattered.length,thetas_scattered);
                log("The original theta values are:\n", thetas_scattered);
                //Sort and identify unique elements
                phis = phis_scattered.dup.sort().uniq().array();
                thetas = thetas_scattered.dup.sort().uniq().array();
                log("The ", phis.length," unique phi values are:\n",phis);
                log("The ", thetas.length," unique theta values are:\n",thetas);
                // Allocate noise level matricies
                dBA = new double[][][](phis.length,thetas.length,group_ids.length);
                EAA = new double[][][](phis.length,thetas.length,group_ids.length);
            }
            else if(nchannels != nchannels0)
            {
                fatal("Database error: number of channels is not constant.");
            }
            // Fetch scattered acoustic data for each sphere
            require_array(group_id,"dBA",dBAs_scattered.length,dBAs_scattered);
            require_array(group_id,"EAA",EAAs_scattered.length,EAAs_scattered);
            if(same_grid)
            {
                //Search through scattered data for grided values
                foreach(m, phi_target; phis)
                {
                    foreach(n, theta_target; thetas)
                    {
                        foreach(c, phi_val; phis_scattered)
                        {
                            if(approxEqual(phi_target,phi_val) && approxEqual(theta_target,thetas_scattered[c]))
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
                fatal("Not implemented");
            }
        }
        // Find range of conditions
        advance_ratio_max = advance_ratios.maxElement;
        advance_ratio_min = advance_ratios.minElement;
        flight_path_angle_max = flight_path_angles.maxElement;
        flight_path_angle_min = flight_path_angles.minElement;
        thrust_coefficient_max = thrust_coefficients.maxElement;
        thrust_coefficient_min = thrust_coefficients.minElement;
        log("Operating condition ranges are:\n",
        "μ  = [", advance_ratio_min, " , " , advance_ratio_max, "]\n",
        "γ  = [", flight_path_angle_min, " , " , flight_path_angle_max, "]\n",
        "CT = [", thrust_coefficient_min, " , " , thrust_coefficient_max, "]\n");
        // Normalize operating conditions
        advance_ratios.each!((ref x) => x = (x - advance_ratio_min) / (advance_ratio_max-advance_ratio_min));
        log("Normalized μ   = ", advance_ratios);
        flight_path_angles.each!((ref x) => x = (x - flight_path_angle_min) / (flight_path_angle_max-flight_path_angle_min));
        log("Normalized γ   = ", flight_path_angles);
        thrust_coefficients.each!((ref x) => x = (x - thrust_coefficient_min) / (thrust_coefficient_max-thrust_coefficient_min));
        log("Normalized CT  = ", thrust_coefficients);
    }
}
