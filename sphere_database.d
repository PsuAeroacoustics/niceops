module sphere_database;

import std.math;
import std.numeric : findRoot;
import std.array : array, staticArray;
import std.algorithm.sorting : sort;
import std.algorithm.iteration : uniq, map, each, sum;
import std.algorithm.searching : minElement, maxElement, until;
import std.range : zip;
import std.string : toStringz;
import std.experimental.logger;
import std.traits : isFloatingPoint;
import netcdf;

// Declare "C" functions from interpolator
extern (C)
{
    int set_interpolator3(size_t n, const double* x1, const double* x2, const double* x3);
    int barycentric(double x1, double x2, double x3, double* weights, size_t* indicies);
    double bilinear_interp(double q11, double q12, double q21, double q22,
            double x1, double x2, double y1, double y2, double x, double y);
}

unittest
{
    //Test bilinear interpolation
    auto x1 = -1.0;
    auto x2 = 1.0;
    auto y1 = -1.0;
    auto y2 = 1.0;
    auto q11 = 0.0;
    auto q22 = 1.0;
    auto q12 = 0.5;
    auto q21 = 0.5;
    auto midpoint = bilinear_interp(q11, q12, q21, q22, x1, x2, y1, y2, 0.5, 0.5);
    assert(midpoint <= 1.0);
    assert(midpoint >= -1.0);

    x1 = 1.0;
    x2 = 1.0;
    y1 = -1.0;
    y2 = 1.0;
    q11 = 0.0;
    q22 = 1.0;
    q12 = 1.0;
    q21 = 0.0;
    midpoint = bilinear_interp(q11, q12, q21, q22, x1, x2, y1, y2, 0.0, 0.0);
    assert(approxEqual(midpoint, 0.5));
}

private int get_scalar(T)(int ncid, string variable_name, ref T value)
{
    int varid;
    auto status = nc_inq_varid(ncid, variable_name.toStringz, &varid);
    if (NC_NOERR != status)
    {
        log("Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        status = nc_get_var(ncid, varid, &value);
        if (NC_NOERR != status)
        {
            fatal("Found but could not read variable ", variable_name, "in ncid ", ncid, ".");
        }
    }
    return status;
}

private void require_scalar(T)(int ncid, string variable_name, ref T value)
{
    if (NC_NOERR != get_scalar(ncid, variable_name, value))
    {
        fatal(variable_name, " required but not found.");
    }
}

private int get_array(T)(int ncid, string variable_name, size_t countp, ref T[] values)
{
    int varid;
    auto status = nc_inq_varid(ncid, variable_name.toStringz, &varid);
    if (NC_NOERR != status)
    {
        info("Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        size_t startp = 0;
        values.length = countp;
        status = nc_get_vara(ncid, varid, &startp, &countp, values.ptr);
        if (NC_NOERR != status)
        {
            fatal("Found but could not read variable ", variable_name, "in ncid ", ncid, ".");
        }
    }
    return status;
}

private void require_array(T)(int ncid, string variable_name, size_t countp, ref T[] values)
{
    if (NC_NOERR != get_array(ncid, variable_name, countp, values))
    {
        fatal(variable_name, " required but not found.");
    }
}

private size_t[2] bracket(T)(T[] a, T val)
{
    // Grab the range up to target value
    auto b = a.until!((a, b) => a > b)(val);
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

unittest
{
    import std.range : iota;

    auto test_range = iota(-90.0, 180.0, 5.0).array();
    auto target = 10.0;
    auto brack = bracket(test_range, target);
    // A matching value in the range
    assert(test_range[brack[0]] <= target && test_range[brack[1]] >= target);
    // A bracketed value
    target = 2.345;
    brack = bracket(test_range, target);
    assert(test_range[brack[0]] <= target && test_range[brack[1]] >= target);
    // at lower limit
    target = 0.0;
    brack = bracket(test_range, target);
    assert(test_range[brack[0]] <= target && test_range[brack[1]] >= target);
    // wrap around upper limit
    target = 180.0;
    brack = bracket(test_range, target);
    assert(test_range[brack[0]] <= target && brack[1] == 0);
}

private T clamp(T)(T x)
{
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

unittest
{
    assert(clamp(0.5) == 0.5);
    assert(clamp(-10.) == 0.);
    assert(clamp(2.3) == 1.0);
}

private double barycentric_interp(in double[][][] metric, size_t i_phi,
        size_t i_theta, in double[4] weights, in size_t[4] indicies)
{
    return weights[0] * metric[i_phi][i_theta][indicies[0]] + weights[1]
        * metric[i_phi][i_theta][indicies[1]] + weights[2]
        * metric[i_phi][i_theta][indicies[2]] + weights[3] * metric[i_phi][i_theta][indicies[3]];
}

private T dot(T)(in T[] x, in T[] y)
{
    return zip(x, y).map!"a[0]*a[1]".sum();
}

unittest
{
    assert(dot([1, 2, 3], [-1, 2, 1]) == 6);
}

private T norm(T)(in T[] x)
{
    return sqrt(x.dot(x));
}

unittest
{
    assert(approxEqual(norm([1., 2., 3.]), sqrt(1. ^^ 2 + 2. ^^ 2 + 3. ^^ 2)));
}

private T[] normalize(T)(in T[] x)
{
    return map!(a => a / x.norm())(x).array();
}

unittest
{
    assert(approxEqual(norm([1., 2., 3., 4.].normalize()), 1.0));
}

private T[3] cross(T)(in T[] x, in T[] y)
{
    assert(x.length == 3 && y.length == 3);
    return [
        x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0]
    ].staticArray;
}

unittest
{
    auto x = cross([1.0, 2.0, 3.0], [3.0, 1.0, 2.0]);
    assert(approxEqual(x[0], 1.0));
    assert(approxEqual(x[1], 7.0));
    assert(approxEqual(x[2], -5.0));
}

private T acosc(T)(T x)
{
    if (x > 1.0)
    {
        x = 1.0;
    }
    else if (x < -1.0)
    {
        x = -1.0;
    }
    return acos(x);
}

unittest
{
    assert(acosc(0.5) == acos(0.5));
    assert(acosc(10.0) == acos(1.0));
    assert(acosc(-10.0) == acos(-1.0));
}

private double[2] lookup_angles(double[3] source, double[3] velocity,
        double[3] normal, double[3] observer)
{
    double[3] r, vn, rp;
    r[] = observer[] - source[];
    r = r.normalize();
    auto n = normal.normalize();
    // Elevation angle (positive above TPP)
    auto elv = asin(fabs(n.dot(r)));
    //Project velocity vector onto TPP to get a "nose" vector
    auto v = velocity.normalize();
    vn = v.dot(n) * n[];
    vn = v[] - vn[];
    vn = normalize(vn);
    // Sphere theta angle is relative to projection of velocity onto TPP
    auto theta = acosc(vn.dot(r));
    //Define starboard as the cross product of the "nose"" and the TPP normal
    auto pn = vn.cross(n);
    //Project the propagation vector onto the plane normal to the "nose"
    rp = r.dot(vn) * vn[];
    rp = r[] - rp[];
    rp = normalize(rp);
    //Then, we can get sphere phi by computing the angle of the propagation vector in the plane
    auto phi = acosc(fabs(n.dot(rp)));
    // Sign is still uncertain, so check to see where we are in relation to starboard
    auto phi_p = acosc(pn.dot(rp) / pn.norm());
    //If we're more than 90 degrees away, phi is defined negative
    if (phi_p > 0.5 * PI)
    {
        phi = -phi;
    }
    // Convert to degrees and return phi and theta packed in an array
    phi *= 180.0 / PI;
    theta *= 180.0 / PI;
    return [phi, theta].staticArray;
}

unittest
{
    immutable r = 123.456;
    double[3] source, observer, velocity1, velocity2, velocity3, velocity4,
        normal1, normal2, normal3, normal4, normal5;
    source = [0, 0, r];
    observer = [r, 0, 0];
    velocity1 = [1, 0, 0];
    velocity2 = [0, 1, 0];
    velocity3 = [0, -1, 0];
    velocity4 = [-1, 0, 0];
    normal1 = [0, 0, 1];
    normal2 = [sqrt(2.), 0, sqrt(2.)];
    normal3 = [-sqrt(2.), 0, sqrt(2.)];
    normal4 = [0, sqrt(2.), sqrt(2.)];
    normal5 = [0, -sqrt(2.), sqrt(2.)];

    //Check various headings only
    auto angles = lookup_angles(source, velocity1, normal1, observer);
    assert(approxEqual(angles[0], 0.0) && approxEqual(angles[1], 45.0));
    angles = lookup_angles(source, velocity2, normal1, observer);
    assert(approxEqual(angles[0], 45.0) && approxEqual(angles[1], 90.0));
    angles = lookup_angles(source, velocity3, normal1, observer);
    assert(approxEqual(angles[0], -45.0) && approxEqual(angles[1], 90.0));
    angles = lookup_angles(source, velocity4, normal1, observer);
    assert(approxEqual(angles[0], 0.0) && approxEqual(angles[1], 135.0));

    //Check various TPP tilt / heading combinations
    angles = lookup_angles(source, velocity1, normal2, observer);
    assert(approxEqual(angles[1], 0.0));
    angles = lookup_angles(source, velocity2, normal2, observer);
    assert(approxEqual(angles[0], 90.0) && approxEqual(angles[1], 90.0));
    angles = lookup_angles(source, velocity3, normal2, observer);
    assert(approxEqual(angles[0], -90.0) && approxEqual(angles[1], 90.0));
    angles = lookup_angles(source, velocity1, normal3, observer);
    assert(approxEqual(angles[0], 0.0) && approxEqual(angles[1], 90.0));
    angles = lookup_angles(source, velocity2, normal3, observer);
    assert(approxEqual(angles[0], 0.0) && approxEqual(angles[1], 90.0));
    angles = lookup_angles(source, velocity3, normal3, observer);
    assert(approxEqual(angles[0], 0.0) && approxEqual(angles[1], 90.0));
    angles = lookup_angles(source, velocity2, normal2, observer);
    assert(approxEqual(angles[0], 90.0) && approxEqual(angles[1], 90.0));
    angles = lookup_angles(source, velocity1, normal4, observer);
    assert(approxEqual(angles[0], -45.0) && approxEqual(angles[1], 45.0));
    angles = lookup_angles(source, velocity1, normal5, observer);
    assert(approxEqual(angles[0], 45.0) && approxEqual(angles[1], 45.0));
}

class Database
{
    double[] phis;
    double[] thetas;
    double[][][] dBA;
    double[][][] EAA;
    const double sphere_radius, advance_ratio_min, advance_ratio_max, flight_path_angle_min,
        flight_path_angle_max, thrust_coefficient_min, thrust_coefficient_max, max_level;
    this(string datafile_path)
    {
        // Open netCDF file
        int file_ncid;
        info("Opening netCDF file: " ~ datafile_path);
        if (NC_NOERR != nc_open(datafile_path.toStringz, NC_NOWRITE, &file_ncid))
        {
            fatal("Can't open netCDF file at " ~ datafile_path);
        }
        scope (exit)
        {
            info("Closing netCDF file: " ~ datafile_path);
            nc_close(file_ncid);
        }

        // Check for "same_grid" flag to detirmine if interpolation is required
        auto same_grid = false;
        int same_grid_flag;
        get_scalar(file_ncid, "same_grid", same_grid_flag);
        same_grid = same_grid_flag ? true : false;
        info("same_grid is set to: ", same_grid);

        // TODO: Enable interpolation when same_grid == false
        if (!same_grid)
            fatal("same_grid == false not supported");

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
        info("There are ", numgroups, " groups.");

        // Allocate storage for operating conditions
        double[] advance_ratios, flight_path_angles, thrust_coefficients;
        advance_ratios.reserve(numgroups);
        flight_path_angles.reserve(numgroups);
        thrust_coefficients.reserve(numgroups);

        // Local operating parameters
        double advance_ratio, flight_path_angle, thrust_coefficient, radius,
            scale, dimensional_radius;
        double[] phis_scattered, thetas_scattered, dBAs_scattered, EAAs_scattered;
        size_t nchannels0;
        auto dba_max = -double.infinity;
        // Cycle through sphere groups and load data
        foreach (i, group_id; group_ids)
        {
            require_scalar(group_id, "advance_ratio", advance_ratio);
            require_scalar(group_id, "flight_path_angle", flight_path_angle);
            require_scalar(group_id, "thrust_coefficient", thrust_coefficient);
            require_scalar(group_id, "radius", radius);
            require_scalar(group_id, "rotor_scale", scale);
            advance_ratios ~= advance_ratio;
            flight_path_angles ~= flight_path_angle;
            thrust_coefficients ~= thrust_coefficient;
            trace("Sphere ", group_id, " loaded with \t μ = ", advance_ratios[i],
                    ",\t γ = ", flight_path_angles[i], ",\t CT = ",
                    thrust_coefficients[i], ",\t scale = ", scale, ",\t radius = ", radius);
            if (i == 0)
            {
                dimensional_radius = radius * scale;
            }
            else if (!approxEqual(radius * scale, dimensional_radius))
            {
                fatal("Database error: sphere radius*scale not constant.");
            }
            // Get number of channels
            int channel_id;
            if (NC_NOERR != nc_inq_dimid(group_id, "channels", &channel_id))
            {
                fatal("Unable to find channel dimension.");
            }
            size_t nchannels;
            if (NC_NOERR != nc_inq_dimlen(group_id, channel_id, &nchannels))
            {
                fatal("Unable to find channel dimension.");
            }
            if (i == 0)
            {
                // Record number of channels to check for changes
                nchannels0 = nchannels;
                info("First sphere contains ", nchannels, " channels.");
                // Allocate arrays to store scattered data
                phis_scattered = new double[](nchannels);
                thetas_scattered = new double[](nchannels);
                dBAs_scattered = new double[](nchannels);
                EAAs_scattered = new double[](nchannels);
                // Get scattered phi and theta
                require_array(group_id, "phi", phis_scattered.length, phis_scattered);
                trace("The original phi values are:\n", phis_scattered);
                require_array(group_id, "theta", thetas_scattered.length, thetas_scattered);
                trace("The original theta values are:\n", thetas_scattered);
                //Sort and identify unique elements
                phis = phis_scattered.dup.sort().uniq().array();
                thetas = thetas_scattered.dup.sort().uniq().array();
                info("The ", phis.length, " unique phi values are:\n", phis);
                info("The ", thetas.length, " unique theta values are:\n", thetas);
                // Allocate noise level matricies
                dBA = new double[][][](phis.length, thetas.length, group_ids.length);
                EAA = new double[][][](phis.length, thetas.length, group_ids.length);
            }
            else if (nchannels != nchannels0)
            {
                fatal("Database error: number of channels is not constant.");
            }
            // Fetch scattered acoustic data for each sphere
            require_array(group_id, "dBA", dBAs_scattered.length, dBAs_scattered);
            require_array(group_id, "EAA", EAAs_scattered.length, EAAs_scattered);
            if (same_grid)
            {
                //Search through scattered data for grided values
                foreach (m, phi_target; phis)
                {
                    foreach (n, theta_target; thetas)
                    {
                        foreach (c, phi_val; phis_scattered)
                        {
                            if (approxEqual(phi_target, phi_val)
                                    && approxEqual(theta_target, thetas_scattered[c]))
                            {
                                //Load sphere data into grid
                                dBA[m][n][i] = dBAs_scattered[c];
                                EAA[m][n][i] = EAAs_scattered[c];
                                if (dBAs_scattered[c] > dba_max)
                                {
                                    dba_max = dBAs_scattered[c];
                                }
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
        max_level = dba_max;
        info("Max level found in database is: ", max_level);
        // Set const sphere_radius to dimensional_radius
        sphere_radius = dimensional_radius;
        // Find range of conditions
        advance_ratio_max = advance_ratios.maxElement;
        advance_ratio_min = advance_ratios.minElement;
        flight_path_angle_max = flight_path_angles.maxElement;
        flight_path_angle_min = flight_path_angles.minElement;
        thrust_coefficient_max = thrust_coefficients.maxElement;
        thrust_coefficient_min = thrust_coefficients.minElement;
        info("Operating condition ranges are:\n", "μ  = [", advance_ratio_min,
                " , ", advance_ratio_max, "]\n", "γ  = [", flight_path_angle_min,
                " , ", flight_path_angle_max, "]\n", "CT = [",
                thrust_coefficient_min, " , ", thrust_coefficient_max, "]\n");
        // Normalize operating conditions
        advance_ratios.each!((ref x) => x = (x - advance_ratio_min) / (
                advance_ratio_max - advance_ratio_min));
        trace("Normalized μ   = ", advance_ratios);
        flight_path_angles.each!((ref x) => x = (x - flight_path_angle_min) / (
                flight_path_angle_max - flight_path_angle_min));
        trace("Normalized γ   = ", flight_path_angles);
        thrust_coefficients.each!((ref x) => x = (x - thrust_coefficient_min) / (
                thrust_coefficient_max - thrust_coefficient_min));
        trace("Normalized CT  = ", thrust_coefficients);

        if (set_interpolator3(advance_ratios.length, advance_ratios.ptr,
                flight_path_angles.ptr, thrust_coefficients.ptr))
        {
            fatal("Could not set interpolator.");
        }
    }

    double eval(double range, double phi, double theta, double advance_ratio,
            double flight_path_angle, double thrust_coefficient)
    {
        // Find the sphere grid points bounding the desired phi, theta
        auto phi_range = bracket(phis, phi);
        auto theta_range = bracket(thetas, theta);
        trace("Found phi (", phi, ") in range ", phis[phi_range[0]], " -- ", phis[phi_range[1]]);
        trace("Found theta (", theta, ") in range ", thetas[theta_range[0]],
                " -- ", thetas[theta_range[1]]);

        // Normalize flight condition and limit range to [0.0,1.0]
        advance_ratio = clamp((advance_ratio - advance_ratio_min) / (
                advance_ratio_max - advance_ratio_min));
        flight_path_angle = clamp((flight_path_angle - flight_path_angle_min) / (
                flight_path_angle_max - flight_path_angle_min));
        thrust_coefficient = clamp((thrust_coefficient - thrust_coefficient_min) / (
                thrust_coefficient_max - thrust_coefficient_min));
        trace("Normalized operating condition is: ", "μ  = ", advance_ratio,
                ", γ  = ", flight_path_angle, ", CT = ", thrust_coefficient);

        // Get barycentric weights and indicies
        size_t[4] indicies;
        double[4] weights;

        if (barycentric(advance_ratio, flight_path_angle, thrust_coefficient,
                weights.ptr, indicies.ptr))
        {
            fatal("Invalid barycentric interpolation.");
        }
        else
        {
            trace("Barycentric weights: ", weights, " and indicies ", indicies);
        }
        // Barycentric interpolation between conditions at each gridpoint
        auto d11 = barycentric_interp(dBA, phi_range[0], theta_range[0], weights, indicies);
        auto d12 = barycentric_interp(dBA, phi_range[0], theta_range[1], weights, indicies);
        auto d21 = barycentric_interp(dBA, phi_range[1], theta_range[0], weights, indicies);
        auto d22 = barycentric_interp(dBA, phi_range[1], theta_range[1], weights, indicies);
        auto e11 = barycentric_interp(EAA, phi_range[0], theta_range[0], weights, indicies);
        auto e12 = barycentric_interp(EAA, phi_range[0], theta_range[1], weights, indicies);
        auto e21 = barycentric_interp(EAA, phi_range[1], theta_range[0], weights, indicies);
        auto e22 = barycentric_interp(EAA, phi_range[1], theta_range[1], weights, indicies);
        // Bilinear interolation between gridpoints
        auto p1 = phis[phi_range[0]];
        auto p2 = phis[phi_range[1]];
        auto t1 = thetas[theta_range[0]];
        auto t2 = thetas[theta_range[1]];
        auto sphere_level = bilinear_interp(d11, d12, d21, d22, p1, p2, t1, t2, phi, theta);
        auto excess_factor = bilinear_interp(e11, e12, e21, e22, p1, p2, t1, t2, phi, theta);
        trace("Sphere level (dB): ", sphere_level, "\tExcess Factor (dB/km): ", excess_factor);
        //Apply spherical spreading and excess atmospheric absorption to surface level
        auto excess_attenuation = excess_factor * (range - sphere_radius) / 1000.0;
        auto spreading_loss = 20.0 * log10(range / sphere_radius);
        trace("Spreading Loss (dB): ", spreading_loss,
                "\tExcess Attenuation (dB): ", excess_attenuation);
        assert(isFinite(sphere_level));
        assert(isFinite(excess_factor));
        assert(isFinite(excess_attenuation));
        assert(isFinite(spreading_loss));
        auto ground_level = sphere_level - spreading_loss - excess_attenuation;
        trace("Ground level: ", ground_level);
        return ground_level;
    }

    double level_at_observer(double[3] observer, double[3] source, double[3] velocity,
            double[3] normal, double advance_ratio, double flight_path_angle,
            double thrust_coefficient)
    {
        double[3] r;
        r[] = observer[] - source[];
        auto slant_range = r.norm();
        trace("R = ", slant_range, " for r = ", r);
        auto angles = lookup_angles(source, velocity, normal, observer);
        trace("Looked up [phi,theta] = ", angles);
        assert(isFinite(angles[0]) && isFinite(angles[1]));
        auto level = eval(slant_range, angles[0], angles[1], advance_ratio,
                flight_path_angle, thrust_coefficient);
        assert(isFinite(level));
        return level;
    }

    double level_along_bearing(double range, double bearing, double altitude, double[3] velocity,
            double[3] normal, double advance_ratio, double flight_path_angle,
            double thrust_coefficient)
    {
        bearing *= PI / 180.0;
        // Define ground velocity, and get heading vector from it
        double[2] ground_velocity = [velocity[0], velocity[1]];
        double[2] heading = ground_velocity.normalize();
        // Place observer on ground at heading rotated by bearing
        double[3] observer = [
            range * (cos(bearing) * heading[0] + sin(bearing) * heading[1]),
            range * (sin(-bearing) * heading[0] + cos(bearing) * heading[0]), 0.0
        ];
        // Put source at altitude above ground
        double[3] source = [0.0, 0.0, altitude];
        return level_at_observer(observer, source, velocity, normal,
                advance_ratio, flight_path_angle, thrust_coefficient);
    }

    double zero_range_level(double altitude, double[3] velocity, double[3] normal,
            double advance_ratio, double flight_path_angle, double thrust_coefficient)
    {
        double[3] observer = [0., 0., 0.];
        double[3] source = [0., 0., altitude];
        return level_at_observer(observer, source, velocity, normal,
                advance_ratio, flight_path_angle, thrust_coefficient);
    }

    double find_range(double contour_level, double bearing, double altitude, double[3] velocity, double[3] normal,
            double advance_ratio, double flight_path_angle,
            double thrust_coefficient, double threshold, bool zero_range_check)
    {
        //Check if level directly below sphere (i.e. range = 0) is above the requested contour level.
        auto zero_level = zero_range_level(altitude, velocity, normal,
                advance_ratio, flight_path_angle, thrust_coefficient);
        trace("Zero range level is: ", zero_level);
        //Initialize range to zero
        auto range = 0.0;
        //If level is greater than contour level, then search for the range associated with the contour level.
        if (!zero_range_check || zero_level > contour_level)
        {
            trace("Proceeding, zero range level above contour level.");
            // Range error
            double residual(double x)
            {
                auto level = level_along_bearing(x, bearing, altitude, velocity,
                        normal, advance_ratio, flight_path_angle, thrust_coefficient);
                return level - contour_level;
            }
            // Tolerance function
            bool tolerance(double a, double b)
            {
                return fabs(a - b) <= threshold;
            }
            // Set maximum range bound from max_level assuming spherical spreading only
            auto range_max = sphere_radius * 10. ^^ (0.05 * (max_level - contour_level));
            trace("Max search range set to: ", range_max);
            range = findRoot(&residual, 0.0, range_max, &tolerance);
            trace("Found range: ", range);
            //Verify solution if zero_range_check is false
            if (!zero_range_check)
            {
                auto check_level = level_along_bearing(range, bearing, altitude, velocity,
                        normal, advance_ratio, flight_path_angle, thrust_coefficient);
                trace("Check level: ", check_level);
                if (check_level < (contour_level - threshold))
                {
                    //Bad solution, set range back to zero
                    range = 0.0;
                    trace("Reset bad range to 0.0.");
                }
            }
        }
        return range;
    }
}
