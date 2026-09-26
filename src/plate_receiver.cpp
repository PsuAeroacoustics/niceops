#include "plate_receiver.hpp"

#include <algorithm>
#include <cmath>

#include <netcdf.h>

#include "util/constants.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

std::size_t dimension(int ncid, const char* name, const std::string& path)
{
    int id = 0;
    std::size_t length = 0;
    if (nc_inq_dimid(ncid, name, &id) != NC_NOERR || nc_inq_dimlen(ncid, id, &length) != NC_NOERR)
        input_error("Plate table ", path, " has no '", name, "' dimension.");
    return length;
}

std::vector<double> variable(int ncid, const char* name, std::size_t size, const std::string& path)
{
    int id = 0;
    std::vector<double> values(size);
    if (nc_inq_varid(ncid, name, &id) != NC_NOERR || nc_get_var_double(ncid, id, values.data()) != NC_NOERR)
        input_error("Plate table ", path, " has no readable '", name, "' variable.");
    return values;
}

double attribute(int ncid, const char* name, const std::string& path)
{
    double value = 0.0;
    if (nc_get_att_double(ncid, NC_GLOBAL, name, &value) != NC_NOERR)
        input_error("Plate table ", path, " has no '", name, "' attribute.");
    return value;
}

std::string text_attribute(int ncid, const char* name, const std::string& path)
{
    std::size_t length = 0;
    if (nc_inq_attlen(ncid, NC_GLOBAL, name, &length) != NC_NOERR)
        input_error("Plate table ", path, " has no '", name, "' attribute.");
    std::string value(length, '\0');
    nc_get_att_text(ncid, NC_GLOBAL, name, value.data());
    return value;
}

bool close_to(double a, double b)
{
    return std::abs(a - b) <= 1e-6 * std::max(std::abs(a), std::abs(b));
}

} // namespace

namespace niceops {

PlateReceiver PlateReceiver::load(const std::string& path)
{
    PlateReceiver plate;
    plate.source_ = path;
    int ncid = 0;
    if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
        input_error("Cannot open plate table ", path, ".");

    plate.n_bands_ = dimension(ncid, "band", path);
    plate.n_sub_ = dimension(ncid, "sub", path);
    const std::size_t n_elevations = dimension(ncid, "elevation", path);
    plate.n_azimuths_ = dimension(ncid, "azimuth", path);
    plate.centres_ = variable(ncid, "band_centre", plate.n_bands_, path);
    plate.frequencies_ = variable(ncid, "frequency", plate.n_bands_ * plate.n_sub_, path);
    plate.elevations_ = variable(ncid, "elevation", n_elevations, path);
    const auto azimuths = variable(ncid, "azimuth", plate.n_azimuths_, path);

    if (n_elevations < 2 || !std::is_sorted(plate.elevations_.begin(), plate.elevations_.end()))
        input_error("Plate table ", path, " elevations must ascend.");
    plate.azimuth_step_ = 2.0 * pi / static_cast<double>(plate.n_azimuths_);
    for (std::size_t a = 0; a < plate.n_azimuths_; a++)
    {
        if (std::abs(azimuths[a] - 360.0 * static_cast<double>(a) / plate.n_azimuths_) > 1e-9)
            input_error("Plate table ", path, " azimuths must be uniform from 0 over 360 degrees.");
    }

    plate.ground_model_ = text_attribute(ncid, "ground_model", path);
    plate.flow_resistance_ = attribute(ncid, "flow_resistance", path);
    plate.porosity_rate_ = attribute(ncid, "porosity_rate", path);
    plate.layer_depth_ = attribute(ncid, "layer_depth_m", path);
    plate.thickness_ = attribute(ncid, "thickness_m", path);
    plate.radius_ = attribute(ncid, "radius_m", path);
    plate.mic_offset_ = std::hypot(attribute(ncid, "mic_x_m", path), attribute(ncid, "mic_y_m", path));

    // Reorder [band][sub][elevation][azimuth] per part into
    // [elevation][azimuth][band][sub][part].
    const std::size_t spectrum = plate.n_bands_ * plate.n_sub_;
    const std::size_t nodes = n_elevations * plate.n_azimuths_;
    const char* parts[4] = {"P_d_real", "P_d_imag", "P_r_real", "P_r_imag"};
    auto values = std::make_shared<std::vector<double>>(nodes * spectrum * 4, 0.0);
    for (std::size_t part = 0; part < 4; part++)
    {
        const auto data = variable(ncid, parts[part], spectrum * nodes, path);
        for (std::size_t s = 0; s < spectrum; s++)
            for (std::size_t node = 0; node < nodes; node++)
                (*values)[(node * spectrum + s) * 4 + part] = data[s * nodes + node];
    }
    nc_close(ncid);
    plate.values_ = std::move(values);

    // The fast path's quadratic in log frequency through each band's edges
    // and centre, at the sub-frequencies' positions across the band.
    plate.edges_.resize(plate.n_bands_ + 1);
    for (std::size_t band = 0; band < plate.n_bands_; band++)
        plate.edges_[band] = plate.centres_[band] * std::pow(2.0, -1.0 / 6.0);
    plate.edges_.back() = plate.centres_.back() * std::pow(2.0, 1.0 / 6.0);
    plate.lagrange_.resize(3 * plate.n_sub_);
    for (std::size_t sub = 0; sub < plate.n_sub_; sub++)
    {
        // Position in log frequency across the band, 0 at the lower edge; the
        // sub-frequencies are the band's midpoints, as panam places them.
        const double t = (static_cast<double>(sub) + 0.5) / static_cast<double>(plate.n_sub_);
        const double expected = plate.centres_[0] * std::pow(2.0, (t - 0.5) / 3.0);
        if (!close_to(plate.frequencies_[sub], expected))
            input_error("Plate table ", path, " sub-frequencies are not the band's ", plate.n_sub_,
                    " equal log steps.");
        plate.lagrange_[3 * sub] = 2.0 * (t - 0.5) * (t - 1.0);
        plate.lagrange_[3 * sub + 1] = -4.0 * t * (t - 1.0);
        plate.lagrange_[3 * sub + 2] = 2.0 * t * (t - 0.5);
    }
    return plate;
}

PlateReceiver PlateReceiver::bound(const std::vector<double>& bands, const GroundModel& ground,
        double sound_speed, std::size_t first_band, std::size_t last_band) const
{
    const std::string& path = source_;
    // The table's bands must be the database's, band for band, over the bands
    // the run reads.
    last_band = std::min(last_band, bands.size());
    if (centres_.size() != bands.size())
        input_error("Plate table ", path, " has ", centres_.size(), " bands; the database has ",
                bands.size(), ".");
    for (std::size_t band = first_band; band < last_band; band++)
    {
        if (!close_to(centres_[band], bands[band]))
            input_error("Plate table ", path, " band ", band, " is ", centres_[band],
                    " Hz; the database's is ", bands[band], " Hz.");
    }
    // P_d and P_r depend on the ground the table was solved over; applying
    // them with another Q would be neither model.
    if (parse_impedance_model(ground_model_) != ground.impedance_model
            || !close_to(flow_resistance_, ground.flow_resistance))
        input_error("Plate table ", path, " was computed over ", ground_model_,
                " ground of flow resistance ", flow_resistance_,
                " kPa s/m^2; the run's ground differs. Match --impedance_model and "
                "--flow_resistance to it, or build a table for this ground (panam "
                "axisymmetric_bem.table and write_netcdf).");
    if (ground.impedance_model == ImpedanceModel::variable_porosity
            && std::abs(porosity_rate_ - ground.porosity_rate) > 1e-9)
        input_error("Plate table ", path, " porosity rate is ", porosity_rate_,
                " 1/m; the run's is ", ground.porosity_rate, ".");
    if (ground.impedance_model == ImpedanceModel::delany_bazley_layer
            && !close_to(layer_depth_, ground.layer_depth))
        input_error("Plate table ", path, " layer depth is ", layer_depth_, " m; the run's is ",
                ground.layer_depth, ".");

    PlateReceiver plate = *this;
    plate.ground_ = ground;
    plate.sound_speed_ = sound_speed;
    plate.admittance_ = ground_impedance(ground, frequencies_, sound_speed);
    plate.edge_admittance_ = ground_impedance(ground, edges_, sound_speed);
    plate.centre_admittance_ = ground_impedance(ground, centres_, sound_speed);
    return plate;
}

void PlateReceiver::level(double source_height, double ground_distance, double azimuth,
        double* correction, std::size_t first_band, std::size_t last_band) const
{
    evaluate<false>(source_height, ground_distance, azimuth, correction, first_band, last_band);
}

void PlateReceiver::level_exact(double source_height, double ground_distance, double azimuth,
        double* correction, std::size_t first_band, std::size_t last_band) const
{
    evaluate<true>(source_height, ground_distance, azimuth, correction, first_band, last_band);
}

template <bool exact>
void PlateReceiver::evaluate(double source_height, double ground_distance, double azimuth,
        double* correction, std::size_t first_band, std::size_t last_band) const
{
    last_band = std::min(last_band, n_bands_);
    // Q at the plate's top, as panam's board_level takes it.
    const double hs = source_height;
    const double image_range = std::hypot(ground_distance, hs + thickness_);
    const double cos_grazing = image_range > 0.0 ? (hs + thickness_) / image_range : 1.0;
    const double kappa_per_hz = pi * image_range / sound_speed_;

    // Bilinear stencil.  Elevation is held at the table's ends.
    const double elevation = std::atan2(hs, ground_distance) * 180.0 / pi;
    std::size_t e0 = 0;
    double we = 0.0;
    if (elevation >= elevations_.back())
    {
        e0 = elevations_.size() - 2;
        we = 1.0;
    }
    else if (elevation > elevations_.front())
    {
        e0 = static_cast<std::size_t>(
                std::upper_bound(elevations_.begin(), elevations_.end(), elevation)
                - elevations_.begin() - 1);
        we = (elevation - elevations_[e0]) / (elevations_[e0 + 1] - elevations_[e0]);
    }
    double turns = azimuth / azimuth_step_;
    turns -= std::floor(turns / static_cast<double>(n_azimuths_)) * static_cast<double>(n_azimuths_);
    auto a0 = static_cast<std::size_t>(turns);
    const double wa = turns - static_cast<double>(a0);
    a0 %= n_azimuths_;
    const std::size_t a1 = (a0 + 1) % n_azimuths_;

    const std::size_t spectrum = n_bands_ * n_sub_;
    const double* table = values_->data();
    const double* corner[4] = {
            table + ((e0 * n_azimuths_ + a0) * spectrum) * 4,
            table + ((e0 * n_azimuths_ + a1) * spectrum) * 4,
            table + (((e0 + 1) * n_azimuths_ + a0) * spectrum) * 4,
            table + (((e0 + 1) * n_azimuths_ + a1) * spectrum) * 4};
    const double weight[4] = {(1.0 - we) * (1.0 - wa), (1.0 - we) * wa, we * (1.0 - wa), we * wa};
    constexpr double ten_over_ln10 = 4.34294481903251827651;
    const bool boundary_loss = ground_.boundary_loss_correction;

    // Q at the lower edge of the first band; each band then adds its centre
    // and upper edge, which is the next band's lower edge.
    std::array<double, 2> q_low {};
    if (!exact && first_band < last_band)
        q_low = spherical_reflection_coefficient(cos_grazing, kappa_per_hz * edges_[first_band],
                edge_admittance_[first_band], boundary_loss);

    for (std::size_t band = first_band; band < last_band; band++)
    {
        std::array<double, 2> q_centre {};
        std::array<double, 2> q_high {};
        if (!exact)
        {
            q_centre = spherical_reflection_coefficient(cos_grazing,
                    kappa_per_hz * centres_[band], centre_admittance_[band], boundary_loss);
            q_high = spherical_reflection_coefficient(cos_grazing,
                    kappa_per_hz * edges_[band + 1], edge_admittance_[band + 1], boundary_loss);
        }
        double energy = 0.0;
        for (std::size_t sub = 0; sub < n_sub_; sub++)
        {
            const std::size_t s = band * n_sub_ + sub;
            double p[4] = {0.0, 0.0, 0.0, 0.0};
            for (std::size_t c = 0; c < 4; c++)
                for (std::size_t part = 0; part < 4; part++)
                    p[part] += weight[c] * corner[c][s * 4 + part];
            std::array<double, 2> q;
            if (exact)
            {
                q = spherical_reflection_coefficient(cos_grazing, kappa_per_hz * frequencies_[s],
                        admittance_[s], boundary_loss);
            }
            else
            {
                const double* l = lagrange_.data() + 3 * sub;
                q = {l[0] * q_low[0] + l[1] * q_centre[0] + l[2] * q_high[0],
                        l[0] * q_low[1] + l[1] * q_centre[1] + l[2] * q_high[1]};
            }
            const double re = p[0] + q[0] * p[2] - q[1] * p[3];
            const double im = p[1] + q[0] * p[3] + q[1] * p[2];
            energy += re * re + im * im;
        }
        correction[band] = ten_over_ln10 * std::log(energy / static_cast<double>(n_sub_));
        q_low = q_high;
    }
}

std::vector<double> PlateReceiver::gain_bounds(std::size_t first_band, std::size_t last_band) const
{
    // An interpolated P_d + Q P_r is a convex combination of the corners'
    // P_d + Q P_r, so its size is at most the weighted sum of the corners'
    // A = |P_d| + |Q|max |P_r|; by Jensen the sub-band mean of its square is
    // then at most the largest corner's mean of A^2.  So the maximum over the
    // table's nodes bounds every geometry, held elevations included.
    last_band = std::min(last_band, n_bands_);
    // The fast path's Q is the quadratic through the edges and centre, so its
    // size is at most sum_k |L_k| |Q_k|max; the exact path's is at most its own
    // |Q|max.  Take the larger, so both are bounded.
    const auto q_exact = reflection_coefficient_bounds(ground_, admittance_);
    const auto q_edge = reflection_coefficient_bounds(ground_, edge_admittance_);
    const auto q_centre = reflection_coefficient_bounds(ground_, centre_admittance_);
    std::vector<double> q_max(q_exact.size());
    for (std::size_t band = 0; band < n_bands_; band++)
    {
        for (std::size_t sub = 0; sub < n_sub_; sub++)
        {
            const double* l = lagrange_.data() + 3 * sub;
            const double quadratic = std::abs(l[0]) * q_edge[band]
                    + std::abs(l[1]) * q_centre[band] + std::abs(l[2]) * q_edge[band + 1];
            q_max[band * n_sub_ + sub] = std::max(q_exact[band * n_sub_ + sub], quadratic);
        }
    }
    const std::size_t spectrum = n_bands_ * n_sub_;
    const std::size_t nodes = elevations_.size() * n_azimuths_;
    std::vector<double> bounds(n_bands_, 1.0);
    for (std::size_t band = first_band; band < last_band; band++)
    {
        double largest = 0.0;
        for (std::size_t node = 0; node < nodes; node++)
        {
            double energy = 0.0;
            for (std::size_t sub = 0; sub < n_sub_; sub++)
            {
                const std::size_t s = band * n_sub_ + sub;
                const double* v = values_->data() + (node * spectrum + s) * 4;
                const double a = std::hypot(v[0], v[1]) + q_max[s] * std::hypot(v[2], v[3]);
                energy += a * a;
            }
            largest = std::max(largest, energy / static_cast<double>(n_sub_));
        }
        bounds[band] = largest;
    }
    return bounds;
}

std::string PlateReceiver::description() const
{
    return source_ + " (plate radius " + std::to_string(radius_) + " m, top "
            + std::to_string(thickness_) + " m up, microphone " + std::to_string(mic_offset_)
            + " m off centre; " + std::to_string(n_sub_) + " sub-frequencies per band, "
            + std::to_string(elevations_.size()) + " elevations x " + std::to_string(n_azimuths_)
            + " azimuths)";
}

} // namespace niceops
