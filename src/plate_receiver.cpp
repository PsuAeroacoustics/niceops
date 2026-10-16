#include "plate_receiver.hpp"

#include <algorithm>
#include <cmath>

#include <netcdf.h>

#include "util/constants.hpp"
#include "util/logger.hpp"
#include "util/netcdf.hpp"

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

// The variable, which must run along `dimensions`: its length is then their
// lengths' product.
std::vector<double> variable(int ncid, const char* name, const std::vector<std::string>& dimensions,
        const std::string& path)
{
    int id = 0;
    if (nc_inq_varid(ncid, name, &id) != NC_NOERR)
        input_error("Plate table ", path, " has no readable '", name, "' variable.");
    if (netcdf::dimension_names(ncid, id) != dimensions)
    {
        std::string expected;
        for (const auto& dimension : dimensions)
            expected += (expected.empty() ? "" : ", ") + dimension;
        input_error("Plate table ", path, ": '", name, "' must be dimensioned (", expected, ").");
    }
    auto values = netcdf::read_doubles(ncid, id);
    if (!values)
        input_error("Plate table ", path, " has no readable '", name, "' variable.");
    return std::move(*values);
}

double attribute(int ncid, const char* name, const std::string& path)
{
    double value = 0.0;
    switch (netcdf::single_value(ncid, NC_GLOBAL, name, value))
    {
    case netcdf::Single::missing:
        input_error("Plate table ", path, " has no '", name, "' attribute.");
    case netcdf::Single::not_one_number:
        input_error("Plate table ", path, ": '", name, "' must be one number.");
    case netcdf::Single::ok:
        break;
    }
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
    plate.centers_ = variable(ncid, "band_center", {"band"}, path);
    plate.frequencies_ = variable(ncid, "frequency", {"band", "sub"}, path);
    plate.elevations_ = variable(ncid, "elevation", {"elevation"}, path);
    const auto azimuths = variable(ncid, "azimuth", {"azimuth"}, path);

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
        const auto data = variable(ncid, parts[part], {"band", "sub", "elevation", "azimuth"}, path);
        for (std::size_t s = 0; s < spectrum; s++)
            for (std::size_t node = 0; node < nodes; node++)
                (*values)[(node * spectrum + s) * 4 + part] = data[s * nodes + node];
    }
    nc_close(ncid);
    plate.values_ = std::move(values);

    // The fast path's quadratic in log frequency through each band's edges
    // and center, at the sub-frequencies' positions across the band.
    plate.edges_.resize(plate.n_bands_ + 1);
    for (std::size_t band = 0; band < plate.n_bands_; band++)
        plate.edges_[band] = plate.centers_[band] * std::pow(2.0, -1.0 / 6.0);
    plate.edges_.back() = plate.centers_.back() * std::pow(2.0, 1.0 / 6.0);
    plate.lagrange_.resize(3 * plate.n_sub_);
    for (std::size_t sub = 0; sub < plate.n_sub_; sub++)
    {
        // Position in log frequency across the band, 0 at the lower edge; the
        // sub-frequencies are the band's midpoints, as panam places them.
        const double t = (static_cast<double>(sub) + 0.5) / static_cast<double>(plate.n_sub_);
        const double expected = plate.centers_[0] * std::pow(2.0, (t - 0.5) / 3.0);
        if (!close_to(plate.frequencies_[sub], expected))
            input_error("Plate table ", path, " sub-frequencies are not the band's ", plate.n_sub_,
                    " equal log steps.");
        plate.lagrange_[3 * sub] = 2.0 * (t - 0.5) * (t - 1.0);
        plate.lagrange_[3 * sub + 1] = -4.0 * t * (t - 1.0);
        plate.lagrange_[3 * sub + 2] = 2.0 * t * (t - 0.5);
    }
    return plate;
}

std::string PlateReceiver::ground_mismatch(const GroundModel& ground) const
{
    const std::string& path = source_;
    // P_d and P_r depend on the ground the table was solved over; applying
    // them with another Q would be neither model.
    if (parse_impedance_model(ground_model_) != ground.impedance_model
            || !close_to(flow_resistance_, ground.flow_resistance))
        return detail::format_message("Plate table ", path, " was computed over ", ground_model_,
                " ground of flow resistance ", flow_resistance_,
                " kPa s/m^2; the run's ground differs. Match --impedance_model and "
                "--flow_resistance to it, or build a table for this ground (panam "
                "axisymmetric_bem.table and write_netcdf).");
    if (ground.impedance_model == ImpedanceModel::variable_porosity
            && std::abs(porosity_rate_ - ground.porosity_rate) > 1e-9)
        return detail::format_message("Plate table ", path, " porosity rate is ", porosity_rate_,
                " 1/m; the run's is ", ground.porosity_rate, ".");
    if (ground.impedance_model == ImpedanceModel::delany_bazley_layer
            && !close_to(layer_depth_, ground.layer_depth))
        return detail::format_message("Plate table ", path, " layer depth is ", layer_depth_,
                " m; the run's is ", ground.layer_depth, ".");
    return {};
}

PlateReceiver PlateReceiver::bound(const std::vector<double>& bands, const GroundModel& ground,
        double sound_speed, std::size_t first_band, std::size_t last_band) const
{
    const std::string& path = source_;
    // The table's bands must be the database's, band for band, over the bands
    // the run reads.
    last_band = std::min(last_band, bands.size());
    if (centers_.size() != bands.size())
        input_error("Plate table ", path, " has ", centers_.size(), " bands; the database has ",
                bands.size(), ".");
    for (std::size_t band = first_band; band < last_band; band++)
    {
        if (!close_to(centers_[band], bands[band]))
            input_error("Plate table ", path, " band ", band, " is ", centers_[band],
                    " Hz; the database's is ", bands[band], " Hz.");
    }
    if (const std::string mismatch = ground_mismatch(ground); !mismatch.empty())
        input_error(mismatch);

    PlateReceiver plate = *this;
    plate.ground_ = ground;
    plate.sound_speed_ = sound_speed;
    plate.admittance_ = ground_impedance(ground, frequencies_, sound_speed);
    plate.edge_admittance_ = ground_impedance(ground, edges_, sound_speed);
    plate.center_admittance_ = ground_impedance(ground, centers_, sound_speed);
    return plate;
}

void PlateReceiver::level(double source_height, double ground_distance, double azimuth,
        double* correction, std::size_t first_band, std::size_t last_band) const
{
    evaluate<false>(source_height, ground_distance, azimuth, correction, first_band, last_band);
}

void PlateReceiver::power_ratio(double source_height, double ground_distance, double azimuth,
        double* ratio, std::size_t first_band, std::size_t last_band) const
{
    evaluate<false, true>(source_height, ground_distance, azimuth, ratio, first_band, last_band);
}

void PlateReceiver::level_exact(double source_height, double ground_distance, double azimuth,
        double* correction, std::size_t first_band, std::size_t last_band) const
{
    evaluate<true>(source_height, ground_distance, azimuth, correction, first_band, last_band);
}

void PlateReceiver::tabulate(double r_max, int angle_nodes, int nodes_per_decade)
{
    const double r_min = 1.0;
    r_max = std::max(r_max, 2.0);
    q_n_s_ = std::max(angle_nodes, 2);
    q_log_r_min_ = std::log(r_min);
    q_log_r_step_ = std::log(10.0) / std::max(nodes_per_decade, 1);
    q_n_r_ = std::max(2, static_cast<int>(std::ceil((std::log(r_max) - q_log_r_min_) / q_log_r_step_)) + 1);
    const std::size_t points = 2 * n_bands_ + 1;
    q_table_.assign(static_cast<std::size_t>(q_n_r_) * q_n_s_ * points * 2, 0.0F);
    const bool boundary_loss = ground_.boundary_loss_correction;
    for (int ir = 0; ir < q_n_r_; ir++)
    {
        const double r = std::exp(q_log_r_min_ + ir * q_log_r_step_);
        const double kappa_per_hz = pi * r / sound_speed_;
        for (int is = 0; is < q_n_s_; is++)
        {
            const double s = static_cast<double>(is) / (q_n_s_ - 1);
            const double cos_grazing = s * s;
            float* out = q_table_.data() + ((static_cast<std::size_t>(ir) * q_n_s_ + is) * points) * 2;
            for (std::size_t p = 0; p < points; p++)
            {
                const bool edge = p % 2 == 0;
                const std::size_t band = p / 2;
                const double f = edge ? edges_[band] : centers_[band];
                const auto beta = edge ? edge_admittance_[band] : center_admittance_[band];
                const auto q = spherical_reflection_coefficient(cos_grazing, kappa_per_hz * f, beta,
                        boundary_loss);
                out[2 * p] = static_cast<float>(q[0]);
                out[2 * p + 1] = static_cast<float>(q[1]);
            }
        }
    }
}

template <bool exact, bool as_ratio>
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

    // The fast path's Q at point p (see q_table_): from the table when there
    // is one, else evaluated.
    const bool tabulated = !exact && !q_table_.empty();
    const float* q_row[4] = {nullptr, nullptr, nullptr, nullptr};
    double q_weight[4] = {0.0, 0.0, 0.0, 0.0};
    if (tabulated)
    {
        const std::size_t points = 2 * n_bands_ + 1;
        double fs = std::sqrt(std::min(std::max(cos_grazing, 0.0), 1.0)) * (q_n_s_ - 1);
        double fr = (std::log(std::max(image_range, 1e-300)) - q_log_r_min_) / q_log_r_step_;
        fr = std::min(std::max(fr, 0.0), static_cast<double>(q_n_r_ - 1));
        const int is = std::min(static_cast<int>(fs), q_n_s_ - 2);
        const int ir = std::min(static_cast<int>(fr), q_n_r_ - 2);
        const double ws = fs - is;
        const double wr = fr - ir;
        const float* base = q_table_.data() + ((static_cast<std::size_t>(ir) * q_n_s_ + is) * points) * 2;
        q_row[0] = base;
        q_row[1] = base + points * 2;
        q_row[2] = base + static_cast<std::size_t>(q_n_s_) * points * 2;
        q_row[3] = q_row[2] + points * 2;
        q_weight[0] = (1.0 - wr) * (1.0 - ws);
        q_weight[1] = (1.0 - wr) * ws;
        q_weight[2] = wr * (1.0 - ws);
        q_weight[3] = wr * ws;
    }
    auto q_at = [&](std::size_t point, double f, std::complex<double> beta) -> std::array<double, 2> {
        if (tabulated)
        {
            double re = 0.0, im = 0.0;
            for (int c = 0; c < 4; c++)
            {
                re += q_weight[c] * q_row[c][2 * point];
                im += q_weight[c] * q_row[c][2 * point + 1];
            }
            return {re, im};
        }
        return spherical_reflection_coefficient(cos_grazing, kappa_per_hz * f, beta, boundary_loss);
    };

    // Q at the lower edge of the first band; each band then adds its center
    // and upper edge, which is the next band's lower edge.
    std::array<double, 2> q_low {};
    if (!exact && first_band < last_band)
        q_low = q_at(2 * first_band, edges_[first_band], edge_admittance_[first_band]);

    // P_d and P_r at this direction for every sub-frequency in range, as one
    // weighted sum of the four corners' contiguous runs (vectorizes), rather
    // than corner by corner inside the band loop.
    thread_local std::vector<double> interpolated;
    const std::size_t lo = first_band * n_sub_ * 4;
    const std::size_t hi = last_band * n_sub_ * 4;
    interpolated.resize(spectrum * 4);
    {
        const double w0 = weight[0], w1 = weight[1], w2 = weight[2], w3 = weight[3];
        const double* c0 = corner[0];
        const double* c1 = corner[1];
        const double* c2 = corner[2];
        const double* c3 = corner[3];
        double* out = interpolated.data();
        for (std::size_t i = lo; i < hi; i++)
            out[i] = w0 * c0[i] + w1 * c1[i] + w2 * c2[i] + w3 * c3[i];
    }

    for (std::size_t band = first_band; band < last_band; band++)
    {
        std::array<double, 2> q_center {};
        std::array<double, 2> q_high {};
        if (!exact)
        {
            q_center = q_at(2 * band + 1, centers_[band], center_admittance_[band]);
            q_high = q_at(2 * band + 2, edges_[band + 1], edge_admittance_[band + 1]);
        }
        double energy = 0.0;
        for (std::size_t sub = 0; sub < n_sub_; sub++)
        {
            const std::size_t s = band * n_sub_ + sub;
            const double* p = interpolated.data() + s * 4;
            std::array<double, 2> q;
            if (exact)
            {
                q = spherical_reflection_coefficient(cos_grazing, kappa_per_hz * frequencies_[s],
                        admittance_[s], boundary_loss);
            }
            else
            {
                const double* l = lagrange_.data() + 3 * sub;
                q = {l[0] * q_low[0] + l[1] * q_center[0] + l[2] * q_high[0],
                        l[0] * q_low[1] + l[1] * q_center[1] + l[2] * q_high[1]};
            }
            const double re = p[0] + q[0] * p[2] - q[1] * p[3];
            const double im = p[1] + q[0] * p[3] + q[1] * p[2];
            energy += re * re + im * im;
        }
        const double mean = energy / static_cast<double>(n_sub_);
        correction[band] = as_ratio ? mean : ten_over_ln10 * std::log(mean);
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
    // The fast path's Q is the quadratic through the edges and center, so its
    // size is at most sum_k |L_k| |Q_k|max; the exact path's is at most its own
    // |Q|max.  Take the larger, so both are bounded.
    const auto q_exact = reflection_coefficient_bounds(ground_, admittance_);
    const auto q_edge = reflection_coefficient_bounds(ground_, edge_admittance_);
    const auto q_center = reflection_coefficient_bounds(ground_, center_admittance_);
    std::vector<double> q_max(q_exact.size());
    for (std::size_t band = 0; band < n_bands_; band++)
    {
        for (std::size_t sub = 0; sub < n_sub_; sub++)
        {
            const double* l = lagrange_.data() + 3 * sub;
            const double quadratic = std::abs(l[0]) * q_edge[band]
                    + std::abs(l[1]) * q_center[band] + std::abs(l[2]) * q_edge[band + 1];
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
            + " m off center; " + std::to_string(n_sub_) + " sub-frequencies per band, "
            + std::to_string(elevations_.size()) + " elevations x " + std::to_string(n_azimuths_)
            + " azimuths)";
}

} // namespace niceops
