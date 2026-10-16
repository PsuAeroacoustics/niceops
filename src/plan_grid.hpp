#pragma once

// A uniform grid over the plane, for finding what a segment or a ray may meet
// without trying everything (buildings.hpp): each item is registered in every
// cell its box overlaps, and a query walks the cells its line crosses, in
// order (Amanatides and Woo's traversal), so that a caller looking for the
// first thing it meets can stop at the cell that holds it.  It only prunes:
// the caller still tests each item exactly, and registering by box, a little
// enlarged, can only add candidates.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace niceops {

class PlanGrid
{
public:
    PlanGrid() = default;

    // `boxes` holds each item's x min, x max, y min, y max; cells are `cell`
    // metres square.
    PlanGrid(const std::vector<double>& boxes, double cell)
        : cell_(cell)
    {
        const std::size_t n = boxes.size() / 4;
        if (n == 0)
            return;
        double x1 = -std::numeric_limits<double>::infinity(), y1 = x1;
        x0_ = y0_ = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < n; i++)
        {
            x0_ = std::min(x0_, boxes[4 * i]);
            x1 = std::max(x1, boxes[4 * i + 1]);
            y0_ = std::min(y0_, boxes[4 * i + 2]);
            y1 = std::max(y1, boxes[4 * i + 3]);
        }
        x0_ -= cell_;
        y0_ -= cell_;
        nx_ = static_cast<std::size_t>((x1 - x0_) / cell_) + 2;
        ny_ = static_cast<std::size_t>((y1 - y0_) / cell_) + 2;
        std::vector<std::vector<std::uint32_t>> cells(nx_ * ny_);
        const double slack = 1e-6;
        for (std::size_t i = 0; i < n; i++)
        {
            const std::size_t ia = column(boxes[4 * i] - slack), ib = column(boxes[4 * i + 1] + slack);
            const std::size_t ja = row(boxes[4 * i + 2] - slack), jb = row(boxes[4 * i + 3] + slack);
            for (std::size_t a = ia; a <= ib; a++)
                for (std::size_t b = ja; b <= jb; b++)
                    cells[b * nx_ + a].push_back(static_cast<std::uint32_t>(i));
        }
        start_.push_back(0);
        for (const auto& c : cells)
        {
            items_.insert(items_.end(), c.begin(), c.end());
            start_.push_back(static_cast<std::uint32_t>(items_.size()));
        }
    }

    // The cells the line (x, y) + t (dx, dy), 0 <= t <= t_max, crosses, in
    // order of t: visit(begin, end, t_exit) for each with its items, the line
    // leaving the cell at t_exit; a visit returning false ends the walk.  An
    // item spans several cells, so it can be visited more than once.
    template <typename Visit>
    void walk(double x, double y, double dx, double dy, double t_max, Visit visit) const
    {
        if (start_.empty())
            return;
        const double inf = std::numeric_limits<double>::infinity();
        const double x1 = x0_ + nx_ * cell_, y1 = y0_ + ny_ * cell_;
        // Where the line is within the grid.
        double t_in = 0.0, t_out = t_max;
        auto slab = [&](double p, double d, double lo, double hi) {
            if (d == 0.0)
                return p >= lo && p < hi;
            double ta = (lo - p) / d, tb = (hi - p) / d;
            if (ta > tb)
                std::swap(ta, tb);
            t_in = std::max(t_in, ta);
            t_out = std::min(t_out, tb);
            return true;
        };
        if (!slab(x, dx, x0_, x1) || !slab(y, dy, y0_, y1) || t_in > t_out)
            return;
        std::size_t ix = column(x + t_in * dx), iy = row(y + t_in * dy);
        const int sx = dx > 0.0 ? 1 : -1, sy = dy > 0.0 ? 1 : -1;
        double next_x = dx == 0.0 ? inf : (x0_ + (ix + (dx > 0.0 ? 1 : 0)) * cell_ - x) / dx;
        double next_y = dy == 0.0 ? inf : (y0_ + (iy + (dy > 0.0 ? 1 : 0)) * cell_ - y) / dy;
        const double step_x = dx == 0.0 ? inf : cell_ / std::abs(dx), step_y = dy == 0.0 ? inf : cell_ / std::abs(dy);
        for (;;)
        {
            const double t_exit = std::min({next_x, next_y, t_out});
            const std::size_t c = iy * nx_ + ix;
            if (!visit(items_.data() + start_[c], items_.data() + start_[c + 1], t_exit) || t_exit >= t_out)
                return;
            if (next_x <= next_y)
            {
                if ((sx < 0 && ix == 0) || (sx > 0 && ix + 1 == nx_))
                    return;
                ix += sx;
                next_x += step_x;
            }
            else
            {
                if ((sy < 0 && iy == 0) || (sy > 0 && iy + 1 == ny_))
                    return;
                iy += sy;
                next_y += step_y;
            }
        }
    }

private:
    std::size_t column(double x) const
    {
        const double c = std::floor((x - x0_) / cell_);
        return static_cast<std::size_t>(std::clamp(c, 0.0, static_cast<double>(nx_ - 1)));
    }
    std::size_t row(double y) const
    {
        const double c = std::floor((y - y0_) / cell_);
        return static_cast<std::size_t>(std::clamp(c, 0.0, static_cast<double>(ny_ - 1)));
    }

    double x0_ = 0.0, y0_ = 0.0, cell_ = 1.0;
    std::size_t nx_ = 0, ny_ = 0;
    std::vector<std::uint32_t> start_;
    std::vector<std::uint32_t> items_;
};

} // namespace niceops
