"""Numerical checks for the optional source workflows (requires NumPy)."""
import csv
import io
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from build_hover_candidate import design, fit, high_only_correction, validate
from export_emission_samples import add_emission_time, label_samples
from score_fidelity_frames import compare, error
from regrid_emission_samples import condition_bins


class FidelityToolsTests(unittest.TestCase):
    def test_condition_bins_assign_every_sample_once_with_bounded_drift(self):
        speed = np.array([40., 41., 45., 51.])
        angle = np.array([-8., -7.1, -6.2, -5.9])
        bins = condition_bins(speed, angle, 5., 2.)
        self.assertEqual(bins.shape, (4, 2))
        self.assertTrue(np.all(np.abs(speed - bins[:, 0] * 5.) <= 2.5))
        self.assertTrue(np.all(np.abs(angle - bins[:, 1] * 2.) <= 1.))
        with self.assertRaises(ValueError):
            condition_bins(speed, angle, 0., 2.)

    def test_energy_score_and_common_coverage(self):
        self.assertAlmostEqual(error(np.array([1000., 1000.]), np.array([997., 997.])), 3.)
        fields = [(name, float) for name in ('mic', 'band', 'time', 'measured', 'snr_db',
                                           'elevation', 'range', 'above_floor_db', 'predicted')]
        baseline = np.array([(1, 500, 1, 70, 20, 3, 1000, 5, 70),
                             (1, 500, 2, 70, 20, 3, 1000, 5, 70)], dtype=fields)
        candidate = baseline.copy()
        candidate['predicted'][0] = -1000
        mask, gaps = compare(baseline, candidate)
        np.testing.assert_array_equal(mask, [False, True])
        self.assertEqual(gaps, 1)
        candidate['time'][1] += 0.5
        with self.assertRaises(ValueError):
            compare(baseline, candidate)

    def test_emission_labels_use_the_source_clock(self):
        track = dict(time=np.array([10., 20., 30.]), ground_speed_knots=np.array([40., 60., 80.]),
                     fpa_deg=np.array([-10., -5., 0.]))
        labels = label_samples([15., 25., 15.], track, 200.)
        np.testing.assert_allclose(labels['ground_speed_knots'], [50., 70., 50.])
        np.testing.assert_allclose(labels['flight_path_angle_deg'], [-7.5, -2.5, -7.5])
        np.testing.assert_allclose(labels['advance_ratio'], np.array([50., 70., 50.]) * 1852 / 3600 / 200)
        with self.assertRaises(ValueError):
            label_samples([35.], track, 200.)
        with self.assertRaises(ValueError):
            label_samples([15.], track, 0.)

    def test_source_patch_refuses_an_unknown_version(self):
        with self.assertRaises(ValueError):
            add_emission_time('def changed(): pass')

    def test_fit_separates_card_offsets_from_azimuth(self):
        rows = []
        for card, offset in [('H7', -3.), ('H8', -1.), ('H9', 1.), ('H10', 3.)]:
            for angle in np.arange(-172.5, 180., 15.):
                correction = 4. + 2. * np.cos(np.radians(angle)) + offset
                rows.append((card, 283427, 500., angle, correction, 20))
        data = np.array(rows, dtype=[('card', 'U3'), ('run', int), ('band', float),
                                    ('azimuth', float), ('correction', float), ('count', int)])
        np.testing.assert_allclose(fit(data, 1), [4., 2., 0.], atol=1e-12)

    def test_validation_retains_unpredictable_held_out_offset(self):
        rows = []
        for k, card in enumerate(['H7', 'H8', 'H9', 'H10']):
            for run in [283427 + k, 283437 + k]:
                for angle in np.arange(-172.5, 180., 15.):
                    rows.append((card, run, 500., angle, 2. * k - 3., 20))
        data = np.array(rows, dtype=[('card', 'U3'), ('run', int), ('band', float),
                                    ('azimuth', float), ('correction', float), ('count', int)])
        output = io.StringIO()
        validate(data, 0, 'high', csv.writer(output))
        scores = list(csv.reader(io.StringIO(output.getvalue())))
        headings = [row for row in scores if row[1] == 'heading']
        self.assertTrue(all(float(row[6]) > 1. for row in headings))
        self.assertTrue(all(abs(float(row[7])) < 1e-12 for row in headings))

    def test_high_ring_correction_does_not_reach_the_horizon(self):
        correction = high_only_correction(lambda p, t, f: np.full(p.shape + (f.size,), 7.))
        # phi=0, theta=elevation gives a direction below the horizon.
        phi = np.zeros((1, 5))
        theta = np.array([[0., 12.5, 20., 30., 90.]])
        result = correction(phi, theta, np.array([500.]))
        np.testing.assert_allclose(result[0, :, 0], [0., 0., 7., 0., 0.], atol=1e-12)


if __name__ == '__main__':
    unittest.main()
