import unittest

from elbow_k_Re_calculation import (
    calculate_row,
    mass_imbalance_percent,
    report_scalar,
)


class ElbowSweepCalculationsTests(unittest.TestCase):
    def test_row_uses_volume_average_properties_for_re_and_k(self):
        row = calculate_row(
            velocity=2.0,
            iterations=50,
            cumulative_iterations=100,
            inlet_total_pressure=120.0,
            outlet_total_pressure=100.0,
            mass_flows=[1.0, 0.0, -0.98],
            density=1.2,
            viscosity=1.8e-5,
            hydraulic_diameter=0.1,
        )
        self.assertAlmostEqual(row["Re"], 13333.3333333333)
        self.assertAlmostEqual(row["k"], 20.0 / 2.4)
        self.assertAlmostEqual(row["mass_imbalance_percent"], 100 * 0.02 / 1.98)
        self.assertEqual(row["total_iterations_performed"], 100)

    def test_mass_imbalance_is_zero_for_balanced_flow(self):
        self.assertEqual(mass_imbalance_percent([1.0, 0.0, -1.0]), 0.0)

    def test_report_scalar_reads_fluent_value_and_units_shape(self):
        result = [{"inlet-pressure": [101325.0, "Pa"]}]
        self.assertEqual(report_scalar(result, "inlet-pressure"), 101325.0)


if __name__ == "__main__":
    unittest.main()
