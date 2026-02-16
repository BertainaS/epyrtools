"""
Comprehensive tests for EPyR Tools lineshape functions.

Tests cover all lineshape types with various parameters, edge cases,
and numerical accuracy validation.
"""

import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy import integrate, special

from epyr.lineshapes import (
    Lineshape,
    convspec,
    create_gaussian,
    create_lorentzian,
    create_pseudo_voigt,
    create_voigt,
    gaussian,
    lorentzian,
    lshape,
    pseudo_voigt,
    voigtian,
)


class TestGaussian:
    """Test Gaussian lineshape functions"""

    def test_gaussian_basic(self):
        """Test basic Gaussian functionality"""
        x = np.linspace(-10, 10, 1000)
        y = gaussian(x, 0, 4)

        # Check peak position (allow ±1 point tolerance for numerical differences)
        assert abs(np.argmax(y) - 500) <= 1  # Center at x=0

        # Check normalization (area should be 1)
        area = np.trapz(y, x)
        assert abs(area - 1.0) < 0.01

        # Check FWHM
        half_max = np.max(y) / 2
        indices = np.where(y >= half_max)[0]
        fwhm_measured = x[indices[-1]] - x[indices[0]]
        assert abs(fwhm_measured - 4.0) < 0.1

    def test_gaussian_derivatives(self):
        """Test Gaussian derivatives"""
        x = np.linspace(-10, 10, 1000)

        # First derivative should have zero at center
        dy = gaussian(x, 0, 4, derivative=1)
        center_idx = np.argmin(np.abs(x))
        assert abs(dy[center_idx]) < 1e-10

        # Second derivative should have maximum at center
        d2y = gaussian(x, 0, 4, derivative=2)
        assert np.argmax(np.abs(d2y)) == center_idx

    def test_gaussian_phase_rotation(self):
        """Test Gaussian phase rotation"""
        x = np.linspace(-10, 10, 1000)

        # Absorption mode
        abs_mode = gaussian(x, 0, 4, phase=0)

        # Dispersion mode
        disp_mode = gaussian(x, 0, 4, phase=np.pi / 2)

        # Should be different
        assert not np.allclose(abs_mode, disp_mode)

        # 45-degree phase should be mixture
        mixed = gaussian(x, 0, 4, phase=np.pi / 4)
        assert np.any(mixed != abs_mode) and np.any(mixed != disp_mode)

    def test_gaussian_both_components(self):
        """Test return_both parameter"""
        x = np.linspace(-10, 10, 1000)
        abs_part, disp_part = gaussian(x, 0, 4, return_both=True)

        # Check they're different
        assert not np.allclose(abs_part, disp_part)

        # Absorption should be symmetric
        assert np.allclose(abs_part, abs_part[::-1], atol=1e-10)

        # Dispersion should be antisymmetric
        assert np.allclose(disp_part, -disp_part[::-1], atol=1e-10)


class TestLorentzian:
    """Test Lorentzian lineshape functions"""

    def test_lorentzian_basic(self):
        """Test basic Lorentzian functionality"""
        x = np.linspace(-20, 20, 1000)
        y = lorentzian(x, 0, 4)

        # Check peak position (allow ±1 point tolerance for numerical differences)
        assert abs(np.argmax(y) - 500) <= 1  # Center at x=0

        # Check normalization
        area = np.trapz(y, x)
        assert abs(area - 1.0) < 0.01

        # Check FWHM
        half_max = np.max(y) / 2
        indices = np.where(y >= half_max)[0]
        fwhm_measured = x[indices[-1]] - x[indices[0]]
        assert abs(fwhm_measured - 4.0) < 0.1

    def test_lorentzian_derivatives(self):
        """Test Lorentzian derivatives"""
        x = np.linspace(-20, 20, 1000)

        # First derivative should have zero at center
        dy = lorentzian(x, 0, 4, derivative=1)
        center_idx = np.argmin(np.abs(x))
        assert abs(dy[center_idx]) < 1e-10

        # Second derivative should have minimum at center
        d2y = lorentzian(x, 0, 4, derivative=2)
        assert np.argmin(d2y) == center_idx

    def test_lorentzian_wide_range(self):
        """Test Lorentzian behavior at wide range (fat tails)"""
        x = np.linspace(-100, 100, 1000)
        y = lorentzian(x, 0, 4)

        # Lorentzian should have significant tails
        tail_ratio = y[50] / y[500]  # Far tail vs center
        assert tail_ratio > 1e-6  # Should have fat tails


class TestVoigtian:
    """Test Voigtian (true Voigt) lineshape functions"""

    def test_voigt_basic(self):
        """Test basic Voigt functionality"""
        x = np.linspace(-15, 15, 1000)
        y = voigtian(x, 0, (3, 2))  # Gaussian=3, Lorentzian=2

        # Check peak position (allow ±1 point tolerance for numerical differences)
        assert abs(np.argmax(y) - 500) <= 1  # Center at x=0

        # Check normalization
        area = np.trapz(y, x)
        assert abs(area - 1.0) < 0.05  # Slightly relaxed for convolution

    def test_voigt_limiting_cases(self):
        """Test Voigt limiting cases"""
        x = np.linspace(-15, 15, 1000)

        # Pure Gaussian limit
        voigt_gauss = voigtian(x, 0, (4, 0))
        pure_gauss = gaussian(x, 0, 4)
        assert np.allclose(voigt_gauss, pure_gauss, rtol=0.1)

        # Pure Lorentzian limit
        voigt_lorentz = voigtian(x, 0, (0, 4))
        pure_lorentz = lorentzian(x, 0, 4)
        assert np.allclose(voigt_lorentz, pure_lorentz, rtol=0.1)


class TestPseudoVoigt:
    """Test pseudo-Voigt lineshape functions"""

    def test_pseudo_voigt_basic(self):
        """Test basic pseudo-Voigt functionality"""
        x = np.linspace(-15, 15, 1000)
        y = pseudo_voigt(x, 0, 4, eta=0.5)  # 50/50 mix

        # Check peak position (allow ±1 point tolerance for numerical differences)
        assert abs(np.argmax(y) - 500) <= 1

        # Check normalization
        area = np.trapz(y, x)
        assert abs(area - 1.0) < 0.05

    def test_pseudo_voigt_limiting_cases(self):
        """Test pseudo-Voigt limiting cases"""
        x = np.linspace(-15, 15, 1000)

        # Pure Gaussian (eta=0)
        pv_gauss = pseudo_voigt(x, 0, 4, eta=0)
        pure_gauss = gaussian(x, 0, 4)
        assert np.allclose(pv_gauss, pure_gauss, rtol=0.01)

        # Pure Lorentzian (eta=1)
        pv_lorentz = pseudo_voigt(x, 0, 4, eta=1)
        pure_lorentz = lorentzian(x, 0, 4)
        assert np.allclose(pv_lorentz, pure_lorentz, rtol=0.01)


class TestLineshapeClass:
    """Test unified Lineshape class"""

    def test_lineshape_creation(self):
        """Test Lineshape class creation"""
        # Gaussian
        gauss = Lineshape("gaussian", width=5.0)
        assert gauss.shape_type == "gaussian"
        assert gauss.width == 5.0

        # Lorentzian
        lorentz = Lineshape("lorentzian", width=3.0, derivative=1)
        assert lorentz.derivative == 1

    def test_lineshape_call(self):
        """Test Lineshape __call__ method"""
        x = np.linspace(-10, 10, 100)

        gauss = Lineshape("gaussian", width=4.0)
        y = gauss(x, center=2.0)

        # Should be same as direct function call
        y_direct = gaussian(x, 2.0, 4.0)
        assert np.allclose(y, y_direct)

    def test_lineshape_convenience_methods(self):
        """Test Lineshape convenience methods"""
        x = np.linspace(-10, 10, 100)
        gauss = Lineshape("gaussian", width=4.0)

        # Test absorption/dispersion
        abs_y = gauss.absorption(x, 0)
        disp_y = gauss.dispersion(x, 0)
        assert not np.allclose(abs_y, disp_y)

        # Test derivative
        deriv_y = gauss.derivative(x, 0, order=1)
        assert np.any(deriv_y != abs_y)

    def test_lineshape_parameter_modification(self):
        """Test Lineshape parameter modification methods"""
        gauss = Lineshape("gaussian", width=4.0)

        # Test width modification
        gauss_wide = gauss.set_width(8.0)
        assert gauss_wide.width == 8.0
        assert gauss.width == 4.0  # Original unchanged

        # Test derivative modification
        gauss_deriv = gauss.set_derivative(1)
        assert gauss_deriv.derivative == 1
        assert gauss.derivative == 0  # Original unchanged

    def test_lineshape_info(self):
        """Test Lineshape info method"""
        gauss = Lineshape("gaussian", width=4.0, phase=np.pi / 4)
        info = gauss.info()

        assert info["shape_type"] == "gaussian"
        assert info["width"] == 4.0
        assert info["phase"] == np.pi / 4
        assert not info["is_absorption"]
        assert not info["is_dispersion"]


class TestFactoryFunctions:
    """Test factory functions for creating lineshapes"""

    def test_factory_functions(self):
        """Test all factory functions"""
        # Gaussian
        gauss = create_gaussian(width=4.0)
        assert isinstance(gauss, Lineshape)
        assert gauss.shape_type == "gaussian"

        # Lorentzian
        lorentz = create_lorentzian(width=5.0)
        assert lorentz.shape_type == "lorentzian"

        # Voigt
        voigt = create_voigt(3.0, 2.0)
        assert voigt.shape_type == "voigt"
        assert voigt.width == (3.0, 2.0)

        # Pseudo-Voigt
        pv = create_pseudo_voigt(width=4.0, alpha=0.7)
        assert pv.shape_type == "pseudo_voigt"
        assert pv.alpha == 0.7


class TestConvspec:
    """Test convolution functions"""


def test_convspec_basic(self):
    """Test basic convolution functionality"""
    # Create delta function
    x = np.linspace(-10, 10, 1000)
    spectrum = np.zeros_like(x)
    spectrum[500] = 1.0  # Delta at center

    # Convolve with Gaussian
    step = x[1] - x[0]
    broadened = convspec(spectrum, step, width=2.0, alpha=1.0)

    # Result should be Gaussian-like
    assert np.argmax(broadened) == 500  # Peak at center

    # Should be normalized
    area = np.trapz(broadened, x)
    assert abs(area - 1.0) < 0.1


class TestInputValidation:
    """Test input validation and error handling"""

    def test_gaussian_validation(self):
        """Test Gaussian input validation"""
        x = np.array([1, 2, 3])

        with pytest.raises(ValueError):
            gaussian(x, "invalid", 4)  # Invalid center

        with pytest.raises(ValueError):
            gaussian(x, 0, -1)  # Negative width

        with pytest.raises(ValueError):
            gaussian(x, 0, 4, derivative=-2)  # Invalid derivative

    def test_lorentzian_validation(self):
        """Test Lorentzian input validation"""
        x = np.array([1, 2, 3])

        with pytest.raises(ValueError):
            lorentzian(x, 0, 0)  # Zero width

        with pytest.raises(ValueError):
            lorentzian(x, 0, 4, phase="invalid")  # Invalid phase

    def test_lineshape_validation(self):
        """Test Lineshape class input validation"""
        with pytest.raises(ValueError):
            Lineshape("invalid_type")  # Invalid shape type

        with pytest.raises(ValueError):
            Lineshape("voigt", width=4.0)  # Voigt needs tuple width


class TestNumericalAccuracy:
    """Test numerical accuracy and edge cases"""

    def test_normalization_accuracy(self):
        """Test that all lineshapes are properly normalized"""
        x = np.linspace(-50, 50, 2000)  # Wide range, high resolution

        # Test various widths
        for width in [0.1, 1.0, 5.0, 20.0]:
            # Gaussian
            y_gauss = gaussian(x, 0, width)
            area_gauss = np.trapz(y_gauss, x)
            assert abs(area_gauss - 1.0) < 0.001, f"Gaussian width={width}"

            # Lorentzian
            y_lorentz = lorentzian(x, 0, width)
            area_lorentz = np.trapz(y_lorentz, x)
            assert abs(area_lorentz - 1.0) < 0.001, f"Lorentzian width={width}"

    def test_symmetry_properties(self):
        """Test symmetry properties of lineshapes"""
        x = np.linspace(-20, 20, 1000)

        # Absorption modes should be symmetric
        for func in [gaussian, lorentzian]:
            y = func(x, 0, 4, phase=0)
            assert np.allclose(
                y, y[::-1], atol=1e-10
            ), f"{func.__name__} absorption not symmetric"

        # Dispersion modes should be antisymmetric
        for func in [gaussian, lorentzian]:
            y = func(x, 0, 4, phase=np.pi / 2)
            assert np.allclose(
                y, -y[::-1], atol=1e-10
            ), f"{func.__name__} dispersion not antisymmetric"

    def test_derivative_consistency(self):
        """Test that analytical derivatives match numerical ones"""
        x = np.linspace(-10, 10, 1000)
        dx = x[1] - x[0]

        for func in [gaussian, lorentzian]:
            # Function and analytical first derivative
            y0 = func(x, 0, 4)
            y1_analytical = func(x, 0, 4, derivative=1)

            # Numerical first derivative
            y1_numerical = np.gradient(y0, dx)

            # Should match reasonably well
            correlation = np.corrcoef(y1_analytical, y1_numerical)[0, 1]
            assert correlation > 0.95, f"{func.__name__} derivative correlation too low"


class TestPerformance:
    """Performance and benchmark tests"""

    @pytest.mark.benchmark
    def test_gaussian_performance(self, benchmark):
        """Benchmark Gaussian performance"""
        x = np.linspace(-20, 20, 10000)

        def run_gaussian():
            return gaussian(x, 0, 4)

        result = benchmark(run_gaussian)
        assert len(result) == 10000

    @pytest.mark.benchmark
    def test_lineshape_class_performance(self, benchmark):
        """Benchmark Lineshape class performance"""
        x = np.linspace(-20, 20, 10000)
        gauss = Lineshape("gaussian", width=4.0)

        def run_lineshape():
            return gauss(x, 0)

        result = benchmark(run_lineshape)
        assert len(result) == 10000


# Integration tests
def test_integration_with_epyr():
    """Test integration with main epyr package"""
    import epyr

    # Should be able to import from epyr
    assert hasattr(epyr, "Lineshape")
    assert hasattr(epyr, "gaussian")
    assert hasattr(epyr, "lorentzian")

    # Should be able to use them
    x = np.linspace(-5, 5, 100)
    y = epyr.gaussian(x, 0, 2)
    assert len(y) == 100
    assert np.max(y) > 0


if __name__ == "__main__":
    # Run basic tests
    pytest.main([__file__, "-v"])
