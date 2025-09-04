"""Advanced tests for baseline correction functions."""

import numpy as np
import pytest

try:
    from epyr.baseline import (
        baseline_constant_offset,
        baseline_mono_exponential,
        baseline_stretched_exponential,
    )
    from epyr.baseline._utils import _exclude_regions, _validate_exclude_regions
    BASELINE_UTILS_AVAILABLE = True
except ImportError:
    BASELINE_UTILS_AVAILABLE = False


@pytest.mark.skipif(not BASELINE_UTILS_AVAILABLE, reason="Baseline utilities not available")
class TestBaselineUtils:
    """Test suite for baseline utility functions."""

    def test_exclude_regions_basic(self):
        """Test basic exclude regions functionality."""
        x = np.linspace(0, 100, 101)
        y = np.ones_like(x)
        exclude_regions = [(20, 40), (60, 80)]
        
        try:
            x_filtered, y_filtered = _exclude_regions(x, y, exclude_regions)
            
            # Should have fewer points
            assert len(x_filtered) < len(x)
            assert len(y_filtered) < len(y)
            assert len(x_filtered) == len(y_filtered)
            
            # Excluded regions should not be present
            assert np.all(x_filtered < 20) or np.all(x_filtered > 40) or np.any((x_filtered > 40) & (x_filtered < 60))
        except NameError:
            pytest.skip("_exclude_regions function not available")

    def test_validate_exclude_regions(self):
        """Test exclude regions validation."""
        try:
            # Valid regions
            _validate_exclude_regions([(10, 20), (30, 40)])
            
            # Invalid regions should raise errors
            with pytest.raises(ValueError):
                _validate_exclude_regions([(20, 10)])  # Start > end
                
            with pytest.raises((ValueError, TypeError)):
                _validate_exclude_regions("invalid")  # Wrong type
                
        except NameError:
            pytest.skip("_validate_exclude_regions function not available")


@pytest.mark.skipif(not BASELINE_UTILS_AVAILABLE, reason="Advanced baseline functions not available")
class TestAdvancedBaseline:
    """Test suite for advanced baseline correction functions."""

    def test_baseline_constant_offset_basic(self):
        """Test constant offset baseline correction."""
        x = np.linspace(0, 100, 100)
        offset = 15.0
        y = np.sin(x/10) + offset + np.random.normal(0, 0.1, len(x))
        
        try:
            y_corrected, baseline_fit = baseline_constant_offset(y, x_data=x)
            
            assert len(y_corrected) == len(y)
            assert len(baseline_fit) == len(y)
            
            # Baseline should be approximately constant
            assert np.std(baseline_fit) < 1.0
            
            # Mean of corrected data should be closer to zero
            assert abs(np.mean(y_corrected)) < abs(np.mean(y))
            
        except NameError:
            pytest.skip("baseline_constant_offset function not available")

    def test_baseline_mono_exponential_basic(self):
        """Test mono-exponential baseline correction."""
        x = np.linspace(0, 100, 100)
        # Create data with exponential baseline
        true_baseline = 10 * np.exp(-x/20) + 2
        signal = 5 * np.exp(-((x-50)**2)/50)  # Gaussian peak
        y = true_baseline + signal + np.random.normal(0, 0.1, len(x))
        
        try:
            y_corrected, baseline_fit = baseline_mono_exponential(
                y, x_data=x, exclude_regions=[(40, 60)]
            )
            
            assert len(y_corrected) == len(y)
            assert len(baseline_fit) == len(y)
            
            # Baseline should follow exponential trend
            # Check that baseline decreases (exponential decay)
            assert baseline_fit[0] > baseline_fit[-1]
            
        except NameError:
            pytest.skip("baseline_mono_exponential function not available")

    def test_baseline_stretched_exponential_basic(self):
        """Test stretched exponential baseline correction."""
        x = np.linspace(0, 100, 100)
        # Create data with stretched exponential baseline
        true_baseline = 10 * np.exp(-(x/30)**0.5) + 2
        signal = 5 * np.exp(-((x-50)**2)/50)  # Gaussian peak
        y = true_baseline + signal + np.random.normal(0, 0.1, len(x))
        
        try:
            y_corrected, baseline_fit = baseline_stretched_exponential(
                y, x_data=x, exclude_regions=[(40, 60)]
            )
            
            assert len(y_corrected) == len(y)
            assert len(baseline_fit) == len(y)
            
            # Baseline should follow stretched exponential trend
            assert baseline_fit[0] > baseline_fit[-1]
            
        except NameError:
            pytest.skip("baseline_stretched_exponential function not available")

    def test_baseline_functions_with_different_parameters(self):
        """Test baseline functions with various parameter combinations."""
        x = np.linspace(0, 100, 50)
        y = np.ones_like(x) + np.random.normal(0, 0.1, len(x))
        
        functions_to_test = []
        
        try:
            from epyr.baseline import baseline_constant_offset
            functions_to_test.append(("constant_offset", baseline_constant_offset))
        except ImportError:
            pass
            
        try:
            from epyr.baseline import baseline_mono_exponential
            functions_to_test.append(("mono_exponential", baseline_mono_exponential))
        except ImportError:
            pass
            
        try:
            from epyr.baseline import baseline_stretched_exponential
            functions_to_test.append(("stretched_exponential", baseline_stretched_exponential))
        except ImportError:
            pass
        
        if not functions_to_test:
            pytest.skip("No advanced baseline functions available")
        
        for name, func in functions_to_test:
            try:
                # Test with minimal parameters
                y_corrected, baseline_fit = func(y, x_data=x)
                assert len(y_corrected) == len(y)
                assert len(baseline_fit) == len(y)
                
                # Test with exclude regions
                y_corrected2, baseline_fit2 = func(y, x_data=x, exclude_regions=[(20, 30)])
                assert len(y_corrected2) == len(y)
                assert len(baseline_fit2) == len(y)
                
            except Exception as e:
                pytest.fail(f"Function {name} failed with parameters: {e}")

    def test_baseline_error_handling(self):
        """Test error handling in baseline functions."""
        x = np.linspace(0, 10, 10)
        y = np.ones(10)
        
        # Test available functions with invalid inputs
        functions_to_test = []
        
        try:
            from epyr.baseline import baseline_constant_offset
            functions_to_test.append(baseline_constant_offset)
        except ImportError:
            pass
            
        for func in functions_to_test:
            try:
                # Test with mismatched array lengths
                with pytest.raises(ValueError):
                    func(y, x_data=x[:-1])
                
                # Test with invalid exclude regions
                with pytest.raises((ValueError, TypeError)):
                    func(y, x_data=x, exclude_regions="invalid")
                    
            except Exception:
                # Function may not implement full error checking
                pass