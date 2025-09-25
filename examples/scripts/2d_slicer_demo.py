#!/usr/bin/env python3
"""
Demonstration of plot_2d_slicer function for visualizing 2D EPR data.

This function allows interactive navigation through 2D EPR data
slice by slice with a slider.
"""

import epyr
import matplotlib.pyplot as plt

def demo_2d_slicer():
    """Interactive 2D viewer demonstration."""
    
    print("🎛️ Interactive 2D EPR Viewer Demonstration")
    print("="*60)
    
    # Enable interactive mode (required for widgets)
    try:
        plt.ion()  # Interactive mode
        
        # Load 2D Rabi data
        print("📂 Loading 2D data...")
        x, y, params, path = epyr.eprload('examples/data/Rabi2D_GdCaWO4_13dB_3057G.DTA')
        
        print(f"✅ Data loaded: {y.shape}")
        print(f"   X-axis ({params.get('XAXIS_NAME', 'X')}): {x[0].shape}")
        print(f"   Y-axis ({params.get('YAXIS_NAME', 'Y')}): {x[1].shape}")
        
        print("\n🚀 Launching interactive viewer...")
        print("   - Use the slider to navigate between slices")
        print("   - Overview shows current position (red line)")
        print("   - Y-scale adjusts automatically")
        
        # Horizontal slices (default)
        print("\n1. HORIZONTAL slices (navigation in Y-axis):")
        slicer_h = epyr.plot_2d_slicer(
            x, y, params, 
            title="2D Rabi Data - Horizontal Slices",
            slice_direction='horizontal'
        )
        
        input("\nPress Enter to see vertical slices...")
        
        # Vertical slices
        print("\n2. VERTICAL slices (navigation in X-axis):")
        slicer_v = epyr.plot_2d_slicer(
            x, y, params, 
            title="2D Rabi Data - Vertical Slices",
            slice_direction='vertical'
        )
        
        print("\n✅ Demonstration complete!")
        print("💡 In Jupyter, use %matplotlib widget for optimal interactivity")
        
        return slicer_h, slicer_v
        
    except Exception as e:
        print(f"❌ Error: {e}")
        print("Make sure you're in an interactive environment (Jupyter with %matplotlib widget)")
        return None, None

def demo_with_different_data():
    """Test with different types of 2D data."""
    
    print("\n📊 Testing with different 2D data files:")
    print("="*50)
    
    # List of available 2D files
    test_files = [
        'examples/data/Rabi2D_GdCaWO4_13dB_3057G.DTA',
        'examples/data/Rabi2D_GdCaWO4_6dB_3770G_2.DTA'
    ]
    
    for i, filename in enumerate(test_files, 1):
        try:
            x, y, params, _ = epyr.eprload(filename, plot_if_possible=False)
            if y.ndim == 2:
                print(f"✅ File {i}: {filename}")
                print(f"   Dimensions: {y.shape}")
                print(f"   X-axis: {params.get('XAXIS_NAME', 'Unknown')} ({params.get('XAXIS_UNIT', 'a.u.')})")
                print(f"   Y-axis: {params.get('YAXIS_NAME', 'Unknown')} ({params.get('YAXIS_UNIT', 'a.u.')})")
            else:
                print(f"⚠️ File {i}: 1D data - {y.shape}")
        except Exception as e:
            print(f"❌ File {i}: error - {e}")

if __name__ == "__main__":
    # Main demonstration
    slicer_horizontal, slicer_vertical = demo_2d_slicer()
    
    # Test with different files
    demo_with_different_data()
    
    print("\n🎯 To use in your own scripts:")
    print("""
import epyr
import matplotlib.pyplot as plt

# Enable interactive mode in Jupyter
%matplotlib widget

# Load your 2D data
x, y, params, _ = epyr.eprload("your_2d_file.DTA")

# Visualize with interactive slicer
epyr.plot_2d_slicer(x, y, params, slice_direction='horizontal')
""")