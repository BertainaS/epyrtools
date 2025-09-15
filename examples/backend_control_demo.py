#!/usr/bin/env python3
"""
Backend Control Demo for EPyR Tools

This script demonstrates how to control matplotlib backends when using 
EPyR Tools in Jupyter notebooks.
"""

import epyr
import matplotlib.pyplot as plt
import numpy as np

def demo_backend_control():
    """Demonstrate backend control options."""
    print("🎨 EPyR Tools - Backend Control Demo")
    print("=" * 40)
    
    print(f"\n📊 Current matplotlib backend: {plt.get_backend()}")
    
    print("\n🔧 Available backend control functions:")
    print("   epyr.setup_inline_backend()     # %matplotlib inline")
    print("   epyr.setup_widget_backend()     # %matplotlib widget") 
    print("   epyr.setup_notebook_backend()   # %matplotlib notebook")
    
    print("\n📋 HOW TO USE IN JUPYTER:")
    print("""
    # Option 1: Use EPyR backend functions
    import epyr
    epyr.setup_inline_backend()        # For static plots
    # or
    epyr.setup_widget_backend()        # For interactive plots
    
    # Option 2: Use standard Jupyter magic commands  
    %matplotlib inline                 # Static plots
    %matplotlib widget                 # Interactive plots (requires ipympl)
    %matplotlib notebook               # Interactive plots (older style)
    
    # Then use baseline correction normally
    x, y, params, _ = epyr.eprload("data.dsc")
    corrected, baseline = epyr.baseline_polynomial_1d(x, y, params)
    plt.plot(x, y, label='Original')
    plt.plot(x, corrected, label='Corrected')
    plt.legend()
    plt.show()
    """)
    
    print("\n💡 RECOMMENDATIONS:")
    print("   🔸 For static plots: Use inline backend")
    print("   🔸 For interactive plots: Use widget backend")
    print("   🔸 For interactive baseline selection: Use widget or notebook backend")
    
    print("\n🆘 TROUBLESHOOTING:")
    print("   • If plots don't appear: Try %matplotlib inline")
    print("   • If interactive doesn't work: Install ipympl (pip install ipympl)")
    print("   • If selection windows stuck: Use epyr.baseline.close_selector_window()")
    
    print("\n✅ Demo complete! EPyR Tools now respects your backend choice.")

if __name__ == "__main__":
    demo_backend_control()