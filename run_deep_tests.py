#!/usr/bin/env python
"""
Deep Protocol Test Runner for EPyR Tools
========================================

This script provides a comprehensive test execution framework with multiple
testing protocols and detailed reporting.

Usage:
    python run_deep_tests.py [options]
    
Test Protocols:
    --smoke     : Quick smoke tests (< 1 minute)
    --standard  : Standard comprehensive tests (< 5 minutes)  
    --deep      : Deep protocol with edge cases (< 15 minutes)
    --scientific: Scientific validation tests (< 30 minutes)
    --all       : All tests (full suite)
    
Additional Options:
    --performance : Include performance benchmarks
    --coverage    : Generate coverage report  
    --parallel    : Run tests in parallel
    --verbose     : Detailed output
    --report      : Generate detailed HTML report
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import pytest


class DeepTestRunner:
    """Comprehensive test runner for EPyR Tools."""
    
    def __init__(self):
        self.project_root = Path(__file__).parent
        self.test_dir = self.project_root / "tests"
        self.results = {}
        
    def run_protocol(self, protocol: str, options: Dict[str, bool]) -> Dict[str, any]:
        """Run a specific test protocol."""
        print(f"\n{'='*60}")
        print(f"Running {protocol.upper()} Protocol Tests")
        print(f"{'='*60}")
        
        # Build pytest command
        cmd = ["python", "-m", "pytest"]
        
        # Add protocol marker
        if protocol != "all":
            cmd.extend(["-m", protocol])
            
        # Add test paths based on protocol
        if protocol == "smoke":
            # Run essential tests only
            cmd.extend([
                "tests/test_constants.py",
                "tests/test_lineshapes.py::TestLineshapesBasic",
                "tests/test_deep_protocol.py::TestEPyRCoreModules::test_constants_module[smoke]"
            ])
        elif protocol == "standard":
            # Run comprehensive module tests
            cmd.extend([
                "tests/test_deep_protocol.py",
                "tests/test_comprehensive_suite.py::TestModuleCoverage",
                "-m", "not slow"
            ])
        elif protocol == "deep":
            # Run deep protocol tests
            cmd.extend([
                "tests/test_deep_protocol.py",
                "tests/test_comprehensive_suite.py",
                "-m", "deep or comprehensive"
            ])
        elif protocol == "scientific":
            # Run scientific validation
            cmd.extend([
                "tests/test_deep_protocol.py::TestScientificValidation",
                "tests/test_comprehensive_suite.py::TestNumericalAccuracy",
                "-m", "scientific"
            ])
        else:  # all
            cmd.append("tests/")
            
        # Add options
        if options.get("verbose", False):
            cmd.extend(["-v", "-s"])
        if options.get("coverage", False):
            cmd.extend(["--cov=epyr", "--cov-report=html", "--cov-report=term"])
        if options.get("parallel", False):
            cmd.extend(["-n", "auto"])
        if options.get("performance", False):
            cmd.extend(["-m", "not slow or performance"])
        
        # Add standard options
        cmd.extend([
            "--tb=short",
            "--durations=10",
            "--maxfail=10 if protocol != 'all' else 0",
            f"--timeout=60 if protocol == 'smoke' else 300"
        ])
        
        # Execute tests
        start_time = time.time()
        
        try:
            result = subprocess.run(cmd, cwd=self.project_root, capture_output=True, text=True)
            execution_time = time.time() - start_time
            
            # Parse results
            success = result.returncode == 0
            output = result.stdout + result.stderr
            
            return {
                "protocol": protocol,
                "success": success,
                "execution_time": execution_time,
                "output": output,
                "return_code": result.returncode,
                "command": " ".join(cmd)
            }
            
        except Exception as e:
            return {
                "protocol": protocol,
                "success": False,
                "execution_time": time.time() - start_time,
                "error": str(e),
                "command": " ".join(cmd)
            }
    
    def run_specific_tests(self, test_patterns: List[str], options: Dict[str, bool]) -> Dict[str, any]:
        """Run specific test patterns."""
        print(f"\n{'='*60}")
        print(f"Running Specific Tests: {', '.join(test_patterns)}")
        print(f"{'='*60}")
        
        cmd = ["python", "-m", "pytest"] + test_patterns
        
        if options.get("verbose", False):
            cmd.extend(["-v", "-s"])
        if options.get("coverage", False):
            cmd.extend(["--cov=epyr", "--cov-report=term"])
            
        start_time = time.time()
        result = subprocess.run(cmd, cwd=self.project_root, capture_output=True, text=True)
        execution_time = time.time() - start_time
        
        return {
            "tests": test_patterns,
            "success": result.returncode == 0,
            "execution_time": execution_time,
            "output": result.stdout + result.stderr,
            "return_code": result.returncode
        }
    
    def generate_report(self, results: List[Dict[str, any]], output_file: Optional[str] = None):
        """Generate comprehensive test report."""
        if output_file is None:
            output_file = self.project_root / "test_report.html"
            
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>EPyR Tools Deep Protocol Test Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
        .protocol {{ margin: 20px 0; border: 1px solid #ddd; border-radius: 5px; }}
        .protocol-header {{ background-color: #e8f4f8; padding: 15px; }}
        .success {{ color: green; }}
        .failure {{ color: red; }}
        .output {{ background-color: #f8f8f8; padding: 15px; font-family: monospace; font-size: 12px; white-space: pre-wrap; }}
        .summary {{ background-color: #fff3cd; padding: 15px; border-radius: 5px; margin: 20px 0; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>EPyR Tools Deep Protocol Test Report</h1>
        <p>Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p>Total Protocols Tested: {len(results)}</p>
    </div>
    
    <div class="summary">
        <h2>Summary</h2>
        <ul>
"""
        
        total_time = 0
        successful_protocols = 0
        
        for result in results:
            status = "✅ PASSED" if result["success"] else "❌ FAILED"
            protocol_name = result.get("protocol", "Unknown")
            exec_time = result.get("execution_time", 0)
            total_time += exec_time
            
            if result["success"]:
                successful_protocols += 1
                
            html_content += f"<li>{protocol_name}: <span class=\"{'success' if result['success'] else 'failure'}\">{status}</span> ({exec_time:.2f}s)</li>\n"
        
        html_content += f"""
        </ul>
        <p><strong>Overall Success Rate: {successful_protocols}/{len(results)} ({100*successful_protocols/len(results):.1f}%)</strong></p>
        <p><strong>Total Execution Time: {total_time:.2f} seconds</strong></p>
    </div>
"""
        
        # Add detailed results
        for result in results:
            protocol_name = result.get("protocol", "Unknown")
            status_class = "success" if result["success"] else "failure"
            status_text = "PASSED" if result["success"] else "FAILED"
            
            html_content += f"""
    <div class="protocol">
        <div class="protocol-header">
            <h3>{protocol_name} Protocol - <span class="{status_class}">{status_text}</span></h3>
            <p>Execution Time: {result.get('execution_time', 0):.2f} seconds</p>
            <p>Return Code: {result.get('return_code', 'N/A')}</p>
        </div>
        <div class="output">{result.get('output', 'No output available')}</div>
    </div>
"""
        
        html_content += """
</body>
</html>"""
        
        with open(output_file, 'w') as f:
            f.write(html_content)
            
        print(f"\n📊 Detailed report saved to: {output_file}")
    
    def print_summary(self, results: List[Dict[str, any]]):
        """Print test execution summary."""
        print(f"\n{'='*60}")
        print("TEST EXECUTION SUMMARY")
        print(f"{'='*60}")
        
        total_time = 0
        successful_protocols = 0
        
        for result in results:
            protocol_name = result.get("protocol", "Unknown")
            success = result["success"]
            exec_time = result.get("execution_time", 0)
            total_time += exec_time
            
            status_icon = "✅" if success else "❌"
            status_text = "PASSED" if success else "FAILED"
            
            print(f"{status_icon} {protocol_name:<15} {status_text:<8} ({exec_time:6.2f}s)")
            
            if success:
                successful_protocols += 1
        
        print(f"\n{'='*60}")
        print(f"Overall Success Rate: {successful_protocols}/{len(results)} ({100*successful_protocols/len(results):.1f}%)")
        print(f"Total Execution Time: {total_time:.2f} seconds")
        print(f"{'='*60}")


def main():
    """Main test runner entry point."""
    parser = argparse.ArgumentParser(
        description="Deep Protocol Test Runner for EPyR Tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Protocol selection
    parser.add_argument("--smoke", action="store_true", help="Run smoke tests")
    parser.add_argument("--standard", action="store_true", help="Run standard tests") 
    parser.add_argument("--deep", action="store_true", help="Run deep protocol tests")
    parser.add_argument("--scientific", action="store_true", help="Run scientific validation tests")
    parser.add_argument("--all", action="store_true", help="Run all tests")
    
    # Test options
    parser.add_argument("--performance", action="store_true", help="Include performance tests")
    parser.add_argument("--coverage", action="store_true", help="Generate coverage report")
    parser.add_argument("--parallel", action="store_true", help="Run tests in parallel")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--report", action="store_true", help="Generate HTML report")
    
    # Specific tests
    parser.add_argument("--test", "-t", action="append", help="Run specific test pattern")
    parser.add_argument("--module", "-m", help="Test specific module")
    
    args = parser.parse_args()
    
    runner = DeepTestRunner()
    results = []
    
    # Determine which protocols to run
    protocols_to_run = []
    if args.smoke:
        protocols_to_run.append("smoke")
    if args.standard:
        protocols_to_run.append("standard")
    if args.deep:
        protocols_to_run.append("deep")
    if args.scientific:
        protocols_to_run.append("scientific")
    if args.all:
        protocols_to_run.append("all")
        
    # Default to smoke if nothing specified
    if not protocols_to_run and not args.test and not args.module:
        protocols_to_run.append("smoke")
    
    options = {
        "performance": args.performance,
        "coverage": args.coverage,
        "parallel": args.parallel,
        "verbose": args.verbose
    }
    
    # Run specific tests if requested
    if args.test or args.module:
        test_patterns = []
        if args.test:
            test_patterns.extend(args.test)
        if args.module:
            test_patterns.append(f"tests/test_{args.module}*.py")
            
        result = runner.run_specific_tests(test_patterns, options)
        results.append(result)
    
    # Run protocol tests
    for protocol in protocols_to_run:
        result = runner.run_protocol(protocol, options)
        results.append(result)
    
    # Generate reports
    runner.print_summary(results)
    
    if args.report:
        runner.generate_report(results)
    
    # Exit with appropriate code
    all_successful = all(result["success"] for result in results)
    sys.exit(0 if all_successful else 1)


if __name__ == "__main__":
    main()