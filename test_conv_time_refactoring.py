#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for conv_time refactoring

This script compares the original and refactored versions to ensure:
1. Correctness (same outputs)
2. Performance improvements
3. Type preservation
"""

import datetime as dt
import numpy as np
import pandas as pd
import time
from geodezyx.conv import conv_time as conv_time_orig
from geodezyx.conv import conv_time_refactored as conv_time_new


def test_single_element():
    """Test single element input/output"""
    print("\n=== Testing Single Element I/O ===")

    test_dt = dt.datetime(2020, 1, 1, 12, 30, 45)

    # Test various conversions
    tests = [
        ('dt2mjd', test_dt),
        ('dt2year_decimal', test_dt),
        ('dt2doy', test_dt),
        ('dt2posix', test_dt),
        ('dt2fracday', test_dt),
    ]

    for func_name, input_val in tests:
        func_orig = getattr(conv_time_orig, func_name)
        func_new = getattr(conv_time_new, func_name)

        result_orig = func_orig(input_val)
        result_new = func_new(input_val)

        match = (result_orig == result_new) or (abs(result_orig - result_new) < 1e-10)
        status = "✓" if match else "✗"
        print(f"{status} {func_name}: {result_orig} vs {result_new}")


def test_list_preservation():
    """Test that list input gives list output"""
    print("\n=== Testing Type Preservation ===")

    test_dts_list = [dt.datetime(2020, 1, 1) + dt.timedelta(hours=i) for i in range(5)]
    test_dts_array = np.array(test_dts_list)
    test_dts_tuple = tuple(test_dts_list)

    # Test with list
    result_list = conv_time_new.dt2mjd(test_dts_list)
    print(f"{'✓' if isinstance(result_list, list) else '✗'} List input → List output: {type(result_list)}")

    # Test with numpy array
    result_array = conv_time_new.dt2mjd(test_dts_array)
    print(f"{'✓' if isinstance(result_array, np.ndarray) else '✗'} Array input → Array output: {type(result_array)}")

    # Test with tuple
    result_tuple = conv_time_new.dt2mjd(test_dts_tuple)
    print(f"{'✓' if isinstance(result_tuple, tuple) else '✗'} Tuple input → Tuple output: {type(result_tuple)}")


def test_datetime_types():
    """Test handling of different datetime types"""
    print("\n=== Testing Datetime Type Handling ===")

    test_date = dt.datetime(2020, 1, 1, 12, 0, 0)
    test_dt_numpy = np.datetime64('2020-01-01T12:00:00')
    test_dt_pandas = pd.Timestamp('2020-01-01 12:00:00')

    # Test conversion to MJD
    mjd_native = conv_time_new.dt2mjd(test_date)
    mjd_numpy = conv_time_new.dt2mjd(test_dt_numpy)
    mjd_pandas = conv_time_new.dt2mjd(test_dt_pandas)

    print(f"Native datetime: MJD = {mjd_native}")
    print(f"Numpy datetime64: MJD = {mjd_numpy}")
    print(f"Pandas Timestamp: MJD = {mjd_pandas}")
    print(f"{'✓' if abs(mjd_native - mjd_numpy) < 1e-6 else '✗'} Native vs Numpy match")
    print(f"{'✓' if abs(mjd_native - mjd_pandas) < 1e-6 else '✗'} Native vs Pandas match")


def benchmark_performance():
    """Benchmark performance improvements"""
    print("\n=== Performance Benchmarking ===")

    # Generate test data
    sizes = [10, 100, 1000, 10000]

    for size in sizes:
        test_dts = [dt.datetime(2020, 1, 1) + dt.timedelta(hours=i) for i in range(size)]

        # Benchmark original
        start = time.time()
        try:
            result_orig = conv_time_orig.dt2mjd(test_dts)
            time_orig = time.time() - start
        except Exception as e:
            time_orig = None
            print(f"Original failed for size {size}: {e}")

        # Benchmark refactored
        start = time.time()
        result_new = conv_time_new.dt2mjd(test_dts)
        time_new = time.time() - start

        if time_orig:
            speedup = time_orig / time_new
            print(f"Size {size:5d}: Original {time_orig*1000:6.2f}ms | New {time_new*1000:6.2f}ms | Speedup: {speedup:.1f}x")
        else:
            print(f"Size {size:5d}: New {time_new*1000:6.2f}ms")


def test_correctness():
    """Test correctness against original implementation"""
    print("\n=== Testing Correctness (against original) ===")

    test_dts = [dt.datetime(2020, 1, 1) + dt.timedelta(hours=i) for i in range(100)]

    functions_to_test = [
        'dt2mjd',
        'dt2year_decimal',
        'dt2doy',
        'dt2posix',
        'dt2fracday',
        'dt2secinday',
    ]

    for func_name in functions_to_test:
        func_orig = getattr(conv_time_orig, func_name)
        func_new = getattr(conv_time_new, func_name)

        try:
            result_orig = func_orig(test_dts)
            result_new = func_new(test_dts)

            # Check if results match
            if isinstance(result_orig, (list, tuple, np.ndarray)):
                result_orig = np.array(result_orig)
                result_new = np.array(result_new)
                max_diff = np.max(np.abs(result_orig - result_new))
                match = max_diff < 1e-9
            else:
                match = (result_orig == result_new)

            status = "✓" if match else "✗"
            print(f"{status} {func_name}: Match = {match}")
        except Exception as e:
            print(f"✗ {func_name}: Error - {e}")


def test_edge_cases():
    """Test edge cases"""
    print("\n=== Testing Edge Cases ===")

    # Empty list
    try:
        result = conv_time_new.dt2mjd([])
        print(f"✓ Empty list handled: {result}")
    except Exception as e:
        print(f"✗ Empty list failed: {e}")

    # Single element list
    single_elem_list = [dt.datetime(2020, 1, 1)]
    result = conv_time_new.dt2mjd(single_elem_list)
    print(f"{'✓' if isinstance(result, list) and len(result) == 1 else '✗'} Single element list: {type(result)}, len={len(result)}")

    # Mixed types in list (should work after normalization)
    try:
        mixed_list = [
            dt.datetime(2020, 1, 1),
            np.datetime64('2020-01-02'),
            pd.Timestamp('2020-01-03')
        ]
        result = conv_time_new.dt2mjd(mixed_list)
        print(f"✓ Mixed datetime types: {len(result)} results")
    except Exception as e:
        print(f"✗ Mixed datetime types failed: {e}")


def main():
    """Run all tests"""
    print("="*70)
    print("CONV_TIME REFACTORING TEST SUITE")
    print("="*70)

    test_single_element()
    test_list_preservation()
    test_datetime_types()
    test_correctness()
    test_edge_cases()
    benchmark_performance()

    print("\n" + "="*70)
    print("Testing complete!")
    print("="*70)


if __name__ == "__main__":
    main()

