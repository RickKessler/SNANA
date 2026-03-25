#!/usr/bin/env python3

"""
Test script for memory optimizations in create_covariance.py
Tests the new SymmetricMatrix class and memory-efficient operations
"""

import numpy as np
import time
import tracemalloc
import sys
import os

# Add current directory to path to import optimized classes
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from create_covariance import SymmetricMatrix, BlockProcessor
    print("Successfully imported optimized classes from create_covariance.py")
except ImportError as e:
    print(f"Error importing: {e}")
    sys.exit(1)

def test_symmetric_matrix(size=1000):
    """Test symmetric matrix storage efficiency"""
    print(f"\n=== Testing SymmetricMatrix with size {size}x{size} ===")
    
    # Start memory tracking
    tracemalloc.start()
    
    # Create test data
    np.random.seed(42)
    data = np.random.randn(size, size)
    symmetric_data = (data + data.T) / 2  # Make symmetric
    
    # Test 1: Memory usage comparison
    print("1. Memory usage comparison:")
    
    # Dense storage
    snapshot1 = tracemalloc.take_snapshot()
    dense_matrix = symmetric_data.copy()
    snapshot2 = tracemalloc.take_snapshot()
    
    # Symmetric storage  
    sym_matrix = SymmetricMatrix(size, dtype='float64')
    for i in range(size):
        for j in range(i, size):
            sym_matrix[i, j] = symmetric_data[i, j]
    snapshot3 = tracemalloc.take_snapshot()
    
    # Calculate memory usage
    dense_memory = sum(stat.size for stat in snapshot2.statistics('lineno'))
    sym_memory = sum(stat.size for stat in snapshot3.statistics('lineno'))
    
    print(f"   Dense matrix memory: {dense_memory / 1024**2:.2f} MB")
    print(f"   Symmetric matrix memory: {sym_memory / 1024**2:.2f} MB") 
    print(f"   Memory savings: {(1 - sym_memory/dense_memory)*100:.1f}%")
    
    # Test 2: Access speed
    print("2. Access speed test:")
    
    # Time dense access
    start = time.time()
    sum_dense = 0
    for i in range(min(100, size)):
        for j in range(min(100, size)):
            sum_dense += dense_matrix[i, j]
    dense_time = time.time() - start
    
    # Time symmetric access  
    start = time.time()
    sum_sym = 0
    for i in range(min(100, size)):
        for j in range(min(100, size)):
            sum_sym += sym_matrix[i, j]
    sym_time = time.time() - start
    
    print(f"   Dense access time: {dense_time:.6f} seconds")
    print(f"   Symmetric access time: {sym_time:.6f} seconds")
    print(f"   Values match: {abs(sum_dense - sum_sym) < 1e-10}")
    
    # Test 3: Conversion accuracy
    print("3. Conversion accuracy test:")
    converted = sym_matrix.to_dense()
    max_diff = np.max(np.abs(converted - symmetric_data))
    print(f"   Maximum difference after conversion: {max_diff}")
    print(f"   Conversion accurate: {max_diff < 1e-14}")
    
    # Cleanup
    sym_matrix.cleanup()
    
    return dense_memory, sym_memory

def test_block_processor(size=2000):
    """Test block-based matrix operations"""
    print(f"\n=== Testing BlockProcessor with size {size}x{size} ===")
    
    np.random.seed(42)
    
    # Create positive definite matrix for inversion test
    A = np.random.randn(size, size)
    matrix = np.dot(A, A.T) + np.eye(size) * 0.01
    
    processor = BlockProcessor(size, block_size=500)
    
    # Test block ranges
    print("1. Block range generation:")
    blocks = processor.get_block_ranges()
    print(f"   Generated {len(blocks)} blocks for {size}x{size} matrix")
    print(f"   First few blocks: {blocks[:3]}")
    
    # Test inversion (on smaller matrix to avoid long computation)
    if size <= 1000:
        print("2. Block inversion test:")
        start = time.time()
        
        try:
            inv_matrix = processor.process_block_inversion(matrix.copy(), None)
            inversion_time = time.time() - start
            
            # Verify inversion
            identity_test = np.dot(matrix, inv_matrix)
            max_off_diag = np.max(np.abs(identity_test - np.eye(size)))
            
            print(f"   Inversion time: {inversion_time:.3f} seconds")
            print(f"   Max off-diagonal error: {max_off_diag:.2e}")
            print(f"   Inversion successful: {max_off_diag < 1e-10}")
            
        except Exception as e:
            print(f"   Inversion failed: {e}")
    else:
        print("2. Skipping inversion test for large matrix")

def test_memory_threshold(sizes=[500, 1000, 2000, 5000]):
    """Test memory scaling with different matrix sizes"""
    print(f"\n=== Memory Scaling Test ===")
    
    results = []
    for size in sizes:
        print(f"Testing size {size}x{size}...")
        
        try:
            dense_mem, sym_mem = test_symmetric_matrix(size)
            results.append((size, dense_mem, sym_mem))
        except MemoryError:
            print(f"   MemoryError at size {size}")
            break
        except Exception as e:
            print(f"   Error at size {size}: {e}")
    
    print("\nMemory Usage Summary:")
    print("Size     Dense(MB)  Symmetric(MB)  Savings(%)")
    print("-" * 50)
    for size, dense, sym in results:
        savings = (1 - sym/dense) * 100
        print(f"{size:4d}     {dense/1024**2:8.1f}     {sym/1024**2:11.1f}     {savings:7.1f}")

def main():
    print("Memory Optimization Test Suite")
    print("=" * 40)
    
    # Test different sizes based on available memory
    try:
        # Small test
        test_symmetric_matrix(500)
        
        # Medium test  
        test_symmetric_matrix(2000)
        
        # Block processor test
        test_block_processor(1000)
        
        # Large test (only if memory allows)
        import psutil
        available_gb = psutil.virtual_memory().available / (1024**3)
        if available_gb > 8:  # Only test large sizes if >8GB available
            print(f"\nAvailable memory: {available_gb:.1f} GB")
            test_memory_threshold([1000, 2000, 3000, 5000])
        else:
            print(f"\nLimited memory ({available_gb:.1f} GB) - skipping large tests")
            
    except KeyboardInterrupt:
        print("\nTest interrupted by user")
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n=== Test completed ===")

if __name__ == "__main__":
    main()