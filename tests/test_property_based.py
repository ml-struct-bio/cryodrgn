"""Property-based tests using Hypothesis for cryoDRGN.

This module uses the Hypothesis library to generate comprehensive test cases
and verify mathematical properties and invariants across the codebase.
"""

import pytest
import numpy as np
import torch
from hypothesis import given, strategies as st, settings, assume, note
from hypothesis.extra.numpy import arrays, array_shapes
import tempfile
from pathlib import Path

from cryodrgn import fft, utils, ctf
from cryodrgn.lattice import Lattice
from cryodrgn.source import ImageSource, write_mrc
from cryodrgn.dataset import ImageDataset


# Custom strategies for cryoDRGN-specific types
@st.composite
def valid_box_sizes(draw):
    """Generate valid box sizes (powers of 2, reasonable range)."""
    power = draw(st.integers(min_value=4, max_value=9))  # 16 to 512
    return 2 ** power

@st.composite
def ctf_parameters(draw):
    """Generate valid CTF parameters."""
    return {
        'apix': draw(st.floats(min_value=0.5, max_value=5.0)),
        'defocus_u': draw(st.floats(min_value=0.5, max_value=10.0)),
        'defocus_v': draw(st.floats(min_value=0.5, max_value=10.0)),
        'defocus_angle': draw(st.floats(min_value=0.0, max_value=180.0)),
        'voltage': draw(st.floats(min_value=80.0, max_value=500.0)),
        'cs': draw(st.floats(min_value=0.0, max_value=5.0)),
        'w': draw(st.floats(min_value=0.0, max_value=0.3)),
        'phase_shift': draw(st.floats(min_value=0.0, max_value=180.0)),
    }

@st.composite
def rotation_matrices(draw):
    """Generate valid 3D rotation matrices."""
    from scipy.spatial.transform import Rotation
    # Generate random rotation
    angles = draw(arrays(np.float64, (3,), elements=st.floats(-np.pi, np.pi)))
    return Rotation.from_euler('xyz', angles).as_matrix()

@st.composite 
def pose_data(draw, n_poses=None):
    """Generate valid pose data (rotations and translations)."""
    if n_poses is None:
        n_poses = draw(st.integers(min_value=1, max_value=100))
    
    rotations = draw(arrays(
        np.float64, 
        (n_poses, 3, 3),
        elements=st.floats(-1.0, 1.0)
    ))
    
    # Ensure they are actually rotation matrices
    for i in range(n_poses):
        U, _, Vt = np.linalg.svd(rotations[i])
        rotations[i] = U @ Vt
        if np.linalg.det(rotations[i]) < 0:
            rotations[i][:, 0] *= -1
    
    translations = draw(arrays(
        np.float64,
        (n_poses, 2), 
        elements=st.floats(-50.0, 50.0)
    ))
    
    return rotations, translations


class TestFFTProperties:
    """Property-based tests for FFT operations."""
    
    @given(
        box_size=valid_box_sizes(),
        n_images=st.integers(min_value=1, max_value=10)
    )
    @settings(max_examples=50, deadline=5000)
    def test_fft_ifft_inverse_property(self, box_size, n_images):
        """Test that FFT and IFFT are inverse operations."""
        # Generate random image data
        data = torch.randn(n_images, box_size, box_size, dtype=torch.float32)
        
        # Apply FFT then IFFT
        fft_data = fft.fft2_center(data)
        reconstructed = fft.ifft2_center(fft_data).real
        
        # Should recover original data (within numerical precision)
        assert torch.allclose(data, reconstructed, rtol=1e-5, atol=1e-6), \
            f"FFT-IFFT not inverse for size {box_size}x{box_size}"
    
    @given(
        box_size=valid_box_sizes(),
        n_images=st.integers(min_value=1, max_value=5)
    )
    @settings(max_examples=30, deadline=5000)
    def test_fft_linearity_property(self, box_size, n_images):
        """Test FFT linearity: FFT(a*x + b*y) = a*FFT(x) + b*FFT(y)."""
        x = torch.randn(n_images, box_size, box_size, dtype=torch.float32)
        y = torch.randn(n_images, box_size, box_size, dtype=torch.float32)
        a = torch.randn(1).item()
        b = torch.randn(1).item()
        
        # Linear combination
        linear_combo = a * x + b * y
        fft_combo = fft.fft2_center(linear_combo)
        
        # FFT of individual terms
        fft_x = fft.fft2_center(x)
        fft_y = fft.fft2_center(y)
        expected_combo = a * fft_x + b * fft_y
        
        assert torch.allclose(fft_combo, expected_combo, rtol=1e-5, atol=1e-6), \
            f"FFT linearity violated for size {box_size}x{box_size}"
    
    @given(box_size=valid_box_sizes())
    @settings(max_examples=20, deadline=3000)
    def test_fft_energy_conservation(self, box_size):
        """Test Parseval's theorem: energy is conserved by FFT."""
        data = torch.randn(1, box_size, box_size, dtype=torch.float32)
        fft_data = fft.fft2_center(data)
        
        # Calculate energy in both domains
        time_energy = torch.sum(data.abs() ** 2)
        freq_energy = torch.sum(fft_data.abs() ** 2) / (box_size * box_size)
        
        assert torch.allclose(time_energy, freq_energy, rtol=1e-4, atol=1e-5), \
            f"Energy not conserved in FFT for size {box_size}x{box_size}"


class TestCTFProperties:
    """Property-based tests for CTF computations."""
    
    @given(params=ctf_parameters())
    @settings(max_examples=50, deadline=3000)
    def test_ctf_symmetry_property(self, params):
        """Test that CTF has expected symmetry properties."""
        device = torch.device('cpu')  # Use CPU for deterministic behavior
        box_size = 64
        
        lattice = Lattice(box_size, device=device)
        freqs = lattice.freqs2d / params['apix']
        
        # Compute CTF
        c = ctf.compute_ctf(
            freqs,
            torch.tensor(params['defocus_u']),
            torch.tensor(params['defocus_v']), 
            torch.tensor(params['defocus_angle']),
            torch.tensor(params['voltage']),
            torch.tensor(params['cs']),
            torch.tensor(params['w'])
        )
        
        # CTF should be real-valued
        assert torch.is_floating_point(c), "CTF should be real-valued"
        
        # CTF should be finite everywhere
        assert torch.isfinite(c).all(), f"CTF contains non-finite values for params {params}"
        
        # CTF should have reasonable magnitude
        assert c.abs().max() <= 10.0, f"CTF values too large for params {params}"
    
    @given(
        params=ctf_parameters(),
        box_size=st.sampled_from([32, 64, 128])
    )
    @settings(max_examples=30, deadline=5000) 
    def test_ctf_scaling_property(self, params, box_size):
        """Test CTF behavior under frequency scaling."""
        device = torch.device('cpu')
        
        # Create lattices with different pixel sizes
        lattice1 = Lattice(box_size, device=device)
        freqs1 = lattice1.freqs2d / params['apix']
        
        # Scale frequencies by factor of 2
        freqs2 = freqs1 * 2.0
        
        c1 = ctf.compute_ctf(
            freqs1,
            torch.tensor(params['defocus_u']),
            torch.tensor(params['defocus_v']),
            torch.tensor(params['defocus_angle']),
            torch.tensor(params['voltage']),
            torch.tensor(params['cs']),
            torch.tensor(params['w'])
        )
        
        c2 = ctf.compute_ctf(
            freqs2,
            torch.tensor(params['defocus_u']),
            torch.tensor(params['defocus_v']),
            torch.tensor(params['defocus_angle']),
            torch.tensor(params['voltage']),
            torch.tensor(params['cs']),
            torch.tensor(params['w'])
        )
        
        # Both should be finite and well-behaved
        assert torch.isfinite(c1).all() and torch.isfinite(c2).all(), \
            "CTF should be finite for all valid parameters"


class TestLatticeProperties:
    """Property-based tests for Lattice operations."""
    
    @given(box_size=valid_box_sizes())
    @settings(max_examples=20, deadline=3000)
    def test_lattice_coordinate_properties(self, box_size):
        """Test basic properties of lattice coordinates."""
        device = torch.device('cpu')
        lattice = Lattice(box_size, device=device)
        
        # Check coordinate shapes
        assert lattice.coords.shape == (box_size * box_size, 3), \
            f"Unexpected coordinate shape for box size {box_size}"
        
        # Coordinates should be within expected range
        max_coord = box_size // 2
        assert lattice.coords.abs().max() <= max_coord, \
            f"Coordinates out of range for box size {box_size}"
        
        # Frequency coordinates should be properly scaled
        assert torch.isfinite(lattice.freqs2d).all(), \
            f"Non-finite frequency coordinates for box size {box_size}"
        
        # Check that we have the right number of unique coordinates
        unique_coords = torch.unique(lattice.coords, dim=0)
        assert unique_coords.shape[0] <= box_size * box_size, \
            f"Too many unique coordinates for box size {box_size}"
    
    @given(
        box_size=valid_box_sizes(),
        extent_fraction=st.floats(min_value=0.1, max_value=1.0)
    )
    @settings(max_examples=30, deadline=3000)
    def test_lattice_mask_properties(self, box_size, extent_fraction):
        """Test properties of lattice masks."""
        device = torch.device('cpu')
        extent = int(box_size * extent_fraction / 2)
        lattice = Lattice(box_size, extent=extent, device=device)
        
        mask = lattice.get_circular_mask(extent)
        
        # Mask should be boolean
        assert mask.dtype == torch.bool, "Mask should be boolean"
        
        # Mask should have correct shape
        assert mask.shape == (box_size * box_size,), \
            f"Unexpected mask shape for box size {box_size}"
        
        # Some points should be masked (unless extent is very large)
        if extent < box_size // 2:
            assert (~mask).any(), f"No points masked for extent {extent}, box size {box_size}"
        
        # Number of masked points should be reasonable
        masked_count = mask.sum().item()
        total_points = box_size * box_size
        mask_fraction = masked_count / total_points
        
        assert 0 <= mask_fraction <= 1, \
            f"Invalid mask fraction {mask_fraction} for extent {extent}, box size {box_size}"


class TestImageSourceProperties:
    """Property-based tests for ImageSource operations."""
    
    @given(
        n_particles=st.integers(min_value=1, max_value=100),
        box_size=st.sampled_from([32, 64, 128])
    )
    @settings(max_examples=20, deadline=10000)
    def test_image_source_consistency(self, n_particles, box_size):
        """Test that ImageSource provides consistent data access."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test data
            particles_data = np.random.randn(n_particles, box_size, box_size).astype(np.float32)
            particles_file = temp_path / "test_particles.mrcs"
            write_mrc(str(particles_file), particles_data, Apix=1.0)
            
            # Load with ImageSource
            source = ImageSource.from_file(str(particles_file), lazy=True)
            
            # Test properties
            assert source.n == n_particles, f"Incorrect particle count"
            assert source.D == box_size, f"Incorrect box size"
            
            # Test random access consistency
            test_indices = np.random.choice(n_particles, min(10, n_particles), replace=False)
            
            for idx in test_indices:
                # Access same image multiple times
                img1 = source.images([idx])
                img2 = source.images([idx])
                
                assert torch.allclose(img1, img2, rtol=1e-6, atol=1e-7), \
                    f"Inconsistent data access for particle {idx}"
                
                # Check shape
                assert img1.shape == (1, box_size, box_size), \
                    f"Unexpected image shape for particle {idx}"
    
    @given(
        n_particles=st.integers(min_value=10, max_value=50),
        chunk_size=st.integers(min_value=1, max_value=20)
    )
    @settings(max_examples=15, deadline=10000)
    def test_chunked_access_completeness(self, n_particles, chunk_size):
        """Test that chunked access covers all particles exactly once."""
        assume(chunk_size <= n_particles)  # Avoid trivial cases
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test data
            particles_data = np.random.randn(n_particles, 64, 64).astype(np.float32)
            particles_file = temp_path / "test_particles.mrcs"
            write_mrc(str(particles_file), particles_data, Apix=1.0)
            
            source = ImageSource.from_file(str(particles_file), lazy=True)
            
            # Collect all indices from chunked access
            all_indices = []
            total_images = 0
            
            for indices, chunk in source.chunks(chunksize=chunk_size):
                all_indices.extend(indices)
                total_images += chunk.shape[0]
                
                # Each chunk should have expected size (except possibly the last)
                assert chunk.shape[0] <= chunk_size, \
                    f"Chunk too large: {chunk.shape[0]} > {chunk_size}"
                assert chunk.shape[1:] == (64, 64), \
                    f"Unexpected chunk shape: {chunk.shape}"
            
            # Should cover all particles exactly once
            assert len(all_indices) == n_particles, \
                f"Wrong number of indices: {len(all_indices)} != {n_particles}"
            assert total_images == n_particles, \
                f"Wrong total images: {total_images} != {n_particles}"
            assert set(all_indices) == set(range(n_particles)), \
                f"Missing or duplicate indices in chunked access"


class TestUtilsProperties:
    """Property-based tests for utility functions."""
    
    @given(
        data=st.one_of(
            st.lists(st.floats(allow_nan=False, allow_infinity=False), min_size=1),
            arrays(np.float64, array_shapes(min_dims=1, max_dims=3), 
                  elements=st.floats(allow_nan=False, allow_infinity=False))
        )
    )
    @settings(max_examples=50, deadline=3000)
    def test_pickle_roundtrip_property(self, data):
        """Test that pickle save/load is a perfect roundtrip."""
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as temp_file:
            try:
                # Save data
                utils.save_pkl(data, temp_file.name)
                
                # Load data back
                loaded_data = utils.load_pkl(temp_file.name)
                
                # Should be identical
                if isinstance(data, np.ndarray):
                    np.testing.assert_array_equal(data, loaded_data)
                else:
                    assert data == loaded_data, f"Pickle roundtrip failed"
                
            finally:
                os.unlink(temp_file.name)
    
    @given(
        rotations=arrays(np.float64, (10, 3, 3), elements=st.floats(-2.0, 2.0)),
        translations=arrays(np.float64, (10, 2), elements=st.floats(-100.0, 100.0))
    )
    @settings(max_examples=30, deadline=3000)
    def test_pose_data_properties(self, rotations, translations):
        """Test properties of pose data handling."""
        # Ensure rotations are valid (orthogonal with det=1)
        for i in range(rotations.shape[0]):
            U, _, Vt = np.linalg.svd(rotations[i])
            rotations[i] = U @ Vt
            if np.linalg.det(rotations[i]) < 0:
                rotations[i][:, 0] *= -1
        
        pose_data = (rotations, translations)
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as temp_file:
            try:
                # Save and load pose data
                utils.save_pkl(pose_data, temp_file.name)
                loaded_rotations, loaded_translations = utils.load_pkl(temp_file.name)
                
                # Check properties are preserved
                np.testing.assert_allclose(rotations, loaded_rotations, rtol=1e-10)
                np.testing.assert_allclose(translations, loaded_translations, rtol=1e-10)
                
                # Rotations should still be valid
                for i in range(loaded_rotations.shape[0]):
                    R = loaded_rotations[i]
                    # Should be orthogonal
                    identity = R @ R.T
                    np.testing.assert_allclose(identity, np.eye(3), rtol=1e-10, atol=1e-12)
                    # Should have determinant 1
                    det = np.linalg.det(R)
                    assert abs(det - 1.0) < 1e-10, f"Invalid rotation determinant: {det}"
                    
            finally:
                os.unlink(temp_file.name)


if __name__ == "__main__":
    """Run property-based tests when called as a script."""
    pytest.main([__file__, "-v", "--hypothesis-show-statistics"])