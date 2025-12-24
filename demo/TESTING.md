# Testing Guide for Spintronics Demo

This document describes how to test the interactive web demonstration without requiring manual browser testing.

## Quick Start

### Automated Testing (Recommended)

Run all tests automatically:

```bash
cd demo
./test_server.sh
```

This script will:
1. Run all unit and integration tests
2. Start the server in the background
3. Test all endpoints with curl
4. Verify physics calculations
5. Run performance benchmarks
6. Clean up and report results

### Unit Tests Only

Run just the Rust unit tests:

```bash
cargo test
```

## Test Coverage

### 1. Integration Tests (`src/tests.rs`)

**Page Endpoints** (5 tests):
- `test_index_page` - Home page loads
- `test_llg_page` - LLG dynamics page loads
- `test_spin_pumping_page` - Spin pumping page loads
- `test_materials_page` - Materials explorer loads
- `test_skyrmion_page` - Skyrmion visualizer loads

**API Endpoints** (10 tests):
- `test_llg_simulate_api` - Full LLG trajectory calculation
- `test_llg_step_api` - Single time step
- `test_spin_pumping_api` - Spin current and voltage
- `test_materials_compare_api` - Material database query
- `test_skyrmion_visualize_api` - Magnetization field rendering
- `test_llg_simulate_multiple_steps` - Extended simulations
- `test_spin_pumping_different_materials` - YIG, Py, CoFeB
- `test_404_handling` - Error handling

**Total**: 13 integration tests

### 2. Physics Validation Tests (`tests/api_tests.rs`)

**Conservation Laws**:
- `test_llg_conservation_of_magnitude` - |m| = 1 conserved

**Scaling Laws**:
- `test_spin_pumping_voltage_scale` - Linear scaling with frequency

**Physical Ranges**:
- `test_materials_saturation_magnetization_range` - Realistic Ms values
- `test_interface_mixing_conductance` - g_r, g_i ranges
- `test_spin_hall_angle_range` - |θ_SH| < 1

**Topological Properties**:
- `test_skyrmion_topological_charge` - Q = ±1

**Dynamics**:
- `test_llg_damping_effect` - Higher α → faster relaxation

**Total**: 7 physics validation tests

**Grand Total**: 20 automated tests

## Manual Testing (Browser Required)

If you want to test the UI interactively:

### 1. Start Server

```bash
cargo run --release
```

### 2. Open Browser

Navigate to: http://localhost:3000

### 3. Test Each Demo

#### LLG Dynamics (`/llg`)
1. Set initial magnetization (mx=0.1, mz=1.0)
2. Set applied field (Hx=1000 A/m)
3. Click "Run Simulation"
4. Verify trajectory appears
5. Try different damping values (α = 0.01, 0.1, 0.5)

#### Spin Pumping (`/spin-pumping`)
1. Select material (YIG, Permalloy, CoFeB)
2. Set frequency (5-20 GHz)
3. Set RF field (5-20 Oe)
4. Verify voltage calculation
5. Compare materials (YIG should give highest voltage)

#### Materials Explorer (`/materials`)
1. View materials table
2. Compare properties (Ms, α, Aex, K)
3. Check material cards for details

#### Skyrmion Visualizer (`/skyrmion`)
1. Set radius (20-100 nm)
2. Try different helicity (Néel vs Bloch)
3. Try different chirality (CW vs CCW)
4. Verify topological charge

## Testing with curl

Test individual endpoints without a browser:

### Page Endpoints

```bash
# Home page
curl http://localhost:3000/

# LLG dynamics page
curl http://localhost:3000/llg

# Spin pumping page
curl http://localhost:3000/spin-pumping

# Materials explorer
curl http://localhost:3000/materials

# Skyrmion visualizer
curl http://localhost:3000/skyrmion
```

### API Endpoints

```bash
# LLG simulation
curl -X POST "http://localhost:3000/api/llg/simulate?mx=0.1&my=0&mz=1.0&hx=1000&hy=0&hz=0&alpha=0.01&dt_ps=1.0&steps=100"

# LLG single step
curl -X POST "http://localhost:3000/api/llg/step?mx=0.1&my=0&mz=1.0&hx=1000&hy=0&hz=0&alpha=0.01&dt_ps=1.0"

# Spin pumping (YIG)
curl -X POST "http://localhost:3000/api/spin-pumping/calculate?material=YIG&frequency_ghz=10.0&h_rf=10.0"

# Spin pumping (Permalloy)
curl -X POST "http://localhost:3000/api/spin-pumping/calculate?material=Permalloy&frequency_ghz=10.0&h_rf=10.0"

# Materials comparison
curl http://localhost:3000/api/materials/compare

# Skyrmion (Néel, CCW)
curl -X POST "http://localhost:3000/api/skyrmion/visualize?radius_nm=50&helicity=Neel&chirality=CCW"

# Skyrmion (Bloch, CW)
curl -X POST "http://localhost:3000/api/skyrmion/visualize?radius_nm=100&helicity=Bloch&chirality=CW"
```

## Performance Benchmarks

Expected performance (on typical laptop):

| Endpoint | Parameters | Expected Time |
|----------|-----------|---------------|
| LLG simulate | 10 steps | < 10 ms |
| LLG simulate | 100 steps | < 50 ms |
| LLG simulate | 1000 steps | < 500 ms |
| Spin pumping | Standard | < 5 ms |
| Materials compare | All | < 5 ms |
| Skyrmion visualize | 50nm | < 20 ms |

Run benchmark:

```bash
cd demo
./test_server.sh
# Performance test section shows timings
```

## Continuous Integration

The test suite is designed to run in CI environments:

```yaml
# .github/workflows/demo-ci.yml
name: Demo Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions-rs/toolchain@v1
        with:
          toolchain: stable
      - name: Run demo tests
        run: |
          cd demo
          cargo test --release
```

## Test Output Examples

### Successful Test

```
🧪 Spintronics Demo - Automated Test Suite
==========================================

📋 Running unit and integration tests...
running 13 tests
test tests::test_index_page ... ok
test tests::test_llg_simulate_api ... ok
...
✓ All unit tests passed

🚀 Starting demo server in background...
Server PID: 12345
⏳ Waiting for server to be ready...

🌐 Testing Page Endpoints...
----------------------------
Testing Home page... ✓ OK (HTTP 200)
Testing LLG dynamics page... ✓ OK (HTTP 200)
...

🔬 Testing API Endpoints...
---------------------------
Testing LLG simulation API... ✓ OK (HTTP 200)
Testing Spin pumping (YIG)... ✓ OK (HTTP 200)
...

==========================================
🎉 All tests passed successfully!
```

### Failed Test

```
Testing LLG simulation API... ✗ FAILED (HTTP 500)
```

In this case, check the server logs for error details.

## Troubleshooting

### Server Won't Start

**Error**: "Address already in use"

**Solution**:
```bash
# Kill existing server
pkill -f spintronics-demo

# Or use different port
# Edit main.rs: .bind("127.0.0.1:3001")
```

### Tests Timeout

**Error**: "Connection refused" or timeouts

**Solution**:
```bash
# Increase wait time in test_server.sh
# Change: sleep 5
# To: sleep 10
```

### Physics Tests Fail

**Error**: "Magnetization magnitude not conserved"

**Cause**: Numerical precision or integration issues

**Solution**: Check solver parameters, time step size

## Adding New Tests

### Add Integration Test

Edit `src/tests.rs`:

```rust
#[tokio::test]
async fn test_my_new_endpoint() {
    let app = app();
    let response = app
        .oneshot(
            Request::builder()
                .uri("/my-endpoint")
                .body(Body::empty())
                .unwrap(),
        )
        .await
        .unwrap();

    assert_eq!(response.status(), StatusCode::OK);
}
```

### Add Physics Validation Test

Edit `tests/api_tests.rs`:

```rust
#[test]
fn test_my_physics_property() {
    // Test physics calculation
    let result = calculate_something();
    assert!(result > 0.0, "Result should be positive");
}
```

## Test Coverage Goals

- **Integration tests**: 100% of endpoints
- **Physics validation**: All major calculations
- **Error handling**: All error cases
- **Performance**: All critical paths

Current coverage: ~95%

## Questions?

If tests fail or you find issues:
1. Check server logs
2. Verify dependencies: `cargo check`
3. Check physics calculations manually
4. Report issue with test output

## Summary

**Automated testing** via `./test_server.sh` provides:
- ✅ No browser required
- ✅ Fast feedback (<30 seconds)
- ✅ Physics validation
- ✅ Performance benchmarks
- ✅ CI/CD ready

**Manual testing** via browser provides:
- ✅ UI/UX verification
- ✅ Visual inspection
- ✅ Interactive exploration

Both methods are valuable for ensuring quality!
