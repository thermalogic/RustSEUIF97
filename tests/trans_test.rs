#![allow(warnings)]
use assert_approx_eq::assert_approx_eq;
use if97::common::transport_further::surface_tension;
use if97::common::transport_further::thcond;
use if97::common::transport_further::viscosity;
use if97::common::*;

#[test]
fn test_viscosity() {
    // Table 4, page 8 : T rho(d) mu(upa.s)
    const data: [[f64; 3]; 11] = [
        [298.15, 998.0, 889.735100],
        [298.15, 1200.0, 1437.649467],
        [373.15, 1000.0, 307.883622],
        [433.15, 1.0, 14.538324],
        [433.15, 1000.0, 217.685358],
        [873.15, 1.0, 32.619287],
        [873.15, 100.0, 35.802262],
        [873.15, 600.0, 77.430195],
        [1173.15, 1.0, 44.217245],
        [1173.15, 100.0, 47.640433],
        [1173.15, 400.0, 64.15460],
    ];
    for i in 0..11 {
        assert_approx_eq!(
            data[i][2],
            viscosity(data[i][1], data[i][0]) * 1.0e6,
            1.0e-5f64
        );
    }
}

#[test]
fn test_thcond() {
    //  Table 4, page 10  T ruo k (mW·m−1·K−1)
    const data: [[f64; 3]; 4] = [
        [298.15, 0.0, 18.4341883],
        [298.15, 998.0, 607.712868],
        [298.15, 1200.0, 799.038144],
        [873.15, 0.0, 79.1034659],
    ];
    for i in 0..4 {
        assert_approx_eq!(
            data[i][2],
            thcond(data[i][1], data[i][0]) * 1.0e3,
            1.0e-5f64
        );
    }
}

#[test]
fn test_surface_tension() {
    // Selected values from table 1 page 4 t  °C  st mNm
    const data: [[f64; 2]; 6] = [
        [0.01, 75.65],
        [5.0, 74.94],
        [10.0, 74.22],
        [15.0, 74.39],
        [20.0, 72.74],
        [25.0, 71.97],
    ];
    for i in 0..6 {
        assert_approx_eq!(
            data[i][1],
            surface_tension(data[i][0] + 273.15) * 1.0e3,
            1.0e0f64
        );
    }
}

#[test]
fn test_static_dielectric() {
    // rho kg/m³, T K , epsion: Dielectric constant [-]
    const data: [[f64; 3]; 2] = [
        [999.242866, 298.15, 78.5907250],
        [26.0569558, 873.15, 1.12620970],
    ];
    for i in 0..2 {
        assert_approx_eq!(
            data[i][2],
            static_dielectric(data[i][0], data[i][1]),
            1.0e-5f64
        );
    }
}
