//! Demo implementations for HTMX endpoints

use axum::extract::Query;
use axum::response::{Html, IntoResponse};
use serde::{Deserialize, Serialize};
use spintronics::dynamics::{zeeman_energy, LlgSolver};
use spintronics::prelude::*;

/// LLG simulation request
#[derive(Debug, Deserialize)]
pub struct LlgRequest {
    /// Initial magnetization direction
    #[serde(default = "default_mx")]
    pub mx: f64,
    #[serde(default)]
    pub my: f64,
    #[serde(default = "default_mz")]
    pub mz: f64,
    /// Applied field (A/m)
    #[serde(default = "default_field")]
    pub hx: f64,
    #[serde(default)]
    pub hy: f64,
    #[serde(default)]
    pub hz: f64,
    /// Damping constant
    #[serde(default = "default_alpha")]
    pub alpha: f64,
    /// Time step (ps)
    #[serde(default = "default_dt")]
    pub dt_ps: f64,
    /// Number of steps
    #[serde(default = "default_steps")]
    pub steps: usize,
}

fn default_mx() -> f64 {
    0.1
}
fn default_mz() -> f64 {
    1.0
}
fn default_field() -> f64 {
    1000.0
}
fn default_alpha() -> f64 {
    0.01
}
fn default_dt() -> f64 {
    1.0
}
fn default_steps() -> usize {
    100
}

/// LLG simulation response (for future JSON API)
#[allow(dead_code)]
#[derive(Debug, Serialize)]
pub struct LlgResponse {
    pub trajectory: Vec<MagPoint>,
    pub energy: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct MagPoint {
    pub time: f64,
    pub mx: f64,
    pub my: f64,
    pub mz: f64,
}

/// Run LLG simulation and return trajectory as HTML
pub async fn llg_simulate(Query(params): Query<LlgRequest>) -> impl IntoResponse {
    // Create ferromagnet
    let fm = Ferromagnet::permalloy().with_alpha(params.alpha);

    // Initial magnetization
    let m0 = Vector3::new(params.mx, params.my, params.mz).normalize() * fm.ms;

    // Applied field
    let h_ext = Vector3::new(params.hx, params.hy, params.hz);

    // Time step in seconds
    let dt = params.dt_ps * 1e-12;

    // Run simulation
    let solver = LlgSolver::new(fm.alpha, dt);
    let mut m = m0;
    let mut trajectory = Vec::with_capacity(params.steps);
    let mut energy = Vec::with_capacity(params.steps);

    for i in 0..params.steps {
        let m_norm = m * (1.0 / fm.ms);
        trajectory.push(MagPoint {
            time: i as f64 * params.dt_ps,
            mx: m_norm.x,
            my: m_norm.y,
            mz: m_norm.z,
        });

        // Calculate Zeeman energy
        let e = zeeman_energy(m_norm, h_ext, fm.ms);
        energy.push(e);

        // Step forward
        let h_eff = h_ext; // Simplified: only external field
        m = solver.step_rk4(m, |_m| h_eff);
    }

    // Generate HTML fragment for HTMX
    let html = generate_trajectory_html(&trajectory);
    Html(html)
}

/// Generate single LLG step (for real-time updates)
pub async fn llg_step(Query(params): Query<LlgRequest>) -> impl IntoResponse {
    let fm = Ferromagnet::permalloy().with_alpha(params.alpha);
    let m0 = Vector3::new(params.mx, params.my, params.mz).normalize() * fm.ms;
    let h_ext = Vector3::new(params.hx, params.hy, params.hz);
    let dt = params.dt_ps * 1e-12;

    let solver = LlgSolver::new(fm.alpha, dt);
    let m_new = solver.step_rk4(m0, |_m| h_ext);
    let m_norm = m_new * (1.0 / fm.ms);

    // Return HTML fragment showing current magnetization
    let html = format!(
        r#"<div class="magnetization" id="mag-display">
            <div class="mag-component">
                <span class="label">mx:</span>
                <span class="value">{:.4}</span>
                <div class="bar" style="width: {}%"></div>
            </div>
            <div class="mag-component">
                <span class="label">my:</span>
                <span class="value">{:.4}</span>
                <div class="bar" style="width: {}%"></div>
            </div>
            <div class="mag-component">
                <span class="label">mz:</span>
                <span class="value">{:.4}</span>
                <div class="bar" style="width: {}%"></div>
            </div>
        </div>"#,
        m_norm.x,
        (m_norm.x.abs() * 100.0),
        m_norm.y,
        (m_norm.y.abs() * 100.0),
        m_norm.z,
        (m_norm.z.abs() * 100.0)
    );

    Html(html)
}

fn generate_trajectory_html(trajectory: &[MagPoint]) -> String {
    let mut html = String::from(r##"<div class="trajectory-result">"##);
    html.push_str("<h3>Magnetization Trajectory</h3>");
    html.push_str(r##"<svg width="600" height="400" viewBox="-1.2 -1.2 2.4 2.4" style="border: 1px solid #ccc; background: #f9f9f9;">"##);

    // Draw unit sphere
    html.push_str(
        r##"<circle cx="0" cy="0" r="1" fill="none" stroke="#ddd" stroke-width="0.02"/>"##,
    );

    // Draw trajectory
    html.push_str(r##"<path d="M "##);
    for (i, point) in trajectory.iter().enumerate() {
        if i == 0 {
            html.push_str(&format!("{:.3} {:.3} ", point.mx, -point.mz));
        } else {
            html.push_str(&format!("L {:.3} {:.3} ", point.mx, -point.mz));
        }
    }
    html.push_str(r##"" fill="none" stroke="#0066cc" stroke-width="0.02" opacity="0.8"/>"##);

    // Mark start and end
    if let Some(first) = trajectory.first() {
        html.push_str(&format!(
            r##"<circle cx="{:.3}" cy="{:.3}" r="0.05" fill="green"/>"##,
            first.mx, -first.mz
        ));
    }
    if let Some(last) = trajectory.last() {
        html.push_str(&format!(
            r##"<circle cx="{:.3}" cy="{:.3}" r="0.05" fill="red"/>"##,
            last.mx, -last.mz
        ));
    }

    html.push_str("</svg>");
    html.push_str(&format!("<p>Simulated {} steps</p>", trajectory.len()));
    html.push_str("</div>");
    html
}

/// Spin pumping calculation
#[derive(Debug, Deserialize)]
pub struct SpinPumpingRequest {
    #[serde(default = "default_frequency")]
    pub frequency_ghz: f64,
    #[serde(default = "default_h_rf")]
    pub h_rf: f64,
    #[serde(default = "default_material")]
    pub material: String,
}

fn default_frequency() -> f64 {
    9.0
}
fn default_h_rf() -> f64 {
    100.0
}
fn default_material() -> String {
    "YIG".to_string()
}

pub async fn spin_pumping_calculate(Query(params): Query<SpinPumpingRequest>) -> impl IntoResponse {
    // Select material
    let fm = match params.material.as_str() {
        "YIG" => Ferromagnet::yig(),
        "Permalloy" => Ferromagnet::permalloy(),
        "CoFeB" => Ferromagnet::cofeb(),
        _ => Ferromagnet::yig(),
    };

    // Interface
    let interface = match params.material.as_str() {
        "YIG" => SpinInterface::yig_pt(),
        "Permalloy" => SpinInterface::py_pt(),
        "CoFeB" => SpinInterface::cofeb_pt(),
        _ => SpinInterface::yig_pt(),
    };

    // Frequency in Hz
    let f = params.frequency_ghz * 1e9;
    let omega = 2.0 * std::f64::consts::PI * f;

    // Precession cone angle (approximate)
    let theta = (params.h_rf / 1000.0).min(0.5); // Simplified

    // Simplified spin current calculation (based on mixing conductance theory)
    // j_s ∝ g_r * ω * θ² * sin(θ) for small cone angles
    // This is a simplified approximation of the full spin pumping formula
    use spintronics::constants::{E_CHARGE, HBAR};
    let js_magnitude = (HBAR / (4.0 * std::f64::consts::PI * E_CHARGE))
        * interface.g_r
        * omega
        * theta.sin()
        * theta; // Simplified: ∝ ω θ sin(θ)

    // ISHE conversion to voltage
    let ishe = InverseSpinHall::platinum();
    let js_vec = Vector3::new(0.0, 0.0, js_magnitude);
    let js_pol = Vector3::new(1.0, 0.0, 0.0);
    let voltage = ishe.voltage(js_vec, js_pol, 5e-3); // 5mm width

    // Generate HTML response
    let html = format!(
        r#"<div class="result-card" id="spin-pumping-result">
            <h3>Spin Pumping Results</h3>
            <div class="result-grid">
                <div class="result-item">
                    <span class="label">Frequency:</span>
                    <span class="value">{:.2} GHz</span>
                </div>
                <div class="result-item">
                    <span class="label">RF Field:</span>
                    <span class="value">{:.1} A/m</span>
                </div>
                <div class="result-item">
                    <span class="label">Spin Current:</span>
                    <span class="value">{:.2e} A/m²</span>
                </div>
                <div class="result-item">
                    <span class="label">ISHE Voltage:</span>
                    <span class="value">{:.2} μV</span>
                </div>
                <div class="result-item">
                    <span class="label">Material:</span>
                    <span class="value">{}</span>
                </div>
            </div>
            <div class="material-info">
                <p>M<sub>s</sub> = {:.2e} A/m, α = {:.4}</p>
                <p>g<sub>r</sub> = {:.2e} Ω⁻¹m⁻², θ<sub>SH</sub> = {:.3}</p>
            </div>
        </div>"#,
        params.frequency_ghz,
        params.h_rf,
        js_magnitude,
        voltage * 1e6, // Convert to μV
        params.material,
        fm.ms,
        fm.alpha,
        interface.g_r,
        ishe.theta_sh
    );

    Html(html)
}

/// Materials comparison request
#[derive(Debug, Deserialize)]
pub struct MaterialsRequest {
    #[allow(dead_code)]
    pub materials: Option<String>,
}

pub async fn materials_compare(Query(_params): Query<MaterialsRequest>) -> impl IntoResponse {
    let materials = vec![
        ("YIG", Ferromagnet::yig()),
        ("Permalloy", Ferromagnet::permalloy()),
        ("CoFeB", Ferromagnet::cofeb()),
        ("Iron", Ferromagnet::iron()),
        ("Cobalt", Ferromagnet::cobalt()),
        ("Nickel", Ferromagnet::nickel()),
    ];

    let mut html = String::from(r#"<div class="materials-table" id="materials-comparison">"#);
    html.push_str("<table>");
    html.push_str("<thead><tr>");
    html.push_str("<th>Material</th><th>M<sub>s</sub> (MA/m)</th><th>α</th><th>A<sub>ex</sub> (pJ/m)</th><th>K (kJ/m³)</th>");
    html.push_str("</tr></thead><tbody>");

    for (name, mat) in materials {
        html.push_str(&format!(
            "<tr><td>{}</td><td>{:.2}</td><td>{:.4}</td><td>{:.1}</td><td>{:.1}</td></tr>",
            name,
            mat.ms / 1e6,
            mat.alpha,
            mat.exchange_a * 1e12,
            mat.anisotropy_k / 1e3
        ));
    }

    html.push_str("</tbody></table></div>");
    Html(html)
}

/// Skyrmion visualization request
#[derive(Debug, Deserialize)]
pub struct SkyrmionRequest {
    #[serde(default = "default_radius")]
    pub radius_nm: f64,
    #[serde(default)]
    pub helicity: String,
    #[serde(default)]
    pub chirality: String,
}

fn default_radius() -> f64 {
    50.0
}

pub async fn skyrmion_visualize(Query(params): Query<SkyrmionRequest>) -> impl IntoResponse {
    let radius = params.radius_nm * 1e-9;

    let helicity = match params.helicity.as_str() {
        "Bloch" => Helicity::Bloch,
        _ => Helicity::Neel,
    };

    let chirality = match params.chirality.as_str() {
        "CW" => Chirality::Clockwise,
        _ => Chirality::CounterClockwise,
    };

    let sk = Skyrmion::new((0.0, 0.0), radius, helicity, chirality);

    // Generate SVG visualization
    let size = 300;
    let grid_size = 50;
    let mut html = format!(
        r#"<div class="skyrmion-viz" id="skyrmion-display">
            <svg width="{}" height="{}" viewBox="0 0 {} {}" style="border: 1px solid #ccc;">"#,
        size, size, grid_size, grid_size
    );

    // Draw magnetization field
    for i in 0..grid_size {
        for j in 0..grid_size {
            let x = (i as f64 - grid_size as f64 / 2.0) * radius * 4.0 / grid_size as f64;
            let y = (j as f64 - grid_size as f64 / 2.0) * radius * 4.0 / grid_size as f64;

            let m = sk.magnetization_at(x, y, 10e-9);
            let m_norm = m.normalize();

            // Color based on mz
            let color_val = ((m_norm.z + 1.0) / 2.0 * 255.0) as u8;
            let color = format!("rgb({}, {}, {})", color_val, 100, 255 - color_val);

            html.push_str(&format!(
                r#"<rect x="{}" y="{}" width="1" height="1" fill="{}"/>"#,
                i, j, color
            ));
        }
    }

    html.push_str("</svg>");
    html.push_str(&format!(
        r#"<div class="skyrmion-info">
            <p>Radius: {:.1} nm</p>
            <p>Helicity: {:?}</p>
            <p>Chirality: {:?}</p>
            <p>Topological Charge: {}</p>
        </div></div>"#,
        params.radius_nm, helicity, chirality, sk.topological_charge
    ));

    Html(html)
}
