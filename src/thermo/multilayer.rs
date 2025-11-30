//! Thermal transport in multilayer structures
//!
//! Models heat flow and spin-thermal effects in multilayer thin films,
//! accounting for interface thermal resistance (Kapitza resistance).

/// Single layer in a multilayer stack
#[derive(Debug, Clone)]
pub struct Layer {
    /// Layer name/identifier
    pub name: String,

    /// Thickness \[m\]
    pub thickness: f64,

    /// Thermal conductivity \[W/(m·K)\]
    pub thermal_conductivity: f64,

    /// Volumetric heat capacity [J/(m³·K)]
    pub heat_capacity: f64,

    /// Electrical resistivity [Ω·m]
    pub resistivity: f64,

    /// Is this layer magnetic?
    pub is_magnetic: bool,
}

impl Layer {
    /// Create a YIG layer
    pub fn yig(thickness: f64) -> Self {
        Self {
            name: "YIG".to_string(),
            thickness,
            thermal_conductivity: 6.0, // W/(m·K)
            heat_capacity: 3.0e6,      // J/(m³·K)
            resistivity: 1.0e10,       // Insulator
            is_magnetic: true,
        }
    }

    /// Create a Pt layer
    pub fn pt(thickness: f64) -> Self {
        Self {
            name: "Pt".to_string(),
            thickness,
            thermal_conductivity: 72.0, // W/(m·K)
            heat_capacity: 2.8e6,       // J/(m³·K)
            resistivity: 2.0e-7,        // Ω·m
            is_magnetic: false,
        }
    }

    /// Create a Cu layer
    pub fn cu(thickness: f64) -> Self {
        Self {
            name: "Cu".to_string(),
            thickness,
            thermal_conductivity: 400.0, // W/(m·K)
            heat_capacity: 3.4e6,        // J/(m³·K)
            resistivity: 1.7e-8,
            is_magnetic: false,
        }
    }

    /// Calculate thermal resistance
    pub fn thermal_resistance(&self, area: f64) -> f64 {
        self.thickness / (self.thermal_conductivity * area)
    }
}

/// Thermal boundary between two layers
#[derive(Debug, Clone)]
pub struct ThermalBoundary {
    /// Interface thermal conductance [W/(m²·K)]
    ///
    /// Also called Kapitza conductance. Inverse is Kapitza resistance.
    /// Typical values: 10⁷ - 10⁹ W/(m²·K)
    pub conductance: f64,
}

impl ThermalBoundary {
    /// Create boundary with default conductance
    pub fn new(conductance: f64) -> Self {
        Self { conductance }
    }

    /// Typical metal/insulator interface
    pub fn metal_insulator() -> Self {
        Self {
            conductance: 1.0e8, // W/(m²·K)
        }
    }

    /// Typical metal/metal interface
    pub fn metal_metal() -> Self {
        Self { conductance: 5.0e8 }
    }

    /// Calculate thermal resistance of interface
    pub fn thermal_resistance(&self, area: f64) -> f64 {
        1.0 / (self.conductance * area)
    }

    /// Calculate heat flux across boundary
    ///
    /// Q = G × ΔT
    pub fn heat_flux(&self, delta_t: f64) -> f64 {
        self.conductance * delta_t
    }
}

/// Multilayer thin film stack
#[derive(Debug, Clone)]
pub struct MultilayerStack {
    /// Layers from bottom to top
    pub layers: Vec<Layer>,

    /// Boundaries between layers (n-1 boundaries for n layers)
    pub boundaries: Vec<ThermalBoundary>,

    /// Cross-sectional area \[m²\]
    pub area: f64,
}

impl MultilayerStack {
    /// Create a new multilayer stack
    pub fn new(area: f64) -> Self {
        Self {
            layers: Vec::new(),
            boundaries: Vec::new(),
            area,
        }
    }

    /// Add a layer to the top of the stack
    pub fn add_layer(&mut self, layer: Layer, boundary: Option<ThermalBoundary>) {
        if !self.layers.is_empty() {
            let boundary = boundary.unwrap_or_else(|| {
                // Auto-detect boundary type
                let prev_layer = self.layers.last().unwrap();
                if prev_layer.is_magnetic != layer.is_magnetic {
                    ThermalBoundary::metal_insulator()
                } else {
                    ThermalBoundary::metal_metal()
                }
            });
            self.boundaries.push(boundary);
        }
        self.layers.push(layer);
    }

    /// Calculate total thermal resistance
    pub fn total_thermal_resistance(&self) -> f64 {
        let mut r_total = 0.0;

        // Add layer resistances
        for layer in &self.layers {
            r_total += layer.thermal_resistance(self.area);
        }

        // Add interface resistances
        for boundary in &self.boundaries {
            r_total += boundary.thermal_resistance(self.area);
        }

        r_total
    }

    /// Calculate effective thermal conductivity
    ///
    /// k_eff = L / R_total
    pub fn effective_thermal_conductivity(&self) -> f64 {
        let total_thickness: f64 = self.layers.iter().map(|l| l.thickness).sum();
        let r_total = self.total_thermal_resistance();

        if r_total > 0.0 {
            total_thickness / (r_total * self.area)
        } else {
            0.0
        }
    }

    /// Calculate temperature profile through the stack
    ///
    /// # Arguments
    /// * `t_bottom` - Temperature at bottom surface \[K\]
    /// * `t_top` - Temperature at top surface \[K\]
    ///
    /// # Returns
    /// Vector of temperatures at each interface (including surfaces)
    pub fn temperature_profile(&self, t_bottom: f64, t_top: f64) -> Vec<f64> {
        let mut temperatures = vec![t_bottom];

        let r_total = self.total_thermal_resistance();
        let delta_t = t_top - t_bottom;

        let mut current_t = t_bottom;

        for i in 0..self.layers.len() {
            // Temperature drop across layer
            let r_layer = self.layers[i].thermal_resistance(self.area);
            let dt_layer = delta_t * (r_layer / r_total);
            current_t += dt_layer;
            temperatures.push(current_t);

            // Temperature drop across interface (if not last layer)
            if i < self.boundaries.len() {
                let r_interface = self.boundaries[i].thermal_resistance(self.area);
                let dt_interface = delta_t * (r_interface / r_total);
                current_t += dt_interface;
                temperatures.push(current_t);
            }
        }

        temperatures
    }

    /// Get total thickness
    pub fn total_thickness(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_layer_creation() {
        let yig = Layer::yig(100.0e-9);
        assert_eq!(yig.name, "YIG");
        assert!(yig.is_magnetic);

        let pt = Layer::pt(5.0e-9);
        assert_eq!(pt.name, "Pt");
        assert!(!pt.is_magnetic);
    }

    #[test]
    fn test_multilayer_stack() {
        let mut stack = MultilayerStack::new(1.0e-12); // 1 μm²

        stack.add_layer(Layer::yig(100.0e-9), None);
        stack.add_layer(Layer::pt(5.0e-9), None);

        assert_eq!(stack.layers.len(), 2);
        assert_eq!(stack.boundaries.len(), 1);
        assert!((stack.total_thickness() - 105.0e-9).abs() < 1e-15);
    }

    #[test]
    fn test_thermal_resistance() {
        let mut stack = MultilayerStack::new(1.0e-12);
        stack.add_layer(Layer::pt(10.0e-9), None);

        let r = stack.total_thermal_resistance();
        assert!(r > 0.0);
    }

    #[test]
    fn test_temperature_profile() {
        let mut stack = MultilayerStack::new(1.0e-12);
        stack.add_layer(Layer::yig(100.0e-9), None);
        stack.add_layer(Layer::pt(10.0e-9), None);

        let profile = stack.temperature_profile(300.0, 310.0);

        // First temperature is bottom
        assert!((profile[0] - 300.0).abs() < 1e-10);

        // Temperatures should increase monotonically
        for i in 1..profile.len() {
            assert!(profile[i] >= profile[i - 1]);
        }

        // Last temperature is top
        assert!((profile.last().unwrap() - 310.0).abs() < 1e-6);
    }
}
