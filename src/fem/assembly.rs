//! Finite element matrix assembly
//!
//! Provides both serial and parallel matrix assembly for FEM.
//! Parallel assembly uses scirs2-core's parallel features for
//! multi-threaded element-wise assembly.

use scirs2_core::{IntoParallelRefIterator, ParallelIterator};

use crate::fem::element::TriangularElement;
use crate::fem::mesh::Mesh2D;

/// Simple sparse matrix in triplet (COO) format
#[derive(Debug, Clone)]
pub struct SparseMatrix {
    /// Number of rows
    pub nrows: usize,
    /// Number of columns
    pub ncols: usize,
    /// Row indices
    pub rows: Vec<usize>,
    /// Column indices
    pub cols: Vec<usize>,
    /// Values
    pub vals: Vec<f64>,
}

impl SparseMatrix {
    /// Create new sparse matrix
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Self {
            nrows,
            ncols,
            rows: Vec::new(),
            cols: Vec::new(),
            vals: Vec::new(),
        }
    }

    /// Add entry to matrix
    pub fn add_entry(&mut self, i: usize, j: usize, val: f64) {
        self.rows.push(i);
        self.cols.push(j);
        self.vals.push(val);
    }

    /// Matrix-vector multiplication
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.nrows];
        for ((&i, &j), &val) in self.rows.iter().zip(&self.cols).zip(&self.vals) {
            y[i] += val * x[j];
        }
        y
    }

    /// Matrix-vector multiplication for a single row
    pub fn matvec_row(&self, x: &[f64], row: usize) -> f64 {
        let mut result = 0.0;
        for ((&i, &j), &val) in self.rows.iter().zip(&self.cols).zip(&self.vals) {
            if i == row {
                result += val * x[j];
            }
        }
        result
    }
}

/// Assemble global stiffness matrix for Laplacian operator (serial)
///
/// # Arguments
/// * `mesh` - Finite element mesh
///
/// # Returns
/// Sparse stiffness matrix
pub fn assemble_stiffness_matrix(mesh: &Mesh2D) -> SparseMatrix {
    let n = mesh.n_nodes();
    let mut k = SparseMatrix::new(n, n);

    // Assemble element contributions
    for elem in &mesh.elements {
        // Get element nodes
        let nodes = [
            mesh.nodes[elem.nodes[0]].position,
            mesh.nodes[elem.nodes[1]].position,
            mesh.nodes[elem.nodes[2]].position,
        ];

        let tri_elem = TriangularElement::new(nodes);
        let gradients = tri_elem.shape_gradients();
        let area = tri_elem.area();

        // Local stiffness matrix (3x3)
        for i in 0..3 {
            for j in 0..3 {
                let k_local =
                    area * (gradients[i][0] * gradients[j][0] + gradients[i][1] * gradients[j][1]);

                let global_i = elem.nodes[i];
                let global_j = elem.nodes[j];

                k.add_entry(global_i, global_j, k_local);
            }
        }
    }

    k
}

/// Assemble global stiffness matrix for Laplacian operator (parallel)
///
/// Uses parallel iteration over elements for faster assembly on multi-core systems.
///
/// # Arguments
/// * `mesh` - Finite element mesh
///
/// # Returns
/// Sparse stiffness matrix
pub fn assemble_stiffness_matrix_parallel(mesh: &Mesh2D) -> SparseMatrix {
    let n = mesh.n_nodes();

    // Compute element contributions in parallel
    let element_entries: Vec<Vec<(usize, usize, f64)>> = mesh
        .elements
        .par_iter()
        .map(|elem| {
            let nodes = [
                mesh.nodes[elem.nodes[0]].position,
                mesh.nodes[elem.nodes[1]].position,
                mesh.nodes[elem.nodes[2]].position,
            ];

            let tri_elem = TriangularElement::new(nodes);
            let gradients = tri_elem.shape_gradients();
            let area = tri_elem.area();

            let mut entries = Vec::with_capacity(9);

            for i in 0..3 {
                for j in 0..3 {
                    let k_local = area
                        * (gradients[i][0] * gradients[j][0] + gradients[i][1] * gradients[j][1]);

                    let global_i = elem.nodes[i];
                    let global_j = elem.nodes[j];

                    entries.push((global_i, global_j, k_local));
                }
            }

            entries
        })
        .collect();

    // Assemble into sparse matrix
    let mut k = SparseMatrix::new(n, n);
    for entries in element_entries {
        for (i, j, val) in entries {
            k.add_entry(i, j, val);
        }
    }

    k
}

/// Assemble global mass matrix (serial)
///
/// # Arguments
/// * `mesh` - Finite element mesh
///
/// # Returns
/// Sparse mass matrix
pub fn assemble_mass_matrix(mesh: &Mesh2D) -> SparseMatrix {
    let n = mesh.n_nodes();
    let mut m = SparseMatrix::new(n, n);

    for elem in &mesh.elements {
        let nodes = [
            mesh.nodes[elem.nodes[0]].position,
            mesh.nodes[elem.nodes[1]].position,
            mesh.nodes[elem.nodes[2]].position,
        ];

        let tri_elem = TriangularElement::new(nodes);
        let area = tri_elem.area();

        for i in 0..3 {
            for j in 0..3 {
                let m_local = if i == j { area / 6.0 } else { area / 12.0 };

                m.add_entry(elem.nodes[i], elem.nodes[j], m_local);
            }
        }
    }

    m
}

/// Assemble global mass matrix (parallel)
///
/// Uses parallel iteration for faster assembly on multi-core systems.
///
/// # Arguments
/// * `mesh` - Finite element mesh
///
/// # Returns
/// Sparse mass matrix
pub fn assemble_mass_matrix_parallel(mesh: &Mesh2D) -> SparseMatrix {
    let n = mesh.n_nodes();

    // Compute element contributions in parallel
    let element_entries: Vec<Vec<(usize, usize, f64)>> = mesh
        .elements
        .par_iter()
        .map(|elem| {
            let nodes = [
                mesh.nodes[elem.nodes[0]].position,
                mesh.nodes[elem.nodes[1]].position,
                mesh.nodes[elem.nodes[2]].position,
            ];

            let tri_elem = TriangularElement::new(nodes);
            let area = tri_elem.area();

            let mut entries = Vec::with_capacity(9);

            for i in 0..3 {
                for j in 0..3 {
                    let m_local = if i == j { area / 6.0 } else { area / 12.0 };
                    entries.push((elem.nodes[i], elem.nodes[j], m_local));
                }
            }

            entries
        })
        .collect();

    // Assemble into sparse matrix
    let mut m = SparseMatrix::new(n, n);
    for entries in element_entries {
        for (i, j, val) in entries {
            m.add_entry(i, j, val);
        }
    }

    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_assemble_stiffness() {
        let mesh = Mesh2D::rectangle(1.0, 1.0, 0.5).unwrap();
        let k = assemble_stiffness_matrix(&mesh);

        assert_eq!(k.nrows, mesh.n_nodes());
        assert_eq!(k.ncols, mesh.n_nodes());
        assert!(!k.vals.is_empty());
    }

    #[test]
    fn test_assemble_mass() {
        let mesh = Mesh2D::rectangle(1.0, 1.0, 0.5).unwrap();
        let m = assemble_mass_matrix(&mesh);

        assert_eq!(m.nrows, mesh.n_nodes());
        assert!(!m.vals.is_empty());
    }

    #[test]
    fn test_parallel_stiffness_assembly() {
        let mesh = Mesh2D::rectangle(1.0, 1.0, 0.25).unwrap();

        let k_serial = assemble_stiffness_matrix(&mesh);
        let k_parallel = assemble_stiffness_matrix_parallel(&mesh);

        // Both should have same dimensions
        assert_eq!(k_serial.nrows, k_parallel.nrows);
        assert_eq!(k_serial.ncols, k_parallel.ncols);

        // Should have same number of non-zero entries
        assert_eq!(k_serial.vals.len(), k_parallel.vals.len());
    }

    #[test]
    fn test_parallel_mass_assembly() {
        let mesh = Mesh2D::rectangle(1.0, 1.0, 0.25).unwrap();

        let m_serial = assemble_mass_matrix(&mesh);
        let m_parallel = assemble_mass_matrix_parallel(&mesh);

        // Both should have same dimensions
        assert_eq!(m_serial.nrows, m_parallel.nrows);
        assert_eq!(m_serial.ncols, m_parallel.ncols);

        // Should have same number of non-zero entries
        assert_eq!(m_serial.vals.len(), m_parallel.vals.len());
    }

    #[test]
    fn test_parallel_assembly_works() {
        // Test that parallel assembly produces a valid matrix
        let mesh = Mesh2D::rectangle(100e-9, 50e-9, 10e-9).unwrap();

        let k_parallel = assemble_stiffness_matrix_parallel(&mesh);
        let m_parallel = assemble_mass_matrix_parallel(&mesh);

        // Verify dimensions
        assert_eq!(k_parallel.nrows, mesh.n_nodes());
        assert_eq!(k_parallel.ncols, mesh.n_nodes());
        assert_eq!(m_parallel.nrows, mesh.n_nodes());
        assert_eq!(m_parallel.ncols, mesh.n_nodes());

        // Verify non-empty
        assert!(!k_parallel.vals.is_empty());
        assert!(!m_parallel.vals.is_empty());

        // Test that parallel assembly matches serial for this mesh
        let k_serial = assemble_stiffness_matrix(&mesh);
        assert_eq!(k_serial.vals.len(), k_parallel.vals.len());
    }
}
