//! Mesh generation and data structures for FEM
//!
//! Uses Delaunay triangulation from scirs2-spatial for generating
//! high-quality finite element meshes.

use crate::vector3::Vector3;
use scirs2_core::ndarray::Array2;
use scirs2_spatial::delaunay::Delaunay;

/// Node in the finite element mesh
#[derive(Debug, Clone)]
pub struct Node {
    /// Node ID
    pub id: usize,
    /// Spatial coordinates \[m\]
    pub position: Vector3<f64>,
}

/// 2D triangular element
#[derive(Debug, Clone)]
pub struct Element2D {
    /// Element ID
    pub id: usize,
    /// Node IDs (3 for triangle)
    pub nodes: [usize; 3],
}

/// 3D tetrahedral element
#[derive(Debug, Clone)]
pub struct Element3D {
    /// Element ID
    pub id: usize,
    /// Node IDs (4 for tetrahedron)
    pub nodes: [usize; 4],
}

/// 2D finite element mesh
#[derive(Debug, Clone)]
pub struct Mesh2D {
    /// Nodes in the mesh
    pub nodes: Vec<Node>,
    /// Elements in the mesh
    pub elements: Vec<Element2D>,
    /// Characteristic mesh size \[m\]
    pub h: f64,
}

impl Mesh2D {
    /// Create a rectangular mesh using Delaunay triangulation
    ///
    /// # Arguments
    /// * `width` - Width in x-direction \[m\]
    /// * `height` - Height in y-direction \[m\]
    /// * `h` - Target element size \[m\]
    ///
    /// # Returns
    /// Triangulated mesh of the rectangle
    pub fn rectangle(width: f64, height: f64, h: f64) -> Result<Self, String> {
        // Generate regular grid points
        let nx = ((width / h) as usize).max(2);
        let ny = ((height / h) as usize).max(2);

        // Generate points as Array2
        let n_points = (nx + 1) * (ny + 1);
        let mut points = Array2::zeros((n_points, 2));
        let mut idx = 0;
        for j in 0..=ny {
            for i in 0..=nx {
                let x = (i as f64) * width / (nx as f64);
                let y = (j as f64) * height / (ny as f64);
                points[[idx, 0]] = x;
                points[[idx, 1]] = y;
                idx += 1;
            }
        }

        // Perform Delaunay triangulation
        let delaunay = Delaunay::new(&points)
            .map_err(|e| format!("Delaunay triangulation failed: {:?}", e))?;

        // Extract nodes
        let nodes: Vec<Node> = (0..n_points)
            .map(|id| Node {
                id,
                position: Vector3::new(points[[id, 0]], points[[id, 1]], 0.0),
            })
            .collect();

        // Extract simplices (triangles)
        let simplices = delaunay.simplices();
        let elements: Vec<Element2D> = simplices
            .iter()
            .enumerate()
            .map(|(id, tri)| Element2D {
                id,
                nodes: [tri[0], tri[1], tri[2]],
            })
            .collect();

        Ok(Self { nodes, elements, h })
    }

    /// Get number of nodes
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Get number of elements
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }

    /// Get element area
    pub fn element_area(&self, elem_id: usize) -> f64 {
        let elem = &self.elements[elem_id];
        let p0 = &self.nodes[elem.nodes[0]].position;
        let p1 = &self.nodes[elem.nodes[1]].position;
        let p2 = &self.nodes[elem.nodes[2]].position;

        // Area = 0.5 * |cross product|
        let v1 = *p1 - *p0;
        let v2 = *p2 - *p0;
        0.5 * v1.cross(&v2).magnitude()
    }
}

/// 3D finite element mesh
#[derive(Debug, Clone)]
pub struct Mesh3D {
    /// Nodes in the mesh
    pub nodes: Vec<Node>,
    /// Elements in the mesh
    pub elements: Vec<Element3D>,
    /// Characteristic mesh size \[m\]
    pub h: f64,
}

impl Mesh3D {
    /// Create a cubic mesh
    ///
    /// # Arguments
    /// * `lx, ly, lz` - Dimensions \[m\]
    /// * `h` - Target element size \[m\]
    pub fn cuboid(lx: f64, ly: f64, lz: f64, h: f64) -> Result<Self, String> {
        // Generate regular grid points
        let nx = ((lx / h) as usize).max(2);
        let ny = ((ly / h) as usize).max(2);
        let nz = ((lz / h) as usize).max(2);

        let mut nodes = Vec::new();
        let mut node_id = 0;

        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    let x = (i as f64) * lx / (nx as f64);
                    let y = (j as f64) * ly / (ny as f64);
                    let z = (k as f64) * lz / (nz as f64);
                    nodes.push(Node {
                        id: node_id,
                        position: Vector3::new(x, y, z),
                    });
                    node_id += 1;
                }
            }
        }

        // Create structured tetrahedral mesh (divide each cube into 5 tetrahedra)
        let mut elements = Vec::new();
        let mut elem_id = 0;

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let n000 = i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
                    let n100 = n000 + 1;
                    let n010 = n000 + (nx + 1);
                    let n110 = n010 + 1;
                    let n001 = n000 + (nx + 1) * (ny + 1);
                    let n101 = n001 + 1;
                    let n011 = n001 + (nx + 1);
                    let n111 = n011 + 1;

                    // 5 tetrahedra per cube
                    elements.push(Element3D {
                        id: elem_id,
                        nodes: [n000, n100, n110, n111],
                    });
                    elem_id += 1;

                    elements.push(Element3D {
                        id: elem_id,
                        nodes: [n000, n110, n010, n111],
                    });
                    elem_id += 1;

                    elements.push(Element3D {
                        id: elem_id,
                        nodes: [n000, n010, n011, n111],
                    });
                    elem_id += 1;

                    elements.push(Element3D {
                        id: elem_id,
                        nodes: [n000, n001, n101, n111],
                    });
                    elem_id += 1;

                    elements.push(Element3D {
                        id: elem_id,
                        nodes: [n000, n100, n101, n111],
                    });
                    elem_id += 1;
                }
            }
        }

        Ok(Self { nodes, elements, h })
    }

    /// Get number of nodes
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Get number of elements
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mesh2d_rectangle() {
        let mesh = Mesh2D::rectangle(100e-9, 50e-9, 10e-9).unwrap();

        assert!(mesh.n_nodes() > 0);
        assert!(mesh.n_elements() > 0);
        assert_eq!(mesh.h, 10e-9);
    }

    #[test]
    fn test_mesh3d_cuboid() {
        let mesh = Mesh3D::cuboid(100e-9, 50e-9, 10e-9, 10e-9).unwrap();

        assert!(mesh.n_nodes() > 0);
        assert!(mesh.n_elements() > 0);
    }

    #[test]
    fn test_element_area() {
        let mesh = Mesh2D::rectangle(10.0, 10.0, 2.0).unwrap();

        if !mesh.elements.is_empty() {
            let area = mesh.element_area(0);
            assert!(area > 0.0);
        }
    }
}
