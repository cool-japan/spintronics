//! Reverse-mode automatic differentiation tape.
//!
//! This module implements a Wengert-list (tape) based reverse-mode AD engine.
//! Every scalar operation records itself on the tape; the backward pass then
//! propagates gradients from outputs back to inputs in reverse order
//! (back-propagation).
//!
//! # Design
//!
//! Each elementary operation is stored as a [`TapeOp`]:
//! - `output_idx` — index of the result node in `values`
//! - `inputs[0..n_inputs]` — indices of the operand nodes
//! - `local_grads[0..n_inputs]` — ∂output/∂input_i (local Jacobian entries)
//!
//! The backward pass accumulates: `grads[input_i] += local_grads[i] * grads[output]`.
//!
//! # References
//!
//! - Griewank & Walther, *Evaluating Derivatives*, SIAM 2008
//! - Baydin et al., "Automatic Differentiation in Machine Learning: a Survey",
//!   *JMLR* **18**, 1 (2018)
//! - Linnainmaa, "Taylor expansion of the accumulated rounding error",
//!   *BIT* **16**, 146 (1976)

use std::cell::RefCell;
use std::ops::{Add, Div, Mul, Neg, Sub};

// ─── Internal tape operation record ─────────────────────────────────────────

/// A single elementary operation recorded on the tape.
///
/// For binary ops: `n_inputs == 2`, both slots of `inputs`/`local_grads` used.
/// For unary ops:  `n_inputs == 1`, only slot 0 used.
#[derive(Clone, Debug)]
struct TapeOp {
    /// Index of the output node in `values`/`grads`.
    output_idx: usize,
    /// Indices of the input nodes.
    inputs: [usize; 2],
    /// Number of valid inputs (1 or 2).
    n_inputs: usize,
    /// Local partial derivatives ∂output/∂inputs[i].
    local_grads: [f64; 2],
}

// ─── Tape ────────────────────────────────────────────────────────────────────

/// A reverse-mode automatic differentiation tape.
///
/// The tape records every operation performed on [`Var`] nodes attached to it.
/// After the forward pass, call [`Tape::backward`] to propagate gradients from
/// a scalar loss variable back through all recorded operations.
///
/// Each iteration of an optimisation loop should use a **fresh** `Tape`:
/// gradients accumulate across `backward` calls on the same tape.
pub struct Tape {
    ops: RefCell<Vec<TapeOp>>,
    values: RefCell<Vec<f64>>,
    grads: RefCell<Vec<f64>>,
}

impl Default for Tape {
    fn default() -> Self {
        Self::new()
    }
}

impl Tape {
    /// Create a new, empty tape.
    pub fn new() -> Self {
        Self {
            ops: RefCell::new(Vec::new()),
            values: RefCell::new(Vec::new()),
            grads: RefCell::new(Vec::new()),
        }
    }

    /// Allocate a leaf node (no op recorded).  Returns the new node index.
    fn push_leaf(&self, value: f64) -> usize {
        let mut vals = self.values.borrow_mut();
        let mut grds = self.grads.borrow_mut();
        let idx = vals.len();
        vals.push(value);
        grds.push(0.0);
        idx
    }

    /// Allocate an output node and record the operation on the tape.
    fn push_op(
        &self,
        value: f64,
        inputs: [usize; 2],
        n_inputs: usize,
        local_grads: [f64; 2],
    ) -> usize {
        let idx = self.push_leaf(value);
        self.ops.borrow_mut().push(TapeOp {
            output_idx: idx,
            inputs,
            n_inputs,
            local_grads,
        });
        idx
    }

    // ── Backward pass ────────────────────────────────────────────────────────

    /// Run the backward pass starting from `loss`.
    ///
    /// Sets `grads[loss.idx] = 1.0` and propagates via the chain rule in reverse
    /// topological order through all recorded operations.
    pub fn backward(&self, loss: Var<'_>) {
        // Seed the loss gradient.
        {
            let mut grds = self.grads.borrow_mut();
            grds[loss.idx] = 1.0;
        }
        // Iterate ops in reverse insertion order.
        let ops = self.ops.borrow();
        let mut grds = self.grads.borrow_mut();
        for op in ops.iter().rev() {
            let out_grad = grds[op.output_idx];
            for i in 0..op.n_inputs {
                grds[op.inputs[i]] += op.local_grads[i] * out_grad;
            }
        }
    }

    /// Clear the tape completely (ops, values, and grads).
    pub fn reset(&self) {
        self.ops.borrow_mut().clear();
        self.values.borrow_mut().clear();
        self.grads.borrow_mut().clear();
    }

    /// Returns the number of values (nodes) currently on the tape.
    pub fn values_len(&self) -> usize {
        self.values.borrow().len()
    }

    /// Returns the number of operations currently recorded.
    pub fn ops_len(&self) -> usize {
        self.ops.borrow().len()
    }

    // ── Internal helpers ─────────────────────────────────────────────────────

    fn value_at(&self, idx: usize) -> f64 {
        self.values.borrow()[idx]
    }

    fn grad_at(&self, idx: usize) -> f64 {
        self.grads.borrow()[idx]
    }
}

// ─── Var ─────────────────────────────────────────────────────────────────────

/// A scalar variable tracked on a [`Tape`].
///
/// `Var` is a lightweight handle: it stores only a reference to the tape and
/// the index of its node.  Arithmetic and transcendental operations on `Var`
/// automatically record the corresponding operation on the shared tape.
///
/// # Lifetime
///
/// The `'t` lifetime ties every `Var` to its originating tape, preventing a
/// variable from escaping its tape's scope.
#[derive(Copy, Clone)]
pub struct Var<'t> {
    tape: &'t Tape,
    /// Node index in the tape's value/grad vectors.
    pub(crate) idx: usize,
}

impl<'t> Var<'t> {
    // ── Construction ─────────────────────────────────────────────────────────

    /// Create a new leaf (input) variable on `tape` with the given `value`.
    ///
    /// Leaf variables have no parents; their gradient after `backward` gives
    /// the partial derivative of the loss with respect to this variable.
    pub fn leaf(tape: &'t Tape, value: f64) -> Self {
        let idx = tape.push_leaf(value);
        Self { tape, idx }
    }

    // ── Accessors ─────────────────────────────────────────────────────────────

    /// Return the forward-pass value of this variable.
    pub fn value(&self) -> f64 {
        self.tape.value_at(self.idx)
    }

    /// Return the accumulated gradient of this variable after `tape.backward`.
    ///
    /// Returns `0.0` before `backward` is called.
    pub fn grad(&self) -> f64 {
        self.tape.grad_at(self.idx)
    }

    // ── Unary differentiable operations ──────────────────────────────────────

    /// Differentiable `sin`: returns `sin(self)`, ∂/∂self = cos(self).
    pub fn sin(self) -> Var<'t> {
        let v = self.value();
        let out = self.tape.push_op(v.sin(), [self.idx, 0], 1, [v.cos(), 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable `cos`: returns `cos(self)`, ∂/∂self = −sin(self).
    pub fn cos(self) -> Var<'t> {
        let v = self.value();
        let out = self
            .tape
            .push_op(v.cos(), [self.idx, 0], 1, [-v.sin(), 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable `exp`: returns `exp(self)`, ∂/∂self = exp(self).
    pub fn exp(self) -> Var<'t> {
        let v = self.value();
        let ev = v.exp();
        let out = self.tape.push_op(ev, [self.idx, 0], 1, [ev, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable natural logarithm: returns `ln(self)`, ∂/∂self = 1/self.
    ///
    /// Panics (numerically) if `self ≤ 0`.
    pub fn ln(self) -> Var<'t> {
        let v = self.value();
        let out = self.tape.push_op(v.ln(), [self.idx, 0], 1, [1.0 / v, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable square root: returns `√self`, ∂/∂self = 1/(2√self).
    pub fn sqrt(self) -> Var<'t> {
        let v = self.value();
        let sv = v.sqrt();
        let grad = if sv.abs() < f64::EPSILON {
            0.0
        } else {
            0.5 / sv
        };
        let out = self.tape.push_op(sv, [self.idx, 0], 1, [grad, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable `tanh`: returns `tanh(self)`, ∂/∂self = 1 − tanh²(self).
    pub fn tanh(self) -> Var<'t> {
        let v = self.value();
        let tv = v.tanh();
        let grad = 1.0 - tv * tv;
        let out = self.tape.push_op(tv, [self.idx, 0], 1, [grad, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable integer power: returns `self^n`, ∂/∂self = n·self^(n−1).
    pub fn powi(self, n: i32) -> Var<'t> {
        let v = self.value();
        let grad = (n as f64) * v.powi(n - 1);
        let out = self.tape.push_op(v.powi(n), [self.idx, 0], 1, [grad, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable real power: returns `self^p`, ∂/∂self = p·self^(p−1).
    pub fn powf(self, p: f64) -> Var<'t> {
        let v = self.value();
        let grad = p * v.powf(p - 1.0);
        let out = self.tape.push_op(v.powf(p), [self.idx, 0], 1, [grad, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable absolute value: returns `|self|`, ∂/∂self = sign(self).
    pub fn abs(self) -> Var<'t> {
        let v = self.value();
        let grad = if v > 0.0 {
            1.0
        } else if v < 0.0 {
            -1.0
        } else {
            0.0
        };
        let out = self.tape.push_op(v.abs(), [self.idx, 0], 1, [grad, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }

    /// Differentiable reciprocal: returns `1/self`, ∂/∂self = −1/self².
    pub fn recip(self) -> Var<'t> {
        let v = self.value();
        let grad = -1.0 / (v * v);
        let out = self.tape.push_op(1.0 / v, [self.idx, 0], 1, [grad, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

// ─── Arithmetic trait implementations ────────────────────────────────────────

// --- Add ---

impl<'t> Add<Var<'t>> for Var<'t> {
    type Output = Var<'t>;
    fn add(self, rhs: Var<'t>) -> Var<'t> {
        let val = self.value() + rhs.value();
        let out = self.tape.push_op(val, [self.idx, rhs.idx], 2, [1.0, 1.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

impl<'t> Add<f64> for Var<'t> {
    type Output = Var<'t>;
    fn add(self, rhs: f64) -> Var<'t> {
        let val = self.value() + rhs;
        let out = self.tape.push_op(val, [self.idx, 0], 1, [1.0, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

impl<'t> Add<Var<'t>> for f64 {
    type Output = Var<'t>;
    fn add(self, rhs: Var<'t>) -> Var<'t> {
        rhs + self
    }
}

// --- Sub ---

impl<'t> Sub<Var<'t>> for Var<'t> {
    type Output = Var<'t>;
    fn sub(self, rhs: Var<'t>) -> Var<'t> {
        let val = self.value() - rhs.value();
        let out = self.tape.push_op(val, [self.idx, rhs.idx], 2, [1.0, -1.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

impl<'t> Sub<f64> for Var<'t> {
    type Output = Var<'t>;
    fn sub(self, rhs: f64) -> Var<'t> {
        let val = self.value() - rhs;
        let out = self.tape.push_op(val, [self.idx, 0], 1, [1.0, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

impl<'t> Sub<Var<'t>> for f64 {
    type Output = Var<'t>;
    fn sub(self, rhs: Var<'t>) -> Var<'t> {
        let val = self - rhs.value();
        let out = rhs.tape.push_op(val, [rhs.idx, 0], 1, [-1.0, 0.0]);
        Var {
            tape: rhs.tape,
            idx: out,
        }
    }
}

// --- Mul ---

impl<'t> Mul<Var<'t>> for Var<'t> {
    type Output = Var<'t>;
    fn mul(self, rhs: Var<'t>) -> Var<'t> {
        let (lv, rv) = (self.value(), rhs.value());
        let out = self.tape.push_op(lv * rv, [self.idx, rhs.idx], 2, [rv, lv]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

impl<'t> Mul<f64> for Var<'t> {
    type Output = Var<'t>;
    fn mul(self, rhs: f64) -> Var<'t> {
        let val = self.value() * rhs;
        let out = self.tape.push_op(val, [self.idx, 0], 1, [rhs, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

impl<'t> Mul<Var<'t>> for f64 {
    type Output = Var<'t>;
    fn mul(self, rhs: Var<'t>) -> Var<'t> {
        rhs * self
    }
}

// --- Div ---

impl<'t> Div<Var<'t>> for Var<'t> {
    type Output = Var<'t>;
    fn div(self, rhs: Var<'t>) -> Var<'t> {
        let (lv, rv) = (self.value(), rhs.value());
        let g_rhs = -lv / (rv * rv);
        let out = self
            .tape
            .push_op(lv / rv, [self.idx, rhs.idx], 2, [1.0 / rv, g_rhs]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

impl<'t> Div<f64> for Var<'t> {
    type Output = Var<'t>;
    fn div(self, rhs: f64) -> Var<'t> {
        let val = self.value() / rhs;
        let out = self.tape.push_op(val, [self.idx, 0], 1, [1.0 / rhs, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

impl<'t> Div<Var<'t>> for f64 {
    type Output = Var<'t>;
    fn div(self, rhs: Var<'t>) -> Var<'t> {
        let rv = rhs.value();
        let val = self / rv;
        let grad = -self / (rv * rv);
        let out = rhs.tape.push_op(val, [rhs.idx, 0], 1, [grad, 0.0]);
        Var {
            tape: rhs.tape,
            idx: out,
        }
    }
}

// --- Neg ---

impl<'t> Neg for Var<'t> {
    type Output = Var<'t>;
    fn neg(self) -> Var<'t> {
        let val = -self.value();
        let out = self.tape.push_op(val, [self.idx, 0], 1, [-1.0, 0.0]);
        Var {
            tape: self.tape,
            idx: out,
        }
    }
}

// ─── Free-function utilities ──────────────────────────────────────────────────

/// Compute a numerical gradient using the central finite-difference formula.
///
/// `(f(x + h) − f(x − h)) / (2h)`
///
/// Typical step size: `h = 1e-5` to `1e-7`.
pub fn finite_diff_grad<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    (f(x + h) - f(x - h)) / (2.0 * h)
}

/// Compare the AD gradient of `f` at `x` with a central finite-difference
/// approximation, returning `true` if they agree within `tol`.
///
/// The closure `f` receives a fresh tape and the single leaf variable `x`.
pub fn check_gradient<F>(f: F, x: f64, tol: f64) -> bool
where
    F: for<'a> Fn(&'a Tape, Var<'a>) -> Var<'a>,
{
    // AD gradient
    let tape = Tape::new();
    let xv = Var::leaf(&tape, x);
    let out = f(&tape, xv);
    tape.backward(out);
    let ad_grad = xv.grad();

    // Finite-difference gradient
    let fd_grad = finite_diff_grad(
        |xi| {
            let t2 = Tape::new();
            let xv2 = Var::leaf(&t2, xi);
            f(&t2, xv2).value()
        },
        x,
        1e-6,
    );

    (ad_grad - fd_grad).abs() < tol
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    // Helper: run backward on the given computation and return (value, grad of x)
    fn run_unary<F: Fn(Var<'_>) -> Var<'_>>(x_val: f64, f: F) -> (f64, f64) {
        let tape = Tape::new();
        let x = Var::leaf(&tape, x_val);
        let out = f(x);
        tape.backward(out);
        (out.value(), x.grad())
    }

    // 1. Leaf node: value correct, grad zero before backward
    #[test]
    fn test_leaf_node() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 7.0);
        assert!((x.value() - 7.0).abs() < 1e-15);
        assert!((x.grad() - 0.0).abs() < 1e-15);
    }

    // 2. Addition: both grads = 1
    #[test]
    fn test_add_grads() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 3.0);
        let y = Var::leaf(&tape, 5.0);
        let z = x + y;
        tape.backward(z);
        assert!((z.value() - 8.0).abs() < 1e-15);
        assert!((x.grad() - 1.0).abs() < 1e-15);
        assert!((y.grad() - 1.0).abs() < 1e-15);
    }

    // 3. Multiplication: grads = [y, x]
    #[test]
    fn test_mul_grads() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 3.0);
        let y = Var::leaf(&tape, 4.0);
        let z = x * y;
        tape.backward(z);
        assert!((z.value() - 12.0).abs() < 1e-15);
        assert!((x.grad() - 4.0).abs() < 1e-15); // ∂z/∂x = y = 4
        assert!((y.grad() - 3.0).abs() < 1e-15); // ∂z/∂y = x = 3
    }

    // 4. x * x: grad = 2x
    #[test]
    fn test_square_grad() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 5.0);
        let z = x * x;
        tape.backward(z);
        assert!((z.value() - 25.0).abs() < 1e-15);
        assert!((x.grad() - 10.0).abs() < 1e-14); // 2 * 5 = 10
    }

    // 5. sin: grad = cos(x)
    #[test]
    fn test_sin_grad() {
        let x_val = PI / 4.0;
        let (val, grad) = run_unary(x_val, |x| x.sin());
        assert!((val - x_val.sin()).abs() < 1e-15);
        assert!((grad - x_val.cos()).abs() < 1e-15);
    }

    // 6. cos: grad = -sin(x)
    #[test]
    fn test_cos_grad() {
        let x_val = PI / 3.0;
        let (val, grad) = run_unary(x_val, |x| x.cos());
        assert!((val - x_val.cos()).abs() < 1e-15);
        assert!((grad - (-x_val.sin())).abs() < 1e-15);
    }

    // 7. exp: grad = exp(x)
    #[test]
    fn test_exp_grad() {
        let x_val = 1.5;
        let (val, grad) = run_unary(x_val, |x| x.exp());
        assert!((val - x_val.exp()).abs() < 1e-14);
        assert!((grad - x_val.exp()).abs() < 1e-14);
    }

    // 8. ln: grad = 1/x
    #[test]
    fn test_ln_grad() {
        let x_val = 2.0;
        let (val, grad) = run_unary(x_val, |x| x.ln());
        assert!((val - x_val.ln()).abs() < 1e-15);
        assert!((grad - 1.0 / x_val).abs() < 1e-15);
    }

    // 9. sqrt: grad = 1/(2*sqrt(x))
    #[test]
    fn test_sqrt_grad() {
        let x_val = 9.0;
        let (val, grad) = run_unary(x_val, |x| x.sqrt());
        assert!((val - 3.0).abs() < 1e-14);
        assert!((grad - 1.0 / 6.0).abs() < 1e-14);
    }

    // 10. tanh: grad = 1 - tanh²(x)
    #[test]
    fn test_tanh_grad() {
        let x_val = 0.5;
        let (val, grad) = run_unary(x_val, |x| x.tanh());
        let tv = x_val.tanh();
        assert!((val - tv).abs() < 1e-15);
        assert!((grad - (1.0 - tv * tv)).abs() < 1e-15);
    }

    // 11. Chain rule: sin(x²) — grad = 2x·cos(x²)
    #[test]
    fn test_chain_rule_sin_square() {
        let x_val = 1.2;
        let (val, grad) = run_unary(x_val, |x| (x * x).sin());
        let expected_val = (x_val * x_val).sin();
        let expected_grad = 2.0 * x_val * (x_val * x_val).cos();
        assert!((val - expected_val).abs() < 1e-14);
        assert!((grad - expected_grad).abs() < 1e-13);
    }

    // 12. Product of three variables: z = (x*y)*w
    #[test]
    fn test_triple_product() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 2.0);
        let y = Var::leaf(&tape, 3.0);
        let w = Var::leaf(&tape, 4.0);
        let z = (x * y) * w;
        tape.backward(z);
        // dz/dx = y*w = 12, dz/dy = x*w = 8, dz/dw = x*y = 6
        assert!((z.value() - 24.0).abs() < 1e-14);
        assert!((x.grad() - 12.0).abs() < 1e-13);
        assert!((y.grad() - 8.0).abs() < 1e-13);
        assert!((w.grad() - 6.0).abs() < 1e-13);
    }

    // 13. powi(x, 3): grad = 3x²
    #[test]
    fn test_powi_grad() {
        let x_val = 2.0;
        let (val, grad) = run_unary(x_val, |x| x.powi(3));
        assert!((val - 8.0).abs() < 1e-14);
        assert!((grad - 12.0).abs() < 1e-14); // 3 * 4 = 12
    }

    // 14. Gradient check: finite diff vs AD on a complex function
    #[test]
    fn test_gradient_check_complex() {
        // f(x) = sin(x²) * exp(-x) + ln(x)
        let x_val = 1.5;
        let ok = check_gradient(
            |_tape, x| {
                let x_sq = x * x;
                let s = x_sq.sin();
                let e = (-x).exp();
                let l = x.ln();
                s * e + l
            },
            x_val,
            1e-8,
        );
        assert!(ok, "AD gradient did not match finite difference");
    }

    // 15. f64 scalar arithmetic: constant does not affect grad path
    #[test]
    fn test_scalar_mul_add() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 3.0);
        let z = x * 2.0 + 5.0; // dz/dx = 2
        tape.backward(z);
        assert!((z.value() - 11.0).abs() < 1e-14);
        assert!((x.grad() - 2.0).abs() < 1e-14);
    }

    // 16. Subtraction: grads = [1, -1]
    #[test]
    fn test_sub_grads() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 7.0);
        let y = Var::leaf(&tape, 3.0);
        let z = x - y;
        tape.backward(z);
        assert!((x.grad() - 1.0).abs() < 1e-15);
        assert!((y.grad() - (-1.0)).abs() < 1e-15);
    }

    // 17. Division: grads = [1/y, -x/y²]
    #[test]
    fn test_div_grads() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 6.0);
        let y = Var::leaf(&tape, 3.0);
        let z = x / y;
        tape.backward(z);
        assert!((z.value() - 2.0).abs() < 1e-15);
        assert!((x.grad() - 1.0 / 3.0).abs() < 1e-14); // 1/y
        assert!((y.grad() - (-6.0 / 9.0)).abs() < 1e-14); // -x/y²
    }

    // 18. Negation: grad = -1
    #[test]
    fn test_neg_grad() {
        let (_, grad) = run_unary(4.0, |x| -x);
        assert!((grad - (-1.0)).abs() < 1e-15);
    }

    // 19. Tape state: values_len and ops_len
    #[test]
    fn test_tape_lengths() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 1.0);
        let y = Var::leaf(&tape, 2.0);
        assert_eq!(tape.values_len(), 2);
        assert_eq!(tape.ops_len(), 0);
        let _z = x + y;
        assert_eq!(tape.values_len(), 3);
        assert_eq!(tape.ops_len(), 1);
    }

    // 20. Reset clears tape
    #[test]
    fn test_tape_reset() {
        let tape = Tape::new();
        let x = Var::leaf(&tape, 1.0);
        let y = Var::leaf(&tape, 2.0);
        let _z = x + y;
        tape.reset();
        assert_eq!(tape.values_len(), 0);
        assert_eq!(tape.ops_len(), 0);
    }
}
