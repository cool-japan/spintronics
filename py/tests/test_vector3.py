"""Tests for the Vector3 3D vector binding (src/python/vector.rs)."""
import pytest

import spintronics as sp


class TestConstruction:
    def test_constructor_sets_components(self):
        v = sp.Vector3(1.0, 2.0, 3.0)
        assert v.x == 1.0
        assert v.y == 2.0
        assert v.z == 3.0

    def test_setters_mutate_components(self):
        v = sp.Vector3(0.0, 0.0, 0.0)
        v.x = 1.5
        v.y = -2.5
        v.z = 3.5
        assert v.to_tuple() == (1.5, -2.5, 3.5)

    def test_repr(self):
        v = sp.Vector3(1.0, 0.0, 0.0)
        assert repr(v) == "Vector3(1.000000, 0.000000, 0.000000)"

    def test_repr_formats_six_decimal_places(self):
        v = sp.Vector3(1.23456789, -2.3, 0.0)
        assert repr(v) == "Vector3(1.234568, -2.300000, 0.000000)"


class TestStaticFactories:
    def test_unit_x(self):
        assert sp.Vector3.unit_x().to_tuple() == (1.0, 0.0, 0.0)

    def test_unit_y(self):
        assert sp.Vector3.unit_y().to_tuple() == (0.0, 1.0, 0.0)

    def test_unit_z(self):
        assert sp.Vector3.unit_z().to_tuple() == (0.0, 0.0, 1.0)

    def test_zero(self):
        assert sp.Vector3.zero().to_tuple() == (0.0, 0.0, 0.0)


class TestDotProduct:
    def test_orthogonal_vectors(self):
        assert sp.Vector3(1, 0, 0).dot(sp.Vector3(0, 1, 0)) == 0.0

    def test_parallel_unit_vectors(self):
        assert sp.Vector3(1, 0, 0).dot(sp.Vector3(1, 0, 0)) == 1.0

    def test_general_vectors(self):
        assert sp.Vector3(1, 2, 3).dot(sp.Vector3(4, 5, 6)) == 32.0

    def test_self_dot_equals_magnitude_squared(self):
        v = sp.Vector3(1, 2, 2)
        assert v.dot(v) == pytest.approx(v.magnitude() ** 2)


class TestCrossProduct:
    def test_x_cross_y_is_z(self):
        result = sp.Vector3.unit_x().cross(sp.Vector3.unit_y())
        assert result.to_tuple() == (0.0, 0.0, 1.0)

    def test_y_cross_z_is_x(self):
        result = sp.Vector3.unit_y().cross(sp.Vector3.unit_z())
        assert result.to_tuple() == (1.0, 0.0, 0.0)

    def test_z_cross_x_is_y(self):
        result = sp.Vector3.unit_z().cross(sp.Vector3.unit_x())
        assert result.to_tuple() == (0.0, 1.0, 0.0)

    def test_anticommutativity(self):
        a = sp.Vector3(1, 2, 3)
        b = sp.Vector3(4, 5, 6)
        ab = a.cross(b)
        ba = b.cross(a)
        negated_ba = tuple(-c for c in ba.to_tuple())
        assert ab.to_tuple() == pytest.approx(negated_ba)

    def test_parallel_vectors_give_zero(self):
        a = sp.Vector3(2, 4, 6)
        b = sp.Vector3(1, 2, 3)
        result = a.cross(b)
        assert result.magnitude() == pytest.approx(0.0, abs=1e-12)

    def test_result_perpendicular_to_both_operands(self):
        a = sp.Vector3(1, 0, 0)
        b = sp.Vector3(0, 1, 1)
        result = a.cross(b)
        assert result.dot(a) == pytest.approx(0.0, abs=1e-12)
        assert result.dot(b) == pytest.approx(0.0, abs=1e-12)


class TestMagnitude:
    def test_three_four_five_triangle(self):
        assert sp.Vector3(3, 4, 0).magnitude() == pytest.approx(5.0)

    def test_zero_vector(self):
        assert sp.Vector3.zero().magnitude() == 0.0

    def test_negative_components(self):
        assert sp.Vector3(-3, -4, 0).magnitude() == pytest.approx(5.0)

    def test_unit_vectors_have_unit_magnitude(self):
        assert sp.Vector3.unit_x().magnitude() == 1.0
        assert sp.Vector3.unit_y().magnitude() == 1.0
        assert sp.Vector3.unit_z().magnitude() == 1.0


class TestNormalize:
    def test_normalize_returns_unit_vector(self):
        n = sp.Vector3(3, 4, 0).normalize()
        assert n.to_tuple() == pytest.approx((0.6, 0.8, 0.0))
        assert n.magnitude() == pytest.approx(1.0)

    def test_normalize_preserves_direction(self):
        v = sp.Vector3(1, 2, 2)
        n = v.normalize()
        # Same direction => cross product with the original is ~zero.
        assert n.cross(v).magnitude() == pytest.approx(0.0, abs=1e-10)

    def test_normalize_zero_vector_stays_zero(self):
        # The Rust implementation guards against division by zero and
        # returns the original (zero) vector unchanged.
        n = sp.Vector3.zero().normalize()
        assert n.to_tuple() == (0.0, 0.0, 0.0)


class TestArithmetic:
    def test_add(self):
        result = sp.Vector3(1, 2, 3) + sp.Vector3(4, 5, 6)
        assert result.to_tuple() == (5.0, 7.0, 9.0)

    def test_sub(self):
        result = sp.Vector3(4, 5, 6) - sp.Vector3(1, 2, 3)
        assert result.to_tuple() == (3.0, 3.0, 3.0)

    def test_scalar_mul(self):
        result = sp.Vector3(1, 2, 3) * 2.0
        assert result.to_tuple() == (2.0, 4.0, 6.0)

    def test_scalar_rmul(self):
        result = 2.0 * sp.Vector3(1, 2, 3)
        assert result.to_tuple() == (2.0, 4.0, 6.0)

    def test_mul_and_rmul_agree(self):
        v = sp.Vector3(1.5, -2.5, 3.5)
        assert (v * 3.0).to_tuple() == (3.0 * v).to_tuple()


class TestConversions:
    def test_to_tuple(self):
        assert sp.Vector3(1.0, 2.0, 3.0).to_tuple() == (1.0, 2.0, 3.0)

    def test_to_list(self):
        assert sp.Vector3(1.0, 2.0, 3.0).to_list() == [1.0, 2.0, 3.0]

    def test_to_tuple_and_to_list_agree(self):
        v = sp.Vector3(1.0, 2.0, 3.0)
        assert list(v.to_tuple()) == v.to_list()


class TestIdentitySemantics:
    def test_equality_is_by_identity_not_value(self):
        # Vector3 does not implement value equality (no __richcmp__ defined
        # in src/python/vector.rs), so two distinct instances with
        # identical components are NOT `==`. Use `.to_tuple()` (or
        # component-wise) comparisons for value equality instead.
        a = sp.Vector3(1, 2, 3)
        b = sp.Vector3(1, 2, 3)
        assert not (a == b)
        assert a.to_tuple() == b.to_tuple()
        assert a == a
