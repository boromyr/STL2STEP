"""
test_geometry.py - seconds-long checks of refit.py's pure geometry the batteries only see as "arcs lost".

    python tools/test_geometry.py
"""
import os, sys, math, json, tempfile, unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import refit as R


def cylinder(r, axis=(0.0, 0.0, 1.0), center=(0.0, 0.0, 0.0)):
    a = np.asarray(axis, float)
    return R.Prim(R.AXIAL, np.asarray(center, float), a / np.linalg.norm(a), float(r), 0.0)


class MeetSpline(unittest.TestCase):
    def test_cylinder_and_tilted_plane(self):
        # the meet is an ellipse: (10 cos t, 10 sin t, 0.5 x + 3)
        cyl = cylinder(10.0)
        n = np.array([-0.5, 0.0, 1.0])
        pl = R.Prim(R.PLANE, np.array([0.0, 0.0, 3.0]), n / np.linalg.norm(n))
        e = lambda t: np.array([10 * math.cos(t), 10 * math.sin(t), 5 * math.cos(t) + 3])
        cv = R._meet_spline(cyl, pl, e(0.2), e(1.0), 1e-3)
        self.assertIsNotNone(cv)
        c = cv.c
        Q = cv.point(np.linspace(float(c.FirstParameter()), float(c.LastParameter()), 41))
        self.assertLess(float(np.abs(cyl.dist(Q)).max()), 1e-3)
        self.assertLess(float(np.abs(pl.dist(Q)).max()), 1e-3)

    def test_end_off_the_meet_refused(self):
        # an end 0.5 mm off the cylinder: walked onto both it moves past tol - not a meet between these ends
        cyl = cylinder(10.0)
        pl = R.Prim(R.PLANE, np.array([0.0, 0.0, 3.0]), np.array([0.0, 0.0, 1.0]))
        P0, P1 = np.array([10.0, 0.0, 3.0]), np.array([0.0, 10.5, 3.0])
        self.assertIsNone(R._meet_spline(cyl, pl, P0, P1, 1e-3))


class NearCoaxCircle(unittest.TestCase):
    def test_axes_apart_by_a_hair(self):
        # a cylinder and a cone r 271 on axes 0.012 mm apart: the circle about the axis between, through both ends
        cyl = cylinder(271.0)
        cone = R.Prim(R.AXIAL, np.array([0.012, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]), 271.0, 1.0)
        P0 = np.array([271 * math.cos(0.1), 271 * math.sin(0.1), 0.0])
        P1 = np.array([271 * math.cos(0.16), 271 * math.sin(0.16), 0.0])
        c = R._near_coax_circle(cyl, cone, P0, P1, 1e-3)
        self.assertIsNotNone(c)
        self.assertAlmostEqual(c.r, 271.0, delta=0.02)
        self.assertLess(float(c.dist(np.vstack([P0, P1])).max()), 1e-9)

    def test_crossed_axes_refused(self):
        tilt = np.array([0.0, math.sin(math.radians(5)), math.cos(math.radians(5))])
        P0, P1 = np.array([10.0, 0, 0]), np.array([0.0, 10.0, 0])
        self.assertIsNone(R._near_coax_circle(cylinder(10.0), cylinder(10.0, axis=tilt), P0, P1, 1e-3))


class LoadHints(unittest.TestCase):
    def test_contours_fields(self):
        doc = {"refit_hints": 1, "pass2": [1, 2], "groups": [
            {"type": "auto", "points": [[0, 0, 0], [1, 0, 0]], "areas": [1, 1], "max_dev": 2.0,
             "surface": {"kind": "axial", "center": [0, 0, 0], "axis": [0, 0, 2], "r0": 5.0, "slope": 0.0}},
            {"type": "plane", "points": [[0, 0, 0]], "surface": {"kind": "axial", "axis": [0, 0, 1]}}]}
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as fh:
            json.dump(doc, fh)
        try:
            h = R.load_hints(fh.name)
        finally:
            os.remove(fh.name)
        self.assertEqual(len(h), 2)
        self.assertEqual(h[0].max_dev, 2.0)
        self.assertAlmostEqual(float(np.linalg.norm(h[0].surface.axis)), 1.0)
        self.assertIsNone(h[1].max_dev)
        self.assertIsNone(h[1].surface)     # (a broken surface: dropped, the group kept)


if __name__ == "__main__":
    unittest.main()
