import math
import itertools
import Rhino.Geometry as geo


class TriPanel:
    def __init__(
        self,
        pt_base: geo.Point3d,
        pt_foot: geo.Point3d,
        pt_vertical: geo.Point3d,
        perp_vec: geo.Vector3d,
    ):
        self.pt_base = pt_base
        self.pt_foot = pt_foot
        self.pt_vertical = pt_vertical
        self.foot_length = pt_vertical.DistanceTo(pt_foot)

        self.perp_vec = perp_vec

    @property
    def crv(self) -> geo.Curve:
        return geo.PolylineCurve(
            [self.pt_base, self.pt_foot, self.pt_vertical, self.pt_base]
        )

    def move(self, ratio: float):
        # type: (...) -> TriPanel
        if ratio <= 0:
            return self

        height = self.foot_length * ratio
        pt_vertical = self.pt_vertical + self.perp_vec * height  # type: ignore

        moved_pts = []  # type: list[geo.Point3d]
        for pt in [self.pt_base, self.pt_foot]:
            vec = pt - self.pt_vertical  # type: ignore
            vec.Unitize()

            seg_length = pt.DistanceTo(self.pt_vertical)
            if seg_length > height:
                move_length = math.sqrt(seg_length**2 - height**2)
            else:
                move_length = 0.0
            moved_pts.append(self.pt_vertical + vec * move_length)

        return TriPanel(moved_pts[0], moved_pts[1], pt_vertical, self.perp_vec)


class TriModule:
    def __init__(self, pt_a: geo.Point3d, pt_b: geo.Point3d, pt_c: geo.Point3d):
        self.pt_a = pt_a
        self.pt_b = pt_b
        self.pt_c = pt_c

        self.incenter = self.get_incenter()
        self.perp_vec = self.get_perp_vec()

        self.base_panels = self.get_tri_panels()

    @property
    def crv(self) -> geo.Curve:
        return geo.PolylineCurve([self.pt_a, self.pt_b, self.pt_c, self.pt_a])

    @property
    def crvs(self) -> list[geo.Curve]:
        return [panel.crv for panel in self.base_panels]

    def get_incenter(self):
        # type: (...) -> geo.Point3d
        pt_a, pt_b, pt_c = self.pt_a, self.pt_b, self.pt_c
        ab = pt_a.DistanceTo(pt_b)
        bc = pt_b.DistanceTo(pt_c)
        ca = pt_c.DistanceTo(pt_a)
        incenter = geo.Point3d((pt_a * bc + pt_b * ca + pt_c * ab) / (ab + bc + ca))  # type: ignore
        return incenter

    def get_perp_vec(self):
        # type: (...) -> geo.Vector3d
        vec_ab = self.pt_b - self.pt_a  # type: ignore
        vec_bc = self.pt_c - self.pt_b  # type: ignore
        vec_perp = geo.Vector3d.CrossProduct(vec_ab, vec_bc)
        vec_perp.Unitize()
        return vec_perp

    def get_tri_panels(self):
        # type: (...) -> list[TriPanel]
        pts = [self.pt_a, self.pt_b, self.pt_c]
        incenter = self.incenter

        tri_panels = []  # type: list[TriPanel]
        for bases in itertools.combinations(pts, 2):
            line = geo.Line(*bases)  # type: ignore
            foot = line.ClosestPoint(incenter, True)
            if not foot:
                raise ValueError("Impossible Triangle")
            tri_panels += [
                TriPanel(base, foot, incenter, self.perp_vec) for base in bases
            ]

        return tri_panels

    def crvs_at_ratio(self, ratio: float):
        # type: (...) -> list[geo.Curve]
        return [panel.move(ratio).crv for panel in self.base_panels]


class TriFacade:
    def __init__(self, tri_modules: list[TriModule]):
        self.tri_modules = tri_modules

    def get_breps(self, sunlight_vec: geo.Vector3d):
        # type: (...) -> list[geo.Brep]
        crvs = []  # type: list[geo.Curve]
        standard_vec = -sunlight_vec  # type: ignore
        for tri_module in self.tri_modules:
            cos_val = tri_module.perp_vec * standard_vec
            remap_val = (cos_val + 1) / 2
            crvs += tri_module.crvs_at_ratio(remap_val)
        return geo.Brep.CreatePlanarBreps(crvs, 0.001)


class Triangulizer:
    def __init__(self, mass: geo.Brep):
        self.mass = mass

    def calculate_angle(self, vec: geo.Vector3d):
        # type: (...) -> float
        return math.atan2(vec.Y, vec.X)  # atan2는 각도를 -pi에서 pi 사이로 반환

    def get_contours(self, height: float):
        # type: (...) -> list[geo.Curve]
        max_z = geo.Box(geo.Plane.WorldXY, self.mass).Z.Max  # type: ignore
        divide_count = math.floor(max_z / height)

        contour_crvs = []  # type: list[geo.Curve]

        for i in range(divide_count + 1):
            transform = geo.Transform.Translation(geo.Vector3d.ZAxis * i * height)  # type: ignore
            cutter_plane = geo.Plane.WorldXY
            cutter_plane.Transform(transform)
            crvs = list(geo.Brep.CreateContourCurves(self.mass, cutter_plane))  # type: ignore
            contour_crvs += crvs

        return contour_crvs

    def get_facade(self, divide_count: int, height: float):
        # type: (...) -> TriFacade
        contours = self.get_contours(height)
        pt_center = geo.AreaMassProperties.Compute(contours[0]).Centroid  # type: ignore
        transform = geo.Transform.Rotation(math.pi / divide_count, pt_center)  # type: ignore

        contour_pts = {}  # type: dict[geo.Curve, list[geo.Point3d]]
        for i, contour in enumerate(contours):
            params = contour.DivideByCount(divide_count, True)
            pts = [contour.PointAt(param) for param in params]
            if i % 2 == 1:
                pts = [transform * pt for pt in pts]
            contour_pts[contour] = pts

        tri_modules = []  # type: list[TriModule]
        for bottom_crv, top_crv in zip(contours[:-1], contours[1:]):
            pts = contour_pts[bottom_crv] + contour_pts[top_crv]
            pts.sort(key=lambda pt: self.calculate_angle(pt - pt_center))

            for i, pt_0 in enumerate(pts):
                pt_1 = pts[(i + 1) % len(pts)]
                pt_2 = pts[(i + 2) % len(pts)]

                if (pt_1 - pt_0).Z > 0:  # type: ignore
                    pt_0, pt_1, pt_2 = pt_2, pt_1, pt_0
                tri_modules.append(TriModule(pt_0, pt_1, pt_2))

        return TriFacade(tri_modules)
