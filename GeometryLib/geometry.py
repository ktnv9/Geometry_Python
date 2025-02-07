import math
import unittest


class Point:

    def __init__(self, x, y, z):

        # initialization
        self.x, self.y, self.z = x, y, z

    def __add__(self, other):

        # p1+p2
        return Point(self.x+other.x, self.y+other.y, self.z+other.z)
    
    def __sub__(self, other):

        # p1-p2 --> returns vector between the two points.
        return Vector(self.x-other.x, self.y-other.y, self.z-other.z)

    def __mul__(self, number):

        # p1*2 (mulitply with any number; order is important)
        return Point(self.x * number, self.y * number, self.z * number)
    
    def __truediv__(self, number):

        if not isinstance(number, (int, float)):
            raise TypeError("number must be numeric")

        if number == 0:
            raise ZeroDivisionError("Can not divide by zero")
        
        # p1*2 (divide with any number; order is important)
        return Point(self.x/number, self.y/number, self.z/number)
    
    def __neg__(self):

        # -p1
        return Point(-self.x, -self.y, -self.z)
    
    def __eq__(self, other):

        # p1 == p2
        return (abs(self.x-other.x) < 0.01) and (abs(self.y-other.y) < 0.01) and (abs(self.z-other.z) < 0.01)
    
    def __repr__(self):

        # (x, y, z)
        return f"({self.x}, {self.y}, {self.z})"
    
    def distace_to(self, other):

        # euclidean distance between two points.
        return (other-self).length()
    
    def point_at_distance(self, distace, direction_vector):

        # computes point at a distance from this point along the given direction vector
        return self + (direction_vector * distace)
    
    def directional_distance_to(self, other, direction_vector):

        # computes the directional distance between two points.
        return (other-self).dot_product(direction_vector)
    
    def polar_angle(self, other):
        return math.atan2(other.x-self.x, other.y-self.y)
    


class Vector:
    
    def __init__(self, x, y, z):
        # initialization
        self.x, self.y, self.z = x, y, z
    
    def __add__(self,other):

        # v1+v2
        return Vector(self.x+other.x, self.y+other.y, self.z+other.z)
    
    def __sub__(self, other):

        # v1-v2
        return Vector(self.x-other.x, self.y-other.y, self.z-other.z)
    
    def __mul__(self, number):

        # v*2 (multiply with any number; order is important)
        return Vector(self.x*number, self.y*number, self.z*number)
    
    def __truediv__(self, number):

        # v/2 (divide with any number; order is important)
        return Vector(self.x/number, self.y/number, self.z/number)
    
    def __neg__(self):

        # -v
        return Vector(-self.x, -self.y, -self.z)
    
    def __eq__(self, other):

        # v1 == v2
        return (abs(self.x-other.x) < 0.0001) and (abs(self.y-other.y) < 0.0001) and (abs(self.z-other.z) < 0.0001)

    def __repr__(self):

        # (x, y, z)
        return f"({self.x}, {self.y}, {self.z})"
    
    def length(self):

        # length of the vector
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def normalize(self):

        # unitize/normalize the vector
        vec_length = self.length()
        return self/vec_length 

    def is_unit_vector(self):

        # checks if the vector is a unit vector
        return abs(self.length() - 1) < 0.0001
    
    def dot_product(self, other):

        # v1.v2
        return (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
    
    def cross_product(self, other):

        # v1 x v2 (order is important)
        return Vector(self.y*other.z - self.z*other.y, self.z*other.x - self.x*other.z, self.x*other.y - self.y*other.x)
    
    def parallel(self, other):

        # checks if two vectors are parallel
        return self.angle_between(other) < 1 #1 deg tol
    
    def anti_parallel(self, other):

        # checks if two vectors are anti-parallel
        return abs(self.angle_between(other) - 180) < 1  # 1 deg tol

    def aligned(self, other):

        # checks if two vectors are aligned (parallel or anti-parallel)
        return self.parallel() or self.anti_parallel()
    
    def orthogonal(self, other):

        # checks if two vectors are perpendicular
        return abs(self.angle_between(other) - 90) < 1 # 1 deg tol

    def angle_between(self, other):

        # computes the angle between two vectors
        return math.degrees(math.acos(self.dot_product(other)/(self.length()*other.length())))

    def relative_position(self, other, orthog_plane_normal):
        cross_prod_vector = self.cross_product(other)
        
        rel_position = "RIGHT"
        if orthog_plane_normal.parallel(cross_prod_vector):
            rel_position = "LEFT"
            
        return rel_position
    
    def parallelogram_area(self, other):
        cross_prod_vector = self.cross_product(other)
        return cross_prod_vector.length()


class LineSegment:

    def __init__(self, point1, point2):
        self.point1, self.point2 = point1, point2

    def vector(self):
        return self.point2 - self.point1

    def length(self):
        vector = self.vector()
        return vector.length()

    def mid_point(self):
        return (self.point1 + self.point2)/2
    
    def parallel(self, other):
        return self.vector().parallel(other.vector())

    def anti_parallel(self, other):
        return self.vector().anti_parallel(other.vector())

    def aligned(self,other):
        return self.vector().aligned(other.vector())

    def orthogonal(self, other):
        return self.vector().orthogonal(other.vector())

    def intersection_point(self, other):
        
        # a, b are points of self
        # c, d are points of other
        # intersection occurs when two parametric equations become equal: a + t(b-a) = c + s(d-c)
        # solve for t & S by apply dot product of orthogonal vectors of self & other on both sides.

        a, b = self.point1, self.point2
        c, d = other.point1, other.point2

        vec_ab = self.vector()
        vec_cd = other.vector()
        vec_ca = a - c
        vec_ac = -vec_ca
        
        orthog_vec_ab = (vec_ab.cross_product(vec_cd)).cross_product(vec_ab)
        orthog_vec_cd = (vec_ab.cross_product(vec_cd)).cross_product(vec_cd)
        
        scalar_s = vec_ca.dot_product(orthog_vec_ab)/vec_cd.dot_product(orthog_vec_ab)
        scalar_t = vec_ac.dot_product(orthog_vec_cd)/vec_ab.dot_product(orthog_vec_cd)

        if (0 <= scalar_s <= 1) and (0 <= scalar_t <= 1): 

            return a + vec_ab * scalar_t

        return None

    def intersects(self, other):
        return self.intersection_point(other) is not None

    def cross_product(self, other):
        return self.vector().cross_product(other.vector())

    def angle_between(self, other):
        return self.vector().angle_between(other.vector())

    def contains_point(self, point):
        ls_length = self.length()
        p1_p_ls_length = LineSegment(self.point1, point).length()
        p_p2_ls_length = LineSegment(point, self.point2).length() 
        return abs(ls_length - (p1_p_ls_length+p_p2_ls_length)) < 0.00001

class Ray:
    def __init__(self, start_point, direction_vector):
        self.start_point = start_point
        self.direction_vector = direction_vector

    def parallel(self, other):
        return self.direction_vector.parallel(other.direction_vector)
    
    def anti_parallel(self, other):
        return self.direction_vector.anti_parallel(other.direction_vector)

    def orthogonal(self, other):
        return self.direction_vector.orthogonal(other.direction_vector)

    def aligned(self, other):
        return self.direction_vector.aligned(other.direction_vector)

    def intersection_point(self, other):
        
        # a & b are the start points of self & other.
        # v1, v2 are the respective ray direction vectors.
        # intersection occurs when two parametric equations become equal: a + t(v1) = b + s(v2)
        # solve for t & S by apply dot product of orthogonal vectors of self & other on both sides.

        a = self.start_point
        v1 = self.direction_vector
        b = other.start_point
        v2 = other.direction_vector

        orthog_v1 = (v1.cross_product(v2)).cross_product(v1)
        orthog_v2 = (v1.cross_product(v2)).cross_product(v2)

        scalar_s = (a-b).dot_product(orthog_v1)/v2.dot_product(orthog_v1)
        scalar_t = (b-a).dot_product(orthog_v2)/v1.dot_product(orthog_v2)

        if scalar_s >= 0 and scalar_t >= 0:
            return a + v1 * scalar_t

        return None

    def intersects(self, other):
        return self.intersection_point(other) is not None

    def cross_product(self, other):
        return self.direction_vector.cross_product(other.direction_vecotr)

    def angle_between(self, other):
        return self.direction_vector.angle_between(other.direction_vector)

    def contains_point(self, point):
        vector_p_sp = point - self.start_point
        return self.direction_vector.parallel(vector_p_sp)
        
class Rectangle:

    def __init__(self, width, height, width_axis = Vector(1,0,0), height_axis = Vector(0,1,0), origin_point = Point(0,0,0), origin_tag = "LB"):
        self.width, self.height = width, height
        self.width_axis, self.height_axis = width_axis, height_axis
        self.orign_point, self.origin_tag = origin_point, origin_tag

    def area(self):
        return self.width * self.height

    def perimeter(self):
        return 2*(self.width + self.height)

    def contains_point(self, point):

        # this method works for generic mis-algined rectangles as well.
        lb_corner = self.anchor_point("LB")
        rt_corner = self.anchor_point("RT")

        dist_along_width = lb_corner.directional_distance_to(point, self.width_axis)
        dist_along_height = rt_corner.directional_distance_to(point, self.height_axis)

        return (dist_along_width < self.width) and (dist_along_height < self.height)

    def anchor_point(self, anchor_tag):
        # possible anchor points = Left-Bottom, Right-Bottom, Left-Top, Right-Top, Center, Left-Center, Right-Center, Top-Center, Bottom-Center
        
        # move from origin tag to center; then move from center to desired anchor tag
        steps_to_ctr = self._get_steps(self.origin_tag, "CTR")
        steps_from_ctr = self._get_steps("CTR", anchor_tag)
        steps = steps_to_ctr + steps_from_ctr
        for step in steps:
            distance = steps[0]
            direction = steps[1]
            dest_pt = dest_pt.point_at_distance(distance, direction)
        return dest_pt

    def _get_steps_to_and_fro_from_ctr(self, source_or_dist):
        
        # create a dictionary: key  = source_or_dist; value = [(distance1, direction1), (distance2, direction2)]
        pass

    def intersects_with(self, other):
        
        if self.width_axis.aligned(other.width_axis) and self.height_axis.aligned(other.height_axis):
            pass
        else:
            # apply Separating-axis-theorem
            pass

    def intersecting_area(self, other):
        pass

class Polygon:

    def __init__(self, points):
        self.points = points

    def vertex_count(self):
        return len(self.points)

    def perimeter(self):
        pass

    def area(self):
        pass

    def diameter(self):
        pass

    def centriod(self):
        x_coords = [point.x for point in self.points]
        y_coords = [point.y for point in self.points]
        z_coords = [point.z for point in self.points]
        return Point(sum(x_coords)/len(x_coords), sum(y_coords)/len(y_coords), sum(z_coords)/len(z_coords))

    def is_convex(self):
        
        num_points = len(self.pionts)
        for i in range(num_points+1):

            point1, point2 = self.point
            if i == num_points:
                pass

    def corner_curvatures(self):
        # list of curvatures at each point of the polygon.
        pass

    def corner_angles(self):
        pass

    def is_regular(self):
        pass
        # check if all corner angles are the same.

    def contains_point(self):
        
        # implement ray casting algorithm
        pass

    def intersects_with(self, other):
        pass

    def intersecting_area(self, other):
        pass

    def convex_hull(self):
        pass

    def bounding_box(self):
        x_coords = [point.x for point in self.points]
        y_coords = [point.y for point in self.points]
        min_x, max_x = min(x_coords), max(x_coords)
        min_y, max_y = min(y_coords), max(y_coords)
        return Rectangle(width=max_x-min_x, height=max_y-min_y, origin_point=Point(min_x, min_y, self.points[0].z))

    def bounding_circle(self):
        pass
        # return circle object

class Circle:

    def __init__(self, center = Point(0,0,0), radius=1):
        # default unit circle.
        self.center = center
        self.radius = radius
        
    def perimeter(self):
        # 2*pi*r
        return 2*math.pi*self.radius
    
    def area(self):

        # pi*r^2
        return math.pi * self.radius**2

    def intersects_with(self, other):

        # distance between centers < sum of radii
        return (self.center-other.center).length() < (self.radius + other.radius)
    
    def contains_point(self, point):

        # distance between center and point < radius
        return (point - self.center).length() < self.radius

    
    


                





class TestPoint(unittest.TestCase):

    def test_add(self):
        p1 = Point(1,2,3)
        p2 = Point(4,5,6)
        self.assertEqual(p1 + p2, Point(5,7,9))

    def test_sub(self):
        p1 = Point(1,2,3)
        p2 = Point(4,5,6)
        self.assertEqual(p2 - p1, Point(3,3,3))

    def test_mul(self):
        p = Point(2,3,4)
        self.assertEqual(p * 2, Point(4,6,8))
    
    def test_truediv(self):
        p = Point(4, 6, 8)
        self.assertEqual(p/2, Point(2,3,4))
        with self.assertRaises(ZeroDivisionError):
            p/0
        with self.assertRaises(TypeError):
            p/"a"

    def test_neg(self):
        p = Point(3,4,5)
        self.assertEqual(-p,Point(-3,-4,-5))

    def test_eq(self):
        p1 = Point(3,4,5)
        p2 = Point(3.00001,4.00002,5)
        self.assertEqual(p1 == p2, True)

if __name__ == "__main__":

   unittest.main()
