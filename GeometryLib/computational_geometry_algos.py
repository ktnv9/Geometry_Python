from .geometry import *

def convex_hull(points) -> Polygon:
    
    # given a set of 2D points, find the smallest convex polygon that encloses all the points.

    # applications: collision detection, robotics, image processing(shape analysis), geographical mapping (boundary set of locations)

    # Graham Scan algorithm | RC - O(N log N) | SC - O(N)

        # find the lowest point (mark it as anchor point)
        anchor_point = min(points, key=lambda p:(p.y, p.x)) #sort by y and then x.

        # sort all other points by polar angle with respect to the anchor.
        sorted_points = sorted(points, key=lambda p: (anchor_point.polar_angle(p), (p.x, p.y)))

        # traverse the sorted list and use a stack to maintain the convex hull, discarding the points that make a right turn (non-convex)
        hull = []
        for point in sorted_points:
            hull.append(point)
            
            while len(hull) >= 3:
                vect1 = hull[-3] - hull[-2]
                vect2 = hull[-2] - hull[-1]
                if vect1.relative_position(vect2, Vector(0,1)) == "RIGHT":
                    hull.pop()
            
        return Polygon(hull)