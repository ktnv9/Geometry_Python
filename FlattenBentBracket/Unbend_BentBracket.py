
# problem statement
#   Given geometric information of the bent bracket profile, create the correspondig geometric information of the flattened structure of Bent bracket.


"""
edge geometric data
face geometric data

    create part_edge_graph (peg)
    create part_face_graph (pfg)

assumption
    only the upper layer of the bent bracket is considered
    it consists of sub_peg, sub_pfg
    if we can unbend the sub graphs, it is sufficient. The bottom layer can be adjusted using the thickness of the bent bracket.

    if we remove the bend cyl faces from the ul (upper layer) we will be left with only segments.
    so we can create a bb_graph.

    un-bending algorithm
    Loop until there is only one segment in the graph
        identify the leaf segments of the garph.
        for each leaf segment
            identify the corresp non leaf segment which shares the same bend cyl face.
            if only leaf segments are available, consider one segment as the fixed segment and the other as moving segment.
            flatten the leaf segement
            merge the leaf segement and the bend into the non-leaf segement.

            
    flatten one segment wrt another fixed segement

        flatten the bend edge wise --> compute the corresp vertices using point at distance method
        use the dictionary or some chain type of data structure to keep track of the transformation each vertex, edge, face are taking 

"""


class BentBracket:
    
    def bends():
        pass

    def segments():
        pass

    def edge_graph():
        pass

    def face_graph():
        pass

    def tree():
        pass

    def flatten():
        pass

class Bend:

    id = 2
    icf = "FACE A"
    ocf = "FACE B"

    def __init__(self, id, icf, ocf):
        pass

    def axis():
        pass

    def faces():
        pass

    def profile_faces():
        pass

    def side_faces():
        pass

    def co_faces():
        pass

    def co_edges():
        pass

    def angle():
        pass

class BendSegment:

    def __init__(self, id):
        self.id = id

    def faces():
        pass

    def LCS():
        pass

    def profile_faces():
        pass

    def side_faces():
        pass





