
class alpha_spheres():

    def __init__(points):

        from scipy.spatial import Voronoi
        from molmodmt import neighbors_lists

        unit_length = points.unit

        self.points = points
        self.vertices = Voronoi(self.points).vertices
        self.neighbor_points, self.radii = m3t.neighbors_lists(item_1=self.vertices, item_2=self.points, num_neighbors=4)


