
class alpha_spheres():

    def __init__(self, points=None):

        from scipy.spatial import Voronoi
        from molmodmt import neighbors_lists

        unit_length = points.unit

        self.points = points
        self.vertices = Voronoi(self.points).vertices*unit_length
        self.neighbor_points, self.radii = neighbors_lists(item_1=self.vertices, item_2=self.points, num_neighbors=4)


