
class alpha_spheres_set():

    def __init__(self, points=None, minimum_radius=None, maximum_radius=None):

        # fpocket
        # minimum_radius = 3.2 * angstroms
        # maximum_radius = 7.6 * angstroms

        # alphaspace
        # minimum_radius = 3.2 * angstroms
        # maximum_radius = 5.4 * angstroms

        self.points=None
        self.n_points=None
        self.vertices=None
        self.neighbor_points=None
        self.radii=None
        self.n_alpha_spheres=None
        self.cluster=None

        if points is not None:

            from scipy.spatial import Voronoi
            from molmodmt import neighbors_lists
            from numpy import array, zeros

            unit_length = points.unit

            self.points = points
            self.n_points = points.shape[0]
            self.vertices = []
            self.neighbor_points = []
            self.radii = []
            self.n_alpha_spheres = 0

            tmp_vertices = Voronoi(self.points[:,:]).vertices*unit_length
            tmp_neighbor_points, tmp_radii = neighbors_lists(item_1=tmp_vertices, item_2=self.points[:,:], num_neighbors=4)
            tmp_neighbor_points = tmp_neighbor_points[0,:,:]
            tmp_radii = tmp_radii[0,:,0]

            if minimum_radius is not None:
                mask = (tmp_radii>minimum_radius)
                tmp_vertices = tmp_vertices[mask,:]
                tmp_neighbor_points = tmp_neighbor_points[mask,:]
                tmp_radii = tmp_radii[mask]

            if maximum_radius is not None:
                mask = (tmp_radii<maximum_radius)
                tmp_vertices = tmp_vertices[mask,:]
                tmp_neighbor_points = tmp_neighbor_points[mask,:]
                tmp_radii = tmp_radii[mask]

            self.radii = tmp_radii
            self.vertices = tmp_vertices
            self.neighbor_points = tmp_neighbor_points
            self.n_alpha_spheres = tmp_vertices.shape[0]

    def removing_alpha_spheres(self, alpha_sphere_indices=None):

        from numpy import ones
        mask = ones([self.n_alpha_spheres], dtype=bool)
        mask[alpha_sphere_indices] = False
        self.vertices = self.vertices[mask,:]
        self.neighbor_points = self.neighbor_points[mask,:]
        self.radii = self.radii[mask]
        self.n_alpha_spheres = self.vertices.shape[0]

    def removing_least_populated_clusters(self, clusters, threshold=1):

        to_be_removed=[]
        clusters_to_be_kept=[]
        for cluster in clusters:
            if len(cluster)<=threshold:
                to_be_removed.extend(list(cluster))
            else:
                clusters_to_be_kept.append(list(cluster))

        translation = {}
        new_clusters=[]

        ii=0
        for alpha_sphere_index in range(self.n_alpha_spheres):
            if alpha_sphere_index not in to_be_removed:
                translation[alpha_sphere_index]=ii
                ii+=1

        for cluster in clusters_to_be_kept:
            new_cluster = [translation[ii] for ii in cluster]
            new_clusters.append(set(new_cluster))

        self.removing_alpha_spheres(to_be_removed)

        return new_clusters

