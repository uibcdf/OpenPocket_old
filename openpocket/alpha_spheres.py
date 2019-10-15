
class get_alpha_spheres_set():

    def __init__(self, points=None, minimum_radius=None, maximum_radius=None):

        # fpocket
        # minimum_radius = 3.2 * angstroms
        # maximum_radius = 7.6 * angstroms

        # alphaspace
        # minimum_radius = 3.2 * angstroms
        # maximum_radius = 5.4 * angstroms

        self.points=None
        self.n_frames=None
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
            self.n_frames = points.shape[0]
            self.n_points = points.shape[1]
            self.vertices = []
            self.neighbor_points = []
            self.radii = []
            self.n_alpha_spheres = zeros([self.n_frames], dtype='int')

            for ii in range(self.n_frames):

                tmp_vertices = Voronoi(self.points[ii,:,:]).vertices*unit_length
                tmp_neighbor_points, tmp_radii = neighbors_lists(item_1=tmp_vertices, item_2=self.points[ii,:,:], num_neighbors=4)
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

                self.radii.append(tmp_radii)
                self.vertices.append(tmp_vertices)
                self.neighbor_points.append(tmp_neighbor_points)
                self.n_alpha_spheres[ii] = tmp_vertices.shape[0]

def removing_alpha_spheres(alpha_spheres_set, alpha_sphere_indices=None, frame_index=None):

    from numpy import ones
    mask = ones([alpha_spheres_set.n_alpha_spheres[frame_index]], dtype=bool)
    mask[alpha_sphere_indices] = False
    alpha_spheres_set.vertices[frame_index]=alpha_spheres_set.vertices[frame_index][mask,:]
    alpha_spheres_set.neighbor_points[frame_index]=alpha_spheres_set.neighbor_points[frame_index][mask,:]
    alpha_spheres_set.radii[frame_index]=alpha_spheres_set.radii[frame_index][mask]
    return alpha_spheres_set

def removing_least_populated_clusters(alpha_spheres_set, clusters, threshold=1):

    for frame_index in range(alpha_spheres_set.n_frames):

        to_be_removed=[]
        for cluster in clusters:
            if len(cluster)<=threshold:
                to_be_removed.extend(list(cluster))

    return to_be_removed

