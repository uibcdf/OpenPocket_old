
def alpha_spheres(item, selection='all', frame_indices='all', syntaxis='MolModMT', minimum_radius =
        None, maximum_radius = None):

    """alpha_spheres(item, selection='all', frame_indices='all', syntaxis='MolModMT',
    minimum_radius=None, maximum_radius=None)

    Get the set of alpha spheres.

    Given a molecular model the whole set of alpha spheres is computed for the atoms and frames
    selected.

    Parameters
    ----------

    item: molecular model
        Molecular model in any supported form by MolModMT(see: XXX).

    selection: str, default='heavy atoms'
       Selection criterion given by a string following any of the selection syntaxis parsable by
       MolModMT, or a list of atoms (see: XXX).

    syntaxis: str, default='MolModMT'
       Syntaxis used to write the argument `selection`. The current options supported by MolModMt
       can be found in (see: XXX).

    Returns
    -------

    Class AlphaSpheres
        The output is an object of class AlphaSpheres (see: XXX).

    Examples
    --------

    See Also
    --------

    Notes
    -----

    """

    # fpocket
    # minimum_radius = 3.2 * angstroms
    # maximum_radius = 7.6 * angstroms

    # alphaspace
    # minimum_radius = 3.2 * angstroms
    # maximum_radius = 5.4 * angstroms

    from molmodmt import select, get
    from numpy import where

    atom_indices = select(item, selection=selection, syntaxis=syntaxis)
    coordinates = get(item, target='atom', indices=atom_indices, frame_indices=frame_indices, coordinates=True)

    alpha_spheres_set = AlphaSpheresSet(points=coordinates)

    if minimum_radius is not None:
        indices_to_remove = where(alpha_spheres_set.radii < minimum_radius)
        alpha_spheres_set.remove(indices_to_remove)

            if maximum_radius is not None:
                mask = (tmp_radii<maximum_radius)
                tmp_vertices = tmp_vertices[mask,:]
                tmp_neighbor_points = tmp_neighbor_points[mask,:]
                tmp_radii = tmp_radii[mask]




class AlphaSpheresSet():

    def __init__(self, points=None):

        """AlphaSpheresSet(points=None)

        Class representing a set of Alpha spheres.

        Paragraph with detailed explanation.

        Parameters
        ----------

        points: numpy array
            Float numpy array with shape: [natoms, 3]

        Returns
        -------

        Class

        Examples
        --------

        See Also
        --------

        Notes
        -----

        """

        self.points=None
        self.n_points=None
        self.vertices=None
        self.neighbor_points=None
        self.radii=None
        self.n_alpha_spheres=None

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

            self.vertices = Voronoi(self.points[:,:]).vertices*unit_length
            self.neighbor_points, self.radii = neighbors_lists(item_1=self.vertices, item_2=self.points, num_neighbors=4)
            self.neighbor_points = self.neighbor_points[0,:,:]
            self.radii = self.radii[0,:,0]

    def remove(self, indices=None):

        from numpy import ones
        mask = ones([self.n_alpha_spheres], dtype=bool)
        mask[indices] = False
        self.vertices = self.vertices[mask,:]
        self.neighbor_points = self.neighbor_points[mask,:]
        self.radii = self.radii[mask]
        self.n_alpha_spheres = self.vertices.shape[0]

class AlphaSpheres_OLD():

    def __init__(self, points=None, minimum_radius= 3.2 * angstroms, maximum_radius= 7.6 * angstroms):

        """AlphaSpheres(item, selection='all', output_indices='atom', syntaxis='MolModMT')
        Get the atom indices corresponding to a selection criterion.
        Paragraph with detailed explanation.
        Parameters
        ----------
        item: molecular model
            Molecular model in any supported form (see: :doc:`/Forms`). The object being acted on by the method.
        selection: str, default='all'
           Selection criterion given by means of a string following any of the selection syntaxis parsable by MolModMT.
        output_indices: str, default='atom'
           The output list can correspond to 'atom', 'group', 'component', 'molecule', 'chain' or 'entity'
           indices.
        syntaxis: str, default='MolModMT'
           Syntaxis used to write the argument `selection`. The current options supported by MolModMt
           can be found in :doc:`/Atoms_Selection`.
        Returns
        -------
        Numpy array of integers
            List of indices in agreement with the selection criterion applied over `item`. The nature
            of the indices is chosen with the impot argument 'output_indices': 'atom' (default),
            'group', 'component', 'molecule', 'chain' or 'entity'.
        Examples
        --------
        :doc:`/Atoms_Selection`
        See Also
        --------
        Notes
        -----
        """

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
