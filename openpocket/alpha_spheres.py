
def get_alpha_spheres_set(item, selection='all', frame_indices='all', syntaxis='MolSysMT', minimum_radius =
        None, maximum_radius = None):

    """get_alpha_spheres_set(item, selection='all', frame_indices='all', syntaxis='MolSysMT',
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

    from molsysmt import select, get

    atom_indices = select(item, selection=selection, syntaxis=syntaxis)
    coordinates = get(item, target='atom', indices=atom_indices, frame_indices=frame_indices, coordinates=True)

    alpha_spheres_set = AlphaSpheresSet(points=coordinates[0,:,:])

    if minimum_radius is not None:
        alpha_spheres_set.remove_small_alpha_spheres(minimum_radius)

    if maximum_radius is not None:
        alpha_spheres_set.remove_big_alpha_spheres(minimum_radius)

    return alpha_spheres_set

class AlphaSpheresSet():

    """Set of alpha-spheres

    Object with a set of alpha-spheres and its main attributes such as radius, etc.

    Attributes
    ----------
    points : numpy.ndarray ([n_points,3])
        Array of coordinates in space of all points used to generate the set of alpha-spheres.
    n_points: int
        Number of points in space to generate the set of alpha-spheres.
    centers: xxx
        No recuerdo ahora.
    radii: numpy.ndarray
        Array with the radii of alpha-spheres.
    points_of_alpha_sphere: xxx
        Indices of points in the surface of each alpha-sphere.
    n_alpha_spheres: int
        Number of alpha-spheres in the set.

    Methods
    -------
    remove(indices)
        Removing alpha-spheres from the set
    get_points_in_surfaces(indices)
        Points in contact with a subset of alpha-spheres
    view(indices)
        3D spatial visualization of points and alpha-spheres

    """

    def __init__(self, points=None):

        """Creating a new instance of AlphaSpheresSet

        With an array of three dimensional positions (`points`) the resultant set of alpha-spheres is returned
        as a class AlphaSpheresSet.

        Parameters
        ----------

        points: numpy array
            Float numpy array with shape: [natoms, 3]

        Examples
        --------

        See Also
        --------

        Notes
        -----

        """

        self.points=None
        self.n_points=None
        self.centers=None
        self.points_of_alpha_sphere=None
        self.radii=None
        self.n_alpha_spheres=None
        self.voronoi=None

        if points is not None:

            from scipy.spatial import Voronoi
            from molsysmt import neighbors_lists
            from numpy import array, zeros, sort

            unit_length = points.unit

            self.points = points
            self.n_points = points.shape[0]
            self.centers = []
            self.points_of_alpha_sphere = []
            self.radii = []
            self.n_alpha_spheres = 0

            self.voronoi = Voronoi(self.points[:,:])
            self.centers = self.voronoi.vertices*unit_length
            self.points_of_alpha_sphere, self.radii = neighbors_lists(item_1=self.centers, item_2=self.points, num_neighbors=4)
            self.points_of_alpha_sphere = sort(self.points_of_alpha_sphere[0,:,:])
            self.radii = self.radii[0,:,0]
            self.n_alpha_spheres = self.centers.shape[0]

            ## This code computes points_of_alpha_sphere with out computing distances
            ## This should be cheaper
            #vor = Voronoi(self.points[:,:])
            #n_vertices = vor.vertices.shape[0]
            #n_regions = len(vor.regions)
            #n_points = vor.npoints
            #points_of_alpha_sphere = [[] for ii in range(n_vertices)]
            #region_point={vor.point_region[ii]:ii for ii in range(n_points)}
            #for region_index in range(n_regions):
            #    region=vor.regions[region_index]
            #    if len(region)>0:
            #        point_index=region_point[region_index]
            #        for vertex_index in region:
            #            if vertex_index != -1:
            #                points_of_alpha_sphere[vertex_index].append(point_index)

    def remove(self, indices):

        """Removing alpha-spheres from the set

        The method removes from the set those alpha-spheres specified by the input argument
        `indices`.

        Parameters
        ----------
        indices : numpy.ndarray, list or tuple (dtype:ints)
            List, tuple or numpy.ndarray with the integer numbers corresponding to the alpha-sphere
            indices to be removed from the set.

        Examples
        --------

        """

        from numpy import ones
        mask = ones([self.n_alpha_spheres], dtype=bool)
        mask[indices] = False
        self.centers = self.centers[mask,:]
        self.points_of_alpha_sphere = self.points_of_alpha_sphere[mask,:]
        self.radii = self.radii[mask]
        self.n_alpha_spheres = self.centers.shape[0]

    def remove_small_alpha_spheres(self, minimum_radius):

        """Removing smaller alpha-spheres than a given radius

        The method removes from the set those alpha-spheres with their radii shorter than a minimum
        radius.

        Parameters
        ----------
        minimum_radius : Quantity
            Those spheres with a radius smaller than this minimum radius will be removed from
            the set.

        Examples
        --------

        """

        from numpy import where
        indices_to_remove = where(self.radii < minimum_radius)
        self.remove(indices_to_remove)

    def remove_big_alpha_spheres(self, maximum_radius):

        """Removing alpha-spheres with large radii

        The method removes from the set those alpha-spheres with their radii larger than a maximum
        radius.

        Parameters
        ----------
        maximum_radius : Quantity
            Those spheres with a radius smaller than this minimum radius will be removed from
            the set.

        Examples
        --------

        """

        from numpy import where
        indices_to_remove = where(self.radii > maximum_radius)
        self.remove(indices_to_remove)

    def get_points_of_alpha_spheres(self, indices):

        """Get the points in contact with a subset of alpha-spheres

        The list of point indices accounting for the points in contact with a subset of alpha-spheres is calculated.

        Parameters
        ----------
        indices : numpy.ndarray, list or tuple (dtype:ints)
            List, tuple or numpy.ndarray with the alpha-sphere indices defining the subset.

        Return
        ------
        points_of_alpha_spheres : list
            List of point indices in contact with one or more alpha-spheres of the subset.

        Examples
        --------

        >>> import openpocket as opp
        >>> from simtk import unit
        >>> points = ([[-1.,  2.,  0.],
        >>>            [ 0.,  2.,  1.],
        >>>            [ 1., -2.,  1.],
        >>>            [ 0.,  1.,  1.],
        >>>            [ 0.,  0.,  0.],
        >>>            [-1., -1.,  0.]]) * unit.angstrom
        >>> aspheres = opp.alpha_spheres.AlphaSpheresSet(points)
        >>> aspheres.get_points_of_alpha_spheres([1,3])
        [0,2,3,4,5]

        """

        point_indices = set([])

        for index in indices:
            point_indices = point_indices.union(set(self.points_of_alpha_sphere[index]))

        return list(point_indices)

    def view(self, indices='all'):

        """3D spatial view of alpha-spheres and points

        An NGLview view is returned with alpha-spheres (gray color) and points (red color).

        Parameters
        ----------
        indices : numpy.ndarray, list or tuple (dtype:ints)
            List, tuple or numpy.ndarray with the alpha-sphere indices defining the subset.

        Returns
        -------
        view : nglview
            View object of NGLview.

        Examples
        --------

        >>> import openpocket as opp
        >>> from simtk import unit
        >>> points = ([[-1.,  2.,  0.],
        >>>            [ 0.,  2.,  1.],
        >>>            [ 1., -2.,  1.],
        >>>            [ 0.,  1.,  1.],
        >>>            [ 0.,  0.,  0.],
        >>>            [-1., -1.,  0.]]) * unit.angstrom
        >>> aspheres = opp.alpha_spheres.AlphaSpheresSet(points)
        >>> view = aspheres.view([1,3])
        >>> view


        """

        import nglview as nv

        view = nv.NGLWidget()

        point_indices = []

        if indices=='all':
            indices=range(self.n_alpha_spheres)
            point_indices=range(self.n_points)
        else:
            point_indices=self.get_points_of_alpha_spheres(indices)

        for index in point_indices:
            atom_coordinates = self.points[index,:]._value
            view.shape.add_sphere(list(atom_coordinates), [0.8,0.0,0.0], 0.2)

        for index in indices:
            sphere_coordinates = self.centers[index,:]._value
            sphere_radius = self.radii[index]._value
            view.shape.add_sphere(list(sphere_coordinates), [0.8,0.8,0.8], sphere_radius)

        return view

