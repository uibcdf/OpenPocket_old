
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
        indices_to_remove = where(alpha_spheres_set.radii > maximum_radius)
        alpha_spheres_set.remove(indices_to_remove)

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
    vertices: xxx
        No recuerdo ahora.
    neighbor_points: xxx
        Indices of points in the surface of each alpha-sphere.
    radii: numpy.ndarray
        Array with the radii of alpha-spheres.
    n_alpha_spheres: int
        Number of alpha-spheres in the set.

    Methods
    -------
    remove(indices)
        Removing alpha-spheres from the set

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
        self.vertices = self.vertices[mask,:]
        self.neighbor_points = self.neighbor_points[mask,:]
        self.radii = self.radii[mask]
        self.n_alpha_spheres = self.vertices.shape[0]

