from simtk.unit import angstroms

def clustering_fpocket(alpha_spheres_set, maximum_distance_first_clustering=1.8*angstroms,
        clustering_method=None):

        #clustering_method='a', 's'
    pass

def distance_between_vertices(alpha_spheres_set, maximum_distance=1.8*angstroms):

    from molmodmt.distances import contact_map
    from networkx import from_numpy_matrix, connected_components

    clusters=[]

    for ii in range(alpha_spheres_set.n_frames):
        cmap = contact_map(item_1=alpha_spheres_set.vertices[ii], threshold=maximum_distance)
        graph = from_numpy_matrix(cmap[0])
        clusters.append(list(connected_components(graph)))

    del(graph)
    return clusters

def neighbor_points_connectivity(alpha_spheres_set):

    for frame_index in range(alpha_spheres_set.n_frames):

        points_in_spheres={}

        for sphere_index in range(alpha_spheres_set.n_spheres):
            for neighbor_point_index in alpha_spheres_set.neighbours_points:
                try:
                    points_in_spheres[neighbor_point_index].append(sphere_index)
                except:
                    points_in_spheres[neighbor_point_index]=[sphere_index]

    pass


