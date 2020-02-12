from simtk.unit import angstroms

def clustering_fpocket(alpha_spheres_set, maximum_distance_first_clustering=1.8*angstroms,
        clustering_method=None):

        #clustering_method='a', 's'
    pass

def distance_connectivity(alpha_spheres_set, threshold=1.8*angstroms):

    from molmodmt.distances import contact_map
    from networkx import from_numpy_matrix, connected_components

    cmap = contact_map(item_1=alpha_spheres_set.centers, threshold=threshold)
    graph = from_numpy_matrix(cmap[0])
    clusters= list(connected_components(graph))
    del(graph)
    return clusters

def neighbor_points_connectivity(alpha_spheres_set):

    from numpy import zeros
    from itertools import combinations
    from networkx import from_numpy_matrix, connected_components

    cmap = zeros([alpha_spheres_set.n_alpha_spheres, alpha_spheres_set.n_alpha_spheres], dtype=bool)

    points_in_spheres={}

    for sphere_index in range(alpha_spheres_set.n_alpha_spheres):
        for neighbor_point_index in alpha_spheres_set.neighbor_points[sphere_index]:
            try:
                points_in_spheres[neighbor_point_index].append(sphere_index)
            except:
                points_in_spheres[neighbor_point_index]=[sphere_index]

    for spheres_indices_in_contact in points_in_spheres.values():
        for ii,jj in combinations(spheres_indices_in_contact,2):
            cmap[ii,jj]=True
            cmap[jj,ii]=True

    graph = from_numpy_matrix(cmap)
    clusters= list(connected_components(graph))
    del(graph)
    return clusters

