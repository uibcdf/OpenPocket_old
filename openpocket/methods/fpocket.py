
def fpocket(item=None, minimum_radius=None, maximum_radius=None)

    from openpocket import get_alpha_spheres_set

    alpha_spheres_set = get_alpha_spheres_set(item, minimum_radius=minimum_radius, maximum_radius=maximum_radius)

    return alpha_spheres_set

