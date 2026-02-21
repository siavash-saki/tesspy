"""
Voronoi diagram tessellation functions.
"""

from collections import defaultdict
from collections.abc import Generator

import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon


def voronoi_polygons(sp_voronoi_obj: Voronoi, diameter: float) -> Generator[Polygon, None, None]:
    """
    Convert a scipy Voronoi object into shapely Polygons.

    Parameters
    ----------
    sp_voronoi_obj : scipy.spatial.Voronoi
        Voronoi object created from point coordinates
    diameter : float
        Controls the extent of infinite Voronoi regions. Should be large
        enough to cover the bounding box of the input points.

    Yields
    ------
    shapely.geometry.Polygon
        Voronoi polygons (finite and infinite regions)
    """
    centroid = sp_voronoi_obj.points.mean(axis=0)

    ridge_direction: dict = defaultdict(list)
    for (p, q), rv in zip(sp_voronoi_obj.ridge_points, sp_voronoi_obj.ridge_vertices):
        u, v = sorted(rv)
        if u == -1:
            tangent = sp_voronoi_obj.points[q] - sp_voronoi_obj.points[p]
            normal = np.array([-tangent[1], tangent[0]]) / np.linalg.norm(tangent)
            midpoint = sp_voronoi_obj.points[[p, q]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - centroid, normal)) * normal
            ridge_direction[p, v].append(direction)
            ridge_direction[q, v].append(direction)

    for i, r in enumerate(sp_voronoi_obj.point_region):
        region = sp_voronoi_obj.regions[r]
        if -1 not in region:
            # finite regions
            yield Polygon(sp_voronoi_obj.vertices[region])
            continue

        # infinite regions
        inf_vertex_id = region.index(-1)
        previous_vertex_id = region[(inf_vertex_id - 1) % len(region)]
        next_vertex_id = region[(inf_vertex_id + 1) % len(region)]
        if previous_vertex_id == next_vertex_id:
            dir_j, dir_k = ridge_direction[i, previous_vertex_id]
        else:
            (dir_j,) = ridge_direction[i, previous_vertex_id]
            (dir_k,) = ridge_direction[i, next_vertex_id]

        length = 2 * diameter / np.linalg.norm(dir_j + dir_k)

        finite_part = sp_voronoi_obj.vertices[
            region[inf_vertex_id + 1 :] + region[:inf_vertex_id]
        ]
        extra_edge = [
            sp_voronoi_obj.vertices[previous_vertex_id] + dir_j * length,
            sp_voronoi_obj.vertices[next_vertex_id] + dir_k * length,
        ]
        yield Polygon(np.concatenate((finite_part, extra_edge)))
