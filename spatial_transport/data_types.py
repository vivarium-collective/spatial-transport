"""Schema types for spatial transport.

The `volumetric` type is shared with cdFBA and imported from there rather than
redefined, so both packages apply counts/concentrations updates identically.
"""
from cdFBA.data_types import Volumetric, register_types as register_cdfba_types

from spatial_transport.processes import register_processes

#=================
# Spatial Types
#=================
edge_type = {
    "neighbors": "list[string]",
    "surface_area": "float",
    "position": "list[float]",
    "periodic": "boolean",
}

compartment_type = {
    "Shared Environment": "volumetric",
    "position": "list[float]",
    "boundaries": "list[string]",
}

# advection vectors are recomputed from scratch each step, so they replace the
# stored vector instead of accumulating like a plain list would
advection_vector_type = "overwrite[list[float]]"


def register_types(core):
    core = register_cdfba_types(core)

    core.register_type("edge_type", edge_type)
    core.register_type("compartment", compartment_type)
    core.register_type("advection_vector", advection_vector_type)

    return register_processes(core)
