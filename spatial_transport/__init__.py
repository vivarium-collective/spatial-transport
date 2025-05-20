from spatial_transport.processes import register_processes

def conditional_apply(schema, current, update, key, core):
    if key in update:
        applied = core.apply(
            schema[key],
            current[key],
            update[key])
    else:
        applied = current[key]

    return applied

def volumetric_update(schema, current, update, top_schema, top_state, path, core):
    updated_counts = conditional_apply(schema, current, update, "counts", core)
    updated_volume = conditional_apply(schema, current, update, "volume", core)
    updated_concentrations = {}

    for key, counts in updated_counts.items():
        updated_concentrations[key] = counts / updated_volume

    applied = {
        "counts": updated_counts,
        "concentrations": updated_concentrations,
        "volume": updated_volume,
    }
    # # Temp
    # import ipdb; ipdb.set_trace()
    # # Temp
    return applied

volumetric_type = {
    "concentrations":"map[float]",
    "counts":"map[float]",
    "volume": {"_type": "float", "_default": 1.0},
    "_apply": volumetric_update
}

edge_type = {
    "neighbors": "list[string]",
    "surface_area": "float",
}

compartment_type = {
    "Shared Environment": "volumetric",
    "position": "list[float]",
}

def register_types(core):
    core.register("volumetric", volumetric_type)
    core.register("edge_type", edge_type)
    core.register("compartment", compartment_type)
    return register_processes(core)