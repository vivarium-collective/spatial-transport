from spatial_transport.processes.diffusion import SimpleDiffusion
from spatial_transport.processes.advection import SimpleAdvection, DynamicAdvection, Peristalsis

def register_processes(core):
    core.register_link("SimpleDiffusion", SimpleDiffusion)
    core.register_link("SimpleAdvection", SimpleAdvection)
    core.register_link("DynamicAdvection", DynamicAdvection)
    core.register_link("Peristalsis", Peristalsis)

    return core
