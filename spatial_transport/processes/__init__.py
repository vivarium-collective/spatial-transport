import spatial_transport
from spatial_transport.processes.diffusion import SimpleDiffusion
from spatial_transport.processes.advection import SimpleAdvection, DynamicAdvection, Peristalsis

def register_processes(core):
    core.register_process("SimpleDiffusion", SimpleDiffusion)
    core.register_process("SimpleAdvection", SimpleAdvection)
    core.register_process("DynamicAdvection", DynamicAdvection)
    core.register_process("Peristalsis", Peristalsis)
    return core