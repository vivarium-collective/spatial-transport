import spatial_transport
from spatial_transport.processes.diffusion import SimpleDiffusion
from spatial_transport.processes.advection import SimpleAdvection

def register_processes(core):
    core.register_process("SimpleDiffusion", SimpleDiffusion)
    core.register_process("SimpleAdvection", SimpleAdvection)
    return core