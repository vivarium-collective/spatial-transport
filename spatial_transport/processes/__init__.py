from spatial_transport.processes.diffusion import SimpleDiffusion

def register_processes(core):
    core.register_process("SimpleDiffusion", SimpleDiffusion)
    return core