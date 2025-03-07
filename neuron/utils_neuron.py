from neuron import h
import os
import platform

def init_neuron():
    """
    Load NMODL mechanisms in NEURON
    The current working directory is assumed to be nerve_rec_characterization/neuron
    """

    h.load_file("stdrun.hoc")

    # Windows
    if os.name == 'nt':

        # myelinated fibers
        h.nrn_load_dll("./fibers/mod_mrg02/nrnmech.dll")

        # unmyelinated fibers
        h.nrn_load_dll("./fibers/mod_sundt15/nrnmech.dll")

    # Linux/macOS
    elif os.name == 'posix':

        if platform.machine()=='x86_64':
            # myelinated fibers
            h.nrn_load_dll("./fibers/mod_mrg02/x86_64/.libs/libnrnmech.so")

            # unmyelinated fibers
            h.nrn_load_dll("./fibers/mod_sundt15/x86_64/.libs/libnrnmech.so")

        elif platform.machine()=='arm64':
            # myelinated fibers
            h.nrn_load_dll("./fibers/mod_mrg02/arm64/.libs/libnrnmech.so")

            # unmyelinated fibers
            h.nrn_load_dll("./fibers/mod_sundt15/arm64/.libs/libnrnmech.so")

        else:
            raise Exception("CPU architecture not supported")
    
    else:
        raise Exception("OS not supported")
