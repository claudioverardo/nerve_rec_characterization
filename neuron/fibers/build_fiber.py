from neuron import h
from fibers.unmyelinated_nerve_fiber import UnmyelinatedNerveFiber
from fibers.myelinated_nerve_fiber import MyelinatedNerveFiber

def build_fiber(fiber_params):

    # Params for definition of fibers
    #
    # =======================================
    # UNMYELINATED FIBERS
    # =======================================
    # fiber_params = {
    #     'fiber_model': fiber_model,
    #     'diameter': diameter,
    #     'n_nodes': n_nodes,
    #     'node_length': node_length,
    #     'center_node_location': center_node_location
    # }
    #
    # =======================================
    # MYELINATED FIBERS
    # =======================================
    # fiber_params = {
    #     'fiber_model': fiber_model,
    #     'diameter': diameter,
    #     'n_internodes': n_internodes,
    #     'center_node_location': center_node_location
    # }
    
    fiber_model = fiber_params['fiber_model']
    if not(isinstance(fiber_model, str)):
        fiber_model = get_fiber_model_id(fiber_model)
    
    if not('center_node_location' in fiber_params):
        fiber_params['center_node_location'] = [0, 0, 0]

    if fiber_model == "ascent_sundt15":
        fiber = UnmyelinatedNerveFiber(
            diameter = fiber_params['diameter'],
            n_nodes = fiber_params['n_nodes'],
            node_length = fiber_params['node_length'],
            center_node_location = fiber_params['center_node_location'],
            fiber_model = fiber_model,
            passive_end_nodes = 0
        )

    elif fiber_model == "ascent_mrg":
        fiber = MyelinatedNerveFiber(
            diameter = fiber_params['diameter'],
            n_internodes = fiber_params['n_internodes'],
            center_node_location = fiber_params['center_node_location'],
            fiber_model = fiber_model,
            passive_end_nodes = 1
        )

    return fiber