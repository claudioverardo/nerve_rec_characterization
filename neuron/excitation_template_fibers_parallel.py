from neuron import h
from fibers.build_fiber import build_fiber
from fibers.utils_fiber import check_fiber_myelinated
from utils_neuron import init_neuron

from multiprocessing import Pool, Value, Lock
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os
import sys

def run_excitation_template_fibers(template_type, fiber_params, stim_amp, stim_delay, stim_dur, stim_tstop, dt, dt_export, export_dir, n_internodes_CV=[100,120], min_height_CV=0):
    """
    Generate the excitation template for a fiber.
    A simulation of intracellular stimulation in NEURON is used. 

    NOTE: n_internodes_CV defines the # of the two internodes used to compute the 
          conduction velocity (CV). If template_type == "internode-CV", the variables 
          at the first of these internodes are stored in the exported template.
    
    Export variables (data): template_fiber_d_(#diameter).mat
    ===================================================================================
    Variable                      Units       Data structure
    ===================================================================================
    fiber_t_rec                   ms          1 x n_time_samples
    fiber_v_rec                   mV          fiber_n_compartments x n_time_samples
    fiber_itot_rec                mA/cm^2     fiber_n_compartments x n_time_samples
    ===================================================================================
    
    Export variables (metadata): template_fiber_d_(#diameter)_meta.mat
    ===================================================================================
    Variable                      Units       Data structure
    ===================================================================================
    fiber_template_type           -           string ("internode-CV" or "full-fiber")
    fiber_model                   -           string      
    fiber_diameter                um           
    fiber_length                  mm          
    fiber_morpho                              struct (for myelinated fibers)
        n_internodes              -     
        diam_ranvier              um          
        diam_mysa                 um          
        diam_flut                 um          
        diam_stin                 um          
        length_ranvier            um          
        length_mysa               um          
        length_flut               um            
        length_stin               um        
        length_internode          um        
        n_lamellae                um          
    fiber_morpho                              struct (for unmyelinated fibers)
        n_nodes                   -
        node_length               um     
    fiber_n_compartments (*)      -           
    fiber_area                    um^2        1 x fiber_n_compartments
    fiber_pos_z                   mm          1 x fiber_n_compartments
    fiber_dt                      ms
    fiber_stim_amp                nA
    fiber_stim_delay              ms
    fiber_stim_dur                ms
    fiber_stim_tstop              ms
    fiber_CV                      m/s
    fiber_n_internodes_CV         -           1 x 2
    fiber_n_nodes_CV              -           1 x 2

    (*) if template_type == "internode-CV":
            fiber_n_compartments = 11 (myelinated fiber)
            fiber_n_compartments = 1 (unmyelinated fiber)
        if template_type == "full-fiber":
            fiber_n_compartments = n_internodes (myelinated fiber)
            fiber_n_compartments = n_nodes (unmyelinated fiber)
    ===================================================================================
    """

    # build fiber model
    fiber_model = fiber_params['fiber_model']
    fiber = build_fiber(fiber_params)

    is_fiber_myelinated = check_fiber_myelinated(fiber_model)

    # run simulation of intracellular stimulation
    h.dt = dt
    v_rec, itot_rec, t_rec = fiber.apply_intrastim(
        stim_node=0, # -1
        rec_nodes=np.arange(fiber.n_nodes), 
        amp=stim_amp,      # [nA] # 0.5 nA * 0.5 durata attiva fino a diametro=12 um sicuro
        delay=stim_delay,  # [ms] 
        dur=stim_dur,      # [ms]
        tstop=stim_tstop   # [ms]
    )
    v_rec = np.array(v_rec)
    itot_rec = np.array(itot_rec)
    t_rec = np.array(t_rec)

    # ===================================== DEBUGGING PLOTS =====================================
    # fig = plt.figure()
    # plt.plot(t_rec, v_rec[round(fiber.n_nodes/2),:])
    # plt.show()
    
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # for i in np.arange(0,fiber.n_nodes,11):
    #     ax.plot3D(t_rec, np.ones(len(t_rec))*fiber.nodes_location[i,2], v_rec[i,:], 'black')
    # plt.show()

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # for i in np.arange(0,fiber.n_nodes,11):
    #     ax.plot3D(t_rec, np.ones(len(t_rec))*fiber.nodes_location[i,2], itot_rec[i,:], 'black')
    # plt.show()
    # ===========================================================================================

    # collect results data/meta
    # fiber_area = np.array([node(0.5).area() for node in fiber.nodes])   # 2*pi*r*Δx, where r is the axonic radius and Δx is the segment length
    fiber_area = np.array([np.pi*node.diam*node.L for node in fiber.nodes])   
    fiber_pos_z = fiber.nodes_location[:,2]

    if is_fiber_myelinated:
        n_nodes_CV = [(n_internodes_CV[0]-1)*11, (n_internodes_CV[1]-1)*11]
    else:
        n_nodes_CV = n_internodes_CV

    if template_type == "internode-CV":
        if is_fiber_myelinated:
            fiber_n_compartments = 11
            idcs_nodes_template = n_nodes_CV[0] + np.arange(0,11)
        else:
            fiber_n_compartments = 1
            idcs_nodes_template = np.array(n_nodes_CV[0])

    elif template_type == "full-fiber":
        fiber_n_compartments = fiber.n_nodes
        idcs_nodes_template = np.arange(0,v_rec.shape[0])
        
    else:
        raise Exception("invalid template type")

    # compute conduction velocity (CV)
    peak1_z = fiber_pos_z[n_nodes_CV[0]]
    peak1_t, _ = scipy.signal.find_peaks(v_rec[n_nodes_CV[0],:], height=min_height_CV)
    peak1_t = peak1_t[0]
    peak1_t = peak1_t*dt
    peak2_z = fiber_pos_z[n_nodes_CV[1]]
    peak2_t, _ = scipy.signal.find_peaks(v_rec[n_nodes_CV[1],:], height=min_height_CV)
    peak2_t = peak2_t[0]
    peak2_t = peak2_t*dt
    fiber_CV = (peak2_z-peak1_z) / (peak2_t-peak1_t)

    # keep only nodes to be exported + resampling at the export timestep
    fiber_area = fiber_area[idcs_nodes_template]
    fiber_pos_z = fiber_pos_z[idcs_nodes_template]

    v_rec = v_rec[idcs_nodes_template,0:len(t_rec):round(dt_export/h.dt)]
    itot_rec = itot_rec[idcs_nodes_template,0:len(t_rec):round(dt_export/h.dt)]
    t_rec = t_rec[0:len(t_rec):round(dt_export/h.dt)]

    # create export directories
    if not(os.path.isdir(export_dir)):
        os.mkdir(export_dir)

    export_dir_fiber_model = export_dir+'/'+fiber_model
    if not(os.path.isdir(export_dir_fiber_model)):
        os.mkdir(export_dir_fiber_model)

    # store results
    if is_fiber_myelinated:
        fiber_morpho = {
            'n_internodes': fiber.n_internodes, \
            'diam_ranvier': fiber.diam_ranvier, \
            'diam_mysa': fiber.diam_mysa, \
            'diam_flut': fiber.diam_flut, \
            'diam_stin': fiber.diam_stin, \
            'length_ranvier': fiber.length_ranvier, \
            'length_mysa': fiber.length_mysa, \
            'length_flut': fiber.length_flut, \
            'length_stin': fiber.length_stin, \
            'length_internode': fiber.length_internode, \
            'n_lamellae': fiber.n_lamellae
        }
    else:
        fiber_morpho = {
            'n_nodes': fiber.n_nodes, \
            'node_length': fiber.node_length
        }

    scipy.io.savemat(export_dir_fiber_model+'/template_fiber_d_' + "{:.2f}".format(fiber.diameter) + '_meta.mat', \
    {   
        'fiber_template_type': template_type,
        'fiber_model': fiber_model,
        'fiber_diameter': fiber.diameter,
        'fiber_length': fiber.nodes_location[-1,2]-fiber.nodes_location[0,2],
        'fiber_morpho': fiber_morpho,
        'fiber_n_compartments': fiber_n_compartments,
        'fiber_area': fiber_area,
        'fiber_pos_z' : fiber_pos_z,
        'fiber_dt': dt_export,
        'fiber_stim_amp': stim_amp,
        'fiber_stim_delay': stim_delay,
        'fiber_stim_dur': stim_dur,
        'fiber_stim_tstop': stim_tstop,
        'fiber_CV': fiber_CV,
        'fiber_n_internodes_CV': n_internodes_CV,
        'fiber_n_nodes_CV': n_nodes_CV
    })

    scipy.io.savemat(export_dir_fiber_model+'/template_fiber_d_' + "{:.2f}".format(fiber.diameter) + '.mat', \
    {   
        'fiber_t_rec': t_rec,
        'fiber_v_rec': v_rec,
        'fiber_itot_rec': itot_rec
    })

    fiber.remove_from_neuron()

def run_excitation_template_fibers_parallel(template_type, fiber_params, stim_amp, stim_delay, stim_dur, stim_tstop, dt, dt_export, export_dir, n_internodes_CV, min_height_CV):

    run_excitation_template_fibers(
        template_type = template_type,
        fiber_params = fiber_params,
        stim_amp = stim_amp,
        stim_delay = stim_delay,
        stim_dur = stim_dur,
        stim_tstop = stim_tstop,
        dt = dt,
        dt_export = dt_export,
        export_dir = export_dir,
        n_internodes_CV = n_internodes_CV,
        min_height_CV = min_height_CV)
    
    with lock_ref:

        print(fiber_params['fiber_model'] + ' d = ' + "{:.2f}".format(fiber_params['diameter']) + ' done')
        counter_ref.value += 1

    return

def init(counter, lock):
    global counter_ref
    counter_ref = counter

    global lock_ref
    lock_ref = lock

    init_neuron()


if __name__ == "__main__":
    # Example:
    # python excitation_template_fibers_parallel.py internode-CV ".\template_fibers" 8

    template_type = sys.argv[1]     # 'internode-CV' or 'full-fiber'
    export_dir = sys.argv[2]        # folder where to export the fiber excitation templates
    n_processes = int(sys.argv[3])  # number of CPU cores to use
                                    # WARNING: must be lower than the total number of available cores

    center_node_location = [0,0,0]

    arguments = []

    # ========================== myelinated fibers ========================== 
    n_internodes_CV = [100, 120]
    dt = 0.0025
    dt_export = dt 
    stim_amp_myel = 0.5   # [nA]
    stim_delay_myel = 0.1 # [ms]
    stim_dur_myel = 1     # [ms]
    stim_tstop_myel = 15  # [ms]
    n_internodes = 360    # IMPORTANT! EVEN NUMBER TO HAVE A RANVIER NODE @ z=0
    min_height_CV = -10

    fiber_model = "ascent_mrg"
    diameters = np.arange(1,16.5,0.5)
    arguments += [
        (template_type, {'fiber_model':fiber_model, 'diameter':diameter, 'n_internodes':n_internodes}, stim_amp_myel, stim_delay_myel, stim_dur_myel, stim_tstop_myel, dt, dt_export, export_dir, n_internodes_CV, min_height_CV)
        for diameter in diameters
    ]

    # ========================= unmyelinated fibers =========================
    n_internodes_CV = [200, 400]
    dt = 0.01
    dt_export = dt
    stim_amp_unmyel = 0.5   # [nA]
    stim_delay_unmyel = 0.1 # [ms]
    stim_dur_unmyel = 1     # [ms]
    stim_tstop_unmyel = 60  # [ms]
    n_nodes = 1251          # IMPORTANT! ODD NUMBER TO HAVE A NODE @ z=0
    node_length = 50/6      # [um] (= 8.3333 um default in ASCENT)
    min_height_CV = -10

    fiber_model = "ascent_sundt15"
    diameters = np.arange(0.1,2.01,0.1)
    arguments += [
        (template_type, {'fiber_model':fiber_model, 'diameter':diameter, 'n_nodes':n_nodes, 'node_length':node_length}, stim_amp_unmyel, stim_delay_unmyel, stim_dur_unmyel, stim_tstop_unmyel, dt, dt_export, export_dir, n_internodes_CV, min_height_CV)
        for diameter in diameters
    ]
    
    # =======================================================================
    
    counter = Value("i", 0)
    lock = Lock()

    with Pool(processes=n_processes, initializer=init, initargs=(counter, lock)) as pool:
        pool.starmap(func=run_excitation_template_fibers_parallel, iterable=arguments)
