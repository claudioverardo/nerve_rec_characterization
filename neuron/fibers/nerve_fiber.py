from neuron import h
import numpy as np
from scipy import signal

class NerveFiber:
    """
    Base class to define biophysical models of nerve fibers.
    It includes methods to simulate intracellular stimulations.
    Subclasses may be derived to instantiate desired active mechanisms.

    Subclasses should implement the following properties:
        self.diameter           float                   [um]
        self.nodes              (n_nodes,) h.Section()  
        self.n_nodes            int
        self.celsius            float                   [°C]
        self.v_init             float                   [mV]
        self.nodes_location     (n_nodes, 3) float      [mm]
    """

    def __init__(self, diameter, n_nodes, node_length, center_node_location = [0,0,0]):
        """
        Create a passive cable model of a nerve fiber.
        """

        # morphology
        self.diameter = diameter
        self.n_nodes = n_nodes
        self.node_length = node_length  # [um]
        self.length = self.n_nodes*self.node_length/1000  # [mm]

        # temperature
        self.celsius = 37 # [°C]
        h.celsius = self.celsius

        # equilibrium potential
        self.v_init = -70
        
        # set biophysics
        self.nodes = []
        for i in range(self.n_nodes):
            node = h.Section()
            node.nseg = 1
            node.diam = self.diameter
            node.L = self.node_length
                
            node.insert('pas')
            node.g_pas = 0.0001
            node.e_pas = -70

            node.insert('extracellular')
            node.xg[0] = pow(10,10)
            node.xc[0] = 0

            node.Ra = 100
            node.cm = 1
                
            self.nodes.append(node)

        # self.stimulation_nodes = list(range(self.n_nodes))
        # self.n_stimulation_nodes = self.n_nodes

        # connect sections
        for i in range(self.n_nodes - 1):
            self.nodes[i].connect(self.nodes[i+1])

        # fiber geometry
        self.center_node_location = center_node_location
        nodes_spacing = [node.L for node in self.nodes[0:-1]] 
        nodes_z = np.cumsum(nodes_spacing)/1000 # [mm]
        nodes_z = np.concatenate(([0], nodes_z))
        nodes_z -= (self.length-self.node_length/1000)/2
        nodes_z = nodes_z + center_node_location[2] 
        self.nodes_location = np.array(
            [
                [center_node_location[0], center_node_location[1], node_z] for node_z in nodes_z
            ]
        )  # [mm]

    def __repr__(self):
        line_1 = "Passive nerve fiber with the following features:"
        line_2 = "fiber diameter: " + str(self.diameter)
        line_3 = "number of nodes: " + str(self.n_nodes)
        line_4 = "center node location: " + str(self.center_node_location)

        return line_1 + "\n" + line_2 + "\n" + line_3 + "\n" + line_4

    def init_simulation(self):
        h.finitialize(self.v_init)

    def apply_intrastim(self, stim_node, rec_nodes, amp, delay, dur, tstop):
        """
        Apply intracellular stimulation to the fiber.

        Parameters
        ----------
        stim_node : int
            Stimulation node identifier.
        rec_nodes : (n_recording_nodes,) list of int
            Recording node identifiers.
        amp : float
            Stimulating current amplitude [nA].
        delay : float
            Stimulation start time [ms].
        dur : float
            Stimulation duration [ms].
        tstop : float
            Duration of simulation [ms] (time of start simulation : 0 ms).

        Returns
        -------
        v_rec : (n_recording_nodes,) list of (n_time_samples,) list of float
            Membrane potential (v) time-courses at the recording nodes [mV].
        itot_rec : (n_recording_nodes,) list of (n_time_samples,) list of float
            Transmembrane current (itot) time-courses at the recording nodes [mA/cm^2].
        t_rec : (n_time_samples,) list of float
            Simulation time [ms].
        """

        stim = h.IClamp(self.nodes[stim_node](0.5))
        stim.amp = amp
        stim.delay = delay
        stim.dur = dur

        n_rec = len(rec_nodes)
        v_rec = [None for i in range(n_rec)]
        itot_rec = [None for i in range(n_rec)]
        for i in range(n_rec):
            v_rec[i] = h.Vector().record(self.nodes[rec_nodes[i]](0.5)._ref_v)
            itot_rec[i] = h.Vector().record(self.nodes[rec_nodes[i]](0.5)._ref_i_membrane)
        t_rec = h.Vector().record(h._ref_t)

        # run simulation
        self.init_simulation()
        h.continuerun(tstop)

        for i in range(n_rec):
            v_rec[i] = v_rec[i].to_python()
            itot_rec[i] = itot_rec[i].to_python()
        t_rec = t_rec.to_python()

        return v_rec, itot_rec, t_rec

    def remove_from_neuron(self):
        for k, v in self.__dict__.items():
            del v

    @staticmethod
    def compute_activation(v_rec, height=0, prominence=60):
        """
        Compute activation from a recorded membrane potential.
        """
        n_spikes = signal.find_peaks(v_rec, height=height, prominence=prominence)[0].size

        activation = n_spikes > 0
        return activation