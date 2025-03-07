from neuron import h
import numpy as np
import scipy as sp
from scipy import interpolate
from fibers.nerve_fiber import NerveFiber
from fibers.utils_fiber import get_fiber_model_id

class UnmyelinatedNerveFiber(NerveFiber):
    """
    Class defining a biophysical model of an unmyelinated nerve fiber (Sundt model).
    Python implementation based on the HOC code available in ASCENT.
    Ref: Sundt et al. J Neurophysiol (2015)
    """

    def __init__(self, diameter, n_nodes, node_length, fiber_model = "ascent_sundt15", passive_end_nodes = 0, center_node_location = [0,0,0]):
        
        self.fiber_model = get_fiber_model_id(fiber_model)
        
        # morphology
        self.diameter = diameter
        self.n_nodes = n_nodes
        self.node_length = node_length  # [um]
        self.length = self.n_nodes*self.node_length/1000  # [mm]
        
        # equilibrium potential
        # NB: moved here the original node.v variable
        if(self.fiber_model == get_fiber_model_id("ascent_sundt15")):
            self.v_init = -60
        else:
            self.v_init = -70

        self.passive_end_nodes = passive_end_nodes

        # temperature
        self.celsius = 37 # [°C]
        h.celsius = self.celsius
    
        # set biophysics
        self.nodes = []
        for i in range(self.n_nodes):
            node = h.Section()
            node.nseg = 1
            node.diam = self.diameter
            node.L = self.node_length

            if (i==0 or i==self.n_nodes-1) and self.passive_end_nodes == 1:
                node.insert('pas')
                node.g_pas = 0.0001
                
                # equilibrium potential
                if(self.fiber_model == get_fiber_model_id("ascent_sundt15")):
                    node.e_pas = -60
                else:
                    node.e_pas = -70
                
                node.insert('extracellular')
                node.xg[0] = pow(10,10) #short circuit, no myelin
                node.xc[0] = 0          #short circuit, no myelin

                node.Ra = pow(10,10)
                node.cm = 1

            else:

                # Sundt15 model
                # NOTE: temperature set by h.celsius
                if(self.fiber_model == get_fiber_model_id("ascent_sundt15")):
                    node.insert("nahh")
                    node.gnabar_nahh = 0.04
                    node.mshift_nahh = -6
                    node.hshift_nahh = 6

                    node.insert("borgkdr")
                    node.gkdrbar_borgkdr = 0.04
                    node.ek = -90

                    node.insert("pas")
                    node.g_pas = 1/10000
                    # node.v = -60
                    # node.e_pas = node.v + (node.ina + node.ik)/node.g_pas
                    node.e_pas = self.v_init + (node.ina + node.ik)/node.g_pas

                    node.Ra = 100
                    node.cm = 1
            
                node.insert("extracellular")
                node.xg[0] = pow(10,10)
                node.xc[0] = 0
                
            self.nodes.append(node)

        # self.stimulation_nodes = list(range(self.n_nodes))
        # self.n_stimulation_nodes = self.n_nodes

        # connect sections
        for i in range(self.n_nodes - 1):
            self.nodes[i].connect(self.nodes[i+1])
        
        # fiber geometry (straight fiber)
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
        line_1 = "Unyelinated nerve fiber with the following features:"
        line_2 = "fiber model: " + get_fiber_model_id(self.fiber_model)
        line_3 = "fiber diameter: " + str(self.diameter)
        line_4 = "number of nodes: " + str(self.n_nodes)
        line_5 = "center node location: " + str(self.center_node_location)

        return line_1 + "\n" + line_2 + "\n" + line_3 + "\n" + line_4 + "\n" + line_5