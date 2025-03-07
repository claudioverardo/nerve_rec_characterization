from neuron import h
import numpy as np
from fibers.nerve_fiber import NerveFiber
from fibers.utils_fiber import get_fiber_model_id

class MyelinatedNerveFiber(NerveFiber):
    """
    Class defining a biophysical model of a myelinated nerve fiber (MRG model).
    Python implementation based on the HOC code available in ASCENT.
    Refs: McIntyre, Richardson, and Grill, J Neurophysiol, 2002
          McIntyre, Grill, Sherman, and Thakor, J Neurophysiol, 2004
    """

    def __init__(self, diameter, n_internodes, center_node_location = [0,0,0], fiber_model="ascent_mrg", passive_end_nodes=1):
        
        self.fiber_model = get_fiber_model_id(fiber_model)
        self.passive_end_nodes = passive_end_nodes

        # morphological properties
        self.length_mysa = 3.0
        self.length_ranvier = 1.0
        self.periwidth_mysa = 0.002
        self.periwidth_flut = 0.004
        self.periwidth_stin = 0.004

        # electrical properties
        self.cm = 2  # [uF/cm2]
        self.rhoa = 0.7e6  # [Ohm-cm]
        self.cm_my = 0.1  # [uF/cm2-lamella]
        self.gm_my = 0.001  # [S/cm2-lamella]

        # temperature
        self.celsius = 37 # [°C]
        h.celsius = self.celsius  

        # equilibrium potential
        self.v_init = -80  # [mV]

        # properties depending from self.n_internodes
        self.n_internodes = n_internodes
        self.n_ranvier = n_internodes + 1
        self.n_mysa = n_internodes * 2
        self.n_flut = n_internodes * 2
        self.n_stin = n_internodes * 6

        self.diameter = diameter

        # properties depending from self.diameter
        # (interpolation based on MRG02-04 data of myelinated fibers)
        #
        # ----------------------------------------------------------------------
        # ASCENT nomencl.           OUR nomencl.            note
        # ----------------------------------------------------------------------
        # fiberD                    => diameter             fixed (see above)
        # nodeD                     => diam_ranvier
        # paraD1 / MYSAD            => diam_mysa
        # paraD2 / FLUTD            => diam_flut
        # axonD                     => diam_stin
        # deltaz                    => length_internode
        # nodelength                => length_ranvier       fixed (see above)
        # paralength1 / MYSAlength  => length_mysa          fixed (see above)
        # paralength2 / FLUTlength  => length_flut
        # interlength               => length_stin
        # nl                        => n_lamellae
        
        self.diam_ranvier = 0.01093*diameter**2 + 0.1008*diameter + 1.099
        self.diam_mysa = self.diam_ranvier
        self.diam_flut = 0.02361*diameter**2 + 0.3673*diameter + 0.7122
        self.diam_stin = self.diam_flut
        if diameter >= 5.643:
            self.length_internode = -8.215*diameter**2 + 272.4*diameter - 780.2
        else:
            self.length_internode = 81.08*diameter + 37.84
        self.length_flut = -0.1652*diameter**2 + 6.354*diameter - 0.2862
        self.n_lamellae = -0.4749*diameter**2 + 16.85*diameter - 0.7648

        self.length_stin = (self.length_internode - self.length_ranvier - (2 * self.length_mysa) - (2 * self.length_flut)) / 6

        # periaxonal resistivity computation
        area_ext_ranvier = np.pi * (self.diam_ranvier / 2 + self.periwidth_mysa)**2
        area_int_ranvier = np.pi * (self.diam_ranvier / 2)**2
        periarea_ranvier = area_ext_ranvier - area_int_ranvier

        area_ext_mysa = np.pi * (self.diam_mysa / 2 + self.periwidth_mysa)**2
        area_int_mysa = np.pi * (self.diam_mysa / 2)**2
        periarea_mysa = area_ext_mysa - area_int_mysa

        area_ext_flut = np.pi * (self.diam_flut / 2 + self.periwidth_flut)**2
        area_int_flut = np.pi * (self.diam_flut / 2)**2
        periarea_flut = area_ext_flut - area_int_flut

        area_ext_stin = np.pi * (self.diam_stin / 2 + self.periwidth_stin)**2
        area_int_stin = np.pi * (self.diam_stin / 2)**2
        periarea_stin = area_ext_stin - area_int_stin

        self.R_pa_ranvier = (self.rhoa * 0.01) / periarea_ranvier
        self.R_pa_mysa = (self.rhoa * 0.01) / periarea_mysa
        self.R_pa_flut = (self.rhoa * 0.01) / periarea_flut
        self.R_pa_stin = (self.rhoa * 0.01) / periarea_stin

        # generate sections
        self.Ranvier = [h.Section() for i in range(self.n_ranvier)]
        self.MYSA = [h.Section() for i in range(self.n_mysa)]
        self.FLUT = [h.Section() for i in range(self.n_flut)]
        self.STIN = [h.Section() for i in range(self.n_stin)]

        # set biophysics
        for sec in h.allsec():
            sec.insert('extracellular')
            sec.e_extracellular = 0

        for i, ranvier in enumerate(self.Ranvier):
            ranvier.nseg = 1
            ranvier.diam = self.diam_ranvier
            ranvier.L = self.length_ranvier
            ranvier.Ra = self.rhoa / 1e4
            ranvier.cm = self.cm
            if (i == 0 or i == self.n_ranvier-1) and self.passive_end_nodes == 1:
                # passive end nodes to reduce edge effets
                ranvier.insert('pas')
                ranvier.g_pas = 0.0001
                ranvier.e_pas = -70
                ranvier.xg[0] = self.gm_my / (2 * self.n_lamellae)
                ranvier.xc[0] = self.cm_my / (2 * self.n_lamellae)
            else:
                ranvier.insert('axnode_myel')
                ranvier.xraxial[0] = self.R_pa_ranvier
                ranvier.xg[0] = 1e10
                ranvier.xc[0] = 0

        for mysa in self.MYSA:
            eta_diam = self.diam_mysa / self.diameter
            mysa.nseg = 1
            mysa.diam = self.diameter
            mysa.L = self.length_mysa
            mysa.Ra = (self.rhoa / eta_diam**2) / 1e4
            mysa.cm = self.cm * eta_diam
            mysa.insert('pas')
            mysa.g_pas = 0.001*eta_diam
            mysa.e_pas = self.v_init 
            mysa.xraxial[0] = self.R_pa_mysa
            mysa.xg[0] = self.gm_my / (2 * self.n_lamellae)
            mysa.xc[0] = self.cm_my / (2 * self.n_lamellae)

        for flut in self.FLUT:
            eta_diam = self.diam_flut / self.diameter
            flut.nseg = 1
            flut.diam = self.diameter
            flut.L = self.length_flut
            flut.Ra = (self.rhoa / eta_diam**2) / 1e4
            flut.cm = self.cm * eta_diam
            flut.insert('pas')
            flut.g_pas = 0.0001*eta_diam
            flut.e_pas = self.v_init 
            flut.xraxial[0] = self.R_pa_flut
            flut.xg[0] = self.gm_my / (2 * self.n_lamellae)
            flut.xc[0] = self.cm_my / (2 * self.n_lamellae)

        for stin in self.STIN:
            eta_diam = self.diam_stin / self.diameter
            stin.nseg = 1
            stin.diam = self.diameter
            stin.L = self.length_stin
            stin.Ra = (self.rhoa / eta_diam**2) / 1e4
            stin.cm = self.cm * eta_diam
            flut.insert('pas')
            flut.g_pas = 0.0001*eta_diam
            flut.e_pas = self.v_init 
            stin.xraxial[0] = self.R_pa_stin
            stin.xg[0] = self.gm_my / (2 * self.n_lamellae)
            stin.xc[0] = self.cm_my / (2 * self.n_lamellae)

        # set connectivity
        for i in range(self.n_internodes):
            self.MYSA[2*i].connect(self.Ranvier[i](1), 0)
            self.FLUT[2*i].connect(self.MYSA[2*i](1), 0)

            self.STIN[6*i].connect(self.FLUT[2*i](1), 0)
            self.STIN[6*i + 1].connect(self.STIN[6*i](1), 0)
            self.STIN[6*i + 2].connect(self.STIN[6*i + 1](1), 0)
            self.STIN[6*i + 3].connect(self.STIN[6*i + 2](1), 0)
            self.STIN[6*i + 4].connect(self.STIN[6*i + 3](1), 0)
            self.STIN[6*i + 5].connect(self.STIN[6*i + 4](1), 0)

            self.FLUT[2*i + 1].connect(self.STIN[6*i + 5](1), 0)
            self.MYSA[2*i + 1].connect(self.FLUT[2*i + 1](1), 0)
            self.Ranvier[i+1].connect(self.MYSA[2*i + 1](1), 0)

        self.nodes = []
        for i in range(self.n_internodes):
            self.nodes.append(self.Ranvier[i])
            self.nodes.append(self.MYSA[2*i])
            self.nodes.append(self.FLUT[2*i])
            self.nodes.append(self.STIN[6*i])
            self.nodes.append(self.STIN[6*i + 1])
            self.nodes.append(self.STIN[6*i + 2])
            self.nodes.append(self.STIN[6*i + 3])
            self.nodes.append(self.STIN[6*i + 4])
            self.nodes.append(self.STIN[6*i + 5])
            self.nodes.append(self.FLUT[2*i + 1])
            self.nodes.append(self.MYSA[2*i + 1])
        self.nodes.append(self.Ranvier[self.n_internodes])
        self.n_nodes = len(self.nodes)

        # self.stimulation_nodes = list(range(11 * self.n_internodes + 1))
        # self.n_stimulation_nodes = len(self.stimulation_nodes)

        self.center_node_location = center_node_location
        node_spacing_internode = np.array(
            [
                (self.length_ranvier + self.length_mysa) / 2,
                (self.length_mysa + self.length_flut) / 2,
                (self.length_flut + self.length_stin) / 2,
                self.length_stin,
                self.length_stin,
                self.length_stin,
                self.length_stin,
                self.length_stin,
                (self.length_stin + self.length_flut) / 2,
                (self.length_flut + self.length_mysa) / 2,
                (self.length_mysa + self.length_ranvier) / 2
            ]
        )
        node_spacing_fiber = np.tile(node_spacing_internode, (1, self.n_internodes))
        nodes_z = np.cumsum(node_spacing_fiber)
        nodes_z = np.concatenate(([0], nodes_z))
        nodes_z -= self.length_internode * (self.n_internodes / 2)
        nodes_z = nodes_z / 1000 + center_node_location[2]  # [mm]
        self.nodes_location = np.array(
            [
                [center_node_location[0], center_node_location[1], node_z] for node_z in nodes_z
            ]
        )  # [mm]

    def __repr__(self):
        line_1 = "Myelinated nerve fiber with the following features:"
        line_2 = "fiber diameter: " + str(self.diameter)
        line_3 = "number of internodes: " + str(self.n_internodes)
        line_4 = "center node location: " + str(self.center_node_location)

        return line_1 + "\n" + line_2 + "\n" + line_3 + "\n" + line_4