import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import time
import scipy.signal
import matplotlib
from IPython import display
import time
import os

class Simulation(object):

    def __init__(self, fs=3600.0, f0=50, z=5, res=2e-4, d=11e-3):
        """
        Initialize the main simulation parameters

        :param fs: Sampling frequency [Samples/Sec] (default=3600)
        :param f0: Base frequency [Hz] (default=50)
        :param z: Numper of periods to simulate [-] (default=5)
        :param res: Resolution of the mesh-grid [m] (default=2e-4)
        :param d: Outer diameter of the whole conductor [m] (default=11e-3)
        """
        assert (d <= 0.1), 'Error: Diameteter must be smaller than 10cm'
        self.created = datetime.datetime.now().isoformat()
        self.fs = fs
        self.f0 = f0
        self.z = z
        self.T = 1 / self.f0  # Period [s]
        self.N = self.fs * self.T * self.z  # Number of samples per period [-]
        self.w0 = 2 * np.pi * self.f0  # Angular frequency  [rad/s]
        self.t = np.arange(0, self.z * self.T, (1 / self.fs))  # Time vector [s]
        self.res = res
        self.d = d

        # define min and max of the grid
        self.x_min = -self.d
        self.x_max = self.d
        self.y_min = -self.d
        self.y_max = self.d

        self.__init_mesh()
        self.conductors = {}
        self.sensors = {}

        self.B_meshgrids = []

    def __init_mesh(self):
        """
        Initialize the mesh-grid

        :return:
        """
        x = np.arange(self.x_min, self.x_max, self.res)
        y = np.arange(self.y_min, self.y_max, self.res)
        x, y = np.meshgrid(x, y)
        self.mesh_x = x
        self.mesh_y = y

    def __add_conductor(self, name, pos_x, pos_y, diam):
        """
        Add a conductor / wire tho the simulation

        :param name: Unique name, usually L1 to L3 or N, PE
        :param pos_x: Position on the x-axis
        :param pos_y: Position on the y-axis
        :param diam: Diameter of the conductor / wire
        :return:
        """
        if (pos_x < self.x_min) or (pos_x > self.x_max):
            raise ValueError('Conductor outside the cable!')
        if (pos_y < self.y_min) or (pos_y > self.y_max):
            raise ValueError('Conductor outside the cable!')

        self.conductors[name] = self.Conductor(name, pos_x, pos_y, diam)

    def add_L1(self, pos_x, pos_y, diam):
        """
        User function to add L1 to the simulation

        :param pos_x: Position on the x-axis
        :param pos_y: Position on the y-axis
        :param diam: Diameter of the conductor / wire
        :return:
        """
        self.__add_conductor('L1', pos_x, pos_y, diam)

        # add default parameters
        self._set_current('L1', A1=1, A3=0, A5=0, A7=0, Phi1=0, Phi3=0, Phi5=0, Phi7=0, Shift=0)
        self._set_voltage('L1', A1=230, A3=0, A5=0, A7=0, Phi1=0, Phi3=0, Phi5=0, Phi7=0, Shift=0)

    def add_L2(self, pos_x, pos_y, diam):
        """
        User function to add L2 to the simulation

        :param pos_x: Position on the x-axis
        :param pos_y: Position on the y-axis
        :param diam: Diameter of the conductor / wire
        :return:
        """
        self.__add_conductor('L2', pos_x, pos_y, diam)

        # add default parameters
        self._set_current('L2', A1=1, A3=0, A5=0, A7=0, Phi1=0, Phi3=0, Phi5=0, Phi7=0, Shift=-120)
        self._set_voltage('L2', A1=230, A3=0, A5=0, A7=0, Phi1=0, Phi3=0, Phi5=0, Phi7=0, Shift=-120)

    def add_L3(self, pos_x, pos_y, diam):
        """
        User function to add L3 to the simulation

        :param pos_x: Position on the x-axis
        :param pos_y: Position on the y-axis
        :param diam: Diameter of the conductor / wire
        :return:
        """
        self.__add_conductor('L3', pos_x, pos_y, diam)

        # add default parameters
        self._set_current('L3', A1=1, A3=0, A5=0, A7=0, Phi1=0, Phi3=0, Phi5=0, Phi7=0, Shift=120)
        self._set_voltage('L3', A1=230, A3=0, A5=0, A7=0, Phi1=0, Phi3=0, Phi5=0, Phi7=0, Shift=120)

    def add_N(self, pos_x, pos_y, diam):
        """
        User function to add N (neutral) to the simulation

        :param pos_x: Position on the x-axis
        :param pos_y: Position on the y-axis
        :param diam: Diameter of the conductor / wire
        :return:
        """
        self.__add_conductor('N', pos_x, pos_y, diam)

    def add_PE(self, pos_x, pos_y, diam):
        """
        User function to add PE (earth) to the simulation

        :param pos_x: Position on the x-axis
        :param pos_y: Position on the y-axis
        :param diam: Diameter of the conductor / wire
        :return:
        """
        self.__add_conductor('PE', pos_x, pos_y, diam)

    def add_sensor(self, name, pos_x, pos_y):
        """
        User function to add one sensor to the simulation

        :param name: Unique name
        :param pos_x: Position on the x-axis
        :param pos_y: Position on the y-axis
        :return:
        """
        self.sensors[name] = self.Sensor(name, pos_x, pos_y)

    def add_sensors(self, n, r=None, phi0=0):
        """
        User function to add multiple sensors to the simulation.
        The sensors will be distributed equal around the conductor
        at a distance of r and an initial offset angle of phi0
        referenced to the positive x-axis

        :param n: Number of sensors to place
        :param r: Optional. Distance (in m) from the center at which the sensors will be placed.
        If not provided, sensors will be placed with an offset of 1mm on the cable
        :param phi0: Optional. Initial angle-offset in degrees. Defaults to 0°
        :return:
        """
        self.sensors = {}
        # assert n > 1, 'Error: number of sensors must be greater than 1'

        if r is None:
            r = self.d/2 + 1e-3
        else:
            assert r > self.x_max, "Error, sensor outside the simulation space"

        angles = np.arange(0, 2 * np.pi, 2 * np.pi / n)
        angles = angles + np.deg2rad(phi0)
        for i, phi in enumerate(angles):
            z = (r) * np.exp(1j * (phi))
            self.add_sensor('m' + str(i + 1), np.real(z), np.imag(z))

    def set_L1_I(self, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift):
        """
        User function to set parameters of the current in L1

        :param A1: RMS-Amplitude (in A) of the base frequency (50Hz)
        :param A3: RMS-Amplitude (in A) of the third harmonic (150Hz)
        :param A5: RMS-Amplitude (in A) of the fifth harmonic (250Hz)
        :param A7: RMS-Amplitude (in A) of the seventh harmonic (350Hz)
        :param Phi1: Phase-Shift (in deg) of the base frequency (50Hz)
        :param Phi3: Phase-Shift (in deg) of the third harmonic (150Hz)
        :param Phi5: Phase-Shift (in deg) of the fifth harmonic (250Hz)
        :param Phi7: Phase-Shift (in deg) of the seventh harmonic (350Hz)
        :param Shift: Global-Shift (in deg) of all frequency components
        :return:
        """
        self._set_current('L1', A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift)

    def set_L2_I(self, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift):
        """
        User function to set parameters of the current in L2

        :param A1: RMS-Amplitude (in A) of the base frequency (50Hz)
        :param A3: RMS-Amplitude (in A) of the third harmonic (150Hz)
        :param A5: RMS-Amplitude (in A) of the fifth harmonic (250Hz)
        :param A7: RMS-Amplitude (in A) of the seventh harmonic (350Hz)
        :param Phi1: Phase-Shift (in deg) of the base frequency (50Hz)
        :param Phi3: Phase-Shift (in deg) of the third harmonic (150Hz)
        :param Phi5: Phase-Shift (in deg) of the fifth harmonic (250Hz)
        :param Phi7: Phase-Shift (in deg) of the seventh harmonic (350Hz)
        :param Shift: Global-Shift (in deg) of all frequency components
        :return:
        """
        self._set_current('L2', A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift)

    def set_L3_I(self, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift):
        """
        User function to set parameters of the current in L3

        :param A1: RMS-Amplitude (in A) of the base frequency (50Hz)
        :param A3: RMS-Amplitude (in A) of the third harmonic (150Hz)
        :param A5: RMS-Amplitude (in A) of the fifth harmonic (250Hz)
        :param A7: RMS-Amplitude (in A) of the seventh harmonic (350Hz)
        :param Phi1: Phase-Shift (in deg) of the base frequency (50Hz)
        :param Phi3: Phase-Shift (in deg) of the third harmonic (150Hz)
        :param Phi5: Phase-Shift (in deg) of the fifth harmonic (250Hz)
        :param Phi7: Phase-Shift (in deg) of the seventh harmonic (350Hz)
        :param Shift: Global-Shift (in deg) of all frequency components
        :return:
        """
        self._set_current('L3', A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift)

    def set_L1_V(self, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift):
        """
        User function to set parameters of the voltage in L1

        :param A1: RMS-Amplitude (in V) of the base frequency (50Hz)
        :param A3: RMS-Amplitude (in V) of the third harmonic (150Hz)
        :param A5: RMS-Amplitude (in V) of the fifth harmonic (250Hz)
        :param A7: RMS-Amplitude (in V) of the seventh harmonic (350Hz)
        :param Phi1: Phase-Shift (in deg) of the base frequency (50Hz)
        :param Phi3: Phase-Shift (in deg) of the third harmonic (150Hz)
        :param Phi5: Phase-Shift (in deg) of the fifth harmonic (250Hz)
        :param Phi7: Phase-Shift (in deg) of the seventh harmonic (350Hz)
        :param Shift: Global-Shift (in deg) of all frequency components
        :return:
        """
        self._set_voltage('L1', A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift)

    def set_L2_V(self, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift):
        """
        User function to set parameters of the voltage in L2

        :param A1: RMS-Amplitude (in V) of the base frequency (50Hz)
        :param A3: RMS-Amplitude (in V) of the third harmonic (150Hz)
        :param A5: RMS-Amplitude (in V) of the fifth harmonic (250Hz)
        :param A7: RMS-Amplitude (in V) of the seventh harmonic (350Hz)
        :param Phi1: Phase-Shift (in deg) of the base frequency (50Hz)
        :param Phi3: Phase-Shift (in deg) of the third harmonic (150Hz)
        :param Phi5: Phase-Shift (in deg) of the fifth harmonic (250Hz)
        :param Phi7: Phase-Shift (in deg) of the seventh harmonic (350Hz)
        :param Shift: Global-Shift (in deg) of all frequency components
        :return:
        """
        self._set_voltage('L2', A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift)

    def set_L3_V(self, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift):
        """
        User function to set parameters of the voltage in L3

        :param A1: RMS-Amplitude (in V) of the base frequency (50Hz)
        :param A3: RMS-Amplitude (in V) of the third harmonic (150Hz)
        :param A5: RMS-Amplitude (in V) of the fifth harmonic (250Hz)
        :param A7: RMS-Amplitude (in V) of the seventh harmonic (350Hz)
        :param Phi1: Phase-Shift (in deg) of the base frequency (50Hz)
        :param Phi3: Phase-Shift (in deg) of the third harmonic (150Hz)
        :param Phi5: Phase-Shift (in deg) of the fifth harmonic (250Hz)
        :param Phi7: Phase-Shift (in deg) of the seventh harmonic (350Hz)
        :param Shift: Global-Shift (in deg) of all frequency components
        :return:
        """
        self._set_voltage('L3', A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift)

    def _set_current(self, name, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift):
        """
        Function to set the current of an existing conductor

        :param name: Unique name of an existing conductor
        :param A1: RMS-Amplitude (in A) of the base frequency (50Hz)
        :param A3: RMS-Amplitude (in A) of the third harmonic (150Hz)
        :param A5: RMS-Amplitude (in A) of the fifth harmonic (250Hz)
        :param A7: RMS-Amplitude (in A) of the seventh harmonic (350Hz)
        :param Phi1: Phase-Shift (in deg) of the base frequency (50Hz)
        :param Phi3: Phase-Shift (in deg) of the third harmonic (150Hz)
        :param Phi5: Phase-Shift (in deg) of the fifth harmonic (250Hz)
        :param Phi7: Phase-Shift (in deg) of the seventh harmonic (350Hz)
        :param Shift: Global-Shift (in deg) of all frequency components
        :return:
        """
        if name in self.conductors.keys():
            sig = self.Signal(A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift, self.fs, self.f0, self.w0, self.t)
            self.conductors[name].set_current_sig(sig)
        else:
            raise ValueError('No conductor with this name!')

    def _set_voltage(self, name, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift):
        """
        Function to set the voltage of an existing conductor (referenced to N)

        :param name: Uniqua name of an existing conductor
        :param A1: RMS-Amplitude (in V) of the base frequency (50Hz)
        :param A3: RMS-Amplitude (in V) of the third harmonic (150Hz)
        :param A5: RMS-Amplitude (in V) of the fifth harmonic (250Hz)
        :param A7: RMS-Amplitude (in V) of the seventh harmonic (350Hz)
        :param Phi1: Phase-Shift (in deg) of the base frequency (50Hz)
        :param Phi3: Phase-Shift (in deg) of the third harmonic (150Hz)
        :param Phi5: Phase-Shift (in deg) of the fifth harmonic (250Hz)
        :param Phi7: Phase-Shift (in deg) of the seventh harmonic (350Hz)
        :param Shift: Global-Shift (in deg) of all frequency components
        :return:
        """
        if name in self.conductors.keys():
            sig = self.Signal(A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift, self.fs, self.f0, self.w0, self.t)
            self.conductors[name].set_voltage_sig(sig)
        else:
            raise ValueError('No conductor with this name!')

    def plot_Waveforms(self):
        """
        Function provides a plot of currents and voltages of the available conductors

        :return:
        """
        self.__calc_N()
        fig, ax = plt.subplots(nrows=2, sharex=True, figsize=(12, 6))

        for name, con in self.conductors.items():
            if not (con.V is None):
                ax[0].plot(self.t, con.V, label=name)
        ax[0].set_title('Line Voltages')
        ax[0].set_xlabel('Time [s]')
        ax[0].set_ylabel('Voltage [V]')

        for name, con in self.conductors.items():
            if not (con.I is None):
                ax[1].plot(self.t, con.I, label=name)
        ax[1].set_title('Line Currents')
        ax[1].set_xlabel('Time [s]')
        ax[1].set_ylabel('Current [A]')
        for axes in ax:
            axes.legend()
            axes.grid()
        fig.tight_layout()

    def __calc_N(self):
        """
        Function to calculate the resulting neutral current from available conductors except PE

        :return:
        """
        signals = []
        if 'N' in self.conductors.keys():
            # check if any other conductor is defined and has a signal
            for name, con in self.conductors.items():
                if (name == 'N') or (name == 'PE'):
                    continue
                signals.append(con.I)
            if len(signals) > 0:
                self.conductors['N'].I = -np.sum(signals, axis=0)
        else:
            raise Warning('Neutral is not defined!')

    def plot_Layout(self):
        """
        Function to plot the layout of the simulation including conductor position and diameter and optionally placed
        sensors

        :return:
        """
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.set_xlim([self.x_min, self.x_max])
        ax.set_ylim([self.y_min, self.y_max])
        ax.set_xlabel('X-Axis [m]')
        ax.set_ylabel('Y-Axis [m]')
        ax.add_artist(plt.Circle((0, 0), (self.d / 2), color='k', fill=False))
        if len(self.conductors) > 0:
            for name, con in self.conductors.items():
                ax.add_artist(plt.Circle((con.x, con.y), con.d / 2, color='k', fill=False))
                ax.text(con.x + 1e-3, con.y, con.name, fontsize=10)

        if len(self.sensors) > 0:
            for sensor in self.sensors.values():
                p_x = sensor.x
                p_y = sensor.y
                tmp_circle = plt.Circle((p_x, p_y), 0.2e-3, color='r', fill=True)
                ax.add_artist(tmp_circle)
                ax.text(p_x, p_y, sensor.name, fontsize=12)
        ax.set_aspect('equal')
        ax.grid(True)

        return fig, ax

    def plot_Fields(self, idx, scale='global'):
        """
        Function to plot individual magnetic field maps

        :param idx: Time-Index for which the corresponding magnetic field map should be plotted
        :param scale: Applied scaling factor of the meshgrid. "global" fixes the scale to the global maximum and minimum
        of the entire time series. "local" scales each plot individually. Default: 'global'
        :return:
        """
        fig, ax = self.plot_Layout()

        if len(self.B_meshgrids) > 0:
            mesh = self.B_meshgrids[idx]
            if scale == 'global':
                vmin = np.min(self.B_meshgrids)
                vmax = np.max(self.B_meshgrids)
            if scale == 'local':
                vmin = np.min(mesh)
                vmax = np.ax(mesh)
            im = ax.imshow(mesh,
                           aspect='equal',
                           cmap='jet',
                           norm=matplotlib.colors.LogNorm(),
                           origin='lower',
                           extent=(self.x_min, self.x_max, self.y_min, self.y_max),
                           vmin=vmin,
                           vmax=vmax)
            cb = plt.colorbar(im)
            cb.set_label('Magnetic Flux Density [uT]')

    #         for sensor in self.sensors.values():
    #             p_x = sensor.x
    #             p_y = sensor.y
    #             tmp_circle = plt.Circle((p_x, p_y), 0.5e-3, color='r', fill=True)
    #             ax.add_artist(tmp_circle)
    #             ax.text(p_x, p_y, sensor.name, fontsize=12)
    #         ax.grid(True)
    #         ax.set_aspect('equal')

    def get_Meshgrids(self):
        """
        Returns an array of all calculated meshgrids

        :return: Array of size N
        """
        return self.B_meshgrids

    def get_Measurements(self):
        """
        Returns a pandas DataFrame of the measured magnetic flux density measured for each sensor.
        Columns represent sensors, rows individual samples.

        :return: pandas DataFrame
        """
        tmp = {}
        for name, sensor in self.sensors.items():
            tmp[name] = sensor.B_fields

        return pd.DataFrame(tmp)
    
    def get_Measurements_Exact(self):
        
        B = {}
        sens = []
        for key in self.sensors.keys():
            B[key] = []
            sens.append(self.sensors[key].x)
            sens.append(self.sensors[key].y)
        sens = np.array(sens).reshape(-1, 2)

        for idx in range(len(self.t)):
            x = []
            for key in self.conductors.keys():
                if key == 'PE':
                    continue
                x.append([self.conductors[key].x, self.conductors[key].y, self.conductors[key].I[idx]])
            x = np.array(x).reshape(len(self.conductors)-1, -1)

            tmp = self._get_B_exact(x, sens)
            for i, key, in enumerate(self.sensors.keys()):
                B[key].append(tmp[i])

        return pd.DataFrame(B)


        

            

    def get_Waveforms(self):
        """
        Returns a pandas DataFrame of all currents and voltages

        :return: pandas DataFrame
        """
        self.__calc_N()
        tmp = {}
        tmp['Time'] = self.t

        for name, con in self.conductors.items():
            if not (con.V is None):
                tmp[name + 'V'] = con.V

        for name, con in self.conductors.items():
            if not (con.I is None):
                tmp[name + 'I'] = con.I

        return pd.DataFrame(tmp)

    def get_Features(self):
        """
        Returns a dictionary consisting of all relevant features of the IEEE1459 standard for each current carrying
        conductor
        :return: dict with conductor names as key and features as values
        """
        self.__calc_N()
        features = {}
        for name, con in self.conductors.items():
            if (name == 'N') or (name == 'PE'):
                continue
            features[name] = con.get_Power()
        return features

    def _get_B_exact(self, x, sens):
        """
        param x: location and current of conductors, [[x1, y1, I1], [x2, y2, I2], ...]
        param sens: location of sensors [[x1, y1], [x2, y2]]
        """
        Theta = 0
        u0 = 4 * np.pi * 1e-7
        N = sens.shape[0]
        M = x.shape[0]
        #     result = np.ones((M, N))*np.nan
        result = []

        for n in range(N):  # n=sensor number
            for m in range(M):  # m = conductor number
                den = 2 * np.pi * np.sqrt((sens[n, 0]) ** 2 + (sens[n, 1]) ** 2) * (
                            (sens[n, 0] - x[m, 0]) ** 2 + (sens[n, 1] - x[m, 1]) ** 2)
                nom = u0 * x[m, 2] * (
                            ((sens[n, 0] * np.cos(Theta) - sens[n, 1] * np.sin(Theta)) * (sens[n, 0] - x[m, 0])) + \
                            +((sens[n, 1] * np.cos(Theta) + sens[n, 0] * np.sin(Theta)) * (sens[n, 1] - x[m, 1])))
                #             result[m,n] = nom/den
                result.append(nom / den)
        result = np.array(result).reshape(N, M)
        result = np.array(result).sum(axis=1) * 1e6
        return result

    def _get_B(self, x, y, xi, yi, d, i):
        """
        Calculates the magnetic field strength at a certain location (x-, y-coordiantes)

        Parameters:
        x, y: meshgrid [Meter]
        xi, yi: location of the conductor [Meter]
        d: dimension of conductor in  [Meter]
        i: current through the conductor [Ampere]

        return:
        bx, by: meshgrid
        """

        mu = 4*np.pi*1e-7  # magnetic constant
        r = np.sqrt((x - xi) ** 2 + (y - yi) ** 2)  # distance to conductor
        #     r = r + 1e-4
        mask = (r) < d / 2  # use only coordinates outside of conductor
        r[mask] = d / 2  # set distance inside conductor to big values
        mag = (mu / (2 * np.pi)) * (abs(i) / r)  # Magnitude of the vector B

        # calculate x- and y-component of magnitude
        r = np.sqrt((x - xi) ** 2 + (y - yi) ** 2)  # distance to conductor
        x = (x - xi) / r
        y = (y - yi) / r
        if (i > 0):
            by = mag * x  # By
            bx = mag * -y  # Bx
        else:
            by = mag * -x  # By
            bx = mag * y  # Bx

        return bx, by, mag

    def run(self):
        """
        Start the simulation and calculate currents and magnetic fields
        :return:
        """
        self.StartTime = datetime.datetime.now().isoformat()
        self.__calc_N()
        # reset calculated meshgrids
        self.B_meshgrids = []
        # reset sensors
        for sensor in self.sensors.values():
            sensor.B_fields = []

        for idx in range(len(self.t)):
            bx = []
            by = []
            for name, con in self.conductors.items():
                if name == 'PE':
                    continue
                _bx, _by, _mag = self._get_B(self.mesh_x,
                                             self.mesh_y,
                                             con.x,
                                             con.y,
                                             con.d,
                                             con.I[idx])
                bx.append(_bx)
                by.append(_by)

            # calculate superposition of all magnetic field vectors
            bx = np.sum(bx, axis=0)
            by = np.sum(by, axis=0)

            # calculate magnitude and convert into uT
            mag = np.sqrt(bx ** 2 + by ** 2) * 1e6
            self.B_meshgrids.append(mag)
            zero = len(mag) / 2

            for sensor in self.sensors.values():
                # get sensor coordinates
                x_s = sensor.x
                y_s = sensor.y

                # get index relative to meshgrid
                x_s_idx = int(zero + x_s / self.res)
                y_s_idx = int(zero + y_s / self.res)

                # get B-Vector in coresponding direction
                x_B = bx[y_s_idx, x_s_idx]
                y_B = by[y_s_idx, x_s_idx]

                # calculate magnitude according to sensors angle in uT
                B_hat = ((-y_s * x_B + x_s * y_B) / (np.sqrt(x_s ** 2 + y_s ** 2)) * 1e6)

                sensor.B_fields.append(B_hat)
        self.EndTime = datetime.datetime.now().isoformat()

    class Conductor(object):

        def __init__(self, name, pos_x, pos_y, diam):

            self.name = name
            self.x = pos_x
            self.y = pos_y
            self.d = diam
            self.sig_I = None
            self.sig_V = None
            self._I = None

        def set_current_sig(self, sig):
            """
            Sets the current of this conductor
            :param sig: Signal object
            :return:
            """
            self.sig_I = sig

        def get_current_sig(self):
            """
            Returns current of this conductor
            :return: Signal object
            """
            return self.sig_I

        def set_voltage_sig(self, sig):
            """
            Sets the voltage of this conductor (referenced to N)

            :param sig: Signal object
            :return:
            """
            self.sig_V = sig

        def get_voltage_sig(self):
            """
            Gets the voltage of this conductor
            :return: Signal object
            """
            return self.sig_V

        def update(self, w0, t):
            """
            Function to update important simulation parameters

            :param w0: Angular frequency [rad/s]
            :param t: time vector [s]
            :return:
            """
            if not (self.sig_I is None) and (self.sig_V is None):
                self.sig_I.update(w0, t)
                self.sig_V.update(w0, t)

        def get_Power(self, verbose=False):
            """
            Returns a dict of relevant IEEE1459 features
            :param verbose: True to plot calculated features.
            :return: dict with feature names as keys and features as values
            """
            i = self.I
            v = self.V
            fs = self.sig_I.fs
            f0 = self.sig_I.f0
            N = len(i)
            Shift_90 = int((fs / f0) / 4)
            assert len(i) == len(v)
            # RMS von Strom und Spannung
            Irms = np.sqrt(np.mean(np.square(i[0:int(N / 2)])))

            Urms = np.sqrt(np.mean(np.square(v[0:int(N / 2)])))

            # Wirk- und Blindleistung
            P = np.dot(v[0:int(N / 2)], i[0:int(N / 2)]) / (N / 2)
            Q = np.dot(v[0:int(N / 2)], i[Shift_90:int(N / 2) + Shift_90]) / (N / 2)

            # Scheinleistung CAVEAT: Gilt nur wenn Strom und Spannung Sinusförmig sind!!
            S = Urms * Irms
            D = np.sqrt(S ** 2 - P ** 2 - Q ** 2)
            PF = P / S

            if verbose:
                print("Urms = %0.4f [V]" % Urms)
                print("Irms = %0.4f [A]" % Irms)
                print("P = %0.4f [W]" % P)
                print("Q = %0.4f [VA]" % Q)
                print("S = %0.4f [Var]" % S)
                print("PF = %0.4f [-]" % PF)
                print("D = %0.4f [-]" % D)

            return {'Urms': Urms, 'Irms': Irms, 'P': P, 'Q': Q, 'S': S, 'PF': PF, 'D': D}

        @property
        def I(self):
            """
            Returns the an array of generated current samples
            :return: array
            """
            if not (self.sig_I is None):
                return self.sig_I.calc()
            else:
                return self._I

        @I.setter
        def I(self, I):
            """
            Sets local current samples
            :param I:
            :return:
            """
            self._I = I

        @property
        def V(self):
            """
            Returns an array of voltage samples
            :return:
            """
            if not (self.sig_V is None):
                return self.sig_V.calc()

    class Signal(object):

        def __init__(self, A1, A3, A5, A7, Phi1, Phi3, Phi5, Phi7, Shift, fs, f0, w0, t):
            self.A1 = A1
            self.A3 = A3
            self.A5 = A5
            self.A7 = A7
            self.Phi1 = Phi1
            self.Phi3 = Phi3
            self.Phi5 = Phi5
            self.Phi7 = Phi7
            self.Shift = Shift
            self.fs = fs
            self.f0 = f0
            self.w0 = w0
            self.t = t

            self.calc()

        def update(self, w0, t):
            """
            Function to update important simulation parameters

            :param w0: Angular frequency [rad/s]
            :param t: time vector [s]
            :return:
            """
            self.w0 = w0
            self.t = t
            self.calc()

        def calc(self):
            """
            Calculates the signal time-series according to the provided parameters
            :return: numpy array
            """
            # Angles of the fundamental component and harmonics
            phi1 = np.deg2rad(self.Phi1 + self.Shift)
            phi3 = np.deg2rad(self.Phi3 + self.Shift * 3)
            phi5 = np.deg2rad(self.Phi5 + self.Shift * 5)
            phi7 = np.deg2rad(self.Phi7 + self.Shift * 7)

            # instantaneous values
            H1 = np.sqrt(2) * self.A1 * np.sin(1 * self.w0 * self.t + phi1)
            H3 = np.sqrt(2) * self.A3 * np.sin(3 * self.w0 * self.t + phi3)
            H5 = np.sqrt(2) * self.A5 * np.sin(5 * self.w0 * self.t + phi5)
            H7 = np.sqrt(2) * self.A7 * np.sin(7 * self.w0 * self.t + phi7)

            return H1 + H3 + H5 + H7

    class Sensor(object):

        def __init__(self, name, pos_x, pos_y, angle=None):
            self.name = name
            self.x = pos_x
            self.y = pos_y
            self.angle = angle
            self.B_fields = []