#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 20:36:06 2020

@author: Wolfgang Jan Zucha, wolfgang.zucha@gmail.com
V01, 28.12.2020
V02, 04.01.2020
© Wolfgang Jan Zucha, Wabern, Schweiz
"""

import re
import os
import PySimpleGUI as sg
import numpy as np
import pandas as pd

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
matplotlib.use("TkAgg")



###################################FUNCTIONS###############

class XRDarray():

    def __init__(self):
        self.array = np.empty((0, 0))
        self.samplename = []

    def opensample(self, sample_name):
        """
        This opens the label
        """
        if sample_name[-2:] == "xy":
            if len(self.array) == 0:
                self.array = open_xy(sample_name)
                self.array = pd.DataFrame(self.array)
            else:
                temp = open_xy(sample_name)
                temp = pd.DataFrame(temp)
                self.array = pd.concat([self.array, temp],
                                       ignore_index=True, axis=1)
        else:
            if len(self.array) == 0:
                self.array = xrdmlwolf(sample_name)
                self.array = pd.DataFrame(self.array)
            else:
                temp = xrdmlwolf(sample_name)
                temp = pd.DataFrame(temp)
                self.array = pd.concat([self.array, temp],
                                       ignore_index=True, axis=1)

        self.samplename.append(os.path.basename(sample_name))

def open_xy(data):
    """This function opens .xy file and
    converts them to a np.array."""
    twotheta, intensity = [], []
    with open(data) as f:
        for line in f:
            row = line.split()
            twotheta.append(row[0])
            intensity.append(row[1])
    xyarray = list(zip(twotheta, intensity))
    xyarray = np.asarray(xyarray)
    xyarray = xyarray.astype(np.float)
    return xyarray

def xrdmlwolf(xrdmlfile):
    """
    Converts a xrdml file in a np.array
    -----------
    Input: xrdml file as string it can be from cubix/xpert or from
            empyrian.
    -----------
    Output: np.arry[x, 2]: 2theta and intensites as float64
    """

    temp = open(xrdmlfile)
    temp = temp.read().splitlines()

    targets = [line for line in temp if "<intensities unit" in line]
    targets = str(targets)

    if len(targets) == 2:
        targets = [line for line in temp if "<counts unit=" in line]
        targets = str(targets)

    intensities = [int(s) for s in targets.split() if s.isdigit()]
    intensities = np.asarray(intensities)

    indices = [i for i, s in enumerate(temp)
               if "<positions axis=\"2Theta\" unit=\"deg\">" in s]
    startvalue = temp[indices[0]+1]
    endvalue = temp[indices[0]+2]

    startvalue = str(re.findall(r"[-+]?\d*\.\d+|\d+", startvalue))
    endvalue = str(re.findall(r"[-+]?\d*\.\d+|\d+", endvalue))

    #removes characters from string
    nochar = ["\'", "[", "]"]
    for i in nochar:
        startvalue = startvalue.replace(i, "")
        endvalue = endvalue.replace(i, "")

    startvalue = float(startvalue)
    endvalue = float(endvalue)

    stepsize = (endvalue - startvalue)/len(intensities)
    xyarray = np.zeros((len(intensities), 2))

    xyarray[:, 0] = np.arange(startvalue, endvalue, stepsize)
    xyarray[:, 1] = intensities

    return xyarray
###################################################################
def draw_figure(canvas, figure):
    """
    This function help
    --------
    Input: the GUI canvas, the fig from matplotlib
    -------
    Returns: the Matplotlib one can embedd in the GUI
    """
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side="top", fill="both", expand=0)
    return figure_canvas_agg

def listToString(s):
    """
    This function converts a string to a list
    -----
    Input: a list
    -----
    Returns a string

    """
    # initialize an empty string
    str1 = ""

    # traverse in the string
    for ele in s:
        str1 += ele

    # return string
    return str1


def collapse(layout, key):
    """
    Helper function that creates a Column that can be later made hidden,
    thus appearing "collapsed"
    :param layout: The layout for the section
    :param key: Key used to make this seciton visible / invisible
    :return: A pinned column that can be placed directly into your layout
    :rtype: sg.pin
    """
    return sg.pin(sg.Column(layout, key=key, visible=False), shrink=False)


###################HYDROTALCITE##############################

class Hydrotalcite():

    def __init__(self):

        self.Lattice = [[0, 0, 3],
                        [0, 0, 6],
                        [1,	0, -8],
                        [1,	0, -2],
                        [1,	0, -5],
                        [1,	1, 3],
                        [1,	1, 0],
                        [1,	0, 10],
                        [1,	0, -11]]

        self.Intensities = (np.array([100, 7.8, 28.7, 27.1,
                                      24.3, 8.9, 8.3, 6.6, 4.6]))/100

        self.A = 3.0460
        self.C = 22.7720

        self.array = np.zeros((len(self.Lattice), 3))
        self.array[:, 2] = self.Intensities
        self.hydrotalcite(self.A, self.C)

        self.volumn = self.cellvolumn(self.A, self.C)


    def hydrotalcite(self, a, c):
        """
        Calculates the 2Theta values of the reflexes of the hydrotalicte
        unit cell.
        ----
        Input: a and c in Angstrom of the lattice
        """

        counter = 0
        for h, k, l in self.Lattice:
            d = 4/3*((h**2 + h*k + k**2)/a**2) + (l**2)/c**2
            d = np.sqrt(1/d)

            zweitheta = 2*np.rad2deg(np.arcsin(1.5406/(2*d)))
            self.array[counter, 0] = d
            self.array[counter, 1] = zweitheta

            counter = counter + 1

    def cellvolumn(self, a, c):
        """
        Calculates the cell volumn of a hexagonal unit cell
        ----
        Input: a and c in Angstrom of the lattice
        Return: The cellvolumn in Angstrom^3

        """
        return 0.866*(a**2)*c

##############################GUI##########################################
class GUI():

    def __init__(self):
        # load classes
        self.xrdarray = XRDarray()
        self.hydro = Hydrotalcite()  
        # load class variables
        self.filenames = []
        self.array = np.empty((0, 0))
        self.hydrogui = True # to activate the hydro-widget
        self.linearaxsis = True # to activate the linear/sqrt-widget

        col3 = [[sg.T("Adjust a and c of hydrotalcite")],
                [self.a_slider()],
                [self.c_slider()],
                [sg.T("Cell volume: "+str(self.hydro.volumn)[:6]+
                      u"\u212B"+"^3",
                      key="cell")]
                ]

        col2 = [[self.listbox(self.filenames)],
                [self.remove_button("Remove selected sample")],
                [collapse(col3, "hydro")]
                ]

        layout = [
            [self.menu_bar()],
            [sg.Canvas(key="CANVAS"), sg.Col(col2)]
            ]
        # axis labels
        self.ylab = "Instensity [counts]"
        matplotlib.rc('font', size=6)
        matplotlib.rcParams['font.sans-serif'] = "Arial"
        self.fig = matplotlib.figure.Figure(figsize=(7.5, 4.5), dpi=200)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("2 Theta, Cu Anode")
        self.ax.set_ylabel(self.ylab)

        self.window = sg.Window("XRDViewerZ", layout, finalize=True,
                                location=(0, 0))

        draw_figure(self.window["CANVAS"].TKCanvas, self.fig)

    """
    Elements of the GUI
    """
    def menu_bar(self):
        menu_def = [["File", ["Open","Save Figure", "Save Array",
                              "Clear all"]],
                    ["View", ["linear/sqrt y-axis"]],
                    ["Special", ["Hydrotalcite"]],
                    ["Get Help", ["Help"]]]
        return sg.Menu(menu_def)
    
    def listbox(self, filenames):
        return sg.Listbox(values=filenames, size=(20, 20), key="files",
                          enable_events=True)
    
    def remove_button(self, text):
        return sg.B(button_text=text)


    def a_slider(self):
        return sg.Slider(range=(2, 4), key="a_lat",
                         default_value=self.hydro.A, resolution=0.005,
                         orientation="h", enable_events=True)

    def c_slider(self):
        return sg.Slider(range=(10, 35), key="c_lat",
                         default_value=self.hydro.C, resolution=0.005,
                         orientation="h", enable_events=True)

    """
    Functions of the GUI  
    """

    def open_file(self):
        filenames = sg.popup_get_file('file to open', no_window=True,
                             multiple_files=True, file_types=(
                                 ("All Files", "*.*"),
                                 ("XRDML Files","*.xrdml"),
                                 ("XY Files", "*.xy")))

        self.filenames = list(filenames)
        for i in self.filenames:
            self.xrdarray.opensample(i)

        self.array = self.xrdarray.array.to_numpy() # convert pandas to numpy
        self.ax.clear()
        self.ax.plot(self.array[:, ::2], self.array[:, 1::2], linewidth=0.5)
        self.ax.set_xlabel("2 Theta, Cu Anode")
        self.ax.set_ylabel(self.ylab)
        self.ax.set_xlim(np.min(self.array[:, ::2]),
                         np.max(self.array[:, ::2]))
        self.ax.legend(self.xrdarray.samplename, loc='upper right',
                       frameon=False)
        self.fig.canvas.draw()

        self.window["files"].update(values=self.xrdarray.samplename)

    def save_fig(self):
        path = sg.popup_get_file('message', save_as=True, no_window=True,
                                 initial_folder=os.getcwd())

        self.fig.savefig(path, dpi=300)

    def save_array(self):
        if self.array is not None:
            path = sg.popup_get_file('message', save_as=True, no_window=True,
                                     default_extension="xy",
                                     initial_folder=os.getcwd())
            try:
                np.savetxt(path, self.array, fmt='%1.4f')
            except FileNotFoundError:
                pass

    def clear_all(self):
        self.filenames = []
        self.array = []
        self.xrdarray.array = np.empty((0, 0))
        
        self.xrdarray.samplename = []
        self.window["files"].update(values=[])
        
        self.ax.clear()
        
        self.ax.set_xlabel("2 Theta, Cu Anode")
        self.ax.set_ylabel(self.ylab)
        
        self.fig.canvas.draw()
        
    def sqrtaxis(self):
        if self.linearaxsis:
            self.array[:,1::2] = np.sqrt(self.array[:, 1::2])
            self.ylab = "Sqrt " + self.ylab
            self.ax.clear()
            self.ax.plot(self.array[:, ::2], self.array[:, 1::2],
                         linewidth=0.5)
            self.ax.set_xlabel("2 Theta, Cu Anode")
            self.ax.set_ylabel(self.ylab)
            self.ax.set_xlim(np.min(self.array[:,::2]),
                             np.max(self.array[:,::2]))
            self.ax.legend(self.xrdarray.samplename, loc='upper right',
                           frameon=False)
            self.fig.canvas.draw()
            self.linearaxsis = False
        else:
            self.array[:, 1::2] = self.array[:, 1::2]**2
            self.ylab = "Instensity [counts]"
            self.ax.clear()
            self.ax.plot(self.array[:, ::2], self.array[:, 1::2],
                         linewidth=0.5)
            self.ax.set_xlabel("2 Theta, Cu Anode")
            self.ax.set_ylabel(self.ylab)
            self.ax.set_xlim(np.min(self.array[:, ::2]),
                             np.max(self.array[:, ::2]))
            self.ax.legend(self.xrdarray.samplename, loc='upper right',
               frameon=False)
            self.fig.canvas.draw()
            self.linearaxsis = True

    def helpfunc(self):
        sg.popup("For help or to report bugs \nplease "+ 
                 "contact wolfgang.zucha@gmail.com",
                 title="huhu")

    # Function of the rest

    def highlight(self, values):
        stringremove = listToString(values["files"])
        index = self.xrdarray.samplename.index(stringremove)

        self.ax.clear()
        self.ax.plot(self.array[:, ::2], self.array[:, 1::2], linewidth=0.5)
        self.ax.plot(self.array[:, 2*index], self.array[:, 2*index+1],
                     linewidth=1.5, color="red")
        self.ax.set_xlabel("2 Theta, Cu Anode")
        self.ax.set_ylabel(self.ylab)
        self.ax.set_xlim(np.min(self.array[:, ::2]), 
                         np.max(self.array[:, ::2]))
        self.ax.legend(self.xrdarray.samplename, loc='upper right',
                       frameon=False)
        self.fig.canvas.draw()

    def remove(self, values):
        try:
        # update listbox
            stringremove = listToString(values["files"])
            index = self.xrdarray.samplename.index(stringremove)
            self.xrdarray.samplename.remove(self.xrdarray.samplename[index])
            self.window["files"].update(values=self.xrdarray.samplename)

            # update plot
            self.array = np.delete(self.array, (2*index, 2*index+1), 1)
            if len(self.array[0,:]) > 1:
                self.ax.clear()
                self.ax.plot(self.array[:, ::2], self.array[:, 1::2],
                             linewidth=0.5)
                self.ax.set_xlabel("2 Theta, Cu Anode")
                self.ax.set_ylabel(self.ylab)
                self.ax.set_xlim(np.nanmin(self.array[:, ::2]),
                                 np.max(self.array[:, ::2]))
                self.ax.legend(self.xrdarray.samplename, loc='upper right',
                               frameon=False)
                self.fig.canvas.draw()
            else:
                self.ax.clear()
                self.ax.set_xlabel("2 Theta, Cu Anode")
                self.ax.set_ylabel(self.ylab)
                self.fig.canvas.draw()
        except ValueError:
            pass

    def hydroshow(self):
        if self.hydrogui:
            self.window["hydro"].update(visible=True)
            for i in range(0,len(self.hydro.array)):
                self.ax.axvline(self.hydro.array[i, 1],
                                ymax=self.hydro.array[i, 2],
                                color="red", linewidth=0.5)
                self.fig.canvas.draw()
            self.hydrogui = False
        else:
            self.window["hydro"].update(visible=False)
            self.hydrogui = True


    def hydroplot(self, values):

        self.hydro.hydrotalcite(values["a_lat"], values["c_lat"])
        self.ax.clear()
        if len(self.array) > 1:
            self.ax.plot(self.array[:,::2], self.array[:, 1::2], linewidth=0.5)
            self.ax.set_xlim(np.nanmin(self.array[:, ::2]),
                             np.max(self.array[:, ::2]))
            self.ax.legend(self.xrdarray.samplename, loc='upper right',
               frameon=False)

        for i in range(0,len(self.hydro.array)):
            self.ax.axvline(self.hydro.array[i, 1],
                            ymax=self.hydro.array[i, 2],
                       color="red", linewidth=0.5)

        self.ax.set_xlabel("2 Theta, Cu Anode")
        self.ax.set_ylabel(self.ylab)

        self.fig.canvas.draw()

        cellvol = self.hydro.cellvolumn(values["a_lat"], values["c_lat"])
        cellvol = str(cellvol)
        cellvol = "Cell volume: "+cellvol[:6]+u"\u212B"+"^3"
        self.window["cell"].update(cellvol)

def main():

    g = GUI()
    
    function = {'Open':g.open_file, "Save Figure":g.save_fig, 
                "Save Array":g.save_array, "Clear all":g.clear_all,
                "Help":g.helpfunc, "Hydrotalcite":g.hydroshow,
                "linear/sqrt y-axis":g.sqrtaxis,
                }
    
    # dict with functions that need array
    function2 = {"files":g.highlight, "Remove selected sample":g.remove, 
                 "a_lat":g.hydroplot, "c_lat":g.hydroplot}
    
    while True:
    
        event, values = g.window.read()
    
        if event in (sg.WINDOW_CLOSED, 'Exit'):
            break
        if event in function:
            function[event]()
            
        if event in function2:
            function2[event](values)
    
    g.window.close()

if __name__ == "__main__":
    main()