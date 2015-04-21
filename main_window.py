#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import wx
import matplotlib.pyplot as plt
import numpy as np
import pyramid_sparse as pyr
import pyramid_process, distance
from os import path
import string
from gettext import install
from matplotlib.patches import Rectangle
from threading import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import specter
import pylab
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

BLOCK = False

# Button definitions
ID_START = wx.NewId()
ID_STOP = wx.NewId()

# Define notification event for thread completion
EVT_RESULT_ID = wx.NewId()

def EVT_RESULT(win, func):
    """Define Result Event."""
    win.Connect(-1, -1, EVT_RESULT_ID, func)

class ResultEvent(wx.PyEvent):
    """Simple event to carry arbitrary result data."""
    def __init__(self, data):
        """Init Result Event."""
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_RESULT_ID)
        self.data = data

# Thread class that executes processing
class WorkerThread(Thread):
    """Worker Thread Class."""
    def __init__(self, notify_window, base_folder, size_pyramid, factor, filter_value):
        """Init Worker Thread Class."""
        Thread.__init__(self)
        self._notify_window = notify_window
        self._want_abort = 0
        self.base_folder = base_folder
        self.size_pyramid = size_pyramid
        self.factor = factor
        self.filter = filter_value
        # This starts the thread running on creation, but you could
        # also make the GUI thread responsible for calling this
        wx.CallAfter(self.start)

    def run(self):
        """Run Worker Thread."""
        # This is the code executing in the new thread. Simulation of
        # a long process (well, 10s here) as a simple loop - you will
        # need to structure your processing so that you periodically
        # peek at the abort variable
        if self.filter:
            pyramid = pyr.build_and_filter(self.base_folder, self.size_pyramid, self.factor)
        else:
            pyramid = pyr.build(self.base_folder, self.size_pyramid, self.factor)
        lev = pyr.level(pyramid, 2)
        # Here's where the result would be returned (this is an
        # example fixed result of the number 10, but it could be
        # any Python object)
        wx.PostEvent(self._notify_window, ResultEvent(pyramid))

    def abort(self):
        """abort worker thread."""
        # Method for use by main thread to signal an abort
        self._want_abort = 1

class MainWindow(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.CAPTION | wx.CLOSE_BOX | wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.SYSTEM_MENU | wx.RESIZE_BORDER | wx.CLIP_CHILDREN
        wx.Frame.__init__(self, *args, **kwds)
        self.Bind(wx.EVT_CLOSE, lambda event: self.Destroy())
        
        self.chromosomes = ["None"]
        self.nb_levels = 0
        self.panel = wx.Panel(self)
        self.control_box = wx.StaticBox(self.panel)

        self.button_loader = wx.Button(self.panel, wx.ID_ANY, ("Load Data"))
        self.button_loader.Bind(wx.EVT_BUTTON, self.OnDir)
        
        self.Status_dir = [wx.StaticText(self.panel, wx.ID_ANY, ("No data")), wx.StaticText(self.panel, wx.ID_ANY, ("No data"))]
        
        self.button_viz = wx.Button(self.panel, wx.ID_ANY, ("Visualize"))
        self.button_viz.Bind(wx.EVT_BUTTON, self.OnVisualize)
        self.button_viz.Disable()
        
        self.sync = wx.CheckBox(self.panel, wx.ID_ANY, label="Keep matrices in sync")
        self.sync.SetValue(False)
        self.sync.Bind(wx.EVT_CHECKBOX, self.OnSync)
        
        self.auto_thresh = wx.CheckBox(self.panel, wx.ID_ANY, label="Auto-adjust threshold")
        self.auto_thresh.SetValue(False)
        
        self.browse_label = wx.StaticBox(self.panel, label="Matrix browsing")
        self.browse_X_label = wx.StaticBox(self.panel, label="X")
        self.browse_X_chr_label = wx.StaticBox(self.panel, label="Chromosome/Contig")
        self.browse_X_chr = wx.Choice(self.panel, wx.ID_ANY, choices=self.chromosomes)
        self.browse_X_chr.SetStringSelection("1")
        self.browse_X_start = wx.TextCtrl(self.panel, wx.ID_ANY, "0")
        self.browse_X_start_label = wx.StaticBox(self.panel, label="X start")
        self.browse_X_end = wx.TextCtrl(self.panel, wx.ID_ANY, "0")
        self.browse_X_end_label = wx.StaticBox(self.panel, label="X end")
                
        self.browse_Y_label = wx.StaticBox(self.panel, label="Y")        
        self.browse_Y_chr_label = wx.StaticBox(self.panel, label="Chromosome/Contig")
        self.browse_Y_chr = wx.Choice(self.panel, wx.ID_ANY, choices=self.chromosomes)
        self.browse_X_chr.SetStringSelection("1")
        self.browse_Y_start_label = wx.StaticBox(self.panel, label="Y start")
        self.browse_Y_start = wx.TextCtrl(self.panel, wx.ID_ANY, "0")
        self.browse_Y_end_label = wx.StaticBox(self.panel, label="Y end")
        self.browse_Y_end = wx.TextCtrl(self.panel, wx.ID_ANY, "0")
        
        self.pixel_value_label = wx.StaticBox(self.panel, label="Pixel value:")
        self.pixel_value = wx.StaticBox(self.panel, label="None")
        
        self.pyramid_box = wx.StaticBox(self.panel, label="Pyramid settings")
        
        self.pyramid_label = wx.StaticText(self.panel, wx.ID_ANY, ("Level"))     
        self.pyramid_level = wx.SpinCtrl(self.panel, wx.ID_ANY, "0",min=0,max=self.nb_levels)
        self.pyramid_level.Disable()
        
        self.matrix_mode_label = wx.StaticText(self.panel, wx.ID_ANY, ("Matrix mode"))
        self.matrix_mode = wx.Choice(self.panel, wx.ID_ANY, choices=[_("Raw"), _("Normalized"), _("Correlation"), _("Convolution"), ("Correlation/Convolution")])
        self.matrix_mode.SetStringSelection("Raw")
        self.matrix_mode.Bind(wx.EVT_CHOICE, self.OnSelectMode)
        self.matrix_mode.Disable()
        
        self.vmax_label = wx.StaticText(self.panel, wx.ID_ANY, ("Saturation threshold"))
        self.vmax = wx.TextCtrl(self.panel, wx.ID_ANY, "99%")
        self.vmax.Disable()
        
        self.marker_label = wx.StaticText(self.panel, wx.ID_ANY, ("Marker size"))
        self.marker_size = wx.TextCtrl(self.panel, wx.ID_ANY, "0.0000000000000000005")
        self.marker_size.Disable()

        self.norm_box = wx.StaticBox(self.panel, label="Normalization settings")
        
        self.norm_label = wx.StaticText(self.panel, wx.ID_ANY, ("Norm"))
        self.norm = wx.Choice(self.panel, wx.ID_ANY, choices=[_("matrix-wise"),_("fragment-wise"), _("SCN"), _("mirnylib"), _("sparsity")])
        self.norm.SetStringSelection("mirnylib")
        self.norm.Bind(wx.EVT_CHOICE, self.OnSelectMode)
        self.norm.Disable()
        
        self.norm_iteration_label = wx.StaticText(self.panel, wx.ID_ANY, ("Norm iteration number"))
        self.norm_iteration = wx.SpinCtrl(self.panel, wx.ID_ANY, "3",min=0,max=999)
        self.norm_iteration.Disable()
        
        self.norm_exposant_label = wx.StaticText(self.panel, wx.ID_ANY, ("Normalization power (or log)"))
        self.norm_exposant = wx.TextCtrl(self.panel, wx.ID_ANY, "1")
        self.norm_exposant.Disable()

        self.norm_order_label = wx.StaticText(self.panel, wx.ID_ANY, ("Normalization order"))
        self.norm_order = wx.TextCtrl(self.panel, wx.ID_ANY, "1")
        self.norm_order.Disable()
        
        self.convolution_box = wx.StaticBox(self.panel, label="Convolution settings")
        
        self.convolution_label = wx.StaticText(self.panel, wx.ID_ANY, ("Convolution kernel sigma"))
        self.convolution = wx.TextCtrl(self.panel, wx.ID_ANY, "1.5")
        self.convolution.Disable()
        
        self.convolution_gaussian_label = wx.StaticText(self.panel, wx.ID_ANY, ("Convolution iterations"))
        self.convolution_gaussian = wx.SpinCtrl(self.panel, wx.ID_ANY, "1",min=0,max=999)
        self.convolution_gaussian.Disable()
        
        self.current_plot_id = 0

        self.color_code_label = wx.StaticText(self.panel, wx.ID_ANY, ("Color code"))
        self.color_code = wx.Choice(self.panel, wx.ID_ANY, choices=pyramid_process.cmaps)
        self.color_code.SetStringSelection('jet')
        self.color_code.Disable()
        
        self.button_reset = wx.Button(self.panel, wx.ID_ANY, ("Reset"))
        self.button_reset.Bind(wx.EVT_BUTTON, self.OnReset)
        
        self.zoom = wx.Button(self.panel, wx.ID_ANY, ("Zoom"))
        self.zoom.Bind(wx.EVT_BUTTON, self.OnZoom)
        
        self.data_mode=["reads","reads"]
        
        #main arrays - one for each canvas
        self.data = [{'raw':None, 'processed':None, 'buffer':None},{'raw':None, 'processed':None, 'buffer':None}]
        self.annotation = [{},{}]
        self.zoom_coordinates = [[None,None],[None,None]]
        #levels chosen
        self.level = [None,None]
        #canvas size
        self.dpi = [110,110]
        #these are used to bind matplotlib events
        self.clickid = [0,0]
        self.selectid = [0,0]
        self.moveid = [0,0]
        #figure and canvas stuff
        self.fig = [Figure((3.0, 7.0), dpi=self.dpi[0],facecolor='1'),Figure((3.0, 7.0), dpi=self.dpi[1])]
        self.canvas = [FigCanvas(self.panel, -1, self.fig[0]),FigCanvas(self.panel, -1, self.fig[1])]
        for fig in self.fig:
            fig.set_tight_layout(True)
        
        #matplotlib event bindings - one for each canvas
        self.clickid[0] = self.canvas[0].mpl_connect('button_press_event', lambda event: self.OnClickCanvas(event, 0))
        self.selectid[0] = self.canvas[0].mpl_connect('button_release_event', lambda event: self.OnSelectCanvas(event, 0))
        self.moveid[0] = self.canvas[0].mpl_connect('motion_notify_event', lambda event: self.OnMouseCanvas(event, 0))
        self.moveid[0] = self.canvas[0].mpl_connect('axes_leave_event', lambda event: self.OnLeaveCanvas(event, 0))

        self.clickid[1] = self.canvas[1].mpl_connect('button_press_event', lambda event: self.OnClickCanvas(event, 1))
        self.selectid[1] = self.canvas[1].mpl_connect('button_release_event', lambda event: self.OnSelectCanvas(event, 1))
        self.moveid[1] = self.canvas[1].mpl_connect('motion_notify_event', lambda event: self.OnMouseCanvas(event, 1))
        self.moveid[1] = self.canvas[1].mpl_connect('axes_leave_event', lambda event: self.OnLeaveCanvas(event, 1))
        
        #adding plots to each canvas
        self.axes = [self.fig[0].add_subplot(111),self.fig[1].add_subplot(111)]
        self.plot_data = None        
        
        #whatever your mouse hovers
        self.current_x_coordinate = 0
        self.current_y_coordinate = 0
        
        #when selecting a square, the place you clicked your mouse button - one for each canvas
        self.first_x_coordinate = [0,0]
        self.first_y_coordinate = [0,0]
        
        #when selecting a square, the place you release your mouse button - one for each canvas
        self.second_x_coordinate = [0,0]
        self.second_y_coordinate = [0,0]
        
        #if you're currently selecting
        self.select_mode = False
        
        #selection rectangles - one for each canvas
        self.select_square = [Rectangle((0,0),0,0),Rectangle((0,0),0,0)]
        
        #whatever chromosome or contig your mouse is hovering
        self.current_chr = [None,None]
        
        #plots for each canvas
        self.base_axes = [self.fig[0].add_subplot(111),self.fig[1].add_subplot(111)]
        
        #this is used to keep track of your original data after zooming
        self.previous_axes = [None,None]
        
        self.fasta_file = ""
        
        self.statusbar = self.CreateStatusBar()
        
        self.__set_properties()
        self.__do_layout()
        self.create_menu()
        
        #each pyramid contains every level of binning
        self.pyramid = [None,None]
        
    def create_menu(self):
        self.menubar = wx.MenuBar()

        menu_file = wx.Menu()
        m_load = menu_file.Append(-1,"&Load data\tL", "Load pyramid data")
        self.Bind(wx.EVT_MENU, self.OnDir, m_load)
        m_manual_load = menu_file.Append(-1, "&Load manually\tM", "Load a single matrix")
        self.Bind(wx.EVT_MENU, self.OnManualLoad, m_manual_load)
        m_pdb_load = menu_file.Append(-1, "&Load PDB\tP", "Load a PDB file")
        self.Bind(wx.EVT_MENU, self.OnPDBLoad, m_pdb_load)
        
        menu_file.AppendSeparator()
        
        m_expt = menu_file.Append(-1, "&Save plot\tS", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.OnSavePlot, m_expt)
        m_export = menu_file.Append(-1, "&Export matrix\tE", "Export matrix to txt")
        self.Bind(wx.EVT_MENU, self.OnExport, m_export)

        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, lambda event: self.Destroy(), m_exit)

        self.menubar.Append(menu_file, "&File")
        
        menu_view = wx.Menu()
        m_viz = menu_view.Append(-1, "&Visualize plot\tEnter", "Visualize plot")
        self.Bind(wx.EVT_MENU, self.OnVisualize, m_viz)
        m_newplot = menu_view.Append(-1, "&New plot\tCtrl-Enter", "New plot")
        self.Bind(wx.EVT_MENU, self.OnNewPlot, m_newplot)
        m_binkb = menu_view.Append(-1, "&Kb binning\tAlt-B", "Re-bin by base pairs")
        self.Bind(wx.EVT_MENU, self.OnRebinKb, m_binkb)
        m_zoom = menu_view.Append(-1, "&Zoom\tZ", "Zoom to coordinates")
        self.Bind(wx.EVT_MENU, self.OnZoom, m_zoom)
        
        self.menubar.Append(menu_view, "&View")
        
        menu_plot = wx.Menu()
        m_distance = menu_plot.Append(-1, "&Distance law\tCtrl-D", "Plot distance law")
        self.Bind(wx.EVT_MENU, self.OnDistance, m_distance)
        m_sparsity = menu_plot.Append(-1, "&Contact sparsity\tCtrl-G", "Plot genome sparsity")
        self.Bind(wx.EVT_MENU, self.OnGenomicSparsity, m_sparsity)
        
        self.menubar.Append(menu_plot, "&Plot")
        
        menu_structure = wx.Menu()
        m_shrec = menu_structure.Append(-1, "ShRec 3D\tCtrl-H", "3D Reconstruction")
        self.Bind(wx.EVT_MENU, self.OnShrec, m_shrec)
        m_pdb = menu_structure.Append(-1, "PDB\tAlt-P", "Export to PDB")
        self.Bind(wx.EVT_MENU, self.OnPDB, m_pdb)
        
        self.menubar.Append(menu_structure, "&Structure")

        self.SetMenuBar(self.menubar)
        
    def OnSavePlot(self, event):
        file_choices = "PNG (*.png)|*.png"
        id_canvas = self.current_plot_id
        dlg = wx.FileDialog(
            self,
            message="Save plot as...",
            defaultDir="~/",
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)

        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas[id_canvas].print_figure(path, dpi=self.dpi[id_canvas]*2)
            print("Saved to %s" % path)
    
    def EraseSelection(self, id_canvas):
        self.select_square[id_canvas].set_width(0)
        self.select_square[id_canvas].set_height(0)
    
    def OnZoom(self, event=None):
        
        id_canvas = self.current_plot_id
        sync = self.sync.GetValue()
        
        if np.abs(self.select_square[id_canvas].get_width())>0.0 and np.abs(self.select_square[id_canvas].get_height())>0.0:
            self.axes[id_canvas]=self.previous_axes[id_canvas]
            self.canvas[id_canvas].draw()
            x1,y1,x2,y2 = tuple([coord for element in self.zoom_coordinates[id_canvas] for coord in element])
            x1,x2 = min(x1,x2),max(x1,x2)
            y1,y2 = min(y1,y2),max(y1,y2)
            self.axes[id_canvas].set_xlim(x1,x2)
            self.axes[id_canvas].set_ylim(y2,y1)
            print "I zoom to the square of coordinates "+str((x1,y1))+" and "+str((x2,y2))
            self.EraseSelection(id_canvas)
            self.canvas[id_canvas].draw()
            
            if sync: #yeah this is ugly there's certainly a nicer way to do it
                self.axes[not id_canvas]=self.previous_axes[not id_canvas]
                self.canvas[not id_canvas].draw()
                self.axes[not id_canvas].set_xlim(x1,x2)
                self.axes[not id_canvas].set_ylim(y1,y2)
                self.EraseSelection(not id_canvas)
                self.canvas[not id_canvas].draw()
                
        else:
            self.zoom_frame = ZoomWindow(None, wx.ID_ANY, "'")
            self.zoom_frame.main_window = self
            app.SetTopWindow(self.zoom_frame)
            self.zoom_frame.Show()
    
    def OnSync(self, event):
        if not (self.data[0] is None) and (not self.data[1] is None):
            if self.sync.GetValue():
                self.OnSelectCanvas(None, 0)
        else:
            self.data_part[0]=self.data[0]
            self.data_part[1]=self.data[1]

    def OnRebinKb(self, event):
        id_canvas = self.current_plot_id
        x1,y1,x2,y2 = tuple([coord for element in self.zoom_coordinates[id_canvas] for coord in element])
        x1,x2 = min(x1,x2),max(x1,x2)
        y1,y2 = min(y1,y2),max(y1,y2)
        M = pyramid_process.square(self.data[id_canvas]['buffer'][x1:x2,y1:y2])
        N = pyramid_process.rebin_kb(M,self.level[id_canvas].S_o_A_frags)
        self.data[id_canvas]['buffer'] = N
        self.display(N)
        
    def OnGenomicSparsity(self, event):
        id_canvas = self.current_plot_id
        s = self.data[id_canvas]['buffer']
        s_total = np.sqrt(np.sum(s))
        sum_sparsity = [np.sum(R) for R in s]/s_total
        try:
            sparsity = self.level[id_canvas].sparsity
            print sparsity
        except AttributeError:
            "warning: no pyramid data"

        try:
            plt.plot(sum_sparsity)
            plt.plot(sparsity)
        except:
            print "Data in left plot is corrupted or absent, skipping"
        plt.show()
        
    def OnShrec(self, event, pdb=False):
        
        id_canvas = self.current_plot_id
        x1,y1,x2,y2 = tuple([coord for element in self.zoom_coordinates[id_canvas] for coord in element])
        x1,x2 = min(x1,x2),max(x1,x2)
        y1,y2 = min(y1,y2),max(y1,y2)
        M = pyramid_process.square(self.data[id_canvas]['buffer'][x1:x2,y1:y2])
        N = pyramid_process.trim_for_sparsity(M)
        N = specter.remove_disconnected_matrix(N)
        
        s = specter.contact_to_coordinates(N,n=3)
        X,Y,Z = specter.coordinates_to_pdb(s)       
        
        if self.data_mode[id_canvas]=="reads":
            contigs = np.array(self.level[id_canvas].S_o_A_frags['id_c']).astype(int)
            print contigs.shape
            m = min(M.shape)
            n = len(self.chromosomes)
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            index = []
            md = np.diag(M)
            nd = np.diag(N)
            print md,nd
            for j in range(min(nd.shape)):
                k = 0
                while md[j+k] != nd[j] and k<2*m:
                    k+=1
                index.append(k+j)
            indices = np.array(index)
            colors = plt.get_cmap(self.color_code.GetStringSelection())(np.linspace(0, 1, n))
            
            for i in range(1,n+1):
                ax.scatter(X[contigs[indices]==i],Y[contigs[indices]==i],Z[contigs[indices]==i], color=colors[i-1])
            plt.show()            
                
        else:
            
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(X,Y,Z)
            plt.show(block=BLOCK)
        
        if pdb:
            n = min(len(X),len(Y),len(Z))
            

            if self.data_mode[id_canvas] == "reads":
                letters = (string.ascii_uppercase+string.ascii_lowercase+string.digits+string.punctuation)*int(n/94 +1)
                length = np.abs(self.level[id_canvas].S_o_A_frags['len_bp'][index])
                print length
                length = np.sqrt(length)
                print length
                lmax = np.max(length)
                length *= 10.0/lmax
                print length
                length[length<0.1] = 0.1
                print length
                length = np.around(length,2)
                print length
                
                
            else:
                letters = "A"*(n+1)
                length = [0.00]*(n+1)
            
            filename = self.Status_dir[id_canvas].GetLabel()+".pdb"
            with open(filename, 'w') as f:
                for i in range(1,n):
                    line = "ATOM" #1-4 "ATOM"
                    line += "  " #5-6 unused
                    line += str(i).rjust(5) #7-11 atom serial number
                    line += " " #12 unused
                    line += "OW".rjust(4) #13-16 atom name
                    line += " " #17 alternate location indicator
                    line += "SOL" #18-20 residue name
                    line += " " #21 unused
                    line += letters[contigs[i]-1] #22 chain identifier
                    line += str(i).rjust(4) #23-26 residue sequence number
                    line += " " #27 code for insertion of residues
                    line += "   " #28-30 unused
                    line += str(X[i]).rjust(8) #31-38 X orthogonal Å coordinate
                    line += str(Y[i]).rjust(8) #39-46 Y orthogonal Å coordinate
                    line += str(Z[i]).rjust(8) #47-54 Z orthogonal Å coordinate
                    line += "1.00".rjust(6) #55-60 Occupancy
                    line += str(length[i-1]).rjust(6) #61-66 Temperature factor
                    line += "      " #67-72 unused
                    line += "    " #73-76 segment identifier 
                    line += "O".rjust(2) #77-78 element symbol
                    line += "\n"
                    f.write(line)
            
    def OnPDB(self, event):
        self.OnShrec(event,pdb=True)        
    
    def OnDistance(self, event):
        id_canvas = self.current_plot_id
        plt.xscale("log")
        plt.yscale("log")
        x1,y1,x2,y2 = tuple([coord for element in self.zoom_coordinates[id_canvas] for coord in element])
        x1,x2 = min(x1,x2),max(x1,x2)
        y1,y2 = min(y1,y2),max(y1,y2)
        try:
            plt.plot(distance.distance_diagonal_law((self.data[0]['buffer'][x1:x2,y1:y2])))
        except:
            print "Data in left plot is corrupted or absent, skipping"
        try:
            plt.plot(distance.distance_diagonal_law((self.data[1]['buffer'][x1:x2,y1:y2])))
        except:
            print "Data in right plot is corrupted or absent, skipping"
        plt.show()
        
    def OnExport(self, event):
        id_canvas = self.current_plot_id
        file_choices = "txt (*.txt)|*.txt| gz (*.gz)|*.gz"
        x1,y1,x2,y2 = tuple([coord for element in self.zoom_coordinates[id_canvas] for coord in element])
        x1,x2 = min(x1,x2),max(x1,x2)
        y1,y2 = min(y1,y2),max(y1,y2)
        s = self.data[id_canvas]['buffer'][x1:x2,y1:y2]
        if np.any(s):
            dlg = wx.FileDialog(
                self,
                message="Save matrix as...",
                defaultDir="~/",
                defaultFile="matrix.txt",
                wildcard=file_choices,
                style=wx.SAVE)
    
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                np.savetxt(path, s)
                print("Saved to %s" % path)
        else:
            wx.MessageBox("Loaded data is corrupt or absent", "Pyramid data error", wx.OK|wx.ICON_ERROR)
            
    def acquire_level(self, id_canvas=None):
        
        if id_canvas is None:
            id_canvas = self.current_plot_id
        
        pyramid_level = self.pyramid_level.GetValue()
        
        if self.data_mode[id_canvas] == "reads":
            print "Loading reads."
            #acquiring data
            if self.pyramid[id_canvas] is None:
                wx.MessageBox("Loaded data is corrupt or absent", "Pyramid data error", wx.OK|wx.ICON_ERROR)
                return False
            else:
                self.level[id_canvas] = self.pyramid[id_canvas].get_level(pyramid_level)
            print "Loaded level "+str(pyramid_level)
            self.chromosomes = tuple([str(i) for i in set(self.level[id_canvas].S_o_A_frags['id_c'])])
            self.browse_X_chr.SetItems(self.chromosomes)
            self.browse_Y_chr.SetItems(self.chromosomes)
            s = self.level[id_canvas].sparse_mat_csr
            s = np.array(s.todense(),dtype=np.float64)
            s_temp = s.astype(float)
            diag = np.diag(s.diagonal())
            transpose = s.T
            s = s_temp + transpose - diag
            self.data[id_canvas]['raw'] = np.copy(s)
            print "Converted to full matrix"   
            
        elif self.data_mode[id_canvas] == "manual":
            print "j'ai chargé la pyramide manuellement"
            if pyramid_level <= len(self.pyramid[id_canvas])-1:
                s = self.pyramid[id_canvas][pyramid_level]
            else:
                s = np.copy(self.pyramid[id_canvas][len(self.pyramid[id_canvas])-1])
                for i in range(pyramid_level-len(self.pyramid[id_canvas])+1):
                   s = pyramid_process.manual_bin(s)
                   self.pyramid[id_canvas].append(s)
                
            self.data[id_canvas]['raw'] = s
                    
        elif self.data_mode[id_canvas] == "model":
            print "j'ai chargé une matrice modèle"
        
        else:
            pass
            
        self.data[id_canvas]['buffer'] = self.data[id_canvas]['processed'] = self.data[id_canvas]['raw']
        shape = self.data[id_canvas]['buffer'].shape      
        
        self.zoom_coordinates[id_canvas] = [[0,0],list(shape)]
        return True
            
    def display(self, M, id_canvas=None):

        if id_canvas is None:
            id_canvas = self.current_plot_id
            
        colormap = self.color_code.GetStringSelection()
        if self.auto_thresh.GetValue():
            m = np.mean(M)
            std = np.std(M)
            self.vmax.SetValue(str(m+std))
            
        threshold = self.vmax.GetValue()
        if threshold[-1] == "%":
            threshold = np.percentile(M,np.float64(threshold[:-1]))
        markersize = float(self.marker_size.GetValue())        
        
        if id_canvas >= 0:
            print self.zoom_coordinates
            print id_canvas
            x1,y1,x2,y2 = tuple([coord for element in self.zoom_coordinates[id_canvas] for coord in element])
            x1,x2 = min(x1,x2),max(x1,x2)
            y1,y2 = min(y1,y2),max(y1,y2)
            print (x1,y1,x2,y2)
            s = M[x1:x2,y1:y2]            
            print "j'affiche la matrice "+str(id_canvas)
            
            self.axes[id_canvas].clear()
            self.axes[id_canvas].spy(s, precision=0, markersize=markersize)
            self.base_axes[id_canvas].spy(s, precision=0, markersize=markersize, markeredgewidth=0)
            canvas_image = self.axes[id_canvas].imshow(s,vmin=0,vmax=threshold,interpolation='nearest')
            canvas_image.set_cmap(colormap)
            self.canvas[id_canvas].draw()

        else:
            plt.figure()
            plt.spy(M, precision=0, markersize=float(self.marker_size.GetValue()))
            try:
                plot_image = plt.imshow(M, vmin=0,vmax=threshold,interpolation='nearest')
                plot_image.set_cmap(colormap)
            except MemoryError:
                print "Warning: memory overload"
            except TypeError:
                print "Type Error, skipping imshow"
            
            plt.show(block=BLOCK)

    def process(self, s):
        
        if s is None:
            print "Can't process empty data!"
            return
            
        print "Matrix size: "+str(s.shape)
        matrix_name = self.matrix_mode.GetStringSelection().lower()
        normalization = self.norm.GetStringSelection()
        iterations = self.norm_iteration.GetValue()        
        norm_exposant=float(self.norm_exposant.GetValue())
        convolution_gaussian_number = self.convolution_gaussian.GetValue()
        convolution_sigma = float(self.convolution.GetValue())
        order = self.norm_order.GetValue()
        M = pyramid_process.process(s,matrix_name,normalization,order,iterations,norm_exposant,convolution_gaussian_number,convolution_sigma)
        print "Processing "+matrix_name.lower()+" matrix"
        return M
        
    def OnVisualize(self, event, id_canvas=None):
        if id_canvas is None:
            id_canvas = self.current_plot_id        
        self.visualize(id_canvas)
        if self.sync.GetValue():
            self.visualize(id_canvas=not id_canvas)
        
    def visualize(self,id_canvas=None):
        if id_canvas is None:
            id_canvas = self.current_plot_id
        
        self.acquire_level(id_canvas=id_canvas)
        M = self.data[id_canvas]['raw']
        P = self.process(M)
        
        if P is None:
            print "Can't visualize empty data!"
            return
            
        self.data[id_canvas]['processed'] = np.copy(P)
        self.data[id_canvas]['buffer'] = P
        self.display(P, id_canvas=id_canvas)
        
    def OnNewPlot(self, event):
        id_canvas = self.current_plot_id
        
        self.acquire_level(id_canvas=id_canvas)
        M = self.data[id_canvas]['raw']
        P = self.process(M)
        self.display(P, id_canvas=-1)
        
    def OnReset(self, event, id_canvas=None):
        if id_canvas is None:
            id_canvas = self.current_plot_id  
        self.reset(id_canvas)
        if self.sync.GetValue():
            self.reset(id_canvas=not id_canvas)

    def reset(self, id_canvas=None):
        if id_canvas is None:
            id_canvas = self.current_plot_id        
        P = self.data[id_canvas]['buffer']
        original_shape = self.data[id_canvas]['buffer'].shape
        self.zoom_coordinates[id_canvas] = [[0,0],list(original_shape)]
        self.display(P,id_canvas=id_canvas)
        
    def OnDir(self, event):
        # Args below are: parent, question, dialog title, default answer
        dd = wx.DirDialog(None, "Select directory to open", bank_path, 0, (10, 10), wx.Size(400, 300))

        # This function returns the button pressed to close the dialog
        ret = dd.ShowModal()

        self.data_directory = ""
        if ret == wx.ID_OK:
            dd.Destroy()
            self.data_directory = dd.GetPath()
            print('You selected: %s\n' % self.data_directory)
            self.frame_2 = LoaderWindow(None, wx.ID_ANY, "'")
            self.frame_2.data_directory = self.data_directory
            self.frame_2.main_window = self
            app.SetTopWindow(self.frame_2)
            self.frame_2.Show()
        else:
            print('You clicked cancel')
            dd.Destroy()
        
    def OnManualLoad(self, event):
        dd = wx.FileDialog(None, "Select matrix file to open", bank_path)
        # This function returns the button pressed to close the dialog
        ret = dd.ShowModal()
        self.data_path = ""
        if ret == wx.ID_OK:
            self.data_path = dd.GetPath()
            print('You selected: %s\n' % self.data_path)
        else:
            print('You clicked cancel')
        dd.Destroy()
        
        if self.data_path:
            print "Loading "+self.data_path
            id_canvas = self.current_plot_id
            matrix = np.loadtxt(self.data_path)
            self.pyramid[id_canvas] = []
            self.pyramid[id_canvas].append(matrix)
            self.data[id_canvas]['raw'] = matrix
            self.data[id_canvas]['trated'] = matrix
            self.data[id_canvas]['buffer'] = matrix
            
            n = np.min(self.data[id_canvas]['buffer'].shape)
            ideal_size = 800
            
            folder = self.data_path.split('/')[-1][:15]
            
            self.Status_dir[id_canvas].SetLabel(folder)
            self.Status_dir[id_canvas].SetForegroundColour((0,130,0))
            font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)
            self.Status_dir[id_canvas].SetFont(font)
            self.button_viz.Enable()
            self.pyramid_level.Enable()
            self.pyramid_level.SetRange(0,int(n/3))
            self.pyramid_level.SetValue(int(n/ideal_size))
            self.vmax.Enable()
            self.marker_size.Enable()
            self.matrix_mode.Enable()
            self.color_code.Enable()
            self.data_mode[id_canvas] = "manual"
            
    def OnPDBLoad(self, event):
        dd = wx.FileDialog(None, "Select PDB file to open", bank_path)
        # This function returns the button pressed to close the dialog
        ret = dd.ShowModal()
        self.data_path = ""
        if ret == wx.ID_OK:
            self.data_path = dd.GetPath()
            print('You selected: %s\n' % self.data_path)
        else:
            print('You clicked cancel')
        dd.Destroy()
        
        if self.data_path:
            print "Loading "+self.data_path
            id_canvas = self.current_plot_id
            matrix = specter.pdb_to_contact(self.data_path)
            self.pyramid[id_canvas] = matrix
            self.data[id_canvas]['raw'] = matrix
            self.data[id_canvas]['trated'] = matrix
            self.data[id_canvas]['buffer'] = matrix
            print matrix.shape
            
            n = np.min(self.data[id_canvas]['buffer'].shape)
            ideal_size = 800
            
            folder = self.data_path.split('/')[-1][:15]
            
            self.Status_dir[id_canvas].SetLabel(folder)
            self.Status_dir[id_canvas].SetForegroundColour((0,130,0))
            font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)
            self.Status_dir[id_canvas].SetFont(font)
            self.button_viz.Enable()
            self.pyramid_level.Enable()
            self.pyramid_level.SetRange(0,int(n/3))
            self.pyramid_level.SetValue(int(n/ideal_size))
            self.vmax.Enable()
            self.marker_size.Enable()
            self.matrix_mode.Enable()
            self.color_code.Enable()
            self.data_mode[id_canvas] = "manual"

    def OnSelectMode(self, event):
        matrix_name = self.matrix_mode.GetStringSelection().lower()
        norm_name = self.norm.GetStringSelection().lower()
        print "Selecting "+norm_name
      
        if "norm" in matrix_name or "corr" in matrix_name or "conv" in matrix_name:
            self.norm_exposant.Enable()
            self.norm.Enable()
            self.norm_iteration.Enable()
            
            if not norm_name in ["scn", "mirnylib"]:
                self.norm_iteration.Disable()
            else:
                self.norm_iteration.Enable()
            
            if not "wise" in norm_name:
                self.norm_order.Disable()
            else:
                self.norm_order.Enable()
                
            if "convolution" in matrix_name:
                self.convolution.Enable()
                self.convolution_gaussian.Enable()
            else:
                self.convolution.Disable()
                self.convolution_gaussian.Disable()
        else:
            self.convolution.Disable()
            self.convolution_gaussian.Disable()
            self.norm_exposant.Disable()
            self.norm.Disable()
            self.norm_iteration.Disable()
        
    def OnClickCanvas(self, event, id_canvas):
        
        if event.xdata and event.ydata:
            self.first_x_coordinate[id_canvas] = event.xdata
            self.first_y_coordinate[id_canvas] = event.ydata
        self.current_plot_id = id_canvas
        self.fig[id_canvas].set_facecolor('1')
        self.EraseSelection(id_canvas)
        self.canvas[id_canvas].draw()
        self.fig[not id_canvas].set_facecolor('0.75')
        self.canvas[not id_canvas].draw()
        self.previous_axes = self.axes
        font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        default_font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        self.Status_dir[id_canvas].SetFont(font)
        self.Status_dir[not id_canvas].SetFont(default_font)        
        if self.first_x_coordinate[id_canvas] and self.first_y_coordinate[id_canvas]:
            self.select_mode=True
    
    def OnSelectCanvas(self, event, id_canvas):
        
        if not event.xdata is None:
            self.second_x_coordinate[id_canvas] = event.xdata
        if not event.ydata is None:
            self.second_y_coordinate[id_canvas] = event.ydata
        self.select_mode = False
        
        if self.second_x_coordinate[id_canvas] != self.first_x_coordinate[id_canvas] and self.first_y_coordinate[id_canvas] != self.second_y_coordinate[id_canvas]:
            
            coordinates_status = "Display "+str(id_canvas+1)+", x1 = "+str(self.first_x_coordinate[id_canvas])+", y1 = "+str(self.first_y_coordinate[id_canvas])+\
                                 "x2 = "+str(self.second_x_coordinate[id_canvas])+", y2 = "+str(self.second_y_coordinate[id_canvas])
                                 
            self.statusbar.SetStatusText(coordinates_status)
            
            first_x_ordered, second_x_ordered = int(min(self.first_x_coordinate[id_canvas],self.second_x_coordinate[id_canvas])),\
                                                int(max(self.first_x_coordinate[id_canvas],self.second_x_coordinate[id_canvas]))
            second_y_ordered, first_y_ordered = int(min(self.first_y_coordinate[id_canvas],self.second_y_coordinate[id_canvas])),\
                                                int(max(self.first_y_coordinate[id_canvas],self.second_y_coordinate[id_canvas]))
            
            self.zoom_coordinates[id_canvas] = [[first_x_ordered,first_y_ordered],[second_x_ordered,second_y_ordered]]
            if self.sync.GetValue():
                self.zoom_coordinates[not id_canvas] = [[first_x_ordered,first_y_ordered],[second_x_ordered,second_y_ordered]]
            
    def OnMouseCanvas(self, event, id_canvas):
        #pixel coords
        if not event.xdata is None:
            self.current_x_coordinate = event.xdata
        if not event.ydata is None:    
            self.current_y_coordinate = event.ydata
        
        coordinates_status = "Display "+str(id_canvas+1)+", x = "+str(self.current_x_coordinate)+", y = "+str(self.current_y_coordinate)
        self.statusbar.SetStatusText(coordinates_status)
        
        #genome coords
        if (not self.data[id_canvas]['buffer'] is None) and self.current_x_coordinate and self.current_y_coordinate:
            self.pixel_value.SetLabel(str(float(self.data[id_canvas]['buffer'][self.current_x_coordinate,self.current_y_coordinate])))
            
            if self.data_mode[id_canvas] == "reads":
                self.browse_X_chr.SetStringSelection(str(self.level[id_canvas].S_o_A_frags['id_c'][int(self.current_x_coordinate)]))
                self.browse_Y_chr.SetStringSelection(str(self.level[id_canvas].S_o_A_frags['id_c'][int(self.current_y_coordinate)]))
                
                self.browse_X_start.SetValue(str(self.level[id_canvas].S_o_A_frags['start_bp'][int(self.current_x_coordinate)]))
                self.browse_Y_start.SetValue(str(self.level[id_canvas].S_o_A_frags['start_bp'][int(self.current_y_coordinate)]))
        
                self.browse_X_end.SetValue(str(self.level[id_canvas].S_o_A_frags['start_bp'][int(self.current_x_coordinate)]\
                                                         +self.level[id_canvas].S_o_A_frags['len_bp'][int(self.current_x_coordinate)]))
                self.browse_Y_end.SetValue(str(self.level[id_canvas].S_o_A_frags['start_bp'][int(self.current_y_coordinate)]\
                                                         +self.level[id_canvas].S_o_A_frags['len_bp'][int(self.current_y_coordinate)]))
            
        if self.select_mode and self.previous_axes and id_canvas == self.current_plot_id:
            self.select_square[id_canvas].set_width(self.current_x_coordinate - self.first_x_coordinate[id_canvas])
            self.select_square[id_canvas].set_height(self.current_y_coordinate - self.first_y_coordinate[id_canvas])
            self.select_square[id_canvas].set_xy((self.first_x_coordinate[id_canvas], self.first_y_coordinate[id_canvas]))
            self.select_square[id_canvas].set_linestyle('dashed')
            self.axes[id_canvas] = self.previous_axes[id_canvas]
            self.axes[id_canvas].add_patch(self.select_square[id_canvas])
            self.canvas[id_canvas].draw()
            
    def OnLeaveCanvas(self, event, id_canvas):
        self.select_mode = False
        
    def __set_properties(self):
        self.SetTitle(_("Pyramid Toolbox"))
        self.SetSize((350, 120))
        self.button_loader.SetMinSize((150, 27))
        self.Status_dir[0].SetMinSize((149, 17))
        self.Status_dir[1].SetMinSize((149, 17))
        self.button_viz.SetMinSize((149, 27))

    def __do_layout(self):
        
        flag_local = wx.ALIGN_CENTER_VERTICAL
        flag_sizer = wx.VERTICAL
        
        global_sizer = wx.GridBagSizer(5,5)
        
        button_sizer = wx.GridBagSizer(5,5)
                  
        control_sizer = wx.StaticBoxSizer(self.control_box, flag_sizer)
        
        local_control_sizer = wx.GridBagSizer(5,5)
        local_control_sizer.Add(self.button_loader, (0,0), flag = flag_local)
        local_control_sizer.Add(self.Status_dir[0], (2,0), flag = flag_local)
        local_control_sizer.Add(self.Status_dir[1], (2,1), flag = flag_local)
        local_control_sizer.Add(self.button_viz, (0,1), flag = flag_local)
        local_control_sizer.Add(self.button_reset, (3,1), flag = flag_local)        
        local_control_sizer.Add(self.zoom, (3,0), flag=flag_local)
        local_control_sizer.Add(self.sync, (4,0), flag = flag_local)
        local_control_sizer.Add(self.auto_thresh, (4,1), flag = flag_local)
        control_sizer.Add(local_control_sizer, flag=wx.LEFT|wx.TOP, border=5)
        
#         local_view_sizer = wx.GridBagSizer(5,5)

#         local_view_sizer.Add(self.set_selection_mode,(1,0), (1,3), flag = flag_local) 
                
        button_sizer.Add(control_sizer,(0,0))
#         button_sizer.Add(local_view_sizer,(1,0))
                
        browse_sizer = wx.StaticBoxSizer(self.browse_label, flag_sizer)
        
        local_browse_sizer = wx.GridBagSizer(5,10)
        
        browse_X_sizer = wx.StaticBoxSizer(self.browse_X_label, flag_sizer)
        local_browse_X_sizer = wx.GridBagSizer(5,10)
        local_browse_X_sizer.Add(self.browse_X_chr_label,(0,0),flag=flag_local)
        local_browse_X_sizer.Add(self.browse_X_chr,(0,1),flag=flag_local)
        local_browse_X_sizer.Add(self.browse_X_start_label,(1,0),flag=flag_local)
        local_browse_X_sizer.Add(self.browse_X_start,(1,1),flag=flag_local)
        local_browse_X_sizer.Add(self.browse_X_end_label,(2,0),flag=flag_local)       
        local_browse_X_sizer.Add(self.browse_X_end,(2,1),flag=flag_local) 
        browse_X_sizer.Add(local_browse_X_sizer,flag=wx.LEFT|wx.TOP, border=5)
        
        browse_Y_sizer = wx.StaticBoxSizer(self.browse_Y_label, flag_sizer)
        local_browse_Y_sizer = wx.GridBagSizer(5,10)
        local_browse_Y_sizer.Add(self.browse_Y_chr_label,(0,0),flag=flag_local)
        local_browse_Y_sizer.Add(self.browse_Y_chr,(0,1),flag=flag_local)
        local_browse_Y_sizer.Add(self.browse_Y_start_label,(1,0),flag=flag_local)
        local_browse_Y_sizer.Add(self.browse_Y_start,(1,1),flag=flag_local)
        local_browse_Y_sizer.Add(self.browse_Y_end_label,(2,0),flag=flag_local)
        local_browse_Y_sizer.Add(self.browse_Y_end,(2,1),flag=flag_local)
        browse_Y_sizer.Add(local_browse_Y_sizer,flag=wx.LEFT|wx.TOP, border=5)
        
        local_pixel_sizer = wx.GridBagSizer(5,10) 
        local_pixel_sizer.Add(self.pixel_value_label,(0,0),flag=flag_local)
        local_pixel_sizer.Add(self.pixel_value,(0,1),flag=flag_local)
        local_browse_sizer.Add(browse_X_sizer,(0,0))
        local_browse_sizer.Add(browse_Y_sizer,(0,1))
        local_browse_sizer.Add(local_pixel_sizer,(1,0))
        browse_sizer.Add(local_browse_sizer,flag=wx.LEFT|wx.TOP, border=5)
        
        button_sizer.Add(browse_sizer,(0,1))
        
        
        pyramid_sizer = wx.StaticBoxSizer(self.pyramid_box, flag_sizer)
        
        local_pyramid_sizer = wx.GridBagSizer(2,5)
        local_pyramid_sizer.Add(self.pyramid_label, (0,0), flag = flag_local)
        local_pyramid_sizer.Add(self.pyramid_level, (0,1), flag = flag_local)
        local_pyramid_sizer.Add(self.vmax_label, (1,0), flag = flag_local)
        local_pyramid_sizer.Add(self.vmax, (1,1), flag = flag_local)
        local_pyramid_sizer.Add(self.marker_label, (2,0), flag = flag_local)
        local_pyramid_sizer.Add(self.marker_size, (2,1), flag = flag_local)
        local_pyramid_sizer.Add(self.matrix_mode_label, (3,0), flag = flag_local)
        local_pyramid_sizer.Add(self.matrix_mode, (3,1), flag = flag_local)
        local_pyramid_sizer.Add(self.color_code_label, (4,0), flag = flag_local)
        local_pyramid_sizer.Add(self.color_code, (4,1), flag = flag_local)
        pyramid_sizer.Add(local_pyramid_sizer, flag=wx.LEFT|wx.TOP, border=5)
        
        button_sizer.Add(pyramid_sizer,(0,2))
        
        
        norm_sizer = wx.StaticBoxSizer(self.norm_box, flag_sizer)
        
        local_norm_sizer = wx.GridBagSizer(2,5)
        local_norm_sizer.Add(self.norm_exposant_label, (0,0), flag = flag_local)
        local_norm_sizer.Add(self.norm_exposant, (0,1), flag = flag_local)
        local_norm_sizer.Add(self.norm_label, (1,0), flag = flag_local)
        local_norm_sizer.Add(self.norm, (1,1), flag = flag_local)      
        local_norm_sizer.Add(self.norm_iteration_label, (2,0), flag = flag_local)
        local_norm_sizer.Add(self.norm_iteration, (2,1), flag = flag_local)
        local_norm_sizer.Add(self.norm_order_label, (3,0), flag = flag_local)
        local_norm_sizer.Add(self.norm_order, (3,1), flag = flag_local)
        norm_sizer.Add(local_norm_sizer, flag=wx.LEFT|wx.TOP, border=5)        
        
        button_sizer.Add(norm_sizer, (0,3))
        
        
        convolution_sizer = wx.StaticBoxSizer(self.convolution_box, flag_sizer)

        local_convolution_sizer = wx.GridBagSizer(2,5)
        local_convolution_sizer.Add(self.convolution_label, (0,0), flag = flag_local)
        local_convolution_sizer.Add(self.convolution, (0,1), flag = flag_local)
        local_convolution_sizer.Add(self.convolution_gaussian_label, (1,0), flag = flag_local)
        local_convolution_sizer.Add(self.convolution_gaussian, (1,1), flag = flag_local)
        convolution_sizer.Add(local_convolution_sizer, flag=wx.LEFT|wx.TOP, border=5)
        
        button_sizer.Add(convolution_sizer, (0,4))
        
        global_sizer.Add(button_sizer,(0,0),(1,2))
        
        canvas_sizer = wx.GridBagSizer(2,5)
        canvas_sizer.Add(self.canvas[0],(0,0),flag=wx.EXPAND | wx.ALL)
        canvas_sizer.Add(self.canvas[1],(0,1),flag=wx.EXPAND | wx.ALL)
        canvas_sizer.AddGrowableRow(0)
        canvas_sizer.AddGrowableCol(1)
        canvas_sizer.AddGrowableCol(0)
       
        global_sizer.Add(canvas_sizer,(1,0),(1,4),flag=wx.EXPAND | wx.ALL)
        global_sizer.AddGrowableRow(1)
        
        self.panel.SetSizer(global_sizer)
        global_sizer.Fit(self)
        self.panel.SetMinSize(global_sizer.GetMinSize())
        self.Layout()
# end of class MainWindow

class LoaderWindow(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.CAPTION | wx.CLOSE_BOX | wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.SYSTEM_MENU | wx.RESIZE_BORDER | wx.CLIP_CHILDREN
        wx.Frame.__init__(self, *args, **kwds)
        self.label_size_pyr = wx.StaticText(self, wx.ID_ANY, _("Pyramid size"), style=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE)
        self.text_ctrl_size_pyr = wx.TextCtrl(self, wx.ID_ANY, _("3"))
        self.label_fact_comp = wx.StaticText(self, wx.ID_ANY, _("Sub sampling factor"))
        self.text_ctrl_factor = wx.TextCtrl(self, wx.ID_ANY, _("3"))
        
        self.filter = wx.CheckBox(self, wx.ID_ANY, label="Apply sparsity filter")

        self.button_create_pyr = wx.Button(self, wx.ID_ANY, _("Build Pyramid"))
        self.button_create_pyr.Bind(wx.EVT_BUTTON, self.OnCreatePyr)

        self.button_return = wx.Button(self, wx.ID_ANY, _("Load Pyramid"))
        self.button_return.Bind(wx.EVT_BUTTON, self.OnReturn)

        self.label_status = wx.StaticText(self, wx.ID_ANY, _("Status:"))
        self.ctrl_status_pyr = wx.StaticText(self, wx.ID_ANY, _("data not created"))

        self.data_directory = ""
        self.data_fasta_file = ""
        self.__set_properties()
        self.__do_layout()
        EVT_RESULT(self, self.OnResult)
        self.worker = None

    def __set_properties(self):
        self.SetTitle(_("Pyramid Builder"))
        self.SetSize((400, 300))
        self.label_size_pyr.SetMinSize((144, 17))
        self.text_ctrl_size_pyr.SetMinSize((50, 27))
        self.label_fact_comp.SetMinSize((144, 17))
        self.text_ctrl_factor.SetMinSize((50, 27))
        self.button_create_pyr.SetMinSize((149, 27))
        self.button_return.SetMinSize((149, 27))
        self.label_status.SetMinSize((144, 17))
        self.ctrl_status_pyr.SetMinSize((144, 17))

    def __do_layout(self):
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_2 = wx.FlexGridSizer(6, 2, 0, 0)
        grid_sizer_2.Add(self.label_size_pyr, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 15)
        grid_sizer_2.Add(self.text_ctrl_size_pyr, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 15)
        grid_sizer_2.Add(self.label_fact_comp, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 15)
        grid_sizer_2.Add(self.text_ctrl_factor, 0, wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 15)
        grid_sizer_2.Add(self.button_create_pyr, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 15)
        grid_sizer_2.Add(self.button_return, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_2.Add(self.label_status, 0, wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 15)
        grid_sizer_2.Add(self.ctrl_status_pyr, 0, wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_2.Add(self.filter, 0, wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 15)
        sizer_2.Add(grid_sizer_2, 1, wx.EXPAND, 0)

        self.SetSizer(sizer_2)
        self.Layout()
        self.Centre()
        
    def OnCreatePyr(self, event):
        
        base_folder = self.data_directory
        size_pyramid = int(self.text_ctrl_size_pyr.GetValue())
        factor = int(self.text_ctrl_factor.GetValue())
        filter_value = self.filter.GetValue()
        if not self.worker:
            self.worker = WorkerThread(self, base_folder, size_pyramid, factor, filter_value=filter_value)
        self.worker.main_window = self.main_window
        self.main_window.nb_levels = size_pyramid-1

    def OnReturn(self, event):
        
        id_canvas = self.main_window.current_plot_id
        folder = self.data_directory.split('/')[-2][:15]
        
        self.main_window.pyramid[id_canvas] = self.pyramid
        self.main_window.base_folder = self.data_directory
        self.main_window.fasta_file = self.data_fasta_file

        self.main_window.Status_dir[id_canvas].SetLabel(folder)
        self.main_window.Status_dir[id_canvas].SetForegroundColour((255,0,0))
        font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        self.main_window.Status_dir[id_canvas].SetFont(font)
        self.main_window.button_viz.Enable()
        self.main_window.pyramid_level.Enable()
        self.main_window.pyramid_level.SetRange(0,self.main_window.nb_levels)
        self.main_window.pyramid_level.SetValue(self.main_window.nb_levels)
        self.main_window.vmax.Enable()
        self.main_window.marker_size.Enable()
        self.main_window.matrix_mode.Enable()
        self.main_window.color_code.Enable()
        self.main_window.data_mode[id_canvas] = "reads"
        
        self.main_window.acquire_level(id_canvas)
        self.Destroy()

    def OnResult(self, event):
        """Show Result status."""
        id_canvas = self.main_window.current_plot_id
        if event.data is None:
            # Thread aborted (using our convention of None return)
            self.ctrl_status_pyr.SetLabel('pyramid not created')
        else:
            # Process results here
            self.ctrl_status_pyr.SetLabel('pyramid built!')
            self.ctrl_status_pyr.SetForegroundColour((255,0,0))
            self.pyramid = event.data
        # In either event, the worker is done
        self.worker = None
        
class ZoomWindow(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.CAPTION | wx.CLOSE_BOX | wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.SYSTEM_MENU | wx.RESIZE_BORDER | wx.CLIP_CHILDREN
        wx.Frame.__init__(self, *args, **kwds)
        
        self.X_label = wx.StaticText(self, wx.ID_ANY, _("X coordinates"))
        self.X1 = wx.TextCtrl(self, wx.ID_ANY, _("0"))
        self.X2 = wx.TextCtrl(self, wx.ID_ANY, _("0"))
        
        self.Y_label = wx.StaticText(self, wx.ID_ANY, _("Y coordinates"))
        self.Y1 = wx.TextCtrl(self, wx.ID_ANY, _("0"))
        self.Y2 = wx.TextCtrl(self, wx.ID_ANY, _("0"))

        self.button_return = wx.Button(self, wx.ID_ANY, _("Zoom"))
        self.button_return.Bind(wx.EVT_BUTTON, self.OnReturn)
        
        self.button_cancel = wx.Button(self, wx.ID_ANY, _("Cancel"))
        self.button_cancel.Bind(wx.EVT_BUTTON, lambda: self.Destroy())

        self.__set_properties()
        self.__do_layout()
        
    def __set_properties(self):
        self.SetTitle(_("Manual Zoom"))
        self.SetSize((300, 100))
        
    def __do_layout(self):
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        global_sizer_2 = wx.GridBagSizer(2,5)
        global_sizer_2.Add(self.X_label, (0,0), flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        global_sizer_2.Add(self.X1, (0,1),flag=wx.ALL | wx.ALIGN_LEFT | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        global_sizer_2.Add(self.X2, (0,2), flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        global_sizer_2.Add(self.Y_label, (1,0), flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        global_sizer_2.Add(self.Y1, (1,1),flag=wx.ALL | wx.ALIGN_LEFT | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        global_sizer_2.Add(self.Y2, (1,2), flag=wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        global_sizer_2.Add(self.button_return, (2,1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        global_sizer_2.Add(self.button_cancel, (2,2), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_BOTTOM | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL)
        sizer_2.Add(global_sizer_2, 1, wx.EXPAND, 0)

        self.SetSizer(sizer_2)
        self.Layout()
        self.Centre()
        
    def OnReturn(self, event):
        id_canvas = self.main_window.current_plot_id
        self.main_window.previous_axes = self.main_window.axes
        x1,x2 = min(float(self.X1.GetValue()), float(self.X2.GetValue())), \
                max(float(self.X1.GetValue()), float(self.X2.GetValue()))
        y1,y2 = min(float(self.Y1.GetValue()), float(self.Y2.GetValue())), \
                max(float(self.Y1.GetValue()), float(self.Y2.GetValue()))
        self.main_window.zoom_coordinates[id_canvas] = [[x1,y1],[x2,y2]]
        print "The coordinates are now "+str((x1,y1))+" and "+str((x2,y2))
        self.main_window.select_square[id_canvas].set_width(x2-x1)
        self.main_window.select_square[id_canvas].set_height(y2-y1)
        self.main_window.select_square[id_canvas].set_xy((x1,y1))
        self.Destroy()
        self.main_window.OnZoom()
        
app = wx.App(False)
bank_path = ""
main_window = None
def run(bank_folder):
    global bank_path
    bank_path = path.join(bank_folder, "analysis")  
    install("app")
    main_window = MainWindow(None, wx.ID_ANY, "")
    app.SetTopWindow(main_window)
    main_window.Show()
    app.MainLoop()

if __name__ == "__main__":
    run(bank_path)
