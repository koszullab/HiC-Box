#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import wx
from subprocess import Popen, call
import os
from Bio.Restriction import RestrictionBatch, CommOnly
from Bio import SeqIO
from multiprocessing import cpu_count
import analysis_main
from gettext import install
import main_window
from shutil import rmtree, move
from matplotlib import pyplot as plt
from scipy import weave

## By default toolbox and working directories are the same 
toolbox_directory = os.path.dirname(os.path.abspath(__file__))
working_directory = toolbox_directory

class main_frame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: main_frame.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER
        wx.Frame.__init__(self, *args, **kwds)
        
        self.panel = wx.Panel(self)

        self.title = wx.StaticText(self.panel, wx.ID_ANY, _("HiC-Box"), style=wx.ALIGN_CENTRE)
        
        self.bank_folder_label = wx.StaticText(self.panel, wx.ID_ANY, _("Output folder"))
        self.bank_folder = wx.TextCtrl(self.panel, wx.ID_ANY, "", size = (300,-1))
        self.bank_output_file = wx.Button(self.panel, wx.ID_ANY, _("Browse Output Folder"))
        self.bank_output_file.Bind(wx.EVT_BUTTON, self.OnOpenOutputDialogButton, source=None)

        #Fasta file
        self.genome_label = wx.StaticText(self.panel, wx.ID_ANY, _("Genome"))
        self.genome_name = wx.TextCtrl(self.panel, wx.ID_ANY, "", size = (300,-1))
        self.genome_file = wx.Button(self.panel, wx.ID_ANY, _("Browse FASTA Genome"))
        self.genome_file.Bind(wx.EVT_BUTTON, self.OnOpenFastaFileDialogButton, source=None)
        
        #Fastq files
        self.bank_label = wx.StaticText(self.panel, wx.ID_ANY, _("Reads"))
        self.bank_name = wx.TextCtrl(self.panel, wx.ID_ANY, "",size = (300,-1))
        self.bank_file = wx.Button(self.panel, wx.ID_ANY, _("Browse FASTQ Reads"))
        self.bank_file.Bind(wx.EVT_BUTTON, self.OnOpenFastqFileDialogButton, source=None)
        
        #Toggle all commercially available restriction enzymes
        self.site_enzyme = wx.Choice(self.panel, wx.ID_ANY, choices=sorted(CommOnly.as_string()))
        self.site_enzyme.Bind(wx.EVT_CHOICE, self.OnLoadEnzyme, source=None)
        
        #Restriction sites can also be manually entered
        self.site_label = wx.StaticText(self.panel, wx.ID_ANY, _("Enzyme or site"))
        self.site_name = wx.TextCtrl(self.panel, wx.ID_ANY, "",size=(300,-1))
        
        #Bowtie indexing
        self.Index = wx.Button(self.panel, wx.ID_ANY, _("Index"))
        self.Index.Bind(wx.EVT_BUTTON, self.OnIndex, source=None)
        
        #Fasta parsing -- renaming contigs according to number and removing N-only ones
        self.Parse = wx.Button(self.panel, wx.ID_ANY, _("Parse"))
        self.Parse.Bind(wx.EVT_BUTTON, self.OnParse, source=None)
        
        #Herve toolbox stuff
        self.Align = wx.Button(self.panel, wx.ID_ANY, _("Align"))
        self.Align.Bind(wx.EVT_BUTTON, self.OnAlign, source=None)
        
        #main window
        self.ComputeResults = wx.Button(self.panel, wx.ID_ANY, _("Pyramids"))
        self.ComputeResults.Bind(wx.EVT_BUTTON, self.OnComputeResults, source=None)
        
        self.Distribution = wx.Button(self.panel, wx.ID_ANY, _("Distribution"))
        self.Distribution.Bind(wx.EVT_BUTTON, self.OnDistribution, source=None)
        
        self.Advanced = wx.Button(self.panel, wx.ID_ANY, _("Advanced"))
        self.Advanced.Bind(wx.EVT_BUTTON, self.OnLoadAdvanced, source=None)
        
        self.Temp = wx.CheckBox(self.panel, wx.ID_ANY, label="Keep temporary files")
        
        ### DEFAULT ADVANCED PARAMETERS ###
        
        self.tag_length=6
        self.looping=False
        self.speed_looping=4
        self.quality_min=30
        self.tot_len_read=700
        self.len_paired_wise_fastq=3
        self.paired_wise_fastq=True
        self.ncpu=cpu_count()
        self.bowtie=os.path.join(working_directory, "bowtie2")

        
        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle(_("Alignment Toolbox"))
        self.title.SetForegroundColour(wx.Colour(0, 0, 0))
        self.title.SetFont(wx.Font(20, wx.TELETYPE, wx.NORMAL, wx.NORMAL, 0, "Arial"))
        self.bank_folder_label.SetForegroundColour(wx.Colour(0, 0, 255))
        self.bank_folder_label.SetFont(wx.Font(13, wx.TELETYPE, wx.NORMAL, wx.NORMAL, 0, "Arial"))
        self.genome_label.SetFont(wx.Font(13, wx.TELETYPE, wx.NORMAL, wx.NORMAL, 0, "Arial"))
        self.bank_label.SetFont(wx.Font(13, wx.TELETYPE, wx.NORMAL, wx.NORMAL, 0, "Arial"))
        self.site_label.SetFont(wx.Font(13, wx.TELETYPE, wx.NORMAL, wx.NORMAL, 0, "Arial"))
        self.site_enzyme.SetStringSelection("HpaII")

    def __do_layout(self):
        global_sizer = wx.GridBagSizer(5,5)
        
        title_sizer = wx.GridBagSizer(5,5)
        title_sizer.AddGrowableCol(0)
        title_sizer.Add(self.title, (1,1), (3,3), flag = wx.ALIGN_CENTER)
        global_sizer.Add(title_sizer, (0,0), flag = wx.ALIGN_CENTER)        
        
        folder_sizer = wx.GridBagSizer(10,10)
        main_selection = wx.GridBagSizer(8,8)
        buttons_sizer = wx.GridBagSizer(5,5)
        
        folder_sizer.Add(self.bank_folder_label, (0,0), flag = wx.ALIGN_CENTER_VERTICAL)
        folder_sizer.Add(self.bank_folder, (0,1), flag = wx.ALIGN_CENTER)
        folder_sizer.Add(self.bank_output_file, (0,2))
        folder_sizer.Add(self.Temp, (0,3))
        global_sizer.Add(folder_sizer, (1,0), flag = wx.ALIGN_CENTER)

        main_selection.Add(self.genome_label, (0,0), flag = wx.ALIGN_CENTER_VERTICAL)
        main_selection.Add(self.genome_name, (0,1))
        main_selection.Add(self.genome_file, (0,2))
        main_selection.Add(self.Index, (0,4))
        main_selection.Add(self.Parse, (0,3))
        main_selection.Add(self.bank_label, (1,0), flag = wx.ALIGN_CENTER_VERTICAL)
        main_selection.Add(self.bank_name, (1,1))
        main_selection.Add(self.bank_file, (1,2))
        main_selection.Add(self.site_label, (2,0), flag = wx.ALIGN_CENTER_VERTICAL)
        main_selection.Add(self.site_name, (2,1))
        main_selection.Add(self.site_enzyme, (2,2))
        main_selection.Add(self.Distribution, (2,3))
        main_selection.Add(self.Advanced, (2,4))
        global_sizer.Add(main_selection, (2,0), flag = wx.EXPAND)

        buttons_sizer.Add(self.Align, (0,0),(2,1), flag=wx.EXPAND)
        buttons_sizer.Add(self.ComputeResults, (0,3),(2,1), flag=wx.EXPAND)       
        buttons_sizer.AddGrowableCol(0)
        buttons_sizer.AddGrowableCol(3)

        global_sizer.Add(buttons_sizer, (3,0), flag = wx.EXPAND)
        self.panel.SetSizer(global_sizer)
        global_sizer.Fit(self)
        self.Layout()
        
    def OnOpenFastaFileDialogButton(self, event):

        filename = ""  # Use  filename as a flag
        dlg = wx.FileDialog(self, message="Choose a FASTA file")
        
        if dlg.ShowModal() == wx.ID_OK:
 
            # get the new filename from the dialog
            filename = dlg.GetPath()
        dlg.Destroy()  # best to do this sooner than later
 
        if filename: # use the file name
            self.genome_name.SetValue(filename)
            print "Loaded "+self.genome_name.GetValue()
        else: 
            print "No file selected"
            
    def OnOpenFastqFileDialogButton(self, event):
    
        filename = ""  # Use  filename as a flag
        dlg = wx.FileDialog(self, message="Choose a FASTQ file", style=wx.FD_MULTIPLE)
        
        if dlg.ShowModal() == wx.ID_OK:
        
            # get the new filename from the dialog
            filename = dlg.GetPaths()
        dlg.Destroy()  # best to do this sooner than later
        
        if filename: # use the file name
            print 'Loaded '+', '.join(filename)
            self.bank_name.SetValue(','.join(filename))
        else: 
            print "No file selected"
    
    def OnOpenOutputDialogButton(self, event):
        foldername = ""  # Use  filename as a flag
        dlg = wx.DirDialog(self, message="Choose output folder")
        
        if dlg.ShowModal() == wx.ID_OK:
        
            # get the new filename from the dialog
            foldername = dlg.GetPath()
        dlg.Destroy()  # best to do this sooner than later
        
        if foldername: # use the file name
            self.bank_folder.SetValue(foldername)
            print "Loaded "+self.bank_folder.GetValue()
        else: 
            print "No file selected"
            
    def OnSaveAs(self, event):
    
        saveFileDialog = wx.FileDialog(self, "Save XYZ file", "", "",
                                       "XYZ files (*.xyz)|*.xyz", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
    
        if saveFileDialog.ShowModal() == wx.ID_CANCEL:
            return
        
        output_stream = wx.OutputStream(saveFileDialog.GetPath())
    
        if not output_stream.IsOk():
            wx.LogError("Cannot save current contents in file '%s'."%saveFileDialog.GetPath())
            return

    def OnLoadAdvanced(self, event):
        
        self.AdvSettings = advanced(None,wx.ID_ANY,"")
        app.SetTopWindow(self.AdvSettings)
        self.AdvSettings.Show()
    
    def OnLoadEnzyme(self,event):
        Enzyme = self.site_enzyme.GetStringSelection()
        Site_Restriction = RestrictionBatch([Enzyme]).get(Enzyme).site
        self.site_name.SetValue(Site_Restriction)
        print "Loaded restriction site "+str(Site_Restriction)+" for "+str(Enzyme)
    
    def OnIndex(self,event):
        print "Indexing..."
        path_to_fasta=self.genome_name.GetValue()
        genome_name = path_to_fasta.split('/')[-1].split('.')[0]
        bowtie_path=self.bowtie
        print bowtie_path+" "+path_to_fasta
        if os.path.exists(path_to_fasta) and os.path.exists(bowtie_path):
            if not(os.path.exists(os.path.join(working_directory,"index"))):
                os.mkdir(os.path.join(working_directory,"index"))
            call([os.path.join(working_directory,"index.sh"),bowtie_path,path_to_fasta,genome_name])
        else:
            missing = "bowtie folder" if not os.path.exists(bowtie_path) else "fasta file"
            print "Error: "+missing+" not found"
    
    def OnComputeResults(self,event):
        print "Pyramid Builder loaded"
        main_window.run(self.bank_folder.GetValue())
        
    def OnAlign(self,event):
        name=self.genome_name.GetValue().split('/')[-1].split('.')[0] 
        name_bank=self.bank_name.GetValue().split('/')[-1].split('.')[0] 
        bank_folder=self.bank_folder.GetValue()
        enzyme=self.site_name.GetValue()
        genome_fasta=self.genome_name.GetValue()
        genome_fastq=self.bank_name.GetValue()
        bowtie2=self.bowtie
        
        print "Settings used:"
        print("tag_length="+str(self.tag_length))
        print("looping="+str(self.looping))
        print("speed_looping="+str(self.speed_looping))
        print("quality_min="+str(self.quality_min))
        print("tot_len_read="+str(self.tot_len_read))
        print("len_paired_wise_fastq="+str(self.len_paired_wise_fastq))
        print("paired_wise_fastq="+str(self.paired_wise_fastq))
        print("ncpu="+str(self.ncpu))
        
        analysis_result = analysis_main.analyze(name=name, 
                 name_bank=name_bank, 
                 enzyme=enzyme,
                 genome_fasta=genome_fasta,
                 genome_index=os.path.join(working_directory,"index",str(name)),
                 genome_fastq=genome_fastq,
                 tag_length=self.tag_length, 
                 looping=self.looping, 
                 speed_looping=self.speed_looping, 
                 quality_min=self.quality_min, 
                 tot_len_read=self.tot_len_read,
                 len_paired_wise_fastq=self.len_paired_wise_fastq,
                 paired_wise_fastq=self.paired_wise_fastq,
                 bowtie2=bowtie2,
                 bank_folder=bank_folder,
                 ncpu=self.ncpu)
        
        if analysis_result and not self.Temp.GetValue():
            rmtree(os.path.join(bank_folder,'tmp0'))
            sams = [f for f in os.listdir(bank_folder) if ((os.path.isfile(os.path.join(bank_folder,f))) and f.split('.')[-1]=="sam")]
            for filename in sams:               
                os.remove(os.path.join(bank_folder,filename))
                print "Removing "+filename
        
    def OnParse(self, event):
        fasta_file = self.genome_name.GetValue()
        genome_parsed = []
        genome = SeqIO.parse(fasta_file, "fasta")
        fasta_parsed = fasta_file+".parsed"
        with open(fasta_parsed,"w") as fasta_handle:
            contig_id = 1
            for record in genome:
                print "Renaming contig "+str(contig_id)
                if record.seq.count('N') != len(record):
                    record.id = str(contig_id)
                    record.description = ""
                    genome_parsed.append(record)
                    contig_id += 1 
                else:
                    print "I found a N-only contig: "+str(contig_id)+", I'm ignoring it"
            SeqIO.write(genome_parsed,fasta_handle,"fasta")
        print "Fasta file fully parsed"
        self.genome_name.SetValue(fasta_parsed)
        print "I renamed the current file in fasta field"

        
    def OnDistribution(self, event):
        fasta_file = self.genome_name.GetValue()
        enzyme = RestrictionBatch([self.site_enzyme.GetStringSelection()]).get(self.site_enzyme.GetStringSelection())
        genome = ""
        if fasta_file and enzyme:
            with open(fasta_file, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta") :
                    genome += record.seq
           
            plt.hist([len(i) for i in enzyme.catalyse(genome)],alpha=.3, bins=1000)
            plt.show()
#  end of class main_frame

class advanced(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER
        wx.Frame.__init__(self, *args, **kwds)
        
        self.panel = wx.Panel(self)
        self.title = wx.StaticText(self.panel, wx.ID_ANY, _("Advanced Settings"), style=wx.ALIGN_CENTRE)
        
        self.bowtie_label = wx.StaticText(self.panel, wx.ID_ANY, _("Bowtie folder"))
        self.bowtie_folder = wx.TextCtrl(self.panel, wx.ID_ANY, "", size = (300,-1))
        self.bowtie_button = wx.Button(self.panel, wx.ID_ANY, _("Browse Bowtie Folder"))
        self.bowtie_button.Bind(wx.EVT_BUTTON, self.OnOpenBowtieButton, source=None)

        ## Advanced parameter values are copied from main frame so as to preserve previous modifications ##
        self.tag_length = wx.StaticText(self.panel, wx.ID_ANY, _("Tag length"))
        self.tag_length_field = wx.SpinCtrl(self.panel, wx.ID_ANY, "", min=0, max=100)
        self.tag_length_field.SetValue(main_frame.tag_length)
        self.looping = wx.CheckBox(self.panel, wx.ID_ANY, _("Looping"))
        self.looping.SetValue(main_frame.looping)
        self.speed_looping = wx.StaticText(self.panel, wx.ID_ANY, _("Speed looping"))
        self.speed_field = wx.SpinCtrl(self.panel, wx.ID_ANY, "", min=0, max=100)
        self.speed_field.SetValue(main_frame.speed_looping)
        self.quality_min = wx.StaticText(self.panel, wx.ID_ANY, _("Quality min"))
        self.quality_min_field = wx.SpinCtrl(self.panel, wx.ID_ANY, "", min=0, max=1000)
        self.quality_min_field.SetValue(main_frame.quality_min)
        self.len_paired_wise_fastq = wx.StaticText(self.panel, wx.ID_ANY, _("Length paired wise FASTQ"))
        self.len_paired_wise_fastq_field = wx.SpinCtrl(self.panel, wx.ID_ANY, "", min=0, max=1000)
        self.len_paired_wise_fastq_field.SetValue(main_frame.len_paired_wise_fastq)

        self.paired_wise_fastq_field = wx.CheckBox(self.panel, wx.ID_ANY, _("Paired wise FASTQ"))
        self.paired_wise_fastq_field.SetValue(main_frame.paired_wise_fastq)
        self.paired_wise_fastq_field.Bind(wx.EVT_CHECKBOX, self.OnTogglePaired, source=None)
        self.total_len_read = wx.StaticText(self.panel, wx.ID_ANY, _("Total reads length"))
        self.total_len_read_field = wx.SpinCtrl(self.panel, wx.ID_ANY, "", min=0, max=10000)
        self.total_len_read_field.SetValue(main_frame.tot_len_read)
        
        self.ncpu = wx.StaticText(self.panel, wx.ID_ANY, _("Core number"))
        self.ncpu_field = wx.SpinCtrl(self.panel, wx.ID_ANY, "", min=0, max=100)
        self.ncpu_field.SetValue(main_frame.ncpu)
        
        self.Apply = wx.Button(self.panel, wx.ID_ANY, _("Apply"))
        self.Apply.Bind(wx.EVT_BUTTON, self.OnApply, source=None)
        
        self.Cancel = wx.Button(self.panel, wx.ID_ANY, _("Cancel"))
        self.Cancel.Bind(wx.EVT_BUTTON, self.OnCancel, source=None)

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle(_("Advanced settings"))
        self.title.SetForegroundColour(wx.Colour(0, 0, 0))
        self.title.SetFont(wx.Font(20, wx.TELETYPE, wx.NORMAL, wx.NORMAL, 0, "Arial"))

    def __do_layout(self):
        global_advanced_sizer = wx.GridBagSizer(5, 15)
        
        title_sizer = wx.GridBagSizer(5,5)
        title_sizer.Add(self.title, (0,0), (3,3), flag = wx.ALIGN_CENTER)
        global_advanced_sizer.Add(title_sizer, (0,0), (1,2), flag = wx.ALIGN_CENTER)
        
        local_bowtie_sizer = wx.GridBagSizer(5,5)
        local_bowtie_sizer.Add(self.bowtie_label,(0,0), flag = wx.ALIGN_CENTER_VERTICAL)
        local_bowtie_sizer.Add(self.bowtie_folder,(0,1), flag = wx.ALIGN_CENTER_VERTICAL)
        local_bowtie_sizer.Add(self.bowtie_button,(0,2), flag = wx.ALIGN_CENTER_VERTICAL)        
        global_advanced_sizer.Add(local_bowtie_sizer, (1,0), (1,2), flag = wx.ALIGN_CENTER_VERTICAL)
        
        local_reads_sizer = wx.GridBagSizer(5,5)
        local_reads_sizer.Add(self.tag_length, (0,0), flag = wx.ALIGN_CENTER_VERTICAL)
        local_reads_sizer.Add(self.tag_length_field, (0,1))
        local_reads_sizer.Add(self.quality_min, (1,0), flag = wx.ALIGN_CENTER_VERTICAL)
        local_reads_sizer.Add(self.quality_min_field, (1,1))
        local_reads_sizer.Add(self.ncpu, (2,0), flag = wx.ALIGN_CENTER)
        local_reads_sizer.Add(self.ncpu_field, (2,1))
        local_reads_sizer.Add(self.speed_looping, (3,0))
        local_reads_sizer.Add(self.speed_field, (3,1))
        global_advanced_sizer.Add(local_reads_sizer, (2,1))
        
        local_length_sizer = wx.GridBagSizer(5,5)
        local_length_sizer.Add(self.paired_wise_fastq_field, (0,0), flag = wx.ALIGN_CENTER_VERTICAL)
        local_length_sizer.Add(self.len_paired_wise_fastq, (1,0), flag = wx.ALIGN_CENTER_VERTICAL)
        local_length_sizer.Add(self.len_paired_wise_fastq_field, (1,1))
        local_length_sizer.Add(self.total_len_read, (2,0), flag = wx.ALIGN_CENTER_VERTICAL)
        local_length_sizer.Add(self.total_len_read_field, (2,1))
        local_length_sizer.Add(self.looping, (3,0))
        global_advanced_sizer.Add(local_length_sizer, (2,0))
        
        global_advanced_sizer.Add(self.Apply, (3,0), flag = wx.ALIGN_LEFT)
        global_advanced_sizer.Add(self.Cancel, (3,1), flag = wx.ALIGN_RIGHT)
        
        self.panel.SetSizer(global_advanced_sizer)
        global_advanced_sizer.Fit(self)
        self.Layout()
    
    def OnTogglePaired(self, event):
        if self.paired_wise_fastq_field.GetValue():
            self.len_paired_wise_fastq_field.Enable()
        else:
            self.len_paired_wise_fastq_field.Disable()
    
    def OnOpenBowtieButton(self, event):
        filename = ""  # Use  filename as a flag
        dlg = wx.DirDialog(self, message="Choose output folder", style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        
        if dlg.ShowModal() == wx.ID_OK:
        
            # get the new filename from the dialog
            filename = dlg.GetPath()
        dlg.Destroy()  # best to do this sooner than later
        
        if filename: # use the file name
            self.bowtie_folder.SetValue(filename)
            print "Loaded "+self.bowtie_folder.GetValue()
        else: 
            print "No file selected"
    
    def OnApply(self,event):
        
        main_frame.tag_length=self.tag_length_field.GetValue()
        main_frame.looping=self.looping.GetValue() 
        main_frame.speed_looping=self.speed_field.GetValue()
        main_frame.quality_min=self.quality_min_field.GetValue()
        main_frame.tot_len_read=self.total_len_read_field.GetValue()
        main_frame.len_paired_wise_fastq=self.len_paired_wise_fastq_field.GetValue()
        main_frame.paired_wise_fastq=self.paired_wise_fastq_field.GetValue() 
        main_frame.ncpu=self.ncpu_field.GetValue()
        
        default_bowtie = not bool(self.bowtie_folder.GetValue())
        if not default_bowtie:
            main_frame.bowtie=self.bowtie_folder.GetValue()
        
        print "Loading settings..."
        print("tag_length="+str(self.tag_length_field.GetValue()))
        print("quality_min="+str(self.quality_min_field.GetValue()))
        print("tot_len_read="+str(self.total_len_read_field.GetValue()))
        print("len_paired_wise_fastq="+str(self.len_paired_wise_fastq_field.GetValue()))
        print("paired_wise_fastq="+str(self.paired_wise_fastq_field.GetValue()))
        if default_bowtie:
            if os.path.exists(main_frame.bowtie):
                print("Warning, no bowtie folder specified - using "+main_frame.bowtie+" by default")
            else:
                wx.MessageBox("Could not find bowtie folder", "Bowtie not found", wx.OK|wx.ICON_ERROR)
        else:
            main_frame.bowtie=self.bowtie_folder.GetValue()
            print("bowtie2="+self.bowtie_folder.GetValue())
        print("ncpu="+str(self.ncpu_field.GetValue()))
        
        self.Destroy()
    
    def OnCancel(self, event):
        self.Destroy()
        
# end of class advanced
if __name__ == "__main__":
    install("app")
    app = wx.App(False)
    main_frame = main_frame(None, wx.ID_ANY, "")
    app.SetTopWindow(main_frame)
    main_frame.Show()
    app.MainLoop()
