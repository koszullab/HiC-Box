__author__ = 'hervemn'

import scipy.sparse as scp
import numpy as np
import matplotlib.pyplot as plt
import os
import wx

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f, 1):
            pass
    return i

class data(object):
    """ Import raw contact data
    """
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.contact_file = os.path.join(input_folder, 'abs_fragments_contacts_weighted.txt')
        self.list_frag = os.path.join(input_folder,'fragments_list.txt')
        self.output_folder = output_folder
        self.sparse_matrix_file = os.path.join(self.output_folder, 'sparse_contacts.txt')
        self.nfrags = file_len(self.list_frag) - 1

    def create_sparse_dict(self):

        self.sparse_dict = dict()
        h = open(self.contact_file, "r")
        all_lines = h.readlines()
        n_lines = len(all_lines)
        for i in range(1, n_lines):
            line = all_lines[i]
            dat = line.split()
            mates = [int(dat[0]), int(dat[1])]
            mates.sort()
            f1 = mates[0]
            f2 = mates[1]
            if f1 in self.sparse_dict:
                if f2 in self.sparse_dict[f1]:
                    self.sparse_dict[f1][f2] += 1
                else:
                    self.sparse_dict[f1][f2] = 1
            else:
                self.sparse_dict[f1] = dict()
                self.sparse_dict[f1][f2] = 1

    def dok_to_csr(self,):
        keys = self.sparse_dict.keys()
        keys.sort()
        self.out_r = []
        self.out_c = []
        self.out_d = []

        for r in keys:
            data = self.sparse_dict[r]
            for c in data.keys():
                self.out_r.append(r)
                self.out_c.append(c)
                self.out_d.append(data[c])

        self.n_on_pxls = len(self.out_d)
        self.np_csr = np.zeros((3, self.n_on_pxls), dtype=np.int32)
        self.np_csr[0, :] = self.out_r
        self.np_csr[1, :] = self.out_c
        self.np_csr[2, :] = self.out_d

        self.sparse_mat = scp.csr_matrix((self.np_csr[2, :], self.np_csr[0:2, :]), shape=(self.nfrags, self.nfrags))


# if __name__ == '__main__':
#     app = wx.App()
#     frame = wx.Frame(None, -1, '')
#     frame.SetToolTip(wx.ToolTip('HiC data loader'))
#     frame.SetPosition(wx.Point(0,0))
#     frame.SetSize(wx.Size(300,250))
#     frame.SetTitle('HiC data loader')
#     frame.Show()
#
#     app.MainLoop()