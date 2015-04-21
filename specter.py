#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
from scipy import linalg, diag
from matplotlib import pyplot as plt
from scipy import sparse
import networkx as nx
import pyramid_process
from Bio import PDB
import scipy.spatial.distance
BLOCK = False

def contact_to_distance(M):
    N = np.zeros(M.shape)
    N[M!=0] = 1/M[M!=0]
    s = sparse.csgraph.floyd_warshall(N)
    plt.figure()
    plt.spy(s, precision=0, markersize=0.00000000000005)
    plot_image = plt.imshow(s, vmin=0,vmax=np.percentile(s,99),interpolation='nearest')
    plot_image.set_cmap("jet")
    plt.show(block=BLOCK)
    print s
    return s
    
def distance_to_gram(M):
    print M.shape
    n, m = M.shape
    bary = np.sum(np.triu(M,1))/(n**2)
    d = np.sum(M**2,0)/n - bary
    G = np.array([[(d[i]+d[j]-M[i][j]**2)/2 for i in range(n)] for j in range(m)])
    plt.figure()
    plt.spy(G, precision=0, markersize=0.000000000000005)
    
    plot_image = plt.imshow(G, vmin=0,vmax=np.abs(np.percentile(G,99)),interpolation='nearest')
    plot_image.set_cmap("jet")
    plt.show(block=BLOCK)
    print G
    return G

def gram_to_coordinates(M, n=None):
    if n == None:
        n = min(M.shape)
    
    N=M
    norm = linalg.norm(M,'fro')
    N = M/norm
    print N.shape
    m = min(N.shape)
    eigenValues,eigenVectors = linalg.eigh(N,eigvals=(m-n,m-1))
    idx = eigenValues.argsort()[::-1]
    L = eigenValues[idx]
    E = eigenVectors[:,idx]
    print L
    V = E*np.sqrt(L)
    return V

def contact_to_coordinates(M, n=None):
    if n==None:
        n = min(M.shape)
    return gram_to_coordinates(distance_to_gram(contact_to_distance(M)),n=n)
    
def remove_disconnected(M):
    central_node = np.argmax([np.count_nonzero(r) for r in M])
    total_graph = nx.from_numpy_matrix(M)
    G = total_graph.subgraph(nx.node_connected_component(total_graph,central_node))
    return G
    
def remove_disconnected_matrix(M):
    G = remove_disconnected(M)
    s = nx.to_numpy_matrix(G)
    return np.array(s)
    
def coordinates_to_pdb(V,filter=None):
    X,Y,Z = tuple(np.split(V,3,1))

    f = None
    if filter == "cube":
        f = (np.mean(X)-2*np.std(X)<X)*(X<np.mean(X)+2*np.std(X))* \
            (np.mean(Y)-2*np.std(Y)<Y)*(Y<np.mean(Y)+2*np.std(Y))* \
            (np.mean(Z)-2*np.std(Z)<Z)*(Z<np.mean(Z)+2*np.std(Z))
                
    elif filter == "sphere":
        f = (X-np.mean(X))**2+(Y-np.mean(Y))**2+(Z-np.mean(Z))**2 \
                        < np.std(X)**2 + np.std(Y)**2 + np.std(Z)**2
    
    else:
        f = abs(X) >= 0.0
               
    Xfiltered = X[f]
    Yfiltered = Y[f]
    Zfiltered = Z[f]

    Xmax = np.max(np.abs(Xfiltered))
    Ymax = np.max(np.abs(Yfiltered))
    Zmax = np.max(np.abs(Zfiltered))   
#
    X = Xfiltered*100.0/Xmax
    Y = Yfiltered*100.0/Ymax
    Z = Zfiltered*100.0/Zmax
    
    X = np.around(X,3)
    Y = np.around(Y,3)
    Z = np.around(Z,3)
    return X,Y,Z

def pdb_to_contact(filename):
    return distance_to_contact(pdb_to_distance(filename))

def pdb_to_distance(filename):
    p = PDB.PDBParser()
    structure = p.get_structure('S', filename)
    for chain in structure.get_chains():
        atoms =  [np.array(atom.get_coord()) for atom in structure.get_atoms()]
    print "Atoms retrieved"
    D = scipy.spatial.distance.pdist(atoms, 'euclidean')
    return scipy.spatial.distance.squareform(D)
    
def distance_to_contact(D):
    m = np.max(1/D[D!=0])
    n = min(D.shape)
    M = np.zeros(D.shape)
    M[D!=0] = 1/D[D!=0]
    M[D==0] = m
    return M
    
