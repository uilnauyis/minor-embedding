#!/usr/bin/python
# reads in guest and host adjacency lists and gets dwave_sapi embeddinig
# usage: sapi2_embed.py guest.alist host.alist 

from minorminer import find_embedding
import sys, networkx as nx

def read_graph(infile=sys.stdin):  # CS220 adjacency list format
    n=int(infile.readline().strip())
    G=nx.empty_graph(n,create_using=nx.Graph()) 
    for u in range(n):
        neighbors=infile.readline().split()
        for v in neighbors: G.add_edge(u,int(v))
    return G

assert len(sys.argv)==3

gfile=open(sys.argv[1],'r')
Guest=read_graph(gfile)
S_size = Guest.order()
S = {}
#for i in range(S_size-1):
#    for j in range(i+1,S_size):
#        if j in Guest[i]: S[(i, j)] = 1
for i in range(S_size):
    for j in Guest[i]: S[(i, j)] = 1
    S[(i,i)]=1
#    for j in Guest[i]: if i<j: S[(i, j)] = 1
#print "guest graph: ", S
gfile.close()

hfile=open(sys.argv[2],'r')
Host=read_graph(hfile)
T_size = Host.order()
T = []
for i in range(T_size):
    for j in Host[i]: T.append((i, j))
#print "host graph: ", T
hfile.close()

#Time=1000
#embeddings = find_embedding(S, set(T), timeout=Time) 
embeddings = find_embedding(S, set(T), verbose=1) 

#print "embeddings: ", embeddings
print  embeddings
