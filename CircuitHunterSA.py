#!/usr/bin/env python
# sim_anneal_circuit_hunter.py - Rotate paralogs by simulated annealing to
#determine the minimul number of interaction crossing edges.
# Copyright (C) 2012 Corey M. Hudson <coreymhudson(at)gmail.com>
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

__version__ = '1.0.2'

import networkx as nx
import sys
import random
import copy
import math
import time
from cStringIO import StringIO
from optparse import OptionParser
##############################################################################
#Major data structure group: 0->graph class 1->first paralog column 2->second
#paralog column 3->first singleton column 4->second singleotn column 5 -> free
#column 6 -> frozen column 7 -> singletons
##############################################################################
#May ultimately be better to code this as a struct, but for now this
#implementation gets around the process of passing the structure
##############################################################################
#Implementation will be
#class Group():
 #   def __init__(self, graph, column1, column2, singleton1, singleton2, free,
                 #frozen)
  #      self.graph = graph
   #     self.column1 = column1
    #    self.column2 = column2
     #   self.singleton1 = singleton1
     #   self.singleton2 = singleton2
     #   self.free = free
     #   self.frozen = frozen
    #def
##############################################################################


def int2bin(n, count=24):  # creates all possible 4 node data structures
    return "".join([str((n >> y) & 1) for y in range(count - 1, -1, -1)])


def get_bin_string(G, g1, g2, i, j):  # creates binary strings to search
    bin_string = StringIO()
    if G.get_edge_data(g1[i], g1[j], default=0) == 0:
        bin_string.write('0')
    else:
        bin_string.write('1')
    if G.get_edge_data(g1[i], g2[j], default=0) == 0:
        bin_string.write('0')
    else:
        bin_string.write('1')
    if G.get_edge_data(g1[j], g2[i], default=0) == 0:
        bin_string.write('0')
    else:
        bin_string.write('1')
    if G.get_edge_data(g2[i], g2[j], default=0) == 0:
        bin_string.write('0')
    else:
        bin_string.write('1')
    return bin_string.getvalue()


class Frequencies:
    Name = "Frequencies"

    def __init__(self, arg1):
        self.frequencies = {}
        for i in range(0, arg1):
            self.frequencies[int2bin(i, 4)] = 0

    def increment_frequencies(self, group):
        G = group[0]
        g1 = group[1]
        g2 = group[2]
        for i in range(len(g1) - 1):
            for j in range(i + 1, len(g1)):
                bin_string = ""
                bin_string = get_bin_string(G, g1, g2, i, j)
                self.frequencies[bin_string] += 1

    def print_graphic(self, key, num):
        if key == '0000':
            print "o o\t", num, "\no o\n"
        elif key == '0001':
            print "o o\t", num, "\n  |\no o\n"
        elif key == '0010':
            print "o o\t", num, "\n / \no o\n"
        elif key == '0011':
            print "o o\t", num, "\n /|\no o\n"
        elif key == '0100':
            print "o o\t", num, "\n \ \no o\n"
        elif key == '0101':
            print "o o\t", num, "\n \|\no o\n"
        elif key == '0110':
            print "o o\t", num, "\n X \no o\n"
        elif key == '0111':
            print "o o\t", num, "\n X|\no o\n"
        elif key == '1000':
            print "o o\t", num, "\n|\no o\n"
        elif key == '1001':
            print "o o\t", num, "\n| |\no o\n"
        elif key == '1010':
            print "o o\t", num, "\n|/ \no o\n"
        elif key == '1011':
            print "o o\t", num, "\n|/|\no o\n"
        elif key == '1100':
            print "o o\t", num, "\n|\ \no o\n"
        elif key == '1101':
            print "o o\t", num, "\n|\| \no o\n"
        elif key == '1110':
            print "o o\t", num, "\n|X \no o\n"
        elif key == '1111':
            print "o o\t", num, "\n|X| \no o\n"
        else:
            return "\n"

    def print_groups(self):
        items = self.frequencies.keys()
        items.sort()
        for key in items:
            self.print_graphic(key, self.frequencies[key])


def check_string(bin_string):  # defines the number of each node quad
    if(bin_string == '0010'):
        return 1
    elif(bin_string == '0100'):
        return 1
    elif(bin_string == '0110'):
        return 1
    elif(bin_string == '1110'):
        return 1
    elif(bin_string == '0111'):
        return 1
    else:
        return 0


def freeze_nodes(group):
    G = group[0]
    g1 = group[1]
    g2 = group[2]
    s1 = group[3]
    s2 = group[4]
    free = []
    frozen = group[6]
    singletons = group[7]
    print singletons
    for i in range(0, len(g1)):
        if(i not in frozen):
            k = 0
            j = 0
            f = 0
            while(k == 0):
                if(j == i):
                    j += 1
                elif(j >= len(g1)):
                    k = 1
                else:
                    bin_string = get_bin_string(G, g1, g2, i, j)
                    if(check_string(bin_string) == 1):
                        k = 1
                        f = 1
                    else:
                        k = 0
                    j += 1
            if(f == 0):
                frozen.append(i)
            else:
                free.append(i)
    free.sort()
    frozen.sort()
    group = []
    group.append(G)
    group.append(g1)
    group.append(g2)
    group.append(s1)
    group.append(s2)
    group.append(free)
    group.append(frozen)
    group.append(singletons)
    return group


def ret_k(k):
    #This function returns a degree value of 0 when a node has no neighbors
    #(default is to return {})
    if k == {}:
        return 0
    else:
        return k


def remove_singletons(group):
    #Function removes singletons from data structure
    #Graph class
    G = group[0]
    #Column classes
    group1 = group[1]
    group2 = group[2]
    #Singleton containers
    s1 = group[3]
    s2 = group[4]
    free = group[5]
    frozen = group[6]
    singletons = group[7]
    #Population singletons containers, based on whether the sum of the degree
    #of both genes is 1
    for i in range(len(group1)):
        if i not in frozen:
            k1 = G.degree(group1[i])
            k2 = G.degree(group2[i])
            if (ret_k(k1) + ret_k(k2)) == 1:
                singletons.append(i)
                s1.append(group1[i])
                s2.append(group2[i])
    #Remove singetons from original data structure.
    for i in range(0, len(group1)):
        if i in singletons:
            frozen.append(i)
    #Return the initial data structure minus the singleton arrays
    group = []
    group.append(G)
    group.append(group1)
    group.append(group2)
    group.append(s1)
    group.append(s2)
    group.append(free)
    group.append(frozen)
    group.append(singletons)
    return group


def round_figures(x, n):  # Function to round x to n significant figures.
    return round(x, int(n - math.ceil(math.log10(abs(x)))))


def time_string(seconds):  # Returns the time: used to calculate the cooling
    #schedule
    s = int(round(seconds))
    h, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    return '%4i:%02i:%02i' % (h, m, s)


class Annealer:
    #Performs simulated annealing by calling functions to calculate the energy,
    #make moves on the state. The temperature can be provided manually (anneal)
    #or estimated automatically (auto)

    def __init__(self, energy, move):
        self.energy = energy
        self.move = move

    def anneal(self, state, Tmax, Tmin, steps, updates=0):
        #Minimizes energy by simulated annealing
        #Keywords:
        #state = initial system arrangment
        #Tmax = maximum temperature
        #Tmin = minimum temperature (much > 0)
        #steps = number of steps requested
        #updates = number of updates to print to screen
        step = 0
        start = time.time()

        def update(T, E, acceptance, improvement):
        #Prints the current temperature, energy, acceptance rate, improvement
        #rate, elapsed time and remaining time
        #Acceptance rate is is the percent of moves since the last update that
        #were accepted by the Metropolis algorithm
        #Improvement rate indicates the percentage of moves since the last
        #update that strictly decreased the energy
            elapsed = time.time() - start
            if step == 0:
                print ' Temperature        Energy    Accept   Improve     \
                Elapsed   Remaining'
                print '%12.2f  %12.2f                      %s            ' % \
                (T, E, time_string(elapsed))
            else:
                remain = (steps - step) * (elapsed / step)
                print '%12.2f  %12.2f  %7.2f%%  %7.2f%%  %s  %s' % \
                (T, E, 100.0 * acceptance, 100.0 * improvement,
                 time_string(elapsed), time_string(remain))
        #Precompute the factor for cooling from Tmax to Tmin
        if Tmin <= 0.0:
            print 'Exponential cooling requires a minimum temperature greater\
            than zero.'
            sys.exit()
        Tfactor = -math.log(float(Tmax) / Tmin)
        #Note initial state
        T = Tmax
        E = self.energy(state)
        prevState = copy.deepcopy(state)
        prevEnergy = E
        bestState = copy.deepcopy(state)
        bestEnergy = E
        trials, accepts, improves = 0, 0, 0
        if updates > 0:
            updateWavelength = float(steps) / updates
            update(T, E, None, None)
        #Move to a new state
        while step < steps:
            step += 1
            T = Tmax * math.exp(Tfactor * step / steps)
            self.move(state)
            E = self.energy(state)
            dE = E - prevEnergy
            trials += 1
            if dE > 0.0 and math.exp(-dE / T) < random.random():
                #Restore previous state
                state = copy.deepcopy(prevState)
                E = prevEnergy
            else:
                #Accept new state and compare to best state
                accepts += 1
                if dE < 0.0:
                    improves += 1
                prevState = copy.deepcopy(state)
                prevEnergy = E
                if E < bestEnergy:
                    bestState = copy.deepcopy(state)
                    bestEnergy = E
            if updates > 1:
                if step // updateWavelength > (step - 1) // updateWavelength:
                    update(T, E, float(accepts) / trials, float(improves)\
                           / trials)
                    trials, accepts, improves = 0, 0, 0
        return bestState, bestEnergy

    def auto(self, state, minutes, steps=2000):
        #Minimize the energy state through automatic selection
        #Keyword arguments:
        #state = initial state of system
        #minutes = time spent annealing
        #steps = number of steps to spend on each stage of exploration

        def run(state, T, steps):
            E = self.energy(state)
            prevState = copy.deepcopy(state)
            prevEnergy = E
            accepts, improves = 0, 0
            for steps in range(steps):
                self.move(state)
                E = self.energy(state)
                dE = E - prevEnergy
                if dE > 0.0 and math.exp(-dE / T) < random.random():
                    state = copy.deepcopy(prevState)
                    E = prevEnergy
                else:
                    accepts += 1
                    if dE < 0.0:
                        improves += 1
                    prevState = copy.deepcopy(state)
                    prevEnergy = E
            return state, E, float(accepts) / steps, float(improves) / steps
        step = 0
        start = time.time()
        print 'Attempting automatic simulated annealing ...'
        #Find the initial guess at temperature
        T = 0.0
        E = self.energy(state)
        while T == 0.0:
            step += 1
            self.move(state)
            T = abs(self.energy(state) - E)
        print 'Exploring temperature landscape:'
        print 'Temperagure  Energy  Accept  Improve Elapsed'

        def update(T, E, acceptance, improvement):
            elapsed = time.time() - start
            print '%12.2f  %12.2f  %7.2f%%  %7.2f%%  %s' %\
            (T, E, 100.0 * acceptance, 100.0 * improvement,\
             time_string(elapsed))
        #Search for Tmax a temperature that gives 98% acceptance
        state, E, acceptance, improvement = run(state, T, steps)
        step += steps
        while acceptance > 0.98:
            T = round_figures(T / 1.5, 2)
            state, E, acceptance, improvement = run(state, T, steps)
            step += steps
            update(T, E, acceptance, improvement)
        while acceptance < 0.98:
            T = round_figures(T * 1.5, 2)
            state, E, acceptance, improvement = run(state, T, steps)
            step += steps
            update(T, E, acceptance, improvement)
        Tmax = T
        #Search for Tmin - a temperature that gives 0% improvement
        while improvement > 0.0:
            T = round_figures(T / 1.5, 2)
            state, E, acceptance, improvement = run(state, T, steps)
            step += steps
            update(T, E, acceptance, improvement)
        Tmin = T
        #Calculate anneal duration
        elapsed = time.time() - start
        duration = round_figures(int(60.0 * minutes * step / elapsed), 2)
        # Perform anneal
        print 'Annealing from %.2f to %.2f over %i steps:' % \
        (Tmax, Tmin, duration)
        return self.anneal(state, Tmax, Tmin, duration, 20)


def freeze_largest_degree(group, H):
    #Starting at the largest node optimizes faster
    G = group[0]
    #Column classes
    group1 = group[1]
    group2 = group[2]
    #Singleton containers
    s1 = group[3]
    s2 = group[4]
    free = group[5]
    frozen = group[6]
    singletons = group[7]
    max = 0
    for i in free:
        tmp = H.degree(i)
        if tmp >= max:
            max = tmp
            top = i
    frozen.append(top)
    group = []
    group.append(G)
    group.append(group1)
    group.append(group2)
    group.append(s1)
    group.append(s2)
    group.append(free)
    group.append(frozen)
    group.append(singletons)
    return group


def print_sub_frequencies(H, G, group1, group2, free, singletons):
    sub_frequencies = {}
    for i in range(0, 16):
            sub_frequencies[int2bin(i, 4)] = 0
    for m in free:
        n = H[m]
        for key in n.iterkeys():
            bin_string = get_bin_string(G, group1, group2, m, key)
            sub_frequencies[bin_string] += 1
    for m2 in singletons:
        n = H[m2]
        for key in n.iterkeys():
            bin_string = get_bin_string(G, group1, group2, m2, key)
            sub_frequencies[bin_string] += 1
    for key in sub_frequencies.iterkeys():
        print key, sub_frequencies[key]


def create_ancestral(group, H):  # generate ancestral network state
    G = group[0]
    group1 = group[1]
    group2 = group[2]
    paralog_dict = {}
    l = nx.connected_component_subgraphs(H)
    s1 = []
    s2 = []
    singletons = []
    print "Subgraphs ", len(l)
    for i in range(len(l)):
        print(l[i].nodes())
        free = l[i].nodes()
        frozen = []
        if len(free) > 1:
            for j in range(0, len(group1)):
                if j not in free:
                    frozen.append(j)
            group = []
            group.append(G)
            group.append(group1)
            group.append(group2)
            group.append(s1)
            group.append(s2)
            group.append(free)
            group.append(frozen)
            group.append(singletons)
            group = remove_singletons(group)
            group = freeze_largest_degree(group, H)
            #group, c = annealer.auto(group, 0.001)
            #group, c = annealer.anneal(group, 10000, 0.001, 1000, 50)
            group, c = annealer.anneal(group, 100000, 0.001, 500 * math.log10
                                       (len(free)), 10)
            singletons = group[7]
            for m in singletons:
                n = H[m]
                bin_string = get_bin_string(G, group1, group2, m, n.keys()[0])
                if (bin_string == '0100') or (bin_string == '0010') or \
                (bin_string == '0111') or (bin_string == '1110'):
                    group1[m], group2[m] = group2[m], group1[m]
            print_sub_frequencies(H, G, group1, group2, free, singletons)
            G = group[0]
            group1 = group[1]
            group2 = group[2]
            s1 = []
            s2 = []
            free = []
            frozen = []
            singletons = []
    free = range(0, len(group1))
    group[5] = free
    group, c = annealer.anneal(group, 10000, 0.001, 2000, 50)
    sf = Frequencies(16)
    sf.increment_frequencies(group)
    sf.print_groups()
    #print_sub_frequencies(H, G, group1, group2, free, singletons)
    return group


def self_edges(group):  # Function to remove edges connecting a gene to itself
    G = group[0]
    G1 = nx.Graph()
    G2 = nx.Graph()
    H1 = list(G.edges_iter(group[1]))
    for k in H1:
        G1.add_edges(k[0], k[1])
    H2 = list(G.edges_iter(group[2]))
    for k in H2:
        G2.add_edges(k[0], k[1])
    R = G1.copy()
    R.remove_nodes_from(n for n in G1 if n in G2)
    return R.number_of_edges()


def crossing_edges(group):  # Function counting edges across columns
    ce = Frequencies(16)
    ce.increment_frequencies(group)
    cross = 0
    for key in ce.frequencies:
        if (key == "0100") or (key == "0010") or (key == "1100") or\
        (key == "0011") or (key == "1010") or (key == "1011") or\
        (key == "1101"):
                cross += ce.frequencies[key]
        if (key == "0110") or (key == "1110") or (key == "0111") or\
        (key == "1111"):
                cross += 2 * ce.frequencies[key]
    return cross


def reassign_groups(group):  # Function that rotates gene pairs randomly
    random.seed()
    G = group[0]
    group1 = group[1]
    group2 = group[2]
    s1 = group[3]
    s2 = group[4]
    free = group[5]
    frozen = group[6]
    singletons = group[7]
    a = random.choice(free)
    #a = random.randint(0, len(group1)-1)
    group1[a], group2[a] = group2[a], group1[a]
    group = []
    group.append(G)
    group.append(group1)
    group.append(group2)
    group.append(s1)
    group.append(s2)
    group.append(free)
    group.append(frozen)
    group.append(singletons)


def remove_pair_edges(group):  # Function to remove edges connecting paralogs
    G = group[0]
    group1 = group[1]
    group2 = group[2]
    for i in range(len(group1)):
        if G.get_edge_data(group1[i], group2[i], default=0) != 0:
            G.remove_edge(group1[i], group2[i])
    group = []
    group.append(G)
    group.append(group1)
    group.append(group2)
    return group


if __name__ == '__main__':
    #Optimized compiler
    try:
        import psyco
        psyco.full()
    except ImportError:
        pass
    G = nx.Graph()
    usage = "usage: %prog [options] PARALOG EDEGFILE"
    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--paralog', dest="paralog", \
                      help="a file containing paired paralogs", \
                      metavar="PARALOG")
    parser.add_option('-e', '--edge', dest="edges", \
                      help="a file containing paired interactions", \
                      metavar="EDGEFILE")
    (options, args) = parser.parse_args()
    paralog_file = options.paralog
    edge_file = options.edges
    f1 = open(paralog_file, 'r')
    print "Reading paralog file"
    group1 = []
    group2 = []
    partition = []
    integer_dict = {}
    H = nx.Graph()
    node_int = 0
    for line in f1:
        n1, n2 = line.split()
        G.add_node(n1)
        G.add_node(n2)
        group1.append(n1)
        group2.append(n2)
        H.add_node(node_int)
        integer_dict[n1] = node_int
        integer_dict[n2] = node_int
        node_int += 1
    f1.close()
    f2 = open(edge_file, 'r')
    print "Reading edgefile"
    for line in f2:
        node1, node2 = line.split()
        G.add_edge(node1, node2)
        nn1 = integer_dict[node1]
        nn2 = integer_dict[node2]
        if nn1 != nn2:
            H.add_edge(nn1, nn2)
    f2.close()
    l = nx.connected_component_subgraphs(H)
    print "Subgraphs ", len(l)
    group = []
    group.append(G)
    group.append(group1)
    group.append(group2)
    f = Frequencies(16)
    f.increment_frequencies(group)
    f.print_groups()
    #Remove edges across pairs
    group = remove_pair_edges(group)
    #Print an initial count of the number of crossing edges
    print(crossing_edges(group))
    #Run the simulated annealer
    annealer = Annealer(crossing_edges, reassign_groups)
    #Generate the ancestral state
    group = create_ancestral(group, H)
    print(crossing_edges(group))
