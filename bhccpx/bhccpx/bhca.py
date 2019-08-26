#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# This file is part of the BHC Complexity Toolkit.
#
# The BHC Complexity Toolkit is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The BHC Complexity Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the BHC Complexity Toolkit.  If not, 
# see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
# Copyright 2019, Mark D. Flood
#
# Author: Mark D. Flood
# Last revision: 22-Jun-2019
# -----------------------------------------------------------------------------

import re

import networkx as nx

Q_FULL = 1
Q_HETERO = 2
Q_FULL_COND = 3
Q_HETERO_COND = 4


class BHCanalyzeError(RuntimeError):
    def __init__(self, msg):
        self.msg = msg
        
        
def get_labels(BHC, dimen, missings=dict()):
    """Gets all possible node labels associated with the attibute key (dimen)
    
    Scans the given BHC graph for all node attributes with attribute name
    given by dimen. The return value is a set of all unique attribute values
    found along this dimension, sorted in alphabetical order. 
    
    If a missings dictionary is provided, the function will update 
    missings with new entries indicating which nodes within BHC are 
    missing attribute values along the given dimension. 
    
    Parameters
    ----------
    BHC : networkx.DiGraph
        A directed graph representing a bank holding company
    dimen : str
        Key indicating which node attribute to scan for values
    missing : dict
        Collection of missing attributes identified for each node
    
    Returns
    -------
    labels : set
        A set containing all of the unique labels in BHC for the given dimen
        
    Examples
    --------
    Setting up a BHC graph with known node attributes:
        
    >>> import bhc_testutil as TEST
    >>> BHCa = TEST.BHC_attribDAG()
    
    Check the geographic jurisdiction labels:
        
    >>> get_labels(BHCa, 'GEO_JURISD')
    ['GERMANY', 'UNITED STATES - CA', 'UNITED STATES - NC', 'UNITED STATES - NY']

    Check the entity type labels:
        
    >>> get_labels(BHCa, 'entity_type')
    ['BHC', 'DEO', 'SMB']

    """
    labels_all = nx.get_node_attributes(BHC, dimen)
    labels = set(labels_all.values())
    for v in BHC.nodes(data=True):
        if not(dimen in v[1]):
            if not(v in missings.keys()):
                missings[v] = set()
            missing_vals = missings[v]
            missing_vals.add(dimen)
            missings[v] = missing_vals
    return sorted(labels)



def get_quotient(BHC, dimen, Qtype):
    """Calculate the quotient of BHC, partitioning by a given dimension
    
    In a quotient, the nodes of the original graph (BHC) are partitioned
    into sets that share an identical value for the indicated attribute
    (dimen). For example, nodes might be clustered into groups that each
    share the same entity type. The resulting quotient graph (BHCq) has 
    one node for each such partition or cluster. 
    
    The original graph (BHC) is converted to an undirected graph (edge
    orientations are removed) before quotienting. The resulting quotient
    graph is undirected. 
    
    There are four possible quotient types, depending on how the edges 
    from the original graph (BHC) translate to edges in the quotient 
    graph (BHCq):
        
      1. Full quotient (Q_FULL): Every edge in BHC translates to an edge
         between  the matching partition nodes in BHCq. For example, for the 
         edge (u,v) in BHC, where u has entity type 'NAT' and v has 
         entity type 'IHC', the quotient graph would have a matching
         edge (NAT, IHC).
      2. Heterogeneous quotient (Q_HETERO): A subgraph of the full quotient, 
         in which any parallel edges between two nodes are collapsed into 
         a single edge. 
      3. Full condensed quotient (Q_FULL_COND): A subgraph of the full 
         quotient, in which any self-loop edges are removed. 
      4. Heterogeneous condensed quotient (Q_HETERO_COND): A subgraph of 
         the heterogeneous quotient, in which any self-loop edges are removed. 
    
    Parameters
    ----------
    BHC : networkx.DiGraph
        A directed graph representing a bank holding company
    dimen : str
        Key indicating which node attribute to scan for values
    Qtype : int
        Indicator of the quotient type to calculate, one of:
            
            1. Q_FULL
            2. Q_HETERO
            3. Q_FULL_COND
            4. Q_HETERO_COND
    
    Returns
    -------
    BHCq : undirected networkx.MultiGraph or networkx.Graph (if condensed)
        The calculated quotient graph (undirected)

    Examples
    --------
    Setting up a BHC graph with known node attributes:
        
    >>> import bhc_testutil as TEST
    >>> BHCa = TEST.BHC_attribDAG()
    
    Full quotient (Q_FULL) by the geographic jurisdiction labels:
        
    >>> get_quotient(BHCa, 'GEO_JURISD', Q_FULL).number_of_nodes()
    4
    >>> get_quotient(BHCa, 'GEO_JURISD', Q_FULL).number_of_edges()
    6

    Full quotient (Q_FULL) by the entity type labels:
        
    >>> get_quotient(BHCa, 'entity_type', Q_FULL).number_of_nodes()
    3
    >>> get_quotient(BHCa, 'entity_type', Q_FULL).number_of_edges()
    6

    Heterogenous condensed quotient (Q_HETERO_COND) by the entity type labels:
        
    >>> get_quotient(BHCa, 'entity_type', Q_HETERO_COND).number_of_nodes()
    3
    >>> get_quotient(BHCa, 'entity_type', Q_HETERO_COND).number_of_edges()
    2

    """
    BHCu = BHC.to_undirected()
    # BHCq is the (undirected) quotient graph to be derived from BHC
    BHCq = None
    if (Q_FULL==Qtype):          # Full (multi-edges and self-loops)
        BHCq = nx.MultiGraph()
    elif (Q_HETERO==Qtype):      # Heterogeneous (multi-edges, no self-loops)
        BHCq = nx.MultiGraph()
    elif (Q_FULL_COND==Qtype):   # Condensed (self-loops, but no multi-edges)
        BHCq = nx.Graph()
    elif (Q_HETERO_COND==Qtype): # Cond. heterog. (no multi-edges/self-loops)
        BHCq = nx.Graph()
    else:
        raise ValueError('ERROR: Unrecognized Qtype: '+Qtype)
    # Establish the nodes of the quotient
    entity_types = nx.get_node_attributes(BHCu,'entity_type')
    nx.set_node_attributes(BHCq, entity_types)
    # Copy edges from the BHC to the quotient
    for e in BHCu.edges(data=True):
        parent = e[0]
        child = e[1]
        e_attrs = e[2]
        if ((dimen in BHCu.node[parent]) and (dimen in BHCu.node[child])):
            parent_label = BHCu.node[parent][dimen]
            child_label = BHCu.node[child][dimen]
            BHCq.add_edge(parent_label, child_label, attr_dict=e_attrs)
    # Remove self-loops, as appropriate
    if (2==Qtype or 4==Qtype):
        removals = list(nx.selfloop_edges(BHCq))
        BHCq.remove_edges_from(removals)
    # Wrap up
    #print('Summing up:', BHC.number_of_nodes(), BHC.number_of_edges(), BHCq.number_of_nodes(), BHCq.number_of_edges())
    return BHCq



def node_equals(u, v, G, dimen):
    """Test whether two nodes (u,v) are equal along a given attribute dimension 
    
    Two nodes (u,v) in a graph (G) are equal along an attribute 
    dimension (dimen), if u and v each have a dimen attribute defined in G,
    and if the respective values of the attributes u.dimen and v.dimen 
    are equal according to standard python semantics. 
    
    Parameters
    ----------
    u : node identifier
        A valid NetworkX node identifier
    v : node identifier
        A valid NetworkX node identifier
    G : nx.Graph or nx.DiGraph
        Any NetworkX graph object
    dimen : str
        Name of an attribute attaching to the nodes of G
    
    Returns
    -------
    testval : bool
        Whether u and v are equal in G for the attribute dimension (dimen)

    Examples
    --------
    Setting up a BHC graph with known node attributes:
        
    >>> import bhc_testutil as TEST
    >>> BHCa = TEST.BHC_attribDAG()
    
    Examine the first three nodes, for reference:
        
    >>> BHCa.nodes(data=True)[0]['GEO_JURISD']
    'UNITED STATES - NY'
    >>> BHCa.nodes(data=True)[0]['entity_type']
    'BHC'
    >>> BHCa.nodes(data=True)[1]['GEO_JURISD']
    'UNITED STATES - NC'
    >>> BHCa.nodes(data=True)[1]['entity_type']
    'SMB'
    >>> BHCa.nodes(data=True)[2]['GEO_JURISD']
    'UNITED STATES - NC'
    >>> BHCa.nodes(data=True)[2]['entity_type']
    'SMB'

    Compare two equal nodes:
        
    >>> node_equals(1, 2, BHCa, 'GEO_JURISD')
    True
    >>> node_equals(1, 2, BHCa, 'entity_type')
    True

    Compare two unequal nodes:
        
    >>> node_equals(0, 1, BHCa, 'GEO_JURISD')
    False
    >>> node_equals(0, 1, BHCa, 'entity_type')
    False
    
    """
    testval = False
    if (dimen in G.node[u] and dimen in G.node[v]):
        testval = (G.node[u][dimen] == G.node[v][dimen])
    return testval
    


def contract(BHC, dimen):
    """Contracts a BHC graph by collapsing homogeneous edges
    
    For an edge e in BHC, form a new graph by shrinking e to a point --
    that is, remove the edge e and replace its two endpoints by one new 
    vertex whose neighbors are the neighbors of both endpoints. If the 
    edge e = {v, x} contains a leaf node of BHC, where vertex v has
    degree one (v is the leaf), then contraction of edge e is equivalent 
    to removing edge e and vertex v, leaving vertex x as it was. Let us 
    call the process of contracting a graph along one edge an elementary 
    edge contraction. Contraction of a graph involves
    """
    # BHCdup is a deep copy of the BHC graph
    BHCdup = BHC.copy()
    edges_contracted = set()
    edges_uncontract = set()
    remap_edges = dict()
    for e in BHCdup.edges():
        if e in remap_edges:
            parent = remap_edges[e][0]
            child = remap_edges[e][1]
        else:
            parent = e[0]
            child = e[1]
        if parent==child:  # Assumes parent/child values are ints
            if BHC.has_edge(parent,child):
                BHC.remove_edge(parent,child)
            edges_contracted.add(e)
        elif (BHC.has_node(parent) and BHC.has_node(child) and 
              node_equals(parent, child, BHC, dimen)):
            for coe in BHC.out_edges(child, data=True):
                BHC.add_edge(parent,coe[1],**coe[2])
                remap_edges[(coe[0],coe[1])] = (parent,coe[1])
            for cie in BHC.in_edges(child, data=True):
                BHC.add_edge(cie[1],parent,**cie[2])
                remap_edges[(cie[0],cie[1])] = (cie[0],parent)
            BHC.remove_node(child)
            edges_contracted.add(e)
        else:
            edges_uncontract.add(e)
    return BHC, len(edges_contracted), len(edges_uncontract)

def get_contraction(BHC, dimen):
    #root = BHC.graph['rootnode']
    number_of_edges_contracted  = -1
    BHCcont = BHC.copy()
    while (number_of_edges_contracted != 0):
        BHCcont, number_of_edges_contracted, n_unc = contract(BHCcont, dimen)
    BHCcont.remove_edges_from(BHCcont.selfloop_edges())
    return BHCcont

def get_disjoint_maximal_homogeneous_subgraphs(BHC, dimen):
    BHCu = BHC.to_undirected()
    BHCdisj = BHCu.copy()
    for e in BHCu.edges(data=True):
        parent = e[0]
        child = e[1]
        if ((dimen in BHCu.node[parent]) and (dimen in BHCu.node[child])):
            if not(node_equals(parent, child, BHCu, dimen)):
                BHCdisj.remove_edge(parent, child)
    return BHCdisj



def number_of_components(BHC):
    """Counts the connected components in the undirected projection of BHC.
    
    Parameters
    ----------
    BHC : networkx.DiGraph
        A directed graph representing a bank holding company
    
    Returns
    -------
    int
        Components in the projection of BHC to a simple undirected graph
        
    Examples
    --------
    Using this function

    Start with a simple DAG, all in one component
    
    >>> import bhc_testutil as TEST
    >>> bhc7 =TEST.BHC_simpleDAG()
    >>> number_of_components(bhc7)
    1
    
    Cut an edge, to separate the simple DAG into two pieces
    
    >>> bhc7x =TEST.BHC_simpleDAG()
    >>> bhc7x.remove_edge(0,1)
    >>> number_of_components(bhc7x)
    2

    """
    BHCu = BHC.to_undirected()
    return nx.number_connected_components(BHCu)



def edge_count(BHC):
    """Counts the number of edges in the undirected projection of BHC.
    
    For the given BHC (DiGraph) object, count the edges in its undirected 
    equivalent. In projecting BHC to an undirected graph, any parallel edges
    (i.e., multiedges) are collapsed to a single undirected edge between the
    nodes. Self-loops are not removed in the projection to the undirected
    graph. 
    
    Parameters
    ----------
    BHC : networkx.DiGraph
        A directed graph representing a bank holding company
    
    Returns
    -------
    int
        Count of edges in the projection of BHC to a simple undirected graph
        
    Examples
    --------
    Using this function

    The count for a simple directed tree, with 7 nodes and 6 edges
    
    >>> import bhc_testutil as TEST
    >>> edge_count(TEST.BHC_simpleDAG())
    6

    The count for a directed tree, with 7 nodes and 7 edges (1 reverse edge)
    
    >>> edge_count(TEST.BHC_simpleDAG_plusreverseedge())
    6

    """
    BHCu = BHC.to_undirected()
    return BHCu.number_of_edges()


def cycle_rank(BHC):
    BHCu = BHC.to_undirected()
    b0 = number_of_components(BHCu)
    rv = BHCu.number_of_edges() - BHCu.number_of_nodes() + b0
    return rv

def aggregate_weight(BHC, weights):
    # Return the total weight of a BHC, where weights is a dict
    # providing the weight of individual subsidiaries, keyed by RSSD.
    # If a given subsidiary has no weight in the weights dict, then 
    # that subsidiary is ignored in the aggregate.
    rv = 0
    for n in BHC.nodes(data=True):
        try:
            wgt = weights[n]
            rv = rv + wgt
        except KeyError as err:
            pass
    return rv



def high_holder_impute(BHC):
    """Impute a high-holder node within a BHC, based on network characteristics
    
    Scans the directed graph of a holding-company (BHC) structure to find 
    nodes of degree zero. Three cases are possible:
               
        * Case 0: There is a single entity with indegree==0 in the BHC 
          directed graph. This node is imputed as the high holder. This is
          the typical case. 
        * Case 1: There are multiple entities with indegree==0 in the BHC
          directed graph. Among these, the node with the largest number of
          subsidiaries (other nodes reachable along directed paths) is 
          imputed as the high holder.
        * Case 2: There are no entities with indegree==0 in the BHC
          directed graph. This case is anomolous, as it implies that all
          nodes participate in a directed cycle. The node with the 
          smallest indegree is imputed as the "high holder."
    
    Parameters
    ----------
    BHC : networkx.DiGraph
        A directed graph representing a bank holding company. This should
        typically be a single (undirected) component, but this restriction 
        is not a requirement. Nodes in the BHC should have integer identifiers.
        
    Returns
    -------
    HH : int
        The node imputed as high holder
    case : int
        One of the three possible structural cases:
            * Case 0: Unique node of indegree zero
            * Case 1: Multiple nodes of indegree zero
            * Case 2: No nodes of indegree zero
            
    Examples
    --------
    Consider a graph of two holding companies (1,2,3,4) and (10,20,30), which
    together own a joint venture subgraph (400,500). The directed edges in 
    the graph consist of:
        
        * First BHC: (1,2), (2,3), (3,4)
        * Second BHC: (10,20), (20,30)
        * Joint venture: (400,500), (2,400), (20,400)
    
    Set up a connected DiGraph for this structure:
    
    >>> import networkx as nx
    >>> BHCg = nx.DiGraph()
    >>> BHCg.add_edges_from([(1,2), (2,3), (3,4)])
    >>> BHCg.add_edges_from([(10,20), (20,30)])
    >>> BHCg.add_edges_from([(400,500), (2,400), (20,400)])
    >>> sorted(BHCg.nodes)
    [1, 2, 3, 4, 10, 20, 30, 400, 500]
    >>> sorted(BHCg.edges)
    [(1, 2), (2, 3), (2, 400), (3, 4), (10, 20), (20, 30), (20, 400), (400, 500)]
    
    There are two nodes of indegree==0: 

    >>> BHCg.in_degree(1)
    0
    >>> BHCg.in_degree(10)
    0
    
    More nodes are reachable from node 1 than from node 10:
    
    >>> sorted(reachable_nodes(BHCg, 1))
    [1, 2, 3, 4, 400, 500]
    >>> sorted(reachable_nodes(BHCg, 10))
    [10, 20, 30, 400, 500]
    
    Node 1 is imputed as the high holder, under case==1:

    >>> high_holder_impute(BHCg)
    (1, 1)
    
    """
    HH = None
    # The set of all nodes with indegree==0
    IN0 = set()
    # The indegree counts for every node in BHC
    indegrees = dict()
    for ent in BHC:
        indegrees[ent] = BHC.in_degree(ent)
        if (0==BHC.in_degree(ent)):
            IN0.add(ent)
    case = -1
    if (1==len(IN0)):
        # There is a unique node with zero indegree
        HH = IN0.pop()
        case = 0
    elif (len(IN0) > 0):
        # No unique ancestor, choose node with indegree==0 and most descendants
        HHreach = dict()
        for ent in IN0:
            HHreach[ent] = len(reachable_nodes(BHC, ent))
        IN0_sort = sorted(HHreach.items(), reverse=True, key=lambda x:x[1])
        HH = IN0_sort[0][0]
        case = 1
    else:
        # No nodes with zero indegree!! Pick the node with smallest indegree
        all_sort = sorted(indegrees.items(), key=lambda x:x[1])
        HH = all_sort[0][0]
        case = 2
    return (HH, case)



def reachable_nodes(BHC, rssd):
    """Finds all nodes in BHC that are reachable from rssd (a node in BHC)
    
    Assembles the set of nodes in the BHC directed graphs that are
    reachable (via a directed path) from the rssd node.
    
    Parameters
    ----------
    BHC : networkx.DiGraph
        A directed graph representing a bank holding company. This should
        typically be a single (undirected) component, but this restriction 
        is not a requirement. 
    rssd : int
        The identifier of a node within the BHC. 
        
    Returns
    -------
    reachable : set
        The set of nodes in BHC reachable from rssd by a directed path.
            
    Examples
    --------
    Consider a graph of two holding companies (1,2,3) and (10,20,30), which
    together own a joint venture subgraph (400,500). The directed edges in 
    the graph consist of:
        
        * First BHC: (1,2), (2,3)
        * Second BHC: (10,20), (20,30)
        * Joint venture: (400,500), (2,400), (20,400)
    
    Set up a connected DiGraph for this structure:
    
    >>> import networkx as nx
    >>> BHCg = nx.DiGraph()
    >>> BHCg.add_edges_from([(1,2), (2,3)])
    >>> BHCg.add_edges_from([(10,20), (20,30)])
    >>> BHCg.add_edges_from([(400,500), (2,400), (20,400)])
    >>> sorted(BHCg.nodes)
    [1, 2, 3, 10, 20, 30, 400, 500]
    >>> sorted(BHCg.edges)
    [(1, 2), (2, 3), (2, 400), (10, 20), (20, 30), (20, 400), (400, 500)]

    The set of reachable nodes depends on the starting node:
        
    >>> sorted(reachable_nodes(BHCg, 1))
    [1, 2, 3, 400, 500]
    >>> len(reachable_nodes(BHCg, 1))
    5
    >>> sorted(reachable_nodes(BHCg, 10))
    [10, 20, 30, 400, 500]
    >>> sorted(reachable_nodes(BHCg, 20))
    [20, 30, 400, 500]
    >>> sorted(reachable_nodes(BHCg, 30))
    [30]
    
    """
    reachable = set()
    for node in BHC:
        if nx.algorithms.shortest_paths.generic.has_path(BHC, rssd, node):
            reachable.add(node)
    return reachable
            
        
    
def high_holder_partition(BHC_graph):
    """Separates a directed graph into unambigous BHCs and joint ventures
    
    A fundamental assumption is that a holding company is a directed graph of
    legal entities connected by ownership/control relationships, where there 
    is a unique consolidated reporting entity (the so-called high holder) 
    designated as the locus of overall control and/or reporting authority. 
    
    In many cases, the officially designated high holder is not available, and 
    there can be numerous candidates. This is further complicated by the 
    possibility of joint ventures, which create directed edges emanating from
    two (or more) distinct BHCs to a single subisidiary they jointly control.
    
    This function partitions the nodes of a connected digraph component into:
        
        * One or more BHC graphs, each with a unique high holder that 
          unambiguously owns/controls its subsidiary nodes
        * One of more joint venture (JV) objects, consisting of:
            * The subgraph of BHC_graph comprising the nodes that are jointly
              owned/controlled by two or more high holders
            * The set of "bridge" edges that connect the JV to the BHC
              components of the partitioned graph
    
    The resulting collection of BHCs may differ from the regulatory 
    assignment. For example, consider a case of two holding companies
    (identified as 1 and 2) which each own a 50 percent stake in a bank
    identified as 3. Regulators might declare that the entire group
    is a single BHC with node 2 as the official high holder. In contrast, 
    this function would partition the graph into two distinct holding 
    companies, each with a unique high holder (1 and 2), and a separate 
    joint venture consisting of the subgraph 3. 
    
    This function takes no consideration of NIC entity types in its analysis. 
    
    Parameters
    ----------
    BHC_graph : networkx.DiGraph
        A single connected component that may contain multiple high holders
        and one or more joint ventures
        
    Returns
    -------
    BHCs : dict
        A collection of unambiguous BHC DiGraph objects (values), each with 
        a unique high holder node (keys)
    JVs : dict
        A collection of joint venture (JV) objects, each a tuple comprising:
            * nx.DiGraph of the entities in the JV 
            * set of bridge edges that connect the JV to other BHCs
        These are collected in a dict where the keys are ordered 
        tuples of BHC high holders controlling the JVs, and each dict
        values is a list of the JV objects controlled by the high-holder team.
            
    Examples
    --------
    Consider a graph of two holding companies (1,2,3) and (10,20,30), which
    together own a joint venture subgraph (400,500). The directed edges in 
    the graph consist of:
        * First BHC: (1,2), (2,3)
        * Second BHC: (10,20), (20,30)
        * Joint venture: (400,500), (2,400), (20,400)
    
    Set up a connected DiGraph involving a joint venture:
    
    >>> import networkx as nx
    >>> BHCg = nx.DiGraph()
    >>> BHCg.add_edges_from([(1,2), (2,3)])
    >>> BHCg.add_edges_from([(10,20), (20,30)])
    >>> BHCg.add_edges_from([(400,500), (2,400), (20,400)])
    >>> sorted(BHCg.nodes)
    [1, 2, 3, 10, 20, 30, 400, 500]
    
    Partition the graph into BHCs and JVs:
    
    >>> (BHCs, JVs) = high_holder_partition(BHCg)
    
    Inspect the BHCs:
    
    >>> sorted(BHCs[1].nodes)
    [1, 2, 3]
    >>> sorted(BHCs[10].nodes)
    [10, 20, 30]
    
    There is a single joint venture here. JVs is a dict keyed by the 
    high-holder team(s) controlling each JV:
    
    >>> len(JVs)
    1
    >>> print(type(JVs))
    <class 'dict'>
    
    For each high-holder team, the dict contains a list of that team's 
    joint ventures. In this example, (1,10) is the only high-holder team:
    
    >>> HHteam = (1,10)
    >>> print(type(JVs[HHteam]))
    <class 'list'>
    
    Each element in the list is a tuple of the JV graph and its set of
    bridge edges:
    
    >>> print(type(JVs[HHteam][0]))
    <class 'tuple'>
    >>> print(type(JVs[HHteam][0][0]))
    <class 'networkx.classes.digraph.DiGraph'>
    >>> print(type(JVs[HHteam][0][1]))
    <class 'set'>
    
    For JVs[(1,10)][0], the first and only JV contolled by high holders
    1 and 10, inspect the JV nodes:
        
    >>> sorted(JVs[HHteam][0][0].nodes)
    [400, 500]
    
    Inspect the bridge edges:
        
    >>> sorted(JVs[HHteam][0][1])
    [(2, 400), (20, 400)]

    """
    # The input needs to be a single connected component
    if not(1==len([nx.connected_components(BHC_graph.to_undirected())])):
        raise ValueError('Cannot process multiple components in BHC_graph: '+
                         str(BHC_graph.nodes))
    # Only nodes with indegree==0 (the IN0 set) can possibly be high holders
    IN0 = set()
    for ent in BHC_graph.nodes:
        if (0==BHC_graph.in_degree(ent)):
            IN0.add(ent)
    # For each entity in BHC_graph, find the IN0 entities that control it
    controllers = dict()
    for in0 in IN0:
        in0_descend = reachable_nodes(BHC_graph, in0)
        for ent in in0_descend:
            if not(ent in controllers):
                controllers[ent] = set()
            controllers[ent].add(in0)
    # Find the rogues -- entities that do not descend from *any* high holder
    rogues = set(BHC_graph).difference(set(controllers.keys()))
    if (len(rogues) > 0):
        # This really should never occur; something strange is going on
        raise BHCanalyzeError('Rogue nodes detected: '+str(rogues))
    # For each node in IN0, find the set of uniquely controlled descendants
    BHC0_nodes = {in0: set([in0]) for in0 in IN0}
    JV_nodes = dict()
    for ent in controllers.keys():
        # cont is the set of IN0 elements that own/control this ent
        cont = controllers[ent]
        if (1==len(cont)):
            in0 = min(cont)
            BHC0_nodes[in0].add(ent)
        elif (len(cont) > 1):
            # The ordered tuple of IN0 nodes controlling ent is a unique key
            owners = tuple(sorted(cont))
            if not(owners in JV_nodes):
                JV_nodes[owners] = set()
            JV_nodes[owners].add(ent)
        else:
            # This really should never occur; something strange is going on
            raise BHCanalyzeError('No owner found for entity: '+str(ent))
    # Build a dictionary of the standalone BHC components (excluding the JVs)
    BHCs = dict()
    for in0 in BHC0_nodes:
        BHCs[in0] = BHC_graph.subgraph(BHC0_nodes[in0])
    # Build a dictionary of the joint ventures (JVs)
    JVs = dict()
    for jvk in JV_nodes.keys():
        JVsubd = BHC_graph.subgraph(JV_nodes[jvk])
        JVsubu = JVsubd.to_undirected()
        jvs = list()
        for jvc in nx.connected_components(JVsubu):
            # Include this JV directed subgraph
            jv = BHC_graph.subgraph(jvc)
            # Separately, include the edges that bridge from jv to the BHCs
            innies = set(BHC_graph.in_edges(jvc))
            outies = set(BHC_graph.out_edges(jvc))
            jv_bridges = (innies.union(outies)).difference(set(jv.edges))
            jvs.append((jv, jv_bridges))
        JVs[jvk] = jvs
    return (BHCs, JVs)



def find_highholders(BankSys, rssd):
    """Scans a banking system to find all high holders for a given entity.
    
    The banking system (BankSys) object is a directed graph, with edge
    orientations defined by ownership/control relationships of the NIC. 
    Starting from a given entity (rssd), this function examines the 
    banking system to find all entities that:
        
      * are ancestors of rssd (including possibly rssd itself), and
      * have no immediate parents themselves. 
     
    In other words, the high-holder list contains any nodes in the 
    hierarchy that are "upstream" of rssd, and have no parents. It is
    possible for rssd to be its own high-holder, including the case where
    it is a free-standing institution that participates in no 
    ownership/control relationships.
    
    Parameters
    ----------
    BankSys : networkx.DiGraph
        A directed graph representing the banking system at a point in time
    rssd : int
        Identifier indicating a specific legal entity node in BankSys
    
    Returns
    -------
    HHs : list
        A list of the high holders of rssd in BankSys
    BHC : networkx DiGraph
        The BHC graph containing rssd
        
    Examples
    --------
    Consider a graph of two holding companies (1,2,3) and (10,20,30), which
    together own a joint venture subgraph (400,500). The directed edges in 
    the graph consist of:
        * First BHC: (1,2), (2,3)
        * Second BHC: (10,20), (20,30)
        * Joint venture: (400,500), (2,400), (20,400)
    
    Set up a connected DiGraph involving a joint venture:
    
    >>> import networkx as nx
    >>> BHCg = nx.DiGraph()
    >>> BHCg.add_edges_from([(1,2), (2,3)])
    >>> BHCg.add_edges_from([(10,20), (20,30)])
    >>> BHCg.add_edges_from([(400,500), (2,400), (20,400)])
    >>> sorted(BHCg.nodes)
    [1, 2, 3, 10, 20, 30, 400, 500]
    
    Inspect the BHC for high holders:
    
    >>> (HHs, BHC) = find_highholders(BHCg, 3)
    >>> HHs
    [1]
    >>> (HHs, BHC) = find_highholders(BHCg, 30)
    >>> HHs
    [10]
    >>> (HHs, BHC) = find_highholders(BHCg, 400)
    >>> HHs
    [1, 10]
    >>> sorted(BHC.nodes)
    [1, 2, 3, 10, 20, 30, 400, 500]
    
    """
    BankSysU = BankSys.to_undirected()
    BHCu = nx.algorithms.components.node_connected_component(BankSysU, rssd)
    BHC = BankSys.subgraph(BHCu)
    IN0 = set()
    for ent in BHC:
        if not(0==BHC.in_degree(ent)):
            continue
        reach = nx.algorithms.shortest_paths.generic.has_path(BHC, ent, rssd)
        if (reach):
            IN0.add(ent)
    HHs = sorted(IN0)
    return (HHs, BHC)


           
def check_lei(lei, display_warnings=True):
    """Verify the checksum on a candidate Legal Entity Identifier (LEI).
    
    Under international standard ISO/CD 17442, LEIs may be issued for 
    "legal entities," which can be essentially any unique party (except 
    a natural person) that has legal standing to enter into formal contracts 
    in its jurisdiction. 
    
    The Global Legal Entity Identifier Foundation (GLEIF) maintains certain 
    data integrity rules for LEIs. One important rule is that each LEI 
    shall be a 20-character string that includes two final numeric characters
    that are the checksum under standard ISO/IEC 7064 (mod 97-10). 
    
    Parameters
    ----------
    lei : str
        A candidate LEI 
    display_warnings : bool
        Whether to display warning messages as issues are detected
    
    Returns
    -------
    check : int
        The final checksum modulus (which should be 1), indicating:
         * If the candidate LEI is checksum-valid, a code 1 is returned
         * If syntax flaws prevent calculation of a checksum, -1 is returned
         * If the candidate LEI is syntax-valid but fails the checksum, 
           a positive value other than 1 is returned

    syntax_errcodes : list
        A list of numeric error codes indicating syntax flaws in the 
        candidate LEI. If this list is empty, no syntax errors were
        detected. Possible syntax error codes (in brackets) are:
            
         1. LEI value is too long
         2. LEI value is too short
         3. LEI value is not uppercase
         4. LEI value does not match the official format (regex):
             * ``[0-9A-Z]{18}[0-9]{2}``
             
         5. LEI checkdigits are in the invalid range [00, 01, 99]
        
    Examples
    --------
    Applying this function to variations on a particular LEI. For example,
    the actual (valid) LEI for Citigroup is: ``6SHGI4ZSSLCXXQSBB395``
    
    Check the valid value for the Citigroup LEI, and see a clean result
    
    >>> check_lei('6SHGI4ZSSLCXXQSBB395', display_warnings=False)
    (1, [])

    A simple true/false answer: is it valid? A valid LEI has final modulus = 1
    
    >>> (1==check_lei('6SHGI4ZSSLCXXQSBB395', display_warnings=False)[0])
    True

    Introduce some syntax errors and try again
    
    >>> check_lei('6shgi4zsslcxxqsbb395', display_warnings=False)
    (-1, [3, 4])
    >>> check_lei('6shgi4...bb395', display_warnings=False)
    (-1, [2, 3, 4])

    Try faulty checksum digits (the final two characters) in a syntax-valid LEI
    
    >>> check_lei('6SHGI4ZSSLCXXQSBB322', display_warnings=False)
    (25, [])
    >>> (1==check_lei('6SHGI4ZSSLCXXQSBB322', display_warnings=False)[0])
    False

    """
    check = -1
    syntax_errcodes = []
    if (len(lei) > 20):
        syntax_errcodes.append(1)
        if (display_warnings): 
            print('WARNING: LEI value is too long:', lei, len(lei))
    if (len(lei) < 20):
        syntax_errcodes.append(2)
        if (display_warnings): 
            print('WARNING: LEI value is too short:', lei, len(lei))
    if (lei.upper() != lei):
        syntax_errcodes.append(3)
        if (display_warnings): 
            print('WARNING: LEI value is not uppercase:', lei)
    if (None==re.search(r'^[0-9A-Z]{18}[0-9]{2}$', lei)):
        syntax_errcodes.append(4)
        if (display_warnings): 
            print('WARNING: LEI does not match the official format:', lei)
    if (lei[len(lei)-2:len(lei)] in ['00', '01', '99']):
        syntax_errcodes.append(5)
        if (display_warnings): 
            print('WARNING: LEI checkdigits in invalid range [00,01,99]:', lei)
    # If no syntax flaws are detected, go ahead and perform the checksum
    if (0==len(syntax_errcodes)):
        ASC_A = ord('A')
        buff = ''
        for i in range(len(lei)):
            char_i = lei[i]
            asc_i = ord(char_i)
            if (asc_i >= ASC_A):
                # Convert letters to their numeric index: A=10, B=11, C=12, etc.
                buff = buff + str(asc_i-55)
            else:
                # Convert digits (0-9) to their string equivalent
                buff = buff + str(char_i)
        # A valid (digitized) LEI should have modulus 1:
        check = int(buff) % 97
        if (display_warnings and 1!=check): 
            print('WARNING: LEI fails checksum:', lei, "modulus="+str(check))
    return (check, syntax_errcodes)



if __name__ == "__main__":
    import doctest
    doctest.testmod()