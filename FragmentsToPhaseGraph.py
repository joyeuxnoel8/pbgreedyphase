#!/usr/bin/env python
import networkx as nx
import sys
import pickle
import argparse
import intervaltree
import bisect
import copy

from FragmentsToPhaseGraphUtils import *

ap = argparse.ArgumentParser(description="Make a graph from a .fragments file")
ap.add_argument("fragments", help=".fragmnes file: read n chrom [tuple#1... tuple#n], where tuple is in the format pos,ref,alt,read,pre,suf")

ap.add_argument("--minBlock", help="Minimum prefix or suffix to consider a tuple trustworthy.", type=int, default=5)
ap.add_argument("--out", help="Output file.", default=None)
ap.add_argument("--minWeight", help="Only output edges with this weight.", type=int,default=5)
ap.add_argument("--minAlleleCov", help="Only count vertices with this allele count.", type=int, default=10)
ap.add_argument("--writeAlleleCount", help="Write TSV of allele counts to this file", default=None)
ap.add_argument("--writef", help="Write fragments to a file", default=None)
ap.add_argument("--writeog", help="Write overlap graph.", default=None)
ap.add_argument("--readog", help="Read overlap graph.", default=None)
ap.add_argument("--writefg", help="Write the fragment graph to a file", default=None)
ap.add_argument("--readfg", help="Read the fragment graph from a file", default=None)
ap.add_argument("--minOverlap", help="Minimum fragment overlap", default=3,type=int)
ap.add_argument("--minAlleleFreq", help="Minimum allele frequency to keep after fixing", default=0.20, type=float)
ap.add_argument("--minEdgeWeight", help="Minimum edge weight",
    default=4, type=int)
ap.add_argument("--og", help="Write overlap graph", default=None)
ap.add_argument("--pg", help="Write phase graph", default=None)
ap.add_argument("--vcf", help="Write phased vcf to this file.", default=None)
ap.add_argument("--contig", help="Use this contig in the vcf.  In the tuple \"name length\"", default=None)
args = ap.parse_args()

fragmentFile = open(args.fragments)

if (args.out is not None):
    outFile = open(args.out,'w')

g = nx.Graph()

if (args.vcf is not None):
    if (args.contig is None):
        print "ERROR.  When generating a VCF, the name and length of the original contig must be specified using the option --contig"
        sys.exit(1)

totalCount = {}


if (args.readfg is not None):
    print "reading pickle"
    dataFile = open(args.readfg, 'rb')
    frags = pickle.load(dataFile)
    g = pickle.load(dataFile)
    fc = pickle.load(dataFile)
    positions = pickle.load(dataFile)
    
else:
    frags = [ParseFragLine(line) for line in fragmentFile]
    MergeFragments(frags)
    refCount = {}
    altCount = {}

    allPositions = {}
                
    for f in frags:
        for t in f.tuples:
            allPositions[t.pos] = True
            if (t.allele != 2):
                if (g.has_node(t.GetNode()) == False):
                    g.add_node(t.GetNode())
    
                if (t.allele == 0):
                    IncCount(refCount, t.pos)
                if (t.allele == 1):
                    IncCount(altCount, t.pos)
                
                IncCount(totalCount, t.pos)

    #
    # Create a graph with each vertex an SNV allele, and each edge if
    # two alleles are connected by a read.
    #
    import pdb
                    
    for f in frags:
        for i in range(0,len(f.tuples)-1):
            for j in range(i+1,len(f.tuples)):
                if (Supported(f.tuples[i], args.minBlock) and Supported(f.tuples[j], args.minBlock)):
                    s = f.tuples[i].GetNode()
                    t = f.tuples[j].GetNode()
                    if (s is not None and t is not None):
                        if (g.has_edge(s,t) == False):
                            g.add_edge(s,t, weight=0)
                        g[s][t]['weight'] +=1


    toRemove = []

    for (s,t,w) in g.edges(data='weight'):
        if (w < args.minWeight):
            toRemove.append((s,t))

    g.remove_edges_from(toRemove)
    
    toRemove = []
    for n in g.nodes():
        if (g.degree(n) == 0):
            toRemove.append(n)
            
    g.remove_nodes_from(toRemove)
    
    nodesToRemove = []
    posToRemove = []
    positions = []
    #
    # Generate a list of potential heterozygous positions that can be
    # phased.
    #
    for c in totalCount.keys():
        refKey = str(c) + ".0"
        altKey = str(c) + ".1"
        if (c not in refCount and c in altCount):
            toRemove.append(altKey)
        if (c not in altCount and c in refCount):
            toRemove.append(refKey)
    
        if (c in refCount and c in altCount):
            if (refCount[c] < args.minAlleleCov or altCount[c] < args.minAlleleCov):
                toRemove.append(refKey)
                toRemove.append(altKey)
                posToRemove.append(c)
            else:
                positions.append(int(c))
    

    print "Removing " + str(len(toRemove)) + " nodes from the node count from support issues"
    g.remove_nodes_from(toRemove)
                
    nodes = g.nodes()
    fc = [ (frags[i].NInformative(nodes), i, frags[i].read )for i in range(0,len(frags))]
    
    nfc = len(fc)

    
    
    posToRemove.sort()
    ap = allPositions.keys()
    for p in posToRemove:
        if (p in allPositions):
            del allPositions[p]

    pd = {p:True for p in positions}
    mp = {}
    for a in allPositions.keys():
        if a not in pd:
            mp[a] = True
    notRepresented = mp.keys()
    notRepresented.sort()
    RemoveUninformative(posToRemove, frags)
    RemoveUninformative(notRepresented, frags)    
    



fc.sort(reverse=True)


(ref,alt,refNuc,altNuc) = StoreAlleleCount(frags, positions)

if (args.writeAlleleCount is not None):
    acFile = open(args.writeAlleleCount, 'w')
    for i in range(0,len(ref)):
        acFile.write("{}\t{}\t{}\t{}\t{}\n".format(positions[i], ref[i], alt[i], refNuc[i], altNuc[i]))
    acFile.close()

#
# Process overlaps.
#

if (args.readog is not None):
    ogFile = open(args.readog, 'rb')
    og = pickle.load(ogFile)
    ogFile.close()
else:
    og = BuildOverlapGraph(frags, positions, args.minOverlap)


#
# Use the overlap graph to count the support per tuple.
#
if (args.readfg is None):
    pos = { positions[i]:i for i in range(0,len(positions))}
    AddOverlapSupport(og, frags, pos)

    #
    # not super necessary, but remove fragments that do not
    # overlap any existing snv
    #
    RemoveNotOverlappingFragments(frags, og)

#
# Write out some of the results for faster re-runs in debugging
#
    
if (args.writefg is not None):
    writeFg = open(args.writefg,'wb')
    pickle.dump(frags, writeFg, pickle.HIGHEST_PROTOCOL)
    pickle.dump(g, writeFg, pickle.HIGHEST_PROTOCOL)
    pickle.dump(fc, writeFg, pickle.HIGHEST_PROTOCOL)
    positions.sort()
    pickle.dump(positions, writeFg, pickle.HIGHEST_PROTOCOL)

   
if (args.writeog is not None):
    ogFile = open(args.writeog, 'wb')
    pickle.dump(og, ogFile, pickle.HIGHEST_PROTOCOL)
    ogFile.close()


if (args.writef is not None):
    outf = open(args.writef, 'wb')
    pickle.dump(frags, outf)
    outf.close()


positions.sort()
pg = MakePG(positions, frags)
nx.write_gexf(pg, "before_fixing.gexf")
nFixed = FixFragments(frags, 0.6)

(ref,alt)= StoreFrequency(positions, frags)


nRemoved = FilterHomozygousSites(positions, ref, alt, frags,
                                args.minAlleleFreq)
pg = MakePG(positions, frags)
nx.write_gexf(pg, "after_fixing.gexf")

# Try a second round.
ClearFragmentSupport(frags)
og = BuildOverlapGraph(frags, positions, args.minOverlap)
pos = { positions[i]:i for i in range(0,len(positions))}
AddOverlapSupport(og, frags, pos)
nFixed = FixFragments(frags, 0.6)
pg = MakePG(positions, frags)
nx.write_gexf(pg, "after_fixing.2.gexf")

PrintWeights(positions, pg)
nLowWeight = RemovePairedLowWeightEdges(positions, pg, args.minEdgeWeight)


nx.write_gexf(pg, "after_fixing.3.gexf")

#PrintDegree(positions, pg)

uninformative = FindUninformative(positions, pg)
nRemoved = RemoveUninformative(uninformative, frags)
RemovePositions(positions, uninformative)
pg = MakePG(positions, frags)


uninformative = FindUninformative(positions, pg, strict=True)
nRemoved = RemoveUninformative(uninformative, frags)
RemovePositions(positions, uninformative)
pg = MakePG(positions, frags)

RemovePairedLowWeightEdges(positions, pg)

lowWeight = DetectLowWeightSites(positions, pg, minWeight=10, maxFraction=0.25)
nRemoved = RemoveUninformative(lowWeight, frags)
RemovePositions(positions, lowWeight)
pg = MakePG(positions, frags)

PrintWeights(positions,pg)
# Attempt to phase past uninformative sites

uninformative = FindUninformative(positions, pg, strict=True)
i=0


def GetComp(node, comps):

    for i in range(0,len(comps)):
        if (node in comps[i]):
            return i
    return None

def GetMinWeight(dg, s, t):
    ug = dg.to_undirected()
    p = nx.shortest_path(ug,s,t)
    if (p is None):
        return None
    else:
        minWeight = None
        for i in range(0,len(p)-1):
            w = ug[p[i]][p[i+1]]['weight']
            if ( minWeight == None or w < minWeight):
                minWeight = w
        return minWeight
    
while (i < len(uninformative)):
    j=i+1

    while (j < len(uninformative) and uninformative[j] == uninformative[j-1]):
        j+=1

    ustart = positions.index(uninformative[i])
    uend   = positions.index(uninformative[j-1])

    if (ustart > 0 and uend < len(positions)-1):
        window = 3
        pstart = max(ustart-window, 0)
        pend   = min(uend+window+1,len(positions))

        psubset = positions[pstart:pend]
        for u in uninformative[i:j]:
            psubset.remove(u)

        fpg = MakeFullPG(psubset, frags)

        #
        # Now remove low coverage edges.
        #
        RemoveLowWeightEdges(fpg,3)
        
        comps = [c for c in nx.weakly_connected_components(fpg)]

        lens = [len(c) for  c in comps]

        if (len(comps) == 2):
            #
            # It may be possible to fix the graph.
            #
            if (lens[0] == lens[1] and lens[0] == len(psubset)):
                #
                # The two compoenents must have the same length, and each must have
                # the same as the number of positions in the subset
                #

                # b is for bridge
                bstart = ustart-1
                bend   = uend+1
                if (bstart >= 0 and bend < len(positions)):
                    #
                    # The position before the uninformative and position after uninformative are
                    #  in the component. Find which component they are in.
                    #
                    bstart0 = str(positions[bstart]) + ".0"
                    bstart1 = str(positions[bstart]) + ".1"
                    bend0   = str(positions[bend]) + ".0"
                    bend1   = str(positions[bend]) + ".1"

                    bstart0Comp = GetComp(bstart0, comps)
                    bstart1Comp = GetComp(bstart1, comps)
                    bend0Comp   = GetComp(bend0, comps)
                    bend1Comp   = GetComp(bend1, comps)
                    #
                    # Make sure the endpoints are in different components

                    if (bstart0Comp != bstart1Comp and bend0Comp != bend1Comp):
                        if (bstart0Comp == bend0Comp):
                            # add an edge
                            w = GetMinWeight(fpg, bstart0,bend0)
                            if (w is not None):
                                pg.add_edge(bstart0, bend0, weight =2)
                        elif (bstart0Comp == bend1Comp):
                            w = GetMinWeight(fpg, bstart0,bend1)
                            if (w is not None):
                                pg.add_edge(bstart0, bend1, weight=w)
                        if (bstart1Comp == bend0Comp):
                            w = GetMinWeight(fpg, bstart1,bend0)
                            if (w is not None):
                                pg.add_edge(bstart1, bend0, weight=w)
                        elif (bstart1Comp == bend1Comp):
                            w = GetMinWeight(fpg, bstart1,bend1)
                            if (w is not None):
                                pg.add_edge(bstart1, bend1, weight=w)

                        # Clear out the uninformative vertices
                        toRemove = positions[ustart:uend+1]
                        for u in toRemove:
                            for phase in [".0", ".1"]:
                                node = str(u) + phase
                                for e in pg.edges(node):
                                    pg.remove_edge(*e)
                                pg.remove_node(str(u)+phase)
                                if (pg.has_edge(bstart0, node)):
                                    pb.remove_edge(bstart0,node)
                                if (pg.has_edge(bstart1,node)):
                                    pg.remove_edge(bstart1,node)
                            positions.remove(u)
                    
        
    
    # move to next pos
    i=j

PrintWeights(positions,pg)
nx.write_gexf(pg, "after_fixing.4.gexf")
print "Ended with " + str(len(positions)) + " sites."


if (args.pg is not None):
    WriteNx(pg, args.pg )
        
        
if (args.og):
    WriteNx(og, args.og)

if (args.out is not None):
    WriteNx(g,args.out)

def GetAllele(node):
    return int(node[-1])

if (args.vcf is not None):

    (ref,alt,refNuc, altNuc) = StoreAlleleCount(frags, positions)

    contigFai = open(args.contig + ".fai")
    line = contigFai.readline()
    fai = line.split()
    contig = [fai[0], int(fai[1])]
    vcf = open(args.vcf, 'w')

    vcf.write("##fileformat=VCFv4.1\n")
    vcf.write("##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observations\">\n")
    vcf.write("##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observations\">\n")
    vcf.write("##contig=<ID={},length={}>\n".format(contig[0], contig[1]))
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")

    hap=0
    
    for i in range(0,len(positions)-1):
        n0 = str(positions[i]) + ".0"
        n1 = str(positions[i]) + ".1"
        if (hap==0):
            genotype = "0|1"
        else:
            genotype = "1|0"

        infoStr="AN=2;AO={};RO={}".format(alt[i],ref[i])
        #
        # VCF is 1-delimited, so offset alignment
        vcf.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\tGT\t{}\n".format(contig[0],
                                                                    positions[i]+1, refNuc[i],altNuc[i],infoStr, genotype))
        if pg.out_degree(n0) == 1 and pg.out_degree(n1) == 1:
            curNode = str(positions[i]) + "." + str(hap)
            nextNode = pg.edges(curNode)[0][1]
            hap = GetAllele(nextNode)
        else:
            print "Resetting haplotype at " + str(positions[i])
            hap = 0
            
