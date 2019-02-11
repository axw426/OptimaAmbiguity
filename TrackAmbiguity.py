import HelperFunctions as hf
import Geometry as geo
import random
import math
import numpy as np
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F
import time
import argparse

hf.SetSeed(2022)
np.random.seed(53234)

#beam properties
width=25 #mm #sigma
posX=0.0 #cm
posY=0.0 #cm

#NLoops=10
pitch=100.0
size=20 #cm
xmax=size*5000

#read command line options
parser=argparse.ArgumentParser(description="Script for controlling tracker ambiguity studies")
parser.add_argument("-nLoops","--nLoops",help="Number of times NProtons will be fired, default is 100",type=int,default=100)
parser.add_argument("-nProtons","--nProtons",help="Number of protons per frame, default 5",type=int,default=5)
parser.add_argument("-protonRange","--protonRange",help="Range of protons to loop over: starts at nProton, default is 1",type=int,default=1)
parser.add_argument("-geo","--geo",help="Geometry label- see Geometry.py for options, default is 2ModuleXUV",type=str,default='2ModuleXUV')
parser.add_argument("-saveHits","--saveHits",help="Save hit maps for visualization (bool), default is 0(=false)",type=int,default=0)
parser.add_argument("-poisson","--poisson",help="Apply poisson fluctuations to nProtons, default is 1(=true)",type=int,default=1)
args=parser.parse_args()
parser.print_help()

NLoops=args.nLoops
minProtons=args.nProtons
protonRange=args.protonRange
saveStripMaps=args.saveHits
geoName=args.geo
print args

useHalfStrips=False

TrackerAngles,TrackerZ,stripTolerance,trackTolerance,effTolerance,pitch,beamSpread=geo.init(geoName)
print "Using geometry"+geoName,TrackerAngles,TrackerZ,stripTolerance,trackTolerance,effTolerance,pitch,beamSpread

#Strip tolerance: arises from spacial separation of XU, XV, and UV overlaps, irrelevant for two strip configuration as only one point of intersection
#Efficiency tolerance: accounts for difference between reconstructed and true hit positions- essentially combination of hit efficiency and ambiguity in association of hits to form tracks
#Track tolerance: currently defines a quality cut for the track reconstruction based on total deviation from parallel beams

#gStyle.SetOptStat(0000)
gROOT.SetBatch(True)

for nMeanProton in range(minProtons,minProtons+protonRange):

        totalProtons=0.0
        nCombinedHits=0.0
        efficiency=0.0
        nTrackerHits=[0.0]*len(TrackerAngles)
        trackerEffs=[[]]*len(TrackerAngles)

        #reset output file and setup histograms for loop
        MyFile =TFile("tracking.root","RECREATE");
        for i in range(NLoops):
                if i%1 == 0:
                        print("Processing loop "+(str)(i))

                #get nProtons according to possion distribution
                if args.poisson>0:
                        nProton=np.random.poisson(nMeanProton)
                else:
                        nProton=nMeanProton
                totalProtons+=nProton
                      
                #get proton starting positions and angular distributions
                XY=hf.GetRandomXY(nProton,size,width,posX,posY)
                mXmY=hf.GetDirections(nProton,beamSpread,False)
                
                TrackerHits=[]
                
                for module in range(len(TrackerAngles)):
                                       
                        #get hits for each set of trackers
                        #returns hit objects of form (XCoord,YCoord,[stripX,stripU,stripV],[radial distance to intersections])
                        Strips,meanXY=hf.GetTrackerStripCoOrds(XY,mXmY,pitch,TrackerAngles[module],TrackerZ[module])
                        Hits=hf.FindOverlaps(Strips,pitch,TrackerAngles[module],stripTolerance,useHalfStrips)
                        Hits=hf.MergeAdjacentHits(Hits,stripTolerance,pitch)
                        TrackerHits.append(Hits)
                        nTrackerHits[module]+=len(Hits)
                        trackerEffs[module].append(hf.GetEfficiency(Hits,meanXY,pitch,stripTolerance))

                        #print XY
                        if saveStripMaps:
                                hf.PlotHitMap("Tracker"+(str)(module)+"Hits",TrackerHits[module],XY,Strips,pitch,size,i,TrackerAngles[module])

                        #for hit in Tracker1Hits:
                        #        area=hf.GetPixelArea(hit,Tracker1Angles,pitch)
                        #        print "Area=",area,hit
                
                #reconstruct tracks
                if len(TrackerAngles)>1:
                        RecoTracks=hf.ReconstructTracks(TrackerHits,trackTolerance,pitch)
                else:
                        RecoTracks=TrackerHits[0]
                        
                nCombinedHits+=len(RecoTracks)
                if saveStripMaps and len(TrackerAngles)>1:
                        hf.DrawTrackMap("TrackMap",RecoTracks,XY,xmax)

                efficiency+=hf.GetTrackEfficiency(effTolerance,XY,mXmY,RecoTracks,TrackerZ,pitch)
                
        for module in range(len(TrackerAngles)):
                print "Module ",module,": Hits= ",nTrackerHits[module]," Ambiguity= ",100*(nTrackerHits[module]-totalProtons)/totalProtons,"%"," Efficiency=",100*sum(trackerEffs[module])/len(trackerEffs[module])


        print "Combined: Tracks= ",nCombinedHits," Ambiguity= ",100*(nCombinedHits-totalProtons)/totalProtons,"%", "Efficiency=",100*efficiency/NLoops


