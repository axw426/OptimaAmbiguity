import HelperFunctions as hf
import Geometry as geo
import random
import math
import numpy as np
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F
import time

hf.SetSeed(2022)

#beam properties
width=25 #mm #sigma
posX=0.0 #cm
posY=0.0 #cm

NLoops=10
pitch=100.0
size=20 #cm
xmax=size*5000

saveStripMaps=True
useHalfStrips=False

TrackerAngles,TrackerZ,stripTolerance,trackTolerance,effTolerance,pitch,beamSpread=geo.init("2ModuleXUV")
print TrackerAngles,TrackerZ,stripTolerance,trackTolerance,effTolerance,pitch,beamSpread


#gStyle.SetOptStat(0000)
gROOT.SetBatch(True)

for nMeanProton in range(10,11):

        totalProtons=0.0
        nCombinedHits=0.0
        efficiency=0.0
        nTrackerHits=[0.0]*len(TrackerAngles)

        #reset output file and setup histograms for loop
        MyFile =TFile("tracking.root","RECREATE");
        for i in range(NLoops):
                if i%1 == 0:
                        print("Processing loop "+(str)(i))

                #get nProtons according to possion distribution
                nProton=nMeanProton
                #nProton=np.random.poisson(nMeanProton)     
                totalProtons+=nProton

                #get proton starting positions and angular distributions
                XY=hf.GetRandomXY(nProton,size,width,posX,posY)
                mXmY=hf.GetDirections(nProton,beamSpread,False)
                
                TrackerHits=[]
                
                for module in range(len(TrackerAngles)):
                                       
                        #get hits for each set of trackers
                        #returns hit objects of form (XCoord,YCoord,[stripX,stripU,stripV],[radial distance to intersections])
                        Strips=hf.GetTrackerStripCoOrds(XY,mXmY,pitch,TrackerAngles[module],TrackerZ[module])
                        Hits=hf.FindOverlaps(Strips,pitch,TrackerAngles[module],stripTolerance,useHalfStrips)
                        #Hits=hf.RemoveAdjacentHits(Hits,stripTolerance,pitch)
                        TrackerHits.append(Hits)
                        nTrackerHits[module]+=len(Hits)

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
                print "Module ",module,": Hits= ",nTrackerHits[module]," Ambiguity= ",100*(nTrackerHits[module]-totalProtons)/totalProtons,"%"


        print "Combined: Tracks= ",nCombinedHits," Ambiguity= ",100*(nCombinedHits-totalProtons)/totalProtons,"%"

        print("Net Efficiency= ",100*efficiency/NLoops)

