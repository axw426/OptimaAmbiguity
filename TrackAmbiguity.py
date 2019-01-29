import HelperFunctions as hf
import random
import math
import numpy as np
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F
import time

hf.SetSeed(2022)

NLoops=100
size=10 #cm
xmax=size*5000
stripTolerance=0.7
pitch=100.0

saveStripMaps=True
useHalfStrips=False

beamSpread=3.0 #mrad

sensorThickness=155
interPlaneDistance=12*1000
interModuleDistance=100*1000
phantomGap=80*1000

planesPerModule=3

if planesPerModule==2: 
        totalDistance=4*interPlaneDistance+2*interModuleDistance+phantomGap

        #trackTolerance=2*math.tan(beamSpread/1000.0)*totalDistance +500 #4 trackers
        trackTolerance=10 #2 trackers
        pos=0

        Tracker1Angles=[0,90]
        Tracker1Z=[pos,pos+interPlaneDistance]
        
        pos+=interModuleDistance+interPlaneDistance
        
        Tracker2Angles=[45,135]
        Tracker2Z=[pos,pos+interPlaneDistance]

        pos+=phantomGap+interModuleDistance
        
        Tracker3Angles=[0,90]
        Tracker3Z=[pos,pos+interPlaneDistance]

        pos+=interModuleDistance+interPlaneDistance

        Tracker4Angles=[45,135]
        Tracker4Z=[pos,pos+interPlaneDistance]

elif planesPerModule==3: 
        totalDistance=4*interPlaneDistance+2*interModuleDistance+phantomGap

        #trackTolerance=2*math.tan(beamSpread/1000.0)*totalDistance +500 #4 trackers
        trackTolerance=stripTolerance*pitch #2 trackers
        pos=0

        Tracker1Angles=[0,60,120]
        Tracker1Z=[pos,pos,pos]
        
        pos+=interModuleDistance+2*interPlaneDistance
        
        Tracker2Angles=[0,60,120]
        Tracker2Z=[pos,pos+interPlaneDistance,pos+2*interPlaneDistance]

        pos+=phantomGap+interModuleDistance
        
        Tracker3Angles=[0,60,120]
        Tracker3Z=[pos,pos+interPlaneDistance]

        pos+=interModuleDistance+interPlaneDistance

        Tracker4Angles=[0,60,120]
        Tracker4Z=[pos,pos+interPlaneDistance]        

#gStyle.SetOptStat(0000)
gROOT.SetBatch(True)

for nMeanProton in range(10,11):

        totalProtons=0.0
        nTracker1Hits=0.0
        nTracker2Hits=0.0
        nTracker3Hits=0.0
        nTracker4Hits=0.0
        nCombinedHits=0.0
        
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
                XY=hf.GetRandomXY(nProton,size)
                mXmY=hf.GetDirections(nProton,beamSpread,False)
                
                
                #get hits for each set of trackers
                #returns hit objects of form (XCoord,YCoord,[stripX,stripU,stripV],[radial distance to intersections])
                Strips1=hf.GetTrackerStripCoOrds(XY,mXmY,pitch,Tracker1Angles,Tracker1Z)
                Tracker1Hits=hf.FindOverlaps(Strips1,pitch,Tracker1Angles,stripTolerance,useHalfStrips)
                #Tracker1Hits=hf.RemoveAdjacentHits(Tracker1Hits,1.0,pitch)
                if saveStripMaps:
                        hf.PlotHitMap("Tracker1Hits",Tracker1Hits,XY,Strips1,pitch,size,i,Tracker1Angles)
                nTracker1Hits+=len(Tracker1Hits)
                
                Strips2=hf.GetTrackerStripCoOrds(XY,mXmY,pitch,Tracker2Angles,Tracker2Z)
                Tracker2Hits=hf.FindOverlaps(Strips2,pitch,Tracker2Angles,stripTolerance,useHalfStrips)
                #Tracker2Hits=hf.RemoveAdjacentHits(Tracker2Hits,1.0,pitch)
                if saveStripMaps:
                        hf.PlotHitMap("Tracker2Hits",Tracker2Hits,XY,Strips2,pitch,size,i,Tracker2Angles)
                nTracker2Hits+=len(Tracker2Hits)

                Strips3=hf.GetTrackerStripCoOrds(XY,mXmY,pitch,Tracker3Angles,Tracker3Z)
                Tracker3Hits=hf.FindOverlaps(Strips3,pitch,Tracker3Angles,stripTolerance,useHalfStrips)
                #Tracker3Hits=hf.RemoveAdjacentHits(Tracker3Hits,1.0,pitch)
                if saveStripMaps:
                        hf.PlotHitMap("Tracker3Hits",Tracker3Hits,XY,Strips3,pitch,size,i,Tracker3Angles)
                nTracker3Hits+=len(Tracker3Hits)
                
                Strips4=hf.GetTrackerStripCoOrds(XY,mXmY,pitch,Tracker4Angles,Tracker4Z)
                Tracker4Hits=hf.FindOverlaps(Strips4,pitch,Tracker4Angles,stripTolerance,useHalfStrips)
                #Tracker4Hits=hf.RemoveAdjacentHits(Tracker4Hits,1.0,pitch)
                if saveStripMaps:
                        hf.PlotHitMap("Tracker4Hits",Tracker4Hits,XY,Strips4,pitch,size,i,Tracker4Angles)
                nTracker4Hits+=len(Tracker4Hits)

                #RecoTracks=hf.ReconstructTracks([Tracker1Hits,Tracker2Hits,Tracker3Hits,Tracker4Hits],trackTolerance,pitch)
                RecoTracks=hf.ReconstructTracks([Tracker1Hits,Tracker2Hits],trackTolerance,pitch)
                nCombinedHits+=len(RecoTracks)
                if saveStripMaps:
                        hf.DrawTrackMap("TrackMap",RecoTracks,XY,xmax)


        Track1Ambiguity=100*(nTracker1Hits-totalProtons)/totalProtons
        Track2Ambiguity=100*(nTracker2Hits-totalProtons)/totalProtons
        Track3Ambiguity=100*(nTracker3Hits-totalProtons)/totalProtons
        Track4Ambiguity=100*(nTracker4Hits-totalProtons)/totalProtons
        CombinedAmbiguity=100*(nCombinedHits-totalProtons)/totalProtons

        print("Hits: Tracker1=",nTracker1Hits,"Tracker2=",nTracker2Hits,"Tracker3=",nTracker3Hits,"Tracker4=",nTracker4Hits,"Combined",nCombinedHits)
        print("Ambiguities(%): Tracker1=",Track1Ambiguity,"Tracker2=",Track2Ambiguity,"Tracker3=",Track3Ambiguity,"Tracker4=",Track4Ambiguity,"Combined",CombinedAmbiguity)
