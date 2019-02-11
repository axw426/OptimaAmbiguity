import HelperFunctions as hf
import random
import numpy as np
import math
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F
import time

hf.SetSeed(2022)

NLoops=2000
size=10 #cm
tolerance=0.7
effTolerance=2.0
angles=[0,60,120]
saveStripMaps=False

xmax=size*5000

#gStyle.SetOptStat(0000)
gROOT.SetBatch(True)


inFile="Run01_2mmCCol_160pA_TrackTrig_SDS_all_events_SYNC.txt"
#inFile="Run01_2mmCCol_160pA_TrackTrig_SDS__all_events_RTS1_SYNC_Data.txt"
pitch=90.8
size=10
xmax=size*5000
tolerance=5.01
angles=[0,60,-60]

#Set up histograms
MyFile =TFile("test.root","RECREATE");
binscale=1
myHitMap=TH2F("myHitMap","Hit Locations",(int)(2*xmax/pitch)/binscale,-xmax,xmax,(int)(2*xmax/pitch)/binscale,-xmax,xmax)
myRefinedHitMap=TH2F("myRefinedHitMap","Hit Locations",(int)(2*xmax/pitch)/binscale,-xmax,xmax,(int)(2*xmax/pitch)/binscale,-xmax,xmax)
myRefinedHitMap2=TH2F("myRefinedHitMap2","Hit Locations",(int)(2*xmax/pitch)/binscale,-xmax,xmax,(int)(2*xmax/pitch)/binscale,-xmax,xmax)

#read in hits, grouped by timestamp. Expects rosetta stile formatting: time, v1,v2,v3,v4, x1,x2,x3,x4, u1,u2,u3,u4
hitsByTimestamp=hf.ReadXUVStripCoOrds(inFile)
print("Timestamps provided "+(str)(len(hitsByTimestamp)))

#loop over all timestamps and find overlaps
reconstructedHits=0
for i in range(len(hitsByTimestamp)):
        StripX=hitsByTimestamp[i][0]
        StripU=hitsByTimestamp[i][1]
        StripV=hitsByTimestamp[i][2]
        
        #print(StripX,StripU,StripV)
        #find hits
        allHits=hf.FindOverlaps([StripX,StripU,StripV,StripX,StripU,StripV],pitch,angles,tolerance,False) #strips X,U,V doubled here as find overlaps expects info on which half of a plane a sensor is in...
        reconstructedHits+=len(allHits)
        #add results to hit map
        for hit in allHits:
                myHitMap.Fill(hit[0],hit[1])
                
        #tidy up based on effective pixel area
        allHits=hf.RemoveAdjacentHits(allHits,tolerance,pitch)
        for hit in allHits:
                myRefinedHitMap.Fill(hit[0],hit[1])

        allHits=hf.RemoveAmbiguities(allHits,angles,pitch)
        for hit in allHits:
                myRefinedHitMap2.Fill(hit[0],hit[1])



print ("Reconstructed hits= "+(str)(reconstructedHits))

myHitMap.Write()        
myRefinedHitMap.Write()        
myRefinedHitMap2.Write()        
MyFile.Close()
