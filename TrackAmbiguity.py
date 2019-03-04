import HelperFunctions as hf
import Geometry as geo
import random
import math
import numpy as np
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F, TMultiGraph, TGraphErrors, TLegend
import time
import argparse
from array import array
import copy

hf.SetSeed(2022)
np.random.seed(53234)

#beam properties
width=25 #mm #sigma
posX=0.0 #cm
posY=0.0 #cm

#read command line options
parser=argparse.ArgumentParser(description="Script for controlling tracker ambiguity studies")
parser.add_argument("-l","--nLoops",help="Number of times NProtons will be fired, default is 100",type=int,default=100)
parser.add_argument("-np","--nProtons",help="Number of protons per frame, default 5",type=int,default=5)
parser.add_argument("-r","--protonRange",help="Range of protons to loop over: starts at nProton, default is 1",type=int,default=1)
parser.add_argument("-g","--geo",help="Geometry label- see Geometry.py for options, default is 2ModuleXUV",type=str,default='2ModuleXUV')
parser.add_argument("-s","--saveHits",help="Save hit maps for visualization (bool), default is 0(=false)",type=int,default=0)
parser.add_argument("-p","--poisson",help="Apply poisson fluctuations to nProtons, default is 1(=true)",type=int,default=1)
parser.add_argument("-split","--useHalfStrips",help="Split Strips In Half, default is 0(=false)",type=int,default=0)
args=parser.parse_args()
parser.print_help()

NLoops=args.nLoops
minProtons=args.nProtons
protonRange=args.protonRange
saveStripMaps=args.saveHits
geoName=args.geo
print args
restrictNTracks=False

if args.useHalfStrips>0:
        useHalfStrips=True
else:
        useHalfStrips=False

TrackerAngles,TrackerZ,ZMeans,stripTolerance,trackTolerance,pitch,beamSpread,size=geo.init(geoName)
#print ZMeans
#print "Using geometry"+geoName,TrackerAngles,TrackerZ,stripTolerance,trackTolerance,pitch,beamSpread,size
xmax=size*5000

#Strip tolerance: arises from spacial separation of XU, XV, and UV overlaps, irrelevant for two strip configuration as only one point of intersection, used for efficiency definition
#Track tolerance: currently defines a quality cut for the track reconstruction based on total deviation from parallel beams

#gStyle.SetOptStat(0000)
gROOT.SetBatch(True)

outfilename="tracking.txt"
f= open(outfilename,"w+")
f.close

_Efficiency,_Purity,_NumberOfProtons=array('d'),array('d'),array('d')
_EfficiencyErr,_PurityErr,_NumberOfProtonsErr=array('d'),array('d'),array('d')


for nMeanProton in range(minProtons,minProtons+protonRange):

        print "\nMean number of protons: ",nMeanProton
        
        totalProtons=0.0

        nTracks=[]
        trackEfficiency=[]
        trackAmbiguity=[]
        correctedTrackAmbiguity=[]

        nTrackerHits=[]
        trackerEffs=[]
        ambiguity=[]
        correctedAmbiguity=[]
        for i in range(len(TrackerAngles)):
                nTrackerHits.append([])
                trackerEffs.append([])
                ambiguity.append([])
                correctedAmbiguity.append([])


                
        #reset output file and setup histograms for loop
        MyFile =TFile("tracking.root","RECREATE");
        for i in range(NLoops):
                if i%10 == 0:
                        print("Processing loop "+(str)(i)+" of "+(str)(NLoops))

                #get nProtons according to possion distribution
                if args.poisson>0:
                        nProton=np.random.poisson(nMeanProton)
                else:
                        nProton=nMeanProton
                totalProtons+=nProton

                if nProton==0:
                        continue
                      
                #get proton starting positions and angular distributions
                #print "Building proton paths",nProton
                XY=hf.GetRandomXY(nProton,size,width,posX,posY)
                mXmY=hf.GetDirections(nProton,beamSpread,False)

       
                #reconstruct hits for each tracker module
                TrackerHits=[]
                MaxNStrips=[]
                for module in range(len(TrackerAngles)):


                        #convert simulated protons to strips 
                        Strips,meanXY=hf.GetTrackerStripCoOrds(XY,mXmY,pitch,TrackerAngles[module],TrackerZ[module])

                        #keep track of maximum number of strips in each module
                        MaxNStrips.append(len(max(Strips,key=len)))                        

                        #find overlap of strips and return hit objects of form (XCoord,YCoord,[stripX,stripU,stripV],[radial distance to intersections])
                        Hits=hf.FindOverlaps(Strips,pitch,TrackerAngles[module],stripTolerance,useHalfStrips)

                        #verify that hits are reconstructed within correct area
                        #Hits=hf.CheckInsideDetectorArea(Hits,pitch,size,TrackerAngles[module])
                        #look for obvious cases of fake hits (adjacent to each other) and merge them by averaging x-y cooords
                        #turned off due to large strip tolerance needed for rear trackers
                        #Hits=hf.MergeAdjacentHits(Hits,stripTolerance,pitch)

                        TrackerHits.append(Hits)
                        nTrackerHits[module].append(len(Hits))
                        
                        #calculate efficiency for finding correct hits by checking there is always a reco hit within tolerance of true hit
                        eff=hf.GetEfficiency(Hits,meanXY,pitch,stripTolerance)
                        trackerEffs[module].append(100*eff)
                        nFakeHits=len(Hits)- (nProton*eff)
                        
                        #calculate ambiguity
                        if len(Hits)>0:
                                ambiguity[module].append(100*(float)(len(Hits)-nProton)/len(Hits))
                                correctedAmbiguity[module].append(100*nFakeHits/len(Hits))

                        if saveStripMaps:
                                hf.PlotHitMap("Tracker"+(str)(module)+"Hits",TrackerHits[module],XY,Strips,pitch,size,i,TrackerAngles[module])

                #reconstruct tracks
                MaxNTracks=max(MaxNStrips)
                if len(TrackerAngles)>1:
                        RecoTracks=hf.ReconstructTracks(copy.deepcopy(TrackerHits),trackTolerance,pitch,MaxNTracks,restrictNTracks,TrackerZ,stripTolerance,beamSpread)
                        hf.WriteTracks(outfilename,RecoTracks,ZMeans,i)
                else:
                        RecoTracks=TrackerHits[0]

                nTracks.append(len(RecoTracks))

                #get track efficiency
                eff=hf.GetTrackEfficiency(stripTolerance*pitch,XY,mXmY,RecoTracks,TrackerZ,pitch)
                trackEfficiency.append(100*eff)

                #calculate track ambiguity rate
                nFakeTracks=len(RecoTracks)-(nProton*eff)
                if len(RecoTracks)>0 :
                        trackAmbiguity.append(100*(len(RecoTracks)-nProton)/len(RecoTracks))
                        correctedTrackAmbiguity.append(100*nFakeTracks/len(RecoTracks))

                if saveStripMaps and len(TrackerAngles)>1:
                        hf.DrawTrackMap("TrackMap",RecoTracks,XY,xmax)

                
        for module in range(len(TrackerAngles)):
                print "Module ",module,": Hits= ",sum(nTrackerHits[module])," Ambiguity= ",np.mean(ambiguity[module]),"%"," Efficiency=",np.mean(trackerEffs[module]),"Efficiency corrected ambiguity= ",np.mean(correctedAmbiguity[module]),"%"

        print "Combined: Tracks= ",sum(nTracks)," Ambiguity= ",np.mean(trackAmbiguity),"%"," Efficiency=",np.mean(trackEfficiency),"Efficiency corrected ambiguity= ",np.mean(correctedTrackAmbiguity),"%"
        print "Errors: Eff. Error=",np.std(trackEfficiency)/math.sqrt(NLoops),"Purity Err",np.std(correctedTrackAmbiguity)/math.sqrt(NLoops),"%"
        
        _Efficiency.append(np.mean(trackEfficiency))
        _Purity.append(100-np.mean(correctedTrackAmbiguity))
        _NumberOfProtons.append(nMeanProton)
        _EfficiencyErr.append(np.std(trackEfficiency)/math.sqrt(NLoops))
        _PurityErr.append(np.std(correctedTrackAmbiguity)/math.sqrt(NLoops))
        _NumberOfProtonsErr.append(0.0)
                       
canvas1 = TCanvas( 'c1', "mycanvas", 200, 10, 700, 500 )

effGraph=TGraphErrors(len(_NumberOfProtons),_NumberOfProtons,_Efficiency,_NumberOfProtonsErr,_EfficiencyErr)
effGraph.SetMarkerStyle(2)
effGraph.SetMarkerSize(2)
effGraph.GetXaxis().SetTitle("Mean N Proton")
effGraph.GetYaxis().SetTitle("Efficiency/Purity (%)")
effGraph.SetTitle("Efficiency")

purGraph=TGraphErrors(len(_NumberOfProtons),_NumberOfProtons,_Purity,_NumberOfProtonsErr,_PurityErr)
purGraph.SetMarkerStyle(5)
purGraph.SetLineColor(2)
purGraph.SetMarkerColor(2)
purGraph.SetMarkerSize(2)
purGraph.GetXaxis().SetTitle("Mean N Proton")
purGraph.GetYaxis().SetTitle("Purity")
purGraph.SetTitle("Purity")
MyFile =TFile("effVsPurity_bs"+(str)(beamSpread)+"mrad_geo"+geoName+"_RestrictTracks"+(str)(restrictNTracks)+"_pitch"+(str)(pitch)+".root","RECREATE");

mg=TMultiGraph()
mg.SetTitle("Efficiency and Purity vs nProtons;Mean N Protons; (%)")
mg.Add(effGraph,"lp")
mg.Add(purGraph,"lp")
mg.Draw("a")
effGraph.Write()
purGraph.Write()
mg.Write()


legend=TLegend(0.1,0.7,0.48,0.9)
legend.AddEntry(effGraph,"Efficiency","lp")
legend.AddEntry(purGraph,"Purity","lp")

legend.Draw("SAME")
canvas1.Write("EffVsPurity")
