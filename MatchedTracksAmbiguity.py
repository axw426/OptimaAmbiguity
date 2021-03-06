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
parser.add_argument("-g","--geo",help="Geometry label- see Geometry.py for options, default is 2ModuleXUV",type=str,default='FullSystem')
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

        RearnTracks=[]
        ReartrackEfficiency=[]
        ReartrackAmbiguity=[]
        RearcorrectedTrackAmbiguity=[]

        CombinednTracks=[]
        CombinedtrackEfficiency=[]
        CombinedtrackAmbiguity=[]
        CombinedcorrectedTrackAmbiguity=[]

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
                if i%1 == 0:
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
                mXmY=hf.GetDirections(nProton,beamSpread[0],False)

                ################### FRONT TRACKERS #######################################################

                
                #reconstruct hits for each tracker module
                TrackerHits=[]
                MaxNStrips=[]
                for module in [0,1]:


                        #convert simulated protons to strips 
                        Strips,meanXY=hf.GetTrackerStripCoOrds(XY,mXmY,pitch,TrackerAngles[module],TrackerZ[module])

                        #keep track of maximum number of strips in each module
                        MaxNStrips.append(len(max(Strips,key=len)))                        

                        #find overlap of strips and return hit objects of form (XCoord,YCoord,[stripX,stripU,stripV],[radial distance to intersections])
                        Hits=hf.FindOverlaps(Strips,pitch,TrackerAngles[module],stripTolerance[0],useHalfStrips)
                        TrackerHits.append(Hits)
                        nTrackerHits[module].append(len(Hits))
                        
                        #calculate efficiency for finding correct hits by checking there is always a reco hit within tolerance of true hit
                        eff=hf.GetEfficiency(Hits,meanXY,pitch,stripTolerance[0])
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
                RecoTracks=hf.ReconstructTracks(copy.deepcopy(TrackerHits),trackTolerance[0],pitch,MaxNTracks,restrictNTracks,TrackerZ[0:2],stripTolerance[0],beamSpread[0])

                #hf.WriteTracks(outfilename,RecoTracks,ZMeans,i)
                nTracks.append(len(RecoTracks))
#                print TrackerHits,len(RecoTracks)
                
                #get track efficiency
                eff=hf.GetTrackEfficiency(stripTolerance[0]*pitch,XY,mXmY,RecoTracks,TrackerZ[0:2],pitch)
                trackEfficiency.append(100*eff)

                #calculate track ambiguity rate
                nFakeTracks=len(RecoTracks)-(nProton*eff)
                if len(RecoTracks)>0 :
                        trackAmbiguity.append(100*(len(RecoTracks)-nProton)/len(RecoTracks))
                        correctedTrackAmbiguity.append(100*nFakeTracks/len(RecoTracks))

                if saveStripMaps:
                        hf.DrawTrackMap("TrackMap",RecoTracks,XY,xmax)



                ################### REAR TRACKERS #######################################################

                phantomdepth=30.0 #cm
                energy=230.0 #MeV
                XY_MS,mXmY_MS=hf.ApplyMultipleScattering(XY,mXmY,phantomdepth,energy,(ZMeans[2]+ZMeans[1])/2)
                #XY_MS,mXmY_MS=XY,mXmY
                
                #reconstruct hits for each tracker module
                TrackerHits=[]
                MaxNStrips=[]
                for module in [2,3,4]:


                        #convert simulated protons to strips 
                        Strips,meanXY=hf.GetTrackerStripCoOrds(XY_MS,mXmY_MS,pitch,TrackerAngles[module],TrackerZ[module])

                        #keep track of maximum number of strips in each module
                        MaxNStrips.append(len(max(Strips,key=len)))                        

                        #find overlap of strips and return hit objects of form (XCoord,YCoord,[stripX,stripU,stripV],[radial distance to intersections])
                        Hits=hf.FindOverlaps(Strips,pitch,TrackerAngles[module],stripTolerance[1],useHalfStrips)
                        TrackerHits.append(Hits)
                        nTrackerHits[module].append(len(Hits))
                        
                        #calculate efficiency for finding correct hits by checking there is always a reco hit within tolerance of true hit
                        eff=hf.GetEfficiency(Hits,meanXY,pitch,stripTolerance[1])
                        trackerEffs[module].append(100*eff)
                        nFakeHits=len(Hits)- (nProton*eff)
                        
                        #calculate ambiguity
                        if len(Hits)>0:
                                ambiguity[module].append(100*(float)(len(Hits)-nProton)/len(Hits))
                                correctedAmbiguity[module].append(100*nFakeHits/len(Hits))

                        if saveStripMaps:
                                hf.PlotHitMap("Tracker"+(str)(module)+"Hits",TrackerHits[module],XY_MS,Strips,pitch,size,i,TrackerAngles[module])

                #reconstruct tracks
                RearMaxNTracks=max(MaxNStrips)
                RearRecoTracks=hf.ReconstructTracks(copy.deepcopy(TrackerHits),trackTolerance[1],pitch,MaxNTracks,restrictNTracks,TrackerZ[2:5],stripTolerance[1],beamSpread[1])
                #hf.WriteTracks(outfilename,RecoTracks,ZMeans,i)
                RearnTracks.append(len(RearRecoTracks))

                #get track efficiency
                Reareff=hf.GetTrackEfficiency(stripTolerance[1]*pitch,XY_MS,mXmY_MS,RearRecoTracks,TrackerZ[2:5],pitch)
                ReartrackEfficiency.append(100*Reareff)

                #calculate track ambiguity rate
                RearnFakeTracks=len(RearRecoTracks)-(nProton*Reareff)
                if len(RearRecoTracks)>0 :
                        ReartrackAmbiguity.append(100*(len(RearRecoTracks)-nProton)/len(RearRecoTracks))
                        RearcorrectedTrackAmbiguity.append(100*RearnFakeTracks/len(RearRecoTracks))

                if saveStripMaps:
                        hf.DrawTrackMap("TrackMap",RearRecoTracks,XY_MS,xmax)


                ################################## COMBINE TRACKS #######################################################


                #do matching between front and rear trackers/ maybe do this as PART of the rear module tracking- would boost efficiency.....
                #print len(RecoTracks),len(RearRecoTracks)
                #for track in RecoTracks:
                #        print track
                arbitraryTolerance=10000 
                CombinedRecoTracks=hf.MatchTrack(RecoTracks,RearRecoTracks,ZMeans,arbitraryTolerance)
                CombinednTracks.append(len(CombinedRecoTracks))

                #get track efficiency
                #print len(RecoTracks),len(RearRecoTracks),len(CombinedRecoTracks)
                Combinedeff=hf.GetCombinedEfficiency(stripTolerance,[XY,XY_MS],[mXmY,mXmY_MS],CombinedRecoTracks,ZMeans,pitch)
                CombinedtrackEfficiency.append(100*Combinedeff)

                #calculate track ambiguity rate
                CombinednFakeTracks=len(CombinedRecoTracks)-(nProton*Combinedeff)
                
                
                if len(CombinedRecoTracks)>0 :
                        CombinedtrackAmbiguity.append(100*(len(CombinedRecoTracks)-nProton)/len(CombinedRecoTracks))
                        CombinedcorrectedTrackAmbiguity.append(100*CombinednFakeTracks/len(CombinedRecoTracks))

                
        print "\nResults:\n"
        for module in [0,1]:
                print "Module ",module,": Hits= ",sum(nTrackerHits[module])," Ambiguity= ",np.mean(ambiguity[module]),"%"," Efficiency=",np.mean(trackerEffs[module]),"Purity= ",100-np.mean(correctedAmbiguity[module]),"%"
        print "Front Combined: Tracks= ",sum(nTracks)," Ambiguity= ",np.mean(trackAmbiguity),"%"," Efficiency=",np.mean(trackEfficiency),"Purity= ",100-np.mean(correctedTrackAmbiguity),"%"
        print "Errors: Eff. Error=",np.std(trackEfficiency)/math.sqrt(NLoops),"Purity Err",np.std(correctedTrackAmbiguity)/math.sqrt(NLoops),"%\n"

        for module in [2,3,4]:
                print "Module ",module,": Hits= ",sum(nTrackerHits[module])," Ambiguity= ",np.mean(ambiguity[module]),"%"," Efficiency=",np.mean(trackerEffs[module]),"Purity= ",100-np.mean(correctedAmbiguity[module]),"%"
        print "Rear Combined: Tracks= ",sum(RearnTracks)," Ambiguity= ",np.mean(ReartrackAmbiguity),"%"," Efficiency=",np.mean(ReartrackEfficiency),"Purity= ",100-np.mean(RearcorrectedTrackAmbiguity),"%"
        print "Errors: Eff. Error=",np.std(ReartrackEfficiency)/math.sqrt(NLoops),"Purity Err",np.std(RearcorrectedTrackAmbiguity)/math.sqrt(NLoops),"%\n"

        print "Full Combined: Tracks= ",sum(CombinednTracks)," Ambiguity= ",np.mean(CombinedtrackAmbiguity),"%"," Efficiency=",np.mean(CombinedtrackEfficiency),"Purity= ",100-np.mean(CombinedcorrectedTrackAmbiguity),"%"
        print "Errors: Eff. Error=",np.std(CombinedtrackEfficiency)/math.sqrt(NLoops),"Purity Err",np.std(CombinedcorrectedTrackAmbiguity)/math.sqrt(NLoops),"%\n"

        _Efficiency.append(np.mean(CombinedtrackEfficiency))
        _Purity.append(100-np.mean(CombinedcorrectedTrackAmbiguity))
        _NumberOfProtons.append(nMeanProton)
        _EfficiencyErr.append(np.std(CombinedtrackEfficiency)/math.sqrt(NLoops))
        _PurityErr.append(np.std(CombinedcorrectedTrackAmbiguity)/math.sqrt(NLoops))
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

