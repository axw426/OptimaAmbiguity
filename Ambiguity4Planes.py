import HelperFunctions as hf
import random
import math
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F
import time

hf.SetSeed(2022)

readExternal=False

NLoops=500
size=10 #cm
tolerance=2.5
effTolerance=2.0
rawAngle=45
saveStripMaps=False

xmax=size*5000

#gStyle.SetOptStat(0000)
gROOT.SetBatch(True)

if readExternal==True:

        inFile="Run01_2mmCCol_160pA_TrackTrig_SDS_all_events_SYNC.txt"
        #inFile="Run01_2mmCCol_160pA_TrackTrig_SDS__all_events_RTS1_SYNC_Data.txt"
        pitch=90.8
        size=10
        xmax=size*5000
        tolerance=5.01
        rawAngle=60

        #Set up histograms
        MyFile =TFile("test4planes.root","RECREATE");
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
                allHits=hf.FindOverlaps(StripX,StripU,StripV,pitch,rawAngle,tolerance)
                reconstructedHits+=len(allHits)
                #add results to hit map
                for hit in allHits:
                        myHitMap.Fill(hit[0],hit[1])

                #tidy up based on effective pixel area
                allHits=hf.RemoveAdjacentHits(allHits,tolerance,pitch)
                for hit in allHits:
                        myRefinedHitMap.Fill(hit[0],hit[1])

                allHits=hf.RemoveAmbiguities(allHits,rawAngle,pitch)
                for hit in allHits:
                        myRefinedHitMap2.Fill(hit[0],hit[1])



        print ("Reconstructed hits= "+(str)(reconstructedHits))

        myHitMap.Write()        
        myRefinedHitMap.Write()        
        myRefinedHitMap2.Write()        
        MyFile.Close()
else:

        #setup histograms
        ambiguityRate=TH2F("ambiguityRate","Ambiguitie Rate (%) for uniform beam of size "+(str)(size)+"cmX"+(str)(size)+"cm", 5,10.0,110, 9,1.5,10.5)
        ambiguityRate.GetXaxis().SetTitle("Strip Pitch (#mum)")
        ambiguityRate.GetYaxis().SetTitle("nProtons")
        ambiguityRateRefined=TH2F("ambiguityRateRefined","Ambiguitie Rate (%) for uniform beam of size "+(str)(size)+"cmX"+(str)(size)+"cm", 5,10.0,110, 9,1.5,10.5)
        ambiguityRateRefined.GetXaxis().SetTitle("Strip Pitch (#mum)")
        ambiguityRateRefined.GetYaxis().SetTitle("nProtons")


        for nProton in range(13,14):
                for pitch in range (100,101,20):

                        nRawHits=0.0
                        nRefinedHits=0.0

                        nProton_RawEffCorrected=0.0
                        nProton_RefinedEffCorrected=0.0

                        #reset output file and setup histograms for loop
                        #MyFile =TFile("/disk/moose/bilpa/Optima/aw/Ambiguities/nP"+(str)(nProton)+"_pitch"+(str)(pitch)+"um_Array"+(str)(size)+"cm.root","RECREATE");
                        MyFile =TFile("test4planes.root","RECREATE");
                        hRawHits=TH1F("hRawHits","n Reconstructed Hits (Raw)",20,-4.5,24.5)
                        hRawEfficiency=TH1F("hRawEfficiency","Raw Efficiency",10,5,105)

                        hRefinedHits=TH1F("hRefinedHits","n Reconstructed Hits (Refined)",20,-4.5,24.5)
                        hRefinedEfficiency=TH1F("hRefinedEfficiency","Refined Efficiency",10,5,105)

                        myHitMap=TH2F("myHitMap","Hit Locations",(int)(2*xmax/pitch),-xmax,xmax,(int)(2*xmax/pitch),-xmax,xmax)

                        
                        for i in range(NLoops):
                                if i%100 == 0:
                                        print("Processing loop "+(str)(i))

                                #generate N=10 random points for x-y 
                                XY=hf.GetRandomXY(nProton,size)

                                #convert points into a strip number for each layer 
                                StripX,StripU,StripV,StripY=hf.GetStripCoOrds4Planes(XY,pitch,rawAngle)

                                #cycle through all strip combinations, find intersections of pairs of planes, see if they are within tolerance of each other
                                #returns hit objects of form (XCoord,YCoord,stripX,stripU,stripV,[radial distance to intersections])
                                allHits=hf.FindOverlaps4Planes(StripX,StripU,StripV,StripY,pitch,rawAngle,tolerance)
                                nRawHits+=len(allHits)
                                hRawHits.Fill(len(allHits))
                                if saveStripMaps:
                                        hf.PlotHitMap4Planes("RawHits",allHits,XY,StripX,StripU,StripV,StripY,pitch,size,i,rawAngle)
                                rawEff=hf.GetEfficiency(allHits,XY,pitch,effTolerance)
                                hRawEfficiency.Fill(100*rawEff)
                                nProton_RawEffCorrected+=nProton*rawEff

                                for hit in allHits:
                                        myHitMap.Fill(hit[0],hit[1])

                                #remove duplicate hits (separated by < tolerance) 
                                refinedHits=hf.RemoveAdjacentHits(allHits,tolerance,pitch)
                                nRefinedHits+=len(refinedHits)
                                hRefinedHits.Fill(len(refinedHits))
                                if saveStripMaps:
                                        hf.PlotHitMap4Planes("RefinedHits",refinedHits,XY,StripX,StripU,StripV,StripY,pitch,size,i,rawAngle)
                                refinedEff=hf.GetEfficiency(refinedHits,XY,pitch,effTolerance)
                                hRefinedEfficiency.Fill(100*refinedEff)
                                nProton_RefinedEffCorrected+=nProton*refinedEff




                        rawAmbiguity=100*(nRawHits-(NLoops*nProton))/(NLoops*nProton)
                        refinedAmbiguity=100*(nRefinedHits-(NLoops*nProton))/(NLoops*nProton)
                        #rawAmbiguity_Corrected=100*(nRawHits-(nProton_RawEffCorrected))/(nProton_RawEffCorrected)
                        #refinedAmbiguity_Corrected=100*(nRefinedHits-(nProton_RefinedEffCorrected))/(nProton_RefinedEffCorrected)
                        
                        print "nProton="+(str)(nProton)+", pitch="+(str)(pitch)+", rawAmbiguity="+(str)(rawAmbiguity)+"%"+", refinedAmbiguity="+(str)(refinedAmbiguity)
                       # print "Efficiency: raw="+(str)(hRawEfficiency.GetMean())+", refined="+(str)(hRefinedEfficiency.GetMean())
                       # print "Efficiency corrected: rawAmbiguity="+(str)(rawAmbiguity_Corrected)+"%"+", refinedAmbiguity="+(str)(refinedAmbiguity_Corrected)
                        hRawHits.Write()
                        hRawEfficiency.Write()
                        hRefinedHits.Write()
                        hRefinedEfficiency.Write()
                        myHitMap.Write()
                        MyFile.Close()
                        ambiguityRate.Fill(pitch,nProton,rawAmbiguity)
                        ambiguityRateRefined.Fill(pitch,nProton,refinedAmbiguity)

        mycanvas = TCanvas( "c1", "Ambiguities", 550,195,800,700 )
        ambiguityRate.Draw("COLZText")
        mycanvas.SaveAs("outputAmbiguities4Plane_"+(str)(size)+"cm.C")
        mycanvas.SaveAs("outputAmbiguities4Plane_"+(str)(size)+"cm.png")
        ambiguityRateRefined.Draw("COLZTEXT")
        mycanvas.SaveAs("outputAmbiguitiesRefined4Plane_"+(str)(size)+"cm.C")
        mycanvas.SaveAs("outputAmbiguitiesRefined4Plane_"+(str)(size)+"cm.png")
