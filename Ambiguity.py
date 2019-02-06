import HelperFunctions as hf
import random
import numpy as np
import math
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F, TGraph, TGraphErrors
import time
from array import array

hf.SetSeed(2022)
np.random.seed(53234)
saveStripMaps=False
useHalfStrips=False

NLoops=1000
size=20 #cm

#beam properties
width=25 #mm
posX=0.0 #cm
posY=0.0 #cm

#X-Y
#tolerance=0.0 #no need for tolerance as X-Y guaranteed to overlap
#angles=[0,90]

#X-U-V
tolerance=0.67 #~pitch*2/3 (triangle innit), area=5774.0 for 100um
angles=[0,60,120]

#X-U-V-Y
#tolerance=1.5 #needed for 4 planes, why?
#angles=[0,45,90,135]

effTolerance=2.5


xmax=size*5000

#gStyle.SetOptStat(0000)
gROOT.SetBatch(True)

_aMeanProton,_aAmbiguity=array('d'),array('d')
_aMeanProtonError,_aAmbiguityError=array('d'),array('d')

#setup histograms
ambiguityRate=TH2F("ambiguityRate","Ambiguitie Rate (%) for uniform beam of width "+(str)(width)+"cm", 5,10.0,110, 13,0.5,13.5)
ambiguityRate.GetXaxis().SetTitle("Strip Pitch (#mum)")
ambiguityRate.GetYaxis().SetTitle("nProtons")
ambiguityRateRefined=TH2F("ambiguityRateRefined","Ambiguitie Rate (%) for uniform beam of width "+(str)(width)+"cm", 5,10.0,110, 13,0.5,13.5)
ambiguityRateRefined.GetXaxis().SetTitle("Strip Pitch (#mum)")
ambiguityRateRefined.GetYaxis().SetTitle("nProtons")
hNProtons=TH1F("hNProtons","",30,-0.5,29.5)

for nMeanProton in range(10,11):
        for pitch in range (100,101,20):

                nRawHits=0.0
                nRefinedHits=0.0
                totalProtons=0.0

                nProton_RawEffCorrected=0.0
                nProton_RefinedEffCorrected=0.0

                #reset output file and setup histograms for loop
                #MyFile =TFile("/disk/moose/bilpa/Optima/aw/Ambiguities/nP"+(str)(nProton)+"_pitch"+(str)(pitch)+"um_Array"+(str)(size)+"cm.root","RECREATE");
                MyFile =TFile("test.root","RECREATE");
                hRawHits=TH1F("hRawHits","n Reconstructed Hits (Raw)",100,-0.5,100.5)
                hRawEfficiency=TH1F("hRawEfficiency","Raw Efficiency",20,-2.5,102.5)

                hRefinedHits=TH1F("hRefinedHits","n Reconstructed Hits (Refined)",100,-0.5,100.5)
                hRefinedEfficiency=TH1F("hRefinedEfficiency","Refined Efficiency",20,-2.5,102.5)

                myHitMap=TH2F("myHitMap","Hit Locations",(int)(2*xmax/pitch),-xmax,xmax,(int)(2*xmax/pitch),-xmax,xmax)

                for i in range(NLoops):
                        if i%1 == 0:
                                print("Processing loop "+(str)(i))

                        #nProton=nMeanProton
                        nProton=np.random.poisson(nMeanProton)
                        hNProtons.Fill(nProton)
                        totalProtons+=nProton
                        if nProton==0:
                                continue
                                
                        #generate N=10 random points for x-y 
                        XY=hf.GetRandomXY(nProton,size,width,posX,posY)
                        #print XY

                        #convert points into a strip number for each layer 
                        Strips=hf.GetStripCoOrds(XY,pitch,angles)
                        #print Strips

                        #cycle through all strip combinations, find intersections of pairs of planes, see if they are within tolerance of each other
                        #returns hit objects of form (XCoord,YCoord,[stripX,stripU,stripV],[radial distance to intersections])
                        allHits=hf.FindOverlaps(Strips,pitch,angles,tolerance,useHalfStrips)
                        nRawHits+=len(allHits)
                        hRawHits.Fill(len(allHits))
                        for hit in allHits:
                                myHitMap.Fill(hit[0],hit[1])
                                
                        rawEff=hf.GetEfficiency(allHits,XY,pitch,effTolerance)
                        hRawEfficiency.Fill(100.0*rawEff)
                        nProton_RawEffCorrected+=nProton*rawEff

                        if saveStripMaps:
                                hf.PlotHitMap("RawHits",allHits,XY,Strips,pitch,size,i,angles)


                        #for hit in allHits:
                        #        area=hf.GetPixelArea(hit,angles,pitch)
                        #        if area>1.0:
                        #                print "Area=",area,hit
                                

                        #remove duplicate hits (separated by < tolerance) 
                        refinedHits=hf.RemoveAdjacentHits(allHits,tolerance,pitch)


                        nRefinedHits+=len(refinedHits)
                        hRefinedHits.Fill(len(refinedHits))
                        refinedEff=hf.GetEfficiency(refinedHits,XY,pitch,effTolerance)
                        hRefinedEfficiency.Fill(100*refinedEff)
                        nProton_RefinedEffCorrected+=nProton*refinedEff

                        if saveStripMaps:
                                hf.PlotHitMap("RefinedHits",refinedHits,XY,Strips,pitch,size,i,angles)
                        #if len(refinedHits)<len(allHits):
                        #        print("Found a duplicate")
                        #        print allHits

                rawAmbiguity=100*(nRawHits-totalProtons)/totalProtons
                refinedAmbiguity=100*(nRefinedHits-totalProtons)/totalProtons

                refinedAmbiguityAlt=100.0*((hRefinedHits.GetMean())-nMeanProton)/nMeanProton
                fractionalUncertaintyAlt=hRefinedHits.GetMeanError()/hRefinedHits.GetMean()

                print "nMeanProton="+(str)(totalProtons)+",nReco= "+(str)(nRawHits)+",nReco2= "+(str)(nRefinedHits)+", pitch="+(str)(pitch)+", rawAmbiguity="+(str)(rawAmbiguity)+"%"+", refinedAmbiguity="+(str)(refinedAmbiguity)
                print "Efficiency: raw="+(str)(hRawEfficiency.GetMean())+", refined="+(str)(hRefinedEfficiency.GetMean())
                print "alt ambiguity=",refinedAmbiguityAlt,"+/-",fractionalUncertaintyAlt*refinedAmbiguityAlt
                hRawHits.Write()
                hRawEfficiency.Write()
                hRefinedHits.Write()
                hRefinedEfficiency.Write()
                myHitMap.Write()
                hNProtons.Write()
                MyFile.Close()
                ambiguityRate.Fill(pitch,nMeanProton,rawAmbiguity)
                ambiguityRateRefined.Fill(pitch,nMeanProton,refinedAmbiguity)
                _aMeanProton.append(nMeanProton)
                _aMeanProtonError.append(0)
                _aAmbiguity.append(refinedAmbiguity)
                _aAmbiguityError.append(refinedAmbiguity*fractionalUncertaintyAlt)
                
mycanvas = TCanvas( "c1", "Ambiguities", 550,195,800,700 )
ambiguityRate.Draw("COLZText")
mycanvas.SaveAs("OutputData/outputAmbiguities_"+(str)(width)+"cm.C")
mycanvas.SaveAs("OutputData/outputAmbiguities_"+(str)(width)+"cm.png")
ambiguityRateRefined.Draw("COLZTEXT")
mycanvas.SaveAs("OutputData/outputAmbiguitiesRefined_"+(str)(width)+"cm.C")
mycanvas.SaveAs("OutputData/outputAmbiguitiesRefined_"+(str)(width)+"cm.png")

if len(_aAmbiguity) >1:
        graph=TGraphErrors(len(_aAmbiguity),_aMeanProton,_aAmbiguity,_aMeanProtonError,_aAmbiguityError)
        graph.SetMarkerStyle(5)
        graph.GetXaxis().SetTitle("Mean N Proton")
        graph.GetYaxis().SetTitle("Ambiguity (%)")
        graph.SetTitle("")
        graph.Draw("AP")
        graph.Fit("pol2")
        mycanvas.SaveAs("OutputData/AmbiguitiesGraph.C")

