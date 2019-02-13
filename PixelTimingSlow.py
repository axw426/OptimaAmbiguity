import HelperFunctions as hf
import random
import math
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F
import time

hf.SetSeed(2022)

NLoops=10
size=1.2 #cm
saveStripMaps=False

nProton=13
pitch=50
xmax=size*5000

#gStyle.SetOptStat(0000)
gROOT.SetBatch(True)

histEntries=[]
hitArray=[]
counterArray=[]
counter=0

MyFile =TFile("pixeltimingSlow.root","RECREATE");

myHist=TH1F("myHist","Time between pixel firings",200,-0.0,10000)

while True:
        if len(histEntries)>5000:
                break

                
        #generate N=10 random points for x-y 
        XY=hf.GetRandomXY(nProton,size)
        
        #convert points into a strip number for each layer 
        coords=hf.GetPixelCoOrds(XY,pitch)

        for coord in coords:
                if coord in hitArray:
                        i=hitArray.index(coord)
                        histEntries.append(counter-counterArray[i])
                        print "Found duplicate after "+(str)(counter-counterArray[i])+" loops, total found so far is "+(str)(len(histEntries))
                        counterArray[i]=counter
                        #counter=0
                        #break
                else:
                       hitArray.append(coord) 
                       counterArray.append(counter) 
        
        counter+=1
        #print counter

for time in histEntries:
        myHist.Fill(time)

myHist.SetTitle("n Timestamps before pixel is fired again")
myHist.GetXaxis().SetTitle("N Timestamps")
myHist.Write()
