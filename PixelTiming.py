import HelperFunctions as hf
import random
import math
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F
import time

hf.SetSeed(2022)

NLoops=10
size=10 #cm
saveStripMaps=False

nProton=13
pitch=50
xmax=size*5000

gStyle.SetOptStat(0000)
gROOT.SetBatch(True)

histEntries=[]
hitArray=[]
counter=0
supercounter=0

MyFile =TFile("pixeltiming.root","RECREATE");

myHist=TH1F("myHist","Time between pixel firings",75,-0.5,149.6)

while True:
        if len(histEntries)>5000:
                break

                
        #generate N=10 random points for x-y 
        XY=hf.GetRandomXY(nProton,size)
        
        #convert points into a strip number for each layer 
        coords=hf.GetPixelCoOrds(XY,pitch)

        for coord in coords:
                if coord in  hitArray:
                        histEntries.append(counter)
                        print "Found duplicate after "+(str)(counter)+" loops, total found so far is "+(str)(len(histEntries))
                        counter=0
                        hitArray=[]
                        break
                else:
                       hitArray.append(coord) 
        
        counter+=1
        supercounter+=1
print supercounter

for time in histEntries:
        myHist.Fill(time)

myHist.SetTitle("n Timestamps before pixel is fired again")
myHist.GetXaxis().SetTitle("N Timestamps")
myHist.Write()
