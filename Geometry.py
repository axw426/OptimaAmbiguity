import HelperFunctions as hf
import random
import math
import numpy as np
from copy import deepcopy
from ROOT import gROOT, TCanvas, TF1, TFile, gStyle,TH2F, TH1F
import time

hf.SetSeed(2022)

#stripTolerance=1.5 #4 planes per module...for some reason- shape varies depending on position

beamSpread=3.0 #mrad

sensorThickness=155.0 #um
interPlaneDistance=12.0*1000 #um
interModuleDistance=100.0*1000 #um
phantomGap=80.0*1000 #um

def init(geoname):

        pitch=100.0
        effTolerance=2500 #needs to be thought about...

        TrackerAngles=[]
        TrackerZ=[]

        if  geoname=="1ModuleXUV":

                stripTolerance=0.67 # maximum radial distance for an equilateral triangle is 2/3 height
                effTolerance=stripTolerance
                trackTolerance=stripTolerance*pitch 
                pos=0

                TrackerAngles.append([0,60,120])
                TrackerZ.append([0,0,0])

        
        elif  geoname=="1ModuleXY":

                stripTolerance=0.0
                effTolerance=math.sqrt(2.0)/2.0
                trackTolerance=stripTolerance*pitch 
                pos=0

                TrackerAngles.append([0,90])
                TrackerZ.append([0,0])


        
        elif geoname=="2ModuleXY":

                stripTolerance=0.0
                
                #trackTolerance=2*math.tan(beamSpread/1000.0)*interModuleDistance
                trackTolerance=10
                
                pos=0

                TrackerAngles.append([0,90])
                TrackerZ.append([pos,pos+interPlaneDistance])
                #TrackerZ.append([0,0])
                
                pos+=interModuleDistance+interPlaneDistance
        
                TrackerAngles.append([45,135])
                TrackerZ.append([pos,pos+interPlaneDistance])
                #TrackerZ.append([0,0])
                

        elif  geoname=="2ModuleXUV":

                stripTolerance=0.67

                trackTolerance=stripTolerance*pitch 
                pos=0

                TrackerAngles.append([0,60,120])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

                pos+=interModuleDistance+2*interPlaneDistance
                
                TrackerAngles.append([30,90,150])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

        elif geoname=="4ModuleXY":

                stripTolerance=0.0

                totalDistance=4*interPlaneDistance+2*interModuleDistance+phantomGap

                trackTolerance=2*math.tan(beamSpread/1000.0)*totalDistance 
                pos=0

                TrackerAngles.append([0,90])
                TrackerZ.append([pos,pos+interPlaneDistance])
                
                pos+=interModuleDistance+interPlaneDistance
        
                TrackerAngles.append([45,135])
                TrackerZ.append([pos,pos+interPlaneDistance])
                
                pos+=phantomGap+interPlaneDistance
        
                TrackerAngles.append([0,90])
                TrackerZ.append([pos,pos+interPlaneDistance])

                pos+=interModuleDistance+interPlaneDistance

                TrackerAngles.append([45,135])
                TrackerZ.append([pos,pos+interPlaneDistance])

        elif  geoname=="4ModuleXUV":

                stripTolerance=0.67
                
                totalDistance=4*interPlaneDistance+2*interModuleDistance+phantomGap

                trackTolerance=2*math.tan(beamSpread/1000.0)*totalDistance +500 
                pos=0

                TrackerAngles.append([0,60,120])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

                pos+=interModuleDistance+2*interPlaneDistance
                
                TrackerAngles.append([30,90,150])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

                pos+=phantomGap+2*interPlaneDistance
                
                TrackerAngles.append([0,60,120])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

                pos+=interModuleDistance+2*interPlaneDistance

                TrackerAngles.append([30,90,150])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])        


        return TrackerAngles,TrackerZ,stripTolerance,trackTolerance,effTolerance,pitch,beamSpread
