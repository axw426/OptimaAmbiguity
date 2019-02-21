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
sigma=2.0 #number of standard deviations of beamspread to try and accept

size=20 #cm

sensorThickness=155.0 #um
interPlaneDistance=12.0*1000 #um
interModuleDistance=100.0*1000 #um
phantomGap=80.0*1000 #um

geometryNames=["1ModuleXY","1ModuleXUV","2ModuleXY","2ModuleXUV","4ModuleXY","4ModuleXUV","2ModuleXUV_RT"]


def init(geoname):

        pitch=100.0

        TrackerAngles=[]
        TrackerZ=[]

        if  geoname=="1ModuleXUV":

                stripTolerance=1.0 # maximum radial distance for an equilateral triangle is 2/3 height, plus wiggle room for Z separation + cases with merged hits
                trackTolerance=stripTolerance*pitch 
                pos=0

                TrackerAngles.append([0,60,120])
                TrackerZ.append([0,0,0])

        
        elif  geoname=="1ModuleXY":

                stripTolerance=math.sqrt(2.0)/2.0
                trackTolerance=stripTolerance*pitch 
                pos=0

                TrackerAngles.append([0,90])
                TrackerZ.append([0,0])


        
        elif geoname=="2ModuleXY":

                #stripTolerance=math.sqrt(2.0)/2.0
                #stripTolerance=0.9 #needs to be slightly higher than expected due to separation between planes
                pixelsize=math.sqrt(2.0)/2.0 #distance from centre of pixel to corner

                stripTolerance=2*math.tan(sigma*beamSpread/1000.0)*interPlaneDistance/pitch + pixelsize # separation of planes+ base strip tolerance

                #add error from strip tolerance + deviation from parallel beams
                trackTolerance=2*math.tan(sigma*beamSpread/1000.0)*interModuleDistance + stripTolerance*pitch

                pos=0

                TrackerAngles.append([0,90])
                TrackerZ.append([pos,pos+interPlaneDistance])
                
                pos+=interModuleDistance+interPlaneDistance
        
                TrackerAngles.append([45,135])
                TrackerZ.append([pos,pos+interPlaneDistance])
                

        elif  geoname=="2ModuleXUV":

                
                pixelsize=0.7 #distance from centre of pixel to corner
                
                stripTolerance=2*math.tan(sigma*beamSpread/1000.0)*interPlaneDistance/pitch + pixelsize

                #add error from strip tolerance + deviation from parallel beams
                trackTolerance=2*math.tan(sigma*beamSpread/1000.0)*interModuleDistance + stripTolerance*pitch

                pos=0

                TrackerAngles.append([0,60,120])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

                pos+=interModuleDistance+2*interPlaneDistance
                
                TrackerAngles.append([30,90,150])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

        elif  geoname=="2ModuleXUV_RT":

                
                pixelsize=0.7 #distance from centre of pixel to corner
                
                stripTolerance=2*math.tan(sigma*beamSpread/1000.0)*interPlaneDistance/pitch + pixelsize

                #add error from strip tolerance + deviation from parallel beams
                trackTolerance=2*math.tan(sigma*beamSpread/1000.0)*interModuleDistance + stripTolerance*pitch

                pos=0

                TrackerAngles.append([0,60,120])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

                pos+=interModuleDistance+2*interPlaneDistance
                
                TrackerAngles.append([30,90,150])
                TrackerZ.append([pos,pos+interPlaneDistance,pos+2*interPlaneDistance])

                
        elif geoname=="4ModuleXY":

                #stripTolerance=math.sqrt(2.0)/2.0
                stripTolerance=0.9 #needs to be slightly higher than expected due to separation between planes

                totalDistance=4*interPlaneDistance+2*interModuleDistance+phantomGap
                trackTolerance=2*math.tan(beamSpread/1000.0)*totalDistance + stripTolerance*pitch

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

                stripTolerance=1.0
                
                totalDistance=4*interPlaneDistance+2*interModuleDistance+phantomGap
                trackTolerance=2*math.tan(beamSpread/1000.0)*totalDistance + stripTolerance*pitch

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

        else:
                print "\nWarning! Tried to use geometry: "+geoname+", but no such geometry exists!\n"
                print "Available geometries are:"
                for  i in geometryNames:
                        print i
                quit()

       
        ZMeans=[]
        for z in TrackerZ:
                ZMeans.append(sum(z)/len(z))
        
                
        print "\n##########################   Geometry   ##########################"
        print "Setting up geometry: "+geoname
        print "Total sensor size= ",size," cm"
        print "Pitch= ",pitch
        print "Tracker Angles= ",TrackerAngles
        print "Tracker Positions= ",TrackerZ
        print "Strip Tolerance (maximum separation of strip intersections)= ",pitch*stripTolerance
        print "Track Tolerance (maximum deviation from parallel beam to accept tracking)= ",trackTolerance
        print "###################################################################\n"
        
        return TrackerAngles,TrackerZ, ZMeans,stripTolerance,trackTolerance,pitch,beamSpread,size
