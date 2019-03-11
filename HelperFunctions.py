import random
import math
from ROOT import gROOT, TCanvas, TF1, TF2, TH2F, TFile, TLine, TH1F, TLegend,TGraph,TGraphErrors,TObject
from array import array
import copy
import numpy as np

def SetSeed(seed):
	random.seed(seed)

def GetRandomXY(nProtons,size,width=6,posX=0,posY=0):
	XY=[]
	for x in range(nProtons):
		#X=random.uniform(-size*5000.0,size*5000.0) #10000 um per cm
		#Y=random.uniform(-size*5000.0,size*5000.0) #10000 um per cm
		X=random.gauss(posX*10000.0,width*1000.0)
		Y=random.gauss(posY*10000.0,width*1000.0)
                XY.append([X,Y])
	        #print X,Y
	return XY

def GetDirections(nProtons,beamSpread,cheat):
        mXmY=[]
        # m = tan(theta)
        for x in range(nProtons):
                #mX=math.tan(random.uniform(-beamSpread/1000.0,beamSpread/1000.0))
                #mY=math.tan(random.uniform(-beamSpread/1000.0,beamSpread/1000.0))
                mX=math.tan(random.gauss(0.0,beamSpread/1000.0))
                mY=math.tan(random.gauss(0.0,beamSpread/1000.0))
                #mX=math.tan(beamSpread/1000.0)
                #mY=math.tan(beamSpread/1000.0)
                if cheat==False:
                        mXmY.append([mX,mY])
                else:
                        mXmY.append([0,0])
        return mXmY

def ConvertXYToStrip(coord,pitch,rawAngle):

        rotatedX= coord[0]*math.cos(math.radians(rawAngle)) - coord[1]*math.sin(math.radians(rawAngle))
        if rotatedX >=0.5:
                Strip=int((rotatedX+pitch/2.0)/pitch)
        elif rotatedX <=-0.5:
                Strip=int((rotatedX-pitch/2.0)/pitch)
        else:
                Strip=0

        return Strip
def CheckStripHalf(coord,pitch,rawAngle):

        #assuming centre of sensor is at y=0.0,x=0.0
        y= (int)((coord[0]*math.sin(math.radians(rawAngle)) + coord[1]*math.cos(math.radians(rawAngle))) /pitch)

        return y

def GetTrackerStripCoOrds(XY,mXmY,pitch,angles,Z):
	StripX=[]
	StripU=[]
	StripV=[]
	StripY=[]
        StripHalfX=[]
        StripHalfU=[]
        StripHalfV=[]
        StripHalfY=[]
        
        for coord,direction in zip(XY,mXmY):

                #updated pos=X0+mX*Z
                coordX=[coord[0]+direction[0]*Z[0],coord[1]+direction[1]*Z[0]]
                StripX.append(ConvertXYToStrip(coordX,pitch,angles[0]))
                StripHalfX.append(CheckStripHalf(coordX,pitch,angles[0]))

                coordU=[coord[0]+direction[0]*Z[1],coord[1]+direction[1]*Z[1]]
		StripU.append(ConvertXYToStrip(coordU,pitch,angles[1]))
                StripHalfU.append(CheckStripHalf(coordU,pitch,angles[1]))

        Strips=[StripX,StripU,StripHalfX,StripHalfU]

        if len(angles)>2:
                for coord,direction in zip(XY,mXmY):
                        coordV=[coord[0]+direction[0]*Z[2],coord[1]+direction[1]*Z[2]]
                        StripV.append(ConvertXYToStrip(coordV,pitch,angles[2]))
                        StripHalfV.append(CheckStripHalf(coordV,pitch,angles[2]))

                Strips=[StripX,StripU,StripV,StripHalfX,StripHalfU,StripHalfV]
                        
        if len(angles)>3:
                for coord,direction in zip(XY,mXmY):
                        coordY=[coord[0]+direction[0]*Z[3],coord[1]+direction[1]*Z[3]]
                        StripY.append(ConvertXYToStrip(coordY,pitch,angles[3]))
                        StripHalfY.append(CheckStripHalf(coordY,pitch,angles[3]))
	        Strips= [StripX,StripU,StripV,StripY,StripHalfX,StripHalfU,StripHalfV,StripHalfY]	
 
        meanZ=sum(Z)/len(Z)
        meanXY=[]
        for coord,direction in zip(XY,mXmY):
                meanXY.append([coord[0]+direction[0]*meanZ,coord[1]+direction[1]*meanZ])
                        
                
        return Strips,meanXY

def GetStripCoOrds(XY,pitch,angles):
	StripX=[]
	StripU=[]
	StripV=[]
	StripY=[]
        StripHalfX=[]
        StripHalfU=[]
        StripHalfV=[]
        StripHalfY=[]
        for coord in XY:
                StripX.append(ConvertXYToStrip(coord,pitch,angles[0]))
		StripU.append(ConvertXYToStrip(coord,pitch,angles[1]))
                StripHalfX.append(CheckStripHalf(coord,pitch,angles[0]))
                StripHalfU.append(CheckStripHalf(coord,pitch,angles[1]))

        Strips=[StripX,StripU,StripHalfX,StripHalfU]

        if len(angles)>2:
                for coord in XY:
                        StripV.append(ConvertXYToStrip(coord,pitch,angles[2]))
                        StripHalfV.append(CheckStripHalf(coord,pitch,angles[2]))
                Strips=[StripX,StripU,StripV,StripHalfX,StripHalfU,StripHalfV]
                        
        if len(angles)>3:
                for coord in XY:
                        StripY.append(ConvertXYToStrip(coord,pitch,angles[3]))
                        StripHalfY.append(CheckStripHalf(coord,pitch,angles[3]))
	        Strips= [StripX,StripU,StripV,StripY,StripHalfX,StripHalfU,StripHalfV,StripHalfY]	

        return Strips

def GetPixelCoOrds(XY,pitch):
        coords=[]
	for coord in XY:
		cellX=(int)((coord[0]-(pitch/2.0))/pitch)
		cellY=(int)((coord[1]-(pitch/2.0))/pitch)
                coords.append([cellX,cellY])
                
	return coords


def FindMandC(U,theta,pitch):
     	xPrime=(U*pitch)
	
	m=-1*math.tan(math.radians(90-theta))
	c=xPrime/math.cos(math.radians(90-theta))

        return m,c

def FindIntersect(X,U,theta,pitch):

	#in X frame, x*100 to get to middle of strip
	#U line described by y=mx+c, where: m=-tan(90-Theta),c= X'/cos(90-theta)

	x=(X*pitch)#+pitch/2.0
	m,c=FindMandC(U,-theta,pitch)
	y=(m*x)+c

	return [x,y]
	
def FindUVIntersect(U,V,uAngle,vAngle,pitch):

        if uAngle==0:
                return FindIntersect(U,V,vAngle,pitch)

        if vAngle==0:
                return FindIntersect(V,U,uAngle,pitch)

        #define y=mx+c for U and V in frame of X layer

        m1,c1=FindMandC(U,-uAngle,pitch)
        m2,c2=FindMandC(V,-vAngle,pitch)

	Xintercept=(c2-c1)/(m1-m2)
	Yintercept1=(m1*Xintercept) + c1
	Yintercept2=(m2*Xintercept) + c2

	return [Xintercept,Yintercept1]
	

def CheckProximity(xycoords,tolerance,pitch):	

        passed=True

        xmean=0.0
	ymean=0.0

        for i in range(len(xycoords)):
                xmean+=xycoords[i][0]
                ymean+=xycoords[i][1]
        xmean/=len(xycoords)        
        ymean/=len(xycoords)        

        rValues=[]
        for i in range(len(xycoords)):
                r= math.sqrt( (xycoords[i][1]-ymean)**2  + (xycoords[i][0]-xmean)**2)
                rValues.append(r)
                if r>=tolerance*pitch:
                        passed=False

        return passed,xmean,ymean,rValues

def RestrictLine(l, xmin, ymin, xmax, ymax,m,c):

        #tidies up lines so they don't cross axis
        y_at_xmin = m*xmin+c
        y_at_xmax = m*xmax+c
        if abs(m) >0:
                x_at_ymin = (ymin-c)/m
                x_at_ymax = (ymax-c)/m

                xcoords=[]
                ycoords=[]
        
                if y_at_xmin <= ymax and y_at_xmin >= ymin :
                        #crosses left axis
                        xcoords.append(xmin)
                        ycoords.append(y_at_xmin)
                if y_at_xmax <= ymax and y_at_xmax >= ymin :
                        #crosses right axis
                        xcoords.append(xmax)
                        ycoords.append(y_at_xmax)
                if x_at_ymin <= xmax and x_at_ymin >= xmin :
                        #crosses bottom
                        xcoords.append(x_at_ymin)
                        ycoords.append(ymin)
                if x_at_ymax <= xmax and x_at_ymax >= xmin :
                        #crosses top
                        xcoords.append(x_at_ymax)
                        ycoords.append(ymax)
                if xcoords:        
                        l.SetX1(xcoords[0])
                        l.SetX2(xcoords[1])
                        l.SetY1(ycoords[0])
                        l.SetY2(ycoords[1])

        return l

def CreateTLines(xmax,U,uAngle,pitch):

        if uAngle==0:
                myline=TLine(U*pitch,-xmax,U*pitch,xmax)
                myline.SetLineColor(1)
                myline.SetLineWidth(1)
                return myline
        else:
	        m1,c1=FindMandC(U,-uAngle,pitch)
	        myline= TLine(-xmax,m1*(-xmax)+c1,xmax,m1*(xmax)+c1 )
                myline.SetLineColor(1)
                myline.SetLineWidth(1)

                return RestrictLine(myline, -xmax, -xmax, xmax, xmax,m1,c1)
        

def checkForDuplicates(xmean,ymean,allHits,tolerance):
        
        for hit in allHits:
                r=math.sqrt((xmean-hit[0])**2 + (ymean-hit[1])**2)
                if r<tolerance*pitch:
                        return True
        return False

def CheckHalves(coord,pitch,halves,angles):

        passed=True

        for half,angle in zip(halves,angles):
                #check that y coordinates in planes frame have same sign
                if half*CheckStripHalf(coord,pitch,angle)<0:
                     passed=False

        return passed
             


def FindOverlaps2Planes(Strips,pitch,angles,tolerance,useHalfStrips):

        allHits=[]
        
	#find the Y coordinate where the X and U or V intercept and see if they are close
	for X,xHalf in zip(Strips[0],Strips[2]):
		for U,uHalf in zip(Strips[1],Strips[3]):

			xycoords=FindUVIntersect(X,U,angles[0],angles[1],pitch)
                        xmean=xycoords[0]
                        ymean=xycoords[1]
                        
                        if useHalfStrips and CheckHalves([xmean,ymean],pitch,[xHalf,uHalf],angles)==False:
                                continue
                        
                        allHits.append([xmean,ymean,[X,U],[0]])

	return allHits


def FindOverlaps3Planes(Strips,pitch,angles,tolerance,useHalfStrips):

        allHits=[]
        
	#find the Y coordinate where the X and U or V intercept and see if they are close
	for X,xHalf in zip(Strips[0],Strips[3]):
		for U,uHalf in zip(Strips[1],Strips[4]):
			   for V,vHalf in zip(Strips[2],Strips[5]):
                                xycoords=[]
				xycoords.append(FindUVIntersect(X,U,angles[0],angles[1],pitch))
				xycoords.append(FindUVIntersect(X,V,angles[0],angles[2],pitch))
				xycoords.append(FindUVIntersect(U,V,angles[1],angles[2],pitch))

				passed,xmean,ymean,rValues=CheckProximity(xycoords,tolerance,pitch)
                                
                                if passed==True:
                                        if useHalfStrips and CheckHalves([xmean,ymean],pitch,[xHalf,uHalf,vHalf],angles)==False:
                                                continue
                                        allHits.append([xmean,ymean,[X,U,V],rValues])

	return allHits

def FindOverlaps4Planes(Strips,pitch,angles,tolerance,useHalfStrips):

        allHits=[]
	#find the Y coordinate where the X and U or V intercept and see if they are close
	for X,xHalf in zip(Strips[0],Strips[4]):
		for U,uHalf in zip(Strips[1],Strips[5]):
			for V,vHalf in zip(Strips[2],Strips[6]):
			        for Y,yHalf in zip(Strips[3],Strips[7]):
                                        xycoords=[]

                                        xycoords.append(FindUVIntersect(X,U,angles[0],angles[1],pitch))
				        xycoords.append(FindUVIntersect(X,V,angles[0],angles[2],pitch))
				        xycoords.append(FindUVIntersect(X,Y,angles[0],angles[3],pitch))

                                        xycoords.append(FindUVIntersect(U,V,angles[1],angles[2],pitch))
                                        xycoords.append(FindUVIntersect(U,Y,angles[1],angles[3],pitch))

                                        xycoords.append(FindUVIntersect(V,Y,angles[2],angles[3],pitch))

				        passed,xmean,ymean,rValues=CheckProximity(xycoords,tolerance,pitch)
                                        if passed==True:
                                                if useHalfStrips and CheckHalves([xmean,ymean],pitch,[xHalf,uHalf,vHalf,yHalf],angles)==False:
                                                        continue
                                                allHits.append([xmean,ymean,[X,U,V,Y],rValues])

	return allHits

#cycle through all strip combinations and see which create an overlap- will want to define a central line for each strip and work out if they intersect 
def FindOverlaps(Strips,pitch,angles,tolerance,useHalfStrips):

        nPlanes=len(angles)
        if nPlanes==2:
                return FindOverlaps2Planes(Strips,pitch,angles,tolerance,useHalfStrips) 
        elif nPlanes==3:
                return FindOverlaps3Planes(Strips,pitch,angles,tolerance,useHalfStrips)
        elif nPlanes==4:
                return FindOverlaps4Planes(Strips,pitch,angles,tolerance,useHalfStrips)  


def ReadXUVStripCoOrds(inFile):
        hitsByTimestamp=[]
        lastTimestamp=-1
        xcoords=[]
        ucoords=[]
        vcoords=[]

        hXStrips=TH1F("hXStrips","n Strip Positions X",1024,-0.5,1023.5)
        hUStrips=TH1F("hUStrips","n Strip Positions U",1024,-0.5,1023.5)
        hVStrips=TH1F("hVStrips","n Strip Positions V",1024,-0.5,1023.5)


        with open(inFile,'r') as f:
                line=f.readline()#first line is empty...
                line=f.readline()#second line isn't real data...
                line=f.readline()#third line isn't real data...
                line=f.readline()
                while line:
                        values=line.split()
                        timestamp=values[0]
                        #if timestamp has changed and there was a hit recorded, save data
                        if timestamp!=lastTimestamp and xcoords:

                                xcoords[:]=[-(var-512) for var in xcoords]
                                ucoords[:]=[-(var-512) for var in ucoords]
                                vcoords[:]=[-(var-512) for var in vcoords]
                                hitsByTimestamp.append([xcoords,ucoords,vcoords])

                                lastTimestamp=timestamp
                                ucoords=[]
                                xcoords=[]
                                vcoords=[]

                                
                        #loop over 4 data entries and only add them if there's a hit in all planes- format of data is V,X,U (i.e. middle value is the 0 degrees plane)
                        for i in range (1,5):

                                if values[i] !='-1' and values[i+4] !='-1' and values[i+8] !='-1':
                                        hXStrips.Fill((int)(values[i+4]))
                                        hVStrips.Fill((int)(values[i]))
                                        hUStrips.Fill((int)(values[i+8]))

                                #check value is valid
                                if ( values[i]!='-1' and values[i+4]!='-1' and values[i+8]!='-1'):

                                        #check isn't a duplicate value (two adjacent strips...), lazy for now, just take first of two adjacent
                                        #if(xcoords and ((int)(values[i])-vcoords[-1]<=2 or (int)(values[i+4])-xcoords[-1]<=2 or int(values[i+8])-ucoords[-1]<=2)):
                                        #        continue

                                        # print (values[i],values[i+4],values[i+8])
                                        vcoords.append((int)(values[i]))
                                        xcoords.append((int)(values[i+4]))
                                        ucoords.append((int)(values[i+8]))

                        line=f.readline()

        hXStrips.Write()
        hUStrips.Write()
        hVStrips.Write()
        return hitsByTimestamp

def PlotHitMap(name,allHits,XY,Strips,pitch,size,loop,angles):

        canvas1 = TCanvas( 'c1'+(str)(loop), "mycanvas", 200, 10, 700, 500 )
        xmax=size*5000

        #draw reconstructed points
        rawmeasured=TH2F(name+(str)(loop),"Hit Locations",3000,-xmax,xmax,3000,-xmax,xmax)
        rawmeasured.GetXaxis().SetTitle("X position (#mum)")
        rawmeasured.GetYaxis().SetTitle("Y position (#mum)")
        rawmeasured.SetMarkerColor(4)
        rawmeasured.SetMarkerSize(2)
        rawmeasured.SetMarkerStyle(8)
        for hit in allHits:
               rawmeasured.Fill(hit[0],hit[1])
        rawmeasured.Draw("")

        linearray=[]
        #draw the strips which were hit
        nPlanes=len(angles)

       
        if nPlanes==3:
                for X,U,V in zip(Strips[0],Strips[1],Strips[2]):

                        linearray.append(CreateTLines(xmax,X,angles[0],pitch))
                        linearray.append(CreateTLines(xmax,U,angles[1],pitch))
                        linearray.append(CreateTLines(xmax,V,angles[2],pitch))

                        if len(allHits)==1:

                                linearray.append(CreateTLines(xmax,X-1,angles[0],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,U-1,angles[1],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,V-1,angles[2],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,X+1,angles[0],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,U+1,angles[1],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,V+1,angles[2],pitch))
                                linearray[-1].SetLineColor(3)

        if nPlanes==4:
                for X,U,V,Y in zip(Strips[0],Strips[1],Strips[2],Strips[3]):

                        linearray.append(CreateTLines(xmax,X,angles[0],pitch))
                        linearray.append(CreateTLines(xmax,U,angles[1],pitch))
                        linearray.append(CreateTLines(xmax,V,angles[2],pitch))
                        linearray.append(CreateTLines(xmax,Y,angles[3],pitch))

                        if len(allHits)==1:
                                linearray.append(CreateTLines(xmax,X-1,angles[0],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,U-1,angles[1],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,V-1,angles[2],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,Y-1,angles[3],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,X+1,angles[0],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,U+1,angles[1],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,V+1,angles[2],pitch))
                                linearray[-1].SetLineColor(3)
                                linearray.append(CreateTLines(xmax,Y+1,angles[3],pitch))
                                linearray[-1].SetLineColor(3)
        for line in linearray:
                line.Draw("SAME")
                               

        #draw where protons really hit
        truth=TH2F("truth"+(str)(loop),"Hit Locations",3000,-xmax,xmax,3000,-xmax,xmax)
        truth.SetMarkerSize(2)
        truth.SetMarkerStyle(23)
        truth.SetMarkerColor(2)
        for coords in XY:
                truth.Fill(coords[0],coords[1])
        truth.Draw("SAME")


        canvas1.Write(name+(str)(loop))

def GetSeparation(hit1,hit2):

        return math.sqrt((hit1[0]-hit2[0])**2 + (hit1[1]-hit2[1])**2)

    

def RemoveAdjacentHits(allHits,tolerance,pitch):

        #thou must not edit python lists while iterating over them...
        refinedHits=[]

        for i in allHits:
                passed=True
                #check if hit is close to any existing accepted hits
                for j in refinedHits:
                        if GetSeparation(i,j)<=tolerance*pitch:
                                passed=False
                                break
                               
                if passed==True:
                        refinedHits.append(i)

        return refinedHits


def SumHits(hits):

        #only bothering to correct x-y position, not sure it makes as much sense for other hit parameters
        xmean=0.0
        ymean=0.0

        for hit in hits:
                xmean+=hit[0]
                ymean+=hit[1]
        xmean/=(float)(len(hits))
        ymean/=(float)(len(hits))

        newHit=hits[0]
        newHit[0]=xmean
        newHit[1]=ymean
        newHit.append("modified")

        return newHit

def MergeAdjacentHits(allHits,tolerance,pitch):

        #thou must not edit python lists while iterating over them...
        refinedHits=[]
        remainingHits=copy.deepcopy(allHits)

        while True:
                if len(remainingHits)==0:
                        break

                i=remainingHits[0]
                toBeSummed=[i]
                remainingHits.remove(i)

                currentHits=copy.deepcopy(remainingHits)
                for j in currentHits:
                        if GetSeparation(i,j)<=tolerance*pitch:
                                toBeSummed.append(j)
                                remainingHits.remove(j)
                refinedHits.append(SumHits(toBeSummed))
                                
                        
        return refinedHits

def GetEfficiency(allHits,XY,pitch,tolerance):

        #check if there is a reconstructed hit within tolerance of the true hit positions
        nExpected=len(XY)
        nFound=0.0

        for true in XY:
                passed=False
                for reco in allHits:
                        if math.sqrt((true[0]-reco[0])**2 + (true[1]-reco[1])**2)<pitch*tolerance:
                                #print math.sqrt((true[0]-reco[0])**2 + (true[1]-reco[1])**2)
                                passed=True
                                break
                if passed ==True:
                        nFound+=1.0


        return nFound/(float)(nExpected)


def GetPixelArea(hit,angles,pitch):

        #for perfect overlap in 3 planes --> area =pitch*pitch/tan(60)
        
        #find crossing points for XU,XV,and UV and assume triangular
        rSum=0.0

        Strips=hit[2]
	XU=FindUVIntersect(Strips[0],Strips[1],angles[0],angles[1],pitch)
	XV=FindUVIntersect(Strips[0],Strips[2],angles[0],angles[2],pitch)
	UV=FindUVIntersect(Strips[1],Strips[2],angles[1],angles[2],pitch)

        area=abs(XU[0]*(XV[1]-UV[1]) + XV[0]*(UV[1]-XU[1]) +UV[0]*(XU[1]-XV[1]))/2.0

        return area
        
def CheckPixelAreas(allHits,areaTolerance,pitch,angles):

        acceptedHits=[]
        for hit in allHits:
                if GetPixelArea(hit,angles,pitch)<areaTolerance:
                        acceptedHits.append(hit)
        return acceptedHits


def RemoveAmbiguities(inHits,rawAngle,pitch):

        XAcceptedHits=[]
        XUAcceptedHits=[]
        XUVAcceptedHits=[]

        #check unique X
        for hitA in inHits:
                passed=True
                for hitB in inHits:
                        if hitA[2][0]==hitB[2][0] and GetPixelArea(hitA,rawAngle,pitch)>GetPixelArea(hitB,rawAngle,pitch):
                                passed=False
                                break
                if passed==True:
                        XAcceptedHits.append(hitA)

        #check unique U
        for hitA in XAcceptedHits:
                passed=True
                for hitB in XAcceptedHits:
                        if hitA[2][1]==hitB[2][1] and GetPixelArea(hitA,rawAngle,pitch)>GetPixelArea(hitB,rawAngle,pitch):
                                passed=False
                                break
                if passed==True:
                       XUAcceptedHits.append(hitA)           

        #check uniqueV
        for hitA in XUAcceptedHits:
                passed=True
                for hitB in XUAcceptedHits:
                        if hitA[2][2]==hitB[2][2] and GetPixelArea(hitA,rawAngle,pitch)>GetPixelArea(hitB,rawAngle,pitch):
                                passed=False
                                break
                if passed==True:
                       XUVAcceptedHits.append(hitA)  

        return XUVAcceptedHits

def GetGradient(hitI,hitJ,TrackerZ):
        
        deltaZ=np.mean(TrackerZ[1])-np.mean(TrackerZ[0])

        mX=(hitJ[0]-hitI[0])/deltaZ
        mY=(hitJ[1]-hitI[1])/deltaZ

        cX=hitI[0]-mX*np.mean(TrackerZ[0])
        cY=hitI[1]-mY*np.mean(TrackerZ[0])

        return [cX,mX,cY,mY]
        
def ReconstructTracks2Planes(Hits,tolerance,pitch,TrackerZ):

        RecoTracks=[]
        #loop over all hits pairs and find pair with minimum separation, add as track and remove hits from hit database
        while len(Hits[0])>0 and len(Hits[1])>0:
                minimumSeparation=100000
                bestIPos=-1
                bestJPos=-1
                
                for i in range(len(Hits[0])):
                        for j in range(len(Hits[1])):
                                hitI=Hits[0][i]
                                hitJ=Hits[1][j]
                                r=GetSeparation(hitI,hitJ)
                                if r<minimumSeparation:
                                        bestIPos=i
                                        bestJPos=j
                                        minimumSeparation=r
                                        bestFitParams=GetGradient(hitI,hitJ,TrackerZ)
                if minimumSeparation<tolerance:
                        #Append hit positions to track list
                        RecoTracks.append([Hits[0][bestIPos],Hits[1][bestJPos],bestFitParams])
                        #delete entries from Hits
                        del Hits[0][bestIPos]
                        del Hits[1][bestJPos]
                else:
                        break

        return RecoTracks

def GetChi2(hits,TrackerZ,stripTolerance,pitch):

        
        _X,_Y,_Z=array('d'),array('d'),array('d')
        _XErr,_YErr,_ZErr=array('d'),array('d'),array('d')

        for hit,Z in zip(hits,TrackerZ):

                _X.append(hit[0])
                _XErr.append(stripTolerance*pitch)
                _Y.append(hit[1])
                _YErr.append(stripTolerance*pitch)
                _Z.append(np.mean(Z))
                _ZErr.append(0.0)

        myfit=TF1("myfit","[0]+[1]*x")

        graphX=TGraphErrors(len(_Z),_Z,_X,_ZErr,_XErr)
        graphX.Fit("myfit","Q")
        myfitX=graphX.GetFunction("myfit")
        chi2X=myfitX.GetChisquare()
        cX=myfitX.GetParameter(0)
        mX=myfitX.GetParameter(1)
        
        graphY=TGraphErrors(len(_Z),_Z,_Y,_ZErr,_YErr)
        graphY.Fit("myfit","Q")
        myfitY=graphY.GetFunction("myfit")
        chi2Y=myfitY.GetChisquare()
        cY=myfitY.GetParameter(0)
        mY=myfitY.GetParameter(1)

        #graphX.Write()
        #graphY.Write()
        del graphX,graphY

        #print chi2X,chi2Y

        return chi2X*chi2Y, [cX,mX,cY,mY]


def GetChi2_SingleFit(hits,TrackerZ,stripTolerance,pitch):

        #does a chi2 to check trajectory across 3 hits is linear in both x and y, should probably be a simultaneous fit
        #causes a memory leak in ROOT6??? Need to run in sl6 for now...

        
        _X,_Y=array('d'),array('d')
        _XErr,_YErr=array('d'),array('d')

        for hit in hits:

                _X.append(hit[0])
                _XErr.append(2*pitch) # 2*pitch seems to give sensible chi2--> really will depend on angle of track
                _Y.append(hit[1])
                _YErr.append(2*pitch)

        #X and Y have shared Z and both are linear with Z so will be linear relative to each other
        myfit=TF1("myfit","[0]+[1]*x")

        #find reasonable starting parameters
        startM=(_Y[1]-_Y[0])/(_X[1]-_X[0])
        startC=_Y[0]-startM*_X[0]
        myfit.SetParameters(startC,startM)
        
        graph=TGraphErrors(len(_Y),_X,_Y,_XErr,_YErr)
        graph.Fit("myfit","Q")
        myfitResult=graph.GetFunction("myfit")
        chi2=myfitResult.GetChisquare()
        
        del graph

        return chi2

def CheckProjection(i,j,TrackerZ,stripTolerance,pitch, beamSpread):

        ZSeparation=np.mean(TrackerZ[1])-np.mean(TrackerZ[0])
        sigma=2
        maxSeparation=2*math.tan(sigma*beamSpread/1000.0)*ZSeparation + stripTolerance*pitch
        separation=GetSeparation(i,j)
        if separation<maxSeparation:
                return True
        else:
                return False

def ReconstructTracks3Planes(Hits,tolerance,pitch,TrackerZ,stripTolerance,beamSpread):

        RecoTracks=[]
        while len(Hits[0])>0 and len(Hits[1])>0 and len(Hits[2])>0:
                bestChi2=100000
                bestIPos=-1
                bestJPos=-1
                bestKPos=-1
                bestFitParam=[]
                for i in range(len(Hits[0])):
                        hitI=Hits[0][i]
                        for j in range(len(Hits[1])):
                                hitJ=Hits[1][j]

                                #speed optimization- check if both XUV points are compatable for expected beam spread
                                if CheckProjection(hitI,hitJ,TrackerZ,stripTolerance,pitch,beamSpread)==False:
                                        continue
                                        
                                for k in range(len(Hits[2])):
                                        hitK=Hits[2][k]

                                        if CheckProjection(hitJ,hitK,TrackerZ,stripTolerance,pitch,beamSpread)==False:
                                                continue
                                
                                        chi2,fitParam=GetChi2([hitI,hitJ,hitK],TrackerZ,stripTolerance,pitch)
                                        if chi2<bestChi2:
                                                bestIPos=i
                                                bestJPos=j
                                                bestKPos=k
                                                bestChi2=chi2
                                                bestFitParam=fitParam
                #print bestChi2
                if bestChi2<tolerance:
                        #Append hit positions to track list
                        RecoTracks.append([Hits[0][bestIPos],Hits[1][bestJPos],Hits[2][bestKPos],bestFitParam])
                        #delete entries from Hits
                        del Hits[0][bestIPos]
                        del Hits[1][bestJPos]
                        del Hits[2][bestKPos]
                else:
                        break
        return RecoTracks

def ReconstructTracks3PlanesAlt(Hits,tolerance,pitch,TrackerZ,stripTolerance,beamSpread):

        #alternative version in which hits aren't removed after being accepted into track, instead accept all allowed combinations
        #gives ~100% efficiency, but drastically sacrifices ambiguity- need to refine tracks afterwards
        
        RecoTracks=[]

        for i in range(len(Hits[0])):
                bestIPos=-1
                bestJPos=-1
                bestKPos=-1
                hitI=Hits[0][i]
                for j in range(len(Hits[1])):
                        hitJ=Hits[1][j]
                        
                        #speed optimization- check if both XUV points are compatable for expected beam spread
                        if CheckProjection(hitI,hitJ,TrackerZ,stripTolerance,pitch,beamSpread)==False:
                                continue
                                        
                        for k in range(len(Hits[2])):
                                hitK=Hits[2][k]
                                
                                if CheckProjection(hitJ,hitK,TrackerZ,stripTolerance,pitch,beamSpread)==False:
                                        continue
                                
                                chi2=GetChi2([hitI,hitJ,hitK],TrackerZ,stripTolerance,pitch)

                                if chi2<tolerance:
                                        #Append hit positions to track list
                                        RecoTracks.append([hitI,hitJ,hitK,chi2])
                
        return RecoTracks


def ReconstructTracksAlt(Hits,tolerance,pitch,TrackerZ,stripTolerance,beamSpread):

        #loosely based on more traditional approaches of segment finding, though not much can be done with 3 points....

        #procedure:
        #1) Start with most precise points (first XUV) as seeds
        #2) Use uncertainty in position + expected MS angles to match to hits in second module i.e. check difference in radius ~ module separation*tanTheta + 2*stripTolerance 
        #3) Use vector from first two hits to project to RT and match hits in acceptable area 
        #4) Filter tracks somehow- merge similar tracks that are nearby or that use the same hits

        tracks=[]
        for i in Hits[0]:
                segments=[]
                for j in Hits[1]:
                        if CheckProjection(i,j,TrackerZ,stripTolerance,pitch,beamSpread):
                                segments.append([i,j])
                                #tracks.append([i,j])
                for segment in segments:
                        moduleDist=np.mean(TrackerZ[1])-np.mean(TrackerZ[0])
                        RTDist=np.mean(TrackerZ[2])-np.mean(TrackerZ[1])
                        
                        #find y-gradient between module and extrapolate to RT
                        yM=(segment[1][1]-segment[0][1])/moduleDist
                        yFinal=segment[1][1]+yM*RTDist
                        
                        #find x-gradient between module and extrapolate to RT
                        xM=(segment[1][0]-segment[0][0])/moduleDist
                        xFinal=segment[1][0]+xM*RTDist

                        #loop over RT hits to see if any agree within tolerance
                        for k in Hits[2]:
                                r= math.sqrt( (k[0]-xFinal)**2 + (k[1]-yFinal)**2 )
                                #print r,stripTolerance*pitch
                                if r<2*stripTolerance*pitch:
                                        #print
                                        tracks.append([i,j,k])

                #filter tracks by chi2 if more than one found

        #filter tracks by proximity to each other to check 
                
        return tracks
                       
def SumRadialSeparation(hits):

        rSum=0.0
        xmean=0.0
	ymean=0.0

        for i in hits:
                xmean+=i[0]
                ymean+=i[1]
        xmean/=len(hits)        
        ymean/=len(hits)        

        for i in hits:
                r= math.sqrt( (i[0]-xmean)**2  + (i[1]-ymean)**2)
                rSum+=r
        return rSum

def ReconstructTracks4Planes(Hits,trackTolerance,pitch):

        nTrackerModules=len(Hits)
        RecoTracks=[]
        
        while len(Hits[0])>0 and len(Hits[1])>0 and len(Hits[2])>0 and len(Hits[3])>0:
                bestChi2=100000
                bestIPos=-1
                bestJPos=-1
                bestKPos=-1
                bestLPos=-1

                #loop over all combos from each tracker module: find chi2 to linear fit
                for i in range(len(Hits[0])):
                        for j in range(len(Hits[1])):
                                for k in range(len(Hits[2])):
                                        for l in range(len(Hits[3])):
                                                hitI=Hits[0][i]
                                                hitJ=Hits[1][j]
                                                hitK=Hits[2][k]
                                                hitL=Hits[3][l]
                                                chi2=SumRadialSeparation([hitI,hitJ,hitK,hitL])
                                                if chi2<bestChi2 and chi2>0.0001:
                                                        bestIPos=i
                                                        bestJPos=j
                                                        bestKPos=k
                                                        bestLPos=l
                                                        bestChi2=chi2

                if bestChi2<trackTolerance : 
                        #Append hit positions to track list
                        RecoTracks.append([Hits[0][bestIPos],Hits[1][bestJPos],Hits[2][bestKPos],Hits[3][bestLPos]])
                        #delete entries from Hits
                        del Hits[0][bestIPos]
                        del Hits[1][bestJPos]
                        del Hits[2][bestKPos]
                        del Hits[3][bestLPos]
                else:
                        break

        return RecoTracks

def ReconstructTracks(Hits,trackTolerance,pitch,MaxNTracks,restrictNTracks,TrackerZ,stripTolerance,beamSpread ):

        if len(Hits)==2:
                AllTracks=ReconstructTracks2Planes(Hits,trackTolerance,pitch,TrackerZ)
        elif len(Hits)==3:
                AllTracks=ReconstructTracks3Planes(Hits,trackTolerance,pitch,TrackerZ,stripTolerance,beamSpread) 
        elif len(Hits)==4:
                AllTracks=ReconstructTracks4Planes(Hits,trackTolerance,pitch)
        else:
                print("Can only cope with 2 or 4 tracker modules!")
                quit()

        # can cut on expected number of tracks based on number of strips fired
        if restrictNTracks and len(AllTracks)>MaxNTracks:
                SelectedTracks=AllTracks[0:MaxNTracks]
        else:
                SelectedTracks=AllTracks
                
        return SelectedTracks
        
def DrawTrackMap(name,Tracks,XY,xmax):

        truth=TH2F("Truth","Hit Locations",1000,-xmax,xmax,1000,-xmax,xmax)
        truth.SetMarkerColor(1)
        truth.SetMarkerSize(2)
        truth.SetMarkerStyle(4)
        tracker1=TH2F("Tracker1","Hit Locations",1000,-xmax,xmax,1000,-xmax,xmax)
        tracker1.SetMarkerColor(2)
        tracker1.SetMarkerSize(2)
        tracker1.SetMarkerStyle(2)
        tracker2=TH2F("Tracker2","Hit Locations",1000,-xmax,xmax,1000,-xmax,xmax)
        tracker2.SetMarkerColor(4)
        tracker2.SetMarkerSize(2)
        tracker2.SetMarkerStyle(5)
        tracker3=TH2F("Tracker3","Hit Locations",1000,-xmax,xmax,1000,-xmax,xmax)
        tracker3.SetMarkerColor(4)
        tracker3.SetMarkerSize(2)
        tracker3.SetMarkerStyle(32)
        tracker4=TH2F("Tracker4","Hit Locations",1000,-xmax,xmax,1000,-xmax,xmax)
        tracker4.SetMarkerColor(4)
        tracker4.SetMarkerSize(2)
        tracker4.SetMarkerStyle(26)

        trackerArray=[]
        trackerArray.append(tracker1)
        trackerArray.append(tracker2)
        trackerArray.append(tracker3)
        trackerArray.append(tracker4)
        
        for coord in XY:
                truth.Fill(coord[0],coord[1])
        
        for track in Tracks:
                for i in range(len(track)):
                        trackerArray[i].Fill(track[i][0],track[i][1])

        canvas1 = TCanvas( 'c1', "mycanvas", 200, 10, 700, 500 )
        trackerArray[0].Draw()
        trackerArray[1].Draw("SAME")
        trackerArray[2].Draw("SAME")
        trackerArray[3].Draw("SAME")
        truth.Draw("SAME")

        legend=TLegend(0.1,0.7,0.48,0.9)
        legend.AddEntry(truth,"Truth","p")
        legend.AddEntry(trackerArray[0],"Tracker1","p")
        legend.AddEntry(trackerArray[1],"Tracker2","p")
        legend.AddEntry(trackerArray[2],"Tracker3","p")
        legend.AddEntry(trackerArray[3],"Tracker4","p")
        legend.Draw("SAME")
        canvas1.Write(name)

        

def GetTrackSeparation(trueTrack,recoTrack,effTolerance):

        passed=True
        for trueHit,recoHit in zip(trueTrack,recoTrack):
                #for every point on track check if the true and reconstructed hits agree
                separation=math.sqrt((trueHit[0]-recoHit[0])**2 + (trueHit[1]-recoHit[1])**2)
                if separation>effTolerance:
                        passed=False
                        break
        return passed
                
def GetTrackEfficiency(effTolerance,XY,mXmY,RecoTracks,TrackerZ,pitch):

        #if only one module, tracks are just hits for one module
        if len(TrackerZ)==1:
               return GetEfficiency(RecoTracks,XY,pitch,effTolerance)

        
        TrackerPositions=[]
        for z in TrackerZ:
                TrackerPositions.append(sum(z)/len(z))

        nTrueFound=0.0
        
        #loop over all true protons
        for coord,direction in zip(XY,mXmY):

                #create a track for the true proton 
                trueTrack=[]
                passed=False
                for position in TrackerPositions:
                        trueTrack.append([coord[0]+direction[0]*position,coord[1]+direction[1]*position])
                #loop over reco tracks and see if any match at each module within hit tolerance
                for recoTrack in RecoTracks:
                        if GetTrackSeparation(trueTrack,recoTrack,effTolerance):
                                nTrueFound+=1.0
                                passed=True
                                break

                #if passed==False:
                #        print "missed track positions=",trueTrack
    

        return nTrueFound/(float)(len(XY))


def GetCombinedEfficiency(stripTolerance,coords,directions,CombinedRecoTracks,ZMeans,pitch):

        nTrueFound=0.0
        
        #loop over all true protons
        for frontCoord,frontDir,rearCoord,rearDir in zip(coords[0],directions[0],coords[1],directions[1]):

                #create a track for the true proton 
                trueTrackFront=[]
                trueTrackRear=[]
                
                passed=False
                for position in ZMeans[0:2]:
                        trueTrackFront.append([frontCoord[0]+frontDir[0]*position,frontCoord[1]+frontDir[1]*position])
                for position in ZMeans[2:5]:
                        trueTrackRear.append([rearCoord[0]+rearDir[0]*position,rearCoord[1]+rearDir[1]*position])

                #loop over reco tracks and see if any match at each module within hit tolerance
                for recoTrack in CombinedRecoTracks:
                        
                        if GetTrackSeparation(trueTrackFront,recoTrack[0],stripTolerance[0]*pitch) and GetTrackSeparation(trueTrackRear,recoTrack[1],stripTolerance[1]*pitch):
                                nTrueFound+=1.0
                                passed=True
                                break

                #if passed==False and len(CombinedRecoTracks)>0 :
                #        print "Strip tolerances",stripTolerance
                #        print "True Track=",trueTrackFront,trueTrackRear
                #        for track in CombinedRecoTracks:
                #                print "RecoTrackFront=",track[0][0][0:2],track[0][1][0:2],track[1][0][0:2],track[1][1][0:2],track[1][2][0:2]
                #        print ""
        return nTrueFound/(float)(len(coords[0]))


def WriteTracks(outfilename,RecoTracks,ZMeans,loop):

        outString=(str)(loop)       
        f= open(outfilename,"a+")
        for module in range(len(ZMeans)):
                for track in range(5):
                        if len(RecoTracks)>track:
                                hit=RecoTracks[track][module]
                                x=hit[0]
                                y=hit[1]
                                z=ZMeans[module]
                        else:
                                x=-1
                                y=-1
                                z=-1
                        outString+=" "+(str)(x)+" "+(str)(y)+" "+(str)(z)

                        
        f.write(outString+"\n")
        f.close()

def ReadJohnFile(inFile):
        lastTimestamp=-1

        stripshift=512


        with open(inFile,'r') as f:
                line=f.readline()
                while line:
                        values=line.split()
                        timestamp=values[0]
                        #check if any hits in entry, check if there are same number of hits in each of the trackers- needed for reconstruction
                        if tracker1Hits and len(tracker1Hits)==len(tracker2Hits) and len(tracker1Hits)==len(tracker3Hits) and len(tracker1Hits)==len(tracker4Hits):

                                hitsByTimestamp.append([tracker1Hits,tracker2Hits,tracker3Hits,tracker4Hits])

                                lastTimestamp=timestamp
                                tracker1Hits=[[],[],[]]
                                tracker2Hits=[[],[],[]]
                                tracker3Hits=[[],[],[]]
                                tracker4Hits=[[],[],[]]

                                
                        #loop over 4 data entries and only add them if there's a hit in all planes- format of data is V,X,U (i.e. middle value is the 0 degrees plane)
                        for i in range (1,5):
                                j=i+4
                                k=i+8
                                l=i+12
                                #check value is valid
                                if ( values[i]!='-1' and values[i+4]!='-1' and values[i+8]!='-1'):
                                       tracker1Hits[2].append((int)(values[i])-stripshift)
                                       tracker1Hits[0].append((int)(values[i+4])-stripshift)
                                       tracker1Hits[1].append((int)(values[i+8])-stripshift)

                                #check value is valid
                                if ( values[j]!='-1' and values[j+4]!='-1' and values[j+8]!='-1'):
                                       tracker2Hits[2].append((int)(values[j])-stripshift)
                                       tracker2Hits[0].append((int)(values[j+4])-stripshift)
                                       tracker2Hits[1].append((int)(values[j+8])-stripshift)

                                #check value is valid
                                if ( values[k]!='-1' and values[k+4]!='-1' and values[k+8]!='-1'):
                                       tracker3Hits[2].append((int)(values[k])-stripshift)
                                       tracker3Hits[0].append((int)(values[k+4])-stripshift)
                                       tracker3Hits[1].append((int)(values[k+8])-stripshift)

                                #check value is valid
                                if ( values[l]!='-1' and values[l+4]!='-1' and values[l+8]!='-1'):
                                       tracker4Hits[2].append((int)(values[l])-stripshift)
                                       tracker4Hits[0].append((int)(values[l+4])-stripshift)
                                       tracker4Hits[1].append((int)(values[l+8])-stripshift)

                        line=f.readline()

        hXStrips.Write()
        hUStrips.Write()
        hVStrips.Write()
        return hitsByTimestamp

def CheckInsideDetectorArea(Hits,pitch,size,angles):

        acceptedHits=[]

        maxStrip=size*10000.0/(2*pitch)
        
        for hit in Hits:
                coord=[[hit[0],hit[1]]]
                Strips=GetStripCoOrds(coord,pitch,angles)[0]
                passed=True
                for strip in Strips:
                        if abs(strip)>maxStrip:
                                passed=False
                                break
                if passed==True:
                     acceptedHits.append(hit)
                else:
                        print "Failed acceptance test"
                        print Strips

        return acceptedHits

hDisplacement=TH1F("hDisplacement","Displacment (um)",50,-20000,20000)
hAngle=TH1F("hAngle","Angle (rad)",50,-0.1,0.1)

def GenerateMSAngle(phantomdepth,ebeam):

        #ebeam = kinetic energy in MeV
        
        m=938.272 #MeV
        E=ebeam+m #K.E. plus rest mass
        p=math.sqrt(E**2 -m**2)
        v=p/E
        vp=v*p
        
        #silicon
        #x=155 #um
        #X0=9.37*10000 #radiation length, convert from cm to um

        #water
        x=phantomdepth*10000 #um
        X0=36.08*10000 #radiation length, convert from cm to um
        distance=math.sqrt(x/X0)
        
        #pdg- theta0 is the sigma for gaussian distribution of the angles
        theta0=(13.6/vp)*distance*(1+0.038*math.log(distance)) 
        
        a1=random.gauss(0.0,1.0)
        a2=random.gauss(0.0,1.0)

        #from pdg for MS, way to generate correlated angular and spacial shifts
        spatialShift=a1*x*theta0/math.sqrt(12.0) + a2*x*theta0/2.0 #um
        thetaShift=a2*theta0 #rad

        hDisplacement.Fill(spatialShift)
        hDisplacement.Write("",TObject.kOverwrite)

        hAngle.Fill(thetaShift)
        hAngle.Write("",TObject.kOverwrite)

        #print spatialShift
        return spatialShift,thetaShift

        
def ApplyMultipleScattering(XY,mXmY,phantomdepth,energy,ZPosition):

        XY_MS,mXmY_MS=[],[]
        halfPhantomWidth=phantomdepth*10000.0/2.0

        for xy,mxmy in zip(XY,mXmY):

               spatialShiftX,thetaShiftX= GenerateMSAngle(phantomdepth,energy)
               spatialShiftY,thetaShiftY= GenerateMSAngle(phantomdepth,energy)

               #should we scale theta by sqrt(2) to account for theta0= sqrt(thetaX**2+thetaY**2)

               #small angle approximation- tan(theta) ~ theta so can just add
               mXmY_MS.append([mxmy[0]+thetaShiftX,mxmy[1]+thetaShiftY])

               currentX=xy[0]+(mxmy[0]*(ZPosition+halfPhantomWidth))+spatialShiftX
               currentY=xy[1]+(mxmy[1]*(ZPosition+halfPhantomWidth))+spatialShiftY
               #currentX=xy[0]
               #currentY=xy[1]

               X0=currentX-mXmY_MS[-1][0]*(ZPosition+halfPhantomWidth)
               Y0=currentY-mXmY_MS[-1][1]*(ZPosition+halfPhantomWidth)

               XY_MS.append([X0,Y0])

               
        return XY_MS,mXmY_MS



def MatchTrack(RecoTracks,RearRecoTracks,ZMeans,tolerance):

        #matches tracks from front and rear of patient by calculating difference in projected path at centre of phantom
        combinedTracks=[]
        
        phantomCentre=(ZMeans[2]+ZMeans[1])/2

   
        while len(RecoTracks)>0 and len(RearRecoTracks)>0:
                bestRadius=1000000000
                bestI=10000
                bestJ=10000
                for i in range(len(RecoTracks)):
                        for j in range(len(RearRecoTracks)):

                                front=RecoTracks[i]
                                rear=RearRecoTracks[j]
                                
                                frontFit=front[-1] #[cX,mX,cY,mY]
                                rearFit=rear[-1]

                                xFront=frontFit[1]*phantomCentre+frontFit[0]
                                yFront=frontFit[3]*phantomCentre+frontFit[2]

                                xRear=rearFit[1]*phantomCentre+rearFit[0]
                                yRear=rearFit[3]*phantomCentre+rearFit[2]

                                radius=math.sqrt( (xFront-xRear)**2 + (yFront-yRear)**2 )

                                if radius<bestRadius:
                                        bestRadius=radius
                                        bestI=i
                                        bestJ=j

                                
                if bestRadius<tolerance:
                        combinedTracks.append([RecoTracks[bestI],RearRecoTracks[bestJ]])
                        #delete best positions tracks
                        del RecoTracks[bestI]
                        del RearRecoTracks[bestJ]
                else:
                        break
                        
        return combinedTracks


#100% efficienct version for debugging
#def MatchTrack(RecoTracks,RearRecoTracks,ZMeans,tolerance):
#
#        #matches tracks from front and rear of patient by calculating difference in projected path at centre of phantom
#        combinedTracks=[]
#        
#        phantomCentre=(ZMeans[2]+ZMeans[1])/2
#        for i in range(len(RecoTracks)):
#                for j in range(len(RearRecoTracks)):
#
#                        front=RecoTracks[i]
#                        rear=RearRecoTracks[j]
#
#                        frontFit=front[-1] #[cX,mX,cY,mY]
#                        rearFit=rear[-1]
#
#                        xFront=frontFit[1]*phantomCentre+frontFit[0]
#                        yFront=frontFit[3]*phantomCentre+frontFit[2]
#
#                        xRear=rearFit[1]*phantomCentre+rearFit[0]
#                        yRear=rearFit[3]*phantomCentre+rearFit[2]
#
#                        radius=math.sqrt( (xFront-xRear)**2 + (yFront-yRear)**2 )
#
#                        if radius<tolerance: 
#                                combinedTracks.append([RecoTracks[i],RearRecoTracks[j]])
#                        #else:
#                        #        print "RejectedTrack=",RecoTracks[i][0][0:2],RecoTracks[i][1][0:2],RearRecoTracks[j][0][0:2],RearRecoTracks[j][1][0:2],RearRecoTracks[j][2][0:2]
#                        #        print "Radius was ",radius
#                        #        print xFront,yFront,xRear,yRear
#                        
#        return combinedTracks
#
