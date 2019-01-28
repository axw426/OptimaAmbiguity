import random
import math
from ROOT import gROOT, TCanvas, TF1, TF2, TH2F, TFile, TLine, TH1F


def SetSeed(seed):
	random.seed(seed)

def GetRandomXY(nProtons,size):
	XY=[]
	for x in range(nProtons):
		X=random.uniform(-size*5000.0,size*5000.0) #10000 um per cm
		Y=random.uniform(-size*5000.0,size*5000.0) #10000 um per cm
	#	X=random.gauss(0.0,6000)
	#	Y=random.gauss(0.0,6000)
                XY.append([X,Y])
	        #print X,Y
	return XY

def ConvertXYToStrip(coord,pitch,rawAngle):

        return (int)(((coord[0]-50)*math.cos(math.radians(rawAngle)) - (coord[1]-50)*math.sin(math.radians(rawAngle))) /pitch)

def CheckStripHalf(coord,pitch,rawAngle):

        y= (int)(((coord[0]-50)*math.sin(math.radians(rawAngle)) + (coord[1]-50)*math.cos(math.radians(rawAngle))) /pitch)

        return y
        #if y>=0:
        #        return 1.0
        #else:
        #        return -1.0
#convert points into a strip number for each layer assuming pitch of 100um
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
               # print("3 plane tracker")
                for coord in XY:
                        StripV.append(ConvertXYToStrip(coord,pitch,angles[2]))
                        StripHalfV.append(CheckStripHalf(coord,pitch,angles[2]))
                Strips=[StripX,StripU,StripV,StripHalfX,StripHalfU,StripHalfV]
                        
        if len(angles)>3:
               # print("4 plane tracker")
                for coord in XY:
                        StripY.append(ConvertXYToStrip(coord,pitch,angles[3]))
                        StripHalfY.append(CheckStripHalf(coord,pitch,angles[3]))
	        Strips= [StripX,StripU,StripV,StripY,StripHalfX,StripHalfU,StripHalfV,StripHalfY]	

        return Strips
#convert points into a strip number for each layer assuming pitch of 100um
def GetPixelCoOrds(XY,pitch):
        coords=[]
	for coord in XY:
		cellX=(int)((coord[0]-(pitch/2.0))/pitch)
		cellY=(int)((coord[1]-(pitch/2.0))/pitch)
                coords.append([cellX,cellY])
                
	return coords


def FindMandC(U,theta,pitch):
     	xPrime=(U*pitch)#+pitch/2.0
	
	m=-1*math.tan(math.radians(90-theta))
	c=xPrime/math.cos(math.radians(90-theta))

        return m,c

def FindIntersect(X,U,theta,pitch):

	#in X frame, x*100 to get to middle of strip
	#U line described by y=mx+c, where: m=-tan(90-Theta),c= X'/cos(90-theta)

	x=(X*pitch)#+pitch/2.0
	m,c=FindMandC(U,theta,pitch)
	y=(m*x)+c

        #print("M=",m,"c=",c,"x=",x,"X'=",xPrime,"y=",y)
	
	return [x,y]
	
def FindUVIntersect(U,V,uAngle,vAngle,pitch):

        if uAngle==0:
                return FindIntersect(U,V,-vAngle,pitch)

        if vAngle==0:
                return FindIntersect(V,U,-uAngle,pitch)

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
        #print(len(xycoords))
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
        rawmeasured=TH2F(name+(str)(loop),"Hit Locations",1000,-xmax,xmax,1000,-xmax,xmax)
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
                        #draw line for X plane
                        x=(X*pitch)#+pitch/2.0
                        myline=TLine(x,-xmax,x,xmax)
                        myline.SetLineColor(1)
                        myline.SetLineWidth(1)
                        linearray.append(myline)
                        linearray.append(CreateTLines(xmax,U,angles[1],pitch))
                        linearray.append(CreateTLines(xmax,V,angles[2],pitch))

        if nPlanes==4:
                for X,U,V,Y in zip(Strips[0],Strips[1],Strips[2],Strips[3]):
                        #draw line for X plane
                        x=(X*pitch)#+pitch/2.0
                        myline=TLine(x,-xmax,x,xmax)
                        myline.SetLineColor(1)
                        myline.SetLineWidth(1)
                        linearray.append(myline)
                        linearray.append(CreateTLines(xmax,U,angles[1],pitch))
                        linearray.append(CreateTLines(xmax,V,angles[2],pitch))
                        linearray.append(CreateTLines(xmax,Y,angles[3],pitch))
                        
        for line in linearray:
                line.Draw("SAME")
                               

        #draw where protons really hit
        truth=TH2F("truth"+(str)(loop),"Hit Locations",1000,-xmax,xmax,1000,-xmax,xmax)
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
                                
                if passed==True:
                        refinedHits.append(i)

        return refinedHits

def GetEfficiency(allHits,XY,pitch,tolerance):

        #check if there is a reconstructed hit within tolerance of the true hit positions
        nExpected=len(XY)
        nFound=0.0

        for true in XY:
                passed=False
                for reco in allHits:
                        if math.sqrt((true[0]-reco[0])**2 + (true[1]-reco[1])**2)<pitch*tolerance:
                                passed=True
                if passed ==True:
                        nFound+=1.0


        return nFound/(float)(nExpected)

def GetPixelArea(hit,angles,pitch):

        #find crossing points for XU,XV,and UV and assume triangular
        rSum=0.0

        Strips=hit[2]
	XU=FindUVIntersect(Strips[0],Strips[1],angles[0],angles[1],pitch)
	XV=FindUVIntersect(Strips[0],Strips[2],angles[0],angles[2],pitch)
	UV=FindUVIntersect(Strips[1],Strips[2],angles[1],angles[2],pitch)

        base=XU[1]-XV[1]
        height=XU[0]-UV[0]

        #print(XU[1]-XV[1])
        return 0.5*base*height
        

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
                if passed==True:
                        XAcceptedHits.append(hitA)

        #check unique U
        for hitA in XAcceptedHits:
                passed=True
                for hitB in XAcceptedHits:
                        if hitA[2][1]==hitB[2][1] and GetPixelArea(hitA,rawAngle,pitch)>GetPixelArea(hitB,rawAngle,pitch):
                                passed=False
                if passed==True:
                       XUAcceptedHits.append(hitA)           

        #check uniqueV
        for hitA in XUAcceptedHits:
                passed=True
                for hitB in XUAcceptedHits:
                        if hitA[2][2]==hitB[2][2] and GetPixelArea(hitA,rawAngle,pitch)>GetPixelArea(hitB,rawAngle,pitch):
                                passed=False
                if passed==True:
                       XUVAcceptedHits.append(hitA)  

        return XUVAcceptedHits

def FindRadialSeparation(i,j):

        return math.sqrt((i[0]-j[0])**2 + (i[1]-j[1])**2)

def ReconstructTracks(Hits,tolerance,pitch):

        RecoTracks=[]
        
        while len(Hits[0])>0 and len(Hits[1])>0:
                #print("Lenghts=",len(Hits[0]),len(Hits[1]))
                minimumSeparation=100000
                bestIPos=-1
                bestJPos=-1
                for i in range(len(Hits[0])):
                        for j in range(len(Hits[1])):
                                hitI=Hits[0][i]
                                hitJ=Hits[1][j]
                                r=FindRadialSeparation(hitI,hitJ)
                                if r<minimumSeparation:
                                        bestIPos=i
                                        bestJPos=j
                                        minimumSeparation=r
               # print minimumSeparation
                if minimumSeparation<tolerance*pitch:
                        #Append hit positions to track list
                        RecoTracks.append([Hits[0][bestIPos],Hits[1][bestJPos]])
                        #delete entries from Hits
                        del Hits[0][bestIPos]
                        del Hits[1][bestJPos]
                else:
                        break

        return RecoTracks
