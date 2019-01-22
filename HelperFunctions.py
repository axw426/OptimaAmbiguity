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

#convert points into a strip number for each layer assuming pitch of 100um
def GetStripCoOrds(XY,pitch,rawAngle):
	StripX=[]
	StripU=[]
	StripV=[]
	for coord in XY:
		StripX.append((int)((coord[0]-0.5*pitch)/pitch))
		StripU.append((int)(((coord[0]-0.5*pitch)*math.cos(math.radians(rawAngle)) - (coord[1]-0.5*pitch)*math.sin(math.radians(rawAngle))) /pitch))
		StripV.append((int)(((coord[0]-0.5*pitch)*math.cos(math.radians(-rawAngle)) - (coord[1]-0.5*pitch)*math.sin(math.radians(-rawAngle))) / pitch))
		#print(StripX[-1],StripU[-1],StripV[-1])
	return StripX,StripU,StripV	

#convert points into a strip number for each layer assuming pitch of 100um
def GetPixelCoOrds(XY,pitch):
        coords=[]
	for coord in XY:
		cellX=(int)((coord[0]-(pitch/2.0))/pitch)
		cellY=(int)((coord[1]-(pitch/2.0))/pitch)
                coords.append([cellX,cellY])
                
	return coords


def GetStripCoOrds4Planes(XY,pitch,rawAngle):
	StripX=[]
	StripU=[]
	StripV=[]
	StripY=[]
	for coord in XY:
		StripX.append((int)((coord[0]-50)/pitch))
		StripU.append((int)(((coord[0]-50)*math.cos(math.radians(rawAngle)) - (coord[1]-50)*math.sin(math.radians(rawAngle))) /pitch))
		StripV.append((int)(((coord[0]-50)*math.cos(math.radians(2*rawAngle)) - (coord[1]-50)*math.sin(math.radians(2*rawAngle))) / pitch))
		StripY.append((int)(((coord[0]-50)*math.cos(math.radians(3*rawAngle)) - (coord[1]-50)*math.sin(math.radians(3*rawAngle))) / pitch))
		#print(StripX[-1],StripU[-1],StripV[-1])
	return StripX,StripU,StripV,StripY	


def FindMandC(U,theta,pitch):
     	xPrime=(U*pitch)#+pitch/2.0
	
	m=-1*math.tan(math.radians(90-theta))
	c=xPrime/math.cos(math.radians(90-theta))

        return m,c

def FindIntersect(X,U,theta,pitch):

	#in X frame, x*100+50 to get to middle of strip
	#U line described by y=mx+c, where: m=-tan(90-Theta),c= X'/cos(90-theta)

	x=(X*pitch)#+pitch/2.0
	m,c=FindMandC(U,theta,pitch)
	y=(m*x)+c

        #print("M=",m,"c=",c,"x=",x,"X'=",xPrime,"y=",y)
	
	return [x,y]
	
def FindUVIntersect(U,V,uAngle,vAngle,pitch):

	#define y=mx+c for U and V in frame of X layer
        m1,c1=FindMandC(U,uAngle,pitch)
	
        m2,c2=FindMandC(V,vAngle,pitch)
	
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

#def CheckStripHalf(xmean,ymean,X,U,V,rawAngle,pitch):
 #       passed=True
  #      if ytrue*ymean<0:
   #             passed=False
#
 #       yU=xmean*tan(rawAngle)
  #      if ytrue*(ymean-yU)<0:
  #              passed=False
#
 #       yV=xmean*tan(-rawAngle)
  #      if ytrue*(ymean-yV)<0:
   #             passed=False
#
 #       return passed


#cycle through all strip combinations and see which create an overlap- will want to define a central line for each strip and work out if they intersect 
def FindOverlaps(StripX,StripU,StripV,pitch,rawAngle,tolerance):

        allHits=[]
        
	#find the Y coordinate where the X and U or V intercept and see if they are close
	for X in StripX:
		for U in StripU:
			for V in StripV:
                                xycoords=[]
				xycoords.append(FindIntersect(X,U,-rawAngle,pitch))
				xycoords.append(FindIntersect(X,V,rawAngle,pitch))
				xycoords.append(FindUVIntersect(U,V,-rawAngle,rawAngle,pitch))
				passed,xmean,ymean,rValues=CheckProximity(xycoords,tolerance,pitch)
                                if passed==True:
                                        allHits.append([xmean,ymean,X,U,V,rValues])
                                #else:
                                 #       print(X,U,V,xycoords[0],xycoords[1],xycoords[2])

	return allHits

def FindOverlaps4Planes(StripX,StripU,StripV,StripY,pitch,rawAngle,tolerance):

        allHits=[]
        
	#find the Y coordinate where the X and U or V intercept and see if they are close
	for X in StripX:
		for U in StripU:
			for V in StripV:
			        for Y in StripY:
                                        xycoords=[]
                                        xycoords.append(FindIntersect(X,U,-rawAngle,pitch))
				        xycoords.append(FindIntersect(X,V,-rawAngle*2,pitch))
				        xycoords.append(FindIntersect(X,Y,-rawAngle*3,pitch))
                                        
                                        xycoords.append(FindUVIntersect(U,V,-rawAngle,-rawAngle*2,pitch))
                                        xycoords.append(FindUVIntersect(U,Y,-rawAngle,-rawAngle*3,pitch))
                                        xycoords.append(FindUVIntersect(V,Y,-rawAngle*2,-rawAngle*3,pitch))
                                        
				        passed,xmean,ymean,rValues=CheckProximity(xycoords,tolerance,pitch)
                                        if passed==True:
                                                allHits.append([xmean,ymean,X,U,V,Y,rValues])

	return allHits

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

def PlotHitMap(name,allHits,XY,StripX,StripU,StripV,pitch,size,loop,theta):

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
        for X,U,V in zip(StripX,StripU,StripV):
                #draw line for X plane
                x=(X*pitch)#+pitch/2.0
                myline=TLine(x,-xmax,x,xmax)
                myline.SetLineColor(1)
                myline.SetLineWidth(1)
                linearray.append(myline)
                linearray.append(CreateTLines(xmax,U,theta,pitch))
                linearray.append(CreateTLines(xmax,V,-theta,pitch))

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

def PlotHitMap4Planes(name,allHits,XY,StripX,StripU,StripV,StripY,pitch,size,loop,theta):

        #print("Entering plot hit map 4 planes")
        #print len(StripY)
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
        for X,U,V,Y in zip(StripX,StripU,StripV,StripY):
         #       print (X,U,V,Y)
                #draw line for X plane
                x=(X*pitch)#+pitch/2.0
                myline=TLine(x,-xmax,x,xmax)
                myline.SetLineColor(1)
                myline.SetLineWidth(1)
                linearray.append(myline)
                linearray.append(CreateTLines(xmax,U,theta,pitch))
                linearray.append(CreateTLines(xmax,V,2*theta,pitch))
                linearray.append(CreateTLines(xmax,Y,3*theta,pitch))

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
#        duplicatePositions=[]
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

        nExpected=len(XY)
        nFound=0.0

        for true in XY:
                passed=False
                for hit in allHits:
                        if math.sqrt((true[0]-hit[0])**2 + (true[1]-hit[1])**2)<pitch*tolerance:
                                passed=True
                if passed ==True:
                        nFound+=1


        return nFound/nExpected

def GetMaxR(hit):

        rMin=-1.0
        for r in hit[5]:
              if r>rMin:
                      rMin=r

        return rMin

def GetSumR(hit):

        rSum=0.0
        for r in hit[5]:
                rSum+=r

        return rSum
        
def GetPixelArea(hit,rawAngle,pitch):

        #find crossing points for XU,XV,and UV and assume triangular
        rSum=0.0

        XU=FindIntersect(hit[2],hit[3],-rawAngle,pitch)
	XV=FindIntersect(hit[2],hit[4],rawAngle,pitch)
	UV=FindUVIntersect(hit[3],hit[4],-rawAngle,rawAngle,pitch)

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
                        if hitA[2]==hitB[2] and GetPixelArea(hitA,rawAngle,pitch)>GetPixelArea(hitB,rawAngle,pitch):
                                passed=False
                if passed==True:
                        XAcceptedHits.append(hitA)

        #check unique U
        for hitA in XAcceptedHits:
                passed=True
                for hitB in XAcceptedHits:
                        if hitA[3]==hitB[3] and GetPixelArea(hitA,rawAngle,pitch)>GetPixelArea(hitB,rawAngle,pitch):
                                passed=False
                if passed==True:
                       XUAcceptedHits.append(hitA)           

        #check uniqueV
        for hitA in XUAcceptedHits:
                passed=True
                for hitB in XUAcceptedHits:
                        if hitA[4]==hitB[4] and GetPixelArea(hitA,rawAngle,pitch)>GetPixelArea(hitB,rawAngle,pitch):
                                passed=False
                if passed==True:
                       XUVAcceptedHits.append(hitA)  

        return XUVAcceptedHits
