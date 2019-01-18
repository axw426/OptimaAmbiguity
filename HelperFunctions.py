import random
import math
from ROOT import gROOT, TCanvas, TF1, TF2, TH2F, TFile, TLine, TH1F


def SetSeed(seed):
	random.seed(seed)

def GetRandomXY(nProtons,size):
	XY=[]
	for x in range(nProtons):
		X=random.uniform(-size*5000,size*5000) #10000 um per cm
		Y=random.uniform(-size*5000,size*5000)
                XY.append([X,Y])
	
	return XY

#convert points into a strip number for each layer assuming pitch of 100um
def GetStripCoOrds(XY,pitch,rawAngle):
	StripX=[]
	StripU=[]
	StripV=[]
	for coord in XY:
		StripX.append((int)(coord[0]/pitch))
		StripU.append((int)((coord[0]*math.cos(math.radians(rawAngle)) - coord[1]*math.sin(math.radians(rawAngle))) /pitch))
		StripV.append((int)((coord[0]*math.cos(math.radians(-rawAngle)) - coord[1]*math.sin(math.radians(-rawAngle))) / pitch))
		#print(StripX[-1],StripU[-1],StripV[-1])
	return StripX,StripU,StripV	

def FindMandC(U,theta,pitch):
     	xPrime=(U*pitch)+pitch/2.0
	
	m=-1*math.tan(math.radians(90-theta))
	c=xPrime/math.cos(math.radians(90-theta))

        return m,c

def FindIntersect(X,U,theta,pitch):
	#in X frame, x*100+50 to get to middle of strip
	#U line described by y=mx+c, where: m=-tan(90-Theta),c= X'/cos(90-theta)
	x=(X*pitch)+pitch/2.0
	m,c=FindMandC(U,theta,pitch)
	y=(m*x)+c
	#print("M=",m,"c=",c,"x=",x,"X'=",xPrime,"y=",y)
	
	return [x,y]
	
def FindUVIntersect(U,V,uAngle,vAngle,pitch):

	#define y=mx+c for U and V in frame of X layer
	x_U=(U*pitch)+pitch/2.0
	m1=-1*math.tan(math.radians(90-uAngle))
	c1=x_U/math.cos(math.radians(90-uAngle))
	
	x_V=(V*pitch)+pitch/2.0
	m2=-1*math.tan(math.radians(90-vAngle))
	c2=x_V/math.cos(math.radians(90-vAngle))
	
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

def checkForDuplicates(xmean,ymean,allHits,tolerance):
        
        for hit in allHits:
                r=math.sqrt((xmean-hit[0])**2 + (ymean-hit[1])**2)
                if r<tolerance*pitch:
                        return True
        return False

	
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

	return allHits

def CreateHitMap(myHitMap,StripX,StripU,StripV,pitch,rawAngle,tolerance):

        nHits=0
	#find the Y coordinate where the X and U or V intercept and see if they are close
	for X in StripX:
		for U in StripU:
			for V in StripV:
				x1,y1=FindIntersect(X,U,-rawAngle,pitch)
				x2,y2=FindIntersect(X,V,rawAngle,pitch)
				x3,y3,m1,c1,m2,c2=FindUVIntersect(U,V,-rawAngle,rawAngle,pitch)
				passed,xmean,ymean=CheckProximity(x1,x2,x3,y1,y2,y3,tolerance,pitch)
                                #print x1,y1,xmean,ymean
                                #print m1,m2
                                if passed==True:
                                        myHitMap.Fill(xmean,ymean)
                                        nHits+=1
        return nHits



def ReadXUVStripCoOrds(inFile):
        hitsByTimestamp=[]
        lastTimestamp=-1
        xcoords=[]
        ucoords=[]
        vcoords=[]

        hXStrips=TH1F("hXStrips","n Strip Positions X",1024,-0.5,1023.5)
        hUStrips=TH1F("hUStrips","n Strip Positions U",1024,-0.5,1023.5)
        hVStrips=TH1F("hVStrips","n Strip Positions V",1024,-0.5,1023.5)

        #outf=open("myoutput.txt","w+")

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
                                #outf=open("myoutput.txt","a+")
                                #outf.write((str)(timestamp)+", "+(str)(vcoords)+", "+(str)(xcoords)+", "+(str)(ucoords)+"\n")
                                #outf.close()
                                xcoords[:]=[(var-512) for var in xcoords]
                                ucoords[:]=[(var-512) for var in ucoords]
                                vcoords[:]=[(var-512) for var in vcoords]
                                hitsByTimestamp.append([xcoords,ucoords,vcoords])

                                #reset vectors
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

                                #elif values[i] !='-1' and values[i+4] !='-1':
                                #        hVStrips.Fill((int)(values[i]))
                                #        hXStrips.Fill((int)(values[i+4]))

                                #elif values[i] !='-1' and values[i+8] !='-1':
                                #        hVStrips.Fill((int)(values[i]))
                                #        hUStrips.Fill((int)(values[i+8]))

                                #elif values[i+4] !='-1' and  values[i+8] !='-1':
                                #        hXStrips.Fill((int)(values[i+4]))
                                #        hUStrips.Fill((int)(values[i+8]))


                                #check value is valid
                                if ( values[i]!='-1' and values[i+4]!='-1' and values[i+8]!='-1'):

                                        #check isn't a duplicate value (two adjacent strips...), lazy for now, just take first of two adjacent
                                        if(xcoords and ((int)(values[i])-vcoords[-1]<=2 or (int)(values[i+4])-xcoords[-1]<=2 or int(values[i+8])-ucoords[-1]<=2)):
                                                continue

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
                x=(X*pitch)+pitch/2.0
                myline=TLine(x,-xmax,x,xmax)
                myline.SetLineColor(1)
                myline.SetLineWidth(1)
                linearray.append(myline)
                #myline.Draw("SAME")

                #draw plane for U Plane
                m1,c1=FindMandC(U,-theta,pitch)
	        myline1= TLine(-xmax,m1*(-xmax)+c1,xmax,m1*(xmax)+c1 )
                myline1.SetLineColor(1)
                myline1.SetLineWidth(1)
                linearray.append(RestrictLine(myline1, -xmax, -xmax, xmax, xmax,m1,c1))
                #myline1.Draw("SAME")

                #draw plane for V Plane
                m2,c2=FindMandC(V,theta,pitch)
	        myline2= TLine(-xmax,m2*(-xmax)+c2,xmax,m2*(xmax)+c2 )
                myline2.SetLineColor(1)
                myline2.SetLineWidth(1)
                linearray.append(RestrictLine(myline2, -xmax, -xmax, xmax, xmax,m2,c2))
                #myline2.Draw("SAME")

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
        

def RemoveAmbiguities(inHits):

        XAcceptedHits=[]
        XUAcceptedHits=[]
        XUVAcceptedHits=[]

        #check unique X
        for hitA in inHits:
                passed=True
                for hitB in inHits:
                        if hitA[2]==hitB[2] and GetMaxR(hitA)>GetMaxR(hitB):
                                passed=False
                if passed==True:
                        XAcceptedHits.append(hitA)

        #check unique U
        for hitA in XAcceptedHits:
                passed=True
                for hitB in XAcceptedHits:
                        if hitA[3]==hitB[3] and GetMaxR(hitA)>GetMaxR(hitB):
                                passed=False
                if passed==True:
                       XUAcceptedHits.append(hitA)           

        #check uniqueV
        for hitA in XUAcceptedHits:
                passed=True
                for hitB in XUAcceptedHits:
                        if hitA[4]==hitB[4] and GetMaxR(hitA)>GetMaxR(hitB):
                                passed=False
                if passed==True:
                       XUVAcceptedHits.append(hitA)  

        return XUVAcceptedHits
