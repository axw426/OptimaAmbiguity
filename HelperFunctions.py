import random
import math
from ROOT import gROOT, TCanvas, TF1, TF2, TH2F, TFile, TLine, TH1F, TLegend,TGraph
from array import array

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

def GetDirections(nProtons,beamSpread,cheat):
        mXmY=[]
        # m = tan(theta)
        for x in range(nProtons):
                mX=math.tan(random.uniform(-beamSpread/1000.0,beamSpread/1000.0))
                mY=math.tan(random.uniform(-beamSpread/1000.0,beamSpread/1000.0))
                if cheat==False:
                        mXmY.append([mX,mY])
                else:
                        mXmY.append([0,0])
        return mXmY

def ConvertXYToStrip(coord,pitch,rawAngle):

        return (int)(((coord[0]-50)*math.cos(math.radians(rawAngle)) - (coord[1]-50)*math.sin(math.radians(rawAngle))) /pitch)

def CheckStripHalf(coord,pitch,rawAngle):

        y= (int)(((coord[0]-50)*math.sin(math.radians(rawAngle)) + (coord[1]-50)*math.cos(math.radians(rawAngle))) /pitch)

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
                #print(angles[0],angles[1])
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

        return Strips

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

        #if nPlanes==2:

        
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

def FindRadius(i,j):

        return math.sqrt((i[0]-j[0])**2 + (i[1]-j[1])**2)

def ReconstructTracks2Planes(Hits,tolerance,pitch):

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
                                r=FindRadius(hitI,hitJ)
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

def GetChi2LinearFit(hits):

        n=len(hits)
        X,Y=array('d'),array('d')
        for i in range(n):
                X.append(hits[i][0])
                Y.append(hits[i][1])

        graph=TGraph(4,X,Y)

        graph.Fit("pol1","Q")

    
        return graph.GetFunction("pol1").GetChisquare()/(float)(graph.GetFunction("pol1").GetNDF())

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
                                                #chi2=GetChi2LinearFit([hitI,hitJ,hitK,hitL])
                                                chi2=SumRadialSeparation([hitI,hitJ,hitK,hitL])
                                                if chi2<bestChi2 and chi2>0.0001:
                                                        bestIPos=i
                                                        bestJPos=j
                                                        bestKPos=k
                                                        bestLPos=l
                                                        bestChi2=chi2

                #print("bestChi2=",bestChi2,[Hits[0][bestIPos],Hits[1][bestJPos],Hits[2][bestKPos],Hits[3][bestLPos]])
                if bestChi2<trackTolerance : #NEED TO DECIDE SENSIBLE VALUE
                        #Append hit positions to track list
                        RecoTracks.append([Hits[0][bestIPos],Hits[1][bestJPos],Hits[2][bestKPos],Hits[3][bestLPos]])
                        #delete entries from Hits
                        del Hits[0][bestIPos]
                        del Hits[1][bestJPos]
                        del Hits[2][bestKPos]
                        del Hits[3][bestLPos]
                else:
                       # print("Best Chi2=",bestChi2)
                        break

        return RecoTracks

def ReconstructTracks(Hits,trackTolerance,pitch):

        if len(Hits)==2:
                return ReconstructTracks2Planes(Hits,trackTolerance,pitch)   
        elif len(Hits)==4:
                return ReconstructTracks4Planes(Hits,trackTolerance,pitch)
        else:
                print("Can only cope with 2 or 4 tracker modules!")
              
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
