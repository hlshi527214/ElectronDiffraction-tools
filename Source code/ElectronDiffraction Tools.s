//Default values of some parameters
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:maxhkl",6)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:mind",8.5)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:max peak no",10)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Error d",0.3)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Error Angle",5)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Order",1)

number spotsize=5
if(!getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Spot size",spotsize))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Spot size",spotsize)
	}

number wavelength=0.1  //convert to XRD data based on Bragg's law
if(!getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Wavelength",wavelength))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Wavelength",wavelength)
	}
	
string unitstring="s"
if(!getpersistentStringnote("ElectronDiffraction Tools:Crystal Tools:Default:Unit",unitstring))
	{
		setpersistentStringnote("ElectronDiffraction Tools:Crystal Tools:Default:Unit",unitstring)
	}
	
	
number d1,d2,Angle, h1=0,k1=0,l1=0,h2=0,k2=0,l2=0,u=0,v=0,w=0


//-------------------------------Define some functions-------------------------------

// function: to match intensities of two spectrum
image MatchSpectralIntensities(image frontspectrum, image backspectrum, number roileft, number roiright, number iterations)
{
if(iterations<10) iterations=10

number xsize=frontspectrum.ImageGetDimensionSize(0)

if((roiright-roileft)==0)
{
roileft=0
roiright=xsize-1
}

// do the initial scaling on the basis of average intensities
number frontsum=sum(frontspectrum[0, roileft, 1, roiright])
number backsum=sum(backspectrum[0, roileft, 1, roiright])
number spectralratio=frontsum/backsum
frontspectrum=frontspectrum/spectralratio

number initialmisfit=sum(sqrt((backspectrum[0, roileft, 1, roiright]-frontspectrum[0, roileft, 1, roiright])**2))
number finalmultiplier=1
number misfit, i

// Iterate using a Montecarlo approach to further minimise the difference between the two spectra
for(i=0; i<iterations; i++)
{
number randnumber=gaussianrandom()/100
number testmultiplier=finalmultiplier+randnumber
image intensitytest=frontspectrum[0, roileft, 1, roiright]*testmultiplier

image misfitplot=sqrt((backspectrum[0, roileft, 1, roiright]-intensitytest)**2)
misfit=sum(misfitplot)

if(misfit<initialmisfit)
{
initialmisfit=misfit
finalmultiplier=testmultiplier
}
}

image scaledforespectrum=frontspectrum*finalmultiplier

return scaledforespectrum
}

// Function: to minimise the difference between the front spectrum and the back spectrum
image MatchSpectralChannels(image frontspectrum, image backspectrum, number roileft, number roiright)
{
number scanrange=roiright-roileft
image shiftedfore=imageclone(frontspectrum)
image misfitplot=realimage("",4, scanrange, 1)
number initialmisfit=1e99, i
number finalshift=0

// Scan -/+ half the ROI to look for a match
for(i=-scanrange/2; i<scanrange/2; i++)
{
shiftedfore=offset(frontspectrum, i, 0)
misfitplot=sqrt((backspectrum[0,roileft, 1, roiright]-shiftedfore[0,roileft, 1, roiright])**2)
number misfit=sum(misfitplot)

if(misfit<initialmisfit)
{
finalshift=i
initialmisfit=misfit
}
}

// Apply the shift and return the shifted spectrum
shiftedfore=offset(frontspectrum, finalshift, 0)

return shiftedfore
}


// Function: to perform Gaussian filter of a FT image.
image gaussfilter_FilterIt(image filter,number inixsize,number iniysize,number inix2,number iniy2, number gaura)
{
	image iniwarp:=RealImage("",8,inix2,iniy2)
	iniwarp=warp(filter,icol*inixsize/inix2,irow*iniysize/iniy2)
	compleximage ftimage:=RealFFT(iniwarp)
	ftimage*=exp(-(((icol-inix2/2)/inixsize)**2+((irow-iniy2/2)/iniysize)**2)*gaura**2)

	if(gaura>=0)
	{
		iniwarp=RealIFFT(ftimage)
	}
	else
	{
		iniwarp=iniwarp-RealIFFT(ftimage)
	}

	image final:=RealImage("Filtered image",8,inixsize,iniysize)
	final=warp(iniwarp,icol/inixsize*inix2,irow/iniysize*iniy2)

	return final
}

//Function: to perform Gaussian filter of an image.
image GaussFilter(image ini, number gaura)
{
image res 
number inix2, iniy2
number inixsize,iniysize

Getsize(ini,inixsize,iniysize)

for(inix2=1;inix2<inixsize;inix2*=2)
{}

for(iniy2=1;iniy2<iniysize;iniy2*=2)
{}

// story begins
if(ini.ImageIsDataTypeRGB())
	{
		res:=RGBImage("Filtered image",4,inixsize,iniysize)
		res= rgb(gaussfilter_FilterIt(ini.red(),inixsize,iniysize,inix2,iniy2,gaura),gaussfilter_FilterIt(ini.green(),inixsize,iniysize,inix2,iniy2,gaura),gaussfilter_FilterIt(ini.blue(),inixsize,iniysize,inix2,iniy2,gaura))
	}
	else
	if(ini.IsPackedComplexImage()||ini.ImageIsDataTypeComplex())
	{
		res:=ComplexImage("Filtered image",16,inixsize,iniysize)
		res=complex(gaussfilter_FilterIt(ini.real(),inixsize,iniysize,inix2,iniy2,gaura),gaussfilter_FilterIt(ini.imaginary(),inixsize,iniysize,inix2,iniy2,gaura))
	} 
	else
	{
		res:=RealImage("Filtered image",8,inixsize,iniysize)
		res= gaussfilter_FilterIt(ini,inixsize,iniysize,inix2,iniy2,gaura)

     	SetSurveyTechnique(res, 1 ) 
	}

return res
}

//Function: to calculate FWHM of a profile
number fwhm(image img, number x,number &leftpart,number &rightpart)
{
number width,xsize,ysize,leftwidth,rightwidth
img.GetSize(xsize,ysize)

number height=img.GetPixel(x,0)
for(number i=x;i>0;i--)
{
number intensity=img.GetPixel(i,0)

if(intensity<0.8*height)
{
leftwidth=abs(x-i)
leftpart=i

break
}
}

for(number i=x;i<xsize;i++)
{
number intensity=img.GetPixel(i,0)

if(intensity<0.8*height)
{
rightwidth=abs(x-i)
rightpart=i

break
}
}

width=leftwidth+rightwidth
return width
}

//Function: to create 2D Gaussian image
image Create2DGaussian(number xsize, number ysize, number centrex, number centrey, number intensity, number sigmax, number sigmay)
	{
		image gaussian=realimage("", 4,xsize, ysize)
		if(centrex>xsize || centrex<0 || centrey>ysize || centrey<0) return gaussian
		gaussian=intensity*exp(-(icol-centrex)**2/(2*sigmax**2))*exp(-(irow-centrey)**2/(2*sigmay**2))

		return gaussian
	}


//Function: to compute chi squred x2 
number computechisqrd(image original, image fit, number centrex, number centrey, number centralpcnt)
	{
		number xsize, ysize
		getsize(original, xsize, ysize)
		
		// Compare only subarea of the image - centred on the target centre.
		if(centralpcnt<25) centralpcnt=25
		if(centralpcnt>75) centralpcnt=75
		
		// Keep the size of the central region sensible >=10 pixels		
		number xdim=xsize*(centralpcnt/100)
		if(xdim<10) xdim=10
		
		number ydim=ysize*(centralpcnt/100)
		if(ydim<10) ydim=10
		
		
		// Define the search area centred on the target centre
		number t=round(centrey-(ydim/2))
		number l=round(centrex-(xdim/2))
		number b=round(centrey+(ydim/2))
		number r=round(centrex+(xdim/2))
		
		if(t<0) t=0
		if(l<0) l=0
		if(b>ysize-1) b=ysize-1
		if(r>xsize-1) r=xsize-1
		

		// measure the difference between the original and the fit (gaussian) images within their subareal regions
		// and compute chisqrd		
		image difference=(original[t,l,b,r]-fit[t,l,b,r])**2/(original[t,l,b,r]+0.000001)
		number chisqrd=sum(difference)
		return chisqrd
	}


//Function: to find the parameter sigma for Gaussian fitting by the Monte Carlo method.
number FindSigmaMonteCarlo(image front, number centrex, number centrey, number subareapcnt, number intensity, number initialsigma, number iterations, number &minchisqrd, number &sigma)
	{
		// Source the image size
		number xsize, ysize, chisqrd
		getsize(front, xsize, ysize)
		number i

		// random walk - vary the initial sigma with a random gaussian function and compute the fit for the resulting sigma		
		for(i=0; i<iterations; i++)
			{
				number sigmatest=initialsigma+gaussianrandom()
				if(sigmatest<1) sigmatest=1
				image testgaussian=Create2DGaussian(xsize, ysize, centrex, centrey, intensity, sigmatest, sigmatest)			
				chisqrd=computechisqrd(front, testgaussian, centrex, centrey, subareapcnt)

				if(chisqrd<minchisqrd) 
					{
						minchisqrd=chisqrd
						sigma=sigmatest
						initialsigma=sigmatest
					}	
			}
		return minchisqrd
	}


// Function: to find the ceneter of a SAED pattern by the Monte Carlo method
number FindCentreMonteCarlo(image front, number sigma, number intensity, number iterations, number centrex, number centrey, number subareapcnt, number &minchisqrd, number &fitcentrex, number &fitcentrey)
	{
		number xsize, ysize, chisqrd, newx, newy
		getsize(front, xsize, ysize)
		image testgaussian=front*0
	
		// If no centre information is provided  estimate the centre as the geometric centre
		if(centrex==0) centrex=xsize/2
		if(centrey==0) centrey=ysize/2
		number i, randcentrex, randcentrey
		fitcentrex=centrex
		fitcentrey=centrey
	
		// Do a random walk looking for an improved fit
		for(i=0; i<iterations; i++)
			{
				// vary the centres with a random walk
				randcentrex=gaussianrandom()
				randcentrey=gaussianrandom()
				newx=randcentrex+fitcentrex
				newy=randcentrey+fitcentrey
				
				
				// Ignore any values which take the centre outside the bounds of the image
				if(newx<0) randcentrex=0
				if(newx>xsize-1) randcentrex=0
				if(newy<0) randcentrey=0
				if(newy>ysize-1) randcentrey=0
			
				// test the fit with the new centre
				testgaussian=Create2DGaussian(xsize, ysize, newx,newy, intensity, sigma, sigma)
				chisqrd=computechisqrd(front, testgaussian, newx, newy, subareapcnt)

				if(chisqrd<minchisqrd) 
					{
						minchisqrd=chisqrd
						fitcentrex=newx
						fitcentrey=newy
					}	
			}
		
		return minchisqrd
	}
	
	
//Function: to calculate FWHM of a profile
number FWHMcalc(image img, number maxx)
{
number xpos, ypos, maxval, i,xsize,ysize,t,l,b,r
number halfmax, tenthmax, lefthalfstartpos,righthalfstartpos, tenthstartpos,left,right
img.GetSize(xsize,ysize)

img.GetSelection(t,l,b,r)
maxval=img.GetPixel(maxx,0)
halfmax=maxval/1.7

for(i=maxx; i>0; i--)  //left part
	{
		number thisval=getpixel(img, i,0)
		
		if(thisval<halfmax)
			{
				left=i
				img[0,0,1,left]=thisval
				break
			}
	}	

for(i=maxx; i<xsize; i++)  //right part
	{
		number thisval=getpixel(img, i,0)
		
		if(thisval<halfmax)
			{
				right=i			
				img[0,right,1,xsize]=thisval
				break
			}
	}	

img.setselection(0,left,1,right)

Return right-left
}

//Function: to find peaks of a profile by the top hat method
image tophatfilter(image inputimage, number defaultbrimwidth, number defaulthatwidth, number tophatthreshold, number PeakorTrough)
	{ 
		number i,  brimmean, hatmean, tophatresult, x,no=1
		number positivecounter, negativecounter, counter,nomax,nomin,leftpart,rightpart
		image lowerbrim, upperbrim, hat
		number previoussign=0	

		number maxval=max(inputimage)
		inputimage=inputimage/maxval
		number xsize, ysize
		getsize(inputimage, xsize, ysize)
		
		imagedisplay imgdisp=imagegetimagedisplay(inputimage,0)

		//begin to find peaks
		for (i=1; i<(xsize-(defaulthatwidth+(2*defaultbrimwidth))); i++) 
			{
				lowerbrim=inputimage[0,i,1,i+defaultbrimwidth]
				upperbrim=inputimage[0,i+defaultbrimwidth+defaulthatwidth,1,i+defaultbrimwidth+defaulthatwidth+defaultbrimwidth]
				hat=inputimage[0,i+defaultbrimwidth,1,i+defaultbrimwidth+defaulthatwidth]

				brimmean=(mean(lowerbrim)+mean(upperbrim))/2
				hatmean=mean(hat)
				tophatresult=hatmean-brimmean
				
				number temptest=sqrt(tophatresult*tophatresult) 
		if(temptest<tophatthreshold || (sgn(tophatresult)!=previoussign && temptest>tophatthreshold))
			{
				if (positivecounter>0) // there is a positive peak
					{
						number peakmax=(positivecounter/counter)+defaultbrimwidth+(defaulthatwidth/2)
						
						fwhm(inputimage, peakmax,leftpart,rightpart)
						number peak1=(leftpart+rightpart)/2
						
						if(peak1>peakmax-5&&peak1<peakmax-5)peakmax=peak1

						if (PeakorTrough==0 || PeakorTrough==2) //mark it
							{ 
							roi theroi=createroi()
							theroi.roisetrange(peakmax, peakmax)
							roisetcolor(theroi,1,0,0)
							theroi.roisetvolatile(0)
							imgdisp.imagedisplayaddroi(theroi)						
							
							}
						positivecounter=0 //clear all peak stores
						counter=0
						nomax=nomax+1
					}
					
					
				if (negativecounter>0) // there's a negative peak (trough)
					{
						number peakmin=(negativecounter/counter)+defaultbrimwidth+(defaulthatwidth/2)
												
						if (PeakorTrough==1 || PeakorTrough==2) 
							{ 
								roi theroi=createroi()
								theroi.roisetrange(peakmin, peakmin)

								imgdisp.imagedisplayaddroi(theroi)
								roisetcolor(theroi,0,0,1) // troughs are blue
								theroi.roisetvolatile(0)								
							}

						counter=0
						nomin=nomin+1
						negativecounter=0
					}					
					
			}

				if(tophatresult>tophatthreshold) // indicates a peak detected
					{
						positivecounter=positivecounter+i
						counter=counter+1
						negativecounter=0
					}			

				if((-1*tophatresult)>tophatthreshold)	//indicates a trough detected
					{
						negativecounter=negativecounter+i
						positivecounter=0
						counter=counter+1
					}					

		previoussign=sgn(tophatresult) // store the sign of this tophatresult for comparison with the next	

			}

inputimage=inputimage*maxval

return(inputimage)			
}

//Function: to create damped edges of an image
image createdampededge(number edgewidth, number ximagewidth, number yimagewidth)
	{
		// The left hand edge		
		image left=realimage("",4,ximagewidth, yimagewidth)
		left=((sin(((icol-(edgewidth/2))/edgewidth)*pi())+1)/2)
		left=tert(icol>edgewidth, 1, left)

		// Form the right hand edge by flipping the left hand edge image		
		image right=imageclone(left)
		right=left[ximagewidth-icol, irow]

		// Calculate the top edge
		image top=imageclone(left)
		top=((sin(((irow-(edgewidth/2))/edgewidth)*pi())+1)/2)
		top=tert(irow>edgewidth, 1, top)

		// Form the bottom edge by flippig the top image
		image bottom=imageclone(left)
		bottom=top[icol, yimagewidth-irow]

		// Multiple the four images together , then threshold and scale them so that non-edge areas are 1
		// while edge areas scale down from 1 to zero at the edge
		image productimage=left*right*top*bottom
		number threshval=getpixel(productimage, edgewidth, yimagewidth/2)
		image edgeimage=tert(productimage>threshval, threshval, productimage)
		edgeimage=edgeimage/threshval
		
		return edgeimage
	}
	
	
//Function: to perform a quick Gaussian fitting
void CoarseFit(image image1D, number weighting, number &peakintensity, number &peakchannel, number &sigma, number &fwhm)
	{
		number xsize, ysize
		getsize(image1D, xsize, ysize)
		image1D.SetSelection(0,0.3*xsize,1,0.7*xsize)
		
		number t,d,l, r  //d-down
		image1D.getselection(t,l,d,r)

		number origin, dispersion
		number xmax, ymax, i, channel, intensity
		number c1, c2, c3, c4, c5, yx2, yx1, yx0
		number delta, a, b, c, alpha,fitintensity
		string imgname

		converttofloat(image1D)

		number peakmax=max(image1D, xmax, ymax)
		image1D=image1D/peakmax

		for (i=l; i<r; i++)
			{
				intensity=getpixel(image1D, i, 0)

				c1=c1+(intensity**weighting*i**4)
				c2=c2+(intensity**weighting*i**3)
				c3=c3+(intensity**weighting*i**2)
				c4=c4+(intensity**weighting*i)
				c5=c5+intensity**weighting

				yx2=yx2+(intensity**weighting*log(intensity)*i**2)
				yx1=yx1+(intensity**weighting*log(intensity)*i)
				yx0=yx0+(intensity**weighting*log(intensity))
			}

		delta=-c3**3+2*c2*c3*c4-c1*c4**2-c2**2*c5+c1*c3*c5
		a=(((-c3**2+c2*c4)*yx0)+((c3*c4-c2*c5)*yx1)+((-c4**2+c3*c5)*yx2))/delta
		b=(((c2*c3-c1*c4)*yx0)+((-c3**2+c1*c5)*yx1)+((c3*c4-c2*c5)*yx2))/delta
		c=(((-c2**2+c1*c3)*yx0)+((c2*c3-c1*c4)*yx1)+((-c3**2+c2*c4)*yx2))/delta

		sigma=sqrt((-1/(2*a)))
		peakchannel=b*sigma*sigma
		peakintensity=peakmax*exp(c+peakchannel**2/(2*sigma**2))
		fwhm=2*sigma*sqrt(2*log(2))
		
		image1D.ClearSelection()
}

//Function: to perform Gaussian fitting of an image.
image FitGaussian(image image1D, number weighting, number &peakintensity, number &peakchannel, number &sigma, number &fwhm, number &Rp)
{
number xsize,ysize
image1D.GetSize(xsize,ysize)

number t,d,l, r
image1D.getselection(t,l,d,r)

number origin, dispersion
number xmax, ymax, i, channel, intensity
number c1, c2, c3, c4, c5, yx2, yx1, yx0
number delta, a, b, c, alpha, fitintensity
string imgname

converttofloat(image1D)

number peakmax=max(image1D, xmax, ymax)
image1D=image1D/peakmax

for (i=l; i<r; i++)
{
intensity=getpixel(image1D, i, 0)

c1=c1+(intensity**weighting*i**4)
c2=c2+(intensity**weighting*i**3)
c3=c3+(intensity**weighting*i**2)
c4=c4+(intensity**weighting*i)
c5=c5+intensity**weighting

yx2=yx2+(intensity**weighting*log(intensity)*i**2)
yx1=yx1+(intensity**weighting*log(intensity)*i)
yx0=yx0+(intensity**weighting*log(intensity))
}

delta=-c3**3+2*c2*c3*c4-c1*c4**2-c2**2*c5+c1*c3*c5
a=(((-c3**2+c2*c4)*yx0)+((c3*c4-c2*c5)*yx1)+((-c4**2+c3*c5)*yx2))/delta
b=(((c2*c3-c1*c4)*yx0)+((-c3**2+c1*c5)*yx1)+((c3*c4-c2*c5)*yx2))/delta
c=(((-c2**2+c1*c3)*yx0)+((c2*c3-c1*c4)*yx1)+((-c3**2+c2*c4)*yx2))/delta

sigma=sqrt((-1/(2*a)))
peakchannel=b*sigma*sigma
peakintensity=peakmax*exp(c+peakchannel**2/(2*sigma**2))
fwhm=2*sigma*sqrt(2*log(2))

image1D = image1D*peakmax

image peakfit=image1D*0
image Rpimage=image1D*0   //calculate Rp

for (i=l; i<r; i++)
{
fitintensity=peakintensity*exp(-((i-peakchannel)**2)/(2*sigma**2))
setpixel(peakfit, i, 0,fitintensity)
}

Rp=sum((image1D[t,l,d,r]-peakfit[t,l,d,r])**2/(image1D[t,l,d,r]+0.000001))

return peakfit
}


// Function: to create tags of the daily used crystals
void createcrystaltags(string crystalname, number lattice, number a, number b, number c, number alpha, number beta, number gamma)
{
string lattices
if(!getpersistentstringnote("ElectronDiffraction Tools:Crystal Tools:Crystals",lattices))

	setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":Lattice",lattice)
	setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":a",a)
	setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":b",b)
	setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":c",c)
	setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":alpha",alpha)
	setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":beta",beta)
	setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":gamma",gamma)
}

//The seven crystal systems: cubic=1、tetragonal=2、orthorhombic=3、hexagonal=4、rhombohedral=5,monoclinic=6、triclinic=7
createcrystaltags("Au-Gold", 1, 4.078, 4.078, 4.078, 90, 90, 90) // Gold
createcrystaltags("Si-Silicon", 1, 5.4282, 5.4282, 5.4282, 90, 90, 90)// Silicon
createcrystaltags("TiO2-anatase", 2, 3.7852, 3.7852, 9.5139, 90, 90, 90)// TiO2 - anatase
createcrystaltags("TiO2-brookite", 3, 5.4558, 9.1819, 5.1429, 90, 90, 90)// TiO2 - brookite
createcrystaltags("TiO2-rutile", 2, 4.5933, 4.5933, 2.9592, 90, 90, 90)// TiO2 - rutile
createcrystaltags("R-CaCO3", 4, 4.978, 4.978, 16.978, 90, 90, 120)// CaCO3 - R
createcrystaltags("R-CaO", 1, 4.8049, 4.8049, 4.8049, 90, 90, 90)// CaO


Object CrystUIFrame

//Function: to create a radio dialog to set crystal system
TagGroup MakeCrystTypeRadios()
	{
		Taggroup crystradio_items
		TagGroup crysttyperadio = DLGCreateRadioList(crystradio_items, 0).dlgchangedmethod("crystalchanged")
		crystradio_items.DLGAddRadioItem("C",1).DLGSide("Left")
		crystradio_items.DLGAddRadioItem("T",2).DLGSide("Left")
		crystradio_items.DLGAddRadioItem("O",3).DLGSide("Left")
		crystradio_items.DLGAddRadioItem("H",4).DLGSide("Left")
		crystradio_items.DLGAddRadioItem("R",5).DLGSide("Left")
		crystradio_items.DLGAddRadioItem("M",6).DLGSide("Left")
		crystradio_items.DLGAddRadioItem("A",7).DLGSide("Left")
		crysttyperadio.DLGIdentifier("Radio")

		return crysttyperadio	
	}
	
// Function: to create a pull down dialog
TagGroup MakeElementPullDown()
{
	TagGroup popup_items
	TagGroup popup = DLGCreatePopup(popup_items, 0, "elementselect")

	string crystalname
	
	// read tags
	taggroup ptg=getpersistenttaggroup()
	taggroup crystaltags
	ptg.TaggroupGetTagAsTagGroup("ElectronDiffraction Tools:Crystal Tools:Crystals", crystaltags)

	number count=crystaltags.taggroupcounttags()
	number i, lattice, a, b, c, alpha, beta, gamma

	for(i=0; i<count; i++)
	{
		String label = crystaltags.TagGroupGetTagLabel( i ) 
		
		if(i==0)popup_items.DLGAddPopupItemEntry(label,1)
		
		else
		{
		popup_items.DLGAddPopupItemEntry(label)
		}
		
	}

	popup.DLGIdentifier("elementselect")

	return popup
}

// --------------------------Class: electron diffraction tools-------------------------
Class ElectronDiffracion_Tools: UIFrame
{
//---------------------Define some function in this class-----------------
//Function: to convert rhombhedral to hexagonal
void rhombtohexconversion(object self, number rhombalpha, number arhomb, number &ahex, number &chex)
	{
		ahex=2*arhomb*sin(rhombalpha/2)
		chex=sqrt(3)*sqrt(arhomb**2+(2*arhomb**2*cos(rhombalpha)))
	}

// Function: to convert hexagonal to rhombhedral
void hextorhombconversion(object self, number ahex, number chex, number &rhombalpha, number &arhomb)
	{
		arhomb=(1/3)*sqrt((3*ahex**2)+chex**2)
		number rhombangle=3/(2*sqrt(3+(chex/ahex)**2))

		rhombangle=asin(rhombangle)*2
		rhombalpha=rhombangle*(180/pi())
	}

	
// Function: to change crystal tabs.
	void crysttabchanged(object self, taggroup tg)
	{
		number tabvalue
		tg.dlggetvalue(tabvalue)
		
		if(tabvalue==0) 	self.setelementisenabled("Add_pattern",1) // the Tables button is disabled for all but dspacings
		else self.setelementisenabled("Add_pattern",0)
		
		if(tabvalue==3) 
			{
				//result("\nZone Axis Calculation: (h2,k2,l2) should be anticlockwise relative to (h1,k1,l1)\n")
			}
	}


	// Function: response when changing lattice parameters a, b, and c.
	void latticechanged(object self, taggroup tg)
	{
			number crystalvalue, latticevalue
			latticevalue=tg.dlggetvalue() 

			self.dlggetvalue("crystalradios", crystalvalue)

			if(crystalvalue==2|| crystalvalue==4) // set a=b for tetragonal and hexagonal
				{
					self.dlgvalue("afield", latticevalue)
					self.dlgvalue("bfield", latticevalue)
				}

			if(crystalvalue==1 || crystalvalue==5) //set a=b=c for cubic and rhombhedral
				{
					self.dlgvalue("afield", latticevalue)
					self.dlgvalue("bfield", latticevalue)
					self.dlgvalue("cfield", latticevalue)
				}
				
//Clear all data
	TagGroup TgList1= self.LookUpElement("hklList1")  
	TagGroup TgList2= self.LookUpElement("hklList2")  
	TagGroup uvwList= self.LookUpElement("uvwList") 
	uvwList.DLGRemoveElement(0)  
	
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags))
	{
	number tagno=dspacingtags.taggroupcounttags()

	for(number i=0;i<200;i++)
	{
	TgList1.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	}
	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}	
	
	}


	// Function: response when changing lattice parameters α, β and γ.
	void alphachanged(object self, taggroup tg)
	{
		number alphaval, crystalvalue
		alphaval=tg.dlggetvalue()
		self.dlggetvalue("crystalradios", crystalvalue)

		if(crystalvalue==1 || crystalvalue==2 || crystalvalue==3) // set α= β= γ=90 for cubic, tetragonal and orthogonal 
		{
			self.dlgvalue("alphafield",90)
			self.dlgvalue("betafield",90)
			self.dlgvalue("gammafield",90)
		}

		if(crystalvalue==4) // set α= β=90,  γ=120 for hexagonal 
		{
			self.dlgvalue("alphafield",90)
			self.dlgvalue("betafield",90)
			self.dlgvalue("gammafield",120)
		}

		if(crystalvalue==5) //set α= β= γ for rhombhedral
		{
			if(alphaval!=90 && alphaval<120)
				{
					self.dlgvalue("betafield",alphaval)
					self.dlgvalue("gammafield",alphaval)
				}
			else
				{
					beep()
					self.dlgvalue("alphafield",0)
					self.dlgvalue("betafield",0)
					self.dlgvalue("gammafield",0)
				}
		}

	if(crystalvalue==6) //set α= γ=90 for monoclinic
		{
			self.dlgvalue("alphafield",90)
			self.dlgvalue("gammafield",90)
		}

	// Clear all data			
	TagGroup TgList1= self.LookUpElement("hklList1")
	TagGroup TgList2= self.LookUpElement("hklList2")
	TagGroup uvwList= self.LookUpElement("uvwList")
	uvwList.DLGRemoveElement(0) 
	
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags))
	{
	number tagno=dspacingtags.taggroupcounttags()

	for(number i=0;i<200;i++)
	{
	TgList1.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	}
	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}		
	}

// Function: to response when d1 is changed
	void d1changed(object self, taggroup tg)
	{
	image img:=getfrontimage()	
	Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d1",d1)
	self.dlgvalue("d1field",d1)
	self.validateview()
	}	
	
// Function: to response when d2 is changed
	void d2changed(object self, taggroup tg)
	{
	image img:=getfrontimage()	
	Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d2",d2)	
	self.dlgvalue("d2field",d2)
	self.validateview()
	}	
	
// Function: to response when Angle is changed
	void Anglechanged(object self, taggroup tg)
	{
	image img:=getfrontimage()	
	Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",Angle)	
	self.dlgvalue("Anglefield",Angle)
	self.validateview()
	}		


// Function: to clear all the user-defined data
void clearfields(object self)
		{
		self.dlgvalue("crystalnamefield","")

		self.lookupElement("dspacingfield").dlgtitle("  "+ format(0, "%7.4f"+" A")) 
		self.dlgvalue("dhfield",0)
		self.dlgvalue("dkfield",0)
		self.dlgvalue("dlfield",0)
		
		self.dlgvalue("angh1field",0)
		self.dlgvalue("angk1field",0)
		self.dlgvalue("angl1field",0)

		self.dlgvalue("angh2field",0)
		self.dlgvalue("angk2field",0)
		self.dlgvalue("angl2field",0)
		self.lookupElement("angleresultfield").dlgtitle(""+ format(0, "%6.2f"+"deg.")) 
		
		self.dlgvalue("diru1field",0)
		self.dlgvalue("dirv1field",0)
		self.dlgvalue("dirw1field",0)

		self.dlgvalue("diru2field",0)
		self.dlgvalue("dirv2field",0)
		self.dlgvalue("dirw2field",0)
		self.lookupElement("dirresultfield").dlgtitle(""+ format(0, "%6.2f"+"deg.")) 
		
		self.dlgvalue("zoneh1field",0)
		self.dlgvalue("zonek1field",0)
		self.dlgvalue("zonel1field",0)

		self.dlgvalue("zoneh2field",0)
		self.dlgvalue("zonek2field",0)
		self.dlgvalue("zonel2field",0)
		self.lookupElement("zoneresultfield").dlgtitle( "  0 , 0 , 0 ")
		self.validateview()
	}


//Function: to normalize indices
	number normaliseindices(object self, number &x, number &y, number &z)
	{
		number i, origx, origy, origz, divisor

		origx=x
		origy=y
		origz=z

		for(i=12; i>1;i--)
			{
				if (mod(x,i)==0 && mod(y,i)==0 && mod(z,i)==0)
				{
					x=x/i
					y=y/i
					z=z/i
				}
			}

		number scaledx, scaledy, scaledz
		if(origx!=0) scaledx=origx/x
		else scaledx=1

		if(origy!=0) scaledy=origy/y
		else scaledy=1

		if(origz!=0) scaledz=origz/z
		else scaledz=1

		number temp1=max(scaledx, scaledy)
		divisor=max(temp1, scaledz)

		return divisor
	}

//Function: to delete the saved crystal data
	void crystclearbutton(object self)
	{
	// ALT+Clear: delete the selected crystal
	if(optiondown()) 
		{
			showalert("Select a Crystal file to load.",2)

			taggroup loadedtags=newtaggroup()
			string path, extension, filelabel, version, filename
			
			if(!opendialog(path)) return

			extension=pathextractextension(path,0)
			filename=pathextractfilename(path,0)
			if(extension!="txt") 
				{
					showalert("You need a *.txt file",1)
					return
				}
	

			if(!twobuttondialog("Do you want to replace any existing crystal data with : "+filename,"Add","Cancel"))
				{
					showalert("Cancelled by user.",2)
					return
				}

			//delete all tags
			taggroup ptags=getpersistenttaggroup()
			taggroup EDTags
			ptags.taggroupgettagastaggroup("ElectronDiffraction tools", EDTags)
			EDTags.taggroupdeletetagwithlabel("Crystal Tools:Crystals")

			TagGroupAddLabeledTagGroup( EDTags, "Crystal Tools:Crystals", loadedtags )
			showalert("         Crystals are successfully loaded!\n\nData is available when ElectronDiffraction tools is next launched.",1)
			return
		}

//clear all the user-defined data
		self.dlgvalue("dhfield",0)
		self.dlgvalue("dkfield",0)
		self.dlgvalue("dlfield",0)
  
		self.dlgvalue("afield",0)
		self.dlgvalue("bfield",0)
		self.dlgvalue("cfield",0)

		self.dlgvalue("alphafield",90)
		self.dlgvalue("betafield",90)
		self.dlgvalue("gammafield",90)

		self.dlgvalue("crystalradios",0)

		crystUIFrame.clearfields()
	
	TagGroup TgList1= self.LookUpElement("hklList1") 
	TagGroup TgList2= self.LookUpElement("hklList2") 
	TagGroup uvwList= self.LookUpElement("uvwList")
	uvwList.DLGRemoveElement(0)  
	
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags))
	{
	number tagno=dspacingtags.taggroupcounttags()

	for(number i=0;i<200;i++)
	{
	TgList1.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	}
	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}	
	
	}


//Function: to response when the crystal (pull-down list) is changed
	void crystalchanged(object self, taggroup tg)
	{
		crystUIFrame.clearfields()
		
		//to set lattice parameters for different type of crystal system
		number crystalvalue
		crystalvalue=tg.dlggetvalue() 
				if(crystalvalue==1 || crystalvalue==2 || crystalvalue==3 ) // Cubic, tetragonal, orthogonal
					{
						self.dlgvalue("alphafield", 90)
						self.dlgvalue("betafield", 90)
						self.dlgvalue("gammafield", 90)
					}

				if(crystalvalue==4) // hexagonal
					{
						self.dlgvalue("alphafield", 90)
						self.dlgvalue("betafield", 90)
						self.dlgvalue("gammafield", 120)
					}
					
				if(crystalvalue==5 || crystalvalue==7) // rhombhedral
					{
						self.dlgvalue("alphafield", 0)
						self.dlgvalue("betafield", 0)
						self.dlgvalue("gammafield", 0)
					}

				if(crystalvalue==6) // monoclinic
					{
						self.dlgvalue("alphafield", 90)
						self.dlgvalue("betafield", 0)
						self.dlgvalue("gammafield", 90)
					}

				// if the selected list is none, set default values
				number pulldownvalue
				taggroup pulldown=self.lookupelement("elementpulldown")
				self.dlggetvalue("elementpulldown", pulldownvalue)
				
				string pdname
				dlggetnthlabel(pulldown, pulldownvalue-1, pdname)

				if (pdname=="")
					{
						self.dlgvalue("afield",0)
						self.dlgvalue("bfield",0)
						self.dlgvalue("cfield",0)
					}
			
				self.validateview()

//clear all the list				
	TagGroup TgList1= self.LookUpElement("hklList1")
	TagGroup TgList2= self.LookUpElement("hklList2")
	TagGroup uvwList= self.LookUpElement("uvwList")
	uvwList.DLGRemoveElement(0)  
	
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags))
	{
	number tagno=dspacingtags.taggroupcounttags()

	for(number i=0;i<200;i++)
	{
	TgList1.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	}
	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}					
				
	}


	//*********************Calc. button in d-hkl pannel*******************
	void calculatedspacing(object self)
	{
		number crystalvalue, alpha, beta, gamma, a, b, c, h, k, l, dspacing
		string crystalname

		//read lattice parameters
		self.dlggetvalue("crystalradios",crystalvalue)
		self.dlggetvalue("alphafield",alpha)
		self.dlggetvalue("betafield",beta)
		self.dlggetvalue("gammafield",gamma)
		self.dlggetvalue("afield",a)
		self.dlggetvalue("bfield",b)
		self.dlggetvalue("cfield",c)

		self.dlggetvalue("dhfield",h)
		self.dlggetvalue("dkfield",k)
		self.dlggetvalue("dlfield",l)

		alpha=alpha/(180/pi())
		beta=beta/(180/pi())
		gamma=gamma/(180/pi())

		if(a<=0 || b<=0 || c<=0 || alpha<=0 || beta<=0 || gamma<=0 || (h**2+k**2+l**2)==0)
			{
				beep()
				exit(0)
			}

		// calculate d-spacing for a given hkl.
		if(crystalvalue==1)  
			{
				dspacing=sqrt(a**2/(h**2+k**2+l**2))
				crystalname="Cubic"
			}

		if(crystalvalue==2)
			{
				number temp=((1/a**2)*(h**2+k**2))+((1/c**2)*l**2)
				dspacing=sqrt(1/temp)
				crystalname="Tetragonal"
			}

		if(crystalvalue==3)
			{
				number temp=(((1/a**2)*h**2)+((1/b**2)*k**2)+((1/c**2)*l**2))
				dspacing=sqrt(1/temp)
				crystalname="Orthorhombic"
			}
	
		if(crystalvalue==4) 
			{
				number temp=((4/(3*a**2))*(h**2+(h*k)+k**2)+((1/c**2)*l**2))
				dspacing=sqrt(1/temp)
				crystalname="Hexagonal"
			}

		if(crystalvalue==5) 
			{
				number temp1=((h**2+k**2+l**2)*sin(alpha)**2)+(2*((h*k)+(k*l)+(h*l))*(cos(alpha)**2-cos(alpha)))
				number temp2=a**2*(1-(3*cos(alpha)**2)+(2*cos(alpha)**3))
				number temp3=temp1/temp2
				dspacing=sqrt(1/temp3)
				crystalname="Rhombohedral"
			}	

		if(crystalvalue==6)
			{
				number temp=((1/a**2)*(h**2/sin(beta)**2))+((1/b**2)*k**2)+((1/c**2)*(l**2/sin(beta)**2))-((2*h*l*cos(beta))/(a*c*sin(beta)**2))
				dspacing=sqrt(1/temp)
				crystalname="Monoclinic"
			}


		if(crystalvalue==7)
			{
				number V2=a**2*b**2*c**2*(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+(2*cos(alpha)*cos(beta)*cos(gamma)))
				number s11=b**2*c**2*sin(alpha)**2
				number s22=a**2*c**2*sin(beta)**2
				number s33=a**2*b**2*sin(gamma)**2

				number s12=a*b*c**2*(cos(alpha)*cos(beta)-cos(gamma))
				number s23=a**2*b*c*(cos(beta)*cos(gamma)-cos(alpha))
				number s31=a*b**2*c*(cos(gamma)*cos(alpha)-cos(beta))

				number temp=(1/v2)*((s11*h**2)+(s22*k**2)+(s33*l**2)+(2*s12*h*k)+(2*s23*k*l)+(2*s31*l*h))
				dspacing=sqrt(1/temp)
				crystalname="Triclinic"
			}

		self.lookupElement("dspacingfield").dlgtitle("  "+ format(dspacing, "%7.4f"+" A" )) 
		result("\n("+h+","+k+","+l+"), "+"d="+format(dspacing, "%2.4f")+" A\n")

		self.validateview()
	}


//*****************************Function: Planes pannel*******************************************
	void calculateplaneangles(object self)
	{
		number crystalvalue, alpha, beta, gamma, a, b, c, h1, k1, l1,h2,k2,l2, u,v,w,angle
		number p1, p2, q1, q2, s1, s2, r1, r2,d1,d2

		string crystalname, material

		// read lattice parameters and planes
		self.dlggetvalue("crystalradios",crystalvalue)
		self.dlggetvalue("alphafield",alpha)
		self.dlggetvalue("betafield",beta)
		self.dlggetvalue("gammafield",gamma)
		self.dlggetvalue("afield",a)
		self.dlggetvalue("bfield",b)
		self.dlggetvalue("cfield",c)

		self.dlggetvalue("angh1field",h1)
		self.dlggetvalue("angk1field",k1)
		self.dlggetvalue("angl1field",l1)

		self.dlggetvalue("angh2field",h2)
		self.dlggetvalue("angk2field",k2)
		self.dlggetvalue("angl2field",l2)

		alpha=alpha/(180/pi())
		beta=beta/(180/pi())
		gamma=gamma/(180/pi())

		// begin to calculate d-spacings and interAngle	
		if(crystalvalue==1) //cubic
			{
				number costheta, temp1, temp2
				temp1=(h1*h2)+(k1*k2)+(l1*l2)
				temp2=sqrt((h1**2+k1**2+l1**2)*(h2**2+k2**2+l2**2))

				costheta=temp1/temp2
				angle=acos(costheta)
				angle=angle*(180/pi())
				
				d1=sqrt(a**2/(h1**2+k1**2+l1**2))
				d2=sqrt(a**2/(h2**2+k2**2+l2**2))
				
				crystalname="Cubic"
			}

		if(crystalvalue==2) // Tetragonal
			{
				number costheta, temp1, temp2, temp3, temp4
				temp1=((1/a**2)*((h1*h2)+(k1*k2)))+((1/c**2)*l1*l2)
				temp2=((1/a**2)*(h1**2+k1**2))+((1/c**2)*l1**2)
				temp3=((1/a**2)*(h2**2+k2**2))+((1/c**2)*l2**2)
				temp4=sqrt(temp2*temp3)

				costheta=temp1/temp4
				angle=acos(costheta)
				angle=angle*(180/pi())
				
				number temp=((1/a**2)*(h1**2+k1**2))+((1/c**2)*l1**2)
				d1=sqrt(1/temp)
				
				temp=((1/a**2)*(h2**2+k2**2))+((1/c**2)*l2**2)
				d2=sqrt(1/temp)
				
				crystalname="Tetragonal"
			}

		if(crystalvalue==3) //Orthogonal
			{
				number costheta, temp1, temp2, temp3, temp4
				temp1=((1/a**2)*h1*h2)+((1/b**2)*k1*k2)+((1/c**2)*l1*l2)
				temp2=((1/a**2)*h1**2)+((1/b**2)*k1**2)+((1/c**2)*l1**2)
				temp3=((1/a**2)*h2**2)+((1/b**2)*k2**2)+((1/c**2)*l2**2)
				temp4=sqrt(temp2*temp3)

				costheta=temp1/temp4
				angle=acos(costheta)
				angle=angle*(180/pi())
				
				number temp=(((1/a**2)*h1**2)+((1/b**2)*k1**2)+((1/c**2)*l1**2))
				d1=sqrt(1/temp)
				
				temp=(((1/a**2)*h2**2)+((1/b**2)*k2**2)+((1/c**2)*l2**2))
				d2=sqrt(1/temp)
				
				crystalname="Orthorhombic"
			}

		if(crystalvalue==4) // Hexagonal
			{
				number costheta, temp1, temp2, temp3, temp4
				temp1=(h1*h2)+(k1*k2)+(0.5*((h1*k2)+(k1*h2)))+(((3*a**2)/(4*c**2))*l1*l2)
				temp2=h1**2+k1**2+(h1*k1)+(((3*a**2)/(4*c**2))*l1**2)
				temp3=h2**2+k2**2+(h2*k2)+(((3*a**2)/(4*c**2))*l2**2)
				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4
				
				number temp=((4/(3*a**2))*(h1**2+(h1*k1)+k1**2)+((1/c**2)*l1**2))
				d1=sqrt(1/temp)
				
				temp=((4/(3*a**2))*(h2**2+(h2*k2)+k2**2)+((1/c**2)*l2**2))
				d2=sqrt(1/temp)
				

				angle=acos(costheta)
				angle=angle*(180/pi())
				crystalname="Hexagonal"
			}

		if(crystalvalue==5) // Rhombohedral
			{
				number costheta, temp, temp1, temp2, temp3, temp4, cellvol, d1, d2
				temp=(1/a**2)
				temp3=(1+cos(alpha)-(2*(cos(alpha)**2)))
				temp2=(1+cos(alpha))*((h1**2+k1**2+l1**2)-((1-tan(0.5*alpha)**2)*((h1*k1)+(k1*l1)+(l1*h1))))
				temp4=(temp*temp2)/temp3
				d1=sqrt(1/temp4)

				temp2=(1+cos(alpha))*((h2**2+k2**2+l2**2)-((1-tan(0.5*alpha)**2)*((h2*k2)+(k2*l2)+(l2*h2))))
				temp4=(temp*temp2)/temp3
				d2=sqrt(1/temp4)

				cellvol=a**3*sqrt((1-(3*cos(alpha)**2)+(2*cos(alpha)**3)))

				temp1=(a**4*d1*d2)/cellvol**2
				temp2=sin(alpha)**2*((h1*h2)+(k1*k2)+(l1*l2))
				temp3=(cos(alpha)**2-cos(alpha))*((k1*l2)+(k2*l1)+(l1*h2)+(l2*h1)+(h1*k2)+(h2*k1))
				costheta=temp1*(temp2+temp3)

				angle=acos(costheta)
				angle=angle*(180/pi())
				
				number tem1=((h1**2+k1**2+l1**2)*sin(alpha)**2)+(2*((h1*k1)+(k1*l1)+(h1*l1))*(cos(alpha)**2-cos(alpha)))
				number tem2=a**2*(1-(3*cos(alpha)**2)+(2*cos(alpha)**3))
				number tem3=tem1/tem2
				d1=sqrt(1/tem3)
				
				tem1=((h1**2+k1**2+l1**2)*sin(alpha)**2)+(2*((h1*k1)+(k1*l1)+(h1*l1))*(cos(alpha)**2-cos(alpha)))
				tem2=a**2*(1-(3*cos(alpha)**2)+(2*cos(alpha)**3))
				tem3=tem1/tem2
				d2=sqrt(1/tem3)
				
				crystalname="Rhombohedral"
			}

		if(crystalvalue==6) // Monoclinic
			{
				number costheta, temp1, temp2, temp3, temp4

				temp1=((1/a**2)*h1*h2)+((1/b**2)*k1*k2*sin(beta)**2)+((1/c**2)*l1*l2)-((1/(a*c))*((l1*h2)+(l2*h1))*cos(beta))
				temp2=((1/a**2)*h1**2)+((1/b**2)*k1**2*sin(beta)**2)+((1/c**2)*l1**2)-(((2*h1*l1)/(a*c))*cos(beta))
				temp3=((1/a**2)*h2**2)+((1/b**2)*k2**2*sin(beta)**2)+((1/c**2)*l2**2)-(((2*h2*l2)/(a*c))*cos(beta))	
				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4

				angle=acos(costheta)
				angle=angle*(180/pi())
				
				number temp=((1/a**2)*(h1**2/sin(beta)**2))+((1/b**2)*k1**2)+((1/c**2)*(l1**2/sin(beta)**2))-((2*h1*l1*cos(beta))/(a*c*sin(beta)**2))
				d1=sqrt(1/temp)
				
				temp=((1/a**2)*(h2**2/sin(beta)**2))+((1/b**2)*k2**2)+((1/c**2)*(l2**2/sin(beta)**2))-((2*h2*l2*cos(beta))/(a*c*sin(beta)**2))
				d2=sqrt(1/temp)
				
				crystalname="Monoclinic"
			}

		if(crystalvalue==7) // triclinic
			{
				number costheta, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12

				temp1=(h1*h2*b**2*c**2*sin(alpha)**2)+(k1*k2*a**2*c**2*sin(beta)**2)+(l1*l2*a**2*b**2*sin(gamma)**2)
				temp2=a*b*c**2*((cos(alpha)*cos(beta))-cos(gamma))*((k1*h2)+(h1*k2))
				temp3=a*b**2*c*((cos(gamma)*cos(alpha))-cos(beta))*((h1*l2)+(l1*h2))
				temp4=a**2*b*c*((cos(beta)*cos(gamma))-cos(alpha))*((k1*l2)+(l1*k2))

				number fval=temp1+temp2+temp3+temp4

				temp5=(h1**2*b**2*c**2*sin(alpha)**2)+(k1**2*a**2*c**2*sin(beta)**2)+(l1**2*a**2*b**2*sin(gamma)**2)
				temp6=2*h1*k1*a*b*c**2*((cos(alpha)*cos(beta))-cos(gamma))
				temp7=2*h1*l1*a*b**2*c*((cos(gamma)*cos(alpha))-cos(beta))
				temp8=2*k1*l1*a**2*b*c*((cos(beta)*cos(gamma))-cos(alpha))

				number ahikili=sqrt(temp5+temp6+temp7+temp8)

				temp9=(h2**2*b**2*c**2*sin(alpha)**2)+(k2**2*a**2*c**2*sin(beta)**2)+(l2**2*a**2*b**2*sin(gamma)**2)
				temp10=2*h2*k2*a*b*c**2*((cos(alpha)*cos(beta))-cos(gamma))
				temp11=2*h2*l2*a*b**2*c*((cos(gamma)*cos(alpha))-cos(beta))
				temp12=2*k2*l2*a**2*b*c*((cos(beta)*cos(gamma))-cos(alpha))

				number ah2k2l2=sqrt(temp9+temp10+temp11+temp12)

				costheta=fval/(ahikili*ah2k2l2)
				angle=acos(costheta)
				angle=angle*(180/pi())
				
				number V2=a**2*b**2*c**2*(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+(2*cos(alpha)*cos(beta)*cos(gamma)))
				number s11=b**2*c**2*sin(alpha)**2
				number s22=a**2*c**2*sin(beta)**2
				number s33=a**2*b**2*sin(gamma)**2

				number s12=a*b*c**2*(cos(alpha)*cos(beta)-cos(gamma))
				number s23=a**2*b*c*(cos(beta)*cos(gamma)-cos(alpha))
				number s31=a*b**2*c*(cos(gamma)*cos(alpha)-cos(beta))

				number temp=(1/v2)*((s11*h1**2)+(s22*k1**2)+(s33*l1**2)+(2*s12*h1*k1)+(2*s23*k1*l1)+(2*s31*l1*h1))
				d1=sqrt(1/temp)
				
				V2=a**2*b**2*c**2*(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+(2*cos(alpha)*cos(beta)*cos(gamma)))
				s11=b**2*c**2*sin(alpha)**2
				s22=a**2*c**2*sin(beta)**2
				s33=a**2*b**2*sin(gamma)**2

				s12=a*b*c**2*(cos(alpha)*cos(beta)-cos(gamma))
				s23=a**2*b*c*(cos(beta)*cos(gamma)-cos(alpha))
				s31=a*b**2*c*(cos(gamma)*cos(alpha)-cos(beta))

				temp=(1/v2)*((s11*h2**2)+(s22*k2**2)+(s33*l2**2)+(2*s12*h2*k2)+(2*s23*k2*l2)+(2*s31*l2*h2))
				d2=sqrt(1/temp)
				
				crystalname="Triclinic"
			}

//----------------------and then calculate the zone axis direction-------------------------
		//skip for 000
		if((h1**2+k1**2+l1**2)==0 || (h2**2+k2**2+l2**2)==0)
			{
				beep()
				exit(0)
			}

		u=(k1*l2)-(k2*l1)
		v=(l1*h2)-(l2*h1)
		w=(h1*k2)-(h2*k1)

		// Normalization
		number divisor=CrystUIFrame.NormaliseIndices(u,v,w)
		u=round(u*100)/100
		v=round(v*100)/100
		w=round(w*100)/100
	
		self.lookupElement("angleresultfield").dlgtitle(""+ format(angle, "%6.1f"+" deg." )) 
		self.lookupElement("angleresultfield1").dlgtitle("["+u+","+v+","+w+"]" )

			angle=round(angle*100)/100
			result("\nd1="+format(d1, "%6.4f")+" A, d2="+format(d2, "%6.4f")+" A\n")
			result("("+h1+","+k1+","+l1+") <==> ("+h2+","+k2+","+l2+") ==> Angle="+angle+" deg.\n")
			result("("+h1+","+k1+","+l1+") <==> ("+h2+","+k2+","+l2+") ==> ["+u+","+v+","+w+"]\n")
			self.validateview()
		
}


//***************************Function: directions pannel***************************************
	void calculatedirectionangles(object self)
	{
		number crystalvalue, alpha, beta, gamma, a, b, c, u1,v1,t1,w1, u2,v2,t2,w2, h,k,l, r_1,r_2,angle, ahex, chex, arhomb, rhombalpha, xu1, xu2, xv1, xv2, xt1, xt2, xw1, xw2
		number p1, p2, q1, q2, s1, s2, r1, r2, hu1, hu2, hv1, hv2, hw1, hw2
		string crystalname, material

		self.dlggetvalue("crystalradios",crystalvalue)
		self.dlggetvalue("alphafield",alpha)
		self.dlggetvalue("betafield",beta)
		self.dlggetvalue("gammafield",gamma)
		self.dlggetvalue("afield",a)
		self.dlggetvalue("bfield",b)
		self.dlggetvalue("cfield",c)

		self.dlggetvalue("diru1field",u1)
		self.dlggetvalue("dirv1field",v1)
		self.dlggetvalue("dirw1field",w1)
		
		self.dlggetvalue("diru2field",u2)
		self.dlggetvalue("dirv2field",v2)
		self.dlggetvalue("dirw2field",w2)

		alpha=alpha/(180/pi())
		beta=beta/(180/pi())
		gamma=gamma/(180/pi())

		//calculate the size of a given lattice direction [uvw]
		r_1=sqrt(u1**2*a**2+v1**2*b**2+w1**2*c**2+2*v1*w1*b*c*cos(alpha)+2*w1*u1*a*c*cos(beta)+2*u1*v1*a*b*cos(gamma) )
		r_2=sqrt(u2**2*a**2+v2**2*b**2+w2**2*c**2+2*v2*w2*b*c*cos(alpha)+2*w2*u2*a*c*cos(beta)+2*u2*v2*a*b*cos(gamma) )
			
			
		// calculate the angle of two vectors	
		if(crystalvalue==1) // cubic
			{
				number costheta, temp1, temp2
				temp1=(u1*u2)+(v1*v2)+(w1*w2)
				temp2=sqrt((u1**2+v1**2+w1**2)*(u2**2+v2**2+w2**2))
				costheta=temp1/temp2

				angle=acos(costheta)
				angle=angle*(180/pi())
				crystalname="Cubic"
			}


		if(crystalvalue==2) //tetragonal
			{
				number costheta, temp1, temp2, temp3, temp4

				temp1=a**2*((u1*u2)+(v1*v2))+(c**2*w1*w2)
				temp2=(a**2*(u1**2+v1**2))+(c**2*w1**2)
				temp3=(a**2*(u2**2+v2**2))+(c**2*w2**2)
				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4

				angle=acos(costheta)
				angle=angle*(180/pi())
				crystalname="Tetragonal"
			}


		if(crystalvalue==3) // orthorhombic
			{
				number costheta, temp1, temp2, temp3, temp4
				temp1=(a**2*u1*u2)+(b**2*v1*v2)+(c**2*w1*w2)
				temp2=(a**2*u1**2)+(b**2*v1**2)+(c**2*w1**2)
				temp3=(a**2*u2**2)+(b**2*v2**2)+(c**2*w2**2)
				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4

				angle=acos(costheta)
				angle=angle*(180/pi())
				crystalname="Orthorhombic"
			}


		if(crystalvalue==4) // hexagonal
			{
				number costheta, temp1, temp2, temp3, temp4

				//convert four indices to three indices			
				t1=-(u1+v1)
				t2=-(u2+v2)

				xu1=u1-t1
				xu2=u2-t2
				xv1=v1-t1
				xv2=v2-t2

				//w1=xw1,w2=xw2
				temp1=(xu1*xu2)+(xv1*xv2)-(0.5*((xu1*xv2)+(xv1*xu2)))+((c**2/a**2)*w1*w2)
				temp2=xu1**2+xv1**2-(xu1*xv1)+((c**2/a**2)*w1**2)
				temp3=xu2**2+xv2**2-(xu2*xv2)+((c**2/a**2)*w2**2)

				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4

				angle=acos(costheta)
				angle=angle*(180/pi())
				crystalname="Hexagonal"
			}

		if(crystalvalue==5) // Rhombohedral
			{

				number costheta, temp1, temp2, temp3, temp4

				// convertion of indices
				xu1=((2*u1)-v1-w1)/3
				xu2=((2*u2)-v2-w2)/3

				xv1=(u1+v1-(2*w1))/3
				xv2=(u2+v2-(2*w2))/3

				xw1=(u1+v1+w1)/3
				xw2=(u2+v2+w2)/3

				p1=((2*xu1)-xv1)
				q1=((2*xv1)-xu1)

				s1=xw1*3
				normaliseindices(self, p1, q1, s1)
				r1=-(p1+q1)

				p2=((2*xu2)-xv2)
				q2=((2*xv2)-xu2)

				s2=xw2*3
				normaliseindices(self, p2, q2, s2)
				r2=-(p2+q2)

				//convert rhombohedron to hexagonal
				rhombtohexconversion(self, alpha, a, ahex, chex)

				temp1=(xu1*xu2)+(xv1*xv2)-(0.5*((xu1*xv2)+(xv1*xu2)))+((chex**2/ahex**2)*xw1*xw2)
				temp2=xu1**2+xv1**2-(xu1*xv1)+((chex**2/ahex**2)*xw1**2)
				temp3=xu2**2+xv2**2-(xu2*xv2)+((chex**2/ahex**2)*xw2**2)

				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4

				angle=acos(costheta)
				angle=angle*(180/pi())
				crystalname="Rhombohedral"

				normaliseindices(self, xu1, xv1, xw1)
				xt1=-(xu1+xv1)
				normaliseindices(self, xu2, xv2, xw2)
				xt2=-(xu2+xv2)
			}

		if(crystalvalue==6) // monoclinic
			{
				number costheta, temp1, temp2, temp3, temp4
				temp1=(a**2*u1*u2)+(b**2*v1*v2)+(c**2*w1*w2)+(a*c*((w1*u2)+(u1*w2))*cos(beta))
				temp2=(a**2*u1**2)+(b**2*v1**2)+(c**2*w1**2)+(2*a*c*u1*w1*cos(beta))
				temp3=(a**2*u2**2)+(b**2*v2**2)+(c**2*w2**2)+(2*a*c*u2*w2*cos(beta))
				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4

				angle=acos(costheta)
				angle=angle*(180/pi())
				crystalname="Monoclinic"
			}



		if(crystalvalue==7) //triclinic
			{
				number costheta, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12

				temp1=(a**2*u1*u2)+(b**2*v1*v2)+(c**2*w1*w2)
				temp2=b*c*((v1*w2)+(w1*v2))*cos(alpha)
				temp3=a*c*((w1*u2)+(u1*w2))*cos(beta)
				temp4=a*b*((u1*v2)+(v1*u2))*cos(gamma)

				number lval=temp1+temp2+temp3+temp4

				temp5=(a**2*u1**2)+(b**2*V1**2)+(c**2*w1**2)
				temp6=2*b*c*v1*w1*cos(alpha)
				temp7=2*c*a*w1*u1*cos(beta)
				temp8=2*a*b*u1*v1*cos(gamma)

				number iuvw1=sqrt(temp5+temp6+temp7+temp8)

				temp9=(a**2*u2**2)+(b**2*V2**2)+(c**2*w2**2)
				temp10=2*b*c*v2*w2*cos(alpha)
				temp11=2*c*a*w2*u2*cos(beta)
				temp12=2*a*b*u2*v2*cos(gamma)

				number iuvw2=sqrt(temp9+temp10+temp11+temp12)
				costheta=lval/(iuvw1*iuvw2)

				angle=acos(costheta)
				angle=angle*(180/pi())
				crystalname="Triclinic"
			}
			
		angle=round(angle*100)/100
		
//----------------------and then calculate the formed plane-------------------------
		h=(v1*w2)-(v2*w1)
		k=(w1*u2)-(w2*u1)
		l=(u1*v2)-(u2*v1)

		// normalization
		number divisor=CrystUIFrame.NormaliseIndices(h,k,l)
		h=round(h*100)/100
		k=round(k*100)/100
		l=round(l*100)/100
		
		//show results	
		self.lookupElement("dirresultfield").dlgtitle(""+ format(angle, "%6.1f"+" deg." )) 
		self.lookupElement("dirresultfield1").dlgtitle("("+h+","+k+","+l+")" )

		result("\nr1="+format(r_1, "%6.4f")+" A, r2="+format(r_2, "%6.4f")+" A\n")
		result("["+u1+","+v1+","+w1+"] <==> ["+u2+","+v2+","+w2+"] ==> Angle="+angle+" deg.\n")
		result("["+u1+","+v1+","+w1+"] <==> ["+u2+","+v2+","+w2+"] ==> ("+h+","+k+","+l+")\n")

		self.validateview()
}


//***************************Function: plane or direction indices convertion***************************************
void plane2direction(object self)
{
number crystalvalue, alpha, beta, gamma, a, b, c
		self.dlggetvalue("crystalradios",crystalvalue)
		self.dlggetvalue("alphafield",alpha)
		self.dlggetvalue("betafield",beta)
		self.dlggetvalue("gammafield",gamma)
		self.dlggetvalue("afield",a)
		self.dlggetvalue("bfield",b)
		self.dlggetvalue("cfield",c)
		
		alpha=alpha/(180/pi())
		beta=beta/(180/pi())
		gamma=gamma/(180/pi())		

        number h,k,l,u,v,w		
		self.dlggetvalue("zoneh1field",h)
		self.dlggetvalue("zonek1field",k)
		self.dlggetvalue("zonel1field",l)
		
		self.dlggetvalue("zoneh2field",u)
		self.dlggetvalue("zonek2field",v)
		self.dlggetvalue("zonel2field",w)		
		
//define G matrix and its inverse
number a11,a12,a13,a21,a22,a23,a31,a32,a33
a11=a**2
a12=Tert(abs(a*b*cos(gamma))<1e-8,0,a*b*cos(gamma) )
a13=tert(abs(a*c*cos(beta))<1e-8,0,a*c*cos(beta))
a21=Tert(abs(a*b*cos(gamma))<1e-8,0,a*b*cos(gamma))
a22=b**2
a23=Tert(abs(b*c*cos(alpha))<1e-8,0,b*c*cos(alpha))
a31=Tert(abs(a*c*cos(beta))<1e-8,0,a*c*cos(beta))
a32=Tert(abs(b*c*cos(alpha))<1e-8,0,b*c*cos(alpha))
a33=c**2 
image G := [3,3] : { {a11,a12,a13},{a21,a22,a23},{a31,a32,a33} }

number b11,b12,b13,b21,b22,b23,b31,b32,b33,t
image G_1 :=MatrixInverse(G)
b11=G_1.GetPixel(0,0)
b12=G_1.GetPixel(0,1)
b13=G_1.GetPixel(0,2)
b21=G_1.GetPixel(1,0)
b22=G_1.GetPixel(1,1)
b23=G_1.GetPixel(1,2)
b31=G_1.GetPixel(2,0)
b32=G_1.GetPixel(2,1)
b33=G_1.GetPixel(2,2)

Result("\n---------------------------------- The direct metric tensors ----------------------------------\n")
result(a11+"\t\t\t"+a12+"\t\t\t"+a13+"\n"+a21+"\t\t\t"+a22+"\t\t\t"+a23+"\n"+a31+"\t\t\t"+a32+"\t\t\t"+a33+"\n")
result("-------------------------------- The reciprocal metric tensors --------------------------------\n")
result(b11+"\t\t\t"+b12+"\t\t\t"+b13+"\n"+b21+"\t\t\t"+b22+"\t\t\t"+b23+"\n"+b31+"\t\t\t"+b32+"\t\t\t"+b33+"\n")
Result("-----------------------------------------------------------------------------------------------\n")

image direction:=[1,3]:{ {u},{v},{w} }
image plane:=[1,3]:{ {h},{k},{l} }

image hkl:=MatrixMultiply(G,direction)
image uvw:=MatrixMultiply(G_1,plane)

number h1=hkl.GetPixel(0,0)
number k1=hkl.GetPixel(0,1)
number l1=hkl.GetPixel(0,2)

image hkl_1=abs(hkl) 
number x,y
for(number i=1;i<4;i++)
{
if(hkl_1.min(x,y)==0)hkl_1.SetPixel(x,y,1e8)
}

number divisor=hkl_1.min()

h1=h1/divisor
k1=k1/divisor
l1=l1/divisor

number u1=uvw.GetPixel(0,0)
number v1=uvw.GetPixel(0,1)
number w1=uvw.GetPixel(0,2)

image uvw_1=abs(uvw)
for(number i=1;i<4;i++)
{
if(uvw_1.min(x,y)==0)uvw_1.SetPixel(x,y,1e8)
}

divisor=uvw_1.min()

u1=u1/divisor
v1=v1/divisor
w1=w1/divisor


if(h**2+l**2+u**2!=0)result("\nPlane to direction:  ("+h+", "+k+", "+l+") <==>  ["+u1+", "+v1+", "+w1+"]\n")
if(u**2+v**2+w**2!=0)result("\nDirection to plane:  ["+u+", "+v+", "+w+"] <==> ("+h1+", "+k1+", "+l1+") \n")	

dlgtitle(self.lookupelement("zoneresultfield"),""+format(u1,"%4.2f")+", "+format(v1,"%4.2f")+", "+format(w1,"%4.2f"))
dlgtitle(self.lookupelement("zoneresultfield1"),""+format(h1,"%4.2f")+", "+format(k1,"%4.2f")+", "+format(l1,"%4.2f"))
}
	
	//*************************Function: to response when d or hkl is changed in d-hkl pannel*******************************************
	void dhorkchanged(object self, taggroup tg)
	{
		number hvalue, kvalue, ivalue, crystalvalue
		self.dlggetvalue("crystalradios",crystalvalue)        
        
		self.validateview()       
	}
	
	
	//*********************************Function: to response when hkl_1 is changed in Planes pannel***********************************
	void anghork1changed(object self, taggroup tg)
	{
		number hvalue, kvalue, ivalue, crystalvalue
		self.dlggetvalue("crystalradios",crystalvalue)       
        
		self.validateview()       
	}	
	
	
	//*********************************Function: to response when hkl_2 is changed in Planes pannel***********************************
	void anghork2changed(object self, taggroup tg)
	{
		number hvalue, kvalue, ivalue, crystalvalue
		self.dlggetvalue("crystalradios",crystalvalue)       
        
		self.validateview()       
	}	
	
	
		//*********************************Function: to response when uvw_1 is changed in Direction sub-pannel***********************************
	void diru1orv1changed(object self, taggroup tg)
	{
		number uvalue, vvalue, tvalue, crystalvalue
		self.dlggetvalue("crystalradios",crystalvalue)      
        
		self.validateview()       
	}	
	
		//*********************************Function: to response when uvw_2 is changed in Direction sub-pannel***********************************
	void diru2orv2changed(object self, taggroup tg)
	{
		number uvalue, vvalue, tvalue, crystalvalue
		self.dlggetvalue("crystalradios",crystalvalue)      
        
		self.validateview()       
	}	
	
	//*********************************Function: to response when (hkl) is changed in h<->u sub-pannel***********************************
	void zonehork1changed(object self, taggroup tg)
	{
		number hvalue, kvalue, ivalue, crystalvalue
		self.dlggetvalue("crystalradios",crystalvalue)
		
		self.validateview()
	}
	
	//*********************************Function: to response when [uvw] is changed in h<->u sub-pannel***********************************
	void zonehork2changed(object self, taggroup tg)
	{
		number hvalue, kvalue, ivalue, crystalvalue
		self.dlggetvalue("crystalradios",crystalvalue)
		
		self.validateview()
	}	
	
	
	//*********************************Function: to response when Crystal list is changed***********************************
	void elementchanged(object self, taggroup tg)
	{
		number flag, a, b, c, alpha, beta, gamma, pulldownmenuval, i, lattice

		self.clearfields()
		pulldownmenuval=dlggetvalue(tg)

		string pdname
		dlggetnthlabel(tg, pulldownmenuval-1, pdname)

		if(pdname!="") // select a crystal
			{
				self.setelementisenabled("crystalradios",1)
				self.setelementisenabled("alphafield",1)
				self.setelementisenabled("betafield",1)
				self.setelementisenabled("gammafield",1)
				self.setelementisenabled("afield",1)
				self.setelementisenabled("bfield",1)
				self.setelementisenabled("cfield",1)
				self.setelementisenabled("addcrystalbutton",1)
			}
		else  
			{
				self.setelementisenabled("crystalradios",1)
				self.setelementisenabled("alphafield",1)
				self.setelementisenabled("betafield",1)
				self.setelementisenabled("gammafield",1)
				self.setelementisenabled("afield",1)
				self.setelementisenabled("bfield",1)
				self.setelementisenabled("cfield",1)
				self.setelementisenabled("addcrystalbutton",1)

				crystUIFrame.crystclearbutton()
			}

		// ----------------------------read and set lattice parameters of the selected crystal---------------------------
		flag=0

		string crystalname

		taggroup ptg=getpersistenttaggroup()
		taggroup crystaltags
		ptg.TaggroupGetTagAsTagGroup("ElectronDiffraction Tools:Crystal Tools:Crystals", crystaltags)

		number count=crystaltags.taggroupcounttags()

		for(i=0; i<count; i++)
		{
			String label = crystaltags.TagGroupGetTagLabel( i ) 
			if(pdname==label) //find a crystal
			{
				flag=1

				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+label+":Lattice", lattice)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+label+":a", a)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+label+":b", b)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+label+":c", c)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+label+":alpha", alpha)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+label+":beta", beta)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+label+":gamma", gamma)

				break
			}

			if(flag==1) break
		}

		//update the paramters
		if(flag==1) 
			{
				self.dlgvalue("crystalradios",lattice)
				self.validateview()
				
				self.dlgvalue("alphafield",alpha)
				self.validateview()
				self.dlgvalue("betafield",beta)
				self.validateview()
				self.dlgvalue("gammafield",gamma)
				self.validateview()
				
				self.dlgvalue("afield",a)
				self.validateview()
				self.dlgvalue("bfield",b)
				self.validateview()
				self.dlgvalue("cfield",c)
				self.validateview()
			}

		self.validateview()
	}


	// -----------------Function: to response when sub-pannels are switched------------------
	void maincalcroutine(object self)
	{
		number tabvalue
		self.dlggetvalue("crysttabs", tabvalue) // 0=d-spacings, 1=planes

		if(tabvalue==0) // d-hkl sub-pannel
			{
				self.calculatedspacing()
				exit(0)
			}

		if(tabvalue==1) // Planes sub-pannel
			{
				self.calculateplaneangles()
				exit(0)
			}
		if(tabvalue==2) // Directions sub-pannel
			{
				self.calculatedirectionangles()
				exit(0)
			}
		if(tabvalue==4) // h<->u sub-pannel
			{
				self.plane2direction()
				exit(0)
			}
}


//---------------------------------Function: + button, to add or delete crystals-----------------------------
void AddCrystal(object self)
{
If(OptionDown())  //ALT+: to delete the selected crystal
{
TagGroup Tgcrystals= self.LookUpElement("elementpulldown")
number CrystalName=Tgcrystals.dlggetvalue()  //note: it begins from 1 instead 0

taggroup ptg=getpersistenttaggroup()
taggroup crystaltags
ptg.TaggroupGetTagAsTagGroup("ElectronDiffraction Tools:Crystal Tools:Crystals", crystaltags)
String label = crystaltags.TagGroupGetTagLabel( CrystalName-1) 

if(!twobuttondialog("The crystal :    "+label+"    ,  will be removed!","Yes","No")) exit(0)

crystaltags.TagGroupDeleteTagWithIndex(CrystalName-1) 
Result("\nCrystal: "+label+" is removed from the database.\n")

//delete all pull-down list
for(number i=0;i<200;i++)
{
Tgcrystals.DLGRemoveElement(0)
}
	
//update the tags
ptg.TaggroupGetTagAsTagGroup("ElectronDiffraction Tools:Crystal Tools:Crystals", crystaltags)

number count=crystaltags.taggroupcounttags()
for(number i=0; i<count; i++)
{
String label = crystaltags.TagGroupGetTagLabel( i ) 
Tgcrystals.DLGAddPopupItemEntry(label)
}	

}

//Add a crystal
else
{
number lattice, a, b, c, alpha, beta, gamma
string crystalname

self.dlggetvalue("crystalradios",lattice)
self.dlggetvalue("alphafield",alpha)
self.dlggetvalue("betafield",beta)
self.dlggetvalue("gammafield",gamma)
self.dlggetvalue("afield",a)
self.dlggetvalue("bfield",b)
self.dlggetvalue("cfield",c)
self.dlggetvalue("crystalnamefield", crystalname)

		if(a==0 || b==0 || c==0)
		{
			showalert("Lattice parameters can not be zero!",1)
			exit(0)
		}

		if(alpha==0 || beta==0 || gamma==0)
		{
			showalert("Cell angles can not be zero!",1)
			exit(0)
		}

	taggroup ptg=getpersistenttaggroup()
	taggroup crystaltags
	ptg.TaggroupGetTagAsTagGroup("ElectronDiffraction Tools:Crystal Tools:Crystals", crystaltags)

	number count=crystaltags.taggroupcounttags()
	number i, label
	for(i=0; i<count; i++)
	{
		String label = crystaltags.TagGroupGetTagLabel( i ) 

		if(crystalname==label||crystalname==""||crystalname==" "||crystalname=="  ")
			{
				beep()
				if(!twobuttondialog("The crystal :    "+crystalname+"    ,  is already in the database, \nOr the crystal name is   None   !","Overwrite","Rename")) GetString("Crystal name?",crystalname,crystalname)
			}
	}
	
	string latticename
	if(lattice==1) latticename="Cubic"
	if(lattice==2) latticename="Tetragonal"
	if(lattice==3) latticename="Orthorhombic"
	if(lattice==4) latticename="Hexagonal"
	if(lattice==5) latticename="Rhombohedral"
	if(lattice==6) latticename="Monoclinic"
	if(lattice==7) latticename="Triclinic"
	
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":Lattice",lattice)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":a",a)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":b",b)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":c",c)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":alpha",alpha)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":beta",beta)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Crystals:"+crystalname+":gamma",gamma)

Result("\nThe added crystal is: "+CrystalName+" ("+latticename+", "+a+","+b+", "+c+", "+alpha+", "+beta+", "+gamma+")\n")

//update tags
	count=crystaltags.taggroupcounttags()
	number selectedno
	for(i=0; i<count; i++)
	{
		String label = crystaltags.TagGroupGetTagLabel( i ) 

		if(crystalname==label)
			{
			selectedno=i
			result(selectedno+"\n")
			break
			}
	}
	
TagGroup Tgcrystals= self.LookUpElement("elementpulldown")
for(i=0;i<200;i++)
{
Tgcrystals.DLGRemoveElement(0)
}
	
ptg.TaggroupGetTagAsTagGroup("ElectronDiffraction Tools:Crystal Tools:Crystals", crystaltags)
count=crystaltags.taggroupcounttags()
for(i=0; i<count; i++)
{
String label = crystaltags.TagGroupGetTagLabel( i ) 

if(i==selectedno)Tgcrystals.DLGAddPopupItemEntry(label,1)

else
{
Tgcrystals.DLGAddPopupItemEntry(label)
}
}	
}
}


//------------SAED buttion: to calculate FT diffractogram from HRTEM image----------------
void SAEDButton(object self)
{
Result("-----------------------------------------------------------------------------------------\nSAED button: \n1) Get the diffractogram of a HRTEM pattern;\n2) ALT key: Get the filtered diffractogram of a HRTEM pattern;\n3) CTRL key: Get an Autocorrelation SAED pattern\.n-----------------------------------------------------------------------------------------\n")

if(ControlDown())
{
image img:=getfrontimage()

number xsize,ysize,size,n,x,y,t,l,b,r
img.getsize(xsize,ysize)

//Keep the image size 2^n
size=max(xsize,ysize)
n=ceil(log2(size))
size=2**n

l=(size-xsize)/2
t=(size-ysize)/2

image roiImg=realimage("",4,size,size)
roiImg=0
roiImg[t,l,t+ysize,l+xsize]=img
roiImg=roiImg.AutoCorrelate()

x=(size+1)/2-l
y=(size+1)/2-t

image new=roiImg[t,l,t+ysize,l+xsize]
new=new-new.min()
new.ImageCopyCalibrationFrom(img)
new.ShowImage()
SetNumberNote(new,"ElectronDiffraction Tools:Center:X", x)
SetNumberNote(new,"ElectronDiffraction Tools:Center:Y", y)
new.SetOrigin(x,y)

new.SetName("Autocorrelation-"+img.GetName())
new.setstringnote("ElectronDiffraction Tools:File type","SAED")
}

else if(OptionDown())
{
image front:=getfrontimage()
string imgname=getname(front)

number t,l,b,r
front.GetSelection(t,l,b,r)

number x0,y0
x0=0.5*(r+l)
y0=0.5*(t+b)

number width=Minimum(b-t,r-l)
width=width/2
front.SetSelection(y0-width,x0-width,y0+width,x0+width)

image temp:=front[].ImageClone()

number xsize,ysize,xscale,yscale
temp.GetSize(xsize,ysize)
front.GetScale(xscale,yscale)

number no=round(log2(xsize))
number factor=2**no/xsize

xsize=2**no
ysize=xsize
image temp1=realimage("",4,xsize,xsize)
temp1=warp(temp,icol/factor,irow/factor)

xscale=xscale/factor

// Create the damped edges filter
number edgesize=100 // the width in pixels over which the edge is smoothed
if(edgesize>(xsize)/2) edgesize=xsize/3
if(edgesize>(ysize)/2) edgesize=ysize/3

image filter=createdampededge(edgesize, xsize, ysize)
image HanImg=temp1*filter

compleximage FFTimg=realFFT(HanImg)
image magimg=sqrt(real(FFTimg)**2+imaginary(FFTimg)**2)

xscale=1/(xsize*xscale)

//interpolation of the image if the size is so small
if(xsize<1024)
{
number factor=xsize/1024

image new=realimage("",4,1024,1024)
new=warp(magimg,icol*factor,irow*factor)

new.GetSize(xsize,ysize)
new.setname("Diffractogram-"+imgname)
new.showimage()

if(OkCancelDialog("AutoCorrelate?")==1)new=new.AutoCorrelate()

xscale=xscale*factor
new.SetScale(xscale,xscale)
new.SetOrigin((xsize+1)/2,(xsize+1)/2)
new.ImageSetDimensionUnitString(0,"1/nm")
new.ImageSetDimensionUnitString(1,"1/nm")
SetNumberNote(new,"ElectronDiffraction Tools:Center:X", (xsize+1)/2 )
SetNumberNote(new,"ElectronDiffraction Tools:Center:Y", (xsize+1)/2 )

new.setstringnote("ElectronDiffraction Tools:File type","SAED")
}

else
{
magimg.setname("Hanning diffractogram-"+imgname)
magimg.showimage()

if(OkCancelDialog("AutoCorrelate?")==1)magimg=magimg.AutoCorrelate()

number  xscale=1/((r-l)*front.ImageGetDimensionScale(0))
magimg.SetScale(xscale,xscale)
magimg.SetOrigin((r-l+1)/2,(r-l+1)/2)
magimg.ImageSetDimensionUnitString(0,"1/nm")
magimg.ImageSetDimensionUnitString(1,"1/nm")
SetNumberNote(magimg,"ElectronDiffraction Tools:Center:X", (r-l+1)/2 )
SetNumberNote(magimg,"ElectronDiffraction Tools:Center:Y", (r-l+1)/2 )

magimg.setstringnote("ElectronDiffraction Tools:File type","SAED")
}

}

else
{
image front:=getfrontimage()
string imgname=getname(front)

number t,l,b,r
front.GetSelection(t,l,b,r)

number x0,y0
x0=0.5*(r+l)
y0=0.5*(t+b)

number width=Minimum(b-t,r-l)
width=width/2
front.SetSelection(y0-width,x0-width,y0+width,x0+width)

image temp:=front[].ImageClone()

number xsize,ysize,xscale,yscale
temp.GetSize(xsize,ysize)
front.GetScale(xscale,yscale)

number no=round(log2(xsize))
number factor=2**no/xsize

xsize=2**no
ysize=xsize
image temp1=realimage("",4,xsize,xsize)
temp1=warp(temp,icol/factor,irow/factor)

xscale=xscale/factor
xscale=1/(xsize*xscale)

compleximage FFTimg=realFFT(temp1)
image magimg=sqrt(real(FFTimg)**2+imaginary(FFTimg)**2)

//interpolation of the image if the size is so small
if(xsize<1024)
{
number factor=xsize/1024

image new=realimage("",4,1024,1024)
new=warp(magimg,icol*factor,irow*factor)

new.GetSize(xsize,ysize)
new.setname("Diffractogram-"+imgname)
new.showimage()

if(OkCancelDialog("AutoCorrelate?")==1)new=new.AutoCorrelate()

xscale=xscale*factor
new.SetScale(xscale,xscale)
new.SetOrigin((xsize+1)/2,(xsize+1)/2)
new.ImageSetDimensionUnitString(0,"1/nm")
new.ImageSetDimensionUnitString(1,"1/nm")
SetNumberNote(new,"ElectronDiffraction Tools:Center:X", (xsize+1)/2 )
SetNumberNote(new,"ElectronDiffraction Tools:Center:Y", (xsize+1)/2 )

new.setstringnote("ElectronDiffraction Tools:File type","SAED")
}

else
{
magimg.setname("Diffractogram-"+imgname)
magimg.showimage()

if(OkCancelDialog("AutoCorrelate?")==1)magimg=magimg.AutoCorrelate()

number  xscale=1/((r-l)*front.ImageGetDimensionScale(0))
magimg.SetScale(xscale,xscale)
magimg.SetOrigin((r-l+1)/2,(r-l+1)/2)
magimg.ImageSetDimensionUnitString(0,"1/nm")
magimg.ImageSetDimensionUnitString(1,"1/nm")
SetNumberNote(magimg,"ElectronDiffraction Tools:Center:X", (r-l+1)/2 )
SetNumberNote(magimg,"ElectronDiffraction Tools:Center:Y", (r-l+1)/2 )

magimg.setstringnote("ElectronDiffraction Tools:File type","SAED")
}
}
}


/*************************************************************************
Center Button: to find the center of a SAED pattern
*************************************************************************/
void CenterButton(object self)
{
Result("-----------------------------------------------------------------------------------------\nCenter button: \n1) Locate the center of the 2D pattern;\n2) ALT key: to set the center by the oval component or the origin.\n-----------------------------------------------------------------------------------------\n")

if(OptionDown())
{
image img:=getfrontimage()
ImageDisplay imgdisp=img.imagegetimagedisplay(0)
setstringnote(img, "ElectronDiffraction Tools:File type","SAED")

number t,l,b,r,x0,y0,Radius
number annots=imgdisp.componentcountchildrenoftype(6)  
if(annots!=0)  //Oval component
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.ComponentGetBoundingRect(t,l,b,r)

x0=l+(r-l)/2
y0=t+(b-t)/2

if (!GetNumber( "Set New Center X (pix)\nX0= "+x0, x0, x0 ))exit(0)
if (!GetNumber( "Set New Center Y (pix)\nY0= "+y0, y0, y0 ))exit(0)

setorigin(img,x0,y0)
setNumberNote(img,"ElectronDiffraction Tools:Center:X", x0 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", y0)

Result("New center is: ("+x0+", "+y0+")\n")

ROI roi_1 = NewROI()
ROISetRectangle(roi_1, y0-5, x0-5, y0+5, x0+5 ) 
roi_1.roisetvolatile(0)
roi_1.roisetcolor(1,0,0)
imgdisp.ImageDisplayAddROI(roi_1)

setstringnote(img, "ElectronDiffraction Tools:File type","SAED")
}

else if(annots==0)  //input the given values
{
number sizex, sizey, centery, centerx, scalex, scaley,x,y,centerx0,centery0

img.getorigin(centerx0,centery0)

if (!GetNumber( "Set New Center X (pix)\nX0= "+centerx0, centerx0, centerx0 ))exit(0)
if (!GetNumber( "Set New Center Y (pix)\nY0= "+centery0, centery0, centery0 ))exit(0)

setorigin(img,centerx0,centery0)
setNumberNote(img,"ElectronDiffraction Tools:Center:X", centerx0 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", centery0)

Result("New center is: ("+centerx0+", "+centery0+")\n")

ROI roi_1 = NewROI()
ROISetRectangle(roi_1, centery0-5, centerx0-5, centery0+5, centerx0+5 ) 
roi_1.roisetvolatile(0)
roi_1.roisetcolor(1,0,0)
imgdisp.ImageDisplayAddROI(roi_1)

setstringnote(img, "ElectronDiffraction Tools:File type","SAED")
}
}
	
else  //Find it using Gaussian fitting
{	
image  img:=getfrontimage()
ImageDisplay imgdisp = img.ImageGetImageDisplay(0)
ImageDocument  imgdoc=getfrontimagedocument()
documentwindow imgwin=ImageDocumentGetWindow(imgdoc)
string imgname=img.getname()

number centrex, centrey, x0,y0,x1,y1,x2,y2,x3,y3,x,y,xsize,ysize,xscale,yscale, top,left,bottom,right,width,times
img.GetSize(xsize,ysize)

//clear all marks
number i
number annots= imgdisp.ImageDisplayCountROIS()
for (i=0; i<annots; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
 }

annots=imgdisp.componentcountchildrenoftype(13)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(2)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(2,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(5)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(6)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

number radius=xsize/150,smooth

if(radius<1)radius=2

component textannot=newtextannotation(imgdisp,xsize/15, ysize/15,"SPACE to mark the desired center!", 12)
textannot.componentsetfillmode(1)
textannot.componentsetdrawingmode(2) 
textannot.componentsetforegroundcolor(1,1,0)
textannot.componentsetfontfacename("Microsoft Sans Serif")
	
number k=1,m=1,noclick=1,twoclicks=1,n //noclick: number of clicks，twoclick：first click defines the center, second click defines the width
while(2>m)
{
number keypress=getkey()

getwindowsize(img, xsize, ysize)
 getscale(img,xscale,yscale)
 
number mouse_win_x, mouse_win_y
Number img_view_t, img_view_l, img_view_b, img_view_r
Number v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y
Number img_win_t, img_win_l, img_win_b, img_win_r
Number i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y
Number i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y
Number mouse_img_x, mouse_img_y

windowgetmouseposition(imgwin, mouse_win_x, mouse_win_y)
imgdisp.ImageDisplayGetImageRectInView( img_view_t, img_view_l, img_view_b, img_view_r )
imgdoc.ImageDocumentGetViewToWindowTransform( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y )
ObjectTransformTransformRect( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                            , img_view_t, img_view_l, img_view_b, img_view_r\
                            , img_win_t, img_win_l, img_win_b, img_win_r );

imgdisp.ComponentGetChildToViewTransform( i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y )
		
ObjectTransformCompose( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                      , i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y\
                      , i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y )
		
ObjectTransformUntransformPoint( i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y\
                               , mouse_win_x, mouse_win_y, mouse_img_x, mouse_img_y );

x1=mouse_img_x
y1=mouse_img_y

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)

		
if(keypress==32) //Space bar to define it
{
number t,l,b,r,leftpart,rightpart
img.GetSize(xsize,ysize)

Image profile_h:= LiveProfile_ExtractLineProfile( img, 0, y1, xsize, y1, 4 )  //horizental projection
Image profile_v:= LiveProfile_ExtractLineProfile( img,x1,0,x1,ysize , 4 )  //vertical projection

number width_h=fwhm(profile_h,x1,leftpart,rightpart)
x1=(leftpart+rightpart)/2

number width_v=fwhm(profile_v,y1,leftpart,rightpart)
y1=(leftpart+rightpart)/2

radius=1*(width_h+width_v)/2
Component box=NewBoxAnnotation(y1-0.5*radius,x1-0.5*radius,y1+0.5*radius,x1+0.5*radius)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)

SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":X", x1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Y", y1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Width", radius)
                                                                                 
if(noclick==1)                                                           
{                                                                              
radius=1*radius
}

else if(noclick==2)
{
radius=1*radius
}

img.ClearSelection()
noclick=noclick+1
}

if(keypress==29||keypress==30) //right and up key
{
radius=radius+0.2*radius
img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}

if(keypress==28||keypress==31) //left and down key
{
radius=radius-0.2*radius

if(radius<1)radius=1

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}
	
if(keypress>0 && keypress!=32&& keypress!=28&& keypress!=29&& keypress!=30&& keypress!=31||noclick==2) // Any key except space pressed - cancel
{
img.ClearSelection()
m=2
}

}


annots=imgdisp.componentcountchildrenoftype(5)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

for(number z=1;z<2;z++)  //refine the center
{
number centrex,centrey
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":X", centrex )
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":Y", centrey )
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":Width", width )
radius=1*width
width=2*width

img := GetFrontImage()
img.GetSize(xsize,ysize)

image temp=img[centrey-width,centrex-width,centrey+width,centrex+width]
top=centrey-width
left=centrex-width
temp.getsize(xsize,ysize)
temp=( ( icol -(xsize+1)/2 )**2 + ( irow - (ysize+1)/2 )**2 < radius**2 ) ? temp : 0 

number x_0,y_0
times=1
times=1/times
image new := RealImage("target",4,xsize/times,ysize/times)
new = temp.warp(icol*times,irow*times)

number fwhm,leftpart,rightpart

new=new.Gaussfilter(20)

new.GetSize(xsize,ysize)

image xproj=realimage("",4,xsize,1)
image yproj=realimage("",4,xsize,1)

xproj=LiveProfile_ExtractLineProfile( new, 0, (ysize+1)/2 , xsize-1, (ysize+1)/2 ,2 )
yproj= LiveProfile_ExtractLineProfile( new,(xsize+1)/2,0,(xsize+1)/2,ysize-1 , 2 ) 

xproj=xproj-xproj.min()
xproj=10000*xproj/xproj.max()+1

yproj=yproj-yproj.min()
yproj=10000*yproj/yproj.max()+1

number xfwhm=fwhm(xproj,  xsize/2, leftpart,rightpart)
number xpeak=(leftpart+rightpart)/2

number yfwhm=fwhm(yproj,  ysize/2, leftpart,rightpart)
number ypeak=(leftpart+rightpart)/2

number peakintensity, xpeakchannel, ypeakchannel,sigma, ThisRpX,ThisRpY,lastRpX=1,lastRpY=1,weighting=1,loops=20,factor=1

image xpro=realimage("",4,19,1)
image ypro=realimage("",4,19,1)

xproj.SetSelection(0,0.05*xsize,1,0.95*xsize)
for(number j=1;j<20;j++)
{
weighting=j*0.5
image peakfitx=FitGaussian(xproj, weighting, peakintensity, xpeakchannel,sigma, fwhm,ThisRpX)
setpersistentnumbernote("ElectronDiffraction Tools:temp:"+z+":deltaX:"+abs(ThisRpX)+":Peak",xPeakchannel)

xpro.SetPixel(j-1,0,xpeakchannel)
}

yproj.SetSelection(0,0.05*ysize,1,0.95*ysize)
for(number j=1;j<20;j++)
{
weighting=j*0.5
image peakfity=FitGaussian(yproj, weighting, peakintensity, ypeakchannel,sigma, fwhm,ThisRpY)
setpersistentnumbernote("ElectronDiffraction Tools:temp:"+z+":deltaY:"+abs(ThisRpY)+":Peak",yPeakchannel)

ypro.SetPixel(j-1,0,ypeakchannel)
}


xpro.GetSize(xsize,ysize)
x1=left+xpro[0,0.5*xsize,1,xsize].mean()

ypro.GetSize(xsize,ysize)
y1=top+ypro[0,0.5*xsize,1,xsize].mean()

Component box=NewBoxAnnotation(y1-0.5*sigma,x1-0.5*sigma,y1+0.5*sigma,x1+0.5*sigma)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)	
setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot"+z+":X",x1)
setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot"+z+":Y",y1)

deletepersistentnote("ElectronDiffraction Tools:temp")
}

//clear components
annots=imgdisp.componentcountchildrenoftype(5)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x0)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y0)

setNumberNote(img,"ElectronDiffraction Tools:Center:X", x0 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", y0)

Result("New center is: ("+x0+", "+y0+")\n")

ROI roi_1 = NewROI()
ROISetRectangle(roi_1, y0-5, x0-5, y0+5, x0+5 ) 
roi_1.roisetvolatile(0)
roi_1.roisetcolor(1,0,0)
imgdisp.ImageDisplayAddROI(roi_1)

setstringnote(img, "ElectronDiffraction Tools:File type","SAED")
}
}

/*************************************************************************
Profile Button: to get an intensity profile from a SAED pattern
*************************************************************************/
void ProfileButton(object self)
{
number xsize, ysize, centrex, centrey, minivalue,xscale,yscale,xwin,ywin

image frontimage:=getfrontimage()
GetSize( frontimage, xsize, ysize )
frontimage.getscale(xscale,yscale)
frontimage.GetWindowPosition(xwin,ywin)
string imgname=frontimage.getname()

Result("-----------------------------------------------------------------------------------------\nProfile button: \n1) to get an intensity profile and its linear image of the SAED pattern;\n2) ALT key: to get a stacked XRD profile;\n3) Ctrl key: to get a surface plot of a linear image.\n4) Spacebar:  to align the stacked profiles or the linear image.\n -----------------------------------------------------------------------------------------\n")

string unit
GetNumberNote(frontimage,"ElectronDiffraction Tools:Center:X", centrex )
GetNumberNote(frontimage,"ElectronDiffraction Tools:Center:Y", centrey )

if(OptionDown())  //ALT+: to get some stacked profiles
{
image img:=getfrontimage()
imagedisplay imgdisp=img.ImageGetImageDisplay(0)
string imgname=img.GetName()

number t,l,b,r,x0,y0,xsize,ysize
img.GetSize(xsize,ysize)
GetNumberNote(img,"ElectronDiffraction Tools:Center:X", x0 )
GetNumberNote(img,"ElectronDiffraction Tools:Center:Y", y0 )

SetNumberNote(img,"ElectronDiffraction Tools:Temp:Center:X", x0 )
SetNumberNote(img,"ElectronDiffraction Tools:Temp:Center:Y", y0 )

number width=Minimum(x0,xsize-x0,y0,ysize-y0)

image new=img[y0-width,x0-width,y0+width,x0+width]
number samples=2*width
number k=2*Pi()/samples

number centrey,centrex
centrey=width+0.5
centrex=width+0.5

image linearimg := RealImage( "Linear image"+imgname, 4, 2*width, samples )
linearimg = warp(new,	(icol * sin(irow * k))/2+ centrex, (icol * cos(irow * k))/2+ centrey )
linearimg=linearimg-linearimg.min()

getsize(linearimg,xsize,ysize)
number xscale=img.ImageGetDimensionScale(0)/2

string unitstring=img.ImageGetDimensionUnitString(0)
linearimg.imagesetdimensioncalibration(0,0,xscale,"s (1/nm)",0)
unitstring="deg."
imagesetdimensioncalibration( linearimg, 1, 0, 360/ysize, unitstring, 0)

//Number of profiles (note: It should be less than 100 if you wan to import in Mdi Jade.)
number y,profileno,spaceno=1,no=30
GetNumber("The number of profiles (<100)",no,no)
number integralwidth=floor(ysize/no)

SetPersistentNumberNote("ElectronDiffraction Tools:Temp:Slice no",no)
SetPersistentNumberNote("ElectronDiffraction Tools:Temp:Width",width)

//Save a total profile
image profileimage=RealImage("",4,xsize,1) 
profileimage.ShowImage()

profileimage.SetName("Profiles-"+img.GetName())
profileimage.ImageCopyCalibrationFrom(linearimg)
profileimage.setstringnote("ElectronDiffraction Tools:File type","Intensity Profile")
profileimage.SetWindowPosition(xwin,ywin)

imagedisplay Profiledisp=profileimage.ImageGetImageDisplay(0)
Profiledisp.lineplotimagedisplaysetgridon(0)  //gridon=0
Profiledisp.lineplotimagedisplaysetbackgroundon(0)  //bgon=0

//and then subprofile
for(number i=0;i<no;i++)
{
y=(i+0.5)*integralwidth

if(i==0)
{
image profile=RealImage("",4,xsize,1) 
profile[icol,0] += linearimg[0,0,integralwidth,xsize]
profile/=integralwidth

profileimage=profile
lineplotimagedisplaysetslicedrawingstyle(Profiledisp,0,1)
lineplotimagedisplaysetslicecomponentcolor(Profiledisp, 0, 0,1,0,0)

object sliceid=imgdisp.imagedisplaygetsliceidbyindex(0)
Profiledisp.imagedisplaySetslicelabelbyid(sliceid,"0")
}

else if(i==no-1)
{
y=i*integralwidth+(ysize-1-i*integralwidth)/2
integralwidth=(ysize-1-i*integralwidth)

image profile=RealImage("",4,xsize,1) 
profile[icol,0] += linearimg[y-0.5*integralwidth,0,y+0.5*integralwidth,xsize]
profile/=integralwidth

Profiledisp.imagedisplayaddimage(profile, ""+i)
}

else
{

image profile=RealImage("",4,xsize,1) 
profile[icol,0] += linearimg[y-0.5*integralwidth,0,y+0.5*integralwidth,xsize]
profile/=integralwidth

Profiledisp.imagedisplayaddimage(profile, ""+i)
}

}

// and then Y-offset
image spec := GetFrontImage()
imageDisplay disp = spec.ImageGetImageDisplay(0)

number deltaX = 0    //X-offset
number deltaY =0.2  //Y-offset
deltaY =deltaY* spec.max() 

// get current reference slice index and its ID
number refSlice_idx = disp.LinePlotImageDisplayGetSlice()
object slice_ref = disp.ImageDisplayGetSliceIDByIndex( refSlice_idx )
 
number int_offset, int_scale
number pos_offset, pos_scale
number nSlices = disp.LinePlotImageDisplayCountSlices()
for ( number i = 1; i < nSlices; i++ )
{
 object slice_src = disp.ImageDisplayGetSliceIDByIndex( i )
 
 // get current transform factors between slice and reference slice
 disp.LinePlotImageDisplayGetImageToGroupTransform( slice_src, slice_ref, int_offset, int_scale, pos_offset, pos_scale )

 pos_offset = ( i / nSlices ) * deltaX
 int_offset = ( i / nSlices ) * deltaY
 
 // set new transform factors between slice and reference slice
 disp.LinePlotImageDisplaySetImageToGroupTransform( slice_src, slice_ref, int_offset, int_scale, pos_offset, pos_scale )
}
}

else if(ControlDown())  //CTRL+ : to get a surface plot
{
number xsize, ysize, centrex, centrey, minivalue,xscale,yscale

image temp:=getfrontimage()
image frontimage:=temp.ImageClone()
GetSize( frontimage, xsize, ysize )
frontimage.getscale(xscale,yscale)
string imgname=frontimage.getname()

frontimage.ShowImage()

imagedisplay lineardisp=frontimage.ImageGetImageDisplay(0)
lineardisp.ImageDisplaySetCaptionOn(1)

lineardisp.ImageDisplayChangeDisplayType(2)

frontimage.GetSize(xsize,ysize)

number x_axis_x,  x_axis_y,  y_axis_x,  y_axis_y,  z_axis
SurfacePlotImageDisplayGetCubeAxes(lineardisp, x_axis_x,  x_axis_y,  y_axis_x,  y_axis_y,  z_axis)

x_axis_y=0
y_axis_x=0
y_axis_y=ysize/2
SurfacePlotImageDisplaySetCubeAxes(lineardisp, x_axis_x,  x_axis_y,  y_axis_x,  y_axis_y,  z_axis)

lineardisp.SurfacePlotImageDisplaySetShadingOn(1)
frontimage.setcolormode(3)
lineardisp.ImageDisplaySetCaptionOn(1)
frontimage.SetWindowPosition(xwin,ywin)
frontimage.setwindowsize(600,600*3.5/4)

setstringnote(frontimage, "ElectronDiffraction Tools:File type","Surface Plot")
}

else if(SpaceDown()) //alignment of profiles
{
image img:=getfrontimage()

number xsize, ysize,iterations=500
getsize(img, xsize, ysize)
string imgname=getname(img)

string imgtype
getstringnote(img, "ElectronDiffraction Tools:File type",imgtype)


if(imgtype=="Intensity profile")
{
// Checks that an ROI defines the match area
imagedisplay imgdisp=img.imagegetimagedisplay(0)
number norois=imgdisp.imagedisplaycountrois()
if(norois<1)
{
showalert("Ensure a ROI defines the match area.",2)
exit(0)
}

number left, right
roi thisroi=imgdisp.imagedisplaygetroi(0)
thisroi.roigetrange(left, right)

// Checks that two spectra are overlaid
number noslices=imgdisp.imagedisplaycountslices()
if(noslices<2)
{
showalert("Ensure that the two spectra to be matched are overlaid.",2)
exit(0)
}

// Source the two slices and info thereon
image back=img{0}
for(number i=1;i<noslices;i++)
{
image fore=img{i}
image shiftedfore=imageclone(fore)
number initialmisfit=1e99
number finalshift=0

// this function scales the intensity of the front spectrum to match the back
image scaledfore=MatchSpectralIntensities(fore, back, 0, 0, iterations)
img{i}=scaledfore

// Then Fit the intensity matched spectrum to the back
shiftedfore=MatchSpectralChannels(scaledfore, back, left, right)
img{i}=shiftedfore

// Now repeat the intensity matching over just the Region of interest
fore=img{i}
scaledfore=MatchSpectralIntensities(fore, back, left, right, iterations)
img{i}=scaledfore

// And finally repeat the energy matching one last time
fore=img{i}
shiftedfore=MatchSpectralChannels(fore, back, left, right)
img{i}=shiftedfore
}

result("All profiles have been aligned!\n")
}

else if(imgtype=="Linear Image")
{
// main script
image img:=getfrontimage()
image new:=img.imageclone()
new=0

image new1:=img.imageclone()
new1=0
new1.SetName("Aligned-"+img.GetName())
new.SetName("Aligned-"+img.GetName())

number xsize, ysize,iterations=150
getsize(img, xsize, ysize)
string imgname=getname(img)

image profile=realimage("",4,xsize,1)

// Checks that an ROI defines the match area
imagedisplay imgdisp=img.imagegetimagedisplay(0)
number norois=imgdisp.imagedisplaycountrois()
if(norois<1)
{
showalert("Ensure a ROI defines the match area.",2)
exit(0)
}

number t,l,b,r
img.GetSelection(t,l,b,r)

image back:= LiveProfile_ExtractLineProfile( img, 0, 0, xsize-1, 0, 1 )
new[0,0,1,xsize]=back

for(number i=1;i<ysize;i++)
{
image fore:= LiveProfile_ExtractLineProfile( img, 0, i, xsize-1, i, 1 )

image shiftedfore=imageclone(fore)

number initialmisfit=1e99
number finalshift=0

// this function scales the intensity of the front spectrum to match the back
image scaledfore=MatchSpectralIntensities(fore, back, 0, 0, iterations)
new[i,0,i+1,xsize]=scaledfore

// Then Fit the intensity matched spectrum to the back
shiftedfore=MatchSpectralChannels(scaledfore, back, l, r)
new[i,0,i+1,xsize]=shiftedfore

// Now repeat the intensity matching over just the Region of interest
fore:= LiveProfile_ExtractLineProfile( new, 0, i, xsize-1, i, 1 )
scaledfore=MatchSpectralIntensities(fore, back, l, r, iterations)
new[i,0,i+1,xsize]=scaledfore

// And finally repeat the energy matching one last time
fore:= LiveProfile_ExtractLineProfile( new, 0, i, xsize-1, i, 1 )
shiftedfore=MatchSpectralChannels(fore, back, l, r)
new[i,0,i+1,xsize]=shiftedfore
}

new1=new  //linear image
new1.ShowImage()
imagedisplay lineardisp1=new1.ImageGetImageDisplay(0)
lineardisp1.ImageDisplaySetCaptionOn(1)
setwindowposition(new1,xwin,ywin)

profile[icol,0] += new  //intensity profile
profile/=ysize
profile.ImageCopyCalibrationFrom(new)
profile.SetName("Aligned profile-"+img.GetName())
profile.ShowImage()
imagedisplay profdisp=profile.imagegetimagedisplay(0)
profdisp.lineplotimagedisplaysetgridon(0)  //gridon=0
profdisp.lineplotimagedisplaysetbackgroundon(0)  //bgon=0
setwindowposition(profile,xwin,ywin)

setstringnote(profile, "ElectronDiffraction Tools:File type","Intensity Profile")
}

}

else  //Get an intensity profile
{
number xwidth,ywidth,minwidth
xwidth=tert(centrex-0.5*xsize>0,xsize-centrex,centrex) 
ywidth=tert(centrey-0.5*ysize>0,ysize-centrey,centrey) 
minwidth=min(xwidth,ywidth)

number samples=2*minwidth
number k=2*Pi()/samples

image linearimg := RealImage( "Linear image-"+imgname, 4, 2*minwidth, samples )
linearimg = warp( frontimage[centrey-minwidth,centrex-minwidth,centrey+minwidth,centrex+minwidth],	(icol * sin(irow * k))/2+ minwidth+0.5, (icol * cos(irow * k))/2+ minwidth+0.5 )

getsize(linearimg,xsize,ysize)
xscale=frontimage.ImageGetDimensionScale(0)

string unitstring=frontimage.ImageGetDimensionUnitString(0)
linearimg.imagesetdimensioncalibration(0,0,xscale/2,"s (1/nm)",0)
unitstring="deg."
imagesetdimensioncalibration( linearimg, 1, 0, 360/ysize, unitstring, 0)
setstringnote(linearimg, "ElectronDiffraction Tools:File type","Linear Image")
linearimg.ShowImage()

imagedisplay lineardisp=linearimg.ImageGetImageDisplay(0)
lineardisp.ImageDisplaySetCaptionOn(1)
setwindowposition(linearimg,xwin,ywin)

image profile:= RealImage("Intensity profile of "+imgname,4,xsize,1) 
profile[icol,0] += linearimg
profile/=samples

profile.ShowImage()
imagecopycalibrationfrom(profile, linearimg)

xscale=linearimg.ImageGetDimensionScale(0)
profile.imagesetdimensionscale(0,xscale)
profile.imagesetdimensionunitstring(0,"s (1/nm)") 

profile.imagesetdimensionorigin(0,0)
profile.ImagesetDimensionCalibration(1,0,1,"Counts",0)
setwindowposition(profile,xwin,ywin)

imagedisplay imgdisp=profile.imagegetimagedisplay(0)
number profsizex, profsizey, maxx, maxy
getsize(profile, profsizex, profsizey)

profile=profile-profile[0,profsizex-8,1,profsizex].min()

number maxval=max(profile[0,(0.5/xscale),1, profsizex])
imgdisp.LinePlotImageDisplaySetContrastLimits( 0, maxval)  
imgdisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
imgdisp.lineplotimagedisplaysetgridon(0)  //gridon=0
imgdisp.lineplotimagedisplaysetbackgroundon(0)  //bgon=0

setstringnote(profile, "ElectronDiffraction Tools:File type","Intensity Profile")
}
}

/*************************************************************************
Peak Button: to find peaks or get a d-spacing table
*************************************************************************/
void PeaksButton(object self)
{
Result("-----------------------------------------------------------------------------------------\n1) To mark the peak defined by an ROI.\n2) Ctrl key: to locate peaks automatically;\n3) ALT key: to output the d-spacing table.\n-----------------------------------------------------------------------------------------\n")

//ALT+:  to get a d-spacing table
if(OptionDown())
{
image img:=GetFrontImage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

number xsize,ysize
img.GetSize(xsize,ysize)

number  filetype=1
GetNumber("Choose the file type:\n1=*.dsp, d-spacing table for Mdi Jade\n2=*.pks, peak profile for Search Match.",filetype,filetype)

number maxval=img[0,0.1*xsize,1,xsize].max()
if(img.max()<100)img=10000*img/maxval   

string unit,imgname
number l,r,xscale,d,wavelength
xscale=img.ImageGetDimensionScale(0)
unit=img.ImageGetDimensionUnitString(0)
imgname=img.GetName()
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Wavelength",wavelength)

number thetascale=2*180/Pi()*asin(wavelength/2*0.1*xscale)   //for copper X-ray data

number annots= imgdisp.ImageDisplayCountROIS()
for (number i=0; i<annots; i++)
{
ROI currentROI = imgdisp.ImageDisplayGetROI( i )
currentROI.roigetrange(l,r)

number twotheta
number intensity=img.GetPixel(0.5*(l+r),0)

if(unit=="s (1/nm)")
{
d=10/(0.5*(l+r)*xscale )
twotheta=0.5*(l+r)*thetascale   //for copper x-ray data
}

if(unit=="s (1/A)")
{
d=10/(0.5*(l+r)*xscale )
}

if(unit=="Q (1/A)")
{
d=2*Pi()/(0.5*(l+r)*xscale )
}

if(unit=="Two theta (deg.)")
{
d=wavelength/(2*sin(0.5*(l+r)*xscale*0.5*Pi()/180)) 
}

number x=d/100
setnumbernote(img,"ElectronDiffraction Tools:temp:"+x+":d",d)
setnumbernote(img,"ElectronDiffraction Tools:temp:"+x+":Intensity",intensity)
setnumbernote(img,"ElectronDiffraction Tools:temp:"+x+":Two theta",twotheta)
}

if(filetype==1)
{
string filename
If (!SaveAsDialog("Save d-spacing table as *.dsp (d-spacing table for Mdi Jade)", GetApplicationDirectory(2,0) +imgname+ ".dsp", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
		
TagGroup tg,tg1
tg=img.ImageGetTagGroup()
tg.taggroupgettagastaggroup("ElectronDiffraction Tools:temp", tg1)
number no=tg1.taggroupcounttags()

WriteFile(fileID, "\n"+no+"    DSPACE    "+wavelength+"    US     RIR= 0\n")

for(number i=no-1;i>-1;i--)
{
number x,intensity

x=val( tg1.TagGroupGetTaglabel( i ))
getnumbernote(img,"ElectronDiffraction Tools:temp:"+x+":d",d)
getnumbernote(img,"ElectronDiffraction Tools:temp:"+x+":Intensity",intensity)

WriteFile(fileID, d+"       "+intensity+"\n")
}

deletenote(img,"ElectronDiffraction Tools:temp")

img=img*maxval/10000  //复原

CloseFile(fileID)
}

else if(filetype==2)
{
string filename
If (!SaveAsDialog("Save d-spacing table as *.pks (peak file for Search Match)", GetApplicationDirectory(2,0) +imgname+ ".pks", filename)) Exit(0)
result(filename+">>  is saved, the wavelength is "+wavelength+" A\n")
number fileID = CreateFileForWriting(filename)
	
TagGroup tg,tg1
tg=img.ImageGetTagGroup()
tg.taggroupgettagastaggroup("ElectronDiffraction Tools:temp", tg1)
number no=tg1.taggroupcounttags()

number twotheta

for(number i=no-1;i>-1;i--)
{
number x,intensity

x=val( tg1.TagGroupGetTaglabel( i ))
getnumbernote(img,"ElectronDiffraction Tools:temp:"+x+":d",d)
getnumbernote(img,"ElectronDiffraction Tools:temp:"+x+":Intensity",intensity)
getnumbernote(img,"ElectronDiffraction Tools:temp:"+x+":Two theta",twotheta)

WriteFile(fileID, twotheta+"        "+intensity+"\n")
}

deletenote(img,"ElectronDiffraction Tools:temp")

img=img*maxval/10000  //复原

CloseFile(fileID)

}
}

//Find peaks by using top-hat filter
else if(ControlDown())
{
image  img:=getfrontimage()
ImageDisplay imgdisp = img.ImageGetImageDisplay(0)

number t,l,b,r,xsize,ysize
img.GetSize(xsize,ysize)
img.GetSelection(t,l,b,r)
number width=abs(r-l)

if(width<2||width==xsize)GetNumber("Peak width?",width,width)

img.ClearSelection()

number hatwidth=width/2, brimwidth=width/4, hatthreshold=0.0005,PeakorTrough=0  //PeakorTrough=0: peak, PeakorTrough=1: trough, PeakorTrough=2: all
tophatfilter(img, brimwidth, hatwidth, hatthreshold, PeakorTrough)
}

//Add a peak from a ROI
else
{
image  img:=getfrontimage()
ImageDisplay imgdisp = img.ImageGetImageDisplay(0)

number t,l,b,r
img.GetSelection(t,l,b,r)

roi theroi=createroi() 
theroi.roisetrange(0.5*(l+r), 0.5*(l+r))
theroi.roisetvolatile(0)
imgdisp.imagedisplayaddroi(theroi)

img.clearselection()
}
}

/*************************************************************************
Calibration Button: to calibrate patterns
*************************************************************************/
void CalibButton(object self)
{
Result("-----------------------------------------------------------------------------------------\n1) Calibrate the profile by a known peak: Drag a ROI at the desired peak; and then set the known d-spacing (A), the middle location of the ROI will be used to calibrate.\n2) Ctrl key: to calibrate SAED pattern by a known spacing between two spot or a scale bar;\n3)ALT key: to output scale parameters.\n-----------------------------------------------------------------------------------------\n")

if(OptionDown())  //Output the scale information of the selected pattern
{
image img:=GetFrontImage()
number xorigin,xscale,yorigin,yscale,xsize,ysize
string xunitstring,yunitstring

img.GetSize(xsize,ysize)
If(ysize<=1)
{
img.ImageGetDimensionCalibration(0,xorigin,xscale,xunitstring,1)
Result("\nOrigin: "+xorigin+";    Scale: "+xscale+";    Units: "+xunitstring+"\n")

Result("colint x0:=-"+xorigin*xscale+" inc:="+xscale+"\n") 
}

If(ysize>1)
{
img.ImageGetDimensionCalibration(0,xorigin,xscale,xunitstring,1)
img.ImageGetDimensionCalibration(1,yorigin,yscale,yunitstring,1)
Result("\nOrigin: ("+xorigin+", "+yorigin+");    Scale: ("+xscale+", "+yscale+");    Units: ("+xunitstring+", "+yunitstring+")\n")
}
}

//CTRL+:  to calibrate SAED pattern
if(ControlDown())
{

Number flag=1
GetNumber("Choose calibration method? (1=spots,  2=Scale bar)",flag,flag)

// calibrate SAED using two diffraction spots
if(flag==1)
{
image  img:=getfrontimage()
ImageDisplay imgdisp = img.ImageGetImageDisplay(0)
ImageDocument  imgdoc=getfrontimagedocument()
documentwindow imgwin=ImageDocumentGetWindow(imgdoc)
string imgname=img.getname()

number centrex, centrey, x0,y0,x1,y1,x2,y2,x3,y3,x,y,xsize,ysize,xscale,yscale, top,left,bottom,right,width,times
img.GetSize(xsize,ysize)

//clear all marks
number i
number annots= imgdisp.ImageDisplayCountROIS()
for (i=0; i<annots; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
 }

annots=imgdisp.componentcountchildrenoftype(13)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(2)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(2,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(5)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(6)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

number radius=xsize/80,smooth
component textannot=newtextannotation(imgdisp,xsize/15, ysize/15,"SPACE to mark the desired center!", 12)
textannot.componentsetfillmode(1)
textannot.componentsetdrawingmode(2) 
textannot.componentsetforegroundcolor(1,1,0)
textannot.componentsetfontfacename("Microsoft Sans Serif")
	
number k=1,m=1,noclick=1,twoclicks=1,n
while(2>m)
{
number keypress=getkey()

getwindowsize(img, xsize, ysize)
 getscale(img,xscale,yscale)
 
number mouse_win_x, mouse_win_y
Number img_view_t, img_view_l, img_view_b, img_view_r
Number v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y
Number img_win_t, img_win_l, img_win_b, img_win_r
Number i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y
Number i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y
Number mouse_img_x, mouse_img_y

windowgetmouseposition(imgwin, mouse_win_x, mouse_win_y)
imgdisp.ImageDisplayGetImageRectInView( img_view_t, img_view_l, img_view_b, img_view_r )
imgdoc.ImageDocumentGetViewToWindowTransform( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y )
ObjectTransformTransformRect( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                            , img_view_t, img_view_l, img_view_b, img_view_r\
                            , img_win_t, img_win_l, img_win_b, img_win_r );

imgdisp.ComponentGetChildToViewTransform( i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y )
		
ObjectTransformCompose( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                      , i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y\
                      , i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y )
		
ObjectTransformUntransformPoint( i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y\
                               , mouse_win_x, mouse_win_y, mouse_img_x, mouse_img_y );

x1=mouse_img_x
y1=mouse_img_y

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)

		
if(keypress==32) //Space bar to mark it
{
number t,l,b,r
img.GetSelection(t,l,b,r)

number smooth=10
image new=img[].Gaussfilter(smooth)

new.GetSize(xsize,ysize)
image xproj=Realimage("",4,xsize,1)
image yproj=Realimage("",4,ysize,1)

number peakposition, peakintensity, peakchannel, sigma, fwhm,weighting=3.5
xproj[icol,0]+=new
xproj=xproj-xproj.min()+1
CoarseFit(xproj, weighting, peakintensity, peakchannel,sigma, fwhm)
x1=peakchannel+l

yproj[irow,0]+=new
yproj=yproj-yproj.min()+1
CoarseFit(yproj, weighting, peakintensity, peakchannel,sigma, fwhm)
y1=t+peakchannel

Component box=NewBoxAnnotation(y1-0.5*fwhm,x1-0.5*fwhm,y1+0.5*fwhm,x1+0.5*fwhm)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)

SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":X", x1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Y", y1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Width", fwhm)

if(noclick==1)
{
radius=0.5*abs(r-l)
}


img.ClearSelection()
noclick=noclick+1
}

if(keypress==29||keypress==30) //right and up key
{
radius=radius+0.15*radius
img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}

if(keypress==28||keypress==31) //left and down key
{
radius=radius-0.15*radius

if(radius<1)radius=1

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}
	
if(keypress>0 && keypress!=32&& keypress!=28&& keypress!=29&& keypress!=30&& keypress!=31||noclick==3) // Any key except space pressed - cancel
{
img.ClearSelection()
m=2
}

}

annots=imgdisp.componentcountchildrenoftype(5)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

for(number z=1;z<3;z++)  //and then refine them by Gaussian fitting
{
number centrex,centrey
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":X", centrex )
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":Y", centrey )
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":Width", width )
width=2*width

img := GetFrontImage()
img.GetSize(xsize,ysize)

image temp=img[centrey-width,centrex-width,centrey+width,centrex+width]
top=centrey-width
left=centrex-width
temp.getsize(xsize,ysize)

number x_0,y_0
times=1
times=1/times
image new := RealImage("target",4,xsize/times,ysize/times)
new = temp.warp(icol*times,irow*times)

number loop=20,fwhm,leftpart,rightpart
image prox=realimage("",4,loop,1)
image proy=realimage("",4,loop,1)

smooth=20
new=new.Gaussfilter(smooth)

new.GetSize(xsize,ysize)
image xproj=Realimage("",4,xsize,1)
image yproj=Realimage("",4,ysize,1)

xproj[icol,0]+=new
xproj=xproj-xproj.min()+1
number xfwhm=fwhm(xproj,  xsize/2, leftpart,rightpart)
number xpeak=(leftpart+rightpart)/2

yproj[irow,0]+=new
yproj=yproj-yproj.min()+1
number yfwhm=fwhm(yproj,  ysize/2, leftpart,rightpart)
number ypeak=(leftpart+rightpart)/2

number peakintensity, xpeakchannel, ypeakchannel,sigma, ThisRpX,ThisRpY,lastRpX=1,lastRpY=1,weighting=1,loops=20,factor=1

for(number i=18;i>0;i--)
{
number t,l,b,r
t=0
b=1
factor=0.05*i

if(i==18)
{
l=xpeak-factor*xsize/2
r=xpeak+factor*xsize/2
}

else
{
l=xpeakchannel-factor*xsize/2
r=xpeakchannel+factor*xsize/2
}

if(abs(r-l)>0.9*xfwhm&&abs(r-l)<xsize)
{
xproj.SetSelection(t,l,b,r)
for(number j=1;j<20;j++)
{
weighting=j*0.5
image peakfitx=FitGaussian(xproj, weighting, peakintensity, xpeakchannel,sigma, fwhm,ThisRpX)
setpersistentnumbernote("ElectronDiffraction Tools:temp:deltaX:"+abs(ThisRpX)+":Peak",xPeakchannel)
}
}

if(i==18)
{
l=ypeak-factor*ysize/2
r=ypeak+factor*ysize/2
}

else
{
l=ypeakchannel-factor*ysize/2
r=ypeakchannel+factor*ysize/2
}

if(abs(r-l)>0.9*yfwhm&&abs(r-l)<ysize)
{
yproj.SetSelection(t,l,b,r)
for(number j=1;j<20;j++)
{
weighting=j*0.5
image peakfity=FitGaussian(yproj, weighting, peakintensity, ypeakchannel,sigma, fwhm,ThisRpY)
setpersistentnumbernote("ElectronDiffraction Tools:temp:deltaY:"+abs(ThisRpY)+":Peak",yPeakchannel)
}
}

if(i==1)
{
taggroup ptags=getpersistenttaggroup()
taggroup temptags
ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:temp:deltaX", temptags)
number value=val(temptags.TagGroupGetTaglabel(0))
getpersistentnumbernote("ElectronDiffraction Tools:temp:deltaX:"+value+":Peak",xPeakchannel) 
x1=xPeakchannel*times+left

ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:temp:deltaY", temptags)
value=val(temptags.TagGroupGetTaglabel(0))
getpersistentnumbernote("ElectronDiffraction Tools:temp:deltaY:"+value+":Peak",yPeakchannel) 
y1=yPeakchannel*times+top

setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot"+z+":X",x1)
setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot"+z+":Y",y1)

Component box=NewBoxAnnotation(y1-0.5*fwhm,x1-0.5*fwhm,y1+0.5*fwhm,x1+0.5*fwhm)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)	
}
}

deletepersistentnote("ElectronDiffraction Tools:temp")
}

Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)

// Lines are drawn between the points
component lineannotation=newlineannotation(y1,x1,y2,x2)
lineannotation.componentsetfillmode(2) 
lineannotation.componentsetdrawingmode(2)
lineannotation.componentsetforegroundcolor(1,0,1)
imgdisp.componentaddchildatend(lineannotation) 

number order=1
getnumber("Divided by?",order,order)  //order of spots

img.GetScale(xscale,yscale)
result("Befor calibration: "+xscale+" (1/nm)\n")

number unitflag=1,dspacing=1
GetNumber("Spacing of the drawn line (A or 1/A)",dspacing,dspacing)
Getnumber("Choose the unit (1=A, 2=1/A)", unitflag,unitflag)
number r=sqrt((x2-x1)**2+(y2-y1)**2)/order

//calibrate it
if(unitflag==1)xscale=10/(r*dspacing)
if(unitflag==2)xscale=10*dspacing/r

img.setScale(xscale,xscale)
img.SetUnitString("1/nm")
result("After calibration: "+xscale+" (1/nm)\n")

deletenote(img,"ElectronDiffraction Tools:Live profile") 
}

//calibrate SAED using scale bar
if(flag==2)  
{
image  img:=getfrontimage()
ImageDisplay imgdisp = img.ImageGetImageDisplay(0)
ImageDocument  imgdoc=getfrontimagedocument()
documentwindow imgwin=ImageDocumentGetWindow(imgdoc)
string imgname=img.getname()

number centrex, centrey, x0,y0,x1,y1,x2,y2,x3,y3,x,y,xsize,ysize,xscale,yscale, top,left,bottom,right,width,times
img.GetSize(xsize,ysize)

//clear all marks, and measure two end points of the scale bars
number i
number annots= imgdisp.ImageDisplayCountROIS()
for (i=0; i<annots; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
 }

annots=imgdisp.componentcountchildrenoftype(13)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(2)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(2,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(5)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(6)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

number radius=2,smooth
component textannot=newtextannotation(imgdisp,xsize/15, ysize/15,"SPACE to mark the desired center!", 12)
textannot.componentsetfillmode(1)
textannot.componentsetdrawingmode(2) 
textannot.componentsetforegroundcolor(1,1,0)
textannot.componentsetfontfacename("Microsoft Sans Serif")

number k=1,m=1,noclick=1,twoclicks=1,n
while(2>m)
{
number keypress=getkey()

getwindowsize(img, xsize, ysize)
 getscale(img,xscale,yscale)
 
number mouse_win_x, mouse_win_y
Number img_view_t, img_view_l, img_view_b, img_view_r
Number v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y
Number img_win_t, img_win_l, img_win_b, img_win_r
Number i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y
Number i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y
Number mouse_img_x, mouse_img_y

windowgetmouseposition(imgwin, mouse_win_x, mouse_win_y)
imgdisp.ImageDisplayGetImageRectInView( img_view_t, img_view_l, img_view_b, img_view_r )
imgdoc.ImageDocumentGetViewToWindowTransform( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y )
ObjectTransformTransformRect( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                            , img_view_t, img_view_l, img_view_b, img_view_r\
                            , img_win_t, img_win_l, img_win_b, img_win_r );

imgdisp.ComponentGetChildToViewTransform( i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y )
		
ObjectTransformCompose( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                      , i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y\
                      , i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y )
		
ObjectTransformUntransformPoint( i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y\
                               , mouse_win_x, mouse_win_y, mouse_img_x, mouse_img_y );

x1=mouse_img_x
y1=mouse_img_y

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)

		
if(keypress==32) //Space bar to mark it
{
Component box=NewBoxAnnotation(y1-radius,x1-radius,y1+radius,x1+radius)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)

SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":X", x1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Y", y1 )

img.ClearSelection()
noclick=noclick+1
}

if(keypress==29||keypress==30) //right and up key
{
radius=radius+0.15*radius
img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}

if(keypress==28||keypress==31) //left and down key
{
radius=radius-0.15*radius

if(radius<1)radius=1

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}
	
if(keypress>0 && keypress!=32&& keypress!=28&& keypress!=29&& keypress!=30&& keypress!=31||noclick==3) // Any key except space pressed - cancel
{
img.ClearSelection()
m=2
}
}

//begin to calibrate it
Getnumbernote(img,"ElectronDiffraction Tools:Live profile:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Live profile:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Live profile:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Live profile:Spot2:Y",y2)

// Lines are drawn between the points
component lineannotation=newlineannotation(y1,x1,y2,x2)
lineannotation.componentsetfillmode(2) 
lineannotation.componentsetdrawingmode(2)
lineannotation.componentsetforegroundcolor(1,0,1)
imgdisp.componentaddchildatend(lineannotation) 

img.GetScale(xscale,yscale)
result("Befor calibration: "+xscale+" (1/nm)\n")

number unitflag=1,dspacing=5
GetNumber("Spacing of the drawn line (1/nm)",dspacing,dspacing)

number r=sqrt((x2-x1)**2+(y2-y1)**2)
xscale=dspacing/r

img.setScale(xscale,xscale)
img.SetUnitString("1/nm")
result("After calibration: "+xscale+" (1/nm)\n")

deletenote(img,"ElectronDiffraction Tools:Live profile") 
}

}


else
{
image img:=getfrontimage()
Number xsize, ysize,xscale,yscale,t,l,b,r
img.getsize(xsize,ysize)
img.getselection(t,l,b,r)
xscale=img.imagegetdimensionscale(0)

string unit
unit=img.ImageGetDimensionUnitString(0)

if(unit=="s (1/nm)")
{
result("xscale (before) is: "+xscale+" 1/nm"+"\n")

number d=1
getnumber("d-spacing (A) of the peak",d,d)
xscale=10/(0.5*(l+r)*d)
img.imagesetdimensionscale(0,xscale)

img.SetSelection(0,0.5*(l+r),1,0.5*(l+r))

result("xscale (after) is: "+xscale+" 1/nm"+"\n")
}

if(unit=="Q (1/A)")
{
result("xscale (before) is: "+xscale+" 1/A"+"\n")

number d=1
getnumber("d-spacing (A) of the peak",d,d)
xscale=2*Pi()/(0.5*(l+r)*d)
img.imagesetdimensionscale(0,xscale)

img.SetSelection(0,0.5*(l+r),1,0.5*(l+r))

result("xscale (after) is: "+xscale+" 1/A"+"\n")
}
}
}

/*************************************************************************
Save Button
*************************************************************************/
void SaveButton(object self)
{
Result("-----------------------------------------------------------------------------------------\n1) Save profiles as a text file;\n2) CTRL key: Save Linear image or Surface plot as *.dif or *.txt;\n3) Alt key: Save SAED pattern as *.dif or *.txt\n-----------------------------------------------------------------------------------------\n")

//Save a linear image
if(ControlDown())
{
number xsize, ysize, centrex, centrey, minivalue,xscale,yscale

image frontimage:=getfrontimage()
GetSize( frontimage, xsize, ysize )
frontimage.getscale(xscale,yscale)
string imgname=frontimage.getname()

number  savetype=1
GetNumber("Choose the file type:\n1=*.dif (for Mdi Jade)\n2=*.txt (for the image plot in Origin)",savetype,savetype)

string filetype
getstringnote(frontimage, "ElectronDiffraction Tools:File type",filetype)

if(filetype=="Linear Image"||filetype=="Surface Plot")
{
number wavelength,no=30
GetNumber("The number of profiles (<500)",no,no)
number width=floor(ysize/no)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Wavelength",wavelength)

if(savetype==1)   //dif file
{
string filename
If (!SaveAsDialog("Save Linear Image as *.dif file (for Mdi Jade)", GetApplicationDirectory(2,0) +imgname+ ".dif", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
number yval, Q,y,profileno
string datetime=getdate(1)+" "+gettime(1)
xscale=2*180/Pi()*asin(wavelength/2*0.1*xscale)

image temp=frontimage-frontimage.min()

number theta1, theta2
theta1=0*xscale
theta2=xsize*xscale

//save a total profile
image profile=RealImage("",4,xsize,1) 
profile[icol,0] += temp
profile/=ysize

WriteFile(fileID, datetime+"  DIF"+" "+imgname+"\n")
WriteFile(fileID, ""+theta1+"    "+xscale+"    1    US    "+wavelength+"    "+theta2+"    "+xsize+"    <Profile 0"+">\n")
for(number i=0; i<xsize; i++)
	{
		yval=getpixel(profile, i, 0)
		WriteFile(fileID, ""+yval+"    ")
		
		if(i==xsize-1)
		{
		WriteFile(fileID, "\n\n")		
		}		
	}

//and then save subprofiles
for(number i=0;i<no;i++)
{
y=(i+0.5)*width

if(i==no-1)
{
y=i*width+(ysize-1-i*width)/2
width=(ysize-1-i*width)

profile := LiveProfile_ExtractLineProfile(temp,0,y,xsize,y,width)
}

else
{
profile := LiveProfile_ExtractLineProfile(temp,0,y,xsize,y,width)
}

profileno=i+1
WriteFile(fileID, ""+theta1+"    "+xscale+"    1    US    "+wavelength+"    "+theta2+"    "+xsize+" <Profile "+profileno+">\n")
for(number i=0; i<xsize; i++)
	{
		yval=getpixel(profile, i, 0)
		WriteFile(fileID, ""+yval+"    ")
		
		if(i==xsize-1)
		{
		WriteFile(fileID, "\n\n")		
		}		
	}	
}

CloseFile(fileID)	
}


if(savetype==2)   //xyy file
{
string filename
If (!SaveAsDialog("Save Linear image as XYY data (for the image plot in Origin)", GetApplicationDirectory(2,0) +imgname+ ".txt", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
number yval, Q,y,profileno
string datetime=getdate(1)+" "+gettime(1)
xscale=2*180/Pi()*asin(wavelength/2*0.1*xscale)

image temp=frontimage-frontimage.min()

number theta,m
for(number i=0;i<xsize;i++)
{
theta=i*xscale		
WriteFile(fileID, theta+"    ")

	for(number j=0;j<no;j++)
	{
	m=j+1
	yval=temp[j*width,i,m*width,i+1].sum()/width
	WriteFile(fileID, yval+"    ")		
	
	if(m==no)
		{
		WriteFile(fileID, "\n")		
		}
	}
	
}
CloseFile(fileID)	
}
}
}

//Save SAED to stacked profiles
else if(OptionDown())
{
number xsize, ysize, centrex, centrey, minivalue,xscale,yscale

image frontimage:=getfrontimage()
GetSize( frontimage, xsize, ysize )
frontimage.getscale(xscale,yscale)
string imgname=frontimage.getname()

number  savetype=1
GetNumber("Choose the file type:\n1=*.dif (for Mdi Jade)\n2=*.txt (for the image plot in Origin)",savetype,savetype)

string filetype
getstringnote(frontimage, "ElectronDiffraction Tools:File type",filetype)

if(filetype=="SAED")
{
number wavelength,no=30
GetNumber("The number of profiles (<500)",no,no)
number width=floor(ysize/no)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Wavelength",wavelength)

GetNumberNote(frontimage,"ElectronDiffraction Tools:Center:X", centrex )
GetNumberNote(frontimage,"ElectronDiffraction Tools:Center:Y", centrey )

number xwidth,ywidth,minwidth
xwidth=tert(centrex-0.5*xsize>0,xsize-centrex,centrex) 
ywidth=tert(centrey-0.5*ysize>0,ysize-centrey,centrey) 
minwidth=min(xwidth,ywidth)

number samples=2*minwidth
number k=2*Pi()/samples

image linearimg := RealImage( "Linear image-"+imgname, 4, 2*minwidth, samples )
linearimg = warp( frontimage[centrey-minwidth,centrex-minwidth,centrey+minwidth,centrex+minwidth],	(icol * sin(irow * k))/2+ minwidth+0.5, (icol * cos(irow * k))/2+ minwidth+0.5 )
linearimg=linearimg-linearimg.min()

getsize(linearimg,xsize,ysize)
xscale=frontimage.ImageGetDimensionScale(0)/2

if(savetype==1)   //dif file
{
string filename
If (!SaveAsDialog("Save SAED pattern as DIF file (for Mdi Jade)", GetApplicationDirectory(2,0) +imgname+ ".dif", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
number yval, Q,y,profileno
string datetime=getdate(1)+" "+gettime(1)
xscale=2*180/Pi()*asin(wavelength/2*0.1*xscale)

image temp=linearimg-linearimg.min()

number theta1, theta2
theta1=0*xscale
theta2=xsize*xscale

image profile=RealImage("",4,xsize,1) 
profile[icol,0] += temp
profile/=ysize

WriteFile(fileID, datetime+"  DIF"+" "+imgname+"\n")
WriteFile(fileID, ""+theta1+"    "+xscale+"    1    US    "+wavelength+"    "+theta2+"    "+xsize+"    <Profile 0"+">\n")
for(number i=0; i<xsize; i++)
	{
		yval=getpixel(profile, i, 0)
		WriteFile(fileID, ""+yval+"    ")
		
		if(i==xsize-1)
		{
		WriteFile(fileID, "\n\n")		
		}		
	}

for(number i=0;i<no;i++)
{
y=(i+0.5)*width

if(i==no-1)
{
y=i*width+(ysize-1-i*width)/2
width=(ysize-1-i*width)

profile := LiveProfile_ExtractLineProfile(temp,0,y,xsize,y,width)
}

else
{
profile := LiveProfile_ExtractLineProfile(temp,0,y,xsize,y,width)
}

profileno=i+1
WriteFile(fileID, ""+theta1+"    "+xscale+"    1    US    "+wavelength+"    "+theta2+"    "+xsize+" <Profile "+profileno+">\n")
for(number i=0; i<xsize; i++)
	{
		yval=getpixel(profile, i, 0)
		WriteFile(fileID, ""+yval+"    ")
		
		if(i==xsize-1)
		{
		WriteFile(fileID, "\n\n")		
		}		
	}	
}

CloseFile(fileID)	
}


if(savetype==2)   //xyy file
{
string filename
If (!SaveAsDialog("Save SAED pattern as txt file (for the image plot in Origin)", GetApplicationDirectory(2,0) +imgname+ ".txt", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
number yval, Q,y,profileno
string datetime=getdate(1)+" "+gettime(1)
xscale=2*180/Pi()*asin(wavelength/2*0.1*xscale)

image temp=linearimg-linearimg.min()

number theta,m
for(number i=0;i<xsize;i++)
{
theta=i*xscale		
WriteFile(fileID, theta+"    ")

	for(number j=0;j<no;j++)
	{
	m=j+1
	yval=temp[j*width,i,m*width,i+1].sum()/width
	WriteFile(fileID, yval+"    ")		
	
	if(m==no)
		{
		WriteFile(fileID, "\n")		
		}
	}
	
}
CloseFile(fileID)	
}
}

}

//save as xyy data
else
{
image img:=getfrontimage()
imagedisplay imgdisp=img.ImageGetImageDisplay(0)

number xsize,ysize,xscale,yscale
img.GetSize(xsize, ysize )
img.getscale(xscale,yscale)
string imgname=img.getname()
string unitstring=img.ImageGetDimensionUnitString(0)

number  savetype=1
GetNumber("Choose the file type:\n1=*.dif (for Mdi Jade)\n2=*.txt (for the image plot in Origin)",savetype,savetype)

if(savetype==2)
{
string filename
If (!SaveAsDialog("Save profiles as XYY data (for the image plot in Origin)", GetApplicationDirectory(2,0) +imgname+ ".txt", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
Result("The exported profiles:  "+imgname+",     X is "+unitstring+"\n")

number noslices=imgdisp.imagedisplaycountslices()

number Q,yval,profileno=0
for(number j=0; j<xsize; j++)  
	{
		Q=j*xscale		
		WriteFile(fileID, Q+"    ")	
		
		for(number i=0;i<noslices;i++) 
			{
				yval=img{i}.GetPixel(j,0)
				WriteFile(fileID, ""+yval+"    ")	
				
				if(i==noslices-1)
					{
					WriteFile(fileID, "\n")		
					}	
			}		
	}
	
CloseFile(fileID)
}

else if(savetype==1)
{
xscale=2*180/Pi()*asin(wavelength/2*0.1*xscale)
number theta1, theta2,wavelength
theta1=0*xscale
theta2=xsize*xscale
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Wavelength",wavelength)

string datetime=getdate(1)+" "+gettime(1)

string filename
If (!SaveAsDialog("Save profiles as DIF data (data for Mdi Jade)", GetApplicationDirectory(2,0) +imgname+ ".dif", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
Result("The exported profiles:  "+imgname+",     X is "+unitstring+"\n")

WriteFile(fileID, datetime+"  DIF"+" "+imgname+"\n")

number noslices=imgdisp.imagedisplaycountslices()
for(number j=0;j<noslices;j++)
{
image profile=imgdisp{j}
WriteFile(fileID, ""+theta1+"    "+xscale+"    1    US    "+wavelength+"    "+theta2+"    "+xsize+" <Profile "+j+">\n")
for(number i=0; i<xsize; i++)
	{
		number yval=getpixel(profile, i, 0)
		WriteFile(fileID, ""+yval+"    ")
		
		if(i==xsize-1)
		{
		WriteFile(fileID, "\n\n")		
		}		
	}	
}
CloseFile(fileID)	
}
}
}

/*************************************************************************
d_Angle Button: to measure d-spacing and angle
*************************************************************************/
void d_Angle(object self)
{
Result("-----------------------------------------------------------------------------------------\n1) Automatically locate 3 spots by the Space bar.\n2) Ctrl key: to automatically locate 2 spots.\n3) Alt key: to manually find 3 spots.\n-----------------------------------------------------------------------------------------\n")
if(ControlDown())
{
//clear lists and tags
	TagGroup TgList1= self.LookUpElement("hklList1")
	TagGroup TgList2= self.LookUpElement("hklList2")
	TagGroup uvwList= self.LookUpElement("uvwList")
	uvwList.DLGRemoveElement(0) 
	
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags))
	{
	number tagno=dspacingtags.taggroupcounttags()

	for(number i=0;i<200;i++)
	{
	TgList1.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	}
	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}	
	

image  img:=getfrontimage()
ImageDisplay imgdisp = img.ImageGetImageDisplay(0)
ImageDocument  imgdoc=getfrontimagedocument()
documentwindow imgwin=ImageDocumentGetWindow(imgdoc)
string imgname=img.getname()

number centrex, centrey, x0,y0,x1,y1,x2,y2,x3,y3,x,y,xsize,ysize,xscale,yscale, top,left,bottom,right,width,times
img.GetSize(xsize,ysize)

//clear all marks
number i
number annots= imgdisp.ImageDisplayCountROIS()
for (i=0; i<annots; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
 }

annots=imgdisp.componentcountchildrenoftype(13)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(2)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(2,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(5)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(6)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

number radius=xsize/150,smooth

if(radius<1)radius=2

component textannot=newtextannotation(imgdisp,xsize/15, ysize/15,"SPACE to mark the desired center!", 12)
textannot.componentsetfillmode(1)
textannot.componentsetdrawingmode(2) 
textannot.componentsetforegroundcolor(1,1,0)
textannot.componentsetfontfacename("Microsoft Sans Serif")
	
number k=1,m=1,noclick=1,twoclicks=1,n
while(2>m)
{
number keypress=getkey()

getwindowsize(img, xsize, ysize)
 getscale(img,xscale,yscale)
 
number mouse_win_x, mouse_win_y
Number img_view_t, img_view_l, img_view_b, img_view_r
Number v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y
Number img_win_t, img_win_l, img_win_b, img_win_r
Number i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y
Number i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y
Number mouse_img_x, mouse_img_y

windowgetmouseposition(imgwin, mouse_win_x, mouse_win_y)
imgdisp.ImageDisplayGetImageRectInView( img_view_t, img_view_l, img_view_b, img_view_r )
imgdoc.ImageDocumentGetViewToWindowTransform( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y )
ObjectTransformTransformRect( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                            , img_view_t, img_view_l, img_view_b, img_view_r\
                            , img_win_t, img_win_l, img_win_b, img_win_r );

imgdisp.ComponentGetChildToViewTransform( i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y )
		
ObjectTransformCompose( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                      , i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y\
                      , i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y )
		
ObjectTransformUntransformPoint( i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y\
                               , mouse_win_x, mouse_win_y, mouse_img_x, mouse_img_y );

x1=mouse_img_x
y1=mouse_img_y

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)

		
if(keypress==32) //Space bar to mark it
{
number t,l,b,r,leftpart,rightpart
img.GetSize(xsize,ysize)

Image profile_h:= LiveProfile_ExtractLineProfile( img, 0, y1, xsize, y1, 4 )  //horizental projection
Image profile_v:= LiveProfile_ExtractLineProfile( img,x1,0,x1,ysize , 4 )  //vertical projection

number width_h=fwhm(profile_h,x1,leftpart,rightpart)
x1=(leftpart+rightpart)/2

number width_v=fwhm(profile_v,y1,leftpart,rightpart)
y1=(leftpart+rightpart)/2

//radius=1*(width_h+width_v)/2
Component box=NewBoxAnnotation(y1-0.5*radius,x1-0.5*radius,y1+0.5*radius,x1+0.5*radius)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)

SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":X", x1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Y", y1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Width", radius=1*(width_h+width_v)/2)
                                                                                 
if(noclick==1)                                                           
{                                                                              
//radius=1*radius
}

else if(noclick==2)
{
//radius=1*radius
}

img.ClearSelection()
noclick=noclick+1
}

if(keypress==29||keypress==30) //right and up key
{
radius=radius+0.2*radius
img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}

if(keypress==28||keypress==31) //left and down key
{
radius=radius-0.2*radius

if(radius<1)radius=1

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}
	
if(keypress>0 && keypress!=32&& keypress!=28&& keypress!=29&& keypress!=30&& keypress!=31||noclick==3) // Any key except space pressed - cancel
{
img.ClearSelection()
m=2
}

}

annots=imgdisp.componentcountchildrenoftype(5)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

for(number z=1;z<3;z++)  //and then refine these centers by Gaussian fitting
{
number centrex,centrey
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":X", centrex )
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":Y", centrey )
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":Width", width )
radius=0.6*width
width=2*width

img := GetFrontImage()
img.GetSize(xsize,ysize)

image temp=img[centrey-width,centrex-width,centrey+width,centrex+width]
top=centrey-width
left=centrex-width
temp.getsize(xsize,ysize)
temp=( ( icol -(xsize+1)/2 )**2 + ( irow - (ysize+1)/2 )**2 < radius**2 ) ? temp : 0 

number x_0,y_0
times=1
times=1/times
image new := RealImage("target",4,xsize/times,ysize/times)
new = temp.warp(icol*times,irow*times)

number fwhm,leftpart,rightpart

new=new.Gaussfilter(20)

new.GetSize(xsize,ysize)

image xproj=realimage("",4,xsize,1)
image yproj=realimage("",4,xsize,1)

xproj=LiveProfile_ExtractLineProfile( new, 0, (ysize+1)/2 , xsize-1, (ysize+1)/2 ,2 )
yproj= LiveProfile_ExtractLineProfile( new,(xsize+1)/2,0,(xsize+1)/2,ysize-1 , 2 ) 

xproj=xproj-xproj.min()
xproj=10000*xproj/xproj.max()+1

yproj=yproj-yproj.min()
yproj=10000*yproj/yproj.max()+1

number xfwhm=fwhm(xproj,  xsize/2, leftpart,rightpart)
number xpeak=(leftpart+rightpart)/2

number yfwhm=fwhm(yproj,  ysize/2, leftpart,rightpart)
number ypeak=(leftpart+rightpart)/2

number peakintensity, xpeakchannel, ypeakchannel,sigma, ThisRpX,ThisRpY,lastRpX=1,lastRpY=1,weighting=1,loops=20,factor=1

image xpro=realimage("",4,19,1)
image ypro=realimage("",4,19,1)

xproj.SetSelection(0,0.05*xsize,1,0.95*xsize)
for(number j=1;j<20;j++)
{
weighting=j*0.5
image peakfitx=FitGaussian(xproj, weighting, peakintensity, xpeakchannel,sigma, fwhm,ThisRpX)
setpersistentnumbernote("ElectronDiffraction Tools:temp:"+z+":deltaX:"+abs(ThisRpX)+":Peak",xPeakchannel)

xpro.SetPixel(j-1,0,xpeakchannel)
}

yproj.SetSelection(0,0.05*ysize,1,0.95*ysize)
for(number j=1;j<20;j++)
{
weighting=j*0.5
image peakfity=FitGaussian(yproj, weighting, peakintensity, ypeakchannel,sigma, fwhm,ThisRpY)
setpersistentnumbernote("ElectronDiffraction Tools:temp:"+z+":deltaY:"+abs(ThisRpY)+":Peak",yPeakchannel)

ypro.SetPixel(j-1,0,ypeakchannel)
}


xpro.GetSize(xsize,ysize)
x1=left+xpro[0,0.5*xsize,1,xsize].mean()

ypro.GetSize(xsize,ysize)
y1=top+ypro[0,0.5*xsize,1,xsize].mean()

Component box=NewBoxAnnotation(y1-0.5*sigma,x1-0.5*sigma,y1+0.5*sigma,x1+0.5*sigma)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)	
setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot"+z+":X",x1)
setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot"+z+":Y",y1)

deletepersistentnote("ElectronDiffraction Tools:temp")
}

Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)
GetNumberNote(img,"ElectronDiffraction Tools:Center:X", x0 )
GetNumberNote(img,"ElectronDiffraction Tools:Center:Y", y0)

x3=x2
y3=y2

x2=x0
y2=y0

string order="1,1"
getstring("Orders of the diffraction spots (m,n)",order,order)


number Bytes=1
width=len(order)

for(number i=0;i<width;i++)
{
string thischar=mid(order,i,1)
if(thischar==",")Bytes=len(left(order,i))
}

m=val(left(order,Bytes))
n=val(right(order,width-Bytes-1))

x1=x2+(x1-x2)/m
y1=y2+(y1-y2)/m

x3=x2+(x3-x2)/n
y3=y2+(y3-y2)/n

Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)
Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

annots=imgdisp.componentcountchildrenoftype(5)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

// Lines are drawn between the points
component lineannotation=newlineannotation(y1,x1,y2,x2)
lineannotation.componentsetfillmode(2) 
lineannotation.componentsetdrawingmode(2)
lineannotation.componentsetforegroundcolor(1,0,1)
imgdisp.componentaddchildatend(lineannotation) 

lineannotation=newlineannotation(y2,x2,y3,x3)
lineannotation.componentsetfillmode(2) 
lineannotation.componentsetdrawingmode(2)
lineannotation.componentsetforegroundcolor(1,0,1)
imgdisp.componentaddchildatend(lineannotation) 

// Computation of the angle using the sine rule
number dxA, dxB, dxC, dyA, dyB, dyC, hypA, hypB, hypC
dxA=x2-x1
dyA=y2-y1
hypA=sqrt(dxA**2+dyA**2)

dxB=x2-x3
dyB=y2-y3
hypB=sqrt(dxB**2+dyB**2)

dxC=x1-x3
dyC=y1-y3
hypC=sqrt(dxC**2+dyC**2)

number gamma=(hypa**2+hypB**2-hypC**2)/(2*hypA*hypB)
number totalangle=acos(gamma)*(180/pi())
totalangle=round(totalangle*10)/10
string outputstring=""+totalangle+" deg."

// Data is output to the results window and onto the image
img.GetScale(xscale,yscale)
result("R1= "+hypA*xscale+" nm-1,   d1= "+10/(hypA*xscale)+" A"+"\n")
result("R2= "+hypB*xscale+" nm-1,   d2= "+10/(hypB*xscale)+" A"+"\n")
result("Angle= "+outputstring+"\n\n")

self.lookupElement("d1field").dlgtitle(" "+format(10/(hypA*xscale), "%2.4f" ))
self.lookupElement("d2field").dlgtitle(" "+format(10/(hypB*xscale), "%2.4f" ))
self.lookupElement("Anglefield").dlgtitle(" "+format(totalangle, "%6.2f" ))

//Save tags
Setnumbernote(img,"ElectronDiffraction Tools:d-Angle:d1",10/(hypA*xscale))
Setnumbernote(img,"ElectronDiffraction Tools:d-Angle:d2",10/(hypB*xscale))
Setnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",totalangle)

SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot1:X",x1)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot1:Y",y1)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot2:X",x2)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot2:Y",y2)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot3:X",x3)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot3:Y",y3)

SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Scale:Xscale",xscale)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:d-Angle:d1",10/(hypA*xscale))
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:d-Angle:d2",10/(hypB*xscale))
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:d-Angle:Angle",totalangle)

deletenote(img,"ElectronDiffraction Tools:Live profile") 
updateimage(img)
}

else if(OptionDown())
{
//clear lists and tags
	TagGroup TgList1= self.LookUpElement("hklList1") 
	TagGroup TgList2= self.LookUpElement("hklList2") 
	TagGroup uvwList= self.LookUpElement("uvwList")
	uvwList.DLGRemoveElement(0)  
	
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags))
	{
	number tagno=dspacingtags.taggroupcounttags()

	for(number i=0;i<200;i++)
	{
	TgList1.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	}
	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}	
	
image  img:=getfrontimage()
ImageDisplay imgdisp =img.ImageGetImageDisplay(0)
string imgname=img.getname()

number i
number annots= imgdisp.ImageDisplayCountROIS()
for (i=0; i<annots; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
 }

annots=imgdisp.componentcountchildrenoftype(13)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(2)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(2,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(5)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(6)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

if(img.ImageGetDataType()==3)
{
image magimg=sqrt(real(img)**2+imaginary(img)**2)
magimg.showimage()
magimg.ImageCopyCalibrationFrom(img)

number xsize,ysize
magimg.getsize(xsize,ysize)
setNumberNote(magimg,"ElectronDiffraction Tools:Center:X", xsize/2+0.5 )
setNumberNote(magimg,"ElectronDiffraction Tools:Center:Y", ysize/2+0.5 )
}

image   frontimage:=getfrontimage()
ImageDisplay frontimage_disp = frontimage.ImageGetImageDisplay(0)
ImageDocument  frontimage_doc=getfrontimagedocument()
documentwindow frontimage_win=ImageDocumentGetWindow(frontimage_doc)
number x0, y0, x1,y1,x,y,xsize,ysize
frontimage.getsize(xsize,ysize)

number k=1,m=1,n=1,noclick=1,radius=xsize/100,points=3
while(2>m)
{
number keypress=getkey()

 //find the coordinator of the spot
 number xsize, ysize,xscale,yscale,x, y,x0,y0
 getwindowsize(img, xsize, ysize)
 getscale(img,xscale,yscale)
 
// you can use mouse click to get the end point
number n, mouse_win_x, mouse_win_y
Number img_view_t, img_view_l, img_view_b, img_view_r
Number v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y
Number img_win_t, img_win_l, img_win_b, img_win_r
Number i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y
Number i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y
Number mouse_img_x, mouse_img_y

 windowgetmouseposition(frontimage_win, mouse_win_x, mouse_win_y)
 
//display how to operate
imagedisplay imgdisp=frontimage.imagegetimagedisplay(0)
number annotationID

		frontimage_disp.ImageDisplayGetImageRectInView( img_view_t, img_view_l, img_view_b, img_view_r )

		frontimage_doc.ImageDocumentGetViewToWindowTransform( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y )
		
		ObjectTransformTransformRect( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                            , img_view_t, img_view_l, img_view_b, img_view_r\
                            , img_win_t, img_win_l, img_win_b, img_win_r );

		frontimage_disp.ComponentGetChildToViewTransform( i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y )
		
		ObjectTransformCompose( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                      , i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y\
                      , i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y )
		
		ObjectTransformUntransformPoint( i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y\
                               , mouse_win_x, mouse_win_y, mouse_img_x, mouse_img_y );

		x1=mouse_img_x
		y1=mouse_img_y

frontimage.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
		
if(keypress==32) // space bar to mark the centers
{
Component box=NewBoxAnnotation(y1-radius,x1-radius,y1+radius,x1+radius)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
frontimage_disp.componentaddchildatend(box)
		
SetNumberNote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot"+noclick+":X", x1 )
SetNumberNote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot"+noclick+":Y", y1 )

frontimage.ClearSelection()

noclick=noclick+1
}	

if(keypress==29||keypress==30) //right and up key
{
radius=radius+0.1*radius
frontimage.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}

if(keypress==28||keypress==31) //left and down key
{
radius=radius-0.1*radius

if(radius<1)radius=1

frontimage.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}
	
if(keypress>0 && keypress!=32&& keypress!=28&& keypress!=29&& keypress!=30&& keypress!=31||noclick==4) // Any key except space pressed - cancel
{
frontimage.ClearSelection()
m=2
}
}

number x2,y2,x3,y3,xscale,yscale
frontimage.GetScale(xscale,yscale)

Getnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)
Getnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Getnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

string order="1,1"
getstring("Orders of the diffraction spots (m,n)",order,order)


number Bytes=1
number width=len(order)

for(number i=0;i<width;i++)
{
string thischar=mid(order,i,1)
if(thischar==",")Bytes=len(left(order,i))
}

m=val(left(order,Bytes))
n=val(right(order,width-Bytes-1))

x1=x2+(x1-x2)/m
y1=y2+(y1-y2)/m

x3=x2+(x3-x2)/n
y3=y2+(y3-y2)/n

Setnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Setnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Setnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Setnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

annots=imgdisp.componentcountchildrenoftype(5)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

// Lines are drawn between the points
component lineannotation=newlineannotation(y1,x1,y2,x2)
lineannotation.componentsetfillmode(2) 
lineannotation.componentsetdrawingmode(2)
lineannotation.componentsetforegroundcolor(1,0,1)
imgdisp.componentaddchildatend(lineannotation) 

lineannotation=newlineannotation(y2,x2,y3,x3)
lineannotation.componentsetfillmode(2) 
lineannotation.componentsetdrawingmode(2)
lineannotation.componentsetforegroundcolor(1,0,1)
imgdisp.componentaddchildatend(lineannotation) 

// Computation of the angle using the sine rule
number dxA, dxB, dxC, dyA, dyB, dyC, hypA, hypB, hypC
dxA=x2-x1
dyA=y2-y1
hypA=sqrt(dxA**2+dyA**2)

dxB=x2-x3
dyB=y2-y3
hypB=sqrt(dxB**2+dyB**2)

dxC=x1-x3
dyC=y1-y3
hypC=sqrt(dxC**2+dyC**2)

number gamma=(hypa**2+hypB**2-hypC**2)/(2*hypA*hypB)
number totalangle=acos(gamma)*(180/pi())
totalangle=round(totalangle*10)/10
string outputstring=""+totalangle+" deg."

// Data is output to the results window and onto the image
result("R1= "+hypA*xscale+" nm-1,   d1= "+10/(hypA*xscale)+" A"+"\n")
result("R2= "+hypB*xscale+" nm-1,   d2= "+10/(hypB*xscale)+" A"+"\n")
result("Angle= "+outputstring+"\n\n")

self.lookupElement("d1field").dlgtitle(" "+format(10/(hypA*xscale), "%2.4f" ))
self.lookupElement("d2field").dlgtitle(" "+format(10/(hypB*xscale), "%2.4f" ))
self.lookupElement("Anglefield").dlgtitle(" "+format(totalangle, "%6.2f" ))

//save tags
Setnumbernote(frontimage,"ElectronDiffraction Tools:d-Angle:d1",10/(hypA*xscale))
Setnumbernote(frontimage,"ElectronDiffraction Tools:d-Angle:d2",10/(hypB*xscale))
Setnumbernote(frontimage,"ElectronDiffraction Tools:d-Angle:Angle",totalangle)

deletenote(frontimage,"ElectronDiffraction Tools:Live profile") 
updateimage(frontimage)
}

else
{
//clear tags and lists
	TagGroup TgList1= self.LookUpElement("hklList1") 
	TagGroup TgList2= self.LookUpElement("hklList2") 
	TagGroup uvwList= self.LookUpElement("uvwList")
	uvwList.DLGRemoveElement(0) 
	
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags))
	{
	number tagno=dspacingtags.taggroupcounttags()

	for(number i=0;i<200;i++)
	{
	TgList1.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	}
	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}	
	
image  img:=getfrontimage()
ImageDisplay imgdisp = img.ImageGetImageDisplay(0)
ImageDocument  imgdoc=getfrontimagedocument()
documentwindow imgwin=ImageDocumentGetWindow(imgdoc)
string imgname=img.getname()

number centrex, centrey, x0,y0,x1,y1,x2,y2,x3,y3,x,y,xsize,ysize,xscale,yscale, top,left,bottom,right,width,times
img.GetSize(xsize,ysize)

number i
number annots= imgdisp.ImageDisplayCountROIS()
for (i=0; i<annots; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
 }

annots=imgdisp.componentcountchildrenoftype(13)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(2)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(2,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(5)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(6)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

number radius=xsize/150,smooth

if(radius<1)radius=2

component textannot=newtextannotation(imgdisp,xsize/15, ysize/15,"SPACE to mark the desired center!", 12)
textannot.componentsetfillmode(1)
textannot.componentsetdrawingmode(2) 
textannot.componentsetforegroundcolor(1,1,0)
textannot.componentsetfontfacename("Microsoft Sans Serif")
	
number k=1,m=1,noclick=1,twoclicks=1,n
while(2>m)
{
number keypress=getkey()

getwindowsize(img, xsize, ysize)
 getscale(img,xscale,yscale)
 
number mouse_win_x, mouse_win_y
Number img_view_t, img_view_l, img_view_b, img_view_r
Number v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y
Number img_win_t, img_win_l, img_win_b, img_win_r
Number i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y
Number i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y
Number mouse_img_x, mouse_img_y

windowgetmouseposition(imgwin, mouse_win_x, mouse_win_y)
imgdisp.ImageDisplayGetImageRectInView( img_view_t, img_view_l, img_view_b, img_view_r )
imgdoc.ImageDocumentGetViewToWindowTransform( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y )
ObjectTransformTransformRect( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                            , img_view_t, img_view_l, img_view_b, img_view_r\
                            , img_win_t, img_win_l, img_win_b, img_win_r );

imgdisp.ComponentGetChildToViewTransform( i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y )
		
ObjectTransformCompose( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                      , i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y\
                      , i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y )
		
ObjectTransformUntransformPoint( i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y\
                               , mouse_win_x, mouse_win_y, mouse_img_x, mouse_img_y );

x1=mouse_img_x
y1=mouse_img_y

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)

		
if(keypress==32) //Space bar to mark it
{
number t,l,b,r,leftpart,rightpart
img.GetSize(xsize,ysize)

Image profile_h:= LiveProfile_ExtractLineProfile( img, 0, y1, xsize, y1, 4 )  //horizental projection
Image profile_v:= LiveProfile_ExtractLineProfile( img,x1,0,x1,ysize , 4 )  //vertical projection

number width_h=fwhm(profile_h,x1,leftpart,rightpart)
x1=(leftpart+rightpart)/2

number width_v=fwhm(profile_v,y1,leftpart,rightpart)
y1=(leftpart+rightpart)/2

Component box=NewBoxAnnotation(y1-0.5*radius,x1-0.5*radius,y1+0.5*radius,x1+0.5*radius)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)

SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":X", x1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Y", y1 )
SetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Width", (width_h+width_v)/2)
                                                                                 
img.ClearSelection()
noclick=noclick+1
}

if(keypress==29||keypress==30) //right and up key
{
radius=radius+0.2*radius
img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}

if(keypress==28||keypress==31) //left and down key
{
radius=radius-0.2*radius

if(radius<1)radius=1

img.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}
	
if(keypress>0 && keypress!=32&& keypress!=28&& keypress!=29&& keypress!=30&& keypress!=31||noclick==4) // Any key except space pressed - cancel
{
img.ClearSelection()
m=2
}

}

annots=imgdisp.componentcountchildrenoftype(5)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

for(number z=1;z<4;z++)  //and then refine centers by Gaussian fitting
{
number centrex,centrey
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":X", centrex )
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":Y", centrey )
GetNumberNote(img,"ElectronDiffraction Tools:Live profile:Spot"+z+":Width", width )
radius=0.6*width
width=2*width

img := GetFrontImage()
img.GetSize(xsize,ysize)

image temp=img[centrey-width,centrex-width,centrey+width,centrex+width]
top=centrey-width
left=centrex-width
temp.getsize(xsize,ysize)
temp=( ( icol -(xsize+1)/2 )**2 + ( irow - (ysize+1)/2 )**2 < radius**2 ) ? temp : 0 

number x_0,y_0
times=1
times=1/times
image new := RealImage("target",4,xsize/times,ysize/times)
new = temp.warp(icol*times,irow*times)

number fwhm,leftpart,rightpart

new=new.Gaussfilter(20)
new.GetSize(xsize,ysize)

image xproj=realimage("",4,xsize,1)
image yproj=realimage("",4,xsize,1)

xproj=LiveProfile_ExtractLineProfile( new, 0, (ysize+1)/2 , xsize-1, (ysize+1)/2 ,2 )
yproj= LiveProfile_ExtractLineProfile( new,(xsize+1)/2,0,(xsize+1)/2,ysize-1 , 2 ) 

xproj=xproj-xproj.min()
xproj=10000*xproj/xproj.max()+1

yproj=yproj-yproj.min()
yproj=10000*yproj/yproj.max()+1

number xfwhm=fwhm(xproj,  xsize/2, leftpart,rightpart)
number xpeak=(leftpart+rightpart)/2

number yfwhm=fwhm(yproj,  ysize/2, leftpart,rightpart)
number ypeak=(leftpart+rightpart)/2

number peakintensity, xpeakchannel, ypeakchannel,sigma, ThisRpX,ThisRpY,lastRpX=1,lastRpY=1,weighting=1,loops=20,factor=1

image xpro=realimage("",4,19,1)
image ypro=realimage("",4,19,1)

xproj.SetSelection(0,0.05*xsize,1,0.95*xsize)
for(number j=1;j<20;j++)
{
weighting=j*0.5
image peakfitx=FitGaussian(xproj, weighting, peakintensity, xpeakchannel,sigma, fwhm,ThisRpX)
setpersistentnumbernote("ElectronDiffraction Tools:temp:"+z+":deltaX:"+abs(ThisRpX)+":Peak",xPeakchannel)

xpro.SetPixel(j-1,0,xpeakchannel)
}

yproj.SetSelection(0,0.05*ysize,1,0.95*ysize)
for(number j=1;j<20;j++)
{
weighting=j*0.5
image peakfity=FitGaussian(yproj, weighting, peakintensity, ypeakchannel,sigma, fwhm,ThisRpY)
setpersistentnumbernote("ElectronDiffraction Tools:temp:"+z+":deltaY:"+abs(ThisRpY)+":Peak",yPeakchannel)

ypro.SetPixel(j-1,0,ypeakchannel)
}


xpro.GetSize(xsize,ysize)
x1=left+xpro[0,0.5*xsize,1,xsize].mean()

ypro.GetSize(xsize,ysize)
y1=top+ypro[0,0.5*xsize,1,xsize].mean()

Component box=NewBoxAnnotation(y1-0.5*sigma,x1-0.5*sigma,y1+0.5*sigma,x1+0.5*sigma)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
imgdisp.componentaddchildatend(box)	
setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot"+z+":X",x1)
setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot"+z+":Y",y1)

deletepersistentnote("ElectronDiffraction Tools:temp")
}


Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

string order="1,1"
getstring("Orders of the diffraction spots (m,n)",order,order)

number Bytes=1
width=len(order)

for(number i=0;i<width;i++)
{
string thischar=mid(order,i,1)
if(thischar==",")Bytes=len(left(order,i))
}

m=val(left(order,Bytes))
n=val(right(order,width-Bytes-1))

x1=x2+(x1-x2)/m
y1=y2+(y1-y2)/m

x3=x2+(x3-x2)/n
y3=y2+(y3-y2)/n

Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Setnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

annots=imgdisp.componentcountchildrenoftype(5)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

// Lines are drawn between the points
component lineannotation=newlineannotation(y1,x1,y2,x2)
lineannotation.componentsetfillmode(2) 
lineannotation.componentsetdrawingmode(2)
lineannotation.componentsetforegroundcolor(1,0,1)
imgdisp.componentaddchildatend(lineannotation) 

lineannotation=newlineannotation(y2,x2,y3,x3)
lineannotation.componentsetfillmode(2) 
lineannotation.componentsetdrawingmode(2)
lineannotation.componentsetforegroundcolor(1,0,1)
imgdisp.componentaddchildatend(lineannotation) 

// Computation of the angle using the sine rule
number dxA, dxB, dxC, dyA, dyB, dyC, hypA, hypB, hypC
dxA=x2-x1
dyA=y2-y1
hypA=sqrt(dxA**2+dyA**2)

dxB=x2-x3
dyB=y2-y3
hypB=sqrt(dxB**2+dyB**2)

dxC=x1-x3
dyC=y1-y3
hypC=sqrt(dxC**2+dyC**2)

number gamma=(hypa**2+hypB**2-hypC**2)/(2*hypA*hypB)
number totalangle=acos(gamma)*(180/pi())
totalangle=round(totalangle*10)/10
string outputstring=""+totalangle+" deg."

// Data is output to the results window and onto the image
img.GetScale(xscale,yscale)
result("R1= "+hypA*xscale+" nm-1,   d1= "+10/(hypA*xscale)+" A"+"\n")
result("R2= "+hypB*xscale+" nm-1,   d2= "+10/(hypB*xscale)+" A"+"\n")
result("Angle= "+outputstring+"\n\n")

self.lookupElement("d1field").dlgtitle(" "+format(10/(hypA*xscale), "%2.4f" ))
self.lookupElement("d2field").dlgtitle(" "+format(10/(hypB*xscale), "%2.4f" ))
self.lookupElement("Anglefield").dlgtitle(" "+format(totalangle, "%6.2f" ))

//Save tags
Setnumbernote(img,"ElectronDiffraction Tools:d-Angle:d1",10/(hypA*xscale))
Setnumbernote(img,"ElectronDiffraction Tools:d-Angle:d2",10/(hypB*xscale))
Setnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",totalangle)

SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot1:X",x1)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot1:Y",y1)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot2:X",x2)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot2:Y",y2)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot3:X",x3)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Spot3:Y",y3)

SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:Spot coordinates:Scale:Xscale",xscale)
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:d-Angle:d1",10/(hypA*xscale))
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:d-Angle:d2",10/(hypB*xscale))
SetPersistentNumberNote("ElectronDiffraction Tools:Crystal Tools:d-Angle:Angle",totalangle)

deletenote(img,"ElectronDiffraction Tools:Live profile") 
updateimage(img)
}
}


/*************************************************************************
Calc. Button: to perform pattern indexing or extract d-spacing table
*************************************************************************/
void CalcButton(object self)
{
Result("-----------------------------------------------------------------------------------------\n1) Calculate the d-spacing table of the given crystal.\n2) Ctrl key: to extract an extended d-spacing table from two given spots.\n-----------------------------------------------------------------------------------------\n")

// create a masked image
if(ControlDown())	
{
TagGroup sizefield= self.LookUpElement("sizefield")  //spot size
number spotsize=sizefield.dlggetvalue()
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Spot size",spotsize)

TagGroup Rfield= self.LookUpElement("Rfield")  //RGB
number Rvalue=Rfield.dlggetvalue()
TagGroup Gfield= self.LookUpElement("Gfield") 
number Gvalue=Gfield.dlggetvalue()
TagGroup B1field= self.LookUpElement("B1field") 
number B1value=B1field.dlggetvalue()

number d1=5,d2=5,angle=90,xsize,ysize
image img:=getfrontimage()
imagedisplay imgdisp=img.ImageGetImageDisplay(0)

string imgname=img.GetName()

number xscale=img.imagegetdimensionscale(0)
img.getsize(xsize,ysize)

TagGroup tg = img.ImageGetTagGroup()
TagGroup tg1
if(tg.taggroupgettagastaggroup("ElectronDiffraction Tools:d-Angle", tg1))
{
getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d1",d1)
getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d2",d2)
getnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",angle)
}

getnumber("d1 (A)",d1,d1)
getnumber("d2 (A)",d2,d2)
getnumber("Theta (deg.)",angle,angle)

number annots=imgdisp.componentcountchildrenoftype(9) 
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(9,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(13)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

number angleval,x1,y1,x2,y2,x3,y3,x1_r,y1_r
Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",angleval)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

result("--------------------------------------------------------------------------------------------------\n")
result("The input data: d1="+format(d1, "%2.4f" )+" A, d2="+format(d2, "%2.4f" )+" A, theta="+format(Angle, "%3.3f" )+" deg.\n")

number r1_cal=(10/d1)/xscale
number r1_exp=sqrt((x1-x2)**2+(y1-y2)**2)

x1_r=x2+1*r1_cal*(x1-x2)/r1_exp
y1_r=y2+1*r1_cal*(y1-y2)/r1_exp

//d and angle correction
angle=(180-angle)*pi()/180
number r3_cal=(10/d2)/xscale
number x3_r=(r3_cal/r1_cal)*(((x1-x2)*cos(angle)-(y1-y2)*sin(angle)))+x2
number y3_r=(r3_cal/r1_cal)*(((x1-x2)*sin(angle)+(y1-y2)*cos(angle)))+y2

//----------------------begin to calculate a d-spacing table-------------------------------
number xij,yij,dij
for(number i=-50;i<=50;i++)
{
for(number j=-50;j<50;j++)
{
if(i==0&&j==0)continue

xij=i*(x1_r-x2)+x2+j*(x3_r-x2)
yij=i*(y1_r-y2)+y2+j*(y3_r-y2)

if(xij>0&&xij<xsize&&yij>0&&yij<ysize)
{
dij=10/(sqrt((xij-x2)**2+(yij-y2)**2 )*xscale)

if(xij>=x2)SetstringNote(img,"ElectronDiffraction Tools:Temp:"+format(dij,"%2.4f"), format(dij,"%2.4f"))
}
}
}

//------------------get masked and unmasked images----------------------------
number radius=1
getnumber("Spot size of the array mask ? ",radius,radius)

setNumberNote(img,"ElectronDiffraction Tools:Center:X", x2 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", y2)
img.setorigin(x2,y2)

component mask=newcomponent(9,0,0,0,0) //Array mask
mask.ComponentSetControlPoint(9,x1_r-x2,y1_r-y2,1)//1st spot
mask.ComponentSetControlPoint(10,x3_r-x2,y3_r-y2,1)//2nd spot
mask.ComponentSetControlPoint(4,radius,radius,1)//radius
imgdisp.componentaddchildatend(mask)

number hasmask,edge=3
compleximage maskimage=createmaskfromannotations(imgdisp, edge,0, hasmask)

image maskimg=sqrt(real(maskimage)*real(maskimage)+imaginary(maskimage)*imaginary(maskimage))
maskimg=img.GetPixel(x2,y2)*maskimg/maskimg.max()
maskimg=maskimg-maskimg.min()
maskimg.ShowImage()
maskimg.imagecopycalibrationfrom(img)

setNumberNote(maskimg,"ElectronDiffraction Tools:Center:X", x2 )
setNumberNote(maskimg,"ElectronDiffraction Tools:Center:Y", y2)
maskimg.setorigin(x2,y2)
maskimg.SetName("Masked image-"+img.GetName())
maskimg.setstringnote("ElectronDiffraction Tools:File type","SAED")

//-------------------------------get intensity profiles---------------------
maskimg.GetSize(xsize,ysize)
number samples=ysize
number k=2*Pi()/samples

number minor = min( xsize, ysize )
image linearimg := RealImage( "Linear image-", 4, minor, samples )
linearimg = warp( maskimg,	(icol * sin(irow * k))+ x2, (icol * cos(irow * k))+ y2 )

xscale=maskimg.ImageGetDimensionScale(0)
string unitstring=maskimg.ImageGetDimensionUnitString(0)
linearimg.imagesetdimensioncalibration(0,0,xscale,unitstring,0)
imagesetdimensioncalibration( linearimg, 1, 0, 360/samples, "deg.", 0)

image profile:= RealImage("Intensity profile of "+maskimg.GetName(),4,minor,1) 
profile[icol,0] += linearimg
profile/=samples

for(number i=0;i<0.5*minor;i++)
{
if(profile.GetPixel(i,0)==0)
{
profile[0,0,1,i]=0
break
}
}

//save d-spacing table
profile=10000*profile/profile.max()

profile.ShowImage()
profile.ImageCopyCalibrationFrom(linearimg)
profile.imagesetdimensionorigin(0,0)
profile.ImagesetDimensionCalibration(1,0,1,"Counts",0)
setwindowposition(profile,150, 320)	

profile.setstringnote("ElectronDiffraction Tools:File type","Intensity Profile")

imagedisplay prodisp=profile.ImageGetImageDisplay(0)
prodisp.LinePlotImageDisplaySetContrastLimits( 0, 1.05*profile.max())  
prodisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
prodisp.lineplotimagedisplaysetgridon(0)  //gridon=0
prodisp.lineplotimagedisplaysetbackgroundon(0)  //bgon=0  

tg=img.ImageGetTagGroup()
tg.taggroupgettagastaggroup("ElectronDiffraction Tools:Temp", tg1)
number no=tg1.taggroupcounttags()

string unit
getpersistentStringnote("ElectronDiffraction Tools:Crystal Tools:Default:Unit",unit)
if(unit=="s")
{
xscale=linearimg.ImageGetDimensionScale(0)
profile.imagesetdimensionscale(0,xscale)
profile.imagesetdimensionunitstring(0,"s (1/nm)") 

string filename
If (!SaveAsDialog("Save d-spacing table as   *.dsp   file", GetApplicationDirectory(2,0) +imgname+ ".dsp", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
WriteFile(fileID, "\n"+no+"    DSPACE    0.1     US     RIR= 0\n") 

for(number i=no-1;i>-1;i--)
{
number d =val( tg1.TagGroupGetTaglabel( i ))
number intensity=profile.GetPixel(10/(d*xscale),0)
WriteFile(fileID, d+"        "+intensity+"\n")
}

CloseFile(fileID)
}

if(unit=="Q")
{
xscale=linearimg.ImageGetDimensionScale(0)/10*2*pi()
profile.imagesetdimensionscale(0,xscale)
profile.imagesetdimensionunitstring(0,"Q (1/A)")

string filename
If (!SaveAsDialog("Save d-spacing table as   *.dsp   file", GetApplicationDirectory(2,0) +imgname+ ".dsp", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
WriteFile(fileID, "\n"+no+"    DSPACE    0.1     US     RIR= 0\n")  

for(number i=no-1;i>-1;i--)
{
number d =val( tg1.TagGroupGetTaglabel( i ))
number intensity=profile.GetPixel(10/(d*linearimg.ImageGetDimensionScale(0)),0)
WriteFile(fileID, d+"        "+intensity+"\n")
}

CloseFile(fileID)
}
				
deletenote(img,"ElectronDiffraction Tools:Temp")
}

//
else if(OptionDown())
{
TagGroup sizefield= self.LookUpElement("sizefield")  //spot size
number spotsize=sizefield.dlggetvalue()
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Spot size",spotsize)

TagGroup Rfield= self.LookUpElement("Rfield")  //RGB
number Rvalue=Rfield.dlggetvalue()
TagGroup Gfield= self.LookUpElement("Gfield") 
number Gvalue=Gfield.dlggetvalue()
TagGroup B1field= self.LookUpElement("B1field") 
number B1value=B1field.dlggetvalue()

number d1=5,d2=5,angle=90,xsize,ysize
image img:=getfrontimage()
imagedisplay imgdisp=img.ImageGetImageDisplay(0)

string imgname=img.GetName()
img.getsize(xsize,ysize)

number xscale=img.imagegetdimensionscale(0)

TagGroup tg = img.ImageGetTagGroup()
TagGroup tg1
if(tg.taggroupgettagastaggroup("ElectronDiffraction Tools:d-Angle", tg1))
{
getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d1",d1)
getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d2",d2)
getnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",angle)
}

getnumber("d1 (A)",d1,d1)
getnumber("d2 (A)",d2,d2)
getnumber("Theta (deg.)",angle,angle)

number annots=imgdisp.componentcountchildrenoftype(6) 
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(13)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

number angleval,x1,y1,x2,y2,x3,y3,x1_r,y1_r
Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",angleval)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

result("--------------------------------------------------------------------------------------------------\n")
result("The input data: d1="+format(d1, "%2.4f" )+" A, d2="+format(d2, "%2.4f" )+" A, theta="+format(Angle, "%3.3f" )+" deg.\n")

number r1_cal=(10/d1)/xscale
number r1_exp=sqrt((x1-x2)**2+(y1-y2)**2)

x1_r=x2+1*r1_cal*(x1-x2)/r1_exp
y1_r=y2+1*r1_cal*(y1-y2)/r1_exp
Component spot1=newovalannotation(y1_r-spotsize, x1_r-spotsize, y1_r+spotsize, x1_r+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
imgdisp.componentaddchildatend(spot1)

//d and angle correction
angle=(180-angle)*pi()/180 
number r3_cal=(10/d2)/xscale
number x3_r=(r3_cal/r1_cal)*(((x1-x2)*cos(angle)-(y1-y2)*sin(angle)))+x2
number y3_r=(r3_cal/r1_cal)*(((x1-x2)*sin(angle)+(y1-y2)*cos(angle)))+y2
spot1=newovalannotation(y3_r-spotsize, x3_r-spotsize, y3_r+spotsize, x3_r+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
imgdisp.componentaddchildatend(spot1)

number xij,yij,dij
for(number i=-20;i<=20;i++)
{
for(number j=-20;j<20;j++)
{
if(i==0&&j==0)continue //skip for 00

xij=i*(x1_r-x2)+x2+j*(x3_r-x2)
yij=i*(y1_r-y2)+y2+j*(y3_r-y2)

if(xij>0&&xij<xsize&&yij>0&&yij<ysize)
{
dij=10/(sqrt((xij-x2)**2+(yij-y2)**2 )*xscale)

if(xij>=x2)SetstringNote(img,"ElectronDiffraction Tools:Temp:"+format(dij,"%2.4f"), format(dij,"%2.4f"))

spot1=newovalannotation(yij-spotsize, xij-spotsize, yij+spotsize, xij+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
imgdisp.componentaddchildatend(spot1)
}
}
}

//save d-spacing table
documentwindow dspacingtable=newscriptwindow("d-spacing Table-"+imgname,400, 200, 650, 800)
windowshow(dspacingtable)
windowselect(dspacingtable)
		
tg=img.ImageGetTagGroup()
tg.taggroupgettagastaggroup("ElectronDiffraction Tools:Temp", tg1)
number no=tg1.taggroupcounttags()
editorwindowaddtext(dspacingtable,"\n"+no+"    DSPACE    0.1     US     RIR= 0\n")

for(number i=no-1;i>-1;i--)
{
number d =val( tg1.TagGroupGetTaglabel( i ))
editorwindowaddtext(dspacingtable,d+"       100\n")
}

string pathname
if (!saveasdialog("","d-spacing Table of "+imgname+".dsp", pathname)) Exit (0)
editorwindowsavetofile(dspacingtable,pathname)

deletenote(img,"ElectronDiffraction Tools:Temp")
}


else
{
image img:=GetFrontImage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

number annots=imgdisp.componentcountchildrenoftype(6) 
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(13)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

	number crystalvalue, alpha, beta, gamma, a, b, c, crysttabs, h, k, l, dspacing,no=1
	number maxhklval, mindval
	string crystalname

	//clear all list and tags
	TagGroup TgList1= self.LookUpElement("hklList1") 
	TagGroup TgList2= self.LookUpElement("hklList2") 
	TagGroup uvwList= self.LookUpElement("uvwList")
	uvwList.DLGRemoveElement(0)  
	
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags))
	{
	number tagno=dspacingtags.taggroupcounttags()

	for(number i=0;i<200;i++)
	{
	TgList1.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	}
	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}

//------------------1) calculate d-spacing table	------------------
	// read lattice parameters
	self.dlggetvalue("crystalradios",crystalvalue)
	self.dlggetvalue("alphafield",alpha)
	self.dlggetvalue("betafield",beta)
	self.dlggetvalue("gammafield",gamma)
	self.dlggetvalue("afield",a)
	self.dlggetvalue("bfield",b)
	self.dlggetvalue("cfield",c)
	
		getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:maxhkl",maxhklval)
		getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:mind",mindval)

		alpha=alpha/(180/pi())
		beta=beta/(180/pi())
		gamma=gamma/(180/pi())
		
		taggroup mcptags
	    ptags=getpersistenttaggroup()

		if(a<=0 || b<=0 || c<=0 || alpha<=0 || beta<=0 || gamma<=0 )
			{
				beep()
				exit(0)
			}

	// begin to calculate d-spacing
	for(h=maxhklval; h>=-maxhklval; h--)
	{
		for(k=maxhklval; k>=-maxhklval; k--)
		{
			for(l=maxhklval; l>=-maxhklval; l--)
				{
					if(h**2+k**2+l**2==0) continue

					if(crystalvalue==1) // Cubic
					{
						dspacing=sqrt(a**2/(h**2+k**2+l**2))
						crystalname="Cubic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1
					}

					if(crystalvalue==2) // Tetragonal
					{
						number temp=((1/a**2)*(h**2+k**2))+((1/c**2)*l**2)
						dspacing=sqrt(1/temp)
						crystalname="Tetragonal"

						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==3) //Orthorhombic
					{
						number temp=(((1/a**2)*h**2)+((1/b**2)*k**2)+((1/c**2)*l**2))
						dspacing=sqrt(1/temp)
						crystalname="Orthorhombic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}
	
					if(crystalvalue==4) //Hexagonal
					{

						number temp=((4/(3*a**2))*(h**2+(h*k)+k**2)+((1/c**2)*l**2))
						dspacing=sqrt(1/temp)
						crystalname="Hexagonal"

						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==5) //Rhombohedral
					{
						number temp1=((h**2+k**2+l**2)*sin(alpha)**2)+(2*((h*k)+(k*l)+(h*l))*(cos(alpha)**2-cos(alpha)))
						number temp2=a**2*(1-(3*cos(alpha)**2)+(2*cos(alpha)**3))
						number temp3=temp1/temp2
						dspacing=sqrt(1/temp3)
						crystalname="Rhombohedral"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}	

					if(crystalvalue==6) //Monoclinic
					{
						number temp=((1/a**2)*(h**2/sin(beta)**2))+((1/b**2)*k**2)+((1/c**2)*(l**2/sin(beta)**2))-((2*h*l*cos(beta))/(a*c*sin(beta)**2))
						dspacing=sqrt(1/temp)
						crystalname="Monoclinic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==7) //Triclinic
					{
						number V2=a**2*b**2*c**2*(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+(2*cos(alpha)*cos(beta)*cos(gamma)))
						number s11=b**2*c**2*sin(alpha)**2
						number s22=a**2*c**2*sin(beta)**2
						number s33=a**2*b**2*sin(gamma)**2

						number s12=a*b*c**2*(cos(alpha)*cos(beta)-cos(gamma))
						number s23=a**2*b*c*(cos(beta)*cos(gamma)-cos(alpha))
						number s31=a*b**2*c*(cos(gamma)*cos(alpha)-cos(beta))

						number temp=(1/v2)*((s11*h**2)+(s22*k**2)+(s33*l**2)+(2*s12*h*k)+(2*s23*k*l)+(2*s31*l*h))
						dspacing=sqrt(1/temp)
						crystalname="Triclinic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

				}		
		}
	}

//------------------2) search d with the given error------------------
number d1val, d2val,Error_d
Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d1",d1val)  	
Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d1",d2val)  
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Error d",Error_d)
	 
number maxd1=d1val+Error_d
number mind1=d1val-Error_d

ptags=getpersistenttaggroup()	
ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)

number peakno=dspacingtags.taggroupcounttags()

for(number i=1; i<peakno+1; i++)
       {
          number d,h,k,l
		  string mcplabel = dspacingtags.TagGroupGetTaglabel( i-1 )
		  
		  d=val(mcplabel)		  
		  number delta_d1=abs(d1val-d )   
		  		  
		   if(d<maxd1 && d>mind1)  
                {
				 taggroup peaktags,peaks
				 peaktags=getpersistenttaggroup()
				 peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d, peaks)
				 number peaknos=peaks.taggroupcounttags()
				 
				String peaklabel = peaks.TagGroupGetTaglabel(0 )
				number peakvalue=val(peaklabel)

				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":h",h)
                getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":k",k)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":l",l)
					
                Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d1+":delta d",delta_d1) 
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d1+":d",d)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d1+":h",h)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d1+":k",k)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d1+":l",l)
				}
		}

//order the d-spacing table
	 taggroup peaktags,peaks
	 peaktags=getpersistenttaggroup()
	 peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1",peaks)
	 no=peaks.taggroupcounttags()
	 for(number i=0;i<no;i++)
	 {
	 String peaklabel = peaks.TagGroupGetTaglabel(i )
	 number peakvalue=val(peaklabel)

	 number d
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+peakvalue+":d",d)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+peakvalue+":h",h)
     getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+peakvalue+":k",k)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+peakvalue+":l",l)	
	 
	 //update the list1
	 TgList1.DLGAddListItem(format(100*peakvalue/d, "%2.1f")+"|("+h+","+k+","+l+")",0)
	 }
	 }
}

/*************************************************************************
hklList1Changed: to response when hkl list1 is changed
*************************************************************************/
void hklList1Changed(object self, taggroup tg)
{
//clear lists
TagGroup uvwList= self.LookUpElement("uvwList")
uvwList.DLGRemoveElement(0)  

TagGroup sizefield= self.LookUpElement("sizefield")  //spot size
number spotsize=sizefield.dlggetvalue()
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Spot size",spotsize)

TagGroup Rfield= self.LookUpElement("Rfield")  //RGB值
number Rvalue=Rfield.dlggetvalue()
TagGroup Gfield= self.LookUpElement("Gfield")  
number Gvalue=Gfield.dlggetvalue()
TagGroup B1field= self.LookUpElement("B1field")
number B1value=B1field.dlggetvalue()

taggroup peaktags,peaks
peaktags=getpersistenttaggroup()
peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", peaks)

number no=	tg.dlggetvalue() 		 
String peaklabel = peaks.TagGroupGetTaglabel(no)
number delta_d=val(peaklabel)

number d1,h1,k1,l1,d2,h2,k2,l2,deltad1
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d+":d",d1)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d+":h",h1)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d+":k",k1)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d+":l",l1)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List1:"+delta_d+":delta d",deltad1)

	TagGroup TgList2= self.LookUpElement("hklList2") 

	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()

	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", dspacingtags)) 
	{
	number tagno=dspacingtags.taggroupcounttags()
	for(number i=0;i<200;i++)
	{
	TgList2.DLGRemoveElement(0)
	}
	
	deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List2")
	}

number Error_d,Error_angle  
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Error d",Error_d)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Error Angle",Error_angle)	
	 
number d1val,d2val,angleval  
image img:=GetFrontImage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

//display the selected spots on the SAED
number annots=imgdisp.componentcountchildrenoftype(6) 
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(13)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

number xsize,ysize
number xscale=img.imagegetdimensionscale(0)
img.getsize(xsize,ysize)

number x1,y1,x2,y2,x,y
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)

number r1_cal=(10/d1)/xscale
number r1_exp=sqrt((x1-x2)**2+(y1-y2)**2)

x=x2-1*r1_cal*(x2-x1)/r1_exp
y=y2-1*r1_cal*(y2-y1)/r1_exp

Component spot1=newovalannotation(Y-spotsize, X-spotsize, Y+spotsize, X+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
imgdisp.componentaddchildatend(spot1)

string fonttype="Microsoft Sans Serif"
component textannot1=newtextannotation(x1, y1, ""+h1+","+k1+","+l1, 2*spotsize)
textannot1.componentsetdrawingmode(2)
textannot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
textannot1.componentsetfontfacename(fonttype) 
imgdisp.componentaddchildatend(textannot1)

Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d1",d1val)
Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:d2",d2val)
Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",angleval)

number maxd2=d2val+Error_d
number mind2=d2val-Error_d	
number maxAngle=angleval+Error_angle
number minAngle=angleval-Error_angle

//read lattice parameters of the selected crystal
number crystalvalue,alpha,beta,gamma,a,b,c,angle,interAngle
string crystalname
		self.dlggetvalue("alphafield",alpha)
		self.dlggetvalue("betafield",beta)
		self.dlggetvalue("gammafield",gamma)
		self.dlggetvalue("afield",a)
		self.dlggetvalue("bfield",b)
		self.dlggetvalue("cfield",c)

		self.dlggetvalue("crystalradios", crystalvalue)

		// Set the crystal name based on the radio
		if(crystalvalue==1) crystalname="C"
		if(crystalvalue==2) crystalname="T"
		if(crystalvalue==3) crystalname="O"
		if(crystalvalue==4) crystalname="H"
		if(crystalvalue==5) crystalname="R"
		if(crystalvalue==6) crystalname="M"
		if(crystalvalue==7) crystalname="A"
		
 		alpha=alpha/(180/pi())
		beta=beta/(180/pi())
		gamma=gamma/(180/pi())
		angle=angle/(180/pi())

//read tags and match d2
ptags=getpersistenttaggroup()	
dspacingtags 
ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)
number peakno=dspacingtags.taggroupcounttags()

for(number i=1; i<peakno+1; i++)
       {
		  string mcplabel = dspacingtags.TagGroupGetTaglabel( i-1 )		  
		  d2=val(mcplabel)
		  
		  number deltad2=abs(d2val-d2 )   
		  		 
		   if(d2<maxd2 && d2>mind2)
                {
				 taggroup peaktags,peaks
				 peaktags=getpersistenttaggroup()
				 peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d2, peaks)
				 number peaknos=peaks.taggroupcounttags()
				
				 for(number j=1;j<peaknos+1;j++)
				{
				String peaklabel = peaks.TagGroupGetTaglabel(j-1)
				number peakvalue=val(peaklabel)

				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d2+":"+peakvalue+":h",h2)
                getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d2+":"+peakvalue+":k",k2)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d2+":"+peakvalue+":l",l2)
				
					if(h2**2+k2**2+l2**2==0) continue			

		//-----------------------begin to calculate angles-----------------------------			
		if(crystalvalue==1) //Cubic
			{
				number costheta, temp1, temp2
				temp1=(h1*h2)+(k1*k2)+(l1*l2)
				temp2=sqrt((h1**2+k1**2+l1**2)*(h2**2+k2**2+l2**2))

				costheta=temp1/temp2
				interAngle=acos(costheta)
				interAngle=interAngle*(180/pi())
				crystalname="Cubic"
			}

		if(crystalvalue==2) //Tetragonal
			{
				number costheta, temp1, temp2, temp3, temp4
				temp1=((1/a**2)*((h1*h2)+(k1*k2)))+((1/c**2)*l1*l2)
				temp2=((1/a**2)*(h1**2+k1**2))+((1/c**2)*l1**2)
				temp3=((1/a**2)*(h2**2+k2**2))+((1/c**2)*l2**2)
				temp4=sqrt(temp2*temp3)

				costheta=temp1/temp4
				interAngle=acos(costheta)
				interAngle=interAngle*(180/pi())
				crystalname="Tetragonal"
			}

		if(crystalvalue==3) //Orthorhombic
			{
				number costheta, temp1, temp2, temp3, temp4
				temp1=((1/a**2)*h1*h2)+((1/b**2)*k1*k2)+((1/c**2)*l1*l2)
				temp2=((1/a**2)*h1**2)+((1/b**2)*k1**2)+((1/c**2)*l1**2)
				temp3=((1/a**2)*h2**2)+((1/b**2)*k2**2)+((1/c**2)*l2**2)
				temp4=sqrt(temp2*temp3)

				costheta=temp1/temp4
				interAngle=acos(costheta)
				interAngle=interAngle*(180/pi())
				crystalname="Hexagonal"
			}

		if(crystalvalue==4) //六角
			{
				number costheta, temp1, temp2, temp3, temp4
				temp1=(h1*h2)+(k1*k2)+(0.5*((h1*k2)+(k1*h2)))+(((3*a**2)/(4*c**2))*l1*l2)
				temp2=h1**2+k1**2+(h1*k1)+(((3*a**2)/(4*c**2))*l1**2)
				temp3=h2**2+k2**2+(h2*k2)+(((3*a**2)/(4*c**2))*l2**2)
				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4

				interAngle=acos(costheta)
				interAngle=interAngle*(180/pi())
				crystalname="Hexagonal"
			}

		if(crystalvalue==5) //Rhombohedral
			{
				number costheta, temp, temp1, temp2, temp3, temp4, cellvol, d1, d2

				temp=(1/a**2)
				temp3=(1+cos(alpha)-(2*(cos(alpha)**2)))
				temp2=(1+cos(alpha))*((h1**2+k1**2+l1**2)-((1-tan(0.5*alpha)**2)*((h1*k1)+(k1*l1)+(l1*h1))))
				temp4=(temp*temp2)/temp3
				d1=sqrt(1/temp4)

				temp2=(1+cos(alpha))*((h2**2+k2**2+l2**2)-((1-tan(0.5*alpha)**2)*((h2*k2)+(k2*l2)+(l2*h2))))
				temp4=(temp*temp2)/temp3
				d2=sqrt(1/temp4)

				cellvol=a**3*sqrt((1-(3*cos(alpha)**2)+(2*cos(alpha)**3)))

				temp1=(a**4*d1*d2)/cellvol**2
				temp2=sin(alpha)**2*((h1*h2)+(k1*k2)+(l1*l2))
				temp3=(cos(alpha)**2-cos(alpha))*((k1*l2)+(k2*l1)+(l1*h2)+(l2*h1)+(h1*k2)+(h2*k1))
				costheta=temp1*(temp2+temp3)

				interAngle=acos(costheta)
				interAngle=interAngle*(180/pi())
				crystalname="Rhombohedral"
			}

		if(crystalvalue==6) //Monoclinic
			{
				number costheta, temp1, temp2, temp3, temp4

				temp1=((1/a**2)*h1*h2)+((1/b**2)*k1*k2*sin(beta)**2)+((1/c**2)*l1*l2)-((1/(a*c))*((l1*h2)+(l2*h1))*cos(beta))
				temp2=((1/a**2)*h1**2)+((1/b**2)*k1**2*sin(beta)**2)+((1/c**2)*l1**2)-(((2*h1*l1)/(a*c))*cos(beta))
				temp3=((1/a**2)*h2**2)+((1/b**2)*k2**2*sin(beta)**2)+((1/c**2)*l2**2)-(((2*h2*l2)/(a*c))*cos(beta))	
				temp4=sqrt(temp2*temp3)
				costheta=temp1/temp4

				interAngle=acos(costheta)
				interAngle=interAngle*(180/pi())
				crystalname="Monoclinic"
			}

		if(crystalvalue==7) // Triclinic
			{
				number costheta, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12

				temp1=(h1*h2*b**2*c**2*sin(alpha)**2)+(k1*k2*a**2*c**2*sin(beta)**2)+(l1*l2*a**2*b**2*sin(gamma)**2)
				temp2=a*b*c**2*((cos(alpha)*cos(beta))-cos(gamma))*((k1*h2)+(h1*k2))
				temp3=a*b**2*c*((cos(gamma)*cos(alpha))-cos(beta))*((h1*l2)+(l1*h2))
				temp4=a**2*b*c*((cos(beta)*cos(gamma))-cos(alpha))*((k1*l2)+(l1*k2))

				number fval=temp1+temp2+temp3+temp4

				temp5=(h1**2*b**2*c**2*sin(alpha)**2)+(k1**2*a**2*c**2*sin(beta)**2)+(l1**2*a**2*b**2*sin(gamma)**2)
				temp6=2*h1*k1*a*b*c**2*((cos(alpha)*cos(beta))-cos(gamma))
				temp7=2*h1*l1*a*b**2*c*((cos(gamma)*cos(alpha))-cos(beta))
				temp8=2*k1*l1*a**2*b*c*((cos(beta)*cos(gamma))-cos(alpha))

				number ahikili=sqrt(temp5+temp6+temp7+temp8)

				temp9=(h2**2*b**2*c**2*sin(alpha)**2)+(k2**2*a**2*c**2*sin(beta)**2)+(l2**2*a**2*b**2*sin(gamma)**2)
				temp10=2*h2*k2*a*b*c**2*((cos(alpha)*cos(beta))-cos(gamma))
				temp11=2*h2*l2*a*b**2*c*((cos(gamma)*cos(alpha))-cos(beta))
				temp12=2*k2*l2*a**2*b*c*((cos(beta)*cos(gamma))-cos(alpha))

				number ah2k2l2=sqrt(temp9+temp10+temp11+temp12)

				costheta=fval/(ahikili*ah2k2l2)
				interAngle=acos(costheta)
				interAngle=interAngle*(180/pi())
				crystalname="Triclinic"
			}
			
			number deltaAngle
			deltaAngle=abs(Angleval-interAngle)

			//for the matched hkl planes
			if(deltaAngle<=Error_Angle)
                {					
				//calculate zone axis direction	   
				number u=(k1*l2)-(k2*l1)
				number v=(l1*h2)-(l2*h1)
				number w=(h1*k2)-(h2*k1)  
	  
				number divisor=CrystUIFrame.NormaliseIndices(u,v,w)
				u=round(u*100)/100
				v=round(v*100)/100
				w=round(w*100)/100
				
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":d1",d1)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":h1",h1)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":k1",k1)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":l1",l1)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta d1",deltad1)
				
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":d2",d2)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":h2",h2)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":k2",k2)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":l2",l2)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta d2",deltad2)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta Angle",deltaangle)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":u",u)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":v",v)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":w",w)
				setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":Angle",interAngle)
				}
				}
				}		
}

//update the lists
number deltad2,deltaangle
taggroup list2peaks
	 peaktags=getpersistenttaggroup()
	 peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2",list2peaks)
	 number list2no=list2peaks.taggroupcounttags()

	 for(number i=0;i<list2no;i++)
	 {
	 String peaklabel = list2peaks.TagGroupGetTaglabel(i )
	 number peakvalue=val(peaklabel)

	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+peakvalue+":d2",d2)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+peakvalue+":h2",h2)
     getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+peakvalue+":k2",k2)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+peakvalue+":l2",l2)	
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+peakvalue+":delta d2",deltad2)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+peakvalue+":delta Angle",deltaangle)
	 
	 TgList2.DLGAddListItem(format(100*peakvalue/d2, "%2.1f")+"|"+format(100*deltaangle/angleval, "%2.1f") +"|("+h2+","+k2+","+l2+")",0)
	 }
}

/*************************************************************************
hklList2Changed: to reponse when hkl list2 is changed
*************************************************************************/
void hklList2Changed(object self, taggroup tg)
{
number no=	tg.dlggetvalue()
TagGroup TgList2= self.LookUpElement("hklList2")  
TagGroup uvwList= self.LookUpElement("uvwList") 
uvwList.DLGRemoveElement(0)

TagGroup sizefield= self.LookUpElement("sizefield")  //spot size
number spotsize=sizefield.dlggetvalue()
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Spot size",spotsize)

TagGroup Rfield= self.LookUpElement("Rfield")  //RGB
number Rvalue=Rfield.dlggetvalue()
TagGroup Gfield= self.LookUpElement("Gfield") 
number Gvalue=Gfield.dlggetvalue()
TagGroup B1field= self.LookUpElement("B1field") 
number B1value=B1field.dlggetvalue()

taggroup peaktags,peaks
peaktags=getpersistenttaggroup()
peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", peaks)

String peaklabel = peaks.TagGroupGetTaglabel(no)
number deltad2=val(peaklabel)

number u,v,w,d1,h1,k1,l1,d2,h2,k2,l2,angle,deltaangle,deltad1
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":u",u)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":v",v)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":w",w)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":d2",d2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":h2",h2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":k2",k2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":l2",l2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":Angle",angle)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta d2",deltad2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta Angle",deltaangle)	

Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":d1",d1)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":h1",h1)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":k1",k1)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":l1",l1)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta d1",deltad1)

uvwList.DLGAddListItem("["+u + ", "+v+", "+w+"]",0)

//clear marks
image img:=GetFrontImage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

number annots=imgdisp.componentcountchildrenoftype(6)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(13)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

number xsize,ysize
number xscale=img.imagegetdimensionscale(0)
img.getsize(xsize,ysize)

number angleval,x1,y1,x2,y2,x3,y3,x1_r,y1_r
Getnumbernote(img,"ElectronDiffraction Tools:d-Angle:Angle",angleval)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

result("--------------------------------------------------------------------------------------------------\n")
result("d1="+format(d1, "%2.4f" )+"("+format(100*deltad1/d1, "%2.1f" )+"), d2="+format(d2, "%2.4f" )+"("+format(100*deltad2/d2, "%2.1f" )+"|"+format(100*deltaangle/angle, "%2.1f" )+"), Angle="+format(Angle, "%2.1f" )+"|"+format(angleval, "%2.1f" )+", ("+h1+","+k1+","+l1+") <--> ("+h2+","+k2+","+l2+ ") ==> ["+u+","+v+","+w+"]\n")

number r1_cal=(10/d1)/xscale
number r1_exp=sqrt((x1-x2)**2+(y1-y2)**2)

x1_r=x2+1*r1_cal*(x1-x2)/r1_exp
y1_r=y2+1*r1_cal*(y1-y2)/r1_exp
Component spot1=newovalannotation(y1_r-spotsize, x1_r-spotsize, y1_r+spotsize, x1_r+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
imgdisp.componentaddchildatend(spot1)

string fonttype="Microsoft Sans Serif"
component textannot1=newtextannotation(x1_r, y1_r, ""+h1+","+k1+","+l1, 2*spotsize)
textannot1.componentsetdrawingmode(2)
textannot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
textannot1.componentsetfontfacename(fonttype) 
imgdisp.componentaddchildatend(textannot1)

//d and angle corrections
angle=(180-angle)*pi()/180 //角度转换
number r3_cal=(10/d2)/xscale
number x3_r=(r3_cal/r1_cal)*(((x1-x2)*cos(angle)-(y1-y2)*sin(angle)))+x2
number y3_r=(r3_cal/r1_cal)*(((x1-x2)*sin(angle)+(y1-y2)*cos(angle)))+y2

if(sqrt((x3_r-x3)**2+(y3_r-y3)**2)>sqrt((x3-x2)**2+(y3-y2)**2))
{
x3_r=2*x2-x3_r
y3_r=2*y2-y3_r
}

spot1=newovalannotation(y3_r-spotsize, x3_r-spotsize, y3_r+spotsize, x3_r+spotsize)

spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)

fonttype="Microsoft Sans Serif"
textannot1=newtextannotation(x3,y3, ""+h2+","+k2+","+l2, 2*spotsize)
textannot1.componentsetdrawingmode(2)
textannot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
textannot1.componentsetfontfacename(fonttype) 

number h3,k3,l3
number xij,yij
for(number i=-10;i<=10;i++)
{
for(number j=-10;j<10;j++)
{
if(i==0&&j==0)continue //skip for 00

xij=i*(x1_r-x2)+x2+j*(x3_r-x2)
yij=i*(y1_r-y2)+y2+j*(y3_r-y2)

if(xij>0&&xij<xsize&&yij>0&&yij<ysize)
{
spot1=newovalannotation(yij-spotsize, xij-spotsize, yij+spotsize, xij+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
imgdisp.componentaddchildatend(spot1)

h3=i*h1+j*h2
k3=i*k1+j*k2
l3=i*l1+j*l2

fonttype="Microsoft Sans Serif"
textannot1=newtextannotation(xij,yij, ""+h3+","+k3+","+l3, 2*spotsize)
textannot1.componentsetdrawingmode(2)
textannot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
textannot1.componentsetfontfacename(fonttype) 
imgdisp.componentaddchildatend(textannot1)
}
}
}
}

/*************************************************************************
ZoneListChanged: to response when zone list is changed
*************************************************************************/
void ZoneListChanged(object self, taggroup tg)
{
TagGroup TgList2= self.LookUpElement("hklList2") 
number no=	TgList2.dlggetvalue() 

TagGroup sizefield= self.LookUpElement("sizefield")  //spot size
number spotsize=sizefield.dlggetvalue()
TagGroup Rfield= self.LookUpElement("Rfield")  //RGB
number Rvalue=Rfield.dlggetvalue()
TagGroup Gfield= self.LookUpElement("Gfield") 
number Gvalue=Gfield.dlggetvalue()
TagGroup B1field= self.LookUpElement("B1field")
number B1value=B1field.dlggetvalue()

taggroup peaktags,peaks
peaktags=getpersistenttaggroup()
peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List2", peaks)

String peaklabel = peaks.TagGroupGetTaglabel(no)
number deltad2=val(peaklabel)

number u,v,w,d1,h1,k1,l1,d2,h2,k2,l2,angle,deltaangle,deltad1
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":u",u)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":v",v)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":w",w)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":d2",d2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":h2",h2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":k2",k2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":l2",l2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":Angle",angle)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta d2",deltad2)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta Angle",deltaangle)	

Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":d1",d1)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":h1",h1)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":k1",k1)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":l1",l1)
Getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:List2:"+deltad2+":delta d1",deltad1)

//Clear marks
image img:=GetFrontImage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

number annots=imgdisp.componentcountchildrenoftype(6)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

annots=imgdisp.componentcountchildrenoftype(13)
for (number i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

number xsize,ysize
number xscale=img.imagegetdimensionscale(0)
img.getsize(xsize,ysize)

number x1,y1,x2,y2,x3,y3,x1_r,y1_r
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:X",x2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot2:Y",y2)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:X",x3)
Getnumbernote(img,"ElectronDiffraction Tools:Spot coordinates:Spot3:Y",y3)

//the first measured spot
Component spot1=newovalannotation(y1-spotsize, x1-spotsize, y1+spotsize, x1+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
imgdisp.componentaddchildatend(spot1)

string fonttype="Microsoft Sans Serif"
component textannot1=newtextannotation(x1, y1, ""+h1+","+k1+","+l1, 2*spotsize)
textannot1.componentsetdrawingmode(2)
textannot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
textannot1.componentsetfontfacename(fonttype) 
imgdisp.componentaddchildatend(textannot1)

//the corrected spots
number ratio=d2/d1
number r1_cal=sqrt((x1-x2)**2+(y1-y2)**2)
number d1_cal=10/(r1_cal*xscale)
number r3_cal=10/(d1_cal*ratio*xscale)
angle=(180-angle)*pi()/180

number x3_r=(r3_cal/r1_cal)*(((x1-x2)*cos(angle)-(y1-y2)*sin(angle)))+x2
number y3_r=(r3_cal/r1_cal)*(((x1-x2)*sin(angle)+(y1-y2)*cos(angle)))+y2

if(sqrt((x3_r-x3)**2+(y3_r-y3)**2)>sqrt((x3-x2)**2+(y3-y2)**2))
{
x3_r=2*x2-x3_r
y3_r=2*y2-y3_r
}

spot1=newovalannotation(y3_r-spotsize, x3_r-spotsize, y3_r+spotsize, x3_r+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)

fonttype="Microsoft Sans Serif"
textannot1=newtextannotation(x3, y3, ""+h2+","+k2+","+l2, 2*spotsize)
textannot1.componentsetdrawingmode(2)
textannot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
textannot1.componentsetfontfacename(fonttype) 

number xij,yij, h3,k3,l3
for(number i=-10;i<=10;i++)
{
for(number j=-10;j<10;j++)
{
if(i==0&&j==0)continue //skip for 00

xij=i*(x1-x2)+x2+j*(x3_r-x2)
yij=i*(y1-y2)+y2+j*(y3_r-y2)

if(xij>0&&xij<xsize&&yij>0&&yij<ysize)
{
spot1=newovalannotation(yij-spotsize, xij-spotsize, yij+spotsize, xij+spotsize)
spot1.componentsetfillmode(1)
spot1.componentsetdrawingmode(2)
spot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
imgdisp.componentaddchildatend(spot1)

h3=i*h1+j*h2
k3=i*k1+j*k2
l3=i*l1+j*l2

fonttype="Microsoft Sans Serif"
textannot1=newtextannotation(xij,yij, ""+h3+","+k3+","+l3, 2*spotsize)
textannot1.componentsetdrawingmode(2)
textannot1.componentsetforegroundcolor(Rvalue,Gvalue,B1value)
textannot1.componentsetfontfacename(fonttype) 
imgdisp.componentaddchildatend(textannot1)
}

}
}
}

//Function: to response when d-spacing is changed in d-hkl sub-pannel
void dspacingfieldChanged(object self, taggroup tg)
{
	number crystalvalue, alpha, beta, gamma, a, b, c, crysttabs, h, k, l, dspacing,no=1
	number maxhklval, mindval
	string crystalname

	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:Temp", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Temp")
	
	// read lattice parameters
	self.dlggetvalue("crystalradios",crystalvalue)
	self.dlggetvalue("alphafield",alpha)
	self.dlggetvalue("betafield",beta)
	self.dlggetvalue("gammafield",gamma)
	self.dlggetvalue("afield",a)
	self.dlggetvalue("bfield",b)
	self.dlggetvalue("cfield",c)
	
		getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:maxhkl",maxhklval)
		getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:mind",mindval)

		alpha=alpha/(180/pi())
		beta=beta/(180/pi())
		gamma=gamma/(180/pi())
		
		taggroup mcptags
	    ptags=getpersistenttaggroup()

		if(a<=0 || b<=0 || c<=0 || alpha<=0 || beta<=0 || gamma<=0 )
			{
				beep()
				exit(0)
			}

	// begin to calculate d-spacing table
	for(h=maxhklval; h>=-maxhklval; h--)
	{
		for(k=maxhklval; k>=-maxhklval; k--)
		{
			for(l=maxhklval; l>=-maxhklval; l--)
				{
				//skip for 000
					if(h**2+k**2+l**2==0) continue

					if(crystalvalue==1) // Cubic
					{
						dspacing=sqrt(a**2/(h**2+k**2+l**2))
						crystalname="Cubic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1
					
					}

					if(crystalvalue==2) // Tetragonal
					{
						number temp=((1/a**2)*(h**2+k**2))+((1/c**2)*l**2)
						dspacing=sqrt(1/temp)
						crystalname="Tetragonal"

						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==3) //Orthorhombic
					{
						number temp=(((1/a**2)*h**2)+((1/b**2)*k**2)+((1/c**2)*l**2))
						dspacing=sqrt(1/temp)
						crystalname="Orthorhombic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}
	
					if(crystalvalue==4) //Hexagonal
					{

						number temp=((4/(3*a**2))*(h**2+(h*k)+k**2)+((1/c**2)*l**2))
						dspacing=sqrt(1/temp)
						crystalname="Hexagonal"

						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==5) //Rhombohedral
					{
						number temp1=((h**2+k**2+l**2)*sin(alpha)**2)+(2*((h*k)+(k*l)+(h*l))*(cos(alpha)**2-cos(alpha)))
						number temp2=a**2*(1-(3*cos(alpha)**2)+(2*cos(alpha)**3))
						number temp3=temp1/temp2
						dspacing=sqrt(1/temp3)
						crystalname="Rhombohedral"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}	

					if(crystalvalue==6) //Monoclinic
					{
						number temp=((1/a**2)*(h**2/sin(beta)**2))+((1/b**2)*k**2)+((1/c**2)*(l**2/sin(beta)**2))-((2*h*l*cos(beta))/(a*c*sin(beta)**2))
						dspacing=sqrt(1/temp)
						crystalname="Monoclinic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==7) //Triclinic
					{
						number V2=a**2*b**2*c**2*(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+(2*cos(alpha)*cos(beta)*cos(gamma)))
						number s11=b**2*c**2*sin(alpha)**2
						number s22=a**2*c**2*sin(beta)**2
						number s33=a**2*b**2*sin(gamma)**2

						number s12=a*b*c**2*(cos(alpha)*cos(beta)-cos(gamma))
						number s23=a**2*b*c*(cos(beta)*cos(gamma)-cos(alpha))
						number s31=a*b**2*c*(cos(gamma)*cos(alpha)-cos(beta))

						number temp=(1/v2)*((s11*h**2)+(s22*k**2)+(s33*l**2)+(2*s12*h*k)+(2*s23*k*l)+(2*s31*l*h))
						dspacing=sqrt(1/temp)
						crystalname="Triclinic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

				}		
		}
	}

	// ********* find hkl indices with the given error**************
	 number d_exp=tg.dlggetvalue() 

	 taggroup ed
	 ed=self.LookUpElement("errorfield")
	 number Error_d=ed.dlggetvalue() 

number maxd=d_exp+Error_d
number mind=d_exp-Error_d

ptags=getpersistenttaggroup()	
ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)
number peakno=dspacingtags.taggroupcounttags()

for(number i=1; i<peakno+1; i++)
       {
          number d,h,k,l
		  string mcplabel = dspacingtags.TagGroupGetTaglabel( i-1 )
		  
		  d=val(mcplabel)		  
		  number delta_d=abs(d_exp-d )   
		  		  
		   if(d<maxd && d>mind)
                {
				 taggroup peaktags,peaks
				 peaktags=getpersistenttaggroup()
				 peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d, peaks)
				 number peaknos=peaks.taggroupcounttags()
				 
				String peaklabel = peaks.TagGroupGetTaglabel(0 )
				number peakvalue=val(peaklabel)

				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":h",h)
                getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":k",k)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":l",l)
					
                Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":delta d",delta_d) 
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":d",d)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":h",h)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":k",k)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":l",l)
				}
		}

//order and output it
result("-----------------------------------------------------------------------------\n")

	 taggroup peaktags,peaks
	 peaktags=getpersistenttaggroup()
	 peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:Temp",peaks)
	 no=peaks.taggroupcounttags()
	 for(number i=0;i<no;i++)
	 {
	 String peaklabel = peaks.TagGroupGetTaglabel(i )
	 number peakvalue=val(peaklabel)

	 number dval,delta_d
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":d",dval)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":h",h)
     getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":k",k)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":l",l)	
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":delta d",delta_d) 
	 result(format(dval,"%2.3f")+"("+format(delta_d,"%2.2f")+"), ("+h+","+k+","+l+")"+"\n") 
	 }
	 
	 deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Temp")	
}

//Function: to response when error is changed in d-hkl sub-pannel
void errorfieldChanged(object self, taggroup tg)
{
	number crystalvalue, alpha, beta, gamma, a, b, c, crysttabs, h, k, l, dspacing,no=1
	number maxhklval, mindval
	string crystalname

	//clear tags
	taggroup dspacingtags
    taggroup ptags=getpersistenttaggroup()
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:dspacing")		
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:List1", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:List1")	
	if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:Temp", dspacingtags)) deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Temp")
	
	// read lattice parameters of the selected crystal
	self.dlggetvalue("crystalradios",crystalvalue)
	self.dlggetvalue("alphafield",alpha)
	self.dlggetvalue("betafield",beta)
	self.dlggetvalue("gammafield",gamma)
	self.dlggetvalue("afield",a)
	self.dlggetvalue("bfield",b)
	self.dlggetvalue("cfield",c)
	
		getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:maxhkl",maxhklval)
		getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:mind",mindval)

		alpha=alpha/(180/pi())
		beta=beta/(180/pi())
		gamma=gamma/(180/pi())
		
		taggroup mcptags
	    ptags=getpersistenttaggroup()

		if(a<=0 || b<=0 || c<=0 || alpha<=0 || beta<=0 || gamma<=0 )
			{
				beep()
				exit(0)
			}

	// begin to calculate d-spacings
	for(h=maxhklval; h>=-maxhklval; h--)
	{
		for(k=maxhklval; k>=-maxhklval; k--)
		{
			for(l=maxhklval; l>=-maxhklval; l--)
				{
				//skip for 00
					if(h**2+k**2+l**2==0) continue

					if(crystalvalue==1) // Cubic
					{
						dspacing=sqrt(a**2/(h**2+k**2+l**2))
						crystalname="Cubic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1
					}

					if(crystalvalue==2) // Tetragonal
					{
						number temp=((1/a**2)*(h**2+k**2))+((1/c**2)*l**2)
						dspacing=sqrt(1/temp)
						crystalname="Tetragonal"

						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==3) //Orthorhombic
					{
						number temp=(((1/a**2)*h**2)+((1/b**2)*k**2)+((1/c**2)*l**2))
						dspacing=sqrt(1/temp)
						crystalname="Orthorhombic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}
	
					if(crystalvalue==4) //Hexagonal
					{

						number temp=((4/(3*a**2))*(h**2+(h*k)+k**2)+((1/c**2)*l**2))
						dspacing=sqrt(1/temp)
						crystalname="Hexagonal"

						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==5) //Rhombohedral
					{
						number temp1=((h**2+k**2+l**2)*sin(alpha)**2)+(2*((h*k)+(k*l)+(h*l))*(cos(alpha)**2-cos(alpha)))
						number temp2=a**2*(1-(3*cos(alpha)**2)+(2*cos(alpha)**3))
						number temp3=temp1/temp2
						dspacing=sqrt(1/temp3)
						crystalname="Rhombohedral"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}	

					if(crystalvalue==6) //Monoclinic
					{
						number temp=((1/a**2)*(h**2/sin(beta)**2))+((1/b**2)*k**2)+((1/c**2)*(l**2/sin(beta)**2))-((2*h*l*cos(beta))/(a*c*sin(beta)**2))
						dspacing=sqrt(1/temp)
						crystalname="Monoclinic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

					if(crystalvalue==7) //Triclinic
					{
						number V2=a**2*b**2*c**2*(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+(2*cos(alpha)*cos(beta)*cos(gamma)))
						number s11=b**2*c**2*sin(alpha)**2
						number s22=a**2*c**2*sin(beta)**2
						number s33=a**2*b**2*sin(gamma)**2

						number s12=a*b*c**2*(cos(alpha)*cos(beta)-cos(gamma))
						number s23=a**2*b*c*(cos(beta)*cos(gamma)-cos(alpha))
						number s31=a*b**2*c*(cos(gamma)*cos(alpha)-cos(beta))

						number temp=(1/v2)*((s11*h**2)+(s22*k**2)+(s33*l**2)+(2*s12*h*k)+(2*s23*k*l)+(2*s31*l*h))
						dspacing=sqrt(1/temp)
						crystalname="Triclinic"
						
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":h",h)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":k",k)
						setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+dspacing+":"+no+":l",l)

						no=no+1						
					}

				}		
		}
	}

	// *********search hkl planes within the given d-error**************
	 number Error_d=tg.dlggetvalue() 

	 taggroup dfield
	 dfield=self.LookUpElement("dspacingfield")
	 number d_exp=dfield.dlggetvalue() 

number maxd=d_exp+Error_d
number mind=d_exp-Error_d

ptags=getpersistenttaggroup()	
ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing", dspacingtags)
number peakno=dspacingtags.taggroupcounttags()

for(number i=1; i<peakno+1; i++)
       {
          number d,h,k,l
		  string mcplabel = dspacingtags.TagGroupGetTaglabel( i-1 )
		  
		  d=val(mcplabel)		  
		  number delta_d=abs(d_exp-d )  
		  		  
		   if(d<maxd && d>mind)
                {
				 taggroup peaktags,peaks
				 peaktags=getpersistenttaggroup()
				 peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d, peaks)
				 number peaknos=peaks.taggroupcounttags()
				 
				String peaklabel = peaks.TagGroupGetTaglabel(0 )
				number peakvalue=val(peaklabel)

				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":h",h)
                getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":k",k)
				getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:dspacing:"+d+":"+peakvalue+":l",l)
					
                Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":delta d",delta_d) 
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":d",d)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":h",h)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":k",k)
				Setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+delta_d+":l",l)
				}
		}

//order and output them
result("-----------------------------------------------------------------------------\n")
	 taggroup peaktags,peaks
	 peaktags=getpersistenttaggroup()
	 peaktags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:Temp",peaks)
	 no=peaks.taggroupcounttags()
	 for(number i=0;i<no;i++)
	 {
	 String peaklabel = peaks.TagGroupGetTaglabel(i )
	 number peakvalue=val(peaklabel)

	 number dval,delta_d
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":d",dval)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":h",h)
     getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":k",k)
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":l",l)	
	 getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Temp:"+peakvalue+":delta d",delta_d) 
	 result(format(dval,"%2.4f")+"("+format(delta_d,"%2.1f")+"), ("+h+","+k+","+l+")"+"\n") 
	 }
	 
	 deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Temp")	
}

//peak button in lattice sub-pannel: to measure d-spacing of the selected peak
void PeakButton(object self)
{
image img:=getfrontimage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

string imgType,unit
Getstringnote(img, "ElectronDiffraction Tools:File type",imgType)
unit=img.ImageGetDimensionUnitString(0)

number l,r,d,xscale
if(imgType=="Intensity Profile"||imgType=="Linear Image")
{
ROI currentROI = imgdisp.ImageDisplayGetROI( 0 )
currentROI.roigetrange(l,r)

xscale=img.ImageGetDimensionScale(0)

if(unit=="s (1/nm)"||unit=="1/nm")
{
d=10/(0.5*(l+r)*xscale )
}

if(unit=="Q (1/A)")
{
d=2*Pi()/(0.5*(l+r)*xscale )
}

currentROI.roisetlabel("d="+d+" A")
currentROI.roisetvolatile(0)
currentROI.roisetcolor(1,0,0)

Result("d= "+d+" A\n")

self.lookupElement("dfield").dlgvalue(d)
}
}

//Spot button in lattice sub-pannel: to measure the d-spacing of the selected spots
void SpotButton(object self)
{
Result("-----------------------------------------------------------------------------------------\nSpace bar: to get d-spacing of the selected spot\nALT+: to get d-spacing of the selected spot by ROI \n-----------------------------------------------------------------------------------------\n")	

if(OptionDown())
{
image img:=GetFrontImage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

number t,b,l,r
number annots= imgdisp.ImageDisplayCountROIS()
if(annots==1)
{
ROI currentROI = imgdisp.ImageDisplayGetROI( 0 )
currentROI.roigetrectangle(t,l,b,r)

number x0,y0,x1,y1,d,xscale
xscale=img.ImageGetDimensionScale(0)
getNumberNote(img,"ElectronDiffraction Tools:Center:X", x0)
getNumberNote(img,"ElectronDiffraction Tools:Center:Y", y0)
x1=0.5*(l+r)
y1=0.5*(t+b)

string unit=img.ImageGetDimensionUnitString(0)
if(unit=="s (1/nm)"||unit=="1/nm")
{
d=10/(sqrt((x1-x0)**2+(y1-y0)**2 )*xscale )
}

if(unit=="Q (1/A)")
{
d=2*Pi()/(sqrt((x1-x0)**2+(y1-y0)**2 )*xscale )
}

Result("d= "+d+" A\n")
self.lookupElement("dfield").dlgvalue(d)
}

}

else
{
image   frontimage:=getfrontimage()
ImageDisplay frontimage_disp = frontimage.ImageGetImageDisplay(0)
ImageDocument  frontimage_doc=getfrontimagedocument()
documentwindow frontimage_win=ImageDocumentGetWindow(frontimage_doc)
number x0, y0, x1,y1,x,y,xsize,ysize
frontimage.getsize(xsize,ysize)

string imgname=frontimage.getname()

string unit
unit=frontimage.ImageGetDimensionUnitString(0)

number k=1,m=1,noclick=1

number radius=xsize/150
if(radius<=1)radius=1

while(2>m)
{
number keypress=getkey()

 //find the coordinator of the spot
 number xsize, ysize,xscale,yscale,x, y,x0,y0
 getwindowsize(frontimage, xsize, ysize)
 getscale(frontimage,xscale,yscale)
 
// you can use mouse click to get the end point
number n, mouse_win_x, mouse_win_y
Number img_view_t, img_view_l, img_view_b, img_view_r
Number v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y
Number img_win_t, img_win_l, img_win_b, img_win_r
Number i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y
Number i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y
Number mouse_img_x, mouse_img_y

 windowgetmouseposition(frontimage_win, mouse_win_x, mouse_win_y)
 
//display how to operate
imagedisplay imgdisp=frontimage.imagegetimagedisplay(0)
number annotationID

		frontimage_disp.ImageDisplayGetImageRectInView( img_view_t, img_view_l, img_view_b, img_view_r )

		frontimage_doc.ImageDocumentGetViewToWindowTransform( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y )
		
		ObjectTransformTransformRect( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                            , img_view_t, img_view_l, img_view_b, img_view_r\
                            , img_win_t, img_win_l, img_win_b, img_win_r );

		frontimage_disp.ComponentGetChildToViewTransform( i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y )
		
		ObjectTransformCompose( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                      , i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y\
                      , i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y )
		
		ObjectTransformUntransformPoint( i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y\
                               , mouse_win_x, mouse_win_y, mouse_img_x, mouse_img_y );

		x1=mouse_img_x
		y1=mouse_img_y

frontimage.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)		
		
if(keypress==32) // space bar to mark spot
{
SetNumberNote(frontimage,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":X", x1 )
SetNumberNote(frontimage,"ElectronDiffraction Tools:Live profile:Spot"+noclick+":Y", y1 )

noclick=noclick+1
}		

if(keypress==29||keypress==30) //right and up key
{
radius=radius+0.2*radius
frontimage.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}

if(keypress==28||keypress==31) //left and down key
{
radius=radius-0.2*radius

if(radius<1)radius=1

frontimage.SetSelection(y1-radius,x1-radius,y1+radius,x1+radius)
}

if(keypress>0 && keypress!=32&& keypress!=28&& keypress!=29&& keypress!=30&& keypress!=31||noclick==2) // Any key except space pressed - cancel
{
m=2
}
}

//and then refine it
frontimage.GetSize(xsize,ysize)
number hatwidth=12, brimwidth=6, hatthreshold=0.005,dataselected=0,peakchannel,width=xsize/20,t,l,b,r,xx,yy

for(number i=1;i<2;i++)
{
GetNumberNote(frontimage,"ElectronDiffraction Tools:Live profile:Spot"+i+":X", x1 )
GetNumberNote(frontimage,"ElectronDiffraction Tools:Live profile:Spot"+i+":Y", y1 )

image projectionx = RealImage("", 4, 2*width, 1)
image projectiony = RealImage("", 4, 2*width, 1)

projectionx[icol, 0] += medianfilter(frontimage[y1-width,x1-width,y1+width,x1+width],3,0) 
projectionx=projectionx-projectionx.min()

image mirrorx=projectionx*0
mirrorx=projectionx[iwidth-icol,irow] 

number maxx,maxy  
number maxval=projectionx[0,width-0.2*width,1,width+0.2*width].max(maxx,maxy)
maxx=width-0.2*width+maxx
FWHMcalc(projectionx, maxx)
projectionx.GetSelection(t,l,b,r)

tophatfilter(projectionx, brimwidth, hatwidth, hatthreshold, dataselected)
getpersistentnumbernote("ElectronDiffraction Tools:maxpeak:Peak1:X",x)
x=x1-width+x

maxval=mirrorx[0,width-0.2*width,1,width+0.2*width].max(maxx,maxy)
maxx=width-0.2*width+maxx
FWHMcalc(mirrorx, maxx)
mirrorx.GetSelection(t,l,b,r)

tophatfilter(mirrorx, brimwidth, hatwidth, hatthreshold, dataselected)
getpersistentnumbernote("ElectronDiffraction Tools:maxpeak:Peak1:X",xx)
xx=x1-width+(2*width-xx)
x=0.5*(x+xx)

projectiony[irow, 0] += medianfilter(frontimage[y1-width,x1-width,y1+width,x1+width],3,0)
projectiony=projectiony-projectiony.min()

image mirrory=projectiony*0
mirrory=projectiony[iwidth-icol,irow] 

maxval=projectiony[0,width-0.2*width,1,width+0.2*width].max(maxx,maxy)
maxx=width-0.2*width+maxx
FWHMcalc(projectiony, maxx)
projectiony.GetSelection(t,l,b,r)

tophatfilter(projectiony, brimwidth, hatwidth, hatthreshold, dataselected)
getpersistentnumbernote("ElectronDiffraction Tools:maxpeak:Peak1:X",y)
y=y1-width+y

maxval=mirrory[0,width-0.2*width,1,width+0.2*width].max(maxx,maxy)
maxx=width-0.2*width+maxx
FWHMcalc(mirrory, maxx)
mirrory.GetSelection(t,l,b,r)

tophatfilter(mirrory, brimwidth, hatwidth, hatthreshold, dataselected)
getpersistentnumbernote("ElectronDiffraction Tools:maxpeak:Peak1:X",yy)
yy=y1-width+(2*width-yy)
y=0.5*(y+yy)

setnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot"+i+":X",x)
setnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot"+i+":Y",y)

Component box=NewBoxAnnotation(y-1,x-1,y+1,x+1)
box.componentsetfillmode(2)
box.componentsetforegroundcolor(1,0,0) // sets the foreground colour to magenta
frontimage_disp.componentaddchildatend(box)
}

number xscale,yscale,d
frontimage.GetScale(xscale,yscale)

getNumberNote(frontimage,"ElectronDiffraction Tools:Center:X", x0)
getNumberNote(frontimage,"ElectronDiffraction Tools:Center:Y", y0)
Getnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot1:X",x1)
Getnumbernote(frontimage,"ElectronDiffraction Tools:Spot coordinates:Spot1:Y",y1)
deletenote(frontimage,"ElectronDiffraction Tools:Live profile") 
deletenote(frontimage,"ElectronDiffraction Tools:Spot coordinates")

if(unit=="s (1/nm)"||unit=="1/nm")
{
d=10/(sqrt((x1-x0)**2+(y1-y0)**2 )*xscale )
}

if(unit=="Q (1/A)")
{
d=2*Pi()/(sqrt((x1-x0)**2+(y1-y0)**2 )*xscale )
}

Result("d= "+d+" A\n")
self.lookupElement("dfield").dlgvalue(d)
}
}

//Find lattice button in lattice sub-pannel
void FindLatticeButton(object self)
{


}


//dfield changed in lattice sub-pannel
void dfieldChanged(object self, taggroup tg)
{

}

//lattice field changed in lattice sub-pannel
void Latticehfieldchanged(object self, taggroup tg)
{

}


//function: [uvw]>> button in Indices sub-pannel
void RunButton1( object self ) 
{
TagGroup tag=self.LookUpElement("Indfield1") 
string thisline=tag.dlggetstringvalue()

//pick out u,v,w
number width=len(thisline)
for(number i=0;i<width;i++)
{
string thischar=mid(thisline,i,1)

number Bytes,notab=1
if(thischar==",")  
{
Bytes=len(left(thisline,i))/100
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:temp:Tabs:"+Bytes,Bytes)
}
}

	taggroup abtags,mcptags
	taggroup ptags=getpersistenttaggroup()
	ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:Read Files:temp:Tabs", mcptags)
	
	number mcpno=mcptags.taggroupcounttags()
	for(number i=1; i<mcpno+1; i++)
		{
		String mcplabel = mcptags.TagGroupGetTaglabel( i-1 )
		number value=val(mcplabel)*100
		setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab "+i,value)
		}
		
deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Read Files:temp")

number Tab1,tab2,tab3,u,v,w,t
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 1",Tab1) 
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 2",Tab2)  

u=val(left(thisline,tab1))
v=val(mid(thisline,tab1+1,tab2-tab1-1))
w=val(right(thisline,width-tab2-1))
deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Read Files")

//convert three indices to four indices
number u1,v1,t1,w1
u1=(2*u-v)/3
v1=(2*v-u)/3
t1=-(u+v)/3
w1=w

image img=RealImage("",4,1,4)
if(abs(u1)>0)img.SetPixel(0,0,abs(u1))
else img.SetPixel(0,0,1e6)

if(abs(v1)>0)img.SetPixel(0,1,abs(v1))
else img.SetPixel(0,1,1e6)

if(abs(t1)>0)img.SetPixel(0,2,abs(t1))
else img.SetPixel(0,2,1e6)

if(abs(w1)>0)img.SetPixel(0,3,abs(u1))
else img.SetPixel(0,3,1e6)

number minval=img.min()
u1=u1/minval
v1=v1/minval
t1=t1/minval
w1=w1/minval

Result("\n[ "+u+", "+v+", "+w +"] ==> [ "+u1+", "+v1+", "+t1+", "+w1+" ]\n")
}

//function: (hkl)>> button in Indices sub-pannel
void RunButton2( object self ) 
{
TagGroup tag=self.LookUpElement("Indfield2") 
string thisline=tag.dlggetstringvalue()

//pick out hkl
number width=len(thisline)
for(number i=0;i<width;i++)
{
string thischar=mid(thisline,i,1)

number Bytes,notab=1
if(thischar==",") 
{
Bytes=len(left(thisline,i))/100
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:temp:Tabs:"+Bytes,Bytes)
}
}

	taggroup abtags,mcptags
	taggroup ptags=getpersistenttaggroup()
	ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:Read Files:temp:Tabs", mcptags) 
	
	number mcpno=mcptags.taggroupcounttags()
	for(number i=1; i<mcpno+1; i++)
		{
		String mcplabel = mcptags.TagGroupGetTaglabel( i-1 )
		number value=val(mcplabel)*100
		setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab "+i,value)
		}
		
deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Read Files:temp")

number Tab1,tab2,tab3,u,v,w,t
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 1",Tab1) 
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 2",Tab2)  

u=val(left(thisline,tab1))
v=val(mid(thisline,tab1+1,tab2-tab1-1))
w=val(right(thisline,width-tab2-1))
deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Read Files")

//convert hkl to hkil
number u1,v1,t1,w1
u1=u
v1=v
t1=-(u+v)
w1=w

Result("\n( "+u+", "+v+", "+w +") ==> ( "+u1+", "+v1+", "+t1+", "+w1+" )\n")
}

//function: [uvtw]>> button in Indices sub-pannel
void RunButton3( object self ) 
{
TagGroup tag=self.LookUpElement("Indfield3") 
string thisline=tag.dlggetstringvalue()

//pick out uvtw
number width=len(thisline)
for(number i=0;i<width;i++)
{
string thischar=mid(thisline,i,1)

number Bytes,notab=1
if(thischar==",")  
{
Bytes=len(left(thisline,i))/100
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:temp:Tabs:"+Bytes,Bytes)
}
}

	taggroup abtags,mcptags
	taggroup ptags=getpersistenttaggroup()
	ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:Read Files:temp:Tabs", mcptags) 
	
	number mcpno=mcptags.taggroupcounttags()
	for(number i=1; i<mcpno+1; i++)
		{
		String mcplabel = mcptags.TagGroupGetTaglabel( i-1 )
		number value=val(mcplabel)*100
		setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab "+i,value)
		}
		
deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Read Files:temp")

number Tab1,tab2,tab3,u,v,w,t
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 1",Tab1) 
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 2",Tab2)  
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 3",Tab3)  

u=val(left(thisline,tab1))
v=val(mid(thisline,tab1+1,tab2-tab1-1))
t=val(mid(thisline,tab2+1,tab3-tab1-1))
w=val(right(thisline,width-tab3-1))
deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Read Files")

//conver uvtw to uvw
number u1,v1,w1
u1=2*u+v
v1=2*v+u
w1=w

Result("\n[ "+u+", "+v+", "+t+", "+w +"] ==> [ "+u1+", "+v1+", "+w1+" ]\n")
}

//function: [uvtw]>> button in Indices sub-pannel
void RunButton4( object self ) 
{
TagGroup tag=self.LookUpElement("Indfield4") 
string thisline=tag.dlggetstringvalue()

//pick out hkil
number width=len(thisline)
for(number i=0;i<width;i++)
{
string thischar=mid(thisline,i,1)

number Bytes,notab=1
if(thischar==",")  
{
Bytes=len(left(thisline,i))/100
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:temp:Tabs:"+Bytes,Bytes)
}
}

	taggroup abtags,mcptags
	taggroup ptags=getpersistenttaggroup()
	ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:Crystal Tools:Read Files:temp:Tabs", mcptags)
	
	number mcpno=mcptags.taggroupcounttags()
	for(number i=1; i<mcpno+1; i++)
		{
		String mcplabel = mcptags.TagGroupGetTaglabel( i-1 )
		number value=val(mcplabel)*100
		setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab "+i,value)
		}
		
deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Read Files:temp")

number Tab1,tab2,tab3,u,v,w,t
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 1",Tab1) 
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 2",Tab2)  
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:PDF file:Tabs:Tab 3",Tab3)  

u=val(left(thisline,tab1))
v=val(mid(thisline,tab1+1,tab2-tab1-1))
t=val(mid(thisline,tab2+1,tab3-tab1-1))
w=val(right(thisline,width-tab3-1))
deletepersistentnote("ElectronDiffraction Tools:Crystal Tools:Read Files")

//hkil convert to hkl
number u1,v1,w1
u1=u
v1=v
w1=w

Result("\n( "+u+", "+v+", "+t+", "+w +") ==> ( "+u1+", "+v1+", "+w1+" )\n")
}


	Void Action(Object self) 
		{
		TagGroup tg = self.LookUpElement("TabList")
		}
		
	Void TabChangeAction(Object self, TagGroup TabListSelf)
		{
		}
		

//---------------------------------Create crystal tools GUI-------------------------------		
TagGroup CreateCrystalTools(Object self)
	{
	TagGroup CrystalTools = DLGCreateTab("Crystal Tools")
		
//----------------- create "Create Lattice box" -------------------------
	TagGroup CreateLatticeBoxs, CreateLattice
	CreateLatticeBoxs = DLGCreateBox("Lattice",CreateLattice).dlgexternalpadding(2,0).dlginternalpadding(2,0)

//create crystal, pull-down lsit, +, Name...
	TagGroup CrystalLable=DLGCreateLabel("Crystal:").dlganchor("West").dlgexternalpadding(0,0)	
	TagGroup CNPopUp=MakeElementPullDown()
	CNPopUp.dlgidentifier("elementpulldown").dlgchangedmethod("elementchanged").dlgexternalpadding(0,0)
	
	TagGroup CNadd=DLGCreatePushButton("+","AddCrystal")
	
	TagGroup CNLable=DLGCreateLabel("Name:").dlganchor("West").dlgexternalpadding(0,0)
	taggroup crystalnamefield=dlgcreatestringfield("",13).dlgidentifier("crystalnamefield").dlgexternalpadding(0,0).dlganchor("West")	
	TagGroup CN=DLGGroupItems(CNLable,crystalnamefield).DLGtablelayout(2,1,0)
	
	TagGroup Crystalname=DLGGroupItems(CrystalLable,CNPopUp,CNadd,CN).DLGtablelayout(4,1,0)

//create Crystal Type radio list
	TagGroup CTLable=DLGCreateLabel("    Type:")
	Taggroup crystalradios=MakeCrystTypeRadios()  
	crystalradios.dlgidentifier("crystalradios").dlgchangedmethod("crystalchanged")
	
	TagGroup CrystalType=DLGGroupItems(CTLable,crystalradios).DLGtablelayout(2,1,0)
	TagGroup CrystalNT=DLGGroupItems(CrystalName,CrystalType).DLGtablelayout(1,2,0)
	
//create unit cell box
	TagGroup Unitcell=DLGCreateLabel("Unit cell ").dlganchor("West")
	
	TagGroup alabel=DLGCreateLabel("a= ")
	TagGroup afield = DLGCreateRealField(0, 13,7).dlgidentifier("afield").dlgchangedmethod("latticechanged")
    taggroup agroup=dlggroupitems(alabel, afield).dlgtablelayout(2,1,0)
	
    taggroup alphalabel=dlgcreatelabel("α=")
    TagGroup alphafield= DLGCreateRealField(0, 13, 0).dlgidentifier("alphafield").dlgchangedmethod("alphachanged")
    taggroup alphagroup=dlggroupitems(alphalabel, alphafield).dlgtablelayout(2,1,0)	
	taggroup a_alphagroup=dlggroupitems(agroup, alphagroup).dlgtablelayout(1,2,0)
	
	TagGroup blabel=DLGCreateLabel("b= ")
	TagGroup bfield = DLGCreateRealField(0, 13,7).dlgidentifier("bfield")
    taggroup bgroup=dlggroupitems(blabel, bfield).dlgtablelayout(2,1,0)	
	
	taggroup betalabel=dlgcreatelabel("β=")
    TagGroup betafield = DLGCreateRealField(0, 13, 0).dlgidentifier("betafield")
    taggroup betagroup=dlggroupitems(betalabel, betafield).dlgtablelayout(2,1,0)
	taggroup b_betagroup=dlggroupitems(bgroup, betagroup).dlgtablelayout(1,2,0)
	
	TagGroup clabel=DLGCreateLabel("c= ")
	TagGroup cfield = DLGCreateRealField(0, 13,7).dlgidentifier("cfield")
    taggroup cgroup=dlggroupitems(clabel, cfield).dlgtablelayout(2,1,0)
	
	taggroup gammalabel=dlgcreatelabel("γ=")
    TagGroup gammafield = DLGCreateRealField(0, 13, 0).dlgidentifier("gammafield")
    taggroup gammagroup=dlggroupitems(gammalabel, gammafield).dlgtablelayout(2,1,0)
	taggroup c_gammagroup=dlggroupitems(cgroup, gammagroup).dlgtablelayout(1,2,0)
	
	taggroup unitsizegroup=dlggroupitems(a_alphagroup,b_betagroup,c_gammagroup).dlgtablelayout(3,1,0)
	taggroup unitcellgroup=dlggroupitems(unitcell,unitsizegroup).dlgtablelayout(1,2,0)
	taggroup latticegroup=dlggroupitems(CrystalNT,unitcellgroup).dlgtablelayout(1,2,0).dlganchor("West")
	
	CreateLatticeBoxs.DLGAddElement(latticegroup)
	
//----------------------------create indexing box----------------------------
	TagGroup CreateIndexingBoxs, IndexingItems
	CreateIndexingBoxs = DLGCreateBox("Indexing",IndexingItems).dlgexternalpadding(2,0 ).dlginternalpadding(2,0)

	//1) create o-O-o and calc buttons
	TagGroup CalcBoxs, CalcItems
	CalcBoxs = DLGCreateBox("",CalcItems).dlgexternalpadding(0,5).dlginternalpadding(0,2)

	CalcItems.DLGAddElement(DLGCreatePushButton("o-O-o","d_angle").dlginternalpadding(0,0).dlgexternalpadding(4,3))	
	CalcItems.DLGAddElement(DLGCreatePushButton("Calc.","Calcbutton").dlginternalpadding(0,0).dlgexternalpadding(4,3))	

	number spotsize //set spot size 
	getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Spot size",spotsize)
	
	TagGroup sizelabel=DLGCreateLabel("Size:") 
	TagGroup Sizefield = DLGCreateIntegerField(spotsize, 4).dlgidentifier("sizefield")
	Taggroup sizegroup=dlggroupitems(sizelabel, Sizefield).dlgtablelayout(2,1,0)
	CalcItems.DLGAddElement(sizegroup)
	
	TagGroup RGBlabel=DLGCreateLabel("RGB:") //set RGB
	TagGroup Rfield = DLGCreateRealField(0, 4,1).dlgidentifier("Rfield")
	TagGroup Gfield = DLGCreateRealField(0.8, 4,1).dlgidentifier("Gfield")
	TagGroup B1field = DLGCreateRealField(0, 4,1).dlgidentifier("B1field")
	Taggroup RGBgroup=dlggroupitems(RGBlabel, Rfield,Gfield,B1field).dlgtablelayout(4,1,0)
	CalcItems.DLGAddElement(RGBgroup)
	
	CalcBoxs.DLGTableLayOut(4,1,0)

	//2) display d1, d2 and angle
	TagGroup d1label=DLGCreateLabel("d1= ")
	Taggroup d1field=DLGCreateLabel( "" + format(0, "%6.4f") + "  A").DLGIdentifier("d1field").dlgexternalpadding(0,0)
    Taggroup d1group=dlggroupitems(d1label, d1field).dlgtablelayout(2,1,0).dlginternalpadding(2,0).dlgexternalpadding(2,2)
	
	TagGroup d2label=DLGCreateLabel("d2= ")
	TagGroup d2field =DLGCreateLabel( "" + format(0, "%6.4f") + "  A").DLGIdentifier("d2field").dlgexternalpadding(0,0)
    Taggroup d2group=dlggroupitems(d2label, d2field).dlgtablelayout(2,1,0).dlginternalpadding(2,0).dlgexternalpadding(2,2)
	
	TagGroup Anglelabel=DLGCreateLabel(" Angle= ")
	TagGroup  Anglefield=DLGCreateLabel( "" + format(0, "%6.2f") + "  deg.").DLGIdentifier("Anglefield").dlgexternalpadding(0,0)
    Taggroup  Anglegroup=dlggroupitems( Anglelabel,  Anglefield).dlgtablelayout(2,1,0).dlginternalpadding(2,0).dlgexternalpadding(2,2)
	
	taggroup d_Anglegroup=dlggroupitems(d1group,d2group,Anglegroup).dlgtablelayout(3,1,0).dlganchor("West")	
	
	//3) (hkl)1, (hkl)2, [uvw] lists
	TagGroup hklLable1=DLGCreateLabel("(hkl)1")
	TagGroup hkllist1,hkllistItems1
	hkllist1=DLGCreatelist(20,5).DLGIdentifier("hklList1").dlgchangedmethod("hklList1Changed")
	TagGroup hklgroup1=DLGGroupItems(hklLable1,hkllist1).DLGtablelayout(1,2,0)
	
	TagGroup hklLable2=DLGCreateLabel("(hkl)2")
	TagGroup hkllist2,hkllistItems2
	hkllist2=DLGCreatelist(22,5).DLGIdentifier("hklList2").dlgchangedmethod("hklList2Changed")
	TagGroup hklgroup2=DLGGroupItems(hklLable2,hkllist2).DLGtablelayout(1,2,0)
	
	TagGroup uvwLable=DLGCreateLabel("[uvw]")
	TagGroup uvwlist,uvwlistItems
	uvwlist=DLGCreatelist(12,5).DLGIdentifier("uvwList").dlgchangedmethod("ZoneListChanged")
	TagGroup hklgroup3=DLGGroupItems(uvwLable,uvwlist).DLGtablelayout(1,2,0)
	
	TagGroup hklgroup=DLGGroupItems(hklgroup1,hklgroup2,hklgroup3).DLGtablelayout(3,1,0)
	
	TagGroup Analysisgroup=DLGGroupItems(CalcBoxs,d_Anglegroup,hklgroup).DLGtablelayout(1,3,0)
	CreateIndexingBoxs.DLGAddElement(Analysisgroup)

//---------------Create profile box-----------------------	
TagGroup CreateProfileBoxs,CreateProfile
CreateProfileBoxs = DLGCreateBox("Profile",CreateProfile).dlgexternalpadding(2,2 ).dlginternalpadding(6,6)	

CreateProfile.DLGAddElement(DLGCreatePushButton("SAED","SAEDButton").dlginternalpadding(0,0).dlgexternalpadding(0,0))	
CreateProfile.DLGAddElement(DLGCreatePushButton("Center","CenterButton").dlginternalpadding(0,0).dlgexternalpadding(0,0))	
CreateProfile.DLGAddElement(DLGCreatePushButton("Profile","ProfileButton").dlginternalpadding(0,0).dlgexternalpadding(0,0))	
CreateProfile.DLGAddElement(DLGCreatePushButton("Peaks","PeaksButton").dlginternalpadding(0,0).dlgexternalpadding(0,0))	
CreateProfile.DLGAddElement(DLGCreatePushButton("Calib.","CalibButton").dlginternalpadding(0,0).dlgexternalpadding(0,0))	
CreateProfile.DLGAddElement(DLGCreatePushButton("Save","SaveButton").dlginternalpadding(0,0).dlgexternalpadding(0,0))	

CreateProfileBoxs.DLGtablelayout(6,1,0)

	//---------------Create Lattice Calculator box-----------------------
	TagGroup CreateCalculatorBoxs,CreateCalculator
	CreateCalculatorBoxs = DLGCreateBox("Lattice Calculator",CreateCalculator).dlgexternalpadding(2,2 ).dlginternalpadding(2,0)	
	
	//create sub-pannels
	TagGroup cryst_items
	TagGroup crysttablist = DLGCreateTabList(cryst_items, 0).DLGIdentifier("crysttabs").dlgchangedmethod("crysttabchanged")
	crysttablist.dlgexternalpadding(0,0)

	TagGroup dspacingstab, anglestab, directionstab, zonestab,LatticeTab,Indicestab
	cryst_items.DLGAddTab( DLGCreateTab(dspacingstab, "d-hkl") );
	cryst_items.DLGAddTab(DLGCreateTab(anglestab, "Planes") );
	cryst_items.DLGAddTab(DLGCreateTab(directionstab, "Directions") );
	cryst_items.DLGAddTab(DLGCreateTab(Indicestab, "Indices") );
	cryst_items.DLGAddTab(DLGCreateTab(zonestab, "h <-> u") );
	cryst_items.DLGAddTab(DLGCreateTab(LatticeTab, "Lattice") )
	
//---------------1. d-hkl sub-pannel---------------
taggroup dcalcbox_items
taggroup dcalcbox=dlgcreatebox("",dcalcbox_items)
	
taggroup dhlabel=dlgcreatelabel("( h k l )  ")
TagGroup dhfield = DLGCreateRealField(0, 6,2).dlgidentifier("dhfield").dlgchangedmethod("dhorkchanged")
taggroup dhgroup=dlggroupitems(dhlabel, dhfield).dlgtablelayout(2,1,0)
TagGroup dkfield = DLGCreateRealField(0, 6,2).dlgidentifier("dkfield").dlgchangedmethod("dhorkchanged")
TagGroup dlfield = DLGCreateRealField(0, 6,2).dlgidentifier("dlfield")
taggroup dhkilgroup=dlggroupitems(dhgroup, dkfield, dlfield).dlgtablelayout(3,1,0)

taggroup dspacinglabel=dlgcreatelabel("d:")
taggroup dspacingfield=DLGCreateRealField(0, 9,4).DLGIdentifier("dspacingfield").dlgchangedmethod("dspacingfieldChanged")
TagGroup errorfield=DLGCreateRealField(0.1, 5,2).DLGIdentifier("errorfield").dlgchangedmethod("errorfieldChanged")
taggroup dspacinglabel2=dlgcreatelabel(" A")
taggroup dspacinggroup=dlggroupitems(dspacinglabel, dspacingfield,errorfield,dspacinglabel2).dlgtablelayout(4,1,0).dlganchor("South")

taggroup dspacingtabgroup=dlggroupitems(dhkilgroup, dspacinggroup).dlgtablelayout(1,2,0).dlgexternalpadding(2,0)
dcalcbox_items.dlgaddelement(dspacingtabgroup)

dspacingstab.dlgaddelement(dcalcbox)

taggroup crystcalcbutton=dlgcreatepushbutton("Calc", "maincalcroutine")
taggroup crystclearbutton=dlgcreatepushbutton("Clear", "crystclearbutton")
taggroup crystbuttonsgroup=dlggroupitems( crystcalcbutton,crystclearbutton).dlgtablelayout(2,1,0)

dspacingstab.dlgaddelement(crystbuttonsgroup)

//---------------2. Planes sub-pannel ---------------
taggroup anglesbox_items
taggroup anglesbox=dlgcreatebox("", anglesbox_items)

taggroup angh1label=dlgcreatelabel("( h k l )1  ")
TagGroup angh1field = DLGCreateRealField(0, 5,2).dlgidentifier("angh1field").dlgchangedmethod("anghork1changed")
taggroup angh1group=dlggroupitems(angh1label, angh1field).dlgtablelayout(2,1,0)
TagGroup angk1field = DLGCreateRealField(0, 5,2).dlgidentifier("angk1field").dlgchangedmethod("anghork1changed")
TagGroup angl1field = DLGCreateRealField(0, 5,2).dlgidentifier("angl1field")
taggroup anghkil1group=dlggroupitems(angh1group, angk1field, angl1field).dlgtablelayout(3,1,0)

taggroup angh2label=dlgcreatelabel("( h k l )2  ")
TagGroup angh2field = DLGCreateRealField(0, 5,2).dlgidentifier("angh2field").dlgchangedmethod("anghork2changed")
taggroup angh2group=dlggroupitems(angh2label, angh2field).dlgtablelayout(2,1,0)
TagGroup angk2field = DLGCreateRealField(0, 5,2).dlgidentifier("angk2field").dlgchangedmethod("anghork2changed")
TagGroup angl2field = DLGCreateRealField(0, 5,2).dlgidentifier("angl2field")
taggroup anghkil2group=dlggroupitems(angh2group, angk2field, angl2field).dlgtablelayout(3,1,0)

taggroup anghkil12group=dlggroupitems(anghkil1group, anghkil2group).dlgtablelayout(1,2,0)

taggroup angleresultlabel=dlgcreatelabel("Angle:")
taggroup angleresultfield=DLGCreateLabel( "" + format(0, "%6.1f") + " deg.").DLGIdentifier("angleresultfield")
taggroup angleresultgroup1=dlggroupitems(angleresultlabel, angleresultfield).dlgtablelayout(2,1,0)
taggroup angleresultlabel1=dlgcreatelabel("Zone:")
taggroup angleresultfield1=DLGCreateLabel( "                      " ).DLGIdentifier("angleresultfield1")
taggroup angleresultgroup2=dlggroupitems(angleresultlabel1, angleresultfield1).dlgtablelayout(2,1,0).dlganchor("West")
taggroup angleresultgroup=dlggroupitems(angleresultgroup1, angleresultgroup2).dlgtablelayout(1,2,0)

taggroup anglefinalgroup=dlggroupitems(anghkil12group, angleresultgroup).dlgtablelayout(2,1,0)

anglesbox_items.dlgaddelement(anglefinalgroup)
anglestab.dlgaddelement(anglesbox)
anglestab.dlgaddelement(crystbuttonsgroup)

//--------------- 3. Directions sub-pannel  ---------------
taggroup directionsbox_items
taggroup directionsbox=dlgcreatebox("", directionsbox_items)

taggroup diru1label=dlgcreatelabel("[ u v w ]1 ")
TagGroup diru1field = DLGCreateRealField(0, 5,2).dlgidentifier("diru1field").dlgchangedmethod("diru1orv1changed")
taggroup diru1group=dlggroupitems(diru1label, diru1field).dlgtablelayout(2,1,0)
TagGroup dirv1field = DLGCreateRealField(0, 5,2).dlgidentifier("dirv1field").dlgchangedmethod("diru1orv1changed")
TagGroup dirw1field = DLGCreateRealField(0, 5,2).dlgidentifier("dirw1field")
taggroup diruvtw1group=dlggroupitems(diru1group, dirv1field, dirw1field).dlgtablelayout(3,1,0)

taggroup diru2label=dlgcreatelabel("[ u v w ]2 ")
TagGroup diru2field = DLGCreateRealField(0, 5,2).dlgidentifier("diru2field").dlgchangedmethod("diru2orv2changed")
taggroup diru2group=dlggroupitems(diru2label, diru2field).dlgtablelayout(2,1,0)
TagGroup dirv2field = DLGCreateRealField(0, 5,2).dlgidentifier("dirv2field").dlgchangedmethod("diru2orv2changed")
TagGroup dirw2field = DLGCreateRealField(0, 5,2).dlgidentifier("dirw2field")
taggroup diruvtw2group=dlggroupitems(diru2group, dirv2field, dirw2field).dlgtablelayout(3,1,0)

taggroup diruvtw12group=dlggroupitems(diruvtw1group,diruvtw2group).dlgtablelayout(1,2,0)

taggroup dirresultlabel=dlgcreatelabel("Angle:")
taggroup dirresultfield=DLGCreateLabel( "" + format(0, "%6.1f") + " deg.").DLGIdentifier("dirresultfield")
taggroup dirresultgroup1=dlggroupitems(dirresultlabel, dirresultfield).dlgtablelayout(2,1,0).dlganchor("West")
taggroup dirresultlabel1=dlgcreatelabel("Plane:")
taggroup dirresultfield1=DLGCreateLabel( "                      " ).DLGIdentifier("dirresultfield1")
taggroup dirresultgroup2=dlggroupitems(dirresultlabel1, dirresultfield1).dlgtablelayout(2,1,0).dlganchor("West")
taggroup dirresultgroup=dlggroupitems(dirresultgroup1, dirresultgroup2).dlgtablelayout(1,2,0)

taggroup dirfinalgroup=dlggroupitems(diruvtw12group, dirresultgroup).dlgtablelayout(2,1,0)

directionsbox_items.dlgaddelement(dirfinalgroup)
directionstab.dlgaddelement(directionsbox)
directionstab.dlgaddelement(crystbuttonsgroup)

//--------------- 4. Indices sub-pannel   ---------------
taggroup Indicesbox_items
taggroup Indicesbox=dlgcreatebox("", Indicesbox_items)

taggroup Indlabel1=dlgcreatelabel("[ u v w ] ")
TagGroup Indfield1 = DLGCreateStringField("0,0,0", 9).dlgidentifier("Indfield1")
TagGroup runButton1=DLGCreatePushButton(">>","RunButton1")
taggroup Indgroup1=dlggroupitems(Indlabel1, Indfield1,runButton1).dlgtablelayout(3,1,0)

taggroup Indlabel2=dlgcreatelabel(" ( h k l ) ")
TagGroup Indfield2 = DLGCreateStringField("0,0,0", 9).dlgidentifier("Indfield2")
TagGroup runButton2=DLGCreatePushButton(">>","RunButton2")
taggroup Indgroup2=dlggroupitems(Indlabel2, Indfield2,RunButton2).dlgtablelayout(3,1,0)
taggroup Indgroup12=dlggroupitems(Indgroup1, Indgroup2).dlgtablelayout(1,2,0)

taggroup Indlabel3=dlgcreatelabel("[ u v t w ] ")
TagGroup Indfield3 = DLGCreateStringField("0,0,0,0", 9).dlgidentifier("Indfield3")
TagGroup runButton3=DLGCreatePushButton(">>","RunButton3")
taggroup Indgroup3=dlggroupitems(Indlabel3, Indfield3,runButton3).dlgtablelayout(3,1,0)

taggroup Indlabel4=dlgcreatelabel("( h k i l ) ")
TagGroup Indfield4 = DLGCreateStringField("0,0,0,0", 9).dlgidentifier("Indfield4")
TagGroup runButton4=DLGCreatePushButton(">>","RunButton4")
taggroup Indgroup4=dlggroupitems(Indlabel4, Indfield4,runButton4).dlgtablelayout(3,1,0)
taggroup Indgroup34=dlggroupitems(Indgroup3, Indgroup4).dlgtablelayout(1,2,0)

taggroup Indgroup=dlggroupitems(Indgroup12, Indgroup34).dlgtablelayout(2,1,0)

Indicesbox_items.dlgaddelement(Indgroup)
Indicestab.dlgaddelement(Indicesbox)

// ---------------  5. h<->u sub-pannel  ---------------
taggroup zonebox_items
taggroup zonebox=dlgcreatebox("", zonebox_items)

taggroup zoneh1label=dlgcreatelabel("( h k l )")
TagGroup zoneh1field = DLGCreateRealField(0, 5,2).dlgidentifier("zoneh1field").dlgchangedmethod("zonehork1changed")
taggroup zoneh1group=dlggroupitems(zoneh1label, zoneh1field).dlgtablelayout(2,1,0)
TagGroup zonek1field = DLGCreateRealField(0, 5,2).dlgidentifier("zonek1field").dlgchangedmethod("zonehork1changed")
TagGroup zonel1field =DLGCreateRealField(0, 5,2).dlgidentifier("zonel1field")
taggroup zonehkil1group=dlggroupitems(zoneh1group, zonek1field, zonel1field).dlgtablelayout(3,1,0)

taggroup zoneh2label=dlgcreatelabel("[u v w]")
TagGroup zoneh2field = DLGCreateRealField(0, 5,2).dlgidentifier("zoneh2field").dlgchangedmethod("zonehork2changed")
taggroup zoneh2group=dlggroupitems(zoneh2label, zoneh2field).dlgtablelayout(2,1,0)
TagGroup zonek2field = DLGCreateRealField(0, 5,2).dlgidentifier("zonek2field").dlgchangedmethod("zonehork2changed")
TagGroup zonel2field = DLGCreateRealField(0, 5,2).dlgidentifier("zonel2field")
taggroup zonehkil2group=dlggroupitems(zoneh2group, zonek2field, zonel2field).dlgtablelayout(3,1,0)

taggroup zonehkil12group=dlggroupitems(zonehkil1group,zonehkil2group).dlgtablelayout(1,2,0)

taggroup zoneresultlabel=dlgcreatelabel("[u v w]")
taggroup zoneresultfield=DLGCreateLabel("    0,    0,    0 " ).DLGIdentifier("zoneresultfield")
taggroup zoneresultgroup=dlggroupitems(zoneresultlabel, zoneresultfield).dlgtablelayout(2,1,0).dlganchor("West")
taggroup zoneresultlabel1=dlgcreatelabel("( h k l )")
taggroup zoneresultfield1=DLGCreateLabel("    0,    0,    0 " ).DLGIdentifier("zoneresultfield1")
taggroup zoneresultgroup1=dlggroupitems(zoneresultlabel1, zoneresultfield1).dlgtablelayout(2,1,0).dlganchor("West")

taggroup zoneresultgroup2=dlggroupitems(zoneresultgroup, zoneresultgroup1).dlgtablelayout(1,2,0)

taggroup zonefinalgroup=dlggroupitems(zonehkil12group, zoneresultgroup2).dlgtablelayout(2,1,0)

zonebox_items.dlgaddelement(zonefinalgroup)
zonestab.dlgaddelement(zonebox)
zonestab.dlgaddelement(crystbuttonsgroup)

// ---------------  6. Lattice sub-pannel   ---------------  
taggroup Latticebox_items
taggroup Latticebox=dlgcreatebox("", Latticebox_items)

TagGroup peakbutton=DLGCreatePushButton("Peak","PeakButton").dlginternalpadding(2,0).dlgexternalpadding(2,0)
TagGroup Spotbutton=DLGCreatePushButton("Spot","SpotButton").dlginternalpadding(2,0).dlgexternalpadding(0,0)
TagGroup FindLatticeButton=DLGCreatePushButton("Calc.","FindLatticeButton").dlginternalpadding(2,0).dlgexternalpadding(2,0)
taggroup Buttongroup=dlggroupitems(peakbutton, Spotbutton,FindLatticeButton).dlgtablelayout(3,1,0)
Latticebox_items.dlgaddelement(Buttongroup)

taggroup dlabel=dlgcreatelabel("d (A):")
taggroup dfield=DLGCreateRealField(0, 9,5).DLGIdentifier("dfield").dlgchangedmethod("dfieldChanged")
taggroup dgroup=dlggroupitems(dlabel,dfield).dlgtablelayout(2,1,0)

taggroup Latticelabel=dlgcreatelabel("( h k l )")
TagGroup Latticehfield = DLGCreateRealField(0, 5,2).dlgidentifier("Latticehfield").dlgchangedmethod("Latticehfieldchanged")
TagGroup Latticekfield = DLGCreateRealField(0, 5,2).dlgidentifier("Latticekfield").dlgchangedmethod("Latticehfieldchanged")
TagGroup Latticelfield =DLGCreateRealField(0, 5,2).dlgidentifier("Latticelfield").dlgchangedmethod("Latticehfieldchanged")
taggroup Latticehklgroup=dlggroupitems(Latticelabel,Latticehfield,Latticekfield,Latticelfield).dlgtablelayout(4,1,0)

taggroup Latticegroup1=dlggroupitems(dgroup,Latticehklgroup).dlgtablelayout(2,1,0).dlginternalpadding(0,5).dlgexternalpadding(0,5)

Latticebox_items.dlgaddelement(Latticegroup1)
Latticetab.dlgaddelement(Latticebox)

CreateCalculatorBoxs.DLGAddElement(crysttablist)
TagGroup Crystalbox=DLGGroupItems(CreateLatticeBoxs,CreateIndexingBoxs,CreateProfileBoxs,CreateCalculatorBoxs).DLGtablelayout(1,4,0)	
CrystalTools.DLGAddElement(Crystalbox)
	
Return CrystalTools
}	
			
		
		// Create the main UI frame		
	TagGroup CreateDialog(Object self)
		{
		TagGroup Tabs = DLGCreateTabList(0,"TabChangeAction").DLGIdentifier("TabList")
		Tabs.DLGAddElement(self.CreateCrystalTools())
		TagGroup DialogItems = DLGCreateDialog("Tabbed")
		DialogItems.DLGAddElement(Tabs)
		
		DialogItems.DLGAddElement(dlgcreatelabel("H. L. Shi (honglongshi@outlook.com), MUC, v1.1, Mar. 2019")).dlgexternalpadding(0,5)
		Return DialogItems
		}
		
	Object Init(Object self)	
		{
		TagGroup ElectronDiffracion_Tools = self.CreateDialog()
		Return self.super.Init(ElectronDiffracion_Tools)
		}
}

CrystUIFrame = Alloc(ElectronDiffracion_Tools).Init()
CrystUIFrame.Display("ElectronDiffracion Tools")

documentwindow dialogwin=getdocumentwindow(0)
number xpos, ypos
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Window position X",xpos)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Window position Y",ypos)

//Positions of GUI
number screenwidth,screenheight
GetScreenSize(screenwidth,screenheight)

if(xpos>=142&&xpos<screenwidth&& ypos>=24&&ypos<screenheight) 
{
windowsetframeposition(dialogwin, xpos, ypos)
}
	
else
{
windowsetframeposition(dialogwin, 142, 24)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Window position X",142)
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Default:Window position Y",24)		
}