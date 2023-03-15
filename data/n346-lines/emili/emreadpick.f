	SUBROUTINE line_pick(stat,wave,ener,itr1,itr2,owave,oele,
     1  oion,oupeng,oloeng,ouplev,ololev,oajk,oupj,oloj,otts,otr1,otr2,
     2  ousels,nlines,oz1,oz2)

*****	This subroutine reads in all lines/levels in the master line list
*****	and stores their contents in various arrays when stat=0
*****   When stat is non-zero, its value determines what lines from the
*****   master list should be sent to the main program or subroutine.
*****	This was written as of 7/03/01 to act as a bridge between
*****   the full line list and the remainder of the qdemili routine

*	INPUT
*	stat	   setup/or search(with parameters) - INTEGER
*	wave(2)    wavelength range of scan - array of REAL
*	ener(2)	   energy range of upper level of scan - array(REAL)
*	itr1       scan line tr1 identifier - INTEGER
*	itr2       scan line tr2 identifier - INTEGER
*
*	OUTPUT
*	owave(*)   wavelengths of selected lines - array of REAL
*	oele(*)    element description - array of CHAR(3)
*	oion(*)    ionization stage description - array of CHAR(9)
*	oupeng(*)  upper energies of selected lines - array of REAL
*	oloeng(*)  lower energies of selected lines - array of REAL
*       ouplev(*)  upper level description - array of CHAR(11)
*	ololev(*)  upper level description - array of CHAR(11)
*	oajk(*)    transition probabilities - array of REAL
*       oupj(*)    upper level "j" value - array of REAL
*       oloj(*)    lower level "j" value - array of REAL
*	otts(*)    transition types - array of CHAR(2)
*	otr1(*)	   output tr1 descriptor - array of INTEGER
*	otr2(*)	   output tr2 descriptor - array of INTEGER
*	nlines	   number of lines extracted - INTEGER

*	INPUT VARIABLES
	
        REAL wave(2),ener(2)
	INTEGER itr1,itr2,stat

*	OUTPUT VARIABLES
	
	REAL*8 owave(1500)
	REAL oupeng(1500),oloeng(1500),oajk(1500),oupj(1500),oloj(1500)
	INTEGER otr1(1500),otr2(1500),nlines
	CHARACTER*2 otts(1500)
	CHARACTER*3 oele(1500)
	CHARACTER*9 oion(1500)
	CHARACTER*11 ololev(1500),ouplev(1500)
	CHARACTER*33 oz1,oz2
	LOGICAL ousels(1500)

*	INTERNAL VARIABLES

	INTEGER k,m,nl,n2,nsig,wpar,wstore(10000)
	REAL lowopt,upopt,wdif,wzero
	LOGICAL optical			
	PARAMETER (optical=.true.,lowopt=3000.,upopt=11000,wdif=1,
     1  wzero=3000)

	DOUBLE PRECISION rfr_ind,wavinv	
	LOGICAL vacumn,fine_structure,accuracy
	REAL*8 acceptwav
	CHARACTER*2 w2upusels,w2lousels
	PARAMETER (vacumn=.false.,fine_structure=.true.,accuracy=.false.)
	PARAMETER (acceptwav=100)
	INTEGER mline,mlev
	PARAMETER (mlin=300000,mlev=45373)
	INTEGER i,ii,jj,cni_grnd(10),ionnr,cni2(10,2560),ind(2560)
	REAL lev1,lev2,levlim(2560),ionlim(2560)
	CHARACTER*11 tmp5
	CHARACTER*33 tmp6,ground(2560)
	CHARACTER*40 n1,fm

	INTEGER k1,s2(mlev),l(mlev),j2(mlev),par(mlev),cni(10,mlev)
	INTEGER nlev,levwarn(mlev),lev_mrk(mlev),nconfig(mlev)
	DOUBLE PRECISION level(mlev),lev_err(mlev)
	CHARACTER*2 tt(mlev)
	CHARACTER*3 ref(mlev)
	CHARACTER*33 config(mlev)
	CHARACTER*11 term(mlev)
	LOGICAL ground_term(mlev)
	
	REAL*8 wav,dwav
	REAL a_ki,wupeng,wloeng,wjup,wjlo
	INTEGER tr1,tr2,ref1,ref2,ref3,wele,wion,wtrt,tlines,wup,wlo
	CHARACTER*33 awupconfig,awloconfig
	INTEGER wupconfig,wloconfig

	LOGICAL nebular
	PARAMETER (nebular=.false.)
	CHARACTER*11 wuplev,wlolev
	REAL*8 swav(mlin),sdwav(mlin)
	REAL sa_ki(mlin),svwave(mlin)
	REAL svwaveele
	INTEGER str1(mlin),str2(mlin)
	
	INTEGER a1,startval,endval
	LOGICAL search(8)
	
	INTEGER w2ele,w2ion,w2trt,oionindex
	REAL w2upeng,w2loeng,w2jup,w2jlo
	CHARACTER*11 w2uplev,w2lolev
	INTEGER w2upconfig,w2loconfig

        CHARACTER*2 trantype(4)
	DATA trantype/'E1','E1','M1','E2'/ 
	CHARACTER*2 slab
	CHARACTER*6 spcs
	INTEGER nrsp
	PARAMETER (nrsp=30)
	INTEGER eonnum(mlev,20),stonnum(mlev,20)
	INTEGER eoninc(mlev),stoninc(mlev)
 
	REAL n_e,t_e,hc,vlight,k_boltz,instres,natwidth

*       COMMON BLOCKS

	COMMON /departs/eonlev,stonlev,eoninc,stoninc
c	COMMON /nebatt/t_e,n_e,k_boltz,vlight,hc	
	COMMON/SPC_NOT/ slab(nrsp),spcs(nrsp)

	DATA slab/'H ','He','Li','Be','B ','C ','N ','O ',
     1  'F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar',
     2  'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni',
     3  'Cu','Zn'/
	DATA spcs/'I','II','III','IV','V','VI','VII','VIII',
     1  'IX','X','XI','XII','XIII','XIV','XV','XVI','XVII',
     2  'XVIII','XIX','XX','XXI','XXII','XXIII','XXIV','XXV',
     3  'XXVI','XXVII','XXVIII','XXIX','XXX'/       

*	SAVED VARIABLES

	SAVE swav,sdwav,sa_ki,str1,str2,ground_term,level,tlines
	SAVE svwave,ind
     
*	FORMATS

10	format(a10,i5,1x,a33,9i3,i6,2f23.12)
20	format(i4,i4,1x,2a,l2,2i3,i6,i3,1x,a,9i3,i6,2f23.12,2i4,1x,a3)
30	format(f19.15,1pe8.1,e13.4,0p,i12,i7)
40	format(2h(f,i2.2,1h.,i2.2,1h,,i2.2,2hx))
50	format(2h(f,i2.2,1h.,i2.2,1h))

*	FUNCTIONS

	vnr(k)=iand(k,31)
	vl(k)=iand(ishft(k,-5),15)
	vn(k)=iand(ishft(k,-9),127)
	flu(kk)=iand(kk,1)
	flr(kk)=iand(kk,2)
	flc(kk)=iand(kk,4)
	ilo(k)=iand(k,1023)
	ihi(k)=iand(ishft(k,-10),1023)
	ion(k)=iand(k,63)
	elm(k)=iand(ishft(k,-6),63)
	spc(k)=iand(ishft(k,-20),2047)
	nd(kk)=iand(kk,15)
	ttype(kk)=iand(ishft(kk,-4),7)
	tpc(kk)=iand(ishft(kk,-7),15)+1.

        hc=1.24E-4
        vlight=3.0E5 
        k_boltz=1.38E-23


**	IF stat=0 THIS IS A SET-UP CALL

	if (stat.eq.0) then
 
*	SETUP INITIAL PARAMETERS

      	k=1.
	nl=1.
	n2=0.
	wpar=0.
	
*	READ IN THE VARIOUS LEVELS

*	READ IN ELEMENT/ION HEADER

	write(*,'(1x,a,$)') 'READING IN LEVELS: '

100	read(23,10,end=200,err=990) tmp5,i,tmp6,(cni_grnd(m),m=1,10),
     1  lev1,lev2

	nl=nl+1
	jj=ionnr(tmp5)
        if (jj.gt.0) then
           ind(jj)=k-1
	else
           write(*,'('' unknown spectrum: '',A)') tmp5
           go to 991
        end if
	ground(jj)=tmp6
        do m=1,10
           cni2(m,jj) = cni_grnd(m)
        end do
        levlim(jj)=lev1
        ionlim(jj)=lev2
 
*	READIN IN THE LEVELS FOR THAT ELEMENT/ION

	do jj=1,i
           	
	   read(23,20,err=991) nconfig(k),k1,config(k),term(k),
     1     ground_term(k),s2(k),l(k),j2(k),par(k),tt(k),
     2     (cni(k2,k),k2=1,10),level(k),lev_err(k),levwarn(k),
     3     lev_mrk(k),ref(k)

	   if (level(k) .eq. -1.0D0) level(k)=1.d100
           k=k+1
           nl=nl+1
        end do

	go to 100

200	continue

	nlev=k-1
	write(*,*) nlev

*	READING IN LINES

        write(*,'(1x,a,$)') 'READING IN LINES: '

300	read(24,30,end=400,err=992) wav,dwav,a_ki,tr1,tr2

*	PUT WAVELENGTHS IN ANGSTROMS AND AIR if (vacumn=.false.)

	svwaveele=wav
	nsig=nd(tr2)

	if (.not. vacumn) call vac2air(wav)

	write(fm,40) nsig+4,nsig,12-nsig
	write(n1,fm) wav
	write(fm,50) nsig+4,nsig
	read(n1,fm) wav

	wav=wav*1.0E4
	dwav=dwav*1.0E4

*	ENTERING FILTERS

	   ref1=spc(tr1)
	   ref2=ihi(tr1)+ind(ref1)
	   ref3=ilo(tr1)+ind(ref1)
           
           wele=elm(ref1)
	   wion=ion(ref1)
	   wtrt=ttype(tr2)
	   wupeng=level(ref2)
	   awupconfig=config(ref2)
	   awloconfig=config(ref3)

*       "Optical" FILTER (3000-11000 A) only included during reading
*	if optical=.true.

	if (optical) then
    	   if (wav.lt.lowopt) go to 300
	   if (wav.gt.upopt) go to 400
	end if

*	"Nebular" FILTER, if Hydrogen and Helium only permitted
*	"E1" transitions, if other element only "E2" or "M1" transitions
*	with upper energy level with energy less than 25000 cm^-1 allowed

	if (nebular) then
	   
	   if ((wele.lt.3).and.(wtrt.ne.0)) go to 300
	   if ((wele.gt.2).and.(wtrt.lt.2 .or. wtrt.gt.3)) go to 300
	
	end if

***	"Fine-Structure" FILTER, eliminate all fine structure lines
*	from Hydrogen-like ions (i.e. those with ground electron
*	configuration of 1s^1 => H I, He II, Li III, Be IV, B V, etc...)

	if (fine_structure) then

	   if (wele.eq.(wion+1)) then

	      wup=index(awupconfig,'*')
	      wlo=index(awloconfig,'*')
	   
	      if ((wup.eq.0).and.(wlo.eq.0)) go to 300
	
	   end if
	
	end if

***     "Accuracy" FILTER, eliminate all transitions with wavlength
*       uncertainties greater than "acceptwav" (in km/sec)

	if (accuracy) then
	
	   if (dwav*vlight/wav .gt. acceptwav) go to 300

	end if

*	STORING THE WORKING ARRAYS

	n2=n2+1

	svwave(n2)=svwaveele
	swav(n2)=wav
	sdwav(n2)=dwav
	sa_ki(n2)=a_ki
	str1(n2)=tr1
	str2(n2)=tr2

*	READING IN READY WAVELENGTH REFERENCE

	do while (wav.ge.((wpar*wdif)+wzero))
	   wpar=wpar+1
	   if (wav.eq.(((wpar-1)*wdif)+wzero)) wstore(wpar)=n2
	   if (wav.gt.(((wpar-1)*wdif)+wzero)) wstore(wpar)=n2-1
	end do
	
*	RETURN TO READ THE NEXT LINE

	go to 300

400	continue

	nlines=n2
	tlines=n2
	write(*,*) n2

*	CLOSING THE INPUT FILES

	close(23)
	close(24)

	end if

**	IF stat not equal 0 THIS IS A SEARCH

	if (stat .ne. 0) then

*	SET-UP INITIAL PARAMETERS

	fin=0.
	nlines=0.

*	All SEARCHES CARRIED NOT CARRIED OUT BY DEFAULT

	do i=0,7
	   search(i+1)=.false.
	   a1=ishft(128,-i)
	   if (iand(stat,a1).eq.a1) search(i+1)=.true.
c	   write(*,*) search(i+1)
	end do


*	DECODE THE INCOMING LINE INDENTIFIER IF ANYTHING
*	BUT A STRICT WAVELENGTH SEARCH

 	if (stat.ne.128) then      
	   
           ref1=spc(itr1)
	   ref2=ihi(itr1)+ind(ref1)
	   ref3=ilo(itr1)+ind(ref1)
           
           wele=elm(ref1)
	   wion=ion(ref1)
	   wtrt=ttype(itr2)
	   wupeng=level(ref2)
	   wloeng=level(ref3)
	   wuplev=term(ref2)
	   wlolev=term(ref3)
	   wjup=(j2(ref2)-1)/2.
	   wjlo=(j2(ref3)-1)/2.
	   wupconfig=nconfig(ref2)
	   wloconfig=nconfig(ref3)

	end if	

*	Wavelength Search

	startval=1.
	endval=tlines

	if (.not. search(1)) go to 410
	
	startval=wstore(int(((wave(1)-wzero)/wdif)+0.5))

410     continue
	
*	READ IN A CANDIDATE LINE

	iii=1

	do ii=startval,endval

           if (search(1)) then
	      if (swav(ii).lt.wave(1)) go to 490
	      if (swav(ii).gt.wave(2)) go to 500
	   end if
	   
*	DETERMINE THE ATRIBUTES OF THE MASTER LIST LINE

           ref1=spc(str1(ii))
	   ref2=ihi(str1(ii))+ind(ref1)
	   ref3=ilo(str1(ii))+ind(ref1)
	   
           w2ele=elm(ref1)
	   w2ion=ion(ref1)
	   w2trt=ttype(str2(ii))
	   w2upeng=level(ref2)
	   w2loeng=level(ref3)
	   w2uplev=term(ref2)
	   w2lolev=term(ref3)	 
	   w2jup=(j2(ref2)-1)/2.
	   w2jlo=(j2(ref3)-1)/2.
	   w2upconfig=nconfig(ref2)
	   w2loconfig=nconfig(ref3)
           w2upusels=tt(ref2)
           w2lousels=tt(ref3)
	   
*	SKIP THOSE LINES WITH AN UNKNOWN TRANSITION TYPE

	   if (w2trt.gt.3) go to 490

*	Element/Ion Search

	   if (.not. search(2)) go to 420
	   if (w2ele.ne.wele) go to 490
           if (w2ion.ne.wion) go to 490
	   
420	   continue
	
*	Transition Type Search

	   if (.not. search(3)) go to 430
           if (w2trt.gt.wtrt) go to 490

430	   continue

*	Upper Level Search

	   if (.not. search(4)) go to 440
           if (w2uplev.ne.wuplev) go to 490
	   if (w2upconfig.ne.wupconfig) go to 490
	     
440	   continue

*	Lower Level Search

           if (.not. search(5)) go to 450
           if (w2lolev.ne.wlolev) go to 490
	   if (w2loconfig.ne.wloconfig) go to 490
	   	
450	   continue

*       Level Energy Search

	   if (.not. search(6)) go to 460
           if (w2upeng.gt.ener(1)) go to 490
	   if (w2upeng.lt.ener(2)) go to 490
	   	      
460	   continue
	   
*	Upper J Search

	   if (.not. search(7)) go to 470
	   if (w2jup.lt.(wjup-1)) go to 490

470	   continue

*	Lower J Search

	   if (.not. search(8)) go to 480
	   if (w2jlo.lt.(wjlo-1)) go to 490

480	   continue

*	WRITE THE SPECTRUM

           fin=fin+1

	   if (fin.eq.1499) then
	      print *,'Array Saturated'
	      go to 500
	   end if

	   owave(fin)=swav(ii)
	   oele(fin)=slab(w2ele)
	   if (w2trt.gt.1) oele(fin)='['//oele(fin)(1:2)
	   if (w2trt.le.1) oele(fin)=' '//oele(fin)(1:2)

	   oion(fin)=spcs(w2ion+1)
	   if (w2trt.gt.0) then
	      oionindex=index(oion(fin),' ') 
	      oion(fin)(oionindex:oionindex)=']'
           end if

           oupeng(fin)=w2upeng*hc
	   oloeng(fin)=w2loeng*hc
	   ololev(fin)=w2lolev
           ouplev(fin)=w2uplev	   	   
	   oajk(fin)=sa_ki(ii)
	   oupj(fin)=w2jup
	   oloj(fin)=w2jlo
	   otts(fin)=trantype(w2trt+1)
	   otr1(fin)=str1(ii)
	   otr2(fin)=str2(ii)
	   if (itr1.eq.str1(ii) .and. itr2.eq.str2(ii)) then
              oz1=config(ref2)
	      oz2=config(ref3)
	   end if
	   
*       Checking to see if both source and destination levels arrise
*       from LS coupling (only LS coupling level lines are checked
*       by multiplet check)

           ousels(fin)=.false.
           if (w2lousels.eq.'LS' .and. w2upusels.eq.'LS')
     1     ousels(fin)=.true.

c	   if (itr1.eq.1745497179) then
c	   write(*,*) '!'
c	   write(*,*) owave(fin),oele(fin),oion(fin),oupeng(fin),
c     1     oloeng(fin),ololev(fin),ouplev(fin),hc
c	   read(*,*)
c	   end if

490	   continue

*	FINISHED WITH LINE

        end do
	
500	continue

	end if

	nlines=fin
	go to 999

*	ERROR MESSAGES

990	write(*,*) 'WARNING: Line_Pick.f err=990'
	write(*,*) 'Error reading the header of the level table'
	stop

991	write(*,*) 'WARNING: Line_Pick.f err=991'
	write(*,*) 'Error reading a level'
	stop

992	write(*,*) 'WARNING: Line_Pick.f err=992'
	write(*,*) 'Error reading a line'
	stop

*	SUBROUNTINE COMPLETED

999	continue
	return	   
	end	

	INTEGER FUNCTION ionnr(str)

******  This function returns the ionization stage of the level
	
	CHARACTER*(*) str
	INTEGER nrsp,i,i1,i2
	CHARACTER*2 slab
	CHARACTER*6 spcs
	PARAMETER (nrsp=30)
	COMMON/SPC_NOT/ slab(nrsp),spcs(nrsp)

	i1=0
	i2=0
	do i=1,nrsp
	   if (str(1:2).eq.slab(i)) i1=i
	end do
	do i=1,nrsp
	   if (str(4:9).eq.spcs(i)) i2=i
	end do
	if (i1.eq.0 .or. i2.eq.0) then
	   ionnr=0
	else
	   ionnr=64*i1+i2-1
	end if
	return
	end
	
		
	

	
	   

		
	
