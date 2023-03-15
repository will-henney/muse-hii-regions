
	SUBROUTINE multicheck(icid,iwave,ivdif,iajk,iupj,iloj,itr1,itr2,ivc,
     1  num1,num2,out1,out2,out3,outstr)

******	This subroutine checks for the existence of other multiplet
******  members of the incoming line in the master line list, and seeks
******  them within the candidate list.  If other multiplet lines do exist 
******  and a match or matches can be found for them in the candidate list
******  the strongest of all such candidate lines, meeting a certain flux
******  tolerance is considered a match to the expected multiplet line

*	INPUT:
*	icid    - the id of the candidate line currently undergoing
*	          an attempted match - INTEGER
*	iwave   - the wavelength of the potential master line - REAL
*	ivdif   - input master line velocity difference with candidate
*	          line - REAL
*	iajk    - the Einstein coefficient of the potential master
*                 line match - REAL
*	iupj    - the upper "j" value of the potential master line
*	          match - REAL
*	itr1    - the input master line "TR1" specifier - INTEGER
*	itr2    - the input master line "TR2" specifier - INTEGER
*	ivc     - the velocity correction calculated for the 
*	          candidate line undergoing an attempted match - REAL
*
*	OUTPUT:
*	outstr	- the output string stating the results of the check
*	          CHARACTER*80
	
*	INPUT VARIABLES

	INTEGER	icid,itr1,itr2,ncan
	REAL*8	iwave
	REAL ivdif,idwave(2),ivc,iajk,iupj,iloj

*	OUTPUT VARIABLES

	CHARACTER*1 outstr
	INTEGER num1,num2
	REAL*8 out1(50)
	REAL out2(50),out3(50)

*	INTERNAL VARIABLES

	INTEGER ai,si
	PARAMETER (ai=1500,si=1500)
	INTEGER	wtr1(ai),wtr2(ai),wcid(si),wlines,witt
	REAL listbound(2)
	REAL*8 wwave(ai)
	REAL wupeng(ai),wloeng(ai),wajk(ai),wupj(ai),wloj(ai)
        CHARACTER*2 wtts(ai)
	CHARACTER*3 wele(ai)
	CHARACTER*9 wion(ai)
	CHARACTER*11 wuplev(ai),wlolev(ai)
	LOGICAL wusels(ai)

	REAL wavespread,wave(2),ener(2),sum,weight
	REAL clamcor,wvdif,usevel,flxrat,ajkrat,finvel
	PARAMETER (wavespread=400.)
	REAL vtol,ftol(4),atol(2)
	PARAMETER (svtol=1.0)
	DATA atol/0.33,3./
	DATA ftol/1.0E-4,1.0,0.1,10./
	CHARACTER*32 addoutstr
	INTEGER store(100),store2(100),store3(100)
	REAL fstore(100),fstore2(100)

	REAL clam(si),clvel(si),crvel(si),cflux(si)
	INTEGER cnum,cid(si),clab(si)
	REAL t_e,n_e,k_boltz,vlight,hpl
	INTEGER i,j,k,m,o,p,ttype,itt,kk
	REAL*8 expw(100)
	INTEGER ex
	REAL instres,natwidth

*	FORMATS

10	format(i3,'/',i3)
20	format(f9.4,'(',f9.4,')',1x,f7.4,4x)

*	COMMON BLOCKS

	COMMON /canatt/clam,clvel,crvel,cflux,cnum,clab,cid
	COMMON /nebatt/t_e,n_e,k_boltz,vlight,hc
	COMMON /exclusion/expw,ex
	COMMON /instrument/instres,natwidth
	COMMON /listlimits/listbound

*       FUNCTIONS

	ttype(kk)=iand(ishft(kk,-4),7)

**	BEGIN ROUTINE

*	READ IN ALL LINES FROM MASTER LIST MATCHING THE INCOMING
*	LINE'S ATTRIBUTES, WITHIN A WAVELENGTH SPREAD GIVEN BY
*	wave(1) AND wave(2)

	wave(1)=iwave-wavespread
	if (ttype(itr2).gt.0) wave(1)=wave(1)-250 
	if (wave(1).lt.listbound(1)) wave(1)=listbound(1)
	wave(2)=iwave+wavespread
	if (ttype(itr2).gt.0) wave(2)=wave(2)+250
	if (wave(2).gt.listbound(2)) wave(2)=listbound(2)
	ener(1)=9999.99
	ener(2)=9999.99	

	call line_pick(248,wave,ener,itr1,itr2,wwave,wele,
     1  wion,wupeng,wloeng,wuplev,wlolev,wajk,wupj,wloj,wtts,wtr1,wtr2,
     2  wusels,wlines)

c	if (itr1.eq.1747374394) then
c	do iiz=1,wlines
c	   write(*,*) wwave(iiz),wele(iiz),wion(iiz),wupeng(iiz),wloeng(iiz)
c	   read(*,*)
c	end do
c	end if

	num1=wlines-1.
	if (num1.lt.0) num1=0
	sum=0.
	weight=0.
	
*       IF THE MULTIPLET LINES ARE ALL WITHIN EITHER THE SPECTRAL 
*       RESOLUTION OR NATURAL LINE WIDTH OF THE OBJECT, THEN BLEND THEM
*       TOGETHER INTO A STASTICAL WEIGHT WEIGHTED WAVELENGTH, AND EXCLUDE 
*       FROM FUTHER MULTIPLET CHECKS

*          CHECK TO SEE IF ALL THE EXTRACTED LINES ARE WITHIN THE 
*          RESOLUTION/NATURAL WIDTH LIMITS
	   
	if (wlines.le.1) go to 50
        if ((wwave(wlines)-wwave(1))*vlight/wwave(1).gt.instres) 
     1  go to 50
	if ((wwave(wlines)-wwave(1))*vlight/wwave(1).gt.natwidth)
     1  go to 50

*	   CALCULATE WEIGHTED AVERAGE OF MULTIPLET LINES
*	   AND BUILD EXCLUSION LIST

	do ip=1,wlines
           sum=sum+(2*wupj(ip)+1)*wwave(ip)
	   weight=weight+(2*wupj(ip)+1)
	   ex=ex+1
	   expw(ex)=wwave(ip)
	end do

*	   ASSIGN OUTPUT PARAMETERS

	num1=0
        num2=0
	out1(1)=sum/weight
	outstr='*'

*	   BYPASS THE MULTIPLET CHECK

	   go to 400
	   
50	continue

        outstr=' '

*	FOR EACH LINE EXTRACTED FROM THE MASTER LIST WE NOW SEARCH THE
*	CANDIDATE LIST FOR MATCHES

	do k=1,cnum
           wcid(k)=k
	end do

	call sort(-1,clam,wcid,cnum)

	num2=0.
	n=0

	do i=1,wlines

c	   if (itr1.ne.1880803562) go to 300
	   k=0

*	REJECT LINE IF IT IS THE ORIGINAL MULTIPLIT LINE WE ARE
*	TESTING IN THIS CHECK

	   if (itr1.eq.wtr1(i)) go to 300

*       REJECT LINE IF J VALUES DON'T JIVE

c	   if (itr1.eq.1880803562) then
c              print *,'!',wupj(i),'!',wloj(i),'!',iupj,'!',iloj,'!',num1
c	      read(*,*)
c	   end if

	   if (wupj(i).lt.(iupj-1) .or. wloj(i).lt.(iloj-1)) then
c	      if (itr1.eq.1880803562) then 
c		 print *,'Yep'
c		 read(*,*)
d              end if
              num1=num1-1
	      go to 300
           end if
	      
*	CYCLING OVER ALL CANDIDATE LINES

	   do j=1,cnum
        
	      if (wcid(j).eq.icid) go to 100
	      if (clam(wcid(j)).lt.wave(1)) go to 100
	      if (clam(wcid(j)).gt.wave(2)) go to 200

*	   CALCULATE THE VELOCITY SHIFT

	      clamcor=clam(wcid(j))/(1+(ivc/vlight))
	      wvdif=vlight*(clamcor-wwave(i))/wwave(i)                 	      	      	
*	   DETERMINE THE VELOCITY RESIDUAL

	      if (ivdif.le.0) usevel=clvel(wcid(j))
	      if (ivdif.gt.0) usevel=crvel(wcid(j))
	      
*	   REJECTION CRITERIA
              
	      itt=ttype(itr2)
	      witt=ttype(wtr2(i))

	      flxrat=cflux(icid)/cflux(wcid(j))
	      
	      if (itr1.eq.537968665) print *,'yes'
              if (abs(wvdif-ivdif).gt.abs(svtol*usevel)) then	      
c		 if (itr1.eq.537968665) then
c          	 write(28,*) 'Wave Reject: ',iwave,wwave(i),wvdif,
c     1           clam(wcid(j))
c		 end if
	         go to 100
	      end if

*           IF THE TRANSITION COEFFICIENTS ARE AVAILABLE
*           USE THEM IN THE RATIO WITHIN TIGHT CONSTRAINTS,
*           EXCEPT WHEN THE STASTICAL WEIGHTS OF SOURCE
*           LEVEL IS INDETRMINATE (iupj>=49) OR DIFFERING 
*           TRANSITION TYPE: THEN
*           USE STRAIGHT FLUX RATIO REJECTION CRITERIA

	      if ((iajk.ne.0).and.(wajk(i).ne.0) .and. (iupj.lt.49)
     1        .and. (itt.eq.witt)) then
                 ajkrat=iajk*(2*iupj+1)/(wajk(i)*(2*wupj(i)+1))
	         if (flxrat.lt.(atol(1)*ajkrat)) then
c		     if (itr1.eq.537968665) then
c                     write(28,*) 'Ajk Reject 1: ',iwave,wwave(i),ajkrat,
c     1               flxrat,iajk,iupj,wajk(i),wupj(i)
c		     end if
	            go to 100
	         end if
	         if (flxrat.ge.(atol(2)*ajkrat)) then
c                    if (itr1.eq.537968665) then                        
c                    write(28,*) 'Ajk Reject 2: ',iwave,wwave(i),ajkrat,
c     1              flxrat,iajk,iupj,wajk(i),wupj(i)
c		    end if
	            go to 100
	         end if
		 
*          ELSE USE A LOOSER RESTRICITON ON THE SIMPLE RATIO OF
*          THE OBSERVED FLUXES.

              else

*          IF THE PERSPECTIVE MULTIPLET MEMBER IS OF A DIFFERENT 
*          TRANSITION TYPE, INHERENTLY STRONGER TYPE TRANSITION THAN
*          THE PUTATIVE ID (say a "M1" [O III] 4959 vs "E2" 4959 etc...)
*          ASSUME THAT ANY LINE FOUND WITHIN A FACTOR 10^4 HIGHER
*          MIGHT BE THAT LINE

		 if (itt.gt.witt) then
c		    if (flxrat.lt.ftol(1) .or. flkrat.ge.ftol(2)) then
		    if (flxrat.lt.ftol(1)) then
c		       if (itr1.eq.537968665) then
c			  write(28,*) 'ITT Reject: ',iwave,wwave(i),flxrat,
c     1                    clam(wcid(j))
c		       end if
		       go to 100
		    end if
		 else

*          ELSE USE A TIGHTER RATIO

		    if (flxrat.lt.ftol(3) .or. flxrat.ge.ftol(4)) then 
c		    if (flxrat.lt.ftol(3)) then
c		    if (itr1.eq.537968665) then
c     	               write(28,*) 'Flux Reject: ',iwave,wwave(i),flxrat,
c     1		       clam(wcid(j))
c		    end if
		       go to 100
		    end if
                 end if
	      
              end if
		 
*	   CREDIT A "FOUND" LINE

	      k=k+1

	      store(k)=wcid(j)
	      fstore(k)=cflux(wcid(j))

100	      continue
	   
	   end do

200	   continue
	   	   
	   if (k.ge.1) then
	
	      n=n+1

              call sort(1,fstore,store,k)
	
	      fstore2(n)=cflux(store(1))
	      store2(n)=store(1)
	      store3(n)=i
	      
	      num2=num2+1

	   end if

300	   continue

	end do

        if (num2.ge.1) then

           call sort(1,fstore2,store2,num2)
	   call sort(1,fstore2,store3,num2)

	   do n=1,num2

	      if (n.gt.50) then
	         print *,'Saturated'
	         stop
	      end if

	      finvel=clam(store2(n))/(1+(ivc/vlight))
	      finvel=(finvel-wwave(store3(n)))*vlight/
     1        wwave(store3(n))
	      write(addoutstr,20) wwave(store3(n)),clam(store2(n)),
     1        finvel
	      out1(n)=wwave(store3(n))
	      out2(n)=finvel
	      out3(n)=cflux(store2(n))/cflux(icid)
	   end do
	end if

400	continue
	
**	END ROUTINE

	return
	end

	
	      
	   
		
	
	

	   
	   
	   
	

	  
	
 
