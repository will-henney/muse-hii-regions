
	PROGRAM EMILI
******	EMission Line IDing program
******  Ongoing production of BDS
******  Main program began 2/6/01
******  VERSION 3: 6/25/02
******  Atomic Line List v2.04 - P.A.M. van Hoof, Queen's University
******                                            Belfast

* Filler Variables

	character*3 empty3
	character*9 empty9
	character*11 empty11
	character*2 empty2
	character*1 header
	integer dummy
	real blank
	real blankarray(2)
	integer cero	

* Loop Variables

	integer i,j,c,d,k,m,n	

* Storage Arrays
	
	integer kept(1500),reject(1500)

* Misc Variables

	integer er,idex,limit1,limit2,limit3,lastpts,dif
	integer outeleindex,outionindex
	integer stg
        character*4 outele
	real vrc,z
	real maxri,omaxri
	real outprob(2),habun(2)
	real limri
	real kmri(1500),rmri(1500),problimit
	logical matched
	real usevel

* External Procedures

	external alphasort,sort,line_pick,ri,ionize,abundance
	external multicheck,matchlist,vac2air
	external redshift,idi

* Control variables
	
	real wave(2),ener(2),invel(2)
	real incvel(2)
	character*1 outstr,outstr2,outstr21,outstr31,outstr3
	character*1 outstr11(4)
	data (outstr11(i),i=1,4)/'A','B','C','D'/

* Master Line Variables

	real*8 mlam(1500)
	real mjlow(1500),mjhi(1500),majk(1500)
        real mupeng(1500),mloweng(1500),mri(1500),mvrc(1500)
	real kmlam(1500)
	real mvdif(1500)
	real mieng(1500),mprob(1500)
	real*8 mout1(50),moutt1(1500,50)
	real mout2(50),mout3(50)
 	real moutt2(1500,50),moutt3(1500,50)
	character*1 moutstr(1500),mclass(1500)
	character*2 mtt(1500)
	character*3 mele(1500)
	character*9 mion(1500)
	character*11 mup(1500),mlow(1500)
	character*73 mstr
	integer mnum,multi1(1500),multi2(1500)
	integer mistat(1500),mtr1(1500),mtr2(1500)
        logical musels(1500)
	real mab(2)
	integer mm2(100),mm3(100),mm5(100),mm7(100),mm8(1500)
	integer mpts(100)
	real mm1(100),mm4(100),mm6(100)
	integer meleval,mionval
	character*7 mwave

* Candidate Line Variables

	real clam(1500),cdlwave(1500),cdrwave(1500),cflux(1500),cfwhm(1500)
	real clvel(1500),crvel(1500),clamcor(1500),csn(1500),cided(1500)
	integer cion(1500),cele(1500)
	real cincr(500,100),c_idwav(500,100),c_idflx(500,100)
	real c_idvel(500,100),c_idtr1(500,100),c_idtr2(500,100)
	integer cnum,count
	integer cid(1500),clab(1500)
	character*7 cwave

	real listbound(2)

* Physical Variables and Parameters

	real n_e,t_e,hc,vlight,k_boltz,instres,natwidth
	real*8 expw(100)
	integer ex

* Common Blocks
 
	common empty2,empty3,empty9,empty11
        common /canatt/clam,clvel,crvel,cflux,cnum,clab,cid
	common /nebatt/t_e,n_e,k_boltz,vlight,hc	
	common /exclusion/expw,ex
	common /instrument/instres,natwidth
	common /listlimits/listbound

* Initial Values and Constants

        parameter (nsigma=5.)
	parameter (maxsvel=200)
	parameter (sq2pi=2.506628275)	
	parameter (ritol=0.001)
	parameter (hbetaflux=1.)
	parameter (svtol=1.0)
	parameter (nmdisplay=3)
 	parameter (cero=0)	
	parameter (matchtol=0.01)

* Output Formats

50	format('Observed Line: ',f8.2,2x,1pe7.1,0p,2x,'S/N: ',f8.2,
     1  1x,'FWHM: ',f5.1)
51	format(a1,1x,f8.2,' | ',a1,1x,f9.3,a1,a3,1x,a5,i12,1x,i6,2x,
     1  1pe7.1,0p,f6.1,1x,i1,'/',i1,1x,i1,a1,$)
52      format(1x,f8.3,f6.1,$)
53	format(i3,f11.2,f11.2,i3)
55	format(a1,1x,f3.0,1x,a1,1x,f7.2,' | ',a1,1x,f7.2,a1,a3,1x,a4,1pe7.1,
     1  0p,f6.1,2x,i1,'/',i1,$)
56      format(f8.2,2x,1pe7.1,0p,1x,' | ',1x,$)
57      format(a,1x,a,1x,f8.2,', ',$)

****  BEGIN MAIN PROGRAM

*  setting constant parameters transfered in common blocks

	hc=1.24E-4
	vlight=3.0E5
	k_boltz=1.38E-23
	
*  greeting

        write(*,*)
        write(*,*) '             :) Welcome to EMILI :) '
        write(*,*) '            Ongoing Production of BDS '
        write(*,*) '       2000-2002 Michigan State University'
        write(*,*) '          Using: Atomic Line List v2.04 '
        write(*,*) 'Property of: P. van Hoof - Queens University-'//
     1  'Belfast'
	write(*,*)                           

*  open all files

	call openall
	
*  set up initial ionization table

	call ionize(1,dummy,dummy,blank,dummy)

*  set up the master list records lines/levels

	call line_pick(0,blankarray,blankarray,dummy,dummy,mlam,
     1  mele,mion,mupeng,mloweng,mup,mlow,majk,mjhi,mjlow,mtt,
     2  mtr1,mtr2,musels,mnum)

*  set up initial abundence and velocity correction grids
*  read in observing parameters
	
        call matchlist

*  calculate the hbeta normalization factor

	call abundance(ishft(1,26),habun)
	call ri(habun,'H','I',1,'E1',8.41E+6,12.75,99.,999,t_e,n_e,hbetarinorm)        
	print *,'Hbeta FLux: ',hbetarinorm


*  readying to read in all candidate lines

	er=1
	d=0

	write(*,*) 'SET-UP: Complete'
	write(*,*) 'Press RETURN to begin run'
	read(*,*)

1000	c=c+1
	read(25,'(a1)',end=2000,err=9900) header
	if (header.eq.'#') go to 1000
	d=d+1
	backspace(25)
	
	read(25,*,err=9900) clam(d),cdlwave(d),cdrwave(d),
     1  cflux(d),cfwhm(d),csn(d)
	count=d
	cid(d)=count
	clab(d)=c

*	Normalize with respect to H beta Flux

	cflux(d)=cflux(d)/hbetaflux
	clvel(d)=cdlwave(d)*vlight/clam(d)
	crvel(d)=cdrwave(d)*vlight/clam(d)
	
	go to 1000
	
2000	continue

	cnum=d

*   Sort the lines in order of wavelength

	call sort(-1,clam,cid,cnum)
c	print *,clam(cid(1))
	
*   Set the list boundries

	listbound(1)=clam(cid(1))
	listbound(2)=clam(cid(cnum))

*   Read in all lines from master list with wavelengths
*   between the limits wave(1)-wave(2)
		
	do j=1,cnum
	   
	   wave(1)=clam(cid(j))*(1-maxsvel/vlight)
	   wave(2)=clam(cid(j))*(1+maxsvel/vlight)
	   ener(1)=99.
	   ener(2)=99.
	   write(*,50) clam(cid(j)),cflux(cid(j)),csn(cid(j)),
     1     cfwhm(cid(j))
	   write(28,50) clam(cid(j)),cflux(cid(j)),csn(cid(j)),
     1     cfwhm(cid(j))

*   Ready the summary list

	   write(27,56) clam(cid(j)),cflux(cid(j))

*   Initialize the exclusion array

	   do ip=1,100
	      expw(ip)=0.
	   end do

	   ex=0.

*   Read in the lines, store elements in "m" series array	
	
	   call line_pick(128,wave,ener,dummy,dummy,mlam,
     1     mele,mion,mupeng,mloweng,mup,mlow,majk,mjhi,mjlow,mtt,
     2     mtr1,mtr2,musels,mnum)
	   
	   m=0
	   maxri=9.0E-29
	   omaxri=9.0E-29
	   
	   do k=1,mnum	      

*      Read in the ionization energy of a particular
*      master line

	      call ionize(2,mtr1(k),mtr2(k),mieng(k),mistat(k))

*      Read in the abundence info for the master line

              call abundance(mtr1(k),mab)

*      Read in the velocity correction to the velocity
*      diference between master and candidate line

	      call redshift(mtr1(k),mtr2(k),mloweng(k),vrc)

              mvrc(k)=vrc

*      Obtain the residual velocity difference
*      after correction

	      clamcor(k)=clam(cid(j))/(1+(vrc/vlight))
	      mvdif(k)=vlight*(clamcor(k)-mlam(k))/mlam(k)

	      if (mvdif(k).le.0) then
                 mprob(k)=mvdif(k)/clvel(cid(j))
	      else
	         mprob(k)=mvdif(k)/crvel(cid(j))
	      end if

*	Calculate relative intensity and normalize to H Beta

	      if (mprob(k).le.(2*nsigma)) then

	         call ri(mab,mele(k),mion(k),mistat(k),mtt(k),majk(k),
     1	         mupeng(k),mloweng(k),mtr1(k),t_e,n_e,mri(k))

	         mri(k)=mri(k)/hbetarinorm

*	      Accept master line as potential match if within
*             wavelength tolerance
	         if (mprob(k).le.nsigma) then 
	            m=m+1
	            kept(m)=k
	            if (mri(k).gt.maxri) maxri=mri(k)
	         else
*             Find strongest line outside of search range
c		 if (mtr1(k).eq.673189889) print *,k,mri(k)
		    if (mri(k).gt.omaxri) then 
		      omaxri=mri(k)
		      stg=k
		    end if
	         end if
              end if
    	   end do

c	   print *,stg,mnum,k
c	   read(*,*)

*	Retain only those master lines with r.i. within "witol" of the 	 
*	maximum calculated for any master line within the probability
*	range

	   limri=ritol*maxri
	    
	   limit1=m
	   m=0
	   n=0


	   do k=1,limit1

	      if ((mri(kept(k))/maxri).gt.ritol) then
	         m=m+1
	         idex=kept(k)
	         kept(m)=idex
		 kmri(m)=mri(idex)
              else
	         n=n+1
	         idex=kept(k)
	         reject(n)=idex
	         rmri(m)=mri(idex)
	      end if
	   end do   
    	
           limit2=m

*	Sort the prespective master list match lines by wavelength

	   do ip=1,limit2
	      kmlam(ip)=mvdif(kept(ip))
	   end do

	   call sort(1,kmlam,kept,limit2)

*	   Call Multiplet Check

	   limit3=m+1
	   kept(limit3)=stg
	   
	   m=0

   	   do k=1,limit3

*            Check to see that the line is a pure LS coupling line

	     if (.not. musels(kept(k))) then
                moutstr(kept(k))='$'
		multi1(kept(k))=0
		multi2(kept(k))=0
		go to 2990
             end if

*	     Check the exclusion array:

	     do ip=1,ex
	        if (expw(ip).eq.mlam(kept(k))) go to 3000
	     end do

*	     Call the actual multiplet check:

	     call multicheck(cid(j),mlam(kept(k)),mvdif(kept(k)),
     1       majk(kept(k)),mjhi(kept(k)),mjlow(kept(k)),
     2	     mtr1(kept(k)),mtr2(kept(k)),
     3       mvrc(kept(k)),multi1(kept(k)),multi2(kept(k)),mout1,
     4       mout2,mout3,outstr)

	     moutstr(kept(k))=outstr

*	     assign values of multiplet check matches:

	     do ip=1,multi2(kept(k))
	        moutt1(kept(k),ip)=mout1(ip)
	        moutt2(kept(k),ip)=mout2(ip)
	        moutt3(kept(k),ip)=mout3(ip)
	     end do

*	     if blended line, calculate new values:

	     if (moutstr(kept(k)).eq.'*') then
                mlam(kept(k))=mout1(1)
                mvdif(kept(k))=(clamcor(kept(k))-mlam(kept(k)))*vlight/
     1          mlam(kept(k))
             end if

*	     assign final output set:

2990	      continue

	      if (kept(k).ne.stg) then
                 m=m+1
	         idex=kept(k)
	         kept(m)=idex
	      end if

3000	      continue

	   end do

*	Sort the master lines by velocity difference

	   limit2=m

	   do ip=1,limit2
	      kmlam(ip)=mvdif(kept(ip))
	   end do

	   call sort(1,kmlam,kept,limit2)

*	Call the IDI routine

	    invel(1)=clvel(cid(j))
	    invel(2)=crvel(cid(j))

	    
	    do ip=1,limit2
	       if (kept(ip).ne.stg) then
	          mm1(ip)=mri(kept(ip))
	          mm2(ip)=multi1(kept(ip))
	          mm3(ip)=multi2(kept(ip))
	          mm4(ip)=mvdif(kept(ip))
	          mm5(ip)=kept(ip)
	       end if
	    end do

	    call idi(mm1,mm2,mm3,mm4,invel,limit2,mm5,mm6,mm7)

	    dif=0
	    lastpts=99.

	    do ip=1,limit2
               mpts(mm7(ip))=nint(mm6(ip))
               if (mpts(mm7(ip)).ne.lastpts) then 
                  lastpts=mpts(mm7(ip))
	          dif=ip
	       end if
	       if (dif.le.4) then              
	          mclass(mm7(ip))=outstr11(dif)   
                  mm8(mm7(ip))=dif    
	       else
                  mclass(mm7(ip))=' '
	          mm8(mm7(ip))=0
	       end if
	    end do

	mclass(stg)='<'
	mpts(stg)=0.

	limit3=m+1
	kept(limit3)=stg
	   
*	Temporary Output

	do k=1,limit3

	   outstr3=' '
	   if (kept(k).eq.stg) outstr3='>'

	   if (mvdif(kept(k)).gt.0) then
	      usevel=crvel(cid(j))
	      if (mvdif(kept(k)).lt.(svtol*usevel)) then
                 outstr2='+'
	      else
                 outstr2=' '
	      end if
	   else
	      usevel=clvel(cid(j))
              if (mvdif(kept(k)).gt.(svtol*usevel)) then
                 outstr2='+'
	      else
                 outstr2=' '
	     end if
           end if
   
*	      write the output:
	 
	      write (*,51) outstr2,
     3        clamcor(kept(k)),outstr3,mlam(kept(k)),
     1        moutstr(kept(k)),mele(kept(k)),mion(kept(k)),
     2        mtr1(kept(k)),mtr2(kept(k)),mri(kept(k)),
     2        mvdif(kept(k)),multi1(kept(k)),multi2(kept(k)),
     4        mpts(kept(k)),mclass(kept(k))
	      write (28,51) outstr2,
     3        clamcor(kept(k)),outstr3,mlam(kept(k)),
     1        moutstr(kept(k)),mele(kept(k)),mion(kept(k)),
     2        mtr1(kept(k)),mtr2(kept(k)),mri(kept(k)),
     2        mvdif(kept(k)),multi1(kept(k)),multi2(kept(k)),
     4	      mpts(kept(k)),mclass(kept(k))
	      
	      do ip=1,nmdisplay
	         if (ip.gt.multi2(kept(k))) go to 4000
                 write(*,52) moutt1(kept(k),ip),moutt2(kept(k),ip)
                 write(28,52) moutt1(kept(k),ip),moutt2(kept(k),ip)
	      end do

*	      write the quick list

4000	      continue

	      outele='    '
	      if (mclass(kept(k)).eq.'A') then
	      if (mele(kept(k))(1:1).eq.' ') then
                 outele(1:2)=mele(kept(k))(2:3)
              else
                 outele=mele(kept(k))
              end if
	      outeleindex=index(outele(2:4),' ')
	      outionindex=index(mion(kept(k)),' ')-1
              write(27,57) outele(1:outeleindex),mion(kept(k))(1:
     1        outionindex),mlam(kept(k))
	      end if

5000	      continue

	      write(*,*)
	      write(28,*)

	   end do

6000	   continue

	   write(27,*)
	   write(28,*)   
	   write(*,*)


	end do
            
*       TEMPORARY ENDING

9000	write(*,*) 'The end'
	close(21)
	close(22)
	close(23)
	close(25)
	close(27)
	close(28)
	stop

9900	write(*,*) 'WARNING:EMILI.F main: '
	write(*,*) 'Error in reading input line list at line: ',c
	stop
   
	end

	
