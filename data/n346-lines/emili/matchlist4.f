
	SUBROUTINE matchlist

******	This routine determines a rough ionization distribution of the 
******  observed nebulae using matched lines provided by the user (read
******  in by routine "matchlist.f") with known labratory wavelengths.
******  Specific lines are searched for in this input list, which establish
******  the ICF for each element, based upon abundances read in from the
******  data file "abun.dat" which are with respect to hydrogen.  This
******  routine builds an array of ICF versus ion (in "iabun" array)
******  which can be accessed by the main program
*
*	CODING COMPLETED: 8/17

	IMPLICIT NONE

*	VARIABLES:

        INTEGER i,ii,iy,j,k,m,n,o,npar,nrsp,numbins,numratio
	INTEGER nion,ind
        CHARACTER*1 header
        PARAMETER (npar=50,nrsp=30,nion=10,numbins=5,numratio=2)	
	LOGICAL obsparameters_found,list_done,use_def_rvcor
	REAL t_e,n_e,width,instres
	CHARACTER*3 matchele(npar)
	CHARACTER*6 matchion(npar)
	CHARACTER*10 response1,response2
	REAL k_boltz,vlight,hc
	INTEGER prevbin,binnum(numbins)
	REAL binsum(numbins)
	DATA (binnum(i),i=1,numbins)/numbins*0./
	DATA (binsum(i),i=1,numbins)/numbins*0./
        REAL wavedif,vtol,bins(numbins)
	PARAMETER (vtol=1.)
	DATA (bins(i),i=1,numbins)/0.0,13.6,24.7,54.5,100.0/
	REAL*8 swave(npar)
	REAL obswave(npar),labwave(npar),obsflx(npar)
	REAL ioneng(nrsp,nion),ivcor(nrsp,nrsp+1)
	REAL rvcor
	INTEGER mindif,iions(nrsp)
	REAL def_rvcor(numbins),def_icf(numbins)
	DATA (def_rvcor(i),i=1,numbins)/0.,0.,0.,0.,0./
	DATA (def_icf(i),i=1,numbins)/0.05,0.35,0.35,0.2,0.1/
	INTEGER str1,str2,sele,sion,sttype,srecord(npar)
	INTEGER elm,ion,ttype,bin,strength,nummatched,telenum,tionnum
	INTEGER sid(npar,npar),sided,slist(npar,npar,npar)
	INTEGER trecord(npar),tlist(npar),tbin,tstrength
	REAL tflx(npar),ix(numbins),ixr(numratio),iabun(nrsp,nrsp+1)
	LOGICAL vacumn,line31,line32,line4,use_def_icf,use_match
	LOGICAL manicf,manvel
        DATA vacumn/.false./
	REAL min_val(numbins),thresh(numbins,numbins),val(numbins,numbins)
	REAL max_val(numbins)
	DATA (min_val(i),i=1,numbins)/5*1.0E-3/
	DATA (thresh(1,i),i=1,numbins)/1.0E-4,1.0E-2,3*1.0/
	DATA (val(1,i),i=1,numbins)/1.0E-3,1.0E-2,1.0E-1,2*1.0/
	DATA (max_val(i),i=1,numbins)/4*0.5,0.3/
  	REAL flx31,flx32,flx4,remain,abun,y,y_def,habun,uabun(nrsp)
	REAL vel(numbins),icf(numbins)
        PARAMETER (y_def=0.10)
	CHARACTER*2 tele,element
	CHARACTER*2 slab(nrsp)
	CHARACTER*6 spcs(nrsp),tion
	REAL depeleamount(30)
	INTEGER depele(30),depelehi(30),depelelow(30),depelenumber
	
*	COMMON BLOCKS:

	COMMON /nebatt/t_e,n_e,k_boltz,vlight,hc
	COMMON /spc_not/slab,spcs
	COMMON /ionization/ioneng,iions
	COMMON /instrument/width,instres	
	COMMON /abun/iabun
	COMMON /velcor/ivcor
        COMMON /calculates/use_def_rvcor,use_def_icf,use_match
	COMMON /calculates2/manicf,manvel
        COMMON /specified/vel,icf
	COMMON /depleted/depeleamount,depele,depelehi,depelelow,
     1  depelenumber

*	FUNCTIONS:

	elm(k)=iand(ishft(k,-26),63)
	ion(k)=iand(ishft(k,-20),63)
	ttype(k)=iand(ishft(k,-4),7)
	bin(k)=iand(ishft(k,-3),7)
	strength(k)=iand(k,7)
	
*	FORMATS:

54	format(f13.11,i10,i4,i4)
55	format(a3,25x,a6,'- ',a6,2x,f9.2)

**	BEGINING SUB-ROUTINE

*       IF A MATCHED LIST IS NOT USED, SKIP AHEAD TO ABUNDANCE FILE READ

        if (.not. use_match) go to 195 

*       OTHERWISE PLUNGE AHEAD

	write (*,*) 'READING IN: User Match List'

*	SETTING INITIAL LOOP PARAMETERS

	i=0
	j=0

*	READ IN THE MATCH LINE ATTRIBUTES

100     j=j+1     
	read(26,'(a1)',end=110,err=900) header
	if (header.eq.'Z') go to 100
	i=i+1
	backspace(26)
	read(26,*,err=900) obswave(i),labwave(i),matchele(i),matchion(i),
     1  obsflx(i)
 	
	go to 100
	
110	continue
	close(26)
	   
	nummatched=i

*	READ IN THE SPECIAL ASSIGNMENT LIST

	write(*,*) 'READING IN: Match Database'
	do i=1,nrsp
	   do j=1,nion
	       sid(i,j)=0.
	   end do
	end do

	i=0
	j=0

130	i=i+1
	read(29,54,end=140,err=900) swave(i),str1,str2,srecord(i)

*	CONVERT TO AIR IF REQUIRED AND INTO ANGSTROMS
	
	if (.not. vacumn) call vac2air(swave(i))
	swave(i)=swave(i)*real(1.0E4)
	
*	DETERMINE ELEMENTAL TYPE AND ION

	sele=elm(str1)
	sion=ion(str1)+1
	sttype=ttype(str2)
	if (sttype.lt.1) sion=sion+1
	
*	ASSIGN TO ARRAY

	sid(sele,sion)=sid(sele,sion)+1
	sided=sid(sele,sion)
	slist(sele,sion,sided)=i
	
	go to 130

*	READ IN THE LIST OF MATCHED LINES

140	continue
	close(29)

	do j=1,nummatched

*	DETERMINE THE SOURCE ION

	  k=index(matchele(j),'[')
	  m=index(matchion(j),']')

	  if (k.ne.0) tele=matchele(j)(2:3)
	  if (k.eq.0) tele=matchele(j)(1:2)
	
	  do n=1,nrsp
	     if (tele.eq.slab(n)) go to 150
	  end do

	  write(*,*) 'ERROR: Cannot understand ',matchele(j)
	  write(*,*) 'on line ',j
	  stop

150	  telenum=n
	  
	  if (m.ne.0) then
             do n=1,6
	        if (matchion(j)(n:n).ne.']') tion(n:n)=matchion(j)(n:n)
	        if (matchion(j)(n:n).eq.']') tion(n:n)=' '
	     end do
	  else
	     tion=matchion(j)
	  end if

	  do n=1,nrsp
	     if (tion.eq.spcs(n)) go to 160
	  end do

155	  write(*,*) 'WARNING: Matchlist3.f err=155' 
	  write(*,*) 'Cannot understand ',matchion(j)	  
	  write(*,*) 'on line ',j
	  stop

160	  tionnum=n	  
	  if ((k.eq.0).and.(m.eq.0)) tionnum=tionnum+1
	  
*	DETERMINE REDSHIFT

	  rvcor=(obswave(j)-labwave(j))*vlight/labwave(j)

*	REFERENCE THE SPECIAL LIST TO LOOK FOR A MATCH

	  do n=1,sid(telenum,tionnum)
	     
	     wavedif=(labwave(j)-swave(slist(telenum,tionnum,n)))*vlight
	     wavedif=abs(wavedif/swave(slist(telenum,tionnum,n)))
             if (wavedif.lt.vtol) go to 170
          end do
	 
	  trecord(j)=999
	  go to 180

170	  trecord(j)=srecord(slist(telenum,tionnum,n))
	  	  
180	  continue

*	DETERMINE WHICH BIN THE LINE BELONGS TO

	   do k=1,(numbins-1)
              if ((ioneng(telenum,tionnum).gt.bins(k)).and.
     1        (ioneng(telenum,tionnum).le.bins(k+1))) go to 190
	   end do

190	   binnum(k)=binnum(k)+1
	   binsum(k)=binsum(k)+rvcor
	end do

*	READ IN THE ABUNDANCE DATA

195	i=0
200	read(22,*) header
	if (header.eq.'!') go to 200
	backspace(22)	
	read(22,*,end=220,err=901) element,abun
 	i=i+1

	if (i.gt.nrsp) go to 220       
	do j=1,nrsp
	   if (element.eq.'H ') habun=abun
	   if (element.eq.'He') y=abun
	   if (element.eq.slab(j)) go to 210
	end do

205	write(*,*) 'WARNING: Matchlist3.f err=205'
	write(*,*) 'Element Not Recognized: ',element
 	stop

210	uabun(j)=abun
	go to 200
220	continue

	if (habun.eq.0) then
	   write(*,*) 'H abundance not specified'
	   write(*,*) 'Please Enter the H abundance: '
	   read(*,*) habun
	end if
	if (habun.ne.1.) then
	   do j=1,nrsp
	      uabun(j)=uabun(j)/habun
	   end do
	   y=y/habun
        end if

	if (y.eq.0) y=y_def

*       IF MATCH LIST NOT BEING EMPLOYED SKIP AHEAD TO ASSIGNMENT
*       OF FACTORS

	if (.not. use_match) go to 255

*	LOOKING TO POPULATE BIN "A"

*	EXAMINE ALL RECORDS TO FIND THE STRONGEST LINE (IF ANY)
*	THAT EXISTS IN THE MATCH LIST FOR THIS BIN
	
   	i=0
	
        do n=1,nummatched
           if (trecord(n).eq.999) go to 230
	   tbin=bin(trecord(n))
           if (tbin.ne.1) go to 230
	   i=i+1
	   tlist(i)=i
	   tflx(i)=obsflx(n)
230	   continue
	end do

        call sort(1,tflx,tlist,i)
	
        if (i.eq.0) then
           ix(1)=min_val(1)
	else
           if (tflx(tlist(1)).lt.thresh(1,1)) ix(1)=val(1,1)
           if (tflx(tlist(1)).gt.thresh(1,1)) ix(1)=val(1,2)
	   if (tflx(tlist(1)).gt.thresh(1,2)) ix(1)=val(1,3)
	end if
	  
*	LOOKING TO POPULATE BIN "E"

*	EXAMINE ALL RECORDS TO FIND THE STRONGEST LINE (IF ANY)
*	THAT EXISTS IN THE MATCH LIST FOR THIS BIN	   

*	RESETTING THE "tflx" AND "tlist" ARRAYS

	do n=1,i
	   tflx(n)=0.
	   tlist(n)=0.
	end do
	i=0

	do n=1,nummatched
           if (trecord(n).eq.999) go to 240
	   tbin=bin(trecord(n))
           if (tbin.ne.5) go to 240
	   i=i+1
	   tlist(i)=i
	   tflx(i)=obsflx(n)
240	   continue
	end do

        call sort(1,tflx,tlist,i)

	ix(5)=min_val(5)
	if (i.ne.0) ix(5)=ix(5)+tflx(tlist(1))
	if (ix(5).gt.max_val(5)) ix(5)=max_val(5)
         

*   	LOOKING TO POPULATE BINS "B","C", AND "D"	   

*	SETTING THE "tflx" AND "tlist" ARRAYS
 
	line31=.false.
	line32=.false.
	line4=.false.
	flx31=0.
	flx32=0.
	flx4=0.
       
	do n=1,nummatched
	   if (trecord(n).eq.999) go to 250
	   tbin=bin(trecord(n))
	   tstrength=strength(trecord(n))
	   if ((tbin.ne.3).and.(tbin.ne.4)) go to 250
	   if ((tbin.eq.3).and.(tstrength.eq.1)) then
              line31=.true.
	      flx31=obsflx(n)
	   end if
	   if ((tbin.eq.3).and.(tstrength.eq.2)) then
              line32=.true.
	      flx32=obsflx(n)
	   end if
	   if ((tbin.eq.4).and.(tstrength.eq.1)) then
	      line4=.true.
	      flx4=obsflx(n)
	   end if
250	   continue
	end do

	ixr(1)=0.
	ixr(2)=0.

	if (line31) ixr(1)=flx31*0.74/y
        if (line32) ixr(1)=ixr(1)+flx32*2.03/y
	if (line31.and.line32) ixr(1)=ixr(1)/2.
	if ((.not.line31).and.(.not.line32)) ixr(1)=999
	if (line4) then
	   if (line31) ixr(2)=0.11*flx4/flx31
	   if (line32) ixr(2)=ixr(2)+0.04*flx4/flx32
	   if (line31.and.line32) ixr(2)=ixr(2)/2.
	   if ((.not.line31).and.(.not.line32)) ixr(2)=999
	else
	   ixr(2)=999
	end if

	remain=1.-ix(1)-ix(5)
	if ((ixr(1).eq.999).and.(ixr(2).eq.999)) then
*       Either everything is fairly low ionization
	   ix(4)=min_val(4)
	   ix(3)=min_val(3)
	   ix(2)=remain-ix(3)-ix(4)
*       Or everthing is rather high ionization
	   if (ix(5).eq.max_val(5)) ix(2)=min_val(2)
	end if
	if ((ixr(1).ne.999).and.(ixr(2).eq.999)) then
*       Probably everything is lower ionization bin 1-3
	   ix(4)=min_val(4)
	   ix(2)=(remain-ix(4))/(1+ixr(1))
	   if (ix(2).lt.min_val(2)) ix(2)=min_val(2)
	   ix(3)=ixr(1)*ix(2)
	   if (ix(3).lt.min_val(3)) ix(3)=min_val(3)
	end if
	if ((ixr(1).ne.999).and.(ixr(2).ne.999)) then
*       Split everything a bit more evenly
	   ix(2)=remain/(1+ixr(1)+ixr(1)*ixr(2))
	   if (ix(2).lt.min_val(2)) ix(2)=min_val(2)
	   ix(3)=ix(2)*ixr(1)
	   if (ix(3).lt.min_val(3)) ix(3)=min_val(3)
	   ix(4)=ix(3)*ixr(2)
	   if (ix(4).lt.min_val(4)) ix(4)=min_val(4)
	end if
	if ((ixr(1).eq.999).and.(ixr(2).ne.999)) then
*       Bit higher ionization
	   ix(2)=min_val(2)
	   ix(4)=(remain-ix(2))/((1/ixr(2))+1)
	   if (ix(4).lt.min_val(4)) ix(4)=min_val(4)
	   ix(3)=ix(4)/ixr(2)
	   if (ix(3).lt.min_val(3)) ix(3)=min_val(3)
	end if

*	IF use_def_icf IS SPECIFIED THEN USE DEFAULT VALUES
*       GIVEN IN VARIABLE DECLARATION ABOVE

 255	continue

	if (use_def_icf) then
           write(*,*) 'WARNING: Using Default ICF'
           do ii=1,numbins
              ix(ii)=def_icf(ii)
	   end do
        end if

*       ELSE USE THE MANUALLY SPECIFIED
*       VALUES PROVIDED IN OPENALL, STORED IN "icf" ARRAY
*       IF NO MATCHED FILE IS SPECIFIED

	if (manicf) then
	   write(*,*) 'WARNING: Using Manually Specified ICF'
	   do ii=1,numbins
	      ix(ii)=icf(ii)
	   end do
	end if

*	RECORD TO OUTPUT FILE

	write(28,*) 'ICF Values: Bin/%'
	write(28,*) 'ix 1: ',ix(1)
	write(28,*) 'ix 2: ',ix(2)	
	write(28,*) 'ix 3: ',ix(3)
	write(28,*) 'ix 4: ',ix(4)
	write(28,*) 'ix 5: ',ix(5)
	write(28,*)

*       IF use_def_rvcor IS SPECIFIED THEN USE DEFAULT VALUES
*       GIVEN IN VARIABLE DECLARATION

        if (use_def_rvcor) then
           write(*,*) 'WARNING: Using Default Velocity Matrix'
           do ii=1,numbins
              binsum(ii)=def_rvcor(ii)
	   end do
	   go to 258
	end if

*       ELSE USE THE MANUALLY SPECIFIED
*       VALUES PROVIDED IN OPENALL, STORED IN "vel" ARRAY

	if (manvel) then
	   write(*,*) 'WARNING: Using Manually Specified Velocity Matrix'
	   do ii=1,numbins
	      binsum(ii)=vel(ii)
	   end do
	   go to 258
	end if

*       ELSE PLUNGE AHEAD AND CALCULATE VELOCITY MATRIX

*	DETERMINE AVERAGE REDSHIFT IN EACH BIN

	do i=1,numbins
	   if (binnum(i).ne.0) then
              binsum(i)=binsum(i)/binnum(i)
	   else
	      binsum(i)=999
	   end if
	end do

*	IF NO VALUES IN A PARTICULAR BIN THEN
*	ASSIGN THE VALUE IN THE NEAREST BIN

	do i=1,numbins
	   if (binsum(i).eq.999) then
	      mindif=99
              do j=1,numbins
	         if ((binsum(j).ne.999).and.(i.ne.j)) then
                   if (abs(i-j).lt.mindif) then
                      binsum(i)=binsum(j)
                      mindif=abs(i-j)
	           end if
	         end if
	      end do
	   end if
	end do

*	IF STILL NO VALUES IN A PARTICULAR BIN OR BY DESIGN
*      (use_def_rvcor=.true.) THEN ASSIGN A DEFAULT VALUE

	do i=1,numbins
	   if ((binsum(i).eq.999).or.(use_def_rvcor)) binsum(i)=def_rvcor(i)
	end do

 258	continue

	write(*,*) 'SETTING UP: Abundance Matrix'
	write(*,*) 'SETTING UP: Velocity Distribution Matrix '
c	if (use_def_rvcor) write(*,*) '      --- DEFAULT SETTING'
	write(*,*)

*	FLUSH OUT THE IVCOR(i,j) MATRIX, SET FLAGS

	do i=1,nrsp
	   do j=1,(nrsp+1)
	      ivcor(i,j)=-9999
	   end do
	end do

*	READ IN IONIZATION POTENTIALS

	do i=1,nrsp
           do j=1,(nrsp+1)
	      do k=1,(numbins-1)
                 if ((ioneng(i,j).gt.bins(k)).and.(ioneng(i,j).le.
     1	         bins(k+1))) go to 260
	      end do
260	      continue
	      ivcor(i,j)=binsum(k)
              if (j.eq.1) prevbin=1
	      if (j.gt.iions(i)) then
	         iabun(i,j)=ix(numbins)*uabun(i)
	         ivcor(i,j)=binsum(numbins)
c	         if (i.eq.12 .and. j.eq.10) print *,'!! ',ix(numbins),
c     1           uabun(i),binsum(numbins)
	         k=prevbin
	      else
	         if (k.eq.prevbin) then
	            iabun(i,j)=ix(k)*uabun(i)
	         else
	            iabun(i,j)=(ix(prevbin)+ix(k))*uabun(i)/2
	         end if
	      end if
	      prevbin=k
	   end do
	end do

*	SPECIAL CASE OF COMPLETELY IONIZED HYDROGEN AND HELIUM
	
	iabun(1,2)=ix(2)*uabun(1)
	iabun(2,2)=ix(3)*uabun(2)
	iabun(2,3)=ix(4)*uabun(2)

*       IF DEPELETIONS EXIST FROM OPENALL, APPLY THEM HERE

	if (depelenumber.gt.0) write(28,*) 'Elements Depeleted: 
     1    Ionization States    Amount'
	do i=1,depelenumber
	   do j=depelelow(i),depelehi(i)
	      iabun(depele(i),j)=iabun(depele(i),j)/depeleamount(i)
	   end do   
	   write(28,55) slab(depele(i)),spcs(depelelow(i)),
     1        spcs(depelehi(i)),depeleamount(i)
	end do
	if (depelenumber.gt.0) write(28,*)
	
*	RECORD TO OUTPUT FILE

	write(28,*) 'Velocity Structure: Bin/Vel(km/s)'
	write(28,*) 'irvcor 1: ',binsum(1)
	write(28,*) 'irvcor 2: ',binsum(2)	
	write(28,*) 'irvcor 3: ',binsum(3)
	write(28,*) 'irvcor 4: ',binsum(4)
	write(28,*) 'irvcor 5: ',binsum(5)
	write(28,*)

**	END SUBROUTINE

	return

*	ERROR INDICATED IN INPUT LINE LIST
900	write(*,*) 'WARNING: Matchlist.f  err=900'   
	write(*,*) 'Error in the input list at line: ',j
	stop 
    
*	ERROR IN THE SPECIAL ASSIGNMENT LIST
901	write(*,*) 'WARNING: Matchlist.f: err=901'	
	write(*,*) 'Error Reading in Special Assignment List'
	write(*,*) 'Error occured in line: ',i
	stop

*	ERROR IN ABUNDANCE DATA
902	write(*,*) 'WARNING Matchlist.f: err=902'
	write(*,*) 'Error Reading in the Abundance Table'
	write(*,*) 'Error occured in line: ',i
	stop 
	
*	ERROR IN READING PREVIOUSLY UNFOUND IONIZATION STAGE
903	write(*,*) 'WARNING Matchlist.f: err=903'
	write(*,*) 'Error in Expected Format of Manual Input'
	stop
	
  	end	 
	
