	
	real wave(2),ener(2),obss(1500),flxs(1500)
	logical usels(1500)
	integer itr1,itr2,stat
	real*8 wav(1500)
	character*1 header,mark,marks(1500)
	character*3 ele(1500),wele(50),multia,multi(50)
	character*9 ion(1500),wion(50)
	real upeng(1500),loeng(1500),ieng
	character*11 uplev(1500),lolev(1500)
	real ajk(1500),upj(1500),loj(1500)
	real lastobs,lastfwhm,lastsn,lastflux
	real sn,fwhm,flux,obs
	character*2 tts(1500)
	integer tr1(1500),tr2(1500),istage,nlines
	integer i,j,k,m,in1,in2,elenum,argind
	character*40 ozz,ozz2
	character*30 wlolev(50),wuplev(50)
	real wlab(50),vel(50),wupj(50),wloj(50)
	character*2 ida,wid(50)
	character*70 line,cmdarg
	character*20 infile,outfile

*       READING IN THE COMMAND LINE ARGUEMENTS

        if (iargc().ne.2) then
	   write(*,*) 'EMREAD.F ERROR: Need to specifiy an output file'
	   stop
	end if

        call getarg(1,cmdarg)
        argind=index(cmdarg,' ')-1
        read(cmdarg(1:argind),'(a20)',err=9996) infile
	call getarg(2,cmdarg)
        argind=index(cmdarg,' ')-1
	read(cmdarg(1:argind),'(a20)',err=9997) outfile

 500	continue
	
*       OPENING ALL FILES

	open(23,file='slev_list')
	open(24,file='sline_list')
        open(11,file=infile,err=9998)
	open(12,file=outfile,err=9999)
	open(13,file='intermediate',err=9999)

	stat=0.

	i=0
 700	i=i+1
	read(13,'(f10.3,f10.4)',end=800) obss(i),flxs(i)
	go to 700
 800	continue
	
*       SET UP THE LINE LIST

	call line_pick(stat,wave,ener,itr1,itr2,wav,ele,ion,upeng,loeng,
     1  uplev,lolev,ajk,upj,loj,tts,tr1,tr2,usels,nlines,ozz,ozz2)

	write(*,*) 'Routine to set-up complete'
	
*	READING IN EMILI FULL OUTPUT LIST

1000	continue

	do i=1,25
	   read(11,'(a1)') header
	end do

 80	format(f9.3,1x,f5.1,1x,f8.4,1x,f8.1,$)
 81	format(2x,a2,1x,a3,1x,a5,f9.3,1x,a1,1x,f5.1,2x,a25,1x,a25,1x,f4.1,
	1     1x,f4.1,2x,a3)
 82	format(35x,a2,1x,a3,1x,a5,f9.3,1x,a1,1x,f5.1,2x,a25,1x,a25,1x,f4.1,
	1      1x,f4.1,2x,a3)
 83	format(15x,f8.2,1pe9.2,0p,6x,f9.2,6x,f6.1)
 84	format(16x,f9.3,a1,10x,i11,i7,16x,a3,1x,a2)

	j=0
	m=0

 2000	continue

	read(11,'(a1)',end=9000,err=9001) header
	
	if (header.eq.'O') then
	   if (m.ne.0) then
	      write(12,80) obs,fwhm,flx,sn
	      if (j.gt.0) then
		 do k=1,j
		    if (k.eq.1) write(12,81) wid(k),wele(k),wion(k),
     1              wlab(k),marks(k),vel(k),wlolev(k),wuplev(k),wloj(k),
     2              wupj(k),multi(k)
	            if (k.ne.1) write(12,82) wid(k),wele(k),wion(k),
     1              wlab(k),marks(k),vel(k),wlolev(k),wuplev(k),loj(k),
     2              wupj(k),multi(k)
	         end do
	      else
		 write(12,*)
              end if
	   end if
	   m=m+1
	   backspace(11)
	   read(11,83) obs,flx,sn,fwhm
	   obs=obss(m)
	   flx=flxs(m)
           j=0
	else if (header.eq.'*') then
	   backspace(11)
	   read(11,84) wave(1),mark,itr1,itr2,multia,ida
	   stat=192
	   ener(1)=99.
	   ener(2)=99.
	   wave(1)=wave(1)-1.
	   wave(2)=wave(1)+2.

	call line_pick(stat,wave,ener,itr1,itr2,wav,ele,ion,upeng,loeng,
     1  uplev,lolev,ajk,upj,loj,tts,tr1,tr2,usels,nlines,ozz,ozz2)

	   j=j+1
	   
	   marks(j)=' '
	   wid(j)=ida
	   multi(j)=multia
	   if (mark.eq.'*') marks(j)=mark

           do i=1,nlines
	      if (tr1(i).eq.itr1 .and. tr2(i).eq.itr2) then
c		 print *,ele(i),ion(i),'!',ozz,'!'
c		 print *,ele(i),ion(i),'!',ozz2,'!'
c		 print *,ele(i),ion(i),'!',uplev(i),'!'
c		 print *,ele(i),ion(i),'!',lolev(i),'!'
	         it1=index(ozz,' ')
	         it2=index(ozz2,' ')
	         iz1=index(uplev(i),' ')-1
	         iz2=index(lolev(i),' ')-1
	         wuplev(j)=uplev(i)(1:iz1)//' '//ozz(1:it1)
	         wlolev(j)=lolev(i)(1:iz2)//' '//ozz2(1:it2)
c		 print *,ele(i),ion(i),'!',wuplev(j),'!'
c		 print *,ele(i),ion(i),'!',wlolev(j),'!'
	         wele(j)=ele(i)
		 wion(j)=ion(i)
		 wlab(j)=wav(i)
		 vel(j)=(wlab(j)-obs)*300000/wlab(j)
		 wloj(j)=loj(i)
		 wupj(j)=upj(i)
		 if (.not. usels(i)) multi(j)='   '
	      end if
	   end do
	end if 
	go to 2000

 9000	continue
	write(12,80) obs,fwhm,flx,sn
	if (j.gt.0) then
	   do k=1,j
	      if (k.eq.1) write(12,81) wid(k),wele(k),wion(k),
     1        wlab(k),marks(k),vel(k),wlolev(k),wuplev(k),wloj(k),
     2        wupj(k),multi(k)    
	      if (k.ne.1) write(12,82) wid(k),wele(k),wion(k),
     1        wlab(k),marks(k),vel(k),lolev(k),uplev(k),loj(k),wupj(k),
     2        multi(k)
	   end do
	else
	   write(12,*)
	end if
	write(*,*) 'Program Completed'
	close(12)
	close(11)
	stop

 9001	write(*,*) 'Error in Reading'
	stop

 9996	write(*,*) 'Error1'
	stop

 9997	write(*,*) 'Error2'
	stop

 9998	write(*,*) 'Error3'
	stop

 9999	write(*,*) 'Error4'
	end
	
	
	
