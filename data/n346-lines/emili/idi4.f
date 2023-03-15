
	SUBROUTINE idi(iflux,plet1,plet2,veldif,velsig,num,inids,pts,outids)

******	Routine to rank a series of putative IDs from the 
******  master line list by a parameter known as the "ID" Index
******  which is calculated here based upon the results of consistency
******  checks, as well as the wavelength agreement with the candidate
******  line attempting to be ID'd and the flux strength.
*
*	INPUT:
*	iflux	- array of flux strengths for the putative IDs
*	  	- ARRAY OF REAL
*	plet1   - array of the number of multiplet lines associated 
*                 with the putative line ID that
*                 should also be seen in the spectra 
*	        - ARRAY OF INTEGER
*	plet2   - array of the number of multiplet lines associated
*	 	  with each putative line ID that were seen in the
*		  spectra - ARRAY OF INTEGER
*	veldif  - array of the velocity difference between the
*		  candidate line and the putative line ID
*		- ARRAY OF REAL
*	velsig	- the sigma of the candidate lines wavelength
*		  measurement in km/s - ARRAY OF REAL(2)
*	num	- the number of putative line IDs for the candidate
*	          line - INTEGER
*	inids   - array of the ID # for each of the putative IDs
*	        - ARRAY OF INTEGER
*
*	OUTPUT:
*	pts	- array of ICI index points for each putative ID
*	        - ARRAY OF INTEGER
*	outids	- array of the ID # for each of the putative IDs
*	 	  ranked according to the ICI index

	IMPLICIT NONE

*	INPUT VARIABLES:

	INTEGER num
	REAL iflux(num),veldif(num),velsig(2)
	INTEGER inids(num),plet1(num),plet2(num)
	 
*	OUTPUT VARIABLES:

	INTEGER outids(num)
	REAL pts(num)

*	OTHER VARIABLES:

	INTEGER i,j,k,work(num)
	REAL usevel(num),topflx,wpts(num)

*	PARAMETER VARIABLES:
	
	REAL tol1,tol2,tol3
	PARAMETER (tol1=1)
	REAL ttol1,ttol2,ttol3
	PARAMETER (ttol1=10.)
	INTEGER t1,t2
	REAL ptsa(5)
	DATA (ptsa(i),i=1,5)/0,1,2,3,4/	

**	BEGIN SUB-ROUTINE

*	RANK THE FLUXES

	do i=1,num
c	   print *,'ID ',inids(i),iflux(i),veldif(i),plet1(i),plet2(i)
c	   print *,velsig(1),velsig(2)
c	   read(*,*)
	   work(i)=i	
	   wpts(i)=0.
	end do

	call sort(1,iflux,work,num)
	topflx=iflux(work(1))	

	j=1
	do i=1,num
	   if (veldif(work(i)).le.0) usevel(work(i))=velsig(1)
	   if (veldif(work(i)).gt.0) usevel(work(i))=velsig(2)
	   if (i.eq.1) then
	      if (num.ne.1) then
	         if ((iflux(work(1))/iflux(work(2))).lt.ttol1) then
                    wpts(work(i))=ptsa(2)
	         else
	            wpts(work(i))=ptsa(1)
	         end if
	      else
	         wpts(work(i))=ptsa(1)
	      end if
	   else 
	      if (alog10(topflx/iflux(work(i))).gt.j) j=j+1
c       	write(*,*) alog10(topflx/iflux(work(i))),j,topflx,iflux(work(i))
              wpts(work(i))=ptsa(j+1)
	   end if
c	   write(*,*) 'IDI: ',work(i),wpts(work(i)),iflux(work(i))
c	   read(*,*)
	end do

*	RANK THE WAVELENGTH DIFFERENCES

	do i=1,num
	   work(i)=i	
	end do

	call sort(1,veldif,work,num)

	do i=1,num
	   if (abs(veldif(work(i))).gt.abs(0.5*usevel(work(i)))) 
     1     wpts(work(i))=wpts(work(i))+1
	   if (abs(veldif(work(i))).gt.abs(1.0*usevel(work(i)))) 
     1     wpts(work(i))=wpts(work(i))+1
           if (abs(veldif(work(i))).gt.abs(1.5*usevel(work(i)))) 
     1     wpts(work(i))=wpts(work(i))+1
	   if (abs(veldif(work(i))).gt.abs(2.0*usevel(work(i))))
     1     wpts(work(i))=wpts(work(i))+1
c           write(*,*) 'IDI: ',work(i),wpts(work(i)),veldif(work(i))
c	   read(*,*)	
        end do

*	RANK THE MULTIPLET STATISTICS

	do i=1,num
	   t1=plet1(i)
	   t2=plet2(i)
	   if (t1.eq.0) wpts(i)=wpts(i)+1
	   if ((t1.eq.2).and.(t2.eq.1)) wpts(i)=wpts(i)+1
           if ((t1.eq.1).and.(t2.eq.0)) wpts(i)=wpts(i)+2
	   if ((t1.gt.2).and.(t2.eq.1)) wpts(i)=wpts(i)+2
	   if ((t1.gt.1).and.(t2.eq.0)) wpts(i)=wpts(i)+3
c           write(*,*) 'IDI: ',i,wpts(i),t1,t2
c	   read(*,*)
	end do

*	SORT THE POINTS

	do i=1,num
	   work(i)=i
	end do

	call sort(-1,wpts,work,num)

*	PREPARE OUTPUT

	do i=1,num
	   outids(i)=inids(work(i))
	   pts(i)=wpts(work(i))
c	   write(*,*) 'IDI :',work(i),outids(i),pts(i)
c	   read(*,*) 
	end do

	return

	end
	
	    
	   	
			

	      
