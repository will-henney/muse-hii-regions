	SUBROUTINE ionize(stat,itr1,itr2,pioneng,ionstate)

******	This routine reads in the ionization energy of all elements (stat=1)
******  When queried (stat=2) the routine returns the ionization energy
******  (pioneng) of the paritcular line.  If it a forbidden line
******  it returns the ionization energy of the level itself, while
******  if its a recombination line, it returns the next level's ionization
******  energy.  It also returns the numerical equiavalent of the ionization
******  state of the ion description (ion)
*
*	INPUT:
*	stat	- status call to routine: 1-initial read 2-query - INTEGER
*	itr1    - input "tr1" specifier of the master list line - INTEGER
*	itr2    - input "tr2" specifier of the master list line - INTEGER
*	
*	OUTPUT:
*	pioneng - the ionization energy of the parent species of
*	          master list line - REAL
*	ionstate - the numeric conversion of the master list lines
*	           ionization state - INTEGER

*	CODING COMPLETED: 8/17

	IMPLICIT NONE

*	INPUT VARIABLES:

        INTEGER stat,itr1,itr2
	
*	OUTPUT VARIABLES:

        INTEGER ionstate
	REAL pioneng

*	OTHER VARIABLES:
	
	INTEGER nrsp,inele,nion,inion,intrt
	PARAMETER (nrsp=30,nion=10)
	INTEGER iions(nrsp)
	REAL ioneng(nrsp,nion)
	INTEGER i,j,k,l,z,ier
	CHARACTER*1 header
	CHARACTER*2 ele(nrsp)
 	CHARACTER*72 bigheader       
	INTEGER elm,ion,ttype

*	COMMON BLOCKS:

        COMMON /ionization/ioneng,iions

*	FUNCTIONS:

        elm(k)=iand(ishft(k,-26),63)
	ion(k)=iand(ishft(k,-20),63)
	ttype(k)=iand(ishft(k,-4),7)

*	FORMATS:

50	format(a2,8(f7.2))

*	BEGIN SUBROUTINE:

**	BEGIN SETUP "stat=1" CALL

	if (stat.eq.1) then

	i=0
	j=0	

*	READ IN THE IONIZATION TABLE

100	i=i+1
	ier=1
	read(21,*,end=110,err=900) header
	if (header.eq.'#') go to 100
        backspace(21)
	ier=2
	read(21,'(a72)',err=900) bigheader
	z=index(bigheader,']')
	z=int((z-3)/7)	
	j=j+1
	backspace(21)
	ier=3
	read(21,50, err=900) ele(j),(ioneng(j,l),l=1,z)
	iions(j)=z
	if (j.eq.nrsp) go to 110
	go to 100

110	continue
	close(21)

*	ASSIGN PADDED VALUES TO COMPLETELY IONIZED H AND HE

	ioneng(1,2)=24.69	
	ioneng(2,3)=99.
	iions(1)=iions(1)+1
	iions(2)=iions(2)+1

	write(*,*) 'SETTING UP: Ionization Energy Matrix'

**	END SETUP "stat=1" CALL

c	do i=1,nrsp
c	   write(*,*) i,ele(i),(ioneng(i,j),j=1,numions(i))
c	   read(*,*)
c	end do

	end if    	   	   

**	BEGIN "stat=2" REFERENCE CALL
      
	if (stat.eq.2) then

*	DECOMPOSE INCOMING "tr1" SPECIFIER

	inele=elm(itr1)
	inion=ion(itr1)+1
        intrt=ttype(itr2)

*	ASSIGN OUTPUT

	ionstate=inion
	pioneng=ioneng(inele,ionstate)
	if (intrt.lt.1) pioneng=ioneng(inele,ionstate+1)
	if (itr1.eq. 814750725) print *,'Ionize: ',inele,inion,pioneng
c	read(*,*)

**	END "stat=2" REFERENCE CALL

	end if

	return

*	ERROR IN READING THE ION DATA TABLE
900	write(*,*) 'Ionize.f  err=900 specific: ',ier
	write(*,*) 'Error reading ion table at line: ',i
	stop

	end

   
	   

	   
	      	
