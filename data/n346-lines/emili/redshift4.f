
	SUBROUTINE redshift(itr1,itr2,loleveng,orvcor)

******* This subroutine sets up the velocity look up
*******	table to get trial corrections to the observed
******* wavelenght of a candidate line depending upon
******* the specified ionization table.  Lines ID'd in
******* the main program are added to improve the table.
*
*	INPUT:
*	itr1	- incoming "tr1" specifier from the master
*	          line list for a putative ID - INTEGER
*	itr2	- incoming "tr2" specifier from the master
*	          line list for a putative ID - INTEGER
*
*	OUTPUT:
*	orvcor  - the velocity correction for the appropriate
*	          for the particular source ion's (specified
*	          by the "tr1" and "tr2" specifiers) ionization
*	          potential - REAL

	IMPLICIT NONE

*	INPUT VARIABLES:

	INTEGER itr1,itr2
	REAL loleveng

*	OUTPUT VARIABLES:

	REAL orvcor

*	OTHER VARIABLES

	INTEGER k,inele,inion,intrt,elm,ion,ttype,nrsp,nion
	PARAMETER (nrsp=30,nion=10)
	REAL ivcor(nrsp,nrsp+1)

*	COMMON BLOCKS:

	COMMON /velcor/ivcor

*	FUNCTIONS:

	elm(k)=iand(ishft(k,-26),63)
	ion(k)=iand(ishft(k,-20),63)	
	ttype(k)=iand(ishft(k,-4),7)

*	DECOMPOSE INCOMING "tr1" and "tr2" SPECIFIERS

	inele=elm(itr1)
	inion=ion(itr1)+1
	intrt=ttype(itr2)

*	IF SOURCE TRANSITION IS RECOMBINATION (dipole , non-intercombination,
*	ttype(itr2)=0) THEN BUMP UP TO NEXT STAGE

	if (intrt.lt.1) inion=inion+1

*       IF SOURCE TRANSITION IS INTERCOMBINATION WITH AN ENERGY > 0.372 EV
*       IN IT'S DESTINATION LEVEL, ASSUME ALSO THAT IS RECOMBINATION

	if (intrt.eq.1 .and. loleveng.gt.0.372) inion=inion+1

*	READ IN AND ASSIGN THE VELOCITY CORRECTION

	orvcor=ivcor(inele,inion)
	if (orvcor.eq.-9999.) then
	    write(*,*) 'Redshift.f  WARNING: '
	    write(*,*) 'Cannot find a velocity correction for element: ',inele
	    write(*,*) 'Ionization Stage: ',inion
	    write(*,*) 'Please enter a velocity correction for this ion'
	    read(*,*,err=900) orvcor
	end if

**	END SUBROUTINE

	return

*	ERROR IN READING MANUAL INPUT
900	write(*,*) 'WARNING: Redshift.f: err=900'
	write(*,*) 'Error in Expected Format of Manual Input'
	stop
	
	end	 

	   
	
	
	
