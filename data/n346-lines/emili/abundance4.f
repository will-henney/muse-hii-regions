
	SUBROUTINE abundance(itr1,oabun)

******	This routine determines a rough ionization distribution of the 
******  observed nebulae using matched lines provided by the user (read
******  in by routine "matchlist.f") with known labratory wavelengths.
******  Specific lines are searched for in this input list, which establish
******  the ICF for each element, based upon abundances read in from the
******  data file "abun.dat" which are with respect to hydrogen.  This
******  routine builds an array of ICF versus ion (in "iabun" array)
******  which can be accessed by the main program
*
*	INPUT:
*	itr1	-the "tr1" specifier of the ion producing a line from the
*                master line list - INTEGER
*
*	OUTPUT:
*	obaun   -output abundances for the master list ion and next higher
*                ionization stage for submission to relative intensity
*	         calculations - ARRAY OF REAL(2)

*	CODING COMPLETED: 8/17

	IMPLICIT NONE

*	INPUT VARIABLES:
	
	INTEGER itr1

*	OUTPUT VARIABLES:

	REAL oabun(2)

*	OTHER VARIABLES:

	INTEGER j,k,inele,inion,elm,ion,nrsp,nion
	PARAMETER (nrsp=30)
	REAL iabun(nrsp,nrsp+1)

*	COMMON BLOCKS:

	COMMON /abun/iabun

*	FUNCTIONS:

	elm(k)=iand(ishft(k,-26),63)
	ion(k)=iand(ishft(k,-20),63)	

*	DECOMPOSE INCOMING "tr1" SPECIFIER

	inele=elm(itr1)
	inion=ion(itr1)+1

c	write(*,*) inele,inion,(iabun(inele,j),j=1,10)

*	LOOK UP THE ABUNDANCE OF THE ION AND ITS NEXT STAGE OF IONIZATION
*	ASSIGN THOSE VALUES

	oabun(1)=iabun(inele,inion)
	oabun(2)=iabun(inele,inion+1)

	if ((oabun(1).eq.0).or.(oabun(2).eq.0)) then
	    write(*,*) 'WARNING: Abundance.f '
	    write(*,*) 'Cannot find an abundance for this element: ',inele
	    write(*,*) 'Ionization Stage: ',inion
	    write(*,*) 'Please enter the abundances wrt to H for this ion'
            write(*,*) 'and the next higher ionization stage:'
	    read(*,*,err=900) oabun(1),oabun(2)
	end if

**	END SUBROUTINE

	return

*	ERROR IN READING IN MANUAL INPUT
900	write(*,*) 'WARNING: Abundance.f: err=900'
	write(*,*) 'Error in Expected Format of Manual Input'
	stop
	
	end	 
		

	

	    
	    


	
	
	
	
	



	
