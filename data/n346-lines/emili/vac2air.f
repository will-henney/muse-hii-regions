	SUBROUTINE vac2air(wave)
	
******	Convert vacumn wavelength to air
******  Provided by P. van Hoof
******  From Peck et al. (1972)
******  Valid from 1850 - 17000 Angstroms
*
*	INPUT:
*	wave	- input vacumn wavelength
*
*	OUTPUT:
*	wave	- output vacumn wavelength


*	INPUT VARIABLES:

	REAL*8 wave,inwave

*	OTHER VARIABLES:
	
	REAL*8 wavnr,ref_ind

**      BEGIN FUNCTION

	wavnr=1.D0/wave

c	print *,'I am here'

* 	CALCULATE THE REFRACTION INDEX

	ref_ind = 1.d0 + 8060.51d-8 + 2480990d-8/(1.32274d2-wavnr**2)
c	write(*,*) wave,wavnr,ref_ind
	ref_ind = ref_ind + 17455.7d-8/(3.932957d1 - wavnr**2)
c	write(*,*) ref_ind

*	CORRECT THE VACUMN WAVELENGTH TO AIR

	wave=wave/ref_ind

** 	END SUBROUTINE

	return
	end
