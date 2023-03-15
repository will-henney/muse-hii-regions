	
	SUBROUTINE ri(rabun,rele,rion,ristat,rtt,rajk,rupeng,rloeng,
     1  tr1,t_e,n_e,riv)

*******	This subroutine calculates the predicted relative intesity
******* of a line relative relative to a hyrdrogen recombination
******* line
*
*	INPUT:
*	rabun 	-abundence of lines parent ion relative to hydrogen
*	rele	-parent element description as drawn from master
*	         list (includes intercombination info)
*	rion	-parent element ion description as drawn from
*	         master list
*	ristat  -ionization state from which line arrises
*	rtt     -transition type originating line
*	rajk    -Einstein coefficient for spontaneous decay
*	rupeng  -energy of upper level of transition
*       rloeng  -energy of lower level of transition
*       n_e     -electron density
*	t_e     -electron temperature
*	
*	OUTPUT:
*	riv	-relative intensity

	real rabun(2)
	real rajk,rupeng,rloeng,n_e,t_e
	character*3 rele
	character*9 rion
	integer ristat,tr1
	character*2 rtt
	real riv
	real lajk,kk,lriv1,lriv2,dd,cc
	logical isinter
	integer itr1

	isinter=.false.

*   Determine if the line is an intercombination line

        if ((index(rele,'[').eq.0).and.(index(rion,']').ne.0))
     1  isinter=.true.

*   Assign it a rough norm. value for the einstein coefficient,
*   and calculate the "K" factor

*	If it is electric diploe transition (tt=E1)
*	and not a intercombination line

	if ((rtt.eq.'E1').and.(.not. isinter)) then
           lajk=1.
	   kk=1.0E-14
	   cc=1.

*	If it is an intercombination line...

	else if (isinter) then 
           lajk=1.0E-4
	   kk=1.0E-9
	   if (rloeng.le.0.372) cc=1.
	   if (rloeng.gt.0.372) cc=0.03

*	If it is a electric quad or magnetic dipole

	else
	   lajk=1.0E-7
	   kk=1.0E-6
	   cc=1.
	end if

*	If the ion is neutral de-value the collisonal term
*	as neutral ions have coulomb repulsion working against excitation
*	Call the pre-factor "dd"

	dd=(ristat-1)*1.0E5
	
	if (ristat.eq.1) dd=1.0e4

*   Calculate relative intesity

        lriv1=rabun(1)*cc*dd*exp(-0.8*rupeng/t_e)/(1+kk*n_e)
        lriv2=rabun(2)*lajk*(ristat**1.7)/t_e**0.7
c	if (tr1.eq.470838288 .or. tr1.eq.472020079) then
c	if (tr1.eq.738199553) then
c	   print *,tr1,rabun(1),ristat,cc,dd,rupeng,t_e,kk,n_e
c	   print *,lriv1
c	   print *,rabun(2),lajk
c	   print *,lriv2
c	   read(*,*)
c	end if

*   Send to main program

        riv=lriv1+lriv2

        return

	end


