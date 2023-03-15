	subroutine alphasort(stat,ins,ids,num)
*	sorts by criteria in 'in' array
*	and returns the list of corresponding ids
	
	integer num,stat
	character*2 ins(num)
	integer ids(num)
	character*2 in(10500)
	integer swapno
	character*2 swapchar 
	
	do i=1,num
	   in(i)=ins(i)
	end do

	do i=1,num
	   do j=(i+1),num
	      if (in(i).gt.in(j)) then       
	         swapchar=in(i)
	  	 swapno=ids(i)
       		 in(i)=in(j)
	      	 ids(i)=ids(j)
	         in(j)=swapchar
	         ids(j)=swapno
	      end if
  	   end do
	end do
 
	return
	end
	            
