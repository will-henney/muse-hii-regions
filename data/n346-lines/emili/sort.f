	subroutine sort(stat,ins,ids,num)
*	sorts by criteria in 'in' array
*	and returns the list of corresponding ids
	
	integer num,stat,i,j
	integer ids(num)
	real ins(num)
	real in(1500)
	integer swapno
	real swap,its
	
	do i=1,num
	   in(i)=ins(i)*real(stat)
	end do

	its=in(1)

	do i=1,num
	   do j=(i+1),num
	      if (in(i).lt.in(j)) then
                 swap=in(i)
	  	 swapno=ids(i)
       		 in(i)=in(j)
	      	 ids(i)=ids(j)
	         in(j)=swap
	         ids(j)=swapno
	      end if
  	   end do
	end do
   
	return
	end
	            
