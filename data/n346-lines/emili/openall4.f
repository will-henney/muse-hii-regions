        SUBROUTINE openall

******  This subroutine reads in a command list from
******  a file specified on the command line used to
******  start EMILI, sets certain parameters, checks
******  that certain conditions necessary for the operation
******  of EMILI are met, then opens the necessary files

*       CODE COMPLETED 5/15/02

        IMPLICIT NONE
      
        CHARACTER*20 commandlist
        CHARACTER*1 dummy
        CHARACTER*60 commandline
        INTEGER line,numele,element
        CHARACTER*20 linele(6)
        INTEGER lgt(6),tempint
        LOGICAL badline
        INTEGER ii,iii,ind,jj,kk
        REAL realarg(3),temparg,sum
        CHARACTER*30 attlabel(8)
        CHARACTER*1 att(8),res
        CHARACTER*20 attarg(5)
        INTEGER depele(30),depelelow(30),depelehi(30),depelenumber
        REAL depeleamount(30)
        INTEGER fatalerrs,nfatalerrs
        INTEGER ferrlinenum(15),nferrlinenum(30)
        INTEGER attlabellen(8)
        CHARACTER*60 ferrline(15),nferrline(30)
        CHARACTER*65 ferrmessage(15),nferrmessage(15)
        INTEGER attarglen(5)
        REAL t_e,n_e,instres,width,kboltz,vlight,hc
        REAL vel(5),icf(4)
        CHARACTER*2 slab(30),tempele
        CHARACTER*6 spcs(30)
        LOGICAL fatal,non_fatal,velset,calcvel,icfset,calcicf
        LOGICAL use_def_icf,use_def_rvcor,use_match
        LOGICAL attset(8)
        LOGICAL manvel,manicf

        DATA fatal,non_fatal/.false.,.false./
        DATA velset,icfset,calcvel,calcicf/4*.false./
        DATA (attset(ii),ii=1,8,1)/8*.false./
        DATA (attlabel(ii),ii=1,8,1)/'Input Line List', 
     1  'Input Matched List','Results List','Short Results List',
     2  'Abundance Table','Electron Temp','Electron Density',
     3  'Inst. Resolution'/
        DATA (attlabellen(ii),ii=1,8,1)/15,18,12,18,15,13,16,16/
        DATA (att(ii),ii=1,8,1)/'L','M','O','D','A','T','N','I'/
        DATA fatalerrs,nfatalerrs/0.,0./
        DATA jj,kk/0.,0./

*       COMMON BLOCKS:

        COMMON /nebatt/t_e,n_e,kboltz,vlight,hc
        COMMON /spc_not/slab,spcs
        COMMON /instrument/width,instres
        COMMON /calculates/use_def_rvcor,use_def_icf,use_match
        COMMON /calculates2/manicf,manvel
        COMMON /specified/vel,icf
        COMMON /depleted/depeleamount,depele,depelehi,depelelow,
     1  depelenumber

***   BEGIN SUBROUTINE

***   Setting Initial Switches in Common Blocks

      use_def_rvcor=.true.
      use_def_icf=.true.
      use_match=.true.
      manvel=.false.
      manicf=.false.

***   Obtaining the command list name from the command line

        call getarg(1,commandlist)
        open(20,file=commandlist,status='OLD',err=9902)

***   Reading in the command list, skipping those lines with a 'Z'

        write(*,*) 'READING IN: Command List'
        write(*,*)

 1000   continue
        line=line+1
        read(20,'(a1)',end=3000,err=9903) dummy
        if (dummy.eq.'Z') go to 1000
        backspace(20)
        read(20,'(a60)',err=9903) commandline
*       splitting the components of the lines
        call split(commandline,linele,lgt,numele,badline)
*       if "split" detect a bad line trigger an error message
        if (badline) go to 9000
*       junction, split off depending number of elements on each line

        go to(2010,2020,2030,2030,2030,2040)numele
        go to 9000
 
 2000   continue

***   ONE item, velocity structure or ICF switch

 2010   if (linele(1)(1:3).ne.'vel' .and. linele(1)(1:3).ne.'icf')
     1  go to 9001
        if (lgt(1).ne.4) go to 9001
        if (index(linele(1),'+').eq.0 .and. index(linele(1),'-').eq.
     1  0) go to 9001
        if (lgt(1).ne.4) go to 9001
        if (linele(1)(1:3).eq.'vel') then
           if (velset) go to 9002 
           if (linele(1)(4:4).eq.'+') then
              calcvel=.true.
              use_def_rvcor=.false.
           end if
           if (linele(1)(4:4).eq.'-') then
              calcvel=.false.
           end if
           velset=.true.
        else
           if (icfset) go to 9003
           if (linele(1)(4:4).eq.'+') then 
              calcicf=.true.
              use_def_icf=.false.
           end if 
           if (linele(1)(4:4).eq.'-') then 
              calcicf=.false.
           end if
           icfset=.true.
        end if         
        go to 1000

***   TWO items, file specifiers, or nebular/instrument attributes


 2020   continue
        if (lgt(1).ne.1) go to 9001
        ind=0
        do ii=1,8
           if (linele(1)(1:lgt(1)).eq.att(ii)) ind=ii   
        end do
        if (ind.eq.0) go to 9001
        if (attset(ind)) go to 9004                            
        if (ind.ge.6) then            
           read(linele(2)(1:lgt(2)),*,err=9005) realarg(ind-5)
        else
           attarg(ind)=linele(2)
           attarglen(ind)=lgt(2)
        end if
        attset(ind)=.true.
        go to 1000

***   THREE, FOUR, or FIVE items, ions to depelete and amount

 2030   if (linele(1)(1:lgt(1)).ne.'deplete') go to 9001
        if (lgt(1).ne.7) go to 9001
        read(linele(3)(1:lgt(3)),*,err=9005) temparg
        write(tempele,'(a2)') linele(2)(1:lgt(2))
        do ii=1,30
           if (tempele.eq.slab(ii)) then
              element=ii
              go to 2035
           end if
        end do
        go to 9010
 2035   continue
        do ii=1,kk
           if (element.eq.depele(ii)) go to 9011
        end do
 2037   continue
        if (numele.ge.4 .and. numele.le.5) then
           read(linele(4)(1:lgt(4)),*,err=9005) tempint
           if (tempint.gt.(element+1)) go to 9012
           if (numele.eq.5) then
              read(linele(5)(1:lgt(5)),*,err=9005) tempint
              if (tempint.gt.(element+1)) go to 9012
           end if
           kk=kk+1
           read(linele(4)(1:lgt(4)),*) depelelow(kk)
           if (numele.eq.5) then
              read(linele(5)(1:lgt(5)),*) depelehi(kk)
           else
              depelehi(kk)=depelelow(kk)
           end if
        else
           kk=kk+1
           depelelow(kk)=1
           depelehi(kk)=element+1
        end if
        if (depelelow(kk).gt.depelehi(kk)) then
           tempint=depelelow(kk)
           depelelow(kk)=depelehi(kk)
           depelehi(kk)=tempint
        end if
        depele(kk)=element
        depeleamount(kk)=temparg
        if (depeleamount(kk).lt.0) depeleamount(kk)=
     1  1./abs(depeleamount(kk))
        go to 1000
      
***   SIX elements only, manually specified ICF or velocity struct.

 2040   continue
        if (linele(1)(1:lgt(1)).ne.'vel' .and. linele(1)(1:lgt(1))
     1  .ne.'icf') go to 9001
        do ii=1,5
           read(linele(ii+1)(1:lgt(ii+1)),*,err=9005) temparg
        end do
        if (linele(1)(1:lgt(1)).eq.'vel') then
           if (velset) go to 9002
           do ii=1,5
              read(linele(ii+1)(1:lgt(ii+1)),*,err=9005) vel(ii)
           end do
           calcvel=.false.
           velset=.true.
           manvel=.true.
           use_def_rvcor=.false.
        else
           if (icfset) go to 9003
           sum=0.
           do ii=1,5
              read(linele(ii+1)(1:lgt(ii+1)),*,err=9005) icf(ii)
              sum=sum+icf(ii)
           end do
c           if (sum.ne.1) go to 9006   
           calcicf=.false.
           icfset=.true.
           manicf=.true.
           use_def_icf=.false.
        end if
        go to 1000

***   Now checking for the validity of the commands issued.
***   Also setting default values for some parameters.


 3000   continue
        jj=jj+1
        if (jj.gt.8) go to 4000
        if (.not. attset(jj)) then 
           if (jj.eq.1) then
              go to 9007
           end if
           if (jj.eq.2) then
              use_match=.false.
              if (calcvel .or. calcicf) go to 9008
              attarg(2)='none'
              attarglen(2)=4
           end if
           if (jj.eq.3) then
              attarg(3)=attarg(1)(1:attarglen(1))//'.out'
              attarglen(3)=attarglen(1)+4
           end if
           if (jj.eq.4) then
              attarg(4)=attarg(1)(1:attarglen(1))//'.dat'
              attarglen(4)=attarglen(1)+4
           end if
           if (jj.eq.5) then
              attarg(5)='abun.dat'
              attarglen(5)=8
           end if
           if (jj.gt.5) go to 9009
        end if
        go to 3000
      
***   Here is where we note fatal/non-fatal errors

        
 4000   continue
        write(*,*) 'Error Status: '
        write(*,*) '--------------------------------'
        if (.not. fatal .and. .not. non_fatal) then
           write(*,*)'Command List Had NO Fatal/Non-Fatal Errors'
        end if
        if (fatal) then
           write(*,*)'Command List Had The Following Fatal Errors:'
           write(*,*)
           do ii=1,fatalerrs
              write(*,'(i3,a2,$)') ii,') ' 
              if (ferrlinenum(ii).ne.99) then
                 write(*,'(a6,$)') 'Line: '
                 write(*,'(i3,a3,a)') ferrlinenum(ii),'-> ',ferrline(ii)
              end if  
              write(*,*) ferrmessage(ii)
              write(*,*)
           end do
           write(*,*)'These Must Be Corrected Before Program Will Run'
           write(*,*)
           if (.not. non_fatal) go to 9904
        end if
        if (non_fatal) then
           write(*,*)'Command List Had The Following Non-Fatal'//
     1     ' Errors:'
           write(*,*)
           do ii=1,nfatalerrs
              write(*,'(i3,a2,$)') ii,') '
              if (nferrlinenum(ii).ne.99) then
                write(*,'(a6,$)') 'Line: '
                write(*,'(i3,2x,a)') nferrlinenum(ii),nferrline(ii)
              end if
              write(*,*) nferrmessage(ii)
              write(*,*)
           end do
           if (fatal) go to 9904
           write(*,*)
           write(*,*)'You Wish To Check The List Before Proceeding'
           write(*,'(a,$)')'Continue Without Changes? (y/n) '
           read(*,*) res
           if (res.ne.'y' .and. res.ne.'Y') go to 9905
        end if

***   Here is where we sum up the Commands to EMILI

 5000   continue
        write(*,*)
        write(*,*) 'Arguements submitted to EMILI: '
        write(*,*) '--------------------------------'
        do ii=1,8
           if (ii.lt.6) then
              write(*,'(1x,a,a2,a)') attlabel(ii)(1:attlabellen(ii)),
     1        ': ',attarg(ii)(1:attarglen(ii))
           else
              write(*,'(1x,a,a2,f7.0)') attlabel(ii)(1:attlabellen(ii)),
     1        ': ',realarg(ii-5)
            endif
        write(*,*)
        end do

        write(*,'(1x,a,$)') 'Velocity Structure: ' 
        if (calcvel) write (*,*) 'To be calculated'
        if (.not. calcvel .and. use_def_rvcor) write(*,*) 'Turned off'
        if (.not. calcvel .and. .not. use_def_rvcor) then
           write(*,*) 'Manually specified'
           write(*,'(5f7.2)')(vel(iii), iii=1,5,1)
        end if
        write(*,*)
        write(*,'(1x,a,$)') 'ICF Values: '
        if (calcicf) write(*,*) 'To be calculated'
        if (.not. calcicf .and. use_def_icf) then
           write(*,*) 'Default values'
        end if
        if (.not. calcicf .and. .not. use_def_icf) then
           write(*,*) 'Manually specified'
           write(*,'(5f7.3)')(icf(iii), iii=1,5,1)
        end if
        write(*,*)

        write(*,'(1x,a,$)') 'Elements Depleted: '
        if (kk.eq.0) write(*,*) 'None'
        if (kk.ne.0) then
           write(*,*)
           do ii=1,kk
             if (depelehi(ii).le.30) then
                write(*,'(1x,a2,2x,a6,a2,a6,2x,f9.2)') slab(depele(ii)),
     1          spcs(depelelow(ii)),'- ',spcs(depelehi(ii)),
     2          depeleamount(ii)
             else
                write(*,'(1x,a2,2x,a6,a2,a6,2x,f9.2)') slab(depele(ii)),
     1          spcs(depelelow(ii)),'- ',spcs(30),
     2          depeleamount(ii)
             end if
           end do
        end if
        depelenumber=kk

        write(*,*) '--------------------------------'
        write(*,*)
        write(*,*) 'Press RETURN to begin set-up'
        read(*,*)


***   Now Opening All Files

6000    continue

*     Ionization Energy Table

        open(21,file='ion.dat',status='OLD',err=9906)

*     Abundance Table

        open(22,file=attarg(5),status='OLD',err=9907)
 
*     The Level List

        open(23,file='slev_list',status='OLD',err=9908)

*     The Line List

        open(24,file='sline_list',status='OLD',err=9909)

*     The INPUT Line List

        open(25,file=attarg(1),status='OLD',err=9910)

*     The MATCHED Line List

        if (use_match)
     1    open(26,file=attarg(2),status='OLD',err=9911)

*     The EMILI Summary List

        open(27,file=attarg(4),err=9913)

*     The EMILI Full Output List

        open(28,file=attarg(3),err=9912)

*     Opening The Special Values List

        open(29,file='multi.dat',status='OLD',err=9914)
        
        write(*,*) 'FILES: Opened For Reading'

***   Assigning Nebular/Instrumental Attributes

        t_e=realarg(1)*1.0E-4
        n_e=realarg(2)
        instres=realarg(3)
        width=realarg(3)

        write(*,*) 'ATTRIBUTES: Assigned'

***   Writing Info To EMILI Output File

        write(28,*) 'EMILI Output File'
        write(28,*) '------------------------------------'
        do ii=1,8
           if (ii.le.5) then
              write(28,'(1x,a,a2,a)') attlabel(ii)(1:attlabellen(ii)),
     1        ': ',attarg(ii)(1:attarglen(ii))
           else
             write(28,'(1x,a,a2,f7.0)') attlabel(ii)(1:attlabellen(ii)),
     1       ': ',realarg(ii-5)
           end if
        end do
        write(28,*)

        go to 9999

***   Generating Internal Error Messages

*     Too Many/Long Arguements in Command Line or Comand Line
*     Too Long

 9000   fatalerrs=fatalerrs+1
        ferrlinenum(fatalerrs)=line
        ferrline(fatalerrs)=commandline
        ferrmessage(fatalerrs)='Too Many/Long Arguements in '//
     1    'Command Line or Command Line Too Long'
        fatal=.true.
        go to 1000

*     Non Standard Format of Command Line

 9001   fatalerrs=fatalerrs+1
        ferrlinenum(fatalerrs)=line
        ferrline(fatalerrs)=commandline
        ferrmessage(fatalerrs)='Something Is Wrong With This Line'
        fatal=.true.
        go to 1000

*     Velocity Structure Command Specified Multiple Times

 9002   nfatalerrs=nfatalerrs+1
        nferrlinenum(nfatalerrs)=line
        nferrline(nfatalerrs)=commandline
        nferrmessage(nfatalerrs)='Velocity Structure Command Specified
     1  Multiple Times'
        non_fatal=.true.
        go to 1000

*     ICF Factors Specified Multiple Times

 9003   nfatalerrs=nfatalerrs+1
        nferrlinenum(nfatalerrs)=line
        nferrline(nfatalerrs)=commandline
        nferrmessage(nfatalerrs)='ICF Factors Specified Multiple Times'
        non_fatal=.true.
        go to 1000

*     A Attribute or File Name Has Been Specified Multiple Times

 9004   nfatalerrs=nfatalerrs+1
        nferrlinenum(nfatalerrs)=line
        nferrline(nfatalerrs)=commandline
        nferrmessage(nfatalerrs)=attlabel(ind)(1:attlabellen(ind))
     1  //' Specified Multiple Times'
        non_fatal=.true.
        go to 1000

*     Numeric Argument Expected in Command Line

 9005   fatalerrs=fatalerrs+1
        ferrlinenum(fatalerrs)=line
        ferrline(fatalerrs)=commandline
        ferrmessage(fatalerrs)='Numeric Argument Expected in'//
     1  ' Command Line'
        fatal=.true.
        go to 1000

*     ICF Values Don't Sum to 1

 9006   nfatalerrs=nfatalerrs+1
        nferrlinenum(nfatalerrs)=line
        nferrline(nfatalerrs)=commandline
        nferrmessage(nfatalerrs)='ICF Values Do Not Sum to 1'
        non_fatal=.true.
        go to 1000

*     Input Line List Not Specified

 9007   fatalerrs=fatalerrs+1
        ferrlinenum(fatalerrs)=99
        ferrline(fatalerrs)='Z'
        ferrmessage(fatalerrs)='Input Line List Not Specified'
        fatal=.true.
        go to 3000

*     Matched List Missing: Neccessary to Calculate ICF or Vel. 

 9008   fatalerrs=fatalerrs+1
        ferrlinenum(fatalerrs)=99
        ferrline(fatalerrs)='Z'
        ferrmessage(fatalerrs)='Matched List Missing: Neccessary to'// 
     1  ' Calculate ICF or Vel.' 
        fatal=.true.
        go to 3000

*     Numeric Values For Temp, Density, or Inst. Res. Not Specified

 9009   fatalerrs=fatalerrs+1
        ferrlinenum(fatalerrs)=99
        ferrline(fatalerrs)='Z'
        ferrmessage(fatalerrs)=attlabel(jj)(1:attlabellen(jj))//
     1  ' Not Specified'
        fatal=.true.
        go to 3000

*     Element Not Recognized
      
 9010   nfatalerrs=nfatalerrs+1
        nferrlinenum(nfatalerrs)=line
        nferrline(nfatalerrs)=commandline
        ferrmessage(nfatalerrs)='Element Not Recognized'
        non_fatal=.true.
        go to 1000

*     Element Depeleted More Than Once

 9011   nfatalerrs=nfatalerrs+1
        nferrlinenum(nfatalerrs)=line
        nferrline(nfatalerrs)=commandline
        nferrmessage(nfatalerrs)='Element Depleted More Than Once'
        non_fatal=.true.
        go to 2037

*     Ionization Stage Out of Bounds

 9012   fatalerrs=nfatalerrs+1
        ferrlinenum(nfatalerrs)=line
        ferrline(nfatalerrs)=commandline
        ferrmessage(nfatalerrs)='Ionization Stage Out of Bounds'
        fatal=.true.
        go to 1000

***   Program Faults

 9901   write(*,*) 'WARNING: Openall.f err=9901'
        write(*,*) 'General Fault Reading Command List'
        stop

 9902   write(*,*) 'WARNING: Openall.f err=9902'
        write(*,*) 'Cannot Open Command List File'
        stop

 9903   write(*,*) 'WARNING: Openall.f err=9903'
        write(*,*) 'Cannot Read a Particular Line From Command List'
        stop

 9904   write(*,*) 'WARNING: Openall.f err=9904'
        write(*,*) 'Fatal Error Stop'
        stop

 9905   write(*,*) 'WARNING: Openall.f err=9905'
        write(*,*) 'User Chose to Terminate Program'
        stop

 9906   write(*,*) 'WARNING: Openall.f err=9906'
        write(*,*) 'Cannot Open Ionization Energy Data'
        stop

 9907   write(*,*) 'WARNING: Openall.f err=9907'
        write(*,*) 'Cannot Open Abundance Table'
        stop
 
 9908   write(*,*) 'WARNING: Openall.f err=9908'
        write(*,*) 'Cannot Open Master Level List'
        stop
      
 9909   write(*,*) 'WARNING: Openall.f err=9909'
        write(*,*) 'Cannot Open Master Transition List'
        stop

 9910   write(*,*) 'WARNING: Openall.f err=9910'
        write(*,*) 'Cannot Open Input Line List'
        stop

 9911   write(*,*) 'WARNING: Openall.f err=9911'
        write(*,*) 'Cannot Open Matched Line List'
        stop

 9912   write(*,*) 'WARNING: Openall.f err=9912'
        write(*,*) 'Cannot Open EMILI Output List'
        stop

 9913   write(*,*) 'WARNING: Openall.f err=9913'
        write(*,*) 'Cannot Open the EMILI Short Output List'
        stop

 9914   write(*,*) 'WARNING: Openall.f err=9914'
        write(*,'(a,$)') 'Cannot Open the EMILI Special Values List: '
        write(*,*) 'multi.dat'
        stop

****  END SUBROUTINE

 9999   continue
        return
        end

*******************************************************************
     
      SUBROUTINE split(instring,outele,lenele,numele,badline)

******  Splits a input string "instring" into elements which
******  are denoted by whitespace in between.  Returns these
******  elements in an array called "outele", recording
******  the number of characters in each element as "lenele",
******  and the number of such elements as numele
*
*
*       INPUT:
*       instring - the line in which to split - CHARACTER*60       
*
*       OUTPUT:
*       outele   - the elements on the line, split by whitespace
*                  ARRAY(6) of CHARACTER*10
*       lenele   - the lenght of each such element
*                  ARRAY(6) of INTEGER
*       numele   - the number of elements (max of 10) - INTEGER
*       badline  - if more than 6 elements, this is a bad line
*                  trigger an error message in OPENALL.F
*                  - LOGICAL

      IMPLICIT NONE

*     INPUT VARIABLES

      CHARACTER*60  instring
    
*     OUTPUT VARIABLES

      CHARACTER*20  outele(6)
      INTEGER       lenele(6),numele
      LOGICAL       badline

*     OTHER VARIABLES
      INTEGER       ii,i1,i2,jj
      CHARACTER*20  woutele(6)
      INTEGER       wlenele(6)
    

***   Set Initial Values

      ii=1
      jj=0
      badline=.false.

***   If Too Long, trip a flag in OPENALL 

      if (len(instring).gt.60) go to 4000

***   Commence Looping Through Line List

 1000 continue
      do while (instring(ii:ii).eq.' ')
         ii=ii+1
         if (ii.gt.60) go to 3000
      end do
      i1=ii
      do while (instring(ii:ii).ne.' ')
         ii=ii+1
         if (ii.gt.60) then
            i2=60
            go to 2000
         end if
      end do
      i2=ii-1
 2000 jj=jj+1
      if (jj.gt.6) go to 4000
      if (i2-i1+1 .gt. 20) go to 4000
      woutele(jj)=instring(i1:i2)
      wlenele(jj)=i2-i1+1
      if (ii.lt.60) go to 1000
      
 3000 continue
      do ii=1,jj
         outele(ii)=woutele(ii)
         lenele(ii)=wlenele(ii)
      end do
      numele=jj
       
      go to 5000

 4000 badline=.true.
 5000 continue
      end






         
      
      

      

     
      
             
               
               

      
            
      
      
      
