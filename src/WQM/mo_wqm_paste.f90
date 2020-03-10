!> \file mo_wqm_paste.f90

!> \brief paste and append variables for water quality model

!> \details 

!> \authors Xiaoqiang Yang --Modified from corresponding mHM routines
!> \date Jun 2016

MODULE mo_wqm_paste


  USE mo_kind, ONLY: i4, sp, dp


  IMPLICIT NONE


  PUBLIC :: paste_conc   ! 
  PUBLIC :: append_3d

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         paste_conc

  !     PURPOSE
  !>        \brief 

  !>        \details 

  !     CALLING SEQUENCE
  !        

  !     INTENT(IN)
  !>        \param[in] 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !        
  !        

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoiang Yang
  !>        \date Jun 2016

  subroutine paste_conc(mat1, mat2, fill_value)

    use mo_wqm_global_variables,  only: maxcols  !maximum number of observed data type
    implicit none
    real(dp), dimension(:,:,:), allocatable, intent(inout) :: mat1
    real(dp), dimension(:,:),                intent(in)    :: mat2
    real(dp), optional,                      intent(in)    :: fill_value

    !local
    integer(i4)             :: t1,t2,t3   !three dims for target matrix (mat1)
    integer(i4)             :: s1,s2,s3   !s1 and s3 two dims for source matrix (mat2)
    real(dp), dimension(:,:,:), allocatable  :: tmp
	
    s1 = size(mat2,1)  !rows
    s2 = 1_i4          !one station per call 
    s3 = size(mat2,2)  !columns, should equals "maxcols"
    if (s3 .ne. maxcols) then
       print*, '***ERROR: paste_conc:: 2nd dimension of source matrix IS NOT CORRECT!'
       STOP
    end if
    if (allocated(mat1)) then
       t1 = size(mat1,1)
       t2 = size(mat1,2)
       t3 = size(mat1,3)
       if ((t1 .ne. s1) .and. .not.present(fill_value)) then
          print*, '***ERROR: paste_conc:: rows of target matrix and source matrix are unequal : (',t1,')  and (',s1,')'
          STOP
       end if       
	   !save mat1 to temporal variable
       allocate(tmp(t1,t2,t3))
       tmp = mat1
       deallocate(mat1)
       
       if (t1 .eq. s1) then
          allocate(mat1(t1,t2+s2,t3))
          mat1(1:t1,1:t2,1:t3) = tmp
          mat1(1:s1,t2+s2,1:t3) = mat2(1:s1,1:s3)
       end if
       if (t1 .gt. s1) then
          allocate(mat1(t1,t2+s2,t3))
          mat1(1:t1,1:t2,1:t3) = tmp
          mat1(1:s1,t2+s2,1:t3) = mat2(1:s1,1:s3)
          mat1(s1+1:t1,t2+s2,1:t3) = fill_value
       end if
       if (t1 .lt. s1) then
          allocate(mat1(s1,t2+s2,t3))
          mat1(1:t1,1:t2,1:t3) = tmp
          mat1(t1+1:s1,1:t2,1:t3) = fill_value
          mat1(1:s1,t2+s2,1:t3) = mat2(1:s1,1:s3)
       end if
   
    else
       t2 = 0_i4   !initial number of gauges
       t1 = s1     !number of obervations
       t3 = s3     !number of measured data types
       allocate(mat1(s1, s2, s3))
       mat1(1:s1, t2+ s2, 1:s3) = mat2(1:s1, 1:s3)
    end if


  end subroutine paste_conc
!
!------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         append_3d

  !     PURPOSE
  !>        \brief 

  !>        \details 

  !     CALLING SEQUENCE
  !        

  !     INTENT(IN)
  !>        \param[in] 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !        
  !        

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  subroutine append_3d(mat1,mat2, fill_value)
  
  implicit none
  real(dp), dimension(:,:,:), allocatable, intent(inout)  :: mat1   !target matrix
  real(dp), dimension(:,:,:), intent(in)     :: mat2                !source matrix
  real(dp), optional,         intent(in)     :: fill_value          !optional nodata value

  !local variables
  integer(i4)       :: t1,t2,t3                    !dimensions of mat1
  integer(i4)       :: s1,s2,s3                    !dimensions of mat2
  real(dp), dimension(:,:,:), allocatable   :: tmp

  s1 = size(mat2,1) 
  s2 = size(mat2,2)
  s3 = size(mat2,3) 

  if (allocated(mat1)) then
     t1 = size(mat1,1)
     t2 = size(mat1,2)
     t3 = size(mat1,3)
     if ((t2 .ne. s2) .and. .not. present(fill_value)) then
        print*, '***ERROR: append_3d:: columns of target and source matrixs are unequal.'
        STOP
     end if
     !save mat1
     allocate(tmp(t1,t2,t3))
     tmp = mat1
     deallocate(mat1)

     if (t2 .eq. s2) then
        allocate(mat1(t1+s1, t2, t3) )
        mat1(1:t1,:,:) = tmp(1:t1,:,:)
        mat1(t1+1: t1+s1,:,:) = mat2(:,:,:)
     end if
     if (t2 .gt. s2) then
        allocate(mat1(t1+s1, t2, t3) )
        mat1(1:t1,1:t2,1:t3) = tmp
        mat1(t1+1: t1+s1, 1:s2,:) = mat2(:,:,:)
        mat1(t1+1: t1+s1, s2+1:t2,:) = fill_value
     end if
     if (t2 .lt. s2) then
        allocate(mat1(t1+s1,s2,t3))
        mat1(1:t1,1:t2,1:t3) = tmp
        mat1(1:t1,t2+1:s2,:) = fill_value
        mat1(t1+1:t1+s1,1:s2,:) = mat2(:,:,:)
     end if
		
  else
     t1 = 0_i4
     allocate (mat1(s1,s2,s3))
     mat1= mat2

  end if


  end subroutine append_3d


END MODULE mo_wqm_paste
