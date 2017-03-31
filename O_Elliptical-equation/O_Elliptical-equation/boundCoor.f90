module boundCoor   !determine boundary coordinates  
    use controlDara
    implicit none
    real(kind=8)::x_in(imax),y_in(imax)    !the x,y coordinates of inside boundary
    real(kind=8)::x_out(imax),y_out(imax)  !the x,y coordinates of outside boundary
contains

subroutine solverBoundary
    call inside
    call outside

end subroutine

subroutine inside               !determine the x,y coordinates of inside boundary
    implicit none
    integer,parameter::fileid=10  !the file id
    character(len=80)::filename   !file name
    integer::error                !the return value of verifying wheather the file is readed correctly 
    logical::alive                !the return value of inquiring wheather the file  exsits
    integer::i   
    
    !write(*,*)  "Filename:"         
    !read(*,"(A80)") filename   

    filename="naca0012.dat" 
    
    inquire(file=filename,exist=alive)   !verify wheather the file is readed correctly 
    
    if(.not.alive) then                  !if not, prit error information and stop the program
        write(*,*) trim(filename),"Does't exist!"   
        stop                             
    end if
    
    open(fileid,file=filename)
    do i=1,imax
        read(fileid,*,iostat=error) x_in(i),y_in(i)  
        write(*,*)  x_in(i),y_in(i)  
        if(error/=0) then
           write(*,*) "Read Error!"
           stop
        end if 
        
        x_in(i)= x_in(i)-Center
    end do
    
    close(fileid)
   
end subroutine

subroutine outside          !determine the x,y coordinates of outside boundary
    implicit none          
    integer::i
    !distribute the x,y coordinates by spliting angles averagely
    !alternative way: split arc length averagely
    do i=1,imax  
        x_out(i)= R_out*cos(-2*pi/(imax-1)*(i-1))
        y_out(i)= R_out*sin(-2*pi/(imax-1)*(i-1))
    end do
    
end subroutine

end module



    
         
        