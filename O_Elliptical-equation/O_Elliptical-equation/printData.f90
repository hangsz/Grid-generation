module printData    !print data
    use ellipse
    implicit none
    
contains

subroutine output    !output needed information
    implicit none
    integer::i,j
    integer,parameter::fileid=7              !the  file  id 
    
    if(judgeSource=='y') then
        open(fileid,file="NACA0012-Y.dat")
        write(fileid,*)  'VARIABLES="X","Y"'
        write(fileid,99)  "zone I=",imax,", J=",jmax
        write(fileid,100) ((x(i,j),y(i,j),i=1,imax),j=1,jmax)
        close(fileid) 
    else
        open(fileid,file="NACA0012-N.dat")
        write(fileid,*)  'VARIABLES="X","Y"'
        write(fileid,99)  "zone I=",imax,", J=",jmax
        write(fileid,100) ((x(i,j),y(i,j),i=1,imax),j=1,jmax)
        close(fileid) 
    end if
    
99  format(A7,I3,1X,A4,I3)
100 format(2X,F11.8,2X,F11.8)  
      
end subroutine

end module
