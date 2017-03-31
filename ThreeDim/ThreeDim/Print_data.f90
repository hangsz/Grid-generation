module Print_data
    use Solution
    implicit none
    
contains
subroutine output
    implicit none
    integer::i,j,k
    integer,parameter::fileid=7
    
    !**************************************
    !输出内边界x，y坐标
    !open(fileid,file="x_inner.txt")
    !write(fileid,"(I3,3X,F11.8,3X,F11.8)") (i,x_inner(i),y_inner(i),i=1,imax) 
    
    !输出外边界x，y坐标
    !open(fileid,file="x_outer.txt")
    !write(fileid,"(I3,3X,F11.8,3X,F11.8)")  (i,x_outer(i),y_outer(i),i=1,imax)
    !write(fileid,"(F11.8,3X,F11.8,3X,F11.8,3X,F11.8)")  (x_outer(i,jmax),y_outer(i,jmax),x_outer(imax-i+1,jmax),y_outer(imax-i+1,jmax),i=1,(imax+1)/2)
    !write(*,"(F5.2,3X,F5.2)")(x_outer(i,jmax),y_outer(i,jmax),i=1,imax)
    !******************************************
    
    !输出物理域所有点的x,y初始化值和p,q的初始化值
    
    !输出x,y,z坐标
    open(fileid,file="X_Y_Z.dat")
    write(fileid,*)  'VARIABLES="X","Y","Z"'
    write(fileid,"(2X,A7,I2,A3,I2,A3,I2)")  "zone I=",imax,",J=",jmax,",K=",kmax-1
    write(fileid,"(2X,F11.8,2X,F11.8,2X,F11.8)") (((x(i,j,k),y(i,j,k),z(i,j,k),i=1,imax),j=1,jmax),k=1,kmax-1)
   ! write(fileid,"(F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8)") ((x(i,j),y(i,j),i=1,imax),j=1,jmax)
    close(fileid) 
    
     open(fileid,file="Lk.dat")
     write(fileid,"(2X,F11.8)") (lk_init(k),k=2,kmax-1)
   ! write(fileid,"(F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8,2X,F11.8)") ((x(i,j),y(i,j),i=1,imax),j=1,jmax)
    close(fileid) 
    
    !输出源项值
    !open(fileid,file="p_q.txt")
    !write(fileid,"(I3,2X,I3,2X,F12.8,2X,F12.8)") ((i,j,p(i,j),q(i,j),i=1,imax),j=1,jmax)
    !close(fileid) 
    !**************************************************************************
    
    !输出队xi,eta的导数
    open(fileid,file="r_xi.dat")
    write(fileid,"(I3,2X,I3,2X,I3,2X,F12.8,2X,F12.8,2X,F12.8)") (((i,j,k,x_xi(i,j,k),y_xi(i,j,k),z_eta(i,j,k),i=1,imax),j=1,jmax),k=2,kmax-1)
    close(fileid)
    
    !输出队eta的导数
    open(fileid,file="r_eta.dat")
    write(fileid,"(I3,2X,I3,2X,I3,2X,F12.8,2X,F12.8,2X,F12.8)") (((i,j,k,x_eta(i,j,k),y_eta(i,j,k),z_eta(i,j,k),i=1,imax),j=1,jmax),k=2,kmax-1)

    close(fileid)
   
    
    !输出对zeta的导数
    open(fileid,file="r_zeta.dat")
    write(fileid,"(I3,2X,I3,2X,I3,2X,F12.8,2X,F12.8,2X,F12.8)") (((i,j,k,x_zeta(i,j,k),y_zeta(i,j,k),z_zeta(i,j,k),i=1,imax),j=1,jmax),k=2,kmax-1)
    close(fileid) 
      
    !输出对xieta的导数
    ! open(fileid,file="x_xieta.txt")
    ! write(fileid,"(I3,2X,I3,2X,F12.8,2X,F12.8,,2X,F12.8)") (((i,j,x_xieta(i,j,k),y_xieta(i,j,k),z_xieta(i,j,k),i=1,imax),j=1,jmax),k=2,kmax-1)
    !close(fileid) 
      
    !输出几个系数
    open(fileid,file="Coeff.txt")
    write(fileid,"(I3,2X,I3,2X,I3,2X,F12.8,2X,F12.8,2X,F12.8,2X,F12.8)") (((i,j,k,alpha1(i,j,k),alpha2(i,j,k),alpha3(i,j,k),beta12(i,j,k),i=1,imax),j=1,jmax),k=2,kmax-1)
    close(fileid) 
    
    !******************************************
    
end subroutine

end module
