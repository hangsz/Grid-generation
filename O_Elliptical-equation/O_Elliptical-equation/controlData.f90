module controlDara          !the control parameters 
    !*******************************************************************
    !constant parameters
    real(kind=8),parameter::pi=3.14159265
    !------------------------------------------------------------
    !parameters used to control the grid
    integer,parameter::imax=103,jmax=100        !maximum  grid points in i  and j  directions    
    real(kind=8),parameter::Center=0.50        !the x coordinate of aerofoil
    real(kind=8),parameter::R_out=6.0        !the radius of the computianal  domain
    
    !-------------------------------------------------------------
    !parameters about iteration 
    integer::iterOut=100                 !the times of outer iteration
    integer::iterMax=3000                !the times of inner iteration
    real(kind=8),parameter::eps=1.0E-7   !the convergence precision of x,y coordinates 
    
    real(kind=8),parameter::omg=1.0      !the relaxation  factor(0-2)
    
    
    !-------------------------------------------------------------
    !parameters used to control the source term
    real(kind=8),parameter::a=0.3,b=0.3,c=0.3,d=0.3        !the factors of  exponential interpolation
    !real(kind=8),parameter::a=0.45,b=0.45,c=0.45,d=0.45  
    !real(kind=8),parameter::a=0.7,b=0.7,c=0.7,d=0.7  
    real(8),parameter::betaSource=0.8
    real(kind=8),parameter::sigma=0.3                      !the correction factor of source term (0-0.5)
    !the  four parameters about the restrain conditions of boundaries 
    real(kind=8)::dr_i=1.0E-5,thetar_i=pi/2.0                !inside boundary        
    real(kind=8)::dr_o=0.8,thetar_o=pi/2.0                !outside boundary
    real(kind=8),parameter::eps1=1.0E-7,eps2=1.0E-4        !the convergence precisions of the four parameters
    !--------------------------------------------------------------
       
end module