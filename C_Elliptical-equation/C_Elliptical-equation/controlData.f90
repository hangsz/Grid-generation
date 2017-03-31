module controlData       !the control parameters 
    !------------------------------------------------------------------------
    !constant parameters
    real(8),parameter::pi=3.14159265
    
    !-------------------------------------------------------------------
    !parameters used to control the grid
    integer,parameter::imax=183,jmax=60 , &             !maximum  grid points in i  and j  directions    
                       aerofoil_num=103 ,&                !the number of grid points on the aerofoil
                       trail_num = 5                    !not include the point on the airfoil
    real(8),parameter::Center=1.008400 , &                !the x coordinate of aerofoil
                       x_left=-10.0, x_right=5.0 , &     !the domain of x coordinate
                       y_bottom=-10.0, y_top=10.0 , &     !the domain of y coordinate
                       r=10.0                            !the radius of the circular bead on outside boundary 
    
    !-------------------------------------------------------------------
    !parameters about iteration 
    integer,parameter::iter_out=10 , &                  !the times of outer iteration
                       itermax=3000                     !the times of inner iteration
                     
    real(8),parameter::eps=1.0E-10 , &                   !the convergence precision of x,y coordinates 
                       omg=1.0                          !the relaxation  factor(0-2)
    
    
    !------------------------------------------------------------------
    !parameters used to control the source term 
    !real(8),parameter::a=0.45,b=0.45,c=0.45,d=0.45  
    !real(8),parameter::a=0.7,b=0.7,c=0.7,d=0.7  
    real(8),parameter::a=0.3,b=0.3,c=0.3,d=0.3 , &        !the factors of  exponential interpolation
                       sigma=0.3  !the correction factor of source term (0-0.5)
    real(8),parameter::betaSource=0.8
    
    !the  four parameters about the restrain conditions of boundaries 
    real(8),parameter::dr_i=1.0E-5,thetar_i=pi/2 , &                !inside boundary        
                       dr_o=0.1,thetar_o=pi/2 , &                !outside boundary
                       eps1=1.0E-6,eps2=1.0E-4        !the convergence precisions of the four parameters
    !-------------------------------------------------------------------
end module