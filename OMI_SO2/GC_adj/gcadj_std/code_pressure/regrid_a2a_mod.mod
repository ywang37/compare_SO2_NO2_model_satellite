  �!  7   k820309    `          17.0        ֺ@^                                                                                                           
       regrid_a2a_mod.F90 REGRID_A2A_MOD              DO_REGRID_A2A DO_REGRID_DKH MAP_A2A #         @     @                                                #IM    #JM    #LON1    #Q1    #IN    #LON2    #Q2              
                                                      
                                                     
                                                     
    p           5 � p        r    n                                       1     5 � p        r    n                                      1                                    
                                                    
      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                
                                                     
                                                    
     p           5 � p        r    n                                       1     5 � p        r    n                                      1                                    D                                                   
 !      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                      #         @     @                            	                 	   #IM 
   #JM    #SIN1    #Q1    #JN    #SIN2    #Q2    #IG    #IV              
  @                              
                     
                                                     
                                                    
    p            5 � p        r    n                                           15 � p        r          5 � p        r    n                                      15 � p        r                                    
                                                    
      p        5 � p        r 
   p          5 � p        r 
     5 � p        r        5 � p        r 
     5 � p        r                                
                                                     
                                                    
    p            5 � p        r    n                                           15 � p        r          5 � p        r    n                                      15 � p        r                                    D                                                   
       p        5 � p        r 
   p          5 � p        r 
     5 � p        r        5 � p        r 
     5 � p        r                                
                                                      
                                             #         @     @                                                #IM    #JM    #FILENAME    #LON_EDGES    #LAT_SINES              
                                                       
                                                       
  @                                                  1          D @                                                  
 %    p           5 � p        r    n                                       1     5 � p        r    n                                      1                                    D @                                                  
 &    p           5 � p        r    n                                       1     5 � p        r    n                                      1                           #         @                                                       #FILENAME    #IM    #JM    #INGRID    #OUTGRID    #IS_MASS    #NETCDF               
  @                                                  1           
  @                                                   
  @                                                  
                                                     
      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                D @                                  03             
     p �         p �         p [           p �         p [                                   
                                                       
 @                                           #         @                                   !                    #FILENAME "   #IM #   #JM $   #INGRID %   #OUTGRID &   #IS_MASS '   #NETCDF (             
  @                              "                    1           
  @                              #                     
  @                              $                    
                                 %                    
 	     p        5 � p        r #   p          5 � p        r #     5 � p        r $       5 � p        r #     5 � p        r $                               D                               &     03             
 
    p �         p �         p [           p �         p [                                   
                                  '                     
 @                               (           #         @                                  )                    #IM *   #JM +   #LON1 ,   #SIN1 -   #Q1 .   #IN /   #JN 0   #LON2 1   #SIN2 2   #Q2 3   #IG 4   #IV 5             
  @                              *                     
  @                              +                    
  @                              ,                    
    p           5 � p        r *   n                                       1     5 � p        r *   n                                      1                                    
  @                              -                    
    p           5 � p        r +   n                                       1     5 � p        r +   n                                      1                                    
  @                             .                    
      p        5 � p        r *   p          5 � p        r *     5 � p        r +       5 � p        r *     5 � p        r +                               
  @                              /                     
  @                               0                    
  @                              1                    
    p           5 � p        r /   n                                       1     5 � p        r /   n                                      1                                    
  @                              2                    
    p           5 � p        r 0   n                                       1     5 � p        r 0   n                                      1                                    D @                             3                    
       p        5 � p        r /   p          5 � p        r /     5 � p        r 0       5 � p        r /     5 � p        r 0                               
  @                              4                     
  @                               5              �   *      fn#fn $   �   4   b   uapp(REGRID_A2A_MOD    �   �       XMAP    �  @   a   XMAP%IM    �  @   a   XMAP%JM      &  a   XMAP%LON1    (  $  a   XMAP%Q1    L  @   a   XMAP%IN    �  &  a   XMAP%LON2    �  $  a   XMAP%Q2    �  �       YMAP    j  @   a   YMAP%IM    �  @   a   YMAP%JM    �  f  a   YMAP%SIN1    P	  $  a   YMAP%Q1    t
  @   a   YMAP%JN    �
  f  a   YMAP%SIN2      $  a   YMAP%Q2    >  @   a   YMAP%IG    ~  @   a   YMAP%IV     �  �       READ_INPUT_GRID #   B  @   a   READ_INPUT_GRID%IM #   �  @   a   READ_INPUT_GRID%JM )   �  L   a   READ_INPUT_GRID%FILENAME *     &  a   READ_INPUT_GRID%LON_EDGES *   4  &  a   READ_INPUT_GRID%LAT_SINES    Z  �       DO_REGRID_A2A '   �  L   a   DO_REGRID_A2A%FILENAME !   >  @   a   DO_REGRID_A2A%IM !   ~  @   a   DO_REGRID_A2A%JM %   �  $  a   DO_REGRID_A2A%INGRID &   �  �   a   DO_REGRID_A2A%OUTGRID &   �  @   a   DO_REGRID_A2A%IS_MASS %   �  @   a   DO_REGRID_A2A%NETCDF      �       DO_REGRID_DKH '   �  L   a   DO_REGRID_DKH%FILENAME !   �  @   a   DO_REGRID_DKH%IM !   :  @   a   DO_REGRID_DKH%JM %   z  $  a   DO_REGRID_DKH%INGRID &   �  �   a   DO_REGRID_DKH%OUTGRID &   R  @   a   DO_REGRID_DKH%IS_MASS %   �  @   a   DO_REGRID_DKH%NETCDF    �  �       MAP_A2A    �  @   a   MAP_A2A%IM    �  @   a   MAP_A2A%JM      &  a   MAP_A2A%LON1    (  &  a   MAP_A2A%SIN1    N  $  a   MAP_A2A%Q1    r  @   a   MAP_A2A%IN    �  @   a   MAP_A2A%JN    �  &  a   MAP_A2A%LON2      &  a   MAP_A2A%SIN2    >   $  a   MAP_A2A%Q2    b!  @   a   MAP_A2A%IG    �!  @   a   MAP_A2A%IV 