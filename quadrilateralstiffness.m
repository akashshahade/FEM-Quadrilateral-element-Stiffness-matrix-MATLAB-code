
clc;clear;

% Material properties
 
E = 1; mu = 0.25;
t = 1;  % thickness

% Calculate D matrix
q = E/(1 - (mu^2)) ;
D = q*[1 mu 0; mu 1 0; 0 0 (1 - mu)*0.5];   
    %node coordinates

xn1 = 2; xn2 = 5; xn3 = 4; xn4 = 1;
yn1 = 1; yn2 = 2; yn3 = 6; yn4 = 4;

%--------------------------------------------------------------------------------------------%

% Gauss Point 1
    
    zeta1 = -0.57735; eta1 = -0.57735;
   
    
    j111 = 0.25*((xn3 - xn4) - (xn1 - xn2)) + eta1*0.25*((xn3 - xn4) + (xn1 - xn2)) ;
    j121 = 0.25*((yn3 - yn4) - (yn1 - yn2)) + eta1*0.25*((yn3 - yn4) + (yn1 - yn2)) ;
    j211 = 0.25*((xn3 - xn2) + (xn4 - xn1)) + zeta1*0.25*((xn3 - xn2) - (xn4 - xn1)) ;
    j221 = 0.25*((yn3 - yn2) + (yn4 - yn1)) + zeta1*0.25*((yn3 - yn2) - (yn4 - yn1)) ;
   
    jacob1 = [j111 j121; j211 j221];
    
    %determinent of jacobian
    
    dj1 = det(jacob1);
    dji1 = 1/dj1 ;
    
    % [A] matrix
    
    a1 = dji1*[jacob1(2,2) -jacob1(1,2) 0 0;0 0 -jacob1(2,1) jacob1(1,1); -jacob1(2,1) jacob1(1,1) jacob1(2,2) -jacob1(1,2)];
    
    %[G] matrix
    
    gm1 = [-(1-eta1) 0 (1-eta1) 0 (1+eta1) 0 -(1+eta1) 0;...
        -(1-zeta1) 0 -(1+zeta1) 0 (1+zeta1) 0 (1-zeta1) 0;...
        0 -(1-eta1) 0 (1-eta1) 0 (1+eta1) 0 -(1+eta1);...
        0 -(1-zeta1) 0 -(1+zeta1) 0 (1+zeta1) 0 (1-zeta1)];
    gmat1 = 0.25*gm1;
    
    %[B] matrix calculation
    
    bmat1 = a1*gmat1;
    
    % elemental [k] calculation for 1st gauss point
    
    keg1 = bmat1.'*D*bmat1*dj1 ;
    
  %--------------------------------------------------------------------------------------------%  
    
    % Gauss Point 2
    
    zeta2 = -0.57735; eta2 = 0.57735;
   
    
    j112 = 0.25*((xn3 - xn4) - (xn1 - xn2)) + eta2*0.25*((xn3 - xn4) + (xn1 - xn2)) ;
    j122 = 0.25*((yn3 - yn4) - (yn1 - yn2)) + eta2*0.25*((yn3 - yn4) + (yn1 - yn2)) ;
    j212 = 0.25*((xn3 - xn2) + (xn4 - xn1)) + zeta2*0.25*((xn3 - xn2) - (xn4 - xn1)) ;
    j222 = 0.25*((yn3 - yn2) + (yn4 - yn1)) + zeta2*0.25*((yn3 - yn2) - (yn4 - yn1)) ;
   
    jacob2 = [j112 j122; j212 j222];
    
    %determinent of jacobian
    
    dj2 = det(jacob2);
    dji2 = 1/dj2 ;
    
    % [A] matrix
    
    a2 = dji2*[jacob2(2,2) -jacob2(1,2) 0 0;0 0 -jacob2(2,1) jacob2(1,1); -jacob2(2,1) jacob2(1,1) jacob2(2,2) -jacob2(1,2)];
    
    %[G] matrix
    
    gm2 = [-(1-eta2) 0 (1-eta2) 0 (1+eta2) 0 -(1+eta2) 0;...
        -(1-zeta2) 0 -(1+zeta2) 0 (1+zeta2) 0 (1-zeta2) 0;...
        0 -(1-eta2) 0 (1-eta2) 0 (1+eta2) 0 -(1+eta2);...
        0 -(1-zeta2) 0 -(1+zeta2) 0 (1+zeta2) 0 (1-zeta2)];
    gmat2 = 0.25*gm2;
    
    %[B] matrix calculation
    
    bmat2 = a2*gmat2;
    
    % elemental [k] calculation for 2nd gauss point
    
    keg2 = bmat2.'*D*bmat2*dj2 ;
    
   %--------------------------------------------------------------------------------------------% 
    
    % Gauss Point 3
    
    zeta3 = 0.57735; eta3 = 0.57735;
   
    
    j113 = 0.25*((xn3 - xn4) - (xn1 - xn2)) + eta3*0.25*((xn3 - xn4) + (xn1 - xn2)) ;
    j123 = 0.25*((yn3 - yn4) - (yn1 - yn2)) + eta3*0.25*((yn3 - yn4) + (yn1 - yn2)) ;
    j213 = 0.25*((xn3 - xn2) + (xn4 - xn1)) + zeta3*0.25*((xn3 - xn2) - (xn4 - xn1)) ;
    j223 = 0.25*((yn3 - yn2) + (yn4 - yn1)) + zeta3*0.25*((yn3 - yn2) - (yn4 - yn1)) ;
   
    jacob3 = [j113 j123; j213 j223];
    
    %determinent of jacobian
    
    dj3 = det(jacob3);
    dji3 = 1/dj3 ;
    
    % [A] matrix
    
    a3 = dji3*[jacob3(2,2) -jacob3(1,2) 0 0;0 0 -jacob3(2,1) jacob3(1,1); -jacob3(2,1) jacob3(1,1) jacob3(2,2) -jacob3(1,2)];
    
    %[G] matrix
    
    gm3 = [-(1-eta3) 0 (1-eta3) 0 (1+eta3) 0 -(1+eta3) 0;...
        -(1-zeta3) 0 -(1+zeta3) 0 (1+zeta3) 0 (1-zeta3) 0;...
        0 -(1-eta3) 0 (1-eta3) 0 (1+eta3) 0 -(1+eta3);...
        0 -(1-zeta3) 0 -(1+zeta3) 0 (1+zeta3) 0 (1-zeta3)];
    gmat3 = 0.25*gm3;
    
    %[B] matrix calculation
    
    bmat3 = a3*gmat3;
    
    % elemental [k] calculation for 1st gauss point
    
    keg3 = bmat3.'*D*bmat3*dj3 ;
    
    
    %--------------------------------------------------------------------------------------------%
    
    % Gauss Point 4
    
    zeta4 = 0.57735; eta4 = -0.57735;
   
    
    j114 = 0.25*((xn3 - xn4) - (xn1 - xn2)) + eta4*0.25*((xn3 - xn4) + (xn1 - xn2)) ;
    j124 = 0.25*((yn3 - yn4) - (yn1 - yn2)) + eta4*0.25*((yn3 - yn4) + (yn1 - yn2)) ;
    j214 = 0.25*((xn3 - xn2) + (xn4 - xn1)) + zeta4*0.25*((xn3 - xn2) - (xn4 - xn1)) ;
    j224 = 0.25*((yn3 - yn2) + (yn4 - yn1)) + zeta4*0.25*((yn3 - yn2) - (yn4 - yn1)) ;
   
    jacob4 = [j114 j124; j214 j224];
    
    %determinent of jacobian
    
    dj4 = det(jacob4);
    dji4 = 1/dj4 ;
    
    % [A] matrix
    
    a4 = dji4*[jacob4(2,2) -jacob4(1,2) 0 0;0 0 -jacob4(2,1) jacob4(1,1); -jacob4(2,1) jacob4(1,1) jacob4(2,2) -jacob4(1,2)];
    
    %[G] matrix
    
    gm4 = [-(1-eta4) 0 (1-eta4) 0 (1+eta4) 0 -(1+eta4) 0;...
        -(1-zeta4) 0 -(1+zeta4) 0 (1+zeta4) 0 (1-zeta4) 0;...
        0 -(1-eta4) 0 (1-eta4) 0 (1+eta4) 0 -(1+eta4);...
        0 -(1-zeta4) 0 -(1+zeta4) 0 (1+zeta4) 0 (1-zeta4)];
    gmat4 = 0.25*gm4;
    
    %[B] matrix calculation
    
    bmat4 = a4*gmat4;
    
    % elemental [k] calculation for 1st gauss point
    
    keg4 = bmat4.'*D*bmat4*dj4 ;
    
    
    %--------------------------------------------------------------------------------------------%
    % Final elemental stiffness matrix
    %--------------------------------------------------------------------------------------------%
    
    kel = keg1 + keg2 + keg3 + keg4 ;
    
    
    
