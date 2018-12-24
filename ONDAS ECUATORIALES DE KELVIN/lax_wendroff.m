function [u_paso_entero, v_paso_entero, h_paso_entero] = lax_wendroff(dx, dy, dt, g, u, v, h, ...
					      SX_n_medio, SY_n_medio);
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%         ESQUEMA NUMERICO DE LAX WENDROF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Recordar cambios de variable
%PUNTOS MEDIOS EN EL TIEMPO Y EN EL ESPACIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh = u.*h;
vh = v.*h;

%h_medio_x_t = h_{i+1/2,j}^{n+1/2}


h_medio_x_t = 0.5.*(h(2:end,:)+h(1:end-1,:)) ...
  -(0.5*dt/dx).*(uh(2:end,:)-uh(1:end-1,:)); %h(2:end-1,:) es h_i,j

%h_medio_y_t = h_{i,j+1/2}^{n+1/2}

h_medio_y_t = 0.5.*(h(:,2:end)+h(:,1:end-1)) ...
  -(0.5*dt/dy).*(vh(:,2:end)-vh(:,1:end-1));

Ux = uh.*u+0.5.*g.*h.^2; 

Uy = uh.*v;


uh_medio_x_t = 0.5.*(uh(2:end,:)+uh(1:end-1,:)) ...
  -(0.5*dt/dx).*(Ux(2:end,:)-Ux(1:end-1,:));


uh_medio_y_t = 0.5.*(uh(:,2:end)+uh(:,1:end-1)) ...
  -(0.5*dt/dy).*(Uy(:,2:end)-Uy(:,1:end-1));

Vx = Uy; %de la definicion del cambio de variables


Vy = vh.*v+0.5.*g.*h.^2;

vh_medio_x_t = 0.5.*(vh(2:end,:)+vh(1:end-1,:)) ...
  -(0.5*dt/dx).*(Vx(2:end,:)-Vx(1:end-1,:));

vh_medio_y_t = 0.5.*(vh(:,2:end)+vh(:,1:end-1)) ...
  -(0.5*dt/dy).*(Vy(:,2:end)-Vy(:,1:end-1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                           USAMOS LOS PUNTOS MEDIOS PARA CALCULAR 
%%%%%%%                          LAS VARIABLES EN EL PROXIMO PASO TEMPORAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h_paso_entero = h(2:end-1,2:end-1) ...
  - (dt/dx).*(uh_medio_x_t(2:end,2:end-1)-uh_medio_x_t(1:end-1,2:end-1)) ...
  - (dt/dy).*(vh_medio_y_t(2:end-1,2:end)-vh_medio_y_t(2:end-1,1:end-1));

Ux_medio_x_t = uh_medio_x_t.*uh_medio_x_t./h_medio_x_t + 0.5.*g.*h_medio_x_t.^2;

Uy_medio_y_t = uh_medio_y_t.*vh_medio_y_t./h_medio_y_t;

uh_paso_entero = uh(2:end-1,2:end-1) ...
  - (dt/dx).*(Ux_medio_x_t(2:end,2:end-1)-Ux_medio_x_t(1:end-1,2:end-1)) ...
  - (dt/dy).*(Uy_medio_y_t(2:end-1,2:end)-Uy_medio_y_t(2:end-1,1:end-1)) ...
  + dt.*SX_n_medio.*0.5.*(h(2:end-1,2:end-1)+h_paso_entero);

Vx_medio_x_t = uh_medio_x_t.*vh_medio_x_t./h_medio_x_t;

Vy_medio_y_t = vh_medio_y_t.*vh_medio_y_t./h_medio_y_t + 0.5.*g.*h_medio_y_t.^2;

vh_paso_entero = vh(2:end-1,2:end-1) ...
  - (dt/dx).*(Vx_medio_x_t(2:end,2:end-1)-Vx_medio_x_t(1:end-1,2:end-1)) ...
  - (dt/dy).*(Vy_medio_y_t(2:end-1,2:end)-Vy_medio_y_t(2:end-1,1:end-1)) ...
  + dt.*SY_n_medio.*0.5.*(h(2:end-1,2:end-1)+h_paso_entero);

u_paso_entero = uh_paso_entero./h_paso_entero;
v_paso_entero = vh_paso_entero./h_paso_entero;

