%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                    Modelo de Shallow Water ----> ONDAS DE GRAVEDAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%               SWE en un CANAL utilizando el esquema numérico de LAX-WENDROFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------------------------------------------------------------
clear all; clc;
%%%%%% CONSTANTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ELEGIR CONDICIONES DE CONTORNO
promediadas = 1; 
no_desliz = 2; %u(0,y) = u(end,y) = v(x,0) = v(x,end) = 0 ---> NO SLIP CONDITION, NO DESLIZAMIENTO

%OPCION ELEGIDA:

conds_contorno = no_desliz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
g    = 9.81; 


dt = 120;   % PASO TEMPORAL DE 1 MINUTO (60 segs)
dt_output = 3600;  % TIEMPO ENTRE CADA OUTPUT DE 1 HORA (3600 segs), ENTRE CADA FRAME
dias_del_exp = 3; %se usa luego más abajo para mostrar los datos
dt_exp = dias_del_exp*24*3600.0;   % CANTIDAD DE DÍAS QUE DURA LA SIMULACIÓN en segundos -- forecast length


nx=254; % Number of zonal gridpoints
ny=150;  % Number of meridional gridpoints

dx=1.0e6; % ESPACIADO DE LA GRILLA EN X, en m ---> LATITUDES
dy=dx;      % ESPACIADO DE LA GRILLA EN Y, en m ---> MERIDIANOS


rango_alturas = [9500 10500]; % rango de alturas a graficar, en metros

nt = fix(dt_exp/dt)+1; % nt = NUMERO DE PASOS TEMPORALES, DE ITERACIONES. Se usa fix para redondear a un número entero.
npasos_entre_outputs = fix(dt_output/dt); % 1hora/4dias = paso entre frame y frame que se muestra
nframes = ceil(nt/npasos_entre_outputs); % NUMERO DE FRAMES QUE SE MUESTRAN

x=(0:nx-1).*dx; % coordenada x, en metros ---> DISTANCIA LATITUDINAL
y=(0:ny-1).*dy; % coordenada y, en metros ---> DISTANCIA MERIDIONAL
[Y,X] = meshgrid(y,x); %Creamos la grilla


%%%%% DEFINICIÓN DE H ---> OROGRAFÍA ELEGIDA, SUPERFICIE PLANA, MONTAÑA, ETC..
 H = zeros(nx, ny); %SUPERFICIE SÓLIDA DE ABAJO ES PLANA

%%% CONDICIONES INICIALES PARA LA ALTURA DEL FLUIDO
std_blob = 10.0.*dy; % Standard deviation of blob (m)
desplazamiento = 9700 + 1000.*exp(-((X-0.9.*mean(x)).^2+(Y-mean(y)).^2)./(2* ...
                                                     std_blob^2)); %(X-centro de la gaussiana en x))


   % CONDICIÓN INICIAL PARA EL CAMPO DE VELOCIDADES DEL VIENTO: INICIALMENTE EN REPOSO
u=zeros(nx, ny);
v=zeros(nx, ny);


% Altura del fluido respecto de la orografía
h = desplazamiento - H;

% Inicializamos las variables donde vamos a guardar la data para graficar
u_ = zeros(nx, ny, nframes);
v_ = zeros(nx, ny, nframes);
h_ = zeros(nx, ny, nframes);
t_ = zeros(1, nframes);

% indice p/guardar los datos que se va a ir actualizando con cada iteración
i_ = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                                  LOOP PRINCIPAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:nt %% iteramos para cada paso temporal
  %cada un número fijo de pasos nos guardamos las variables
  if mod(n-1,npasos_entre_outputs) == 0
    u_(:,:,i_) = u; %en la 1era iteración agarra la cond inicial, u0=0, luego se va actualizando!
    v_(:,:,i_) = v;
    h_(:,:,i_) = h;
    t_(i_) = (n-1).*dt;
    i_ = i_+1;
  end

  % Determinamos los términos fuente/sumidero SX,SY, con f=0 en este caso
  SX_n_medio = - (g/(2*dx)).*(H(3:end,2:end-1)-H(1:end-2,2:end-1));
  SY_n_medio = - (g/(2*dy)).*(H(2:end-1,3:end)-H(2:end-1,1:end-2));

  %CORREMOS EL LAX-WENDROFF PARA MOVERNOS 1 PASO EN EL TIEMPO
  [unext, vnext, hnext] = lax_wendroff(dx, dy, dt, g, u, v, h, ...
                                     SX_n_medio, SY_n_medio);

  % Actualizamos las velocidades y la altura
  u = unext([end 1:end 1],[1 1:end end]);
  v = vnext([end 1:end 1],[1 1:end end]);
  h(:,2:end-1) = hnext([end 1:end 1],:);

  %IMPONEMOS CONDICIONES DE CONTORNO, PISAMOS DONDE QUEREMOS
  % CONDS. DE CONTORNO PROMEDIADAS
  if conds_contorno == 1
     x_borde = ( h(2,:) + h(end-1,:) )./2;
     h(1,:) = x_borde;
     h(end,:) = x_borde;

     y_borde = ( h(:,2) + h(:,end-1) )./2;
     h(:,1) = y_borde;
     h(:,end) = y_borde;
     
    u_borde = ( u(2,:) + u(end-1,:) )./2;
    u(1,:) = u_borde;
    u(end,:) = u_borde;


    v_borde = ( v(:,2) + v(:,end-1) )./2;
    v(:,1) = v_borde;
    v(:,end) = v_borde;
  end
  
  if conds_contorno == 2
      v(:,[1 end]) = 0;
     u([1 end],:) = 0;
  end
      


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%                     ANIMACIÓN                                                           %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,'defaultaxesfontsize',20,...
    'paperpositionmode','auto','color','w');
drawnow

% EJES EN MILES DE km, x e y están en metros
x_1000km = x.*1e-6;
y_1000km = y.*1e-6;


ncol=128;
colormap(jet(ncol));
%colormap(hsv(ncol));

%intervalo para el ploteo entre cada vector velocidad
interval = 10; %6

if mean(rango_alturas) > 1000
  escala_altura = 0.001;
  titulo_altura = 'Altura (km)';
else
  escala_altura = 1;
  titulo_altura = 'Altura (m)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  figure(1)
 fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% h  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %LOOP PARA CADA FRAME DE LA ANIMACIÓN
    
for it = 1:nframes
  clf %comando: Clear current figure window

  % sacamos la altura y la componente de velocidad de cada frame, para cada it
  h = squeeze(h_(:,:,it));
  u = squeeze(u_(:,:,it));
  v = squeeze(v_(:,:,it));


  handle = image(x_1000km, y_1000km, (h'+H').*escala_altura); %image --> grafico de la matriz en colores
  set(handle,'CDataMapping','scaled');
  set(gca,'ydir','normal');
  caxis(rango_alturas.*escala_altura);

  hold on %clave
  contourf(x_1000km, y_1000km, H',[1:1000:8001],'k');
  

  %GRAFICAMOS EL CAMPO DE VELOCIDADES
  quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
         u(3:interval:end, 3:interval:end)',...
         v(3:interval:end, 3:interval:end)','k','linewidth',0.2);

 
  xlabel('X - Distancia latitudinal');
  ylabel('Y - Distancia meridional');
  title(['\bf' titulo_altura]);
  text(0, max(y_1000km), ['Tiempo = ' num2str(t_(it)./3600) ' Horas'],...
       'verticalalignment','bottom','fontsize',12);

  daspect([1 1 1]);
  axis([0 max(x_1000km) 0 max(y_1000km)]);
   colorbar
   drawnow
    F(it)=getframe(gcf);
    s = size(F(it).cdata);
    fprintf('%d %d\n', s(2), s(1))
   %eval(['print -dpng frame',num2str(it,'%02d'),'.png']); %me guardo cada iteración como una imagen ---> frame
   
   
   
end

video = VideoWriter('video_h_ondas_de_g.avi');
open(video)
writeVideo(video,F)
close(video)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vorticidad  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure(2)
  fig2=gcf;
fig2.Units='normalized';
fig2.OuterPosition=[0 0 1 1];

 %LOOP PARA CADA FRAME DE LA ANIMACIÓN
    
for it = 1:nframes
  clf %comando: Clear current figure window

  % sacamos la altura y la componente de velocidad de cada frame, para cada it
  h = squeeze(h_(:,:,it));
  u = squeeze(u_(:,:,it));
  v = squeeze(v_(:,:,it));
    
  
  vorticity = zeros(size(u));
  vorticity(2:end-1,2:end-1) = (1/dy).*(u(2:end-1,1:end-2)-u(2:end-1,3:end)) ...
     + (1/dx).*(v(3:end,2:end-1)-v(1:end-2,2:end-1));


   handle = image(x_1000km, y_1000km, vorticity');
  set(handle,'CDataMapping','scaled');
  set(gca,'ydir','normal');
  caxis([-3 3].*1e-4);

  
  xlabel('X - Distancia latitudinal');
  ylabel('Y - Distancia meridional');
  title('\bfVorticidad relativa(s^{-1})');
  text(0, max(y_1000km), ['Tiempo = ' num2str(t_(it)./3600) ' Horas'],...
       'verticalalignment','bottom','fontsize',12);

  % Other axes properties and plot a colorbar
  daspect([1 1 1]);
  axis([0 max(x_1000km) 0 max(y_1000km)]);
  colorbar
   drawnow
    Fv(it)=getframe(gcf);
        %2 lineas de abajo sirven para chequear el size de los frames si hay
    %problemas para hacer el video
    %sv = size(Fv(it).cdata);
    %fprintf('%d %d\n', sv(2), sv(1))
   %eval(['print -dpng frame',num2str(it,'%02d'),'.png']); %me guardo cada iteración como una imagen ---> frame
   
   
   
end

video_vort = VideoWriter('video_vort_ondas_de_gravedad.avi');
open(video_vort)
writeVideo(video_vort,Fv)
close(video_vort)