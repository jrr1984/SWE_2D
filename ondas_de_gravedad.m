%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                    Modelo de Shallow Water ----> ONDAS DE GRAVEDAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%               SWE en un CANAL utilizando el esquema numérico de LAX-WENDROFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------------------------------------------------------------
clear all; clc;
%%%%%% CONSTANTES
g    = 9.81; 
f = 0.; %fluido no rota
beta = 0;

dt = 60;   % PASO TEMPORAL DE 1 MINUTO (60 segs)
dt_output = 3600;  % TIEMPO ENTRE CADA OUTPUT DE 1 HORA (3600 segs), ENTRE CADA FRAME
dias_del_exp = 4; %se usa luego más abajo para mostrar los datos
dt_exp = dias_del_exp*24*3600.0;   % CANTIDAD DE DÍAS QUE DURA LA SIMULACIÓN en segundos -- forecast length


nx=254; % Number of zonal gridpoints
ny=50;  % Number of meridional gridpoints

dx=100.0e3; % ESPACIADO DE LA GRILLA EN X, en m ---> LATITUDES
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
std_blob = 8.0.*dy; % Standard deviation of blob (m)
   height = 9750 + 1000.*exp(-((X-0.25.*mean(x)).^2+(Y-mean(y)).^2)./(2* ...
                                                     std_blob^2));


% Coriolis parameter as a matrix of values varying in y only
F = f+beta.*(Y-mean(y));

% CONDICIÓN INICIAL PARA EL CAMPO DE VELOCIDADES DEL VIENTO: INICIALMENTE EN REPOSO
u=zeros(nx, ny);
v=zeros(nx, ny);





% Define h as the depth of the fluid (whereas "height" is the height of
% the upper surface)
h = height - H;

% Initialize the 3D arrays where the output data will be stored
u_save = zeros(nx, ny, nframes);
v_save = zeros(nx, ny, nframes);
h_save = zeros(nx, ny, nframes);
t_save = zeros(1, nframes);

% Index to stored data
i_save = 1;
% ------------------------------------------------------------------
% SECTION 3: Main loop

for n = 1:nt
  % Every fixed number of timesteps we store the fields
  if mod(n-1,npasos_entre_outputs) == 0
    max_u = sqrt(max(u(:).*u(:)+v(:).*v(:)));
    disp(['Time = ' num2str((n-1)*dt/3600) ...
	  ' hours (max ' num2str(dias_del_exp*24) ...
		   '); max(|u|) = ' num2str(max_u)]);
    u_save(:,:,i_save) = u;
    v_save(:,:,i_save) = v;
    h_save(:,:,i_save) = h;
    t_save(i_save) = (n-1).*dt;
    i_save = i_save+1;
  end

  % Compute the accelerations
  u_accel = F(2:end-1,2:end-1).*v(2:end-1,2:end-1) ...
              - (g/(2*dx)).*(H(3:end,2:end-1)-H(1:end-2,2:end-1));
  v_accel = -F(2:end-1,2:end-1).*u(2:end-1,2:end-1) ...
              - (g/(2*dy)).*(H(2:end-1,3:end)-H(2:end-1,1:end-2));

  % Call the Lax-Wendroff scheme to move forward one timestep
  [unew, vnew, h_new] = lax_wendroff(dx, dy, dt, g, u, v, h, ...
                                     u_accel, v_accel);

  % Update the wind and height fields, taking care to enforce 
  % boundary conditions 
  u = unew([end 1:end 1],[1 1:end end]);
  v = vnew([end 1:end 1],[1 1:end end]);
  v(:,[1 end]) = 0;
  h(:,2:end-1) = h_new([end 1:end 1],:);

end

disp('Now run "animate" to animate the simulation');



% This script animates the height field and the vorticity produced by
% a shallow water model. It should be called only after shallow_water_model
% has been run.

% Set the size of the figure
set(gcf,'units','inches');
pos=get(gcf,'position');
pos([3 4]) = [10.5 5];
set(gcf,'position',pos)

% Set other figure properties and draw now
set(gcf,'defaultaxesfontsize',12,...
    'paperpositionmode','auto','color','w');
drawnow

% Axis units are thousands of kilometers (x and y are in metres)
x_1000km = x.*1e-6;
y_1000km = y.*1e-6;

% Set colormap to have 64 entries
ncol=64;
colormap(jet(ncol));
% colormap(hclmultseq01);
% colormap(flipud(cbrewer('div','RdBu',32)))

% Interval between arrows in the velocity vector plot
interval = 6;

% Decide whether to show height in metres or km
if mean(rango_alturas) > 1000
  height_scale = 0.001;
  height_title = 'Height (km)';
else
  height_scale = 1;
  height_title = 'Height (m)';
end




 figure(1)
% Loop through the frames of the animation
for it = 1:nframes
  clf

  % Extract the height and velocity components for this frame
  h = squeeze(h_save(:,:,it));
  u = squeeze(u_save(:,:,it));
  v = squeeze(v_save(:,:,it));

  % First plot the height field
 
  % Plot the height field
  handle = image(x_1000km, y_1000km, (h'+H').*height_scale);
  set(handle,'CDataMapping','scaled');
  set(gca,'ydir','normal');
  caxis(rango_alturas.*height_scale);

  % Plot the orography as black contours every 1000 m
  hold on
  warning off
  contour(x_1000km, y_1000km, H',[1:1000:8001],'k');
  warning on

  % Plot the velocity vectors
  quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
         u(3:interval:end, 3:interval:end)',...
         v(3:interval:end, 3:interval:end)','k','linewidth',0.2);

  % Write the axis labels, title and time of the frame
  xlabel('X distance (1000s of km)');
  ylabel('Y distance (1000s of km)');
  title(['\bf' height_title]);
  text(0, max(y_1000km), ['Time = ' num2str(t_save(it)./3600) ' hours'],...
       'verticalalignment','bottom','fontsize',12);

  % Set other axes properties and plot a colorbar
  daspect([1 1 1]);
  axis([0 max(x_1000km) 0 max(y_1000km)]);
  colorbar
  
  % Compute the vorticity
  vorticity = zeros(size(u));
  vorticity(2:end-1,2:end-1) = (1/dy).*(u(2:end-1,1:end-2)-u(2:end-1,3:end)) ...
     + (1/dx).*(v(3:end,2:end-1)-v(1:end-2,2:end-1));

  

  % Other axes properties and plot a colorbar
  daspect([1 1 1]);
  axis([0 max(x_1000km) 0 max(y_1000km)]);
  colorbar

   drawnow
  
   eval(['print -dpng frame',num2str(it,'%02d'),'.png']); %me guardo cada iteración como una imagen ---> frame
   
end


%%%%%%%%%%%%%%%%%%% VIDEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the images
 images    = cell(nframes,1);
 for i=1:length(images)
     images{i} = imread(sprintf('/home/jrr/Descargas/FINALES E1 MINOTTI/lax_w_matlab/frame%02d.png',i));
 end
 % create the video writer with 30 fps
 writerObj = VideoWriter('video.avi');
 writerObj.FrameRate = 10;
   % open the video writer
   open(writerObj);
   % write the frames to the video
    for u=1:10
       % convert the image to a frame
       frame = im2frame(images{u});
       for v=1:10
           writeVideo(writerObj, frame);
       end
   end
   % close the writer object
   close(writerObj);
  % implay('Velocity.avi');
