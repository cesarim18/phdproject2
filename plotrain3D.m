function plotrain3D

nx = 128; ny=128;  m=100+1;

Tfinal = 2.0;
Tfinal = Tfinal*60
%stop

Lx = 128;
Ly = 128;
Lz = 15;

dx = Lx/nx; dy=Ly/ny; dz=Lz/(m-1); dz_u=Lz/m;

x = 0:dx:Lx-dx;% Lx-Lx/nx=(1-1/nx)*Lx=((nx-1)/nx)*Lx
y = 0:dy:Ly-dy;
z = 0:dz:Lz ;
z_u = 0:dz_u:Lz ;

[xGrid, yGrid, zGrid] = meshgrid(y,x,z);
%%
display('Movie: QrYHalf')
WYHalf = load('Run_36/Data2D/WYHalf.dat');
qrYHalf = load('Run_36/Data2D/QrYHalf.dat');
qvYHalf = load('Run_36/Data2D/QvYHalf.dat');
qrZavg = load('Run_36/Data2D/QrZavg.dat');
qrZMax = load('Run_36/Data2D/QrZMax.dat');
movies = size(WYHalf,1)/nx/m
w = reshape(WYHalf,nx,m,movies);
NoStepsMov=Tfinal/200;%(movies-1);
MovAvi5 = VideoWriter('Run_36/Videos/WYHalf.avi');
MovAvi5.Quality = 100;
MovAvi5.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    clf
    figure(1)
    contour(x,z,transpose(w(:,:,mov)),100);
	colorbar;
    
    axis([0 Lx 0 Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['w (x,z) contour at y=Ly/2 in m/s, T = ', num2str(Tmov),' minutes'])
    %print('Frame','-depsc',figure(1))
    
	mov
	open(MovAvi5)
    F = getframe(figure(1));
    writeVideo(MovAvi5,F);
end
%stop
%% %% %%%%%%%Qr y-half
display('Movie: QrYHalf')
size(qrYHalf)
display('file loaded')
movies = size(qrYHalf,1)/nx/m
%stop

NoStepsMov=Tfinal/100;%(movies-1);

qr = reshape(qrYHalf,nx,m,movies);
display('Movie: QvYHalf')
size(qvYHalf)
display('file loaded')
movies = size(qvYHalf,1)/nx/m
%stop

%NoStepsMov=Tfinal/(movies-2);%(movies-1);

qv = reshape(qvYHalf,nx,m,movies);

MovAvi5 = VideoWriter('Run_36/Videos/QrYHalf.avi');
MovAvi5.Quality = 100;
MovAvi5.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    clf
    figure(1)
    contour(x,z,transpose(qr(:,:,mov)),100);
	colorbar;
    
    axis([0 Lx 0 Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['qr (x,z) contour at y=Ly/2 in g/Kg, T = ', num2str(Tmov),' minutes'])
    %print('Frame','-depsc',figure(1))
    
	mov
	open(MovAvi5)
    F = getframe(figure(1));
    writeVideo(MovAvi5,F);
end

close(MovAvi5);
%stop
% 
% %% %%%%%%%Qv y-half
% display('Movie: QvYHalf')
% qvYHalf = load('Run_36/Data2D/QvYHalf.dat');
% size(qvYHalf)
% display('file loaded')
% movies = size(qvYHalf,1)/nx/m
% %stop
% 
% NoStepsMov=Tfinal/(movies-2);%(movies-1);
% 
% qv = reshape(qvYHalf,nx,m,movies);
qt = qv + qr;
%qv_bg = load('Run_36/Data1D/Qt_bg.dat');
%for iz=1:m
%    qv(:,iz,:) = qvYHalf(:,iz,:) - qv_bg(iz,1);
%end

MovAvi5 = VideoWriter('Run_36/Videos/QtYHalf.avi');
MovAvi5.Quality = 100;
MovAvi5.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    clf
    figure(1)
    contour(x,z,transpose(qt(:,:,mov)),100);
	colorbar;
    
    axis([0 Lx 0 Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['qt (x,z) contour at y=Ly/2 in g/Kg, T = ', num2str(Tmov),' minutes'])
    %print('Frame','-depsc',figure(1))
    
	mov
	open(MovAvi5)
    F = getframe(figure(1));
    writeVideo(MovAvi5,F);
end

close(MovAvi5);


MovAvi5 = VideoWriter('Run_36/Videos/QvYHalf.avi');
MovAvi5.Quality = 100;
MovAvi5.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    clf
    figure(1)
    contour(x,z,transpose(qv(:,:,mov)),100);
	colorbar;
    
    axis([0 Lx 0 Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['qv (x,z) contour at y=Ly/2 in g/Kg, T = ', num2str(Tmov),' minutes'])
    %print('Frame','-depsc',figure(1))
    
	mov
	open(MovAvi5)
    F = getframe(figure(1));
    writeVideo(MovAvi5,F);
end

close(MovAvi5);

%% %%%% Qr z-avg
display('Movie: QrZavg')
display('file loaded')
movies = size(qrZavg,1)/nx/ny

qrZavg = reshape(qrZavg,nx,ny,movies);

MovAvi0 = VideoWriter('Run_36/Videos/QrZavg.avi');
MovAvi0.Quality = 100;
MovAvi0.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);

    clf
    figure(1)
    contour(x,y,transpose(qrZavg(:,:,mov)),100); colorbar;
    
    axis([0 Lx 0 Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-avg qr (x,y) contour in g/Kg, T = ', num2str(Tmov),' minutes'])
    
	mov
	open(MovAvi0)
	F = getframe(figure(1));
	writeVideo(MovAvi0,F);
end

close(MovAvi0);


%% %%%%% Qr z-max:
display('Movie: QrZMax')
display('file loaded')
movies = size(qrZMax,1)/nx/ny

qrZMax = reshape(qrZMax,nx,ny,movies);

MovAvi1 = VideoWriter('Run_36/Videos/QrZMax.avi');
MovAvi1.Quality = 100;
MovAvi1.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);

    clf
    figure(1)
    contour(x,y,transpose(qrZMax(:,:,mov)),100); colorbar;
    
    axis([0 Lx 0 Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-max qr (x,y) contour in g/Kg, T = ', num2str(Tmov),' minutes'])
    
	mov
	open(MovAvi1)
    F = getframe(figure(1));
    writeVideo(MovAvi1,F);
end

close(MovAvi1);
stop
% 
%% %%%%%%%U y-half
% display('Movie: UYHalf')
% UYHalf = load('Run_36/Data2D/UYHalf.dat');
% size(UYHalf)
% 
% display('file loaded')
% movies = size(UYHalf,1)/nx/(m+1)
% 
% NoStepsMov=0.5;%Tfinal/movies;%(movies-1);
% 
% UYHalf = reshape(UYHalf,nx,m+1,movies);
% 
% MovAvi5 = VideoWriter('Run_36/Videos/UYHalf.avi');
% MovAvi5.Quality = 100;
% MovAvi5.FrameRate = 1;
% for mov=1:movies-1
%     
%     Tmov=NoStepsMov*(mov-1);
%     clf
%     figure(1)
%     contour(x,z_u,transpose(UYHalf(:,:,mov)),100);
% 	colorbar;
%     %caxis([.1 5])
%     
%     axis([0 Lx 0 Lz])
%     axis on
%     xlabel('x in km');
%     ylabel('z in km');    
%     
%     title(['u (x,z) contour at y=Ly/2 in m/s, T = ', num2str(Tmov),' minutes'])
%     print('Frame','-depsc',figure(1))
%     
% 	mov
% 	open(MovAvi5)
%     F = getframe(figure(1));
%     writeVideo(MovAvi5,F);
% end
% 
% close(MovAvi5);

%% %%%%%%%V y-half
% display('Movie: VYHalf')
% VYHalf = load('Run_36/Data2D/VYHalf.dat');
% size(VYHalf)
% 
% display('file loaded')
% movies = size(VYHalf,1)/nx/(m+1)
% 
% NoStepsMov=0.5;%Tfinal/movies;%(movies-1);
% 
% VYHalf = reshape(VYHalf,nx,m+1,movies);
% 
% MovAvi5 = VideoWriter('Run_36/Videos/VYHalf.avi');
% MovAvi5.Quality = 100;
% MovAvi5.FrameRate = 1;
% for mov=1:movies-1
%     
%     Tmov=NoStepsMov*(mov-1);
%     clf
%     figure(1)
%     contour(x,z_u,transpose(VYHalf(:,:,mov)),100);
% 	colorbar;
%     %caxis([.1 5])
%     
%     axis([0 Lx 0 Lz])
%     axis on
%     xlabel('x in km');
%     ylabel('z in km');    
%     
%     title(['v (x,z) contour at y=Ly/2 in m/s, T = ', num2str(Tmov),' minutes'])
%     print('Frame','-depsc',figure(1))
%     
% 	mov
% 	open(MovAvi5)
%     F = getframe(figure(1));
%     writeVideo(MovAvi5,F);
% end
% 
% close(MovAvi5);

%% %%%%%%% W y-half
display('Movie: WYHalf')
WYHalf = load('Run_36/Data2D/WYHalf.dat');
size(WYHalf)

display('file loaded')
movies = size(WYHalf,1)/nx/m

NoStepsMov=0.5;%Tfinal/movies;%(movies-1);

WYHalf = reshape(WYHalf,nx,m,movies);

MovAvi5 = VideoWriter('Run_36/Videos/WYHalf.avi');
MovAvi5.Quality = 100;
MovAvi5.FrameRate = 1;
for mov=1:movies-1
    
    Tmov=NoStepsMov*(mov-1);
    clf
    figure(1)
    contour(x,z,transpose(WYHalf(:,:,mov)),100);
	colorbar;
    %caxis([.1 5])
    
    axis([0 Lx 0 Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['w (x,z) contour at y=Ly/2 in m/s, T = ', num2str(Tmov),' minutes'])
    print('Frame','-depsc',figure(1))
    
	mov
	open(MovAvi5)
    F = getframe(figure(1));
    writeVideo(MovAvi5,F);
end

close(MovAvi5);




stop

%%

Th = load('Run_36/Data3D/Th_6.dat');
Th = reshape(Th,nx,ny,m);


xslice=10*(Lx-dx);% 10*Lx/2; %0:10*Lx/128:10*Lx;%
yslice=10*(Ly-dy); %0:10*Ly/32:10*(Ly-dy); %10*Lx/2; %0:10*Ly/128:10*Ly; %
zslice=0:10*Lz/100:10*Lz; 

h=slice(10*xGrid,10*yGrid,10*zGrid,3*Th,yslice,xslice,zslice); colorbar;
%caxis([.1*10^(-3) 2.5*10^(-3)])
view(78,54)
axis tight

alpha('color')
set(h,'EdgeColor','none','FaceColor','interp',...
 'FaceAlpha','interp')

stop

Th_y = zeros(nx,m);
figure(1)
Th_y(:,:) = 0.5*(Th(:,ny/2 -1,:)+Th(:,ny/2+1,:));
contour(10*x,10*z,3*transpose(Th_y(:,:)),100); colorbar;

figure(2)
Th_y(:,:) = 0.5*(Th(:,ny/2,:)+Th(:,ny/2,:));
contour(10*x,10*z,3*transpose(Th_y(:,:)),100); colorbar;
%stop


%% %%%%%%%%%% Theta at low levels
display('Movie: Theta at 5 km')
thetaZpt5Km = load('Data2D/ThetaZpt5km_88.dat');
movies = size(thetaZpt5Km,1)/nx/ny
thetaZpt5Km = reshape(thetaZpt5Km,nx,ny,movies);
MovAvi9 = VideoWriter('ThetaZpt5Km2.avi');
MovAvi9.Quality = 100;
MovAvi9.FrameRate = 5;
NoStepsMov = 1;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*1; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,transpose(3*thetaZpt5Km(:,:,mov)),100); colorbar;
    %caxis([-0.001 0.001])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['theta (x,y) contour at z=0.5 km in Kelvin, T = ', num2str(Tmov),' minutes'])
    
	mov
	open(MovAvi9)
    F = getframe(figure(1));
    writeVideo(MovAvi9,F);
end

close(MovAvi9);
stop
%%
% revisar graficas
% u_vort = load('Data2D/U_ZFix_95.dat');
% mov_u = size(u_vort,1)/nx/ny
% u_vort = reshape(u_vort,nx,ny,mov_u);
% v_vort = load('Data2D/V_ZFix_95.dat');
% mov_v = size(v_vort,1)/nx/ny
% v_vort = reshape(v_vort,nx,ny,mov_v);
% vort = load('Data2D/VoZFix_95.dat');
% mov_vo = size(vort,1)/nx/ny
% vort = reshape(vort,nx,ny,mov_vo);
% 
% movies = mov_vo;%min(mov_u,mov_v,mov_vo)
% 
% time=zeros(movies,1);
% NoStepsMov=Tfinal/(movies-1);
% for mov=1:movies
%     time(mov)=NoStepsMov*(mov-1)/60;
% end
% 
% MovAvi3 = VideoWriter('vort95.avi');
% MovAvi3.Quality = 100;
% MovAvi3.FrameRate = 5;
% for mov=1:movies
%     
%     Tmov=NoStepsMov*(mov-1);
%     Tmov=Tmov*1; %In real time (hours)
% 
%     clf
%     figure(1)
%     contour(10*x,10*y,Us*transpose(vort(:,:,mov)),100); colorbar;
%     %caxis([-0.001 0.001])
%     hold on
% 	q = quiver(10*x,10*y,Us*transpose(u_vort(:,:,mov)),Us*transpose(v_vort(:,:,mov)),'b')
%     %q.AutoScaleFactor = 2
%     hold off
%     axis([0 10*Lx 0 10*Ly])
%     %axis on
%     xlabel('x in km');
%     ylabel('y in km');    
%     
%     title(['MultiScale FARE z=7.5 km vorticity (x,y), T = ', num2str(Tmov),' minutes'])
% 	mov
% 	open(MovAvi3)
%     F = getframe(figure(1));
%     writeVideo(MovAvi3,F);
% end
% close(MovAvi3);
% stop



%% %%%%% Qr y-half
display('Movie: QrYHalf')
qrYHalf = load('Data2D/QrYHalf_90.dat');
qvYHalf = load('Data2D/QvYHalf_90.dat');
display('file loaded')
movies = size(qrYHalf,1)/nx/m

NoStepsMov=Tfinal/movies;%movies-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)/60;
end

qrYHalf = reshape(qrYHalf,nx,m,movies);
qvYHalf = reshape(qvYHalf,nx,m,movies);
qtYHalf = qrYHalf + qvYHalf;
MovAvi2 = VideoWriter('QtYHalf_MS.avi');
MovAvi2.Quality = 100;
MovAvi2.FrameRate = 5;

for mov=1:movies
    
    Tmov=2.25*mov;%NoStepsMov*(mov-1);
    %Tmov=Tmov*1; %In real time (hours)

    clf
    figure(1)
    subplot(131)
        contour(10*x,10*z,qs*transpose(qrYHalf(:,:,mov)),100); colorbar;
        %caxis([.1 5])

        %caxis([0 5*10^(-3)])
        %caxis([-2*10^(-4) 6*10^(-4)])

        axis([0 10*Lx 0 10*Lz])
        axis on
        xlabel('x in km');
        ylabel('z in km');    

        %title(['qr (x,z) contour at y=Ly/2 in g/Kg, T = ', num2str(Tmov),' seconds'])
        title(['qr (x,z)']);% contour at y=Ly/2 in g/Kg, T = ', num2str(Tmov),' seconds'])
    subplot(132)
        contour(10*x,10*z,qs*transpose(qvYHalf(:,:,mov)),100); colorbar;
        %caxis([.1 5])

        %caxis([0 5*10^(-3)])
        %caxis([-2*10^(-4) 6*10^(-4)])

        axis([0 10*Lx 0 10*Lz])
        axis on
        xlabel('x in km');
        ylabel('z in km');    

        title(['qv (x,z) T = ', num2str(Tmov),' seconds'])

    subplot(133)
       contour(10*x,10*z,qs*transpose(qtYHalf(:,:,mov)),100); colorbar;
         %caxis([.1 5])

        %caxis([0 5*10^(-3)])
        %caxis([-2*10^(-4) 6*10^(-4)])

        axis([0 10*Lx 0 10*Lz])
        axis on
        xlabel('x in km');
        ylabel('z in km');    

        title(['qt (x,z)']);% contour at y=Ly/2 in g/Kg, T = ', num2str(Tmov),' seconds'])
        
        
    print('Frame','-depsc',figure(1))

   
	mov
	open(MovAvi2)
    F = getframe(figure(1));
    writeVideo(MovAvi2,F);
end

close(MovAvi2);
stop

%%
% revisar graficas
% u_vort = load('Data2D/U_ZFix_95.dat');
% mov_u = size(u_vort,1)/nx/ny
% u_vort = reshape(u_vort,nx,ny,mov_u);
% v_vort = load('Data2D/V_ZFix_95.dat');
% mov_v = size(v_vort,1)/nx/ny
% v_vort = reshape(v_vort,nx,ny,mov_v);
% vort = load('Data2D/VoZFix_95.dat');
% mov_vo = size(vort,1)/nx/ny
% vort = reshape(vort,nx,ny,mov_vo);
% 
% movies = mov_vo;%min(mov_u,mov_v,mov_vo)
% 
% time=zeros(movies,1);
% NoStepsMov=Tfinal/(movies-1);
% for mov=1:movies
%     time(mov)=NoStepsMov*(mov-1)/60;
% end
% 
% MovAvi3 = VideoWriter('vort95.avi');
% MovAvi3.Quality = 100;
% MovAvi3.FrameRate = 5;
% for mov=1:movies
%     
%     Tmov=NoStepsMov*(mov-1);
%     Tmov=Tmov*1; %In real time (hours)
% 
%     clf
%     figure(1)
%     contour(10*x,10*y,Us*transpose(vort(:,:,mov)),100); colorbar;
%     %caxis([-0.001 0.001])
%     hold on
% 	q = quiver(10*x,10*y,Us*transpose(u_vort(:,:,mov)),Us*transpose(v_vort(:,:,mov)),'b')
%     %q.AutoScaleFactor = 2
%     hold off
%     axis([0 10*Lx 0 10*Ly])
%     %axis on
%     xlabel('x in km');
%     ylabel('y in km');    
%     
%     title(['MultiScale FARE z=7.5 km vorticity (x,y), T = ', num2str(Tmov),' minutes'])
% 	mov
% 	open(MovAvi3)
%     F = getframe(figure(1));
%     writeVideo(MovAvi3,F);
% end
% close(MovAvi3);
% stop
% %%
% % revisar graficas
% u_vort = load('Data2D/U_ZFix_97.dat');
% mov_u = size(u_vort,1)/nx/ny
% u_vort = reshape(u_vort,nx,ny,mov_u);
% v_vort = load('Data2D/V_ZFix_97.dat');
% mov_v = size(v_vort,1)/nx/ny
% v_vort = reshape(v_vort,nx,ny,mov_v);
% vort = load('Data2D/VoZFix_97.dat');
% mov_vo = size(vort,1)/nx/ny
% vort = reshape(vort,nx,ny,mov_vo);
% 
% movies = mov_vo;%min(mov_u,mov_v,mov_vo)
% 
% time=zeros(movies,1);
% NoStepsMov=Tfinal/(movies-1);
% for mov=1:movies
%     time(mov)=NoStepsMov*(mov-1)/60;
% end
% 
% MovAvi3 = VideoWriter('vort97.avi');
% MovAvi3.Quality = 100;
% MovAvi3.FrameRate = 5;
% for mov=1:movies
%     
%     Tmov=NoStepsMov*(mov-1);
%     Tmov=Tmov*1; %In real time (hours)
% 
%     clf
%     figure(1)
%     contour(10*x,10*y,Us*transpose(vort(:,:,mov)),100); colorbar;
%     %caxis([-0.001 0.001])
%     hold on
% 	q = quiver(10*x,10*y,Us*transpose(u_vort(:,:,mov)),Us*transpose(v_vort(:,:,mov)),'b')
%     %q.AutoScaleFactor = 2
%     hold off
%     axis([0 10*Lx 0 10*Ly])
%     %axis on
%     xlabel('x in km');
%     ylabel('y in km');    
% 
%      title(['FARE z=7.5 km vorticity (x,y), T = ', num2str(Tmov),' minutes'])
% 
% 	mov
% 	open(MovAvi3)
%     F = getframe(figure(1));
%     writeVideo(MovAvi3,F);
% end
% close(MovAvi3);
% 
% stop

%	
% quiver(x,y,x,y)
% 
% vorticity = load('Data2D/VoZFix_13.dat');
% 
% 
% 
% stop
%% %%%%%%% Qv z-max
display('Movie: QvZMax')
qvZMax = load('Data2D/QvZMax_95.dat');

display('file loaded')
movies = size(qvZMax,1)/nx/ny

NoStepsMov=Tfinal/movies;%(movies-1);
% NoStepsMov=Tfinal/(500-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)/60;
end


qvZMax = reshape(qvZMax,nx,ny,movies);

MovAvi3 = VideoWriter('QvZMax_95.avi');
MovAvi3.Quality = 100;
MovAvi3.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*1; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,qs*transpose(qvZMax(:,:,mov)),100); colorbar;
    %caxis([.1 5])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-max qv (x,y) contour in g/Kg, T = ', num2str(Tmov),' minutes'])
    
	mov
	open(MovAvi3)
    F = getframe(figure(1));
    writeVideo(MovAvi3,F);
end

close(MovAvi3);

stop
%% %%%%%%%% Qv z-avg
display('Movie: QvZavg')
qvZavg = load('Data2D/QvZavg_13.dat');
display('file loaded')
movies = size(qvZavg,1)/nx/ny

NoStepsMov=Tfinal/(movies-1);
% NoStepsMov=Tfinal/(500-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;

end

qvZavg = reshape(qvZavg,nx,ny,movies);

MovAvi4 = VideoWriter('Videos_13/QvZavg.avi');
MovAvi4.Quality = 100;
MovAvi4.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,qs*transpose(qvZavg(:,:,mov)),100); colorbar;
    caxis([0.1 5])
    
    axis([0 10*(Lx-dx) 0 10*(Ly-dy)])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-avg qv (x,y) contour in g/Kg, T = ', num2str(Tmov),' hours'])

    mov
    open(MovAvi4)
    F = getframe(figure(1));
    writeVideo(MovAvi4,F);
end

close(MovAvi4);
stop



%% %%%% Qr z-avg
display('Movie: QrZavg')
qrZavg = load('Data2D/QrZavg_14.dat');
display('file loaded')
movies = size(qrZavg,1)/nx/ny

NoStepsMov=Tfinal/(movies-1);
% NoStepsMov=Tfinal/(500-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;
end

qrZavg = reshape(qrZavg,nx,ny,movies);

MovAvi0 = VideoWriter('Videos_14/QrZavg.avi');
MovAvi0.Quality = 100;
MovAvi0.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,qs*transpose(qrZavg(:,:,mov)),100); colorbar;
    caxis([.1 5])
    
    axis([0 10*Lx 0 10*Ly])

    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-avg qr (x,y) contour in g/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi0)
	F = getframe(figure(1));
	writeVideo(MovAvi0,F);
end

close(MovAvi0);


%% %%%%% Qr z-max:
display('Movie: QrZMax')
qrZMax = load('Data2D/QrZMax_14.dat');
display('file loaded')
movies = size(qrZMax,1)/nx/ny

NoStepsMov=Tfinal/(movies-1);
% NoStepsMov=Tfinal/(500-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;
end

qrZMax = reshape(qrZMax,nx,ny,movies);

MovAvi1 = VideoWriter('Videos_14/QrZMax.avi');
MovAvi1.Quality = 100;
MovAvi1.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,qs*transpose(qrZMax(:,:,mov)),100); colorbar;

    caxis([0.1 5])
    %caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-max qr (x,y) contour in g/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi1)
    F = getframe(figure(1));
    writeVideo(MovAvi1,F);
end

close(MovAvi1);

%% %%%%%%% Qv z-max
display('Movie: QvZMax')
qvZMax = load('Data2D/QvZMax_14.dat');
display('file loaded')
movies = size(qvZMax,1)/nx/ny

NoStepsMov=Tfinal/(movies-1);
% NoStepsMov=Tfinal/(500-1);

time=zeros(movies,1);

for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;
end


qvZMax = reshape(qvZMax,nx,ny,movies);

MovAvi3 = VideoWriter('Videos_14/QvZMax.avi');
MovAvi3.Quality = 100;
MovAvi3.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,qs*transpose(qvZMax(:,:,mov)),100); colorbar;
    caxis([.1 5])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-max qv (x,y) contour in g/Kg, T = ', num2str(Tmov),' hours'])

    
	mov
	open(MovAvi3)
    F = getframe(figure(1));
    writeVideo(MovAvi3,F);
end

close(MovAvi3);


%% %%%%%%%% Qv z-avg
display('Movie: QvZavg')
qvZavg = load('Data2D/QvZavg_14.dat');
display('file loaded')
movies = size(qvZavg,1)/nx/ny

NoStepsMov=Tfinal/(movies-1);
% NoStepsMov=Tfinal/(500-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;
end

qvZavg = reshape(qvZavg,nx,ny,movies);

MovAvi4 = VideoWriter('Videos_14/QvZavg.avi');
MovAvi4.Quality = 100;
MovAvi4.FrameRate = 5;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,qs*transpose(qvZavg(:,:,mov)),100); colorbar;
    caxis([0.1 5])
    
    axis([0 10*(Lx-dx) 0 10*(Ly-dy)])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-avg qv (x,y) contour in g/Kg, T = ', num2str(Tmov),' hours'])

    mov
    open(MovAvi4)
    F = getframe(figure(1));
    writeVideo(MovAvi4,F);
end

close(MovAvi4);


stop


%% %%%%%%%%%% Theta at low levels
display('Movie: Theta at 5 km')
thetaZpt5Km = load('2D_ThetaZpt5km_30.dat');
display('file loaded')
movies = size(thetaZpt5Km,1)/nx/ny


NoStepsMov=Tfinal/(movies-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;
end


thetaZpt5Km = reshape(thetaZpt5Km,nx,ny,movies);

MovAvi6 = VideoWriter('ThetaZpt5Km.avi');
MovAvi6.Quality = 100;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,transpose(3*thetaZpt5Km(:,:,mov)),30); colorbar;
    caxis([-3 0.5])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['theta (x,y) contour at z=0.5 km in Kelvin, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi6)
    F = getframe(figure(1));
    writeVideo(MovAvi6,F);
end

close(MovAvi6);


%% %%%%%% U y-half
display('Movie: UYHalf')
UYHalf = load('2D_UYHalf_30.dat');
display('file loaded')
movies = size(UYHalf,1)/nx/m

NoStepsMov=Tfinal/(movies-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;
end

UYHalf = reshape(UYHalf,nx,m,movies);

MovAvi7 = VideoWriter('UYHalf.avi');
MovAvi7.Quality = 100;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*z,(epsilon^2)*transpose(UYHalf(:,:,mov)),30); colorbar;
    caxis([0 5*10^(-3)])
    caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['U (x,y) contour at y=Ly/2 in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi7)
    F = getframe(figure(1));
    writeVideo(MovAvi7,F);
end

close(MovAvi7);

%% %%%%%% V y-half
display('Movie: VYHalf')
VYHalf = load('2D_VYHalf_30.dat');
display('file loaded')
movies = size(VYHalf,1)/nx/m

NoStepsMov=Tfinal/(movies-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;
end

VYHalf = reshape(VYHalf,nx,m,movies);

MovAvi8 = VideoWriter('VYHalf.avi');
MovAvi8.Quality = 100;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*z,(epsilon^2)*transpose(VYHalf(:,:,mov)),30); colorbar;
    caxis([0 5*10^(-3)])
    caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['V (x,y) contour at y=Ly/2 in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi8)
    F = getframe(figure(1));
    writeVideo(MovAvi8,F);
end

close(MovAvi8);

%% %%%%%% W y-half
display('Movie: WYHalf')
WYHalf = load('2D_WYHalf_30.dat');
display('file loaded')
movies = size(WYHalf,1)/nx/m

NoStepsMov=Tfinal/(movies-1);

time=zeros(movies,1);
for mov=1:movies
    time(mov)=NoStepsMov*(mov-1)*15/60;
end

WYHalf = reshape(WYHalf,nx,m,movies);

MovAvi9 = VideoWriter('WYHalf.avi');
MovAvi9.Quality = 100;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*z,(epsilon^2)*transpose(WYHalf(:,:,mov)),30); colorbar;
    caxis([0 5*10^(-3)])
    caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['W (x,y) contour at y=Ly/2 in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi9)
    F = getframe(figure(1));
    writeVideo(MovAvi9,F);
end

close(MovAvi9);


stop

qr = load('3D_Qr_10.dat');
qr = reshape(qr,nx,ny,m);

qr_y = zeros(nx,m);
figure(1)
qr_y(:,:) = 0.5*(qr(:,ny/2 -1,:)+qr(:,ny/2+1,:));
contour(10*x,10*z,(epsilon^2)*transpose(qr_y(:,:)),100); colorbar;

figure(2)
qr_y(:,:) = 0.5*(qr(:,ny/2,:)+qr(:,ny/2,:));
contour(10*x,10*z,(epsilon^2)*transpose(qr_y(:,:)),100); colorbar;
stop

xslice=10*(Lx-dx);% 10*Lx/2; %0:10*Lx/128:10*Lx;%
yslice=10*(Ly-dy); %0:10*Ly/32:10*(Ly-dy); %10*Lx/2; %0:10*Ly/128:10*Ly; %
zslice=0:10*Lz/100:10*Lz; 

h=slice(10*xGrid,10*yGrid,10*zGrid,(epsilon^2)*qr,yslice,xslice,zslice); colorbar;
%caxis([.1*10^(-3) 2.5*10^(-3)])
view(78,54)
axis tight

alpha('color')
set(h,'EdgeColor','none','FaceColor','interp',...
 'FaceAlpha','interp')

stop

% 
% %%%%%%%%% Qr z-max
% display('Movie: QrZMax')
% qrZMax = load('2D_QrZMax_30.dat');
% qrZMax = reshape(qrZMax,nx,ny,movies);
% 
% MovAvi0 = VideoWriter('QrZMax.avi');
% MovAvi0.Quality = 100;
% for mov=1:movies
%     
%     Tmov=NoStepsMov*(mov-1);
%     Tmov=Tmov*15/60; %In real time (hours)
% 
%     clf
%     figure(1)
%     contour(10*x,10*y,(epsilon^2)*transpose(qrZMax(:,:,mov)),30); colorbar;
%     %caxis([0.3*10^(-3) 10^(-3)])
%     
%     axis([0 10*Lx 0 10*Ly])
%     axis on
%     xlabel('x in km');
%     ylabel('y in km');    
%     
%     title(['z-max qr (x,y) contour in Kg/Kg, T = ', num2str(Tmov),' hours'])
%     
% 	mov
% 	open(MovAvi0)
%     F = getframe(figure(1));
%     writeVideo(MovAvi0,F);
%     pause(pa)
% end
% 
% close(MovAvi0);
% 
% %%%%%%%%%% Qr z-avg
% display('Movie: QrZavg')
% qrZavg = load('2D_QrZavg_30.dat');
% qrZavg = reshape(qrZavg,nx,ny,movies);
% 
% MovAvi1 = VideoWriter('QrZavg.avi');
% MovAvi1.Quality = 100;
% for mov=1:movies
%     
%     Tmov=NoStepsMov*(mov-1);
%     Tmov=Tmov*15/60; %In real time (hours)
% 
%     clf
%     figure(1)
%     contour(10*x,10*y,(epsilon^2)*transpose(qrZavg(:,:,mov)),30); colorbar;
%     caxis([0.2*10^(-3) 0.7*10^(-3)])
%     
%     axis([0 10*(Lx-dx) 0 10*(Ly-dy)])
%     axis on
%     xlabel('x in km');
%     ylabel('y in km');    
%     
%     title(['z-avg qr (x,y) contour in Kg/Kg, T = ', num2str(Tmov),' hours'])
% 
%     mov
%     open(MovAvi1)
%     F = getframe(figure(1));
%     writeVideo(MovAvi1,F);
%     pause(pa)
% end
% 
% close(MovAvi1);

%% %%%%%%% Qr y-half

display('Movie: QrYHalf')
qrYHalf = load('2D_QvYHalf_30.dat');
%size(qrYHalf,1)/nx/m
qrYHalf = reshape(qrYHalf,nx,m,movies);

MovAvi2 = VideoWriter('QrYHalf.avi');
MovAvi2.Quality = 100;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*z,(epsilon^2)*transpose(qrYHalf(:,:,mov)),30); colorbar;
    %caxis([0.3*10^(-3) 10^(-3)])
    
    axis([0 10*Lx 0 10*Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['qr (x,y) contour at y=Ly/2 in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi2)
    F = getframe(figure(1));
    writeVideo(MovAvi2,F);
    %pause(pa)
end

close(MovAvi2);

stop

%% %%%% Qv z-avg
display('Movie: QvZavg')
qvZavg = load('2D_QvZavg_30.dat');
qvZavg = reshape(qvZavg,nx,ny,movies);

MovAvi3 = VideoWriter('QvZavg.avi');
MovAvi3.Quality = 100;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,(epsilon^2)*transpose(qvZavg(:,:,mov)),30); colorbar;
    caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-avg qv (x,y) contour in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi3)
	F = getframe(figure(1));
	writeVideo(MovAvi3,F);
    pause(pa)
end

close(MovAvi3);

%%%%%%% Qv z-max:
display('Movie: QvZMax')
qvZMax = load('2D_QvZMax_30.dat');
qvZMax = reshape(qvZMax,nx,ny,movies);

MovAvi4 = VideoWriter('QvZMax.avi');
MovAvi4.Quality = 100;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*y,(epsilon^2)*transpose(qvZMax(:,:,mov)),30); colorbar;
    %caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Ly])
    axis on
    xlabel('x in km');
    ylabel('y in km');    
    
    title(['z-max qv (x,y) contour in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi4)
    F = getframe(figure(1));
    writeVideo(MovAvi4,F);
    pause(pa)
end

close(MovAvi4);


%%%%%%%% Qv y-half
display('Movie: QvYHalf')
qvYHalf = load('2D_QvYHalf_30.dat');
qvYHalf = reshape(qvYHalf,nx,m,movies);

MovAvi5 = VideoWriter('QvYHalf.avi');
MovAvi5.Quality = 100;
for mov=1:movies
    
    Tmov=NoStepsMov*(mov-1);
    Tmov=Tmov*15/60; %In real time (hours)

    clf
    figure(1)
    contour(10*x,10*z,(epsilon^2)*transpose(qvYHalf(:,:,mov)),30); colorbar;
    caxis([0 5*10^(-3)])
    %caxis([-2*10^(-4) 6*10^(-4)])
    
    axis([0 10*Lx 0 10*Lz])
    axis on
    xlabel('x in km');
    ylabel('z in km');    
    
    title(['qv (x,y) contour at y=Ly/2 in Kg/Kg, T = ', num2str(Tmov),' hours'])
    
	mov
	open(MovAvi5)
    F = getframe(figure(1));
    writeVideo(MovAvi5,F);
    pause(pa)
end

close(MovAvi5);
stop

% %%%%%%%%%% Qr z-avg
% display('Movie: uYHalfHorzmean')
% 
% MovAvi2 = VideoWriter('uYHalfHorzmean.avi');
% MovAvi2.Quality = 100;
% for mov=1:movies
%     clf
%     figure(1)
%     %hold on
%     plot(uYHalfHmean(1:m-1,mov),z(1:m-1),'r .',ubg,z,'b -');
%         
% 	mov
% 	open(MovAvi2)
%     F = getframe(figure(1));
%     writeVideo(MovAvi2,F);
% end
% 
% close(MovAvi2);
% 
% 
% stop
% 
% 
% 
% 
% %%%%Steady squall line
% 
% MovAvi1 = avifile('QrYHalfSteady.avi');
% MovAvi1.fps = 1;
% for mov=it_0:it_1
%     
%     Tmov=NoStepsMov*(mov-1);
%     Tmov=Tmov*15/60; %In real time (hours)
% 
%     clf
%     figure(1)
%     contour(10*x,10*z,(epsilon^2)*transpose(qrYHalfSteady(:,:,mov)),30); colorbar;
%     caxis([0.3*10^(-3) 10^(-3)])
%     
%     axis([0 10*Lx 0 10*Lz])
%     axis on
%     xlabel('x in km');
%     ylabel('z in km');    
%     
%     title(['Steady qr (x,y) contour at y=Ly/2 in Kg/Kg, T = ', num2str(Tmov),' hours'])
%     
%     F = getframe(figure(1));
%     MovAvi1 = addframe(MovAvi1,F);
% end
% 
% MovAvi1 = close(MovAvi1);
% 
% %%%%%%%%%%

%Plots:

%%%%%%(x,t) contour of qv

figure(1)
clf
contour(10*x,time,transpose(qrYHalfZavg),60); colorbar;
hold on
plot(u_squall*time,time+5,'r *')
title('(x,t) contour qr at y=Ly/2, z-averaged')
print('qrYHalfZavg_xtcontour','-depsc',figure(1))

%%%%%%(x,t) contour of qv

figure(1)
clf
contour(10*x,time,transpose(qvYHalfZavg),60); colorbar;
hold on
plot(u_squall*time,time+5,'r *')
title('(x,t) contour qv at y=Ly/2, z-averaged')
print('qvYHalfZavg_xtcontour','-depsc',figure(1))

%%(x,t) contour of theta at low levels

figure(1)
clf
contour(10*x,time,transpose(thetaZpt5KmYavg),60); colorbar;
hold on
plot(u_squall*time,time+5,'r *')
title('(x,t) contour \theta at y=Ly/2, z=')
print('ThetaZpt5KmYavg_xtcontour','-depsc',figure(1))

%Steady squall line

figure(1)
clf
contour(10*x,10*z,(epsilon^2)*transpose(qrYHalfSteady_intT(:,:)),60); colorbar;
print('QrYHalfSteadyTint','-depsc',figure(1))

%%%All variables

mov=9; %round(0.8*60/(15*NoStepsMov)+1);
%mov=40;
Tmov=NoStepsMov*(mov-1)*15/60;

figure(1)
subplot(3,3,1)
contour(10*x,10*z,(epsilon^2)*transpose(qrYHalf(:,:,mov)),30); colorbar('FontSize',8);
axis([0 10*Lx 0 10*Lz])
axis on
xlabel('x in km','FontSize',8);
ylabel('z in km','FontSize',8);
set(gca,'FontSize',8);
title(['q_r in kg/kg, T = ', num2str(Tmov),' hrs'],'FontSize',8)

c=colorbar;
x1=get(gca,'position');
x2=get(c,'Position');
x2(3)=0.01;
set(c,'Position',x2)
set(gca,'position',x1)
set(gca,'FontSize',8)

subplot(3,3,2)
contour(10*x,10*z,(epsilon^2)*transpose(qvYHalf(:,:,mov)),30); colorbar('FontSize',8);
axis([0 10*Lx 0 10*Lz])
axis on
xlabel('x in km','FontSize',8);
ylabel('z in km','FontSize',8);
set(gca,'FontSize',8);
title(['q_v in kg/kg, T = ', num2str(Tmov),' hrs'],'FontSize',8)

c=colorbar;
x1=get(gca,'position');
x2=get(c,'Position');
x2(3)=0.01;
set(c,'Position',x2)
set(gca,'position',x1)
set(gca,'FontSize',8)

subplot(3,3,3)
contour(10*x,10*y,3*transpose(thetaZpt5Km(:,:,mov)),30); colorbar('FontSize',8);
axis([0 10*Lx 0 10*Lz])
axis on
xlabel('x in km','FontSize',8);
ylabel('z in km','FontSize',8);
set(gca,'FontSize',8);
title(['\theta in Kelvin, T = ', num2str(Tmov),' hrs'],'FontSize',8)

c=colorbar;
x1=get(gca,'position');
x2=get(c,'Position');
x2(3)=0.01;
set(c,'Position',x2)
set(gca,'position',x1)
set(gca,'FontSize',8)

subplot(3,3,4)
contour(10*x,10*z,(10/9)*transpose(uYHalf(:,:,mov)),30); colorbar('FontSize',8);
axis([0 10*Lx 0 10*Lz])
axis on
xlabel('x in km','FontSize',8);
ylabel('z in km','FontSize',8);
set(gca,'FontSize',8);
title(['u in m/s, T = ', num2str(Tmov),' hrs'],'FontSize',8)

c=colorbar;
x1=get(gca,'position');
x2=get(c,'Position');
x2(3)=0.01;
set(c,'Position',x2)
set(gca,'position',x1)
set(gca,'FontSize',8)

subplot(3,3,5)
contour(10*x,10*z,(10/9)*transpose(vYHalf(:,:,mov)),30); colorbar('FontSize',8);
axis([0 10*Lx 0 10*Lz])
axis on
xlabel('x in km','FontSize',8);
ylabel('z in km','FontSize',8);
set(gca,'FontSize',8);
title(['v in m/s, T = ', num2str(Tmov),' hrs'],'FontSize',8)

c=colorbar;
x1=get(gca,'position');
x2=get(c,'Position');
x2(3)=0.01;
set(c,'Position',x2)
set(gca,'position',x1)
set(gca,'FontSize',8)

subplot(3,3,6)
contour(10*x,10*z,(10/9)*transpose(wYHalf(:,:,mov)),30); colorbar('FontSize',8);
axis([0 10*Lx 0 10*Lz])
axis on
xlabel('x in km','FontSize',8);
ylabel('z in km','FontSize',8);
set(gca,'FontSize',8);
title(['w in m/s, T = ', num2str(Tmov),' hrs'],'FontSize',8)

c=colorbar;
x1=get(gca,'position');
x2=get(c,'Position');
x2(3)=0.01;
set(c,'Position',x2)
set(gca,'position',x1)
set(gca,'FontSize',8)

print('AllVariablesT24','-depsc',figure(1))

% mov=7;
% xslice=10*(Lx-dx);% 10*Lx/2; %0:10*Lx/128:10*Lx;%
% yslice=10*(Ly-dy); %0:10*Ly/32:10*(Ly-dy); %10*Lx/2; %0:10*Ly/128:10*Ly; %
% zslice=0:10*Lz/20:10*Lz; 
% 
% qr_mov=qr(:,:,:,mov);
% 
% h=slice(10*xGrid,10*yGrid,10*zGrid,(epsilon^2)*qr_mov,yslice,xslice,zslice);
% %caxis([.1*10^(-3) 2.5*10^(-3)])
% view(78,54)
% axis tight
% 
% alpha('color')
% set(h,'EdgeColor','none','FaceColor','interp',...
%  'FaceAlpha','interp')
% 
% stop

%%variables:
% 
% 
%  
% uYHalf = load('2D_UYHalf_30.dat');
% uYHalf = reshape(uYHalf,nx,m,movies);
% 
% vYHalf = load('2D_VYHalf_30.dat');
% vYHalf = reshape(vYHalf,nx,m,movies);
% 
% wYHalf = load('2D_WYHalf_30.dat');
% wYHalf = reshape(wYHalf,nx,m,movies);
% 
% 
% 
% qrYHalfZavg=zeros(nx,movies);
% for mov=1:movies
%     for ix=1:nx
%         qrYHalfZavg(ix,mov)=sum(qrYHalf(ix,:,mov))/m;
%     end
% end
% 
% qvYHalfZavg=zeros(nx,movies);
% for mov=1:movies
%     for ix=1:nx
%         qvYHalfZavg(ix,mov)=sum(qvYHalf(ix,:,mov))/m;
%     end
% end
% 
% thetaZpt5KmYavg=zeros(nx,movies);
% for mov=1:movies
%     for ix=1:nx
%         thetaZpt5KmYavg(ix,mov)=sum(thetaZpt5Km(ix,:,mov))/ny;
%     end
% end
% 
% % %Steady squall line
% it_0=max(round(0/24*100),1);
% it_1=round(2/24*100);
% 
% qrYHalfSteady=zeros(nx,m,movies);
% for mov=1:movies
%     Tmov=NoStepsMov*(mov-1);
%     Tmov=Tmov*15/60; %In real time (hours)
%     
%     for ix=1:nx
%         xi=dx*(ix-0.5);
%         xi_squall=mod(xi+u_squall*Tmov/10,12.8); %1/10 to make it non-dml
%         ix_squall=max(mod(round(xi_squall/dx+0.5),nx),1);
%         for iz=1:m
%             qrYHalfSteady(ix,iz,mov)=qrYHalf(ix_squall,iz,mov);
%         end
%     end
% end
% 
% qrYHalfSteady_intT=zeros(nx,m);
% for ix=1:nx
%     for iz=1:m
%         qrYHalfSteady_intT(ix,iz)=sum(qrYHalfSteady(ix,iz,it_0:it_1))/(it_1-it_0+1);
%     end
% end
% 
% %Horizontal mean
% uYHalfHmean=zeros(m,movies);
% for mov=1:movies
%     for k=1:m
%         uYHalfHmean(k,mov)=sum(uYHalf(:,k,mov))/nx;
%     end
% end
