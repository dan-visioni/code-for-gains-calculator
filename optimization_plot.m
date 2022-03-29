clear all
dir_new = '/Volumes/D_Visioni2/CESM2_MA/sys_ID/';
dir_ma = '/Volumes/D_Visioni2/CESM2_MA/ssp245ma/';
sspnam = '_b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.001_zm_monthly.nc';

vol2 = '/Volumes/D_Visioni2/GeoMIP/ssp245/CESM2-WACCM/';
vol3 = '/Volumes/D_Visioni2/GeoMIP/ssp245/UKESM1-0-LL/';

%% reading values for CESM

loc = {'30N','15N','0N','15S','30S'};
longn = '_INJANN';

vart = 'TREFHT';
T_base_new = ncread([dir_ma vart sspnam],'TREFHT'); nb = 20*12+1;
T_base_mm = T_base_new(:,nb:nb+119);
T_base_new = mean(T_base_mm(:,49:120),2);

for i=1:length(loc)
    AOD_new{i} = ncread([dir_new 'AODVISstdn' longn loc{i} '_12Tg.001_zm_monthly.nc'],'AODVISstdn');
    T_new{i} = ncread([dir_new 'TREFHT' longn loc{i} '_12Tg.001_zm_monthly.nc'],'TREFHT');
end

lat = ncread([dir_new 'AODVISstdn' longn loc{i} '_12Tg.001_zm_monthly.nc'],'lat');
%% reading values for GISS bulk

foldG = '/Volumes/D_Visioni2/GISS/';
locG = {'c4bk','c7bk','c6bk','c1bk'};

for i=1:length(loc)
    AOD_G{i} = ncread([foldG 'GISS_corr_' loc{i} '_zm.nc'],'AOD550nm');
    T_G{i} = ncread([foldG 'GISS_corr_' loc{i} '_zm.nc'],'surfT');
end

AOD_base_G = ncread([foldG 'GISS_ssp245_zm.nc'],'AOD550nm'); nb = 15*12+1;
AOD_base_Gmm = AOD_base_G(:,nb:nb+119);
AOD_base_Gnew = mean(AOD_base_Gmm(:,37:120),2);

T_base_G = ncread([foldG 'GISS_ssp245_zm.nc'],'surfT');nb = 15*12+1;
T_base_Gmm = T_base_G(:,nb:nb+119);
T_base_Gnew = mean(T_base_Gmm(:,49:120),2);

latG =-89:2:89;

%% reading values for GISS modal

foldGm = '/Volumes/D_Visioni2/GISSmodal/';
locG = {'c4bk','c7bk','c6bk','c1bk'};

for i=1:length(loc)
    AOD_Gm{i} = ncread([foldGm 'GISS_modal_corr_' loc{i} '_zm.nc'],'AOD550nm');
    T_Gm{i} = ncread([foldGm 'GISS_modal_corr_' loc{i} '_zm.nc'],'surfT');
end

AOD_base_Gm = ncread([foldGm 'GISS_modal_ssp245_zm.nc'],'AOD550nm'); nb = 20*12+1;
AOD_base_Gmmm = AOD_base_Gm(:,nb:nb+119);
AOD_base_Gnewm = mean(AOD_base_Gmmm(:,37:120),2);

T_base_Gm = ncread([foldGm 'GISS_modal_ssp245_zm.nc'],'surfT');nb = 20*12+1;
T_base_Gmmm = T_base_Gm(:,nb:nb+119);
T_base_Gnewm = mean(T_base_Gmmm(:,49:120),2);


%% reading values for UKESM

foldG = '/Volumes/D_Visioni2/UKESM/';
loc2 = {'30N','15N','0N2','15S','30S'};
for i=1:length(loc)
    AOD_U{i} = ncread([foldG 'UKESM_' loc{i} '_zm.nc'],'AOD550nm');
    T_U{i} = ncread([foldG 'UKESM_' loc{i} '_zm.nc'],'surfT');
end

AOD_base_U = ncread([foldG 'UKESM_ssp245_zm.nc'],'AOD550nm'); nb = 20*12+1;
latU = ncread([foldG 'UKESM_ssp245_zm.nc'],'latitude');
AOD_base_Umm = AOD_base_U(:,nb:nb+119);
AOD_base_Unew = mean(AOD_base_Umm(:,37:120),2);

T_base_U = ncread([foldG 'UKESM_ssp245_zm.nc'],'surfT');
nb = 20*12+1;
T_base_Umm = T_base_U(:,nb:nb+119);
T_base_Unew = mean(T_base_Umm(:,49:120),2);

%% selecting colors (uses Brewermap, can skip)

fc = brewermap(4,'Pastel1');
fc = brighten(fc,-.8);
fc1 = brighten(fc,.8);

ww = cos(lat/180*pi); sw = sum(ww);
wwG = cos(latG/180*pi); swG = sum(wwG);
wwU = cos(latU/180*pi); swU = sum(wwU);



%% averages over selected number of years

years = 7; 
dy = 10-years; dm = dy*12+1; dme = 120;

for i=1:5
    CESM{i} = mean(AOD_new{i}(:,dm:dme),2);
    GISS{i} = mean(AOD_G{i}(:,dm:dme),2)-mean(AOD_base_Gmm(:,dm:dme),2);
    GISSm{i} = mean(AOD_Gm{i}(:,dm:dme),2)-mean(AOD_base_Gmmm(:,dm:dme),2);
    UKESM{i} = mean(AOD_U{i}(:,dm:dme),2)-mean(AOD_base_Umm(:,dm:dme),2);
end

%% determines solution to least-square problem

%% defines legendre polynomials
L0 = ones(length(lat),1);
L1=sind(lat);
L1N = L0+sind(lat);
L1S = L0-sind(lat);
L2 = L0+(1.5*L1.^2-.5);

awt = cosd(lat); 

%% finds solution for L0 using 15N+15S

psi=[CESM{2},CESM{4}]/12;
[q,fiterr0{1,1}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L0,-eye(2),zeros(2,1));
inj0(1,1) = q(1); inj0(1,2) = q(2);

%% finds solution for L1N using 15N+30N

psi=[CESM{1},CESM{2}]/12;
[q,fiterr0{2,1}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L1N,-eye(2),zeros(2,1));
inj1N(1,1) = q(1); inj1N(1,2) = q(2);

%% finds solution for L1S using 15S+30S

psi=[CESM{4},CESM{5}]/12;
[q,fiterr0{3,1}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L1S,-eye(2),zeros(2,1));
inj1S(1,1) = q(1); inj1S(1,2) = q(2);

%% finds solution for L2 using 30N+30S

psi=[CESM{1},CESM{5}]/12;
[q,fiterr0{4,1}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L2,-eye(2),zeros(2,1));
inj2(1,1) = q(1); inj2(1,2) = q(2);

%% repeats for UKESM

latU=double(latU);
awt = cosd(latU); 
L0 = ones(length(latU),1);
L1=sind(latU);
L1N = L0+sind(latU);
L1S = L0-sind(latU);
L2 = L0+(1.5*L1.^2-.5);

psi=[UKESM{2},UKESM{4}]/12;
[q,fiterr0{1,2}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L0,-eye(2),zeros(2,1));
inj0(2,1) = q(1); inj0(2,2) = q(2);

psi=[UKESM{1},UKESM{2}]/12;
[q,fiterr0{2,2}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L1N,-eye(2),zeros(2,1));
inj1N(2,1) = q(1); inj1N(2,2) = q(2);

psi=[UKESM{4},UKESM{5}]/12;
[q,fiterr0{3,2}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L1S,-eye(2),zeros(2,1));
inj1S(2,1) = q(1); inj1S(2,2) = q(2);

psi=[UKESM{1},UKESM{5}]/12;
[q,fiterr0{4,2}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L2,-eye(2),zeros(2,1));
inj2(2,1) = q(1); inj2(2,2) = q(2);

%% repeats for GISS bulk

latG=double(latG');
awt = cosd(latG); 
L0 = ones(length(latG),1);
L1=sind(latG);
L1N = L0+sind(latG);
L1S = L0-sind(latG);
L2 = L0+(1.5*L1.^2-.5);

psi=double([GISS{2},GISS{4}]/12);
[q,fiterr0{1,3}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L0,-eye(2),zeros(2,1));
inj0(3,1) = q(1); inj0(3,2) = q(2);

psi=double([GISS{1},GISS{2}]/12);
[q,fiterr0{2,3}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L1N,-eye(2),zeros(2,1));
inj1N(3,1) = q(1); inj1N(3,2) = q(2);

psi=double([GISS{4},GISS{5}]/12);
[q,fiterr0{3,3}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L1S,-eye(2),zeros(2,1));
inj1S(3,1) = q(1); inj1S(3,2) = q(2);

psi=double([GISS{1},GISS{5}]/12);
[q,fiterr0{4,3}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L2,-eye(2),zeros(2,1));
inj2(3,1) = q(1); inj2(3,2) = q(2);

%% repeats for GISS modal

psi=double([GISSm{2},GISSm{4}]/12);
[q,fiterr0{1,4}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L0,-eye(2),zeros(2,1));
inj0(4,1) = q(1); inj0(4,2) = q(2);

psi=double([GISSm{1},GISSm{2}]/12);
[q,fiterr0{2,4}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L1N,-eye(2),zeros(2,1));
inj1N(4,1) = q(1); inj1N(4,2) = q(2);

psi=double([GISSm{4},GISSm{5}]/12);
[q,fiterr0{3,4}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L1S,-eye(2),zeros(2,1));
inj1S(4,1) = q(1); inj1S(4,2) = q(2);

psi=double([GISSm{1},GISSm{5}]/12);
[q,fiterr0{4,4}]=lsqlin((sqrt(awt)*ones(1,2)).*psi,sqrt(awt).*L2,-eye(2),zeros(2,1));
inj2(4,1) = q(1); inj2(4,2) = q(2);



%%
mods = {'CESM','UKESM','GISS','GISSm'};

%% getting the AOD pattern resulting from the optimization

opt_0_CESM = (CESM{2}*inj0(1,1)+CESM{4}*inj0(1,2))/12;
opt_0_GISS = (GISS{2}*inj0(3,1)+GISS{4}*inj0(3,2))/12;
opt_0_GISSm = (GISSm{2}*inj0(4,1)+GISSm{4}*inj0(4,2))/12;
opt_0_UKESM = (UKESM{2}*inj0(2,1)+UKESM{4}*inj0(2,2))/12;

opt_1N_CESM = (CESM{1}*inj1N(1,1)+CESM{2}*inj1N(1,2))/12;
opt_1N_GISS = (GISS{1}*inj1N(3,1)+GISS{2}*inj1N(3,2))/12;
opt_1N_GISSm = (GISSm{1}*inj1N(4,1)+GISSm{2}*inj1N(4,2))/12;
opt_1N_UKESM = (UKESM{1}*inj1N(2,1)+UKESM{2}*inj1N(2,2))/12;

opt_1S_CESM = (CESM{4}*inj1S(1,1)+CESM{5}*inj1S(1,2))/12;
opt_1S_GISS = (GISS{4}*inj1S(3,1)+GISS{5}*inj1S(3,2))/12;
opt_1S_GISSm = (GISSm{4}*inj1S(4,1)+GISSm{5}*inj1S(4,2))/12;
opt_1S_UKESM = (UKESM{4}*inj1S(2,1)+UKESM{5}*inj1S(2,2))/12;

opt_2_CESM = (CESM{1}*inj2(1,1)+CESM{5}*inj2(1,2))/12;
opt_2_GISS = (GISS{1}*inj2(3,1)+GISS{5}*inj2(3,2))/12;
opt_2_GISSm = (GISSm{1}*inj2(4,1)+GISSm{5}*inj2(4,2))/12;
opt_2_UKESM = (UKESM{1}*inj2(2,1)+UKESM{5}*inj2(2,2))/12;

L0 = ones(length(lat),1);
L1=sind(lat);
L1N = L0+sind(lat);
L1S = L0-sind(lat);
L2 = L0+(1.5*L1.^2-.5);


figure(4)
set(gcf, 'Position',  [200, 200, 1200,1000])

subplot(8,2,[1 5])
plot(lat,L0,'--k','Linewidth',2)
hold on
box on
plot(lat,opt_0_CESM,'Color',fc(1,:),'Linewidth',4)
plot(latU,opt_0_UKESM,'Color',fc(2,:),'Linewidth',4)
plot(latG,opt_0_GISS,'Color',fc(3,:),'Linewidth',4)
plot(latG,opt_0_GISSm,'Color',fc(4,:),'Linewidth',4)
axis([-90 90 0 2])
plot([-15 -15],[0 2],':k','Linewidth',.5)
plot([15 15],[0 2],':k','Linewidth',.5)
set(gca,'Linewidth',2,'Fontsize',18,'XTick',-90:30:90)
ylabel('AOD')
fS = '%.0f';
fS2 = '%.2f';
title('a) L0 optimization with 15N and 15S injections')

subplot(8,2,[7])
hold on
box on
for i=1:4
text(.2,3.6-i*.75,mods{i},'Color',fc(i,:),'Fontsize',16)
text(2,3.6-i*.75,['15N inj.=' num2str(inj0(i,1),fS)],'Color',fc(i,:),'Fontsize',16)
text(1,3.6-i*.75,['15S inj.=' num2str(inj0(i,2),fS)],'Color',fc(i,:),'Fontsize',16)
text(3,3.6-i*.75,['Residual=' num2str(fiterr0{i,1},fS2)],'Color',fc(i,:),'Fontsize',16)
end

axis([.0 4.5 0 3.5])
set(gca,'Linewidth',2,'Fontsize',18,'XTick',10:20,'YTick',10:20)

subplot(8,2,[2 6])
plot(lat,L0+L1,'--k','Linewidth',2)
hold on
box on
plot(lat,opt_1N_CESM,'Color',fc(1,:),'Linewidth',4)
plot(latU,opt_1N_UKESM,'Color',fc(2,:),'Linewidth',4)
plot(latG,opt_1N_GISS,'Color',fc(3,:),'Linewidth',4)
plot(latG,opt_1N_GISSm,'Color',fc(4,:),'Linewidth',4)
axis([-90 90 0 2])
set(gca,'Linewidth',2,'Fontsize',18,'XTick',-90:30:90)
plot([30 30],[0 2],':k','Linewidth',.5)
plot([15 15],[0 2],':k','Linewidth',.5)

ylabel('AOD')
title('b) L1N optimization with 15N and 30N injections')

subplot(8,2,[8])
hold on
box on
for i=1:4
text(.2,3.6-i*.75,mods{i},'Color',fc(i,:),'Fontsize',16)
text(2,3.6-i*.75,['30N inj.=' num2str(inj1N(i,1),fS)],'Color',fc(i,:),'Fontsize',16)
text(1,3.6-i*.75,['15N inj.=' num2str(inj1N(i,2),fS)],'Color',fc(i,:),'Fontsize',16)
text(3,3.6-i*.75,['Residual=' num2str(fiterr0{i,2},fS2)],'Color',fc(i,:),'Fontsize',16)
end
axis([.0 4.5 0 3.5])
set(gca,'Linewidth',2,'Fontsize',18,'XTick',10:20,'YTick',10:20)


subplot(8,2,[9 13])
plot(lat,L0-L1,'--k','Linewidth',2)
hold on
box on
plot(lat,opt_1S_CESM,'Color',fc(1,:),'Linewidth',4)
plot(latU,opt_1S_UKESM,'Color',fc(2,:),'Linewidth',4)
plot(latG,opt_1S_GISS,'Color',fc(3,:),'Linewidth',4)
plot(latG,opt_1S_GISSm,'Color',fc(4,:),'Linewidth',4)
axis([-90 90 0 2])
set(gca,'Linewidth',2,'Fontsize',18,'XTick',-90:30:90)
plot([-30 -30],[0 2],':k','Linewidth',.5)
plot([-15 -15],[0 2],':k','Linewidth',.5)

ylabel('AOD')
xlabel('Latitude')
title('c) L1S optimization with 15S and 30S injections')

subplot(8,2,[15])
hold on
box on
for i=1:4
text(.2,3.6-i*.75,mods{i},'Color',fc(i,:),'Fontsize',16)
text(2,3.6-i*.75,['15S inj.=' num2str(inj1S(i,1),fS)],'Color',fc(i,:),'Fontsize',16)
text(1,3.6-i*.75,['30S inj.=' num2str(inj1S(i,2),fS)],'Color',fc(i,:),'Fontsize',16)
text(3,3.6-i*.75,['Residual=' num2str(fiterr0{i,3},fS2)],'Color',fc(i,:),'Fontsize',16)
end
axis([.0 4.5 0 3.5])
set(gca,'Linewidth',2,'Fontsize',18,'XTick',10:20,'YTick',10:20)


subplot(8,2,[10 14])
plot(lat,L2,'--k','Linewidth',2)
hold on
box on
plot(lat,opt_2_CESM,'Color',fc(1,:),'Linewidth',4)
plot(latU,opt_2_UKESM,'Color',fc(2,:),'Linewidth',4)
plot(latG,opt_2_GISS,'Color',fc(3,:),'Linewidth',4)
plot(latG,opt_2_GISSm,'Color',fc(4,:),'Linewidth',4)
axis([-90 90 0 2])
set(gca,'Linewidth',2,'Fontsize',18,'XTick',-90:30:90)
plot([-30 -30],[0 2],':k','Linewidth',.5)
plot([30 30],[0 2],':k','Linewidth',.5)

ylabel('AOD')
xlabel('Latitude')
title('d) L2 optimization with 30S and 30N injections')

subplot(8,2,[16])
hold on
box on
for i=1:4
text(.2,3.6-i*.75,mods{i},'Color',fc(i,:),'Fontsize',16)
text(2,3.6-i*.75,['30N inj.=' num2str(inj2(i,1),fS)],'Color',fc(i,:),'Fontsize',16)
text(1,3.6-i*.75,['30S inj.=' num2str(inj2(i,2),fS)],'Color',fc(i,:),'Fontsize',16)
text(3,3.6-i*.75,['Residual=' num2str(fiterr0{i,4},fS2)],'Color',fc(i,:),'Fontsize',16)
end
axis([.0 4.5 0 3.5])
set(gca,'Linewidth',2,'Fontsize',18,'XTick',10:20,'YTick',10:20)


set(gcf,'renderer','painters')
print(gcf,'-depsc2',['optimization_figure.eps'])