%% Walker's control gain calculator

% Walker Raymond Lee, Cornell University (wl644@cornell.edu)
% Last updated: March 2022

% This script takes in "matrix runs" of injections at 30N, 15N, 15S, and
% 30S and a "background" run and uses the data to compute feedforward and
% feedback gains for a "GLENS"-type SAI simulation to control T0, T1, and
% T2.

%% User Input: simulation names and variable names

% in this section, edit the names of .nc files to read and the names of
% which variables to pull.

% names of surface temperature, AOD, and latitude variables in your model
temperature_variable = 'surfT';
AOD_variable = 'AOD550nm';
lat_variable = 'lat';

% is the temperature in Celsius (1) as opposed to Kelvin (0)?
temperature_in_celsius = 1;

% is the latitude in radians (1) as opposed to degrees (0)?
latitude_in_radians = 0;

% what calendar does your model use? 360 day or 365 day?
calendar = '365 day';

% names of .nc files to read that contain your model output
background_temperature_filename = 'GISS_modal_ssp245_zm.nc';
inj30N_temperature_filename = 'GISS_modal_4Tg_30N_zm.nc';
inj15N_temperature_filename = 'GISS_modal_4Tg_15N_zm.nc';
inj15S_temperature_filename = 'GISS_modal_4Tg_15S_zm.nc';
inj30S_temperature_filename = 'GISS_modal_4Tg_30S_zm.nc';

inj30N_AOD_filename = 'GISS_modal_4Tg_30N_zm.nc';
inj15N_AOD_filename = 'GISS_modal_4Tg_15N_zm.nc';
inj15S_AOD_filename = 'GISS_modal_4Tg_15S_zm.nc';
inj30S_AOD_filename = 'GISS_modal_4Tg_30S_zm.nc';

% injection rate of your matrix runs, in Tg/yr
injection_rate = 4;

% what years do your background and matrix simulations span?
background_yrs = 2015:2064;
matrix_yrs = 2035:2044;

% what period in your background run are your "target" temperatures from?
reference_start = 2020;
reference_end = 2039;

%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW THIS LINE! %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Import data

% this section imports temperature and AOD data. It also imports the
% latitude grid and computes the first three Legendre polynomials as a
% function of sin(latitude).

% temperature data
backgd_T = double(ncread(background_temperature_filename,temperature_variable));
inj30N_T = double(ncread(inj30N_temperature_filename,temperature_variable));
inj15N_T = double(ncread(inj15N_temperature_filename,temperature_variable));
inj15S_T = double(ncread(inj15S_temperature_filename,temperature_variable));
inj30S_T = double(ncread(inj30S_temperature_filename,temperature_variable));

% convert to K, if necessary
if temperature_in_celsius == 1
    backgd_T = backgd_T + 273.15;
    inj30N_T = inj30N_T + 273.15;
    inj15N_T = inj15N_T + 273.15;
    inj15S_T = inj15S_T + 273.15;
    inj30S_T = inj30S_T + 273.15;
end

% aerosol optical depth data
inj30N_AOD = double(ncread(inj30N_AOD_filename,AOD_variable));
inj15N_AOD = double(ncread(inj15N_AOD_filename,AOD_variable));
inj15S_AOD = double(ncread(inj15S_AOD_filename,AOD_variable));
inj30S_AOD = double(ncread(inj30S_AOD_filename,AOD_variable));

% latitude and Legendre polynomials
lat = double(ncread(background_temperature_filename,lat_variable));

if latitude_in_radians == 1
    x = sin(lat);
else
    x = sind(lat);
end

L0 = 1;
L1 = x;
L2 = 1/2*(3*x.^2-1);

%% Compute T0/T1/T2 and l0/l1/l2

% this section computes temperature and AOD mapped onto the 0th-, 1st-, and
% 2nd-order Legendre polynomials.

backgd_T0 = squeeze(sum(backgd_T.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
backgd_T1 = squeeze(sum(backgd_T.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
backgd_T2 = squeeze(sum(backgd_T.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));

inj30N_T0 = squeeze(sum(inj30N_T.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
inj30N_T1 = squeeze(sum(inj30N_T.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
inj30N_T2 = squeeze(sum(inj30N_T.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));
inj30N_l0 = squeeze(sum(inj30N_AOD.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
inj30N_l1 = squeeze(sum(inj30N_AOD.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
inj30N_l2 = squeeze(sum(inj30N_AOD.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));

inj15N_T0 = squeeze(sum(inj15N_T.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
inj15N_T1 = squeeze(sum(inj15N_T.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
inj15N_T2 = squeeze(sum(inj15N_T.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));
inj15N_l0 = squeeze(sum(inj15N_AOD.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
inj15N_l1 = squeeze(sum(inj15N_AOD.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
inj15N_l2 = squeeze(sum(inj15N_AOD.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));

inj15S_T0 = squeeze(sum(inj15S_T.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
inj15S_T1 = squeeze(sum(inj15S_T.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
inj15S_T2 = squeeze(sum(inj15S_T.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));
inj15S_l0 = squeeze(sum(inj15S_AOD.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
inj15S_l1 = squeeze(sum(inj15S_AOD.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
inj15S_l2 = squeeze(sum(inj15S_AOD.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));

inj30S_T0 = squeeze(sum(inj30S_T.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
inj30S_T1 = squeeze(sum(inj30S_T.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
inj30S_T2 = squeeze(sum(inj30S_T.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));
inj30S_l0 = squeeze(sum(inj30S_AOD.*L0.*cosd(lat))/sum(cosd(lat).*L0.^2));
inj30S_l1 = squeeze(sum(inj30S_AOD.*L1.*cosd(lat))/sum(cosd(lat).*L1.^2));
inj30S_l2 = squeeze(sum(inj30S_AOD.*L2.*cosd(lat))/sum(cosd(lat).*L2.^2));

%% Annual means

% this section uses the "annual average" function to compute the annual
% means of the temperature and AOD variables.

backgd_T0 = annual_average(backgd_T0,[],[],calendar);
backgd_T1 = annual_average(backgd_T1,[],[],calendar);
backgd_T2 = annual_average(backgd_T2,[],[],calendar);

inj30N_T0 = annual_average(inj30N_T0,[],[],calendar);
inj30N_T1 = annual_average(inj30N_T1,[],[],calendar);
inj30N_T2 = annual_average(inj30N_T2,[],[],calendar);
inj30N_l0 = annual_average(inj30N_l0,[],[],calendar);
inj30N_l1 = annual_average(inj30N_l1,[],[],calendar);
inj30N_l2 = annual_average(inj30N_l2,[],[],calendar);

inj15N_T0 = annual_average(inj15N_T0,[],[],calendar);
inj15N_T1 = annual_average(inj15N_T1,[],[],calendar);
inj15N_T2 = annual_average(inj15N_T2,[],[],calendar);
inj15N_l0 = annual_average(inj15N_l0,[],[],calendar);
inj15N_l1 = annual_average(inj15N_l1,[],[],calendar);
inj15N_l2 = annual_average(inj15N_l2,[],[],calendar);

inj15S_T0 = annual_average(inj15S_T0,[],[],calendar);
inj15S_T1 = annual_average(inj15S_T1,[],[],calendar);
inj15S_T2 = annual_average(inj15S_T2,[],[],calendar);
inj15S_l0 = annual_average(inj15S_l0,[],[],calendar);
inj15S_l1 = annual_average(inj15S_l1,[],[],calendar);
inj15S_l2 = annual_average(inj15S_l2,[],[],calendar);

inj30S_T0 = annual_average(inj30S_T0,[],[],calendar);
inj30S_T1 = annual_average(inj30S_T1,[],[],calendar);
inj30S_T2 = annual_average(inj30S_T2,[],[],calendar);
inj30S_l0 = annual_average(inj30S_l0,[],[],calendar);
inj30S_l1 = annual_average(inj30S_l1,[],[],calendar);
inj30S_l2 = annual_average(inj30S_l2,[],[],calendar);

%% test figures

% this section produces a "sanity check" figure showing the l0, l1, and l2
% produced by each of the matrix runs.

figure;
set(gcf,'units','inches','position',[0,0,15,5])
left = 0.07; center = 0.40; right = 0.73;
bottom = 0.15; width = 0.25; height = 0.75;
cb_ofst = 0.02; cb_btm = 0.05; cb_wd = 0.02; cb_ht = 0.75;
pos1 = [left   bottom width height];
pos2 = [center bottom width height];
pos3 = [right  bottom width height];

subplot('Position',pos1)
hold on; grid on
plot(matrix_yrs,inj30N_l0,'b')
plot(matrix_yrs,inj15N_l0,'Color',[0.2,0.2,0.8])
plot(matrix_yrs,inj15S_l0,'Color',[0.8,0.2,0.2])
plot(matrix_yrs,inj30S_l0,'r')
title('$\ell_0$','Interpreter','LaTeX')
legend('30N','15N','15S','30S')
set(gca,'FontSize',15)

subplot('Position',pos2)
hold on; grid on
plot(matrix_yrs,inj30N_l1,'b')
plot(matrix_yrs,inj15N_l1,'Color',[0.2,0.2,0.8])
plot(matrix_yrs,inj15S_l1,'Color',[0.8,0.2,0.2])
plot(matrix_yrs,inj30S_l1,'r')
title('$\ell_1$','Interpreter','LaTeX')
legend('30N','15N','15S','30S')
set(gca,'FontSize',15)

subplot('Position',pos3)
hold on; grid on
plot(matrix_yrs,inj30N_l2,'b')
plot(matrix_yrs,inj15N_l2,'Color',[0.2,0.2,0.8])
plot(matrix_yrs,inj15S_l2,'Color',[0.8,0.2,0.2])
plot(matrix_yrs,inj30S_l2,'r')
title('$\ell_2$','Interpreter','LaTeX')
legend('30N','15N','15S','30S')
set(gca,'FontSize',15)

%% difference between matrix and background

% this section computes the temperature change for each of the matrix runs,
% i.e. how T0, T1, and T2 have changed relative to the background.

inj30N_T0_diff = inj30N_T0 - backgd_T0(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));
inj30N_T1_diff = inj30N_T1 - backgd_T1(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));
inj30N_T2_diff = inj30N_T2 - backgd_T2(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));

inj15N_T0_diff = inj15N_T0 - backgd_T0(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));
inj15N_T1_diff = inj15N_T1 - backgd_T1(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));
inj15N_T2_diff = inj15N_T2 - backgd_T2(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));

inj15S_T0_diff = inj15S_T0 - backgd_T0(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));
inj15S_T1_diff = inj15S_T1 - backgd_T1(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));
inj15S_T2_diff = inj15S_T2 - backgd_T2(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));

inj30S_T0_diff = inj30S_T0 - backgd_T0(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));
inj30S_T1_diff = inj30S_T1 - backgd_T1(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));
inj30S_T2_diff = inj30S_T2 - backgd_T2(background_yrs >= matrix_yrs(1) & background_yrs <= matrix_yrs(end));

%% linear fit

% this section turns all of the temperature change and AOD data into long
% vectors and computes linear fits in order to extract sensitivities (i.e.,
% how much a change in l0 causes a change in T0, and so on)

dl0 = cat(1,inj30N_l0,inj15N_l0,inj15S_l0,inj30S_l0);
dl1 = cat(1,inj30N_l1,inj15N_l1,inj15S_l1,inj30S_l1);
dl2 = cat(1,inj30N_l2,inj15N_l2,inj15S_l2,inj30S_l2);

dT0 = cat(1,inj30N_T0_diff,inj15N_T0_diff,inj15S_T0_diff,inj30S_T0_diff);
dT1 = cat(1,inj30N_T1_diff,inj15N_T1_diff,inj15S_T1_diff,inj30S_T1_diff);
dT2 = cat(1,inj30N_T2_diff,inj15N_T2_diff,inj15S_T2_diff,inj30S_T2_diff);

% T0 ~ l0 + l1 + l2
model1 = fitlm([dl0,dl1,dl2],dT0,'Intercept',false,'VarNames',{'l0','l1','l2','T0'});
% T0 ~ l0
model2 = fitlm(dl0,dT0,'Intercept',false,'VarNames',{'l0','T0'});
% T1 ~ l0 + l1 + l2
model3 = fitlm([dl0,dl1,dl2],dT1,'Intercept',false,'VarNames',{'l0','l1','l2','T1'});
% T1 ~ l0 + l1
model4 = fitlm([dl0,dl1],dT1,'Intercept',false,'VarNames',{'l0','l1','T1'});
% T2 ~ l0 + l1 + l2
model5 = fitlm([dl0,dl1,dl2],dT2,'Intercept',false,'VarNames',{'l0','l1','l2','T2'});

% extract sensitivities
dT0dl0 = table2array(model2.Coefficients(1,1)); dT0dl0_err = table2array(model2.Coefficients(1,2));
dT1dl0 = table2array(model4.Coefficients(1,1)); dT1dl0_err = table2array(model4.Coefficients(1,2));
dT1dl1 = table2array(model4.Coefficients(2,1)); dT1dl1_err = table2array(model4.Coefficients(2,2));
dT2dl0 = table2array(model5.Coefficients(1,1)); dT2dl0_err = table2array(model5.Coefficients(1,2));
dT2dl1 = table2array(model5.Coefficients(2,1)); dT2dl1_err = table2array(model5.Coefficients(2,2));
dT2dl2 = table2array(model5.Coefficients(3,1)); dT2dl2_err = table2array(model5.Coefficients(3,2));

%% changes in target variables

% this section determines how T0, T1, and T2 will change over time in the
% background run; in other words, how much temperature change the SAI run
% will need to offset. This section also plots the expected changes as a
% sanity check.

% temperature targets (also called "reference values")
target_T0 = mean(backgd_T0(background_yrs >= reference_start & background_yrs <= reference_end));
target_T1 = mean(backgd_T1(background_yrs >= reference_start & background_yrs <= reference_end));
target_T2 = mean(backgd_T2(background_yrs >= reference_start & background_yrs <= reference_end));

% how the background simulation deviates from the temperature targets
backgd_T0_diff = backgd_T0(background_yrs >= ceil(mean([reference_start reference_end]))) - target_T0;
backgd_T1_diff = backgd_T1(background_yrs >= ceil(mean([reference_start reference_end]))) - target_T1;
backgd_T2_diff = backgd_T2(background_yrs >= ceil(mean([reference_start reference_end]))) - target_T2;

% line of best fit to the deviations from the temperature targets
T0_fit = polyfit(ceil(mean([reference_start reference_end])):background_yrs(end),backgd_T0_diff,1);
T1_fit = polyfit(ceil(mean([reference_start reference_end])):background_yrs(end),backgd_T1_diff,1);
T2_fit = polyfit(ceil(mean([reference_start reference_end])):background_yrs(end),backgd_T2_diff,1);

figure;
set(gcf,'units','inches','position',[0,0,15,5])

subplot('Position',pos1)
hold on
plot(ceil(mean([reference_start reference_end])):background_yrs(end),backgd_T0_diff)
plot(ceil(mean([reference_start reference_end])):background_yrs(end),polyval(T0_fit,ceil(mean([reference_start reference_end])):background_yrs(end)))
title('fit to change in T_0')
xlabel('year')
ylabel('T_0 change (K)')
set(gca,'FontSize',15)

subplot('Position',pos2)
hold on
plot(ceil(mean([reference_start reference_end])):background_yrs(end),backgd_T1_diff)
plot(ceil(mean([reference_start reference_end])):background_yrs(end),polyval(T1_fit,ceil(mean([reference_start reference_end])):background_yrs(end)))
title('fit to change in T_1')
xlabel('year')
ylabel('T_1 change (K)')
set(gca,'FontSize',15)

subplot('Position',pos3)
hold on
plot(ceil(mean([reference_start reference_end])):background_yrs(end),backgd_T2_diff)
plot(ceil(mean([reference_start reference_end])):background_yrs(end),polyval(T2_fit,ceil(mean([reference_start reference_end])):background_yrs(end)))
title('fit to change in T_2')
xlabel('year')
ylabel('T_2 change (K)')
set(gca,'FontSize',15)

%% feeds

% This section computes the feedforward gains needed to offset the expected
% changes in temperature computed in the previous section.

% compute the l0 feed
l0_feed = -T0_fit(1) / dT0dl0;

% compute the l1 feed
T1_change_due_to_l0 = l0_feed*dT1dl0;
T1_leftover = -T1_fit(1) - T1_change_due_to_l0;
l1_feed = T1_leftover / dT1dl1;

% if the l1 feed is larger than the l0 feed, it must be reduced
if abs(l1_feed) > l0_feed
    l1_feed_constrained = l0_feed*sign(l1_feed_constrained);
else
    l1_feed_constrained = l1_feed;
end

% compute the l2 feed
T2_change_due_to_l0 = l0_feed*dT2dl0;
T2_change_due_to_l1 = l1_feed_constrained*dT2dl1;
T2_leftover = -T2_fit(1) - T2_change_due_to_l0 - T2_change_due_to_l1;
l2_feed = T2_leftover / dT2dl2;

% l2 cannot be negative
if l2_feed < 0
    l2_feed_constrained = 0;
else
    l2_feed_constrained = l2_feed;
end

% |l1| + l2 must not > l0
if l2_feed_constrained > l0_feed - abs(l1_feed_constrained)
    l2_feed_constrained = l0_feed - abs(l1_feed_constrained);
end

%% injections

% finalized l0, l1, l2
l0 = l0_feed;
l1S = max([-l1_feed_constrained,0]);
l1N = max([l1_feed_constrained,0]);
l2 = l2_feed_constrained;

% injection rates assuming original "M" matrix
inj_30N = 20*l1N + 40*l2;
inj_15N = 30*(l0 - l1N - l1S - l2) + 45*l1N;
inj_15S = 30*(l0 - l1N - l1S - l2) + 45*l1S;
inj_30S = 20*l1S + 40*l2;

% injection rate if you just want the l0 feedforward (useful when l1 or l2
% are very small, or T1 and T2 very noisy and you don't trust the
% "expected" changes)
inj_for_just_l0 = 30*l0;

%% new M matrix?

% this function computes the new matrix "M" for this model, which relates
% injections at each of the individual locations to changes in l0, l1, and
% l2.

% how l0 changes based on injections at specific latitudes
dl0d30N = mean(inj30N_l0(6:end))/injection_rate;
dl0d15N = mean(inj15N_l0(6:end))/injection_rate;
dl0d15S = mean(inj15S_l0(6:end))/injection_rate;
dl0d30S = mean(inj30S_l0(6:end))/injection_rate;

% how l1 changes based on injections at specific latitudes
dl1d30N = mean(inj30N_l1(6:end))/injection_rate;
dl1d15N = mean(inj15N_l1(6:end))/injection_rate;
dl1d15S = mean(inj15S_l1(6:end))/injection_rate;
dl1d30S = mean(inj30S_l1(6:end))/injection_rate;

% how l2 changes based on injections at specific latitudes
dl2d30N = mean(inj30N_l2(6:end))/injection_rate;
dl2d15N = mean(inj15N_l2(6:end))/injection_rate;
dl2d15S = mean(inj15S_l2(6:end))/injection_rate;
dl2d30S = mean(inj30S_l2(6:end))/injection_rate;

% use fmincon to find the optimal combinations:
options = optimoptions(@fmincon,'Display','off');
% l0 = 1, constrained by l1 = 0; inj15N, inj15S >= 0
data_l0 = fmincon(@(x)abs([dl0d15S,dl0d15N]*x - 1), [1/(dl0d15S + dl0d15N);1/(dl0d15S + dl0d15N)],[],[],[dl1d15S,dl1d15N],0,[0;0],[Inf;Inf],[],options);
% l1N = 1, constrained by l0 = 1; inj30N, inj15N >= 0;
data_l1N = fmincon(@(x)abs([dl1d15N,dl1d30N]*x - 1), [dl0d15N, dl0d30N; dl1d15N, dl1d30N]\[1;1],[],[],[dl0d15N,dl0d30N],1,[0;0],[Inf;Inf],[],options);
% l1S = 1, constrained by l0 = 1; inj30S, inj15S >= 0;
data_l1S = fmincon(@(x)abs([dl1d30S,dl1d15S]*x + 1), [dl0d30S, dl0d15S; dl1d30S, dl1d15S]\[1;-1],[],[],[dl0d30S,dl0d15S],1,[0;0],[Inf;Inf],[],options);
% l2 = 1, constrained by l0 = 1, l1 = 0; inj30S, inj30N >= 0;
data_l2 = fmincon(@(x)abs([dl2d30S,dl2d30N]*x - 1), [dl0d30S, dl0d30N; dl2d30S, dl2d30N]\[1;1],[],[],[dl0d30S,dl0d30N;dl1d30S,dl1d30N],[1;0],[0;0],[Inf;Inf],[],options);

% new M matrix
M = zeros(4);
M(2:3,1) = data_l0;
M(3:4,2) = data_l1N;
M(1:2,3) = data_l1S;
M([1,4],4) = data_l2;

% new injection rates for feedforwards
F = [1,1,1,1;0,1,0,0;0,0,1,0;0,0,0,1];
new_inj_all = M*inv(F)*[l0;l1S;l1N;l2];
new_inj_for_just_l0 = M*inv(F)*[l0;0;0;0];

%% Print out results

% temperature targets (2020-2039 average)
fprintf('\n')
fprintf('Temperature Targets\n')
fprintf('T0 target: %4.4f\n',target_T0)
fprintf('T1 target: %4.4f (%4.4f)\n',target_T1,target_T1/3)
fprintf('T2 target: %4.4f (%4.4f)\n',target_T2,target_T2/5)

% sensitivities
fprintf('\n')
fprintf('Sensitivities\n')
fprintf('T0 sensitivity to l0: %4.4f +/- %4.4f\n',dT0dl0,dT0dl0_err)
fprintf('T1 sensitivity to l0: %4.4f +/- %4.4f\n',dT1dl0,dT1dl0_err)
fprintf('T1 sensitivity to l1: %4.4f +/- %4.4f\n',dT1dl1,dT1dl1_err)
fprintf('T2 sensitivity to l0: %4.4f +/- %4.4f\n',dT2dl0,dT2dl0_err)
fprintf('T2 sensitivity to l1: %4.4f +/- %4.4f\n',dT2dl1,dT2dl1_err)
fprintf('T2 sensitivity to l2: %4.4f +/- %4.4f\n',dT2dl2,dT2dl2_err)

% behavior of target metrics under SSP
fprintf('\n')
fprintf('Behavior of target metrics under SSP that needs to be offset\n')
fprintf('T0 behavior: %+4.4f per year\n',T0_fit(1));
fprintf('T1 behavior: %+4.4f per year\n',T1_fit(1));
fprintf('T2 behavior: %+4.4f per year\n',T2_fit(1));

% amount of l0/l1/l2 needed (feedforward)
fprintf('\n')
fprintf('amount of l0/l1/l2 needed (feedforward)\n')
fprintf('\n')
fprintf('l0:\n')
fprintf('%4.4f K/yr / %4.4f K/l0 = %4.4f l0/yr\n',-T0_fit(1),dT0dl0,l0_feed)
fprintf('\n')
fprintf('l1:\n')
fprintf('%+4.4f K/yr required in total\n',-T1_fit(1))
fprintf('%4.4f K/l0 * %4.4f l0/yr = %4.4f K/yr from l0\n',dT1dl0,l0_feed,T1_change_due_to_l0)
fprintf('%+4.4f K/yr T1 leftover\n',T1_leftover)
fprintf('%4.4f K/yr / %4.4f K/l1 = %4.4f l1/yr\n',T1_leftover,dT1dl1,l1_feed)
if l1_feed_constrained ~= l1_feed
    fprintf('controller constraint: new l1 = %4.4f\n',l1_feed_constrained);
end
fprintf('\n')
fprintf('l2:\n')
fprintf('%+4.4f K/yr required in total\n',-T2_fit(1))
fprintf('%4.4f K/l0 * %4.4f l0/yr = %4.4f K/yr from l0\n',dT2dl0,l0_feed,T2_change_due_to_l0)
fprintf('%4.4f K/l1 * %4.4f l1/yr = %4.4f K/yr from l1\n',dT2dl1,l1_feed_constrained,T2_change_due_to_l1)
fprintf('%+4.4f K/yr T2 leftover\n',T2_leftover)
fprintf('%4.4f K/yr / %4.4f K/l2 = %4.4f l2/yr\n',T2_leftover,dT2dl2,l2_feed)
if l2_feed_constrained ~= l2_feed
    fprintf('controller constraint: new l2 = %4.4f\n',l2_feed_constrained);
end

% final values of l0, l1, l2
fprintf('\n')
fprintf('Final value of l0 feedforward: %4.4f\n',l0_feed)
fprintf('Final value of l1 feedforward: %4.4f\n',l1_feed_constrained)
fprintf('Final value of l2 feedforward: %4.4f\n',l2_feed_constrained)

% injection rates (feedforward)
fprintf('\n')
fprintf('injection rates for just l0 feed:\n')
fprintf('15N and 15S: %4.4f\n',inj_for_just_l0)
fprintf('\n')
fprintf('injection rates for all feeds:\n')
fprintf('30N: %4.4f\n',inj_30N)
fprintf('15N: %4.4f\n',inj_15N)
fprintf('15S: %4.4f\n',inj_15S)
fprintf('30S: %4.4f\n',inj_30S)

% new M matrix?
fprintf('\n')
fprintf('new injection matrix "M":\n')
M
fprintf('\n')
fprintf('injection rates for just l0 feed with new "M":\n')
fprintf('30N: %4.4f\n',new_inj_for_just_l0(1))
fprintf('15N: %4.4f\n',new_inj_for_just_l0(2))
fprintf('15S: %4.4f\n',new_inj_for_just_l0(3))
fprintf('30S: %4.4f\n',new_inj_for_just_l0(4))
fprintf('\n')
fprintf('injection rates for all feeds with new "M":\n')
fprintf('30N: %4.4f\n',new_inj_all(1))
fprintf('15N: %4.4f\n',new_inj_all(2))
fprintf('15S: %4.4f\n',new_inj_all(3))
fprintf('30S: %4.4f\n',new_inj_all(4))

% feedback gains
fprintf('\n')
fprintf('Feedback gains:\n')
fprintf('\n')
fprintf('l0:\n')
fprintf('GLENS gain: 0.028\n')
fprintf('GLENS sensitivity: -5.2 T0/l0\n')
fprintf('New sensitivity: %4.4f T0/l0\n',dT0dl0)
fprintf('New l0 gain: %4.4f\n',0.028*-5.2/dT0dl0)

fprintf('\n')
fprintf('l1:\n')
fprintf('GLENS gain: 0.13\n')
fprintf('GLENS sensitivity: -4.4 T1/l1\n')
fprintf('New sensitivity: %4.4f T1/l1\n',dT1dl1)
fprintf('New l1 gain: %4.4f\n',0.13*-4.4/dT1dl1)

fprintf('\n')
fprintf('l2:\n')
fprintf('GLENS gain: 0.39\n')
fprintf('GLENS sensitivity: -1.6 T2/l2\n')
fprintf('New sensitivity: %4.4f T2/l2\n',dT2dl2)
fprintf('New l1 gain: %4.4f\n',0.39*-1.6/dT2dl2)

%% Annual_average.m

% The ultimate annual mean calculator!

% Walker Raymond Lee
% Version 2.1 (March 2022)
% Cornell University (wl644@cornell.edu)
% MacMartin Research Group

%% Version history

% 2.1 (Mar '22): now allows for [] inputs, for example if you want to specify a 360-day calendar but don't care about times/years
% 2.0 (Jul '21): can now compute seasonal averages; also should now work on data with multiple entries per month
% 1.2 (Jul '21): patched a bug wherein years with only one sample of data would throw an error
% 1.1 (Jun '21): patched a bug wherein 360-day calendar means computed incorrectly
% 1.0 (Mar '21): original release

%% Documentation

% syntax: [output,years,months,month_names,times_formatted] = annual_average(input,times,start_date,calendar,season)

% ----------------------------------------
% Input arguments
% ----------------------------------------

% input: the matrix you want to take the annual mean of. Can be up to 5
% dimensions.
%
% times (optional): the time stamps for your data. The list of times should
% come directly from an output file and have the format "days since XXX",
% where XXX is the start date. If you include this, the script will figure
% out which dimension corresponds to time and use your time stamps to
% compute annual means; otherwise, it will assume the last dimension is
% time and average over every 12 entries, weighting by the number of days
% in each month.
%
% start date (optional): the reference date for your "times" input. Must be
% in the format [day,month,year]; for example, if the date format is
% "days since 1-1-2010," enter start_date as [1,1,2010]. If you skip this
% input, the "years" output will equal "years since start of simulation"
% (i.e., the script will use a dummy start date of [1,1,1]).
%
% calendar (optional): whether your model's calendar is '365 day' or '360
% day'. If skipped, the default setting is '365 day'.
%
% season (optional): this must be one of the following:
%      'ANN' - average over the whole year (default if left empty)
%      'MAM' - average over the spring (March, April, and May)
%      'JJA' - average over the summer (June, July, and August)
%      'SON' - average over the autumn (September, October, and November)
%      'DJF' - average over the winter (December, January, and February)
% Because consecutive December, January, and February are in different
% calendar years, when computing the DJF average of (for example) 2000,
% this script will look at Dec. of 2000 and Jan/Feb of 2001. As such, the
% average for last year of whatever data you're looking at will read NaN
% (for example if you want to find DJF averages for 2000-2010, the DJF
% average for the year 2010 will return NaN because the code will look for
% Jan and Feb of 2011, which isn't there).

% ----------------------------------------
% Output arguments
% ----------------------------------------

% "output" - the annual means of the input data. If you included times and a start date, time will be the first dimension of the output matrix.
% "years" - the years corresponding to the annual means
% "months" - which month each input time corresponds to (numbers 1-12)
% "month_names" - which month each input time corresponds to (strings)
% "times_formatted" - the input times, reformatted as "year.fraction_of_year" (good for troubleshooting)
%
% (note that in "times_formatted", Dec. 31 of the year 1999 reads as
% 2000.00, but it was counted as part of 1999 during averaging)


%%

function [output,years,months,month_names,times_formatted] = annual_average(input,varargin)

%% step 0: if there is no times/start date/calendar/season input, make up dummies

input_dims = size(input); % dimensions of the input matrix
n_dims = length(input_dims); % number of input dimensions

% check if input matrix is too big
if n_dims > 5
    error('This code only supports input fields with 5 dimensions or fewer. Please input a smaller matrix, or send Walker an angry email')
end

% check number of arguments
if nargin < 0 || nargin > 5
    error('This script needs 1-5 input arguments. Please check your input, see the documentation, or send Walker an angry email');
end

% argument 5: season
try
    season = varargin{4}; % if there are 5 inputs, input 5 = calendar
catch
    season = [];
end

if isempty(season)
    season = 'ANN';
end

if strcmp(season,'ANN') == 0 && strcmp(season,'MAM') == 0 && strcmp(season,'JJA') == 0 && strcmp(season,'SON') == 0 && strcmp(season,'DJF') == 0
    error('The season variable must be "ANN", "MAM", "JJA", "SON", or "DJF". Please check your input, see the documentation, or send Walker an angry email');
end

% argument 4: calendar
try
    calendar = varargin{3};
catch
    calendar = [];
end

if isempty(calendar)
    calendar = '365 day';
end

if strcmp(calendar,'360 day') == 0 && strcmp(calendar,'365 day') == 0
    error('The calendar must be "365 day" or "360 day". Please check your "calendar" input or send Walker an angry email');
end

if strcmp(calendar,'365 day') == 1
    daysinmonth = cumsum([31,28,31,30,31,30,31,31,30,31,30,31]);
elseif strcmp(calendar,'360 day') == 1
    daysinmonth = cumsum([30,30,30,30,30,30,30,30,30,30,30,30]);
end

% argument 3: start date
try
    start_date = varargin{2};
catch
    start_date = [];
end

if isempty(start_date)
    start_date = [1,1,1];
end

if length(start_date) ~= 3
    error('The syntax for "start date" needs to be [day, month, year]. Please check your syntax, or send Walker an angry email');
end

% argument 2: times
try
    times = varargin{1};
catch
    times = [];
end

if isempty(times)
    times = zeros(input_dims(end),1);
    for i = 1:length(times)
        if strcmp(calendar,'365 day') == 1
            times(i) = 365*floor((i-0.1)/12) + daysinmonth(rem(i-1,12)+1);
        elseif strcmp(calendar,'360 day') == 1
            times(i) = 360*floor((i-0.1)/12) + daysinmonth(rem(i-1,12)+1);
        end
    end
end

%% step 0.5: reorganize array such that time is now the first dimension

% number of time steps
n_times = length(times);
% find dimension number corresponding to times
time_dim = find(input_dims == n_times);

if isnan(time_dim)
    error("The length of your time vector doesn't match any of the data dimensions. Please double check, or send Walker an angry email")
end

% move the time dimension to position 1
new_order = 1:n_dims;
new_order(new_order==time_dim) = [];
new_order = [time_dim, new_order];
% change the demension order of the input
input = permute(input,new_order);

%% step 1: convert time data from "days since ___" to "year.days"

% based on calendar input, determine how many days have passed by the start
% of each month
if strcmp(calendar,'365 day') == 1
    days_passed = cumsum([0,31,28,31,30,31,30,31,31,30,31,30,31]);
    days_in_year = 365;
elseif strcmp(calendar,'360 day') == 1
    days_passed = cumsum([0,30,30,30,30,30,30,30,30,30,30,30,30]);
    days_in_year = 360;
else
    error('The calendar must be "365 day" or "360 day". Please check your "calendar" input');
end

months = zeros(length(times),1);

% determine which month each time step belongs to
for i = 1:length(times)
    if times(i) - floor(times(i)/days_in_year)*days_in_year == 0
        months(i) = 12;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(2)
        months(i) = 1;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(3)
        months(i) = 2;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(4)
        months(i) = 3;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(5)
        months(i) = 4;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(6)
        months(i) = 5;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(7)
        months(i) = 6;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(8)
        months(i) = 7;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(9)
        months(i) = 8;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(10)
        months(i) = 9;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(11)
        months(i) = 10;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(12)
        months(i) = 11;
    elseif times(i) - floor(times(i)/days_in_year)*days_in_year <= days_passed(13)
        months(i) = 12;
    end
end

% convert times data to "year.frac_of_year"
times_formatted = start_date(3) + (days_passed(start_date(2))+start_date(1)+times-1)/days_in_year;

%% step 2: compute annual means

% round time down to the year; this means midnight on Dec. 31 will count as
% the correct year instead of carrying over to the next year
yrs_temp = floor(times_formatted-0.01);
data_temp = input;
months_temp = months;

% outputs start out as empty arrays
years = []; output = [];


while isempty(yrs_temp) == 0 % if the rounded list of times is not empty

    % find all the entries for the same year as the 1st entry, and the
    % year immediately after (needed for DJF)
    i = find(yrs_temp == yrs_temp(1));
    j = find(yrs_temp == yrs_temp(1)+1);
    
    % add this year to the list of years
    years = cat(1,years,yrs_temp(1));
    
    % weight by days per month, if using 365-day calendar. Also sort
    % into monthly bins (in case of multiple entries per month)
    this_years_data = data_temp(i,:,:,:,:);
    this_years_months = months_temp(i);
    next_years_data = data_temp(j,:,:,:,:);
    next_years_months = months_temp(j);

    if strcmp(calendar,'365 day') == 1
        this_years_data_sorted(1,:,:,:,:,:) = mean(this_years_data(this_years_months==1,:,:,:,:,:),1)*31/365;
        this_years_data_sorted(2,:,:,:,:,:) = mean(this_years_data(this_years_months==2,:,:,:,:,:),1)*28/365;
        this_years_data_sorted(3,:,:,:,:,:) = mean(this_years_data(this_years_months==3,:,:,:,:,:),1)*31/365;
        this_years_data_sorted(4,:,:,:,:,:) = mean(this_years_data(this_years_months==4,:,:,:,:,:),1)*30/365;
        this_years_data_sorted(5,:,:,:,:,:) = mean(this_years_data(this_years_months==5,:,:,:,:,:),1)*31/365;
        this_years_data_sorted(6,:,:,:,:,:) = mean(this_years_data(this_years_months==6,:,:,:,:,:),1)*30/365;
        this_years_data_sorted(7,:,:,:,:,:) = mean(this_years_data(this_years_months==7,:,:,:,:,:),1)*31/365;
        this_years_data_sorted(8,:,:,:,:,:) = mean(this_years_data(this_years_months==8,:,:,:,:,:),1)*31/365;
        this_years_data_sorted(9,:,:,:,:,:) = mean(this_years_data(this_years_months==9,:,:,:,:,:),1)*30/365;
        this_years_data_sorted(10,:,:,:,:,:) = mean(this_years_data(this_years_months==10,:,:,:,:,:),1)*31/365;
        this_years_data_sorted(11,:,:,:,:,:) = mean(this_years_data(this_years_months==11,:,:,:,:,:),1)*30/365;
        this_years_data_sorted(12,:,:,:,:,:) = mean(this_years_data(this_years_months==12,:,:,:,:,:),1)*31/365;
        next_years_data_sorted(1,:,:,:,:,:) = mean(next_years_data(next_years_months==1,:,:,:,:,:),1)*31/365;
        next_years_data_sorted(2,:,:,:,:,:) = mean(next_years_data(next_years_months==2,:,:,:,:,:),1)*28/365;
    end

    % calculate means for this year
    
    if strcmp(season,'ANN') == 1 % annual average:
        
        if strcmp(calendar,'365 day') == 1
            average = sum(this_years_data_sorted,1);
        elseif strcmp(calendar,'360 day') == 1
            average = mean(this_years_data_sorted,1);
        end
        output = cat(1,output,average);
        
    elseif strcmp(season,'MAM') == 1 % spring average:
        
        if strcmp(calendar,'365 day') == 1
            average = sum(this_years_data_sorted(3:5,:,:,:,:,:),1) * 365/(31+30+31);
        elseif strcmp(calendar,'360 day') == 1
            average = mean(this_years_data_sorted(3:5,:,:,:,:,:),1);
        end
        output = cat(1,output,average);
        
    elseif strcmp(season,'JJA') == 1 % summer average:
        
        if strcmp(calendar,'365 day') == 1
            average = sum(this_years_data_sorted(6:8,:,:,:,:,:),1) * 365/(30+31+31);
        elseif strcmp(calendar,'360 day') == 1
            average = mean(this_years_data_sorted(6:8,:,:,:,:,:),1);
        end
        output = cat(1,output,average);
        
    elseif strcmp(season,'SON') == 1 % autumn average:
        
        if strcmp(calendar,'365 day') == 1
            average = sum(this_years_data_sorted(9:11,:,:,:,:,:),1) * 365/(30+31+30);
        elseif strcmp(calendar,'360 day') == 1
            average = mean(this_years_data_sorted(9:11,:,:,:,:,:),1);
        end
        output = cat(1,output,average);
        
    elseif strcmp(season,'DJF') == 1 % winter average:
        
        if strcmp(calendar,'365 day') == 1
            average = sum([this_years_data_sorted(12,:,:,:,:,:);next_years_data_sorted(1:2,:,:,:,:,:)],1) * 365/(31+31+28);
        elseif strcmp(calendar,'360 day') == 1
            average = mean([this_years_data_sorted(12,:,:,:,:,:);next_years_data_sorted(1:2,:,:,:,:,:)]);
        end
        output = cat(1,output,average);
        
    end
    
    % delete year, month, and data from temporary lists
    months_temp(i) = [];
    yrs_temp(i) = [];
    data_temp(i,:,:,:,:) = [];
    
end

%% bookkeeping: return month names

month_names = cell(length(months),1);

for i = 1:length(months)
    if months(i) == 1
        month_names{i} = 'January';
    elseif months(i) == 2
        month_names{i} = 'February';
    elseif months(i) == 3
        month_names{i} = 'March';
    elseif months(i) == 4
        month_names{i} = 'April';
    elseif months(i) == 5
        month_names{i} = 'May';
    elseif months(i) == 6
        month_names{i} = 'June';
    elseif months(i) == 7
        month_names{i} = 'July';
    elseif months(i) == 8
        month_names{i} = 'August';
    elseif months(i) == 9
        month_names{i} = 'September';
    elseif months(i) == 10
        month_names{i} = 'October';
    elseif months(i) == 11
        month_names{i} = 'November';
    elseif months(i) == 12
        month_names{i} = 'December';
    end
end

end % end of script