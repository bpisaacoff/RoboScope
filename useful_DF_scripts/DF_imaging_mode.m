%% Darkfield_imaging_mode
% written by BPI 5/1/15
%
% This script calculates and plots the corrected darkfield spectrum for
% data input as an image. Either .asc or .tif works.
%
% note that you need the function isAprxEq to be in your Matlab folder. You
% can find it on the server in Matlab\code
%
% NOTE to achieve the highest possible accuracy you should click to choose
% the background as near to the NP as possible, the auto find usually
% chooses an area far away, but this is very slightly less accurate. For
% everyday accuracy needs auto find is fine, but for high accuracy
% situations you should click.
%
% last update 11/11/15 BPI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control Panel
% This is the only section you should need to edit in normal circumstances

% the folder where the data is stored, with the slash at the end Z:\ or X:\
% folder='T:\Lab Members\Ben Isaacoff\Data\6_04_18\';

% the filenames for the images, including the extension (.asc, or .tif)
% DF_im='DF_18.mat';%The darkfield file
BF_im='BF.mat';% The brightfield image

%Boolean whether to use the optional background file specified below
sep_BG=0;% Should be uniform
DF_BG_im='nta1000bg.asc';% The optional darkfiled background file

% Boolean, set to 1 to have the program autofind the GNR and the BG
% and set to 0 to click to choose the areas yourself
auto_find=0;

%Fit to a Lorenztian? Enter zero to not fit and your peak guess if you do
%want to fit
Lor_pg=650;%in nm please

% plotting booleans
plot_raw_DF=0;% plot the image showing the chosen regions and the raw DF spectra
plot_BF=0;% plot the brightfield spectrum

diff_wid=20;%Diffraction width in pixels, EVEN NUMBER PLEASE, this is the width of the region around the GNR

%Change the boolean for whichever ND filter used, only the ones with red
% tape on top. If you use the orange ones (or some other) you need to
% measure those yourself.
which_ND_filter1=0;%ND1
which_ND_filter2=1;%ND2
which_ND_filter3=0;%ND03

srv_lett='T';%The letter that you've mapped the server to

X_lims=[450,800];%The limits of the x-axis that you want to plot

%% Importing & Finding

%import the BF
if strcmp(BF_im((end-3):end),'.asc')
    imptemp=importdata([folder,BF_im]);
    if size(imptemp,2)~=513% check that it's oriented so that the pixels are in the second dimension
        imptemp=imptemp';
    end
    wl=imptemp(:,1);
    imptemp=imptemp(:,2:end);%getting rid of the wl vector
elseif strcmp(BF_im((end-3):end),'.mat')
    load([folder,BF_im],'spectrum');
    imptemp=spectrum;
    if size(imptemp,2)~=513% check that it's oriented so that the pixels are in the second dimension
        imptemp=imptemp';
    end
    wlimp=imptemp(:,1);
    wl=linspace(400,1000,601)';
    imptemp=imptemp(:,2:end);%getting rid of the wl vector
    for ii=1:size(imptemp,2)
    imptemp1(:,ii)=interp1(wlimp,imptemp(:,ii),wl);
    end
imptemp=imptemp1;
elseif strcmp(BF_im((end-3):end),'.tif')
    imptemp=importdata([folder,BF_im]);
    if size(imptemp,2)~=512% check that it's oriented so that the pixels are in the second dimension
        imptemp=imptemp';
    end
end

%Find the vertical edges of the slit.
% %NOTE this code is weird (aka not reccomended), don't worry about how it works
mean_im=mean(imptemp,1);%The mean of all wavelength for each vertical pixel
pdiff=abs(smooth(diff(smooth(mean_im,10)),10));%basically the derivative of the previous line
[~,diff_fpsloc,~,proms]=findpeaks(pdiff/max(pdiff));%find the local maxima
filtered_fps=diff_fpsloc(proms>.02);%filter the peaks by prominence
edge1=100;edge2=380;%pull out the first and last peaks that pass the filter
if ~(isAprxEq(edge1,103,.2) && isAprxEq(edge2,375,.2)) && ~(isAprxEq(edge1,137,.2) && isAprxEq(edge2,409,.2))
    warning('Edge detection likely failed, check the borders')
end
edge1=edge1+10;edge2=edge2-10;%shifting the edges to account for the poor edge finding

%keeping only the data in the slit for the BF
BF=mean(imptemp(:,:),2);%mean of all pixel rows inside the slit
dark_pxs=0.5*(mean(imptemp(:,1:(edge1-50)),2)+mean(imptemp(:,(edge2+50):512),2));%average of the dark pixels
BF=BF-dark_pxs;%subtracting the dark pixels from the BF signal

%import the DF
if strcmp(DF_im((end-3):end),'.asc')
    imptemp=importdata([folder,DF_im]);
    if size(imptemp,2)~=513% check that it's oriented so that the pixels are in the second dimension
        imptemp=imptemp';
    end
    if imptemp(:,1)~=wl; error('Wavelength vectors do not match!');end
    imptemp=imptemp(:,2:end);%getting rid of the wl vector
elseif strcmp(DF_im((end-3):end),'.mat')
    load([folder,DF_im],'spectrum','img');
    imptemp=spectrum;
    if size(imptemp,2)~=513% check that it's oriented so that the pixels are in the second dimension
        imptemp=imptemp';
    end
    wlimp=imptemp(:,1);
    wl=linspace(400,1000,601)';
    imptemp=imptemp(:,2:end);%getting rid of the wl vector
    for ii=1:size(imptemp,2)
    imptemp1(:,ii)=interp1(wlimp,imptemp(:,ii),wl);
    end
imptemp=imptemp1;
elseif strcmp(DF_im((end-3):end),'.tif')
    imptemp=importdata([folder,DF_im]);
    if size(imptemp,2)~=512% check that it's oriented so that the pixels are in the second dimension
        imptemp=imptemp';
    end
end

if auto_find
    % pulling out the actual data
    act_data=imptemp(:,edge1:edge2);
    
    %find the diffraction width section of the wavelength collapsed  data
    %with the max and the min intensities, corresponding to the GNR
    %location and the BG location
    wl_clps_data=mean(act_data,1);%collapsing the data to the pixel dimension
    min_int=inf;max_int=0;
    for ii=1:(length(wl_clps_data)-diff_wid)
        cur_int=mean(wl_clps_data((ii):(diff_wid+ii)));
        if cur_int<min_int
            min_int=cur_int;
            BG_loc=ii+diff_wid/2;
        end
        if cur_int>max_int
            max_int=cur_int;
            GNR_loc=ii+diff_wid/2;
        end
    end
    %shifting the coordinates to the full image coordinates
   GNR_loc=GNR_loc+edge1;BG_loc=BG_loc+edge1;
else
    % This section is for the user to click to input the location of the
    % GNR & the BG
    %This bit modified from Beth's ROI_picker
    
    figure;
    pcolor(wl,1:size(imptemp,2),imptemp')
    shading('flat')
    axis tight
    xlim(X_lims);
    xlabel('wavelength (nm)')
    ylabel('y position (pixel)')
    title('Click first on the particle THEN click on the background')
    
    % Click and choose fiduciaries
    ROIcents = round(ginput(2));
    
    GNR_loc=ROIcents(1,2);
    BG_loc=ROIcents(2,2);
    close gcf %close the current figure
end

%pull out the data and sum over the vertical pixels
DF_raw=mean(imptemp(:,(GNR_loc-diff_wid/2):(GNR_loc+diff_wid/2)),2);
DF_BG=mean(imptemp(:,(BG_loc-diff_wid/2):(BG_loc+diff_wid/2)),2);

if sep_BG
    %import the DF_BG
    
    if strcmp(DF_BG_im((end-3):end),'.asc')
        imptemp1=importdata([folder,DF_BG_im]);
        if size(imptemp1,2)~=513% check that it's oriented so that the pixels are in the second dimension
            imptemp1=imptemp1';
        end
        if imptemp1(:,1)~=wl; error('Wavelength vectors do not match!');end
        imptemp1=imptemp1(:,2:end);%getting rid of the wl vector
       DF_BG=mean(imptemp1(:,edge1:edge2),2);%mean of all pixel rows inside the slit
    elseif strcmp(DF_BG_im((end-3):end),'.mat')
        load([folder,DF_BG_im],'subspect');
        imptemp1=subspect;
        if size(imptemp1,2)~=513% check that it's oriented so that the pixels are in the second dimension
            imptemp1=imptemp1';
        end
        if imptemp1(:,1)~=wl; error('Wavelength vectors do not match!');end
        imptemp1=imptemp1(:,2:end);%getting rid of the wl vector
       DF_BG=mean(imptemp1(:,edge1:edge2),2);%mean of all pixel rows inside the slit
    elseif strcmp(DF_BG_im((end-3):end),'.tif')
        imptemp1=importdata([folder,DF_BG_im]);
        if size(imptemp1,2)~=512% check that it's oriented so that the pixels are in the second dimension
            imptemp1=imptemp1';
        end
       DF_BG=mean(imptemp1(:,edge1:edge2),2);%mean of all pixel rows inside the slit
    end
end

%% Plotting
%show it
if plot_raw_DF
    %showing the boundaries
    raw_max=max(max(imptemp));
    imptemp(:,(GNR_loc-diff_wid/2))=raw_max;imptemp(:,(GNR_loc+diff_wid/2))=raw_max;
    if ~sep_BG
        imptemp(:,(BG_loc-diff_wid/2))=raw_max;imptemp(:,(BG_loc+diff_wid/2))=raw_max;
    end
   imptemp(:,edge1)=raw_max;imptemp(:,edge2)=raw_max;
    figure
    axis image
    hold on
    mesh(wl,1:size(imptemp,2),imptemp')
    xlabel('Wavelength (nm)');
    ylabel('vertical pixel number');
    title('Raw DF data image');
    
    figure;
    plot(wl,DF_raw/max(DF_raw),wl,DF_BG/max(DF_BG),'LineWidth',2)
    plot(wl,DF_raw,wl,DF_BG,'LineWidth',2)
    xlim(X_lims);
    xlabel('Wavelength (nm)');
    ylabel('Scattering')
    title('Raw DF collapsed spectra');
    legend('GNR','BG')
end

if plot_BF
    figure;
    plot(wl,BF,'LineWidth',2)
    xlim(X_lims);
    xlabel('Wavelength (nm)');
    ylabel('Scattering')
    title('Raw BF FVB');
end

%% Import the ND filters
ND_1=ones(601,2);ND_2=ones(601,2);ND_3=ones(601,2);
if which_ND_filter1==1
    imptemp=importdata([srv_lett,':\UVVis Data\Esther\R-ND1.txt']);
    ND_1=imptemp.data(:,[1,2]);
    if ND_1(:,1)~=wl; error('Wavelength vectors do not match');end
end
if which_ND_filter2==1
    imptemp=importdata([srv_lett,':\UVVis Data\Esther\R-ND2.txt']);
    ND_2=imptemp.data(:,[1,2]);
    if ND_2(:,1)~=wl; error('Wavelength vectors do not match');end
end
if which_ND_filter3==1
    imptemp=importdata([srv_lett,':\UVVis Data\Esther\R-ND03.txt']);
    ND_3=imptemp.data(:,[1,2]);
    if ND_3(:,1)~=wl; error('Wavelength vectors do not match');end
end
ND_total=ND_1(:,2).*ND_2(:,2).*ND_3(:,2)/100;

%% DF calculation

DF=(DF_raw-DF_BG)./(BF./ND_total);
%DF=(DF-min(DF(wl>X_lims(1,1)&wl<X_lims(1,2))))/(max(DF(wl>X_lims(1,1)&wl<X_lims(1,2)))-min(DF(wl>X_lims(1,1)&wl<X_lims(1,2))));
figure;
hold on
plot(wl,DF,'LineWidth',2)
xlim(X_lims);
xlabel('Wavelength (nm)');
ylabel('Scattering (a.u.)')
ylim([0,1.05*max(DF(100:500))])
box on
grid on


if Lor_pg~=0
    [~,actLor_pg_pos]=min(abs(wl-Lor_pg));
    actLor_pg=wl(actLor_pg_pos);
    
    lorFun=@(c,xdata) c(1)./(1+((xdata-c(2))/c(3)).^2)+c(4);
    
    fitwind=100;%nm around pg to try to fit
    wl2fit=wl(wl>=(Lor_pg-fitwind) & wl<=(Lor_pg+fitwind));
    DF2fit=DF(wl>=(Lor_pg-fitwind) & wl<=(Lor_pg+fitwind));
    
    fitguess=[DF(actLor_pg_pos),actLor_pg,25,0];
    lowerbound=[0,min(wl2fit),1,-1];
    upperbound=[1e2*fitguess(1),max(wl2fit),100,100];
    [fitted] = lsqcurvefit(lorFun,fitguess,wl2fit,DF2fit,lowerbound,upperbound,optimset('Display','off'));
    
    plot(wl2fit,lorFun(fitted,wl2fit),'LineWidth',2)
    title(['Corrected DF spectra \newline Lorentzian max at ',num2str(fitted(2),4),' nm, FWHM = ',num2str(2*fitted(3),3),' nm'])
    
    ylim([0,1.05*lorFun(fitted,fitted(2))])
    
end
hold off