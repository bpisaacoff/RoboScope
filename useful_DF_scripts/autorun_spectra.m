folder='T:\Lab Members\Ben Isaacoff\Data\6_06_18\';

% DFii=1;

fname=['DF_',num2str(DFii)];

stepnglue([folder,fname],400,1000,40,'coarse');

DF_im=[fname,'.mat'];
DF_imaging_mode


DFii=DFii+1;