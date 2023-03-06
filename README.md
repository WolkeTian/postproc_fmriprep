# postproc_fmriprep
# Calculate functional connectome based outputs from fmriprep (ICA-aroma enabled; fsaverage output space enabled).
# Example:
# bidsout = 'E:\Projects\SD_PANAS\prep_out';
# default method to extract time series from roi/parcel: mean
# cal_fc_from_fmriprep(bidsout); 
# The other method to extract time series from roi/parcel: 1st PC
# cal_fc_from_fmriprep(bidsout, 'pca');


