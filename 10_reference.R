# RRHO reference toptable: one tissue x DKT parcel x MRI feature from 08.
# 10_RRHO_map.R writes processed_data/rrho/<ref_fullname>.csv and 10d_ggseg.R reads it,
# so both take the reference from here. Change it in this file only.
# ref_tissue <- 'posteriorcingulatecortex'
# ref_parcel <- 'entorhinal'
# ref_feature <- 'CurvInd'
# ref_tissue <- 'temporalcortex'
# ref_parcel <- 'parsorbitalis'
# ref_feature <- 'SurfArea_norm'
# ref_tissue <- 'temporalcortex'
# ref_parcel <- 'superiortemporal'
# ref_feature <- 'ThickAvg'
# ref_tissue <- 'Headofcaudatenucleus'
# ref_parcel <- 'insula'
# ref_feature <- 'ThickAvg'
ref_tissue <- 'frontalcortex'
ref_parcel <- 'lateraloccipital'
ref_feature <- 'ThickAvg'
ref_fullname <- paste(ref_tissue, ref_parcel, ref_feature, sep = '_')
