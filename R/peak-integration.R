

# Copyright 2020-2022 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms
# shewhart4hrms is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any 
# later version.
# 
# shewhart4hrms is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along 
# with shewhart4hrms. If not, see <https://www.gnu.org/licenses/>.


#' Pick peaks
#' 
#' @export
peakpicking_BfG_cpp <- function(i, rawData, mz_step, rt_min_scan, rt_max_scan, sn, int_threshold, 
                                NoiseScans, peakwidth_min, peakwidth_max, precursormzTol, maxPeaksPerSignal) {
  
  maxima <- NULL
  
  XIC <- xcms::rawEIC(rawData, mzrange = c(i,i+mz_step))
  XIC <- XIC$intensity
  
  maxima <- peakPickingBfGC(mz = i, mz_step = mz_step, XIC, scantime = rawData@scantime, 
                            min_intensity = int_threshold, sn = sn, noisescans = NoiseScans, 
                            peakwidth_min = peakwidth_min, peakwidth_max = peakwidth_max, 
                            maxPeaksPerSignal = maxPeaksPerSignal)
  #browser(expr = abs(i - 269.17) < 0.1)
  if (nrow(maxima) > 0) {
    for (j in 1:nrow(maxima)) {
      mass_spectrum <- xcms::getScan(rawData, maxima[j,5], mzrange = c(i,i+mz_step))
      exactmass <- 0
      if (nrow(mass_spectrum) > 0) {
        exactmass <- mass_spectrum[which.max(mass_spectrum[,2]),1]
        maxima[j,3] <- max(mass_spectrum[,2])
        if (nrow(mass_spectrum) == 1) maxima[j,3] <- maxima[j,3]-maxima[j,14]
      }
      if ((exactmass == 0) | (exactmass < i+mz_step/4) | (exactmass > i+mz_step/4*3)) exactmass = 0 
      maxima[j,1] <- exactmass
      ms2 <- which((rawData@msnRt > maxima[j,6]) & (rawData@msnRt < maxima[j,7]) & (abs(rawData@msnPrecursorMz-exactmass) <= exactmass/1000000*precursormzTol))
      if (length(ms2) > 0) maxima[j,16] <- ms2[which.min(abs(rawData@msnRt[ms2]-maxima[j,2]))]
    }
  }
  maxima <- maxima[maxima[,1] > 0,,drop = FALSE]
  if (nrow(maxima) > 0) return(maxima)
}

