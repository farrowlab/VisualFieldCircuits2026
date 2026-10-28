function tot = input_totals(T, A)
%INPUT_TOTALS Input totals and baseline densities of one experiment.
%   tot = INPUT_TOTALS(T, A) uses the area table T of one experiment
%   (from count_cells_per_area) and the area indices A (from area_indices).
%   Midbrain summary rows (MBsen, MBmot, MBsta) and IC subdivisions are
%   excluded from all totals and volumes.
%
%   tot.nMBsen, nMBmot, nMBsta   cells in the midbrain summary rows
%   tot.nTotal                   all inputs
%   tot.nSC                      cells in the sSC (injection site)
%   tot.nTotalminusSC            inputs outside the sSC
%   tot.nCTX, nsCTX, nTH         cells in cortex, cortical subplate, thalamus
%   tot.nTotalminuCTX            inputs outside cortex, cortical subplate and sSC
%   tot.nMax                     cells in the area with most cells
%   tot.aSC, tot_vol             sSC volume, total volume
%   tot.aTotalminusSC            volume outside the sSC
%   tot.SCDensity                sSC cell density
%   tot.denTotalminusSC          baseline cell density outside the sSC
%   tot.denSum                   sum of area densities outside the sSC

count = cell2mat(T.GroupCount);
vol   = cell2mat(T.vol_all);

keep = A.notMBorIC;          % rows counted in totals
keepNoSC = keep;
keepNoSC(A.SC) = false;

tot.nMBsen        = sum(count(A.MBsen));
tot.nMBmot        = sum(count(A.MBmot));
tot.nMBsta        = sum(count(A.MBsta));
tot.nTotal        = sum(count(keep));
tot.nSC           = sum(count(A.SC));
tot.nTotalminusSC = tot.nTotal - tot.nSC;
tot.nCTX          = sum(count(A.CTX));
tot.nsCTX         = sum(count(A.sCTX));
tot.nTotalminuCTX = tot.nTotal - tot.nCTX - tot.nsCTX - tot.nSC;
tot.nTH           = sum(count(A.TH));
tot.nMax          = max(count(keep));

tot.aSC             = sum(vol(A.SC));
tot.SCDensity       = tot.nSC / tot.aSC;
tot.tot_vol         = sum(vol(keep));
tot.aTotalminusSC   = tot.tot_vol - tot.aSC;
tot.denTotalminusSC = tot.nTotalminusSC / tot.aTotalminusSC;
tot.denSum          = sum(count(keepNoSC) ./ vol(keepNoSC), 'omitnan');
end
