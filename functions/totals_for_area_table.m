function D = totals_for_area_table(OD, tableid)
%TOTALS_FOR_AREA_TABLE Collect per-experiment totals for build_area_table.
%   D = TOTALS_FOR_AREA_TABLE(OD, tableid) gathers the totals computed by
%   input_totals (stored as 1-by-N cell fields of OD) for the experiments
%   tableid, in the form expected by build_area_table.

D.nTotal          = [OD.nTotal{tableid}];
D.nTotalminusSC   = [OD.nTotalminusSC{tableid}];
D.nMax            = [OD.nMax{tableid}];
D.nTotalminusCTX  = [OD.nTotalminuCTX{tableid}];
D.nSC             = [OD.nSC{tableid}];
D.SCDensity       = [OD.SCDensity{tableid}];
D.denTotalminusSC = [OD.denTotalminusSC{tableid}];
D.denSum          = [OD.denSum{tableid}];
end
