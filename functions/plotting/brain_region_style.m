function st = brain_region_style()
%BRAIN_REGION_STYLE Names and colours of the major brain regions in the figures.
%   st = BRAIN_REGION_STYLE() returns, indexed by the parent group in
%   column 6 of BBB (1 CTX, 2 sCTX, 3 TH, 4 HYP, 5 MB, 6 HB, 7 CB):
%     st.names       region labels
%     st.textColor   colour of the region labels
%     st.shadeColor  colour of the background shading
%     st.shadeAlpha  transparency of the background shading

st.names      = ["CTX", "sCTX", "TH", "HYP", "MB", "HB", "CB"];
st.textColor  = [0.8784 0.5020 0.1490; 0.4314 0.2118 0.0471; 0.8784 0.0353 0.1490; ...
                 0.6392 0.0588 0.4471; 0.196 0.533 0.741; 0.1451 0.7216 0.0588; 0.312 0.145 0.441];
st.shadeColor = [0.992 0.682 0.380; 0.8706 0.7608 0.6784; 0.861 0.271 0.349; ...
                 0.6392 0.1647 0.4824; 0.196 0.533 0.741; 0.5176 0.8196 0.4745; 0.312 0.145 0.441];
st.shadeAlpha = [0.3 0.75 0.22 0.25 0.25 0.25 0.2] - 0.1;
end
