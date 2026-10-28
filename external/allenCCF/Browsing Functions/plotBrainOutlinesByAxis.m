
function plotBrainOutlinesByAxis(av,axH,LineColor, LineWidth)

    
    axes(axH);
    
    %sliceIm = squeeze(max(av,[],v)>1);
    sliceIm = av;
    
    c = contourc(double(sliceIm), [3 3]);
    coordsReg = makeSmoothCoords(c);
    for cidx = 1:numel(coordsReg)
        h = plot(axH, coordsReg(cidx).x,coordsReg(cidx).y, 'Color', LineColor, 'LineWidth', LineWidth); hold on;
        
    end
    
    set(axH, 'YDir', 'reverse'); 
    axis equal


%view([90 -90]); % sagittal should be flipped
set(axH, 'XDir', 'reverse'); 
box off;

end