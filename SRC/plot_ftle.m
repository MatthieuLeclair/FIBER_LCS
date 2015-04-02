function fig = plot_ftle(X, Y, FTLE, ttl, visibility)
   
   DXcm = 10;
   DYcm = DXcm * (max(max(Y)) - min(min(Y))) / (max(max(X)) - min(min(X)));
   MLcm = 1 ;  MRcm = 2 ;
   MBcm = .7;  ; MTcm = .7;
   
   CBMLcm = .5; CBWcm = .3;
   
   xSize = MLcm + DXcm + MRcm;
   ySize = MBcm + DYcm + MTcm;
   
   ML = MLcm/xSize; DX = DXcm/xSize; MR = MRcm/xSize; CBML = CBMLcm/xSize; CBW = CBWcm/xSize;
   MB = MBcm/ySize; DY = DYcm/ySize; MT = MTcm/ySize; 
   
   fig = figure('Visible'      ,  visibility      , ...
                'PaperUnits'   , 'centimeters'    , ...
                'PaperPosition', [0 0 xSize ySize], ...
                'PaperSize'    , [xSize,ySize]           );
   
   curax = axes('Position', [ML ; MB ; DX ; DY] );
   imagesc( X, Y, FTLE )
   set(curax,'YDir','normal')
   
   clb = colorbar('Location' , 'EastOutside');
   set( curax, 'Position', [ML ; MB ; DX ; DY] )
   set( clb  , 'Position', [ML+DX+CBML; MB; CBW; DY] )
   caxis([0 0.4])
   title(ttl);

end