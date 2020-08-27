function BP_mov = BP_filter(mov,f1,f2)
    %highpass
    LP_mov = movmean(mov,f1,3,'Endpoints','shrink');
    %lowpass
    HP_mov = movmean(mov,f2,3,'Endpoints','shrink');
    %take absol
    BP_mov = abs(LP_mov-HP_mov);
end