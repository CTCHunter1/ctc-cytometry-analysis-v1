function shiftXTicks(dy)
xtickpos = get(gca, 'Xtick');
ylimvals = get(gca, 'Ylim');
set(gca, 'XTickLabel', []);

    for ii = 1:length(xtickpos);
    t = text(xtickpos(ii), ylimvals(1), num2str(xtickpos(ii)));
    set(t, 'Units', 'Normal');
    set(t, 'Position', get(t, 'Position') + [0, dy, 0] -[0, .18, 0] );
    %set(t, 'Units', 'pixels');
    t1 = get (t, 'Position');    
    t1ex = get(t, 'Extent');
    set(t, 'Position', [t1(1) - t1ex(3)/4, t1(2), t1(3)]);
    end
    
end