function [histw] = histcounts_w(vv, ww, edges)
  histw = zeros(1,length(edges));
  for i=1:length(vv)
    ind = find(edges < vv(i), 1, 'last' );
    if ~isempty(ind)
      histw(ind) = histw(ind) + ww(i);
    end
  end
  histw = histw./sum(histw);
end