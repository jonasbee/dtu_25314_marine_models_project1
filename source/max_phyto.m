function max_P = max_phyto(P)
    max_val = max(P(end,:));
    [row, col] = find(P(end,:)==max_val);
    max_P = [max_val col];
end