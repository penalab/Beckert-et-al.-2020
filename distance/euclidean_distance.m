function e_dis = euclidean_distance(spt_1, spt_2, bins)

    e_dis = sqrt( sum( (hist(spt_1, bins) - hist(spt_2, bins)) .^2 ) );
    
end