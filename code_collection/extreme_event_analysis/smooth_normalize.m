function f_out = smooth_normalize(f, f_smoothed)
    
   var1 = var(f_smoothed(:)); 
   var2 = var(f(:));
   mean1 = mean(f_smoothed(:));
   mean2 = mean(f(:));
   f_out = (f_smoothed - mean1) / sqrt(var1) * sqrt(var2) + mean2;

end



