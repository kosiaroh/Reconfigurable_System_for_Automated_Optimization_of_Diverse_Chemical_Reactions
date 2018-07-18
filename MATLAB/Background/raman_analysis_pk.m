[node, ramanshift, Y, peak, modulus] = textread(raman_file, 'headerlines', 15) ;
diketone = trapz(peak(1074:1085));
product = trapz(peak(988:1007));
diketone_solv = 1; 
product_solv = 1;
product = product/product_solv;
diketone = diketone/diketone_solv;