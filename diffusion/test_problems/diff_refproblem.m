%diff_refproblem_unix    sets up reference examples  (Unix version)
%   TIFISS scriptfile: LR; 24 January 2019.
%Copyright (c) 2016 by D.J. Silvester, Qifeng Liao
gohome;
clear;
close all;
fprintf('\nspecification of reference diffusion problem.');
fprintf('\nchoose specific example');
fprintf('\n     1  Square domain, unit coefficients');
fprintf('\n     2  Square domain, anisotropic coefficients');
fprintf('\n     3  Circle domain, unit coefficients, analytic solution');
fprintf('\n     4  punched ticket domain, unit coefficients');
fprintf('\n     5  L-shaped domain, constant load function');
fprintf('\n');

sn = default('',1);
%
if sn==1
    !/bin/cp ./diffusion/test_problems/unit_rhs.m            ./diffusion/specific_rhs.m
    !/bin/cp ./diffusion/test_problems/zero_bc.m             ./diffusion/specific_bc.m
    !/bin/cp ./diffusion/test_problems/unit_adiff.m          ./diffusion/specific_adiff.m
    !/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m     ./diffusion/specific_gradcoeff.m
    square_adiff;
elseif sn==2
    !/bin/cp ./diffusion/test_problems/unit_rhs.m            ./diffusion/specific_rhs.m
    !/bin/cp ./diffusion/test_problems/zero_bc.m             ./diffusion/specific_bc.m
    !/bin/cp ./diffusion/test_problems/anisotropic_adiff.m   ./diffusion/specific_adiff.m
    !/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m     ./diffusion/specific_gradcoeff.m
    square_adiff;
elseif sn==3
    !/bin/cp ./diffusion/test_problems/unit_rhs.m            ./diffusion/specific_rhs.m
    !/bin/cp ./diffusion/test_problems/circle_bc.m           ./diffusion/specific_bc.m
    !/bin/cp ./diffusion/test_problems/unit_adiff.m          ./diffusion/specific_adiff.m
    !/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m     ./diffusion/specific_gradcoeff.m
    circle_adiff; 
elseif sn==4
    !/bin/cp ./diffusion/test_problems/unit_rhs.m            ./diffusion/specific_rhs.m
    !/bin/cp ./diffusion/test_problems/zero_bc.m             ./diffusion/specific_bc.m
    !/bin/cp ./diffusion/test_problems/unit_adiff.m          ./diffusion/specific_adiff.m
    !/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m     ./diffusion/specific_gradcoeff.m
    obstacle_adiff;
elseif sn==5
    !/bin/cp ./diffusion/test_problems/unit_rhs.m            ./diffusion/specific_rhs.m
    !/bin/cp ./diffusion/test_problems/zero_bc.m             ./diffusion/specific_bc.m
    !/bin/cp ./diffusion/test_problems/unit_adiff.m          ./diffusion/specific_adiff.m
    !/bin/cp ./diffusion/test_problems/zeros_gradcoeff.m     ./diffusion/specific_gradcoeff.m
    ell_adiff;
else
    error('Oops.. reference problem datafile not found!');
end