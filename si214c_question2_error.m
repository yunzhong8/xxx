function error=si214c_question2_error
gohome,cd datafiles
data = load ('system_adiff.mat');
zzq_p1_x_gal=data.x_gal;
zzq_p1_xy=data.xy;
zzq_p1_m=data.M;
zzq_p1_x = zzq_p1_xy(:,1);
zzq_p1_y = zzq_p1_xy(:,2);
u = @(x, y) sin(4 * pi * x) .* sin(4 * pi * y);
zzq_p1_exact_x_gal=u(zzq_p1_x,zzq_p1_y);
zzq_p1_xy_error = zzq_p1_exact_x_gal - zzq_p1_x_gal; 
error = zzq_p1_xy_error'* zzq_p1_m * zzq_p1_xy_error;
error = sqrt(error);

end
