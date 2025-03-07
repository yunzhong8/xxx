%SQUARE_ADIFF   solve variable diffusion problem in unit square domain
%   TIFISS scriptfile: DJS; 5 March 2017. (LR modification 16 October 2017; DJS 19 December 2017)
% Copyright (c) 2017 D.J. Silvester, Alex Bespalov, Qifeng Liao

  gohome, cd datafiles
  save currentsn sn
  clear variables
  load currentsn

% Define geometry
  pde = 1; domain = 11;

  if sn==1
      square_domain(1,1);
  elseif sn==2 
      square_domain(1,0);
  elseif sn==9
      square_domain(2,1);
  else
      square_domain(1,0);
  end
  load square_grid

  [evt,eboundt] = p1grid(xy,mv,bound,mbound);

% Choose red or bisec3 uniform sub-division for spatial estimation 1/2
  subdivPar = 2;

% Uniform refinement  
  nref = default('refinement level (default 0)',0);
  for i = 1:nref
      [xy,evt,bound,eboundt] = p1_refinement(xy,evt,bound,eboundt);
  end
zzq_p1_xy=xy;zzq_p1_evt = evt; zzq_p1_bound = bound;zzq_p1_ebound=eboundt;
save('zzq.mat','zzq_p1_xy','zzq_p1_evt','zzq_p1_bound','zzq_p1_ebound');
%********************************************finite element caculate question************************** 
  pmethod = default('P1/P2 approximation 1/2? (default P1)',1);
  if pmethod == 1       
      [A,M,f] = femp1_adiff(xy,evt); 
  elseif pmethod == 2
      % check for nonconstant coefficients
      [dcoeffdx,dcoeffdy] = specific_gradcoeff(x,y,size(evt,1));
      if norm(dcoeffdx)+ norm(dcoeffdx)~=0
          error('Oops.. P2 with variable coefficients is not allowed!');
      end
      xyp1 = xy;
      boundp1 = bound;
      [xy,evt,bound] = p2grid(xy,evt,bound,eboundt,0);
      [A,M,f] = femp2_diff(xy,evt,xy,evt);
  end

% Boundary conditions
  [Agal,fgal] = nonzerobc(A,f,xy,bound);

% Save resulting system
  fprintf('system saved in system_adiff.mat ...\n')
  gohome
  cd datafiles
  save system_adiff pmethod Agal M  fgal  xy

% Compute solution
  tic
  fprintf('solving linear system ...  ');
  x_gal = Agal\fgal;   
  fprintf('done\n');
  etoc = toc; fprintf('Galerkin system solved in  %8.3e seconds\n',etoc); 
  save system_adiff x_gal  -append


if pmethod == 1
    fprintf('question2 P1 error');
    p1_error=si214c_question2_error;
    fprintf('\n\n******p1 Error is %.6f******\n',p1_error);
elseif pmethod == 2
end
%*******************************************************use result of
%compution plot figure

% Plot solution - refine grid
  xmin = min(xy(:,1)); xmax = max(xy(:,1));
  ymin = min(xy(:,2)); ymax = max(xy(:,2));
  x = linspace(xmin,xmax,301);
  y = linspace(ymin,ymax,301);
  tsolplot(pmethod,x_gal,evt(:,1:3),xy,x,y,12,[xmin,xmax,ymin,ymax]);
  drawnow;


% Compute a posteriori error estimation 
% -------------------------------------------------------------------

% Validate mesh
  [eex,tve,els] = tedgegen(xy,evt);

  if pmethod == 1      
      fprintf('\n');
      pestim = default('Error estimation: linear/quadratic bubble functions 1/2 (default 1)',1);  
      if pestim == 1
      % Linear midpoint hat functions
          fprintf('<strong>diffpost using 3 edge midpoint linear functions...</strong>\n'); 
          tic; 
          [elerr_p,fe,ae] = diffpost_p1_with_p1(xy,evt,eex,tve,els,eboundt,x_gal,subdivPar);
          toc;
          [err_p,elerr_p] = diffpost_p1_bc(ae,fe,elerr_p,xy,evt,eboundt);
      elseif pestim == 2
          % Quadratic midpoint bubble functions
          fprintf('\n<strong>diffpost using 4 quadratic bubble functions...</strong>\n');
          tic;
          [elerr_p,fe,ae] = diffpost_p1_with_p2(xy,evt,eex,tve,els,eboundt,x_gal); 
          toc;    
          [err_p,elerr_p] = diffpost_p1_bc(ae,fe,elerr_p,xy,evt,eboundt);
      else
          error('Value not valid! Choose between either 1 or 2!');
      end
      fprintf('<strong>Estimated energy error: %10.4e</strong>\n\n',norm(elerr_p,2));
      terrplot(pmethod,x_gal,elerr_p,evt,xy,x,y,101,[xmin,xmax,ymin,ymax]);

  elseif pmethod == 2
      % check performance with that of legacy code
      %tic,[elerr_x_p4,error_total_p4,fe_p4,ae_p4] = ...
      %    diffpost_p2_with_p4_x(xy,evt,eboundt,x_gal); toc
      tic;
      [elerr_p_p4,fe,ae] = diffpost_p2_with_p4(xy,evt,eex,tve,els,eboundt,x_gal); 
      toc
      %fprintf('sanity check %10.4e\n',norm(elerr_p_p4-elerr_x_p4,inf))
      terrplot(pmethod,x_gal,elerr_p_p4,evt,xy,x,y,102,[xmin,xmax,ymin,ymax]);

  end % end if

% end script  
