mystartdefaults

dispersion_relation=true; compute_dos=false;
CohenBergstresser1966
tic

save_path = [pwd(), '\figures'] ##save figures here
if ~exist(save_path, 'dir') ##If doesnt exist then mkdir
    mkdir(save_path);
end

recipunit =  1.0e+10;
ekinunit = ((hbar*recipunit)^2/(2*elm))/qel;
m=14;
semiconductor = material_list(m,1:4);

[qpath,tix,til]=BZpath(BZstep,qs,qe,qs_str,qe_str);

#To calculate reciprocal unit lattice vectors
g = zeros(4,3);
g(1:3,1) = cross(a(:,2),a(:,3))/cell_volume;
g(1:3,2) = cross(a(:,3),a(:,1))/cell_volume;
g(1:3,3) = cross(a(:,1),a(:,2))/cell_volume;
for i=1:3;
  g(4,i) = g(1:3,i)'*g(1:3,i);
end


#Build reciprocal lattice
min_norm = sqrt(min(g(4,:)));
limit = ceil(sqrt(sqrt(cutoff))/min_norm);
fprintf('Cutoff condition requires %d positive steps along reciprocal unit lattice vectors\n',limit);

nodes = (2 * limit + 1)^3;
fprintf('Generate (2* %1d + 1)^3 = %4d reciprocal lattice vectors\n',limit,nodes);

G = zeros(5,nodes);
i = 0;
for h=-limit:limit;
  for k=-limit:limit;
    for l=-limit:limit;
      i+=1;
      G(1:3,i) = h*g(1:3,1) + k*g(1:3,2) + l*g(1:3,3);
      G(4,i) = G(1:3,i)'*G(1:3,i);
      G(5,i) = sqrt(G(4,i));
    endfor
  endfor
end

#Sorting the G vectors according to the norm
[G(5,:),perm] = sort(G(5,:));
G(1:4,:) = G(1:4,perm);

#Finding number of G points such that |G|^2<Gs_max
n = 1;
for i=2:length(G);
  if G(4,i)<=cutoff;
    n++;
  endif
end

ngx = ceil(n/2);
n = 2*ngx-1;   #size of hamiltonian

#Defining the pseudopotential from form factors of the semiconductor
V = zeros(1,n);
a = ls(1,m);               #Lattice spacing

fprintf('\n            G(1)            G(2)            G(3)            G(4)            V\n')

ff(:,1) = [-0.7724,-0.6951,-0.5015,-0.6775,-0.6527,-0.5107,-0.5627,-0.5112,-0.5251,-0.4460,-0.4677,-0.4496,-0.3925,-0.3117]';
for i=1:n;
  Vs = 0;
  Va = 0;
  if G(4,i)<=Gs_max;
    if G(4,i)==0;
      Vs = ff(m,1) * Rydberg;
    endif
    if G(4,i)==3;
      Vs = ff(m,2) * Rydberg;
      Va = ff(m,5) * Rydberg;
    endif
    if G(4,i)==4;
      Va = ff(m,6) * Rydberg;
    endif
    if G(4,i)==8;
      Vs = ff(m,3) * Rydberg;
    endif
    if G(4,i)==11;
      Va = ff(m,7) * Rydberg;
      Vs = ff(m,4) * Rydberg;
    endif
  endif
  V(1,i) = Vs * cos(2*pi*G(1:3,i)'*tau(1:3,1)) - 1i*Va * sin(2*pi*G(1:3,i)'*tau(1:3,2));  #eqn 3.107
  fprintf('%15.6G %15.6G %15.6G %15.6G %15.6G\n',[G(1,i),G(2,i),G(3,i),G(4,i),V(1,i)]);
end


#Defining Hamiltonian matrix
H = zeros(n,n);
#Adding potential terms in Hamiltonian
for i=1:n;
  for j=1:n;
    for k=1:n;
      if norm(G(1:3,i)-G(1:3,j)-G(1:3,k))<tol
        H(i,j) = V(1,k);
        break
      endif
    endfor
  endfor
end

#Adding Kinetic energy terms
for k=1:length(qpath)
  for i=1:n
    H(i,i) = V(1) + ekinunit*((2*pi/a)^2) * ((qpath(1:3,k)-G(1:3,i))' * (qpath(1:3,k)-G(1:3,i)));  #eqn 3.21
  endfor

  if(!ishermitian(H,tol))
    printf('\nHamiltonian matrix not hermitian : fatal error.\n');
    return;
  else
    [v,ev] = eig(H);
    E = real(diag(ev));
    [E,perm] = sort(E);
    v = v(:,perm);

    for g=1:nband
      bandstructure(g,k) = E(g);
    endfor
  endif
end

#num_k_points = 1:length(qpath);
graph = plot(qpath(5,:),bandstructure,'-', 'LineWidth', 1);
hold on;
graph_markers = plot(qpath(5, 1:5:end), bandstructure(:, 1:5:end), 'o', 'MarkerSize', 1.2, 'Color', 'k','MarkerFaceColor', 'w');

xlim([0,max(qpath(5,:))]);
ylim([-4,8]);
set(gca, 'ytick', [-16:1:24], 'FontSize',8,'FontWeight','bold')
set(gca,'xtick',tix);
set(gca,'xticklabel',til,'FontSize',10,'FontWeight','bold');
set(gca, 'XTickLabelRotation', 45);
xlabel ('k [{{2\pi}/{a}}]', 'Fontsize', 12,'FontWeight','bold')
ylabel('E [eV]', 'Fontsize', 12,'FontWeight','bold');
grid();
title(semiconductor, 'FontSize',12,'FontWeight','bold');
##cellarr = strcat(save_path, '\', semiconductor,'_16_bands.png');
##filename = strcat(cellarr, '');
##print(filename, '-dpng', '-r300', '-tight');
waitfor(graph);
toc


