function [] = diffusion_stats()

  fp=[];
  tbc=[];
  dt=[];
  v=[];

  for i=1:50
    # Collision points
    pc_x = load(["diffusion-results/", "collision-position-", num2str(i), "-x.tsv"]);
    pc_y = load(["diffusion-results/", "collision-position-", num2str(i), "-y.tsv"]);
    pc_z = load(["diffusion-results/", "collision-position-", num2str(i), "-z.tsv"]);
    
    # Compute mean free path
    n = length(pc_x);
    
    for j=2:n
      fp(end+1)=sqrt((pc_x(j,2)-pc_x(j-1,2))^2+(pc_z(j,2)-pc_z(j-1,2))^2+(pc_z(j,2)-pc_z(j-1,2))^2);
      dt(end+1)=pc_x(j,1)-pc_x(j-1,1);
      v(end+1)=fp(end)/dt(end); # velocity
    endfor
    
    # Position sampled every 0.1 ns
    p_x = load(["diffusion-results/", "position-", num2str(i), "-x.tsv"]);
    p_y = load(["diffusion-results/", "position-", num2str(i), "-y.tsv"]);
    p_z = load(["diffusion-results/", "position-", num2str(i), "-z.tsv"]);
    
    # Compute mean square displacement
    
    # Collision times
    ct = load(["diffusion-results/", "collisions-", num2str(i), ".tsv"]);
    
    # Compute time between collisions
    n = length(ct);
    
    for j=2:n
      tbc(end+1)=ct(j)-ct(j-1);
    endfor
 
  endfor
  
  # velocity histogram 
  n=length(v);
  k=ceil(2*n^(1/3)); # Rice rule
  [f, x]=hist(v, k); 
  bar(x, f);
  X(:,1)=x;
  X(:,2)=f;
  dlmwrite('calcium-velocity-histogram-25c.tsv',X,"\t");
  mean_free_path=mean(fp)
  mean_time_between_collisions=mean(tbc)

end