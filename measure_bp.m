function [S, N]=measure_bp(Phi_T, phi_old, B_up, B_dn, Proj_k_half, t_bp, t_pop, N_up, N_par, N_sites, Lx, Ly, N_x, N_y)
  %% Implement backwards propagation with the stored values
  % calculate backpropagated bra
  Phi_left=Phi_T;
    for ii=1:t_bp+1
        Phi_left=Proj_k_half*Phi_left;
        for jj=1:N_sites
            Phi_left(jj,1:N_up)=Phi_left(jj,1:N_up)*B_up(t_bp-ii+2,jj);
            Phi_left(jj,N_up+1:N_par)=Phi_left(jj,1+N_up:N_par)*B_dn(t_bp-ii+2,jj);
        end
        Phi_left=Proj_k_half*Phi_left;
        if mod(ii,t_pop)==0
            [Phi_left(:,1:N_up), Trsh]=qr(Phi_left(:,1:N_up),0);
            [Phi_left(:,1+N_up:N_par), Trsh]=qr(Phi_left(:,1+N_up:N_par),0);
        end
    end
    %% Once the bra is back propagated, the observable is calculated the usual way
    % calculate inverse matrices in order to calculate Green's function
    inv_O_up=inv(Phi_left(:,1:N_up)'*phi_old(:,1:N_up));
    inv_O_dn=inv(Phi_left(:,1+N_up:N_par)'*phi_old(:,1+N_up:N_par));
    temp_up=phi_old(:,1:N_up)*inv_O_up;
    temp_dn=phi_old(:,N_up+1:N_par)*inv_O_dn;
    G_up=temp_up*Phi_left(:,1:N_up)';
    G_dn=temp_dn*Phi_left(:,N_up+1:N_par)';
    Id=eye(N_sites);
    G_up_tilde=Id-G_up;
    G_dn_tilde=Id-G_dn;
    % As second try I will measure the spin correlation
    aux1=0;
    aux2=0;
    qx=linspace(-pi,pi,N_x);
    qy=linspace(-pi,pi,N_y);
    ind=0;
    S=zeros(N_x*N_y,1);
    N=zeros(N_x*N_y,1);
    
    Corr=zeros(N_sites,N_sites);
    Charg=zeros(N_sites,N_sites);
    for k=1:N_sites
        for l=1:N_sites
            [Corr(k,l),Charg(k,l)]=measure_corr(k,l,G_up, G_up_tilde, G_dn, G_dn_tilde);
        end
    end
    for i=1:N_x
        for j=1:N_y
            ind=ind+1;
            for k=1:N_sites
                for l=1:N_sites
                    ix1=floor(k/Lx);
                    iy1=mod(k,Ly);
                    if iy1==0
                        iy1=Ly;
                    end
                    ix2=floor(l/Lx);
                    iy2=mod(l,Ly);
                    if iy2==0
                        iy2=Ly;
                    end
                    S(ind)=S(ind)+exp(sqrt(-1)*((ix1-ix2)*qx(i)+(iy1-iy2)*qy(j)))*Corr(k,l);
                    N(ind)=N(ind)+exp(sqrt(-1)*((ix1-ix2)*qx(i)+(iy1-iy2)*qy(j)))*Charg(k,l);
                end
            end
            S(ind)=S(ind)/N_sites;
            N(ind)=N(ind)/N_sites;
        end
    end
end