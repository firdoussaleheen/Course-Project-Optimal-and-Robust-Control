function c = Hmatrix(Ft, k, Ntau, Qf, G, W, tau)

    syms t tau zetaa D H


    x0 = 1;         %% initial conditions, n x 1 vector
    t0 = 0;         %% beginning integration time
    tf = 1;         %% ending integration time
    
            %%% INT(S,v,a,b) is the definite integral of S with respect to v
            %%% from a to b.
            
    phi_tau  = exp(int(Ft, t, t0, tau));    %%% phi in 
    phi_tf   = exp(int(Ft, t, t0, tf));
    Nzeta    = subs(Ntau, zetaa);
    phi_zeta = exp(int(Ft, t, t0, zetaa));
    c    = cell(k);
    
    for l = 1:k
      
        Kx6 = kayx(tau, k, t0, tf);     % Kx6
        Kx2 = subs(Kx6{l+1}, t, tau);   % Kx2
        Kx  = subs(Kx2, tau, zetaa);     % Kx  x
        Kx3 = subs(Kx6{l+1}, t, tf);    % Kx3
        Kx4 = subs(Kx3, tau, tf);       % Kx4 x
        Kx5 = subs(Kx3, tau, zetaa);     % Kx5 x
        Kx1 = subs(Kx2, tau, tf);       % Kx1 x
   % K_tau_zeta
   % K_tf_tf
   % K_tf_zeta
   % K_tau_tf
   % K_tau_tau
        if l==1
                    %%% H = (**24) or (*1.27)
            H = phi_tf'*Qf*phi_tf + int(phi_tau'*Ntau*phi_tau, tau, t0, tf); 
            
        else        %%% H = (*1.28)
            first_int  = phi_tau'*Ntau*Kx*Nzeta*phi_zeta;
            second_int = int(first_int, zetaa, t0, tf);  %%%dzeta
            H = int(second_int, tau, t0, tf) + phi_tf'*Qf*Kx4*Qf*phi_tf + ... %%% dtau
                phi_tf'*Qf*int(Kx5*Nzeta*phi_zeta, zetaa, t0, tf) + ...    
                int(phi_tau'*Ntau*Kx1, tau, t0, tf)*Qf*phi_tf;
        end %% end if/else
        
        D = int(trace(subs(H, t0, zetaa)*G*W*G'), zetaa, t0, tf);     %%% (**26)
      
        c{l,1} = H;
        c{l,2} = D;
    end %% for

end