#altern_psatz.jl
#certify that a polynomial q(A, B) is nonnegative over the L2-specified consistency set
#use a theorem of alternatives to get a simpler program

mutable struct psatz
    mu      #equality constrained multipliers
    quad_x  #x: polynomials (s, tau) from quad_mult    
    quad_u  #u: polynomials (s, tau) from quad_mult   
    blocksize #sizes of PSD matrices involved in the constraint 
end

mutable struct quad_mult
    s       #poly vector inside the SOC
    tau     #poly vector for the SOC radius
    Gram    #Gram matrices in the certificate
    blocksize #size of the Gram matrices
end

function make_mult_quad(n, model, vars, order, SPARSE=true)
    quad_out = 0;
    if SPARSE
        quad_out = make_mult_quad_sparse(n, model, vars, order);        
    else
        quad_out = make_mult_quad_dense(n, model, vars, order);
    end

    return quad_out;
end

function make_mult_quad_dense(n, model, vars, order)
    #make_mult_quad_dense: create the multipliers associated with a single robust SOC constraint:
    #
    #[tau s\\ s I*tau] in SOS^(n+1) [A, B]. Each time index t has its own SOC constraint
    #    
    #
    #Input:
    #   model:  JuMP model which is containing all constraints
    #   vars:   the variables A and B
    #   order:  the order of polynomials to be used (degree/2)
    #
    #Output:    
    #   Gram:   Matrix that should be PSD
    #   s:      vector of polynomials
    #   z:      vector of polynomials

    #generate all monomials
    degree = 2*order;
    n = size(vars.A, 2)
    vars_flat = vec([vars.A vars.B]);
    # vars_c = col_vars(vars);
    mon = reverse(monomials(vars_flat, 0:order));

    #create the block Gram matrix
    len_gram = (n+1)*length(mon);
    nv = length(mon);
    
    #catch the polynomials and store them in arrays
    # Gram=Array{LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}}(undef, n, 1);
    # Zeta = Matrix{Polynomial{true, Float64}}}(undef)
    # s = Array{Polynomial}(undef, n, 1);
    # z = Array{Polynomial}(undef, n, 1);
    # zcorner = Array{Polynomial}(undef, n, 1);
    
    #the sum of the z elements is each entry of zcorner

        #declare the PSD matrix variable 
        Gram_curr = JuMP.@variable(model, [1:len_gram, 1:len_gram], PSD);
        Gram = Gram_curr;

        eye_n = 1* Matrix(I, n+1, n+1)
        mon_kron =  kron(mon, eye_n);        
        #index out the polynomial pieces
        # s[i] = mon'*Gram_s*mon;
        # z[i] = mon'*Gram_z*mon;
        Zeta = mon_kron'*Gram_curr*mon_kron;

        #fetch the multipliers
        s = Zeta[2:end, 1];
        tau = Zeta[1, 1];

        nv = length(mon);

        # @constraint(model, Gram_curr >= 0, PSDCone());
        for i = 1:n
            for j=i:n
                Zeta_curr = Zeta[i, j];                
                # zeta_coeff = coefficients(Zeta_curr);
                if i==j
                    #on-diagonal
                    @constraint(model, coefficients(Zeta_curr - tau)==0);
                else
                    Zeta_opp = Zeta[j, i];
                    #off-diagonal
                    @constraint(model, coefficients(Zeta_opp + Zeta_curr)==0);
                end
            end
        end



        blocksize = [len_gram];
        #add the constraints

    return quad_mult(vec(s), tau, Gram, blocksize)


end

function make_mult_psd_dense(model, vars, order)
    #make_mult_quad_dense: create the multipliers associated with a single robust SOC constraint:
    #
    #[Zeta] in SOS^(n+1) [A, B]. Each time index t has its own SOC constraint
    #    
    #
    #Input:
    #   model:  JuMP model which is containing all constraints
    #   vars:   the variables A and B
    #   order:  the order of polynomials to be used (degree/2)
    #
    #Output:    
    #   Gram:   Matrix that should be PSD
    #   s:      vector of polynomials
    #   z:      vector of polynomials

    #generate all monomials
    degree = 2*order;
    n = size(vars.A, 2)
    vars_flat = vec([vars.A vars.B]);
    # vars_c = col_vars(vars);
    mon = reverse(monomials(vars_flat, 0:order));

    #create the block Gram matrix
    len_gram = (n+1)*length(mon);
    nv = length(mon);
    
    #catch the polynomials and store them in arrays
    # Gram=Array{LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}}(undef, n, 1);
    # Zeta = Matrix{Polynomial{true, Float64}}}(undef)
    # s = Array{Polynomial}(undef, n, 1);
    # z = Array{Polynomial}(undef, n, 1);
    # zcorner = Array{Polynomial}(undef, n, 1);
    
    #the sum of the z elements is each entry of zcorner

        #declare the PSD matrix variable 
        Gram_curr = JuMP.@variable(model, [1:len_gram, 1:len_gram], PSD);
        Gram = Gram_curr;

        eye_n = 1* Matrix(I, n+1, n+1)
        mon_kron =  kron(mon, eye_n);        
        #index out the polynomial pieces
        # s[i] = mon'*Gram_s*mon;
        # z[i] = mon'*Gram_z*mon;
        Zeta = mon_kron'*Gram_curr*mon_kron;

    return quad_mult_dense(Gram, Zeta)


end

function make_mult_quad_sparse(n, model, vars, order)
    #make_mult_quad_sparse: create the multipliers associated with a single robust SOC constraint:
    #
    #[sum_i z_it, s_it; s_it, z_it] in SOS^2[A, B]. Each time index t has its own SOC constraint
    #    
    #
    #Input:
    #   n:      number of dimensions of SOC constraint
    #   model:  JuMP model which is containing all constraints
    #   vars:   the variables A and B
    #   order:  the order of polynomials to be used (degree/2)
    #
    #Output:    
    #   Gram:   Matrix that should be PSD
    #   s:      vector of polynomials
    #   z:      vector of polynomials

    #generate all monomials
    vars_flat = vec([vars.A vars.B]);
    mon = reverse(monomials(vars_flat, 0:order));

    #create the block Gram matrix
    len_gram = 2*length(mon);
    nv = length(mon);
    
    #catch the polynomials and store them in arrays
    Gram=Array{LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}}(undef, n, 1);
    s = Array{Polynomial}(undef, n, 1);
    z = Array{Polynomial}(undef, n, 1);
    zcorner = Array{Polynomial}(undef, n, 1);
    
    #the sum of the z elements is each entry of zcorner
    tau = 0;
    for i = 1:n
        #declare the PSD matrix variable 
        Gram_curr = JuMP.@variable(model, [1:len_gram, 1:len_gram], PSD);
        Gram[i] = Gram_curr;

        #store the Gram matrix pieces
        Gram_corner = Gram_curr[1:nv, 1:nv]
        Gram_s = Gram_curr[1:nv, (nv+1):end];
        Gram_z = Gram_curr[(nv+1):end, (nv+1):end];

        #index out the polynomial pieces
        s[i] = mon'*Gram_s*mon;
        z[i] = mon'*Gram_z*mon;
        zcorner[i] = mon'*Gram_corner'*mon;

        tau = tau + z[i];        
    end

    #now iterate back through the corners and set the corners equal to the sum
    for i = 1:n
        # corner_diff = tau - zcorner[i];
        # corner_coeff = coefficients(corner_diff, monomials(corner_diff));
        # @constraint!(model, corner_coeff==0);
        @constraint(model, coefficients(tau - zcorner[i])==0);
    end

    blocksize = Int64.(ones(n, 1))*len_gram;

    return quad_mult(vec(s), tau, Gram, blocksize)


end

function quad_psatz(q, order, model, data, vars, SPARSE=false)
#enforce that q>=0 over the consistency set described in data
# (assume only x-noise for now)
    #create the multipliers 
    n = size(data.X, 1);
    m = size(data.U, 1);
    T = size(data.U, 2);

        

    vars_flat = vec([vars.A vars.B]);

    blocksize = [];

    #create the residual
    Xnext = data.X[:, 2:end];
    Xprev = data.X[:, 1:end-1];
    U = data.U[:, 1:end-1];
    h0 = Xnext - vars.A*Xprev - vars.B*U;

    #the contribution such that q==psatz_term
    psatz_term = 0;

    #form quad_x (quadratic inequality constraints for dx) 
    if data.epsilon[1] > 0
        quad_x = Array{quad_mult}(undef, T, 1);
        for k = 1:T        
            quad_x[k] = make_mult_quad(n, model, vars, order, SPARSE);  
            blocksize = [blocksize; quad_x[k].blocksize];
        end
    else
        quad_x = zeros(T, 1);
    end

    #form quad_u (quadratic inequality constraints for du)
    if data.epsilon[2] > 0
        quad_u = Array{quad_mult}(undef, T-1, 1);
        for k = 1:T-1        
            quad_u[k] = make_mult_quad(m, model, vars, order, SPARSE);        
            blocksize = [blocksize; quad_u[k].blocksize];
        end
    else
        quad_u = zeros(T-1, 1);
    end

    
    #form mu (process noise multipliers for w)
    mu = Array{quad_mult}(undef, n, T-1);
    if data.epsilon[3] > 0
        #yes process noise    
        for k = 1:(T-1)
            mu[k] = make_mult_quad(n, model, vars, order-1, SPARSE);        
            blocksize = [blocksize; mu[k].blocksize];
            
            psatz_term = psatz_term + sum(vec(mu[k].s).*h0[:, k]);
            psatz_term = psatz_term + mu[k].tau*data.epsilon[3]; #check the sign, might be -tau
        end
        
    else
        #no process noise        
        for k = 1:(T-1)
            s_curr = Array{Polynomial}(undef, n, 1);
            

            # mu[k].tau = 0;
            for j = 1:n
                # mon = reverse(monomials(vars_flat, 0:(2:order-1)));
                # coeff_curr = @variable(model, 1:n, [1:length(mon)])
                # mu[k].s = coeff_curr'*mon;
            

                # mucurr, mucurrc, mucurrb = add_poly!(model, vars_flat, 2*order-1);
                
                s_curr[j] = make_poly(model, vars_flat, 2*order-1);
            end   
            mu[k] = quad_mult(vec(s_curr), 0, [], []);         
            psatz_term = psatz_term + sum(vec(s_curr).*h0[:, k]);
            # v, vc, vb = add_poly!(model, x, 2d)
            # psatz_term = psatz_term + mu[:, k]'*h0[k];
        end
        # psatz_term = sum(mu .* h0);
    end
    
    #now iterate through the constraints (for x)
    s_sig_x = sqrt(data.Sigma[1]);
    s_sig_u = sqrt(data.Sigma[2]);
    s_sig_w = sqrt(data.Sigma[3]); #TODO: figure out how s_sig_w enters
    for k = 1:T
        
        if data.epsilon[1] > 0
            s_term_x = -s_sig_x \ quad_x[k].s;
            mu_con_x = mu_eq_con(T, k, mu, vars.A);

            for i=1:n
                @constraint(model, coefficients(s_term_x[i] - mu_con_x[i])==0);
            end

            psatz_term = psatz_term + quad_x[k].tau*data.epsilon[1];    
        end

        if (data.epsilon[2] > 0) & (k<= T-1)
            s_term_u = -s_sig_u \ quad_u[k].s;
            mu_con_u = vars.B'*mu[k].s;

            for j = 1:m
                @constraint(model, coefficients.(s_term_u[j] - mu_con_u[j])==0);
            end
            psatz_term = psatz_term + quad_u[k].tau*data.epsilon[2];    
        end
                
    end


    #seal up with sigma0 (the psatz)
    mon = reverse(monomials(vars_flat, 0:order));

    #create the block Gram matrix
    len_gram = length(mon);
    #declare the PSD matrix variable 
    Gram0 = JuMP.@variable(model, [1:len_gram, 1:len_gram], PSD);
    sigma0 = mon'*Gram0*mon;

    @constraint(model, coefficients(q-psatz_term - sigma0)==0);
    # model, info_seal = add_psatz!(model, q-psatz_term, vars_flat, [], [],order);
    blocksize = [blocksize; len_gram]
    # @constraint(model, coefficients)
    
    # quad_u = 0
    return psatz(mu, quad_x, quad_u, blocksize);
end

function mu_eq_con(T, k, mu, A)
    #the equality constraint involved in data-consistency
    #(assuming no process noise)
    #process noise will be added later
    mu_con = 0*mu[1].s;
    if k<=T-1
        mu_con = mu_con + A'*mu[k].s;
    end
    if k>=2
        mu_con = mu_con - mu[k-1].s;
    end
    return mu_con;
end