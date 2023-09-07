mutable struct psatz
    mu      #equality constrained multipliers
    quad_x  #x: polynomials (s, tau) from quad_mult    
    quad_u  #u: polynomials (s, tau) from quad_mult    
end

mutable struct quad_mult
    s       #poly vector inside the SOC
    tau     #poly vector for the SOC radius
    Gram    #Gram matrices in the certificate
end

function make_mult_quad(model, vars, order, SPARSE=true)
    quad_out = 0;
    if SPARSE
        quad_out = make_mult_quad_sparse(model, vars, order);        
    else
        quad_out = make_mult_quad_dense(model, vars, order);
    end

    return quad_out;
end

function make_mult_quad_dense(model, vars, order)
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



        #add the constraints

    return quad_mult(s, tau, Gram)


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

function make_mult_quad_sparse(model, vars, order)
    #make_mult_quad_sparse: create the multipliers associated with a single robust SOC constraint:
    #
    #[sum_i z_it, s_it; s_it, z_it] in SOS^2[A, B]. Each time index t has its own SOC constraint
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
    n = size(vars.A, 2)
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

    return quad_mult(s, tau, Gram)


end

function quad_psatz(q, order, model, data, vars, SPARSE=false)
#enforce that q>=0 over the consistency set described in data
# (assume only x-noise for now)
    #create the multipliers 
    n = size(data.X, 1);
    m = size(data.U, 1);
    T = size(data.U, 2);

        

    vars_flat = vec([vars.A vars.B]);

    #create the residual
    Xnext = data.X[:, 2:end];
    Xprev = data.X[:, 1:end-1];
    U = data.U[:, 1:end-1];
    h0 = Xnext - vars.A*Xprev - vars.B*U;

    #the contribution such that q==psatz_term
    psatz_term = 0;

    #form mu (equality constraint multipliers)
    mu = Array{Polynomial}(undef, n, T-1);
    for k = 1:(T-1)
        for j = 1:n
            mucurr, mucurrc, mucurrb = add_poly!(model, vars_flat, 2*order-1);
            mu[j, k] = mucurr;
        end
        # v, vc, vb = add_poly!(model, x, 2d)
        # psatz_term = psatz_term + mu[:, k]'*h0[k];
    end
    psatz_term = sum(mu .* h0);

    #form quad (the quadratic inequality constraints (s, tau))    
    if data.epsilon[1] > 0
        quad_x = Array{quad_mult}(undef, T, 1);
        for k = 1:T        
            quad_x[k] = make_mult_quad(model, vars, order, SPARSE);        
        end
    else
        quad_x = zeros(T, 1);
    end

    #form quad_u (quadratic inequality constraints for u)
    if data.epsilon[2] > 0
        quad_u = Array{quad_mult}(undef, T-1, 1);
        for k = 1:T-1        
            quad_u[k] = make_mult_quad(model, vars, order, SPARSE);        
        end
    else
        quad_u = zeros(T-1, 1);
    end
    
    #now iterate through the constraints (for x)
    s_sig_x = sqrt(data.Sigma[1]);
    s_sig_u = sqrt(data.Sigma[1]);
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
            mu_con_u = vars.B'*mu[:, k];

            for j = 1:m
                @constraint(model, coefficients.(s_term_u[j] - mu_con_u[j])==0);
            end
            psatz_term = psatz_term + quad_u[k].tau*data.epsilon[2];    
        end
                
    end

    #seal up with the psatz constraint
    model = add_psatz!(model, q-psatz_term, vars_flat, [], [],order);
    # @constraint(model, coefficients)
    
    # quad_u = 0
    return psatz(mu, quad_x, quad_u);
end

function mu_eq_con(T, k, mu, A)
    #the equality constraint involved in data-consistency
    #(assuming no process noise)
    #process noise will be added later
    mu_con = 0*mu[:, 1];
    if k<=T-1
        mu_con = mu_con + A'*mu[:, k];
    end
    if k>=2
        mu_con = mu_con - mu[:, k-1];
    end
    return mu_con;
end