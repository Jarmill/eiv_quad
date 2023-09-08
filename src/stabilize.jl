#stabilize ()
function quad_ss(data, order)

    #create the optimizer
    model = Model(optimizer_with_attributes(Mosek.Optimizer));


        
    #solve the SOS problem
    optimize!(model)

    #recover the solution, process the output
    # output =  slice_recover(model, vars, poly, order, opts, info);

    return output
end


function slice_run(order, opts)
    #slice_max: the main routine
    #
    #Inputs:
    #   order:  the polynomial order of the slicing SOS programs
    #   opts:   options (slice_interface.jl/slice_options)
    #Outputs:
    #   slice_out

    #create the optimizer
    model = Model(optimizer_with_attributes(Mosek.Optimizer));

    #create the variables and polynomials
    vars = make_slice_vars(opts);
    poly = make_slice_poly!(model, order, vars, opts);

    #form the SOS program
    info = slice_program!(model, order, poly, vars, opts)

    #form the objective 
    slice_obj = slice_objective!(model, poly, opts);        
    

    #solve the SOS problem
    optimize!(model)
end