using JuMP
using DynamicPolynomials
using eiv_quad
using Revise
model = Model();

n = 2;
m = 1;
order = 1;

@polyvar A[1:n, 1:n];
@polyvar B[1:n, 1:m];


vs = sys_vars(A, B)
# qm = make_mult_quad(model, vs, order);
qm = make_mult_quad(n, model, vs, order, false);
qms = make_mult_quad(n, model, vs, order, true);

# qms3 = make_mult_quad(3, model, vs, order, true);