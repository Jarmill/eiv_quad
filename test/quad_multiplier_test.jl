using eiv_quad
using JuMP

model = Model();

n = 2;
m = 1;
order = 1;

@polyvar A[1:n, 1:n];
@polyvar B[1:n, 1:m];


vs = vars(A, B)
qm = make_mult_quad(model, vs, order);
