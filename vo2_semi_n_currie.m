function result=vo2_semi_n_currie(x)
x=x*10^9;
data=[3.99069767441860e+002	2.72527472527473e+000
4.25872093023256e+002	2.72527472527473e+000
4.56395348837209e+002	2.68131868131868e+000
5.56686046511628e+002	2.51648351648352e+000
6.22093023255814e+002	2.45054945054945e+000
7.04941860465116e+002	2.40659340659341e+000
7.92151162790698e+002	2.38461538461538e+000
9.01162790697674e+002	2.40659340659341e+000
9.66569767441860e+002	2.42857142857143e+000
1.05813953488372e+003	2.47252747252747e+000
1.26308139534884e+003	2.54945054945055e+000
1.68168604651163e+003	2.64835164835165e+000
2.04796511627907e+003	2.68131868131868e+000
2.20058139534884e+003	2.68131868131868e+000
];

result=interp1(data(:,1),data(:,2),x);