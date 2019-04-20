%% LOAD FAKE DATA

clc; close all; clear
load hospital
DATA = dataset2table(hospital(:,[6 4 2 3]));
DATA.Properties.RowNames={};
DATA.BP = DATA.BloodPressure(:,1);
DATA.BP(DATA.Weight<140) = DATA.BP(DATA.Weight<140).*.8;
DATA.BP(DATA.Weight>190) = DATA.BP(DATA.Weight>190).*1.2;
DATA.BP(DATA.Sex=='Male') = DATA.BP(DATA.Sex=='Male').*1.2;
DATA.BP(DATA.Age>40) = DATA.BP(DATA.Age>40).*.9;
DATA.Sex = grp2idx(DATA.Sex).*2-3;
DATA.BP = round(DATA.BP);
DAT = DATA(:,[5 2 3 4]);

clearvars -except DAT

disp(DAT(1:10,:)) 


% Previewing the DATA table we just created...
% 
% 
%   DATA =
%     100×4 table
%     BP     Weight     Sex      Age
%     ___    ______    ______    ___
%     149    176       Male      38 
%     118    163       Male      43 
%     100    131       Female    38 
%      94    133       Female    40 
%      88    119       Female    49 
%     109    142       Female    46 
%     130    142       Female    33 
%     138    180       Male      40 
%     138    183       Male      28 
%     ...
%     123    177       Male      48 



%% =================== Part 3: Cost and Gradient descent ===================

data = [DAT.BP DAT.Weight DAT.Age];
y    = DAT.Sex;


X       = [ones(size(data,1), 1), data];  % Add a column of ones to x
theta   = zeros(size(X,2), 1);            % initialize fit params


% COMPUTE COST
m = length(y); % number of training examples

J = sum((X * theta - y).^2) ./ (2*m);

disp(J);


% SET NUMBER OF L1 NEURONS EQUAL TO NUMBER OF VARIANTS
L1neurons  = size(X,2);
nLabels = numel(unique(y));


% NEURAL NET PARAMETERS
%----------------------------------------------------------------------
lambda = 0.005;                 % .001 - .01 is a good operational window
epsInit = 0.22;                 % random initial theta weights
maxIters  = 50;                 % 20-50 iterations should be sufficient
L2neurons = 3;                  % 10 - 50 neurons should be sufficient
%----------------------------------------------------------------------



% ESTABLISH NEURON ACTIVATION FUNCTION
sig    = @(z) 1 ./ (1 + exp(-z) );     % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));    % sigmoid gradient function
sstf = @(n) 2 ./ (1 + exp(-2*n)) - 1;  % sigmoid symetric transfer func
% subplot(2,1,1), plot(sstf(-5:.1:5))
% subplot(2,1,2), plot(sig(-5:.1:5))





% PREP NEURAL NET SYNAPTIC WEIGHTS (THETAS) & GRADIENT DECENT COST FUNC
%----------------------------------------------------------------------

% INITIALIZE THETA WEIGHTS
initTheta1 = rand(L2neurons, L1neurons) * 2 * epsInit - epsInit;
initTheta2 = rand(nLabels, L2neurons) * 2 * epsInit - epsInit;


% UNROLL THETA WEIGHTS 
initial_Thetas = [initTheta1(:) ; initTheta2(:)];
nn_Thetas = initial_Thetas;


% ESTABLISH COST FUNCTION
J = NNCostF(nn_Thetas, L1neurons, L2neurons, nLabels, TRAINX, TRAINL, lambda);




























%% TUTORIAL
%{

Scripps hospital just invested in some fancy new machine learning tech 
that takes a pictures people as they walk in the front door, and it 
automatically estimates their weight, age, and sex. Those are good bits
of information to know, but what they really want to know is a person's 
blood pressue (BP) as they walk through the door.

So they contact the Malinow Consulting Group and ask if we could take 
this information and somehow provide accurate BP estimates. Set to task,
the first thing we do is acquire a dataset of the BP, Weight, Age, and
Sex for 100 people. We want to see if these three IVs are related to BP. 

MATLAB has a build-in dataset just like this. Here I'll give some code 
to load and modify this dataset so it ultimatly contains just the info 
we are interested in (raw values slightly adjusted for demo reasons)...


----------------------------------------------------------------
     Command Window
----------------------------------------------------------------

    clear
    load hospital
    DATA = dataset2table(hospital(:,[6 4 2 3]))
    DATA.Properties.RowNames={}
    DATA.BP = DATA.BloodPressure(:,1)
    DATA.BP(DATA.Weight<140) = DATA.BP(DATA.Weight<140).*.8
    DATA.BP(DATA.Weight>190) = DATA.BP(DATA.Weight>190).*1.2
    DATA.BP(DATA.Sex=='Male') = DATA.BP(DATA.Sex=='Male').*1.2;
    DATA.BP(DATA.Age>40) = DATA.BP(DATA.Age>40).*.9
    DATA.Sex = grp2idx(DATA.Sex).*2-3;
    DATA.BP = round(DATA.BP)
    DATA = DATA(:,[5 2 3 4])

    >>
----------------------------------------------------------------

Previewing the DATA table we just created...


  DATA =
    100×4 table
    BP     Weight     Sex      Age
    ___    ______    ______    ___
    149    176       Male      38 
    118    163       Male      43 
    100    131       Female    38 
     94    133       Female    40 
     88    119       Female    49 
    109    142       Female    46 
    130    142       Female    33 
    138    180       Male      40 
    138    183       Male      28 
    ...
    123    177       Male      48 



One stats method that could inform us of whether there's any predictive 
value in those IV parameters is a linear regression model. Ultimately 
though, the hospital wants a piece of machine learning software they can 
integrate into their current system, so it can get better over time.

With those marching orders, we start building a simple regression-based
machine learning tool that uses gradient descent to minimize a cost
function.




############################################################
STEP-1: CREATE A TRAINING MATRIX WITH IV DATA
############################################################

To create this matrix is quite simple. We just need to put the IV data
for our 100 participants in a matrix such that participants are in rows
and each IV, (weight, age, sex) is a column.

There is just one weird thing we need to do though. We'll put a column
of all ones as the first column of this matrix, followed by the three
actual IVs. This column of ones is called a 'bias'. Its main purpose
is to allow a regression estimator, and later a neural net, to center
the intercept or move it around as it sees fit. It also let's you be a
little sloppy with your data prep. Really you should always convert your 
raw data to z-scores before model building (or at least mean-center the 
data at zero). But as long as there is a bias column of all ones...
whatever. Anyway, let's just go ahead and make that matrix:


    X  =  [ones(100,1) DATA.Weight DATA.Sex DATA.Age];



The other thing we will need is a matrix of parameter coefficients for
fitting the linear model. Since we the input matrix we just created
has four columns (1 bias + 3 IVs), we'll need four theta coefficients:


    t = [0; 0; 0; 0]



Now let's define a cost function. The primary purpose of a cost function
is to capture some information about the error in the DV estimates.
Something like the variance or standard deviation would probably do just 
fine. However, since those terms don't technically apply here, nor are we
bounded or obligated to use them, we'll simply define a cost function 
in the spirit of representing mean error...

    y = DATA.BP
    n = length(y)


    J = sum((X * t - y).^2) ./ (2*n)

For comparison sake, here is the equation for variance:

    v = sum((y_pred - mean(y)).^2) ./ (n-1)


Where y_pred is the predicted value of y based on the GLM.

Ok we got our trainin matrix 'X', our theta coefficient matrix 't',
and our cost function 'J'. What next?




############################################################
STEP-2: RE-LEARN LINEAR REGRESSION AND MATRIX ALGEBRA
############################################################


Here we're building a regression-based machine learning tool. Naturally,
we need to define a linear model. Something like...


    BP = Weight*w + Age*a + Sex*s


Yeah. That's basically it. Except, we can compute all the BP estimates
in a single step if we use matrix multiplication instead:

   
    BP = X*t


For fun, let's remind ourselves what this matrix multiplication is doing
using a simple example. Say that we have two matrices 'A' and 'B'. 

A: 4x3 matrix
B: 3x1 matrix

Matrix multiplication requires that the number of columns in A are equal
to the number of rows in B. So indeed, we can do A*B; but not B*A. It
is easy to visualize whether a matrix multiplication will be valid
simply by putting the dimensions of the two matrices side-by-side, in
the order you want to perform the multiplication.

    A*B  4x3 3x1    4x[3] <--OK--> [3]x1
    B*A  3x1 4x3    3x[1] <--NO--> [4]x3

If the two inner dims are the same, it will work. If not, it won't.
Thus, by design, we created a 2x1 'theta' matrix above, because our
data matrix 'X' is Nx2... X*theta... Nx2 2x1... valid.

What will the dimensions be for the result?...

    C = A*B

That's also easy to determine when looking at these dimensions 
side-by-side. The output matrix will be the same size as the 
outer dimensions of A and B:

    A*B  4x3 3x1    [4]x3 3x[1]

Thus C is going to be a 4x1 matrix. Now that we've confirmed A*B is a 
valid matrix multiplication, and we should expect output will be 4x1
let's actually run through it by hand. Given...


    A = [0 1 2;
         1 1 1;
         1 5 1;
         2 2 4]

    B = [1;
         2;
         3]


--------------------------------
 Command Window
--------------------------------

 >> C = A*B
 C =
      8
      6
     14
     18

 >>
--------------------------------



Great, it worked, and output a 4x1 matrix. Let's do it by hand! Reminder,
A matrix multiplication the dot product of the rows and columns...


        [0 1 2;     [1;
         1 1 1;  *   2;
         1 5 1;      3]
         2 2 4]


       sum( [0 1 2] .* [1 2 3] )  %= [ 8 ;
       sum( [1 1 1] .* [1 2 3] )  %=   6 ;
       sum( [1 5 1] .* [1 2 3] )  %=  14 ;
       sum( [2 2 4] .* [1 2 3] )  %=  18 ]


Let's do it once more, this time B is a 3x2 matrix...


    A = [0 1 2;
         1 1 1;
         1 5 1;
         2 2 4]

    B = [1  0;
         2 -1;
         3  2]


        [0 1 2;     [1  0;
         1 1 1;  *   2  2;
         1 5 1;      3 -1]
         2 2 4]



Again, we know this is still valid...  4x3 3x2  ...with the main
difference being the dimensions of the output matrix. Instead of 
being 4x1 it will be 4x2. In fact, since all the matrix data is still
the same, sans the additional column in B, the first column of the
output matrix will be exactly the same as before:

       sum( [0 1 2] .* [1 2 3] )  %= [ 8 __ ;
       sum( [1 1 1] .* [1 2 3] )  %=   6 __ ;
       sum( [1 5 1] .* [1 2 3] )  %=  14 __ ;
       sum( [2 2 4] .* [1 2 3] )  %=  18 __ ]

We just need to compute the output for the second column:


       sum( [0 1 2] .* [0 2 -1] )  %= [__  3 ;
       sum( [1 1 1] .* [0 2 -1] )  %=  __  1 ;
       sum( [1 5 1] .* [0 2 -1] )  %=  __ -3 ;
       sum( [2 2 4] .* [0 2 -1] )  %=  __  6 ]


Thus, all together, C = A*B equals...

----------------------------------------------------------------
 Command Window
----------------------------------------------------------------

>> C = A*B
C =
     8     3
     6     1
    14    -3
    18     6

>>
----------------------------------------------------------------


That should clarify the low-level nuts and bolts of matrix multiplication.



Now we can turn to *why* are we doing it? At the top of this doc,
I wrote...


% Right. Currently the theta coefficients are set to zeros:
%
%    theta = [0; 0] 
%
% So in a linear model like...
%
%    y = X*theta
%
% all the predictions for y are going to be zero, if the theta 
% coefficents are all zero. Naturally this is going to result 
% in a huge amount of error...


To expand on that, let's say a hospital just invested in some fancy new
machine learning tech that takes a pictures people as they walk in the
front door, and automatically estimates of their weight, age, and sex. 
However, what they'd really like to know is a person's blood pressue (BP). 
So they contact the Malinow Consulting Group and ask if we could take 
this information and provide accurate BP estimates. The first thing
we could do is take a dataset of known BPs Weights Ages and Sexes and
create a linear model, like so...


    BP = Weight*w + Age*a + Sex*s


Computing the model's predictions for say 100 participants would require
multiplying each of the 100 Weights by 'w', each 100 Ages by 'a', etc.
This is tedious. Thankfully, matrix multiplication lets us do all these
things, simultaniously. Moreover it's likely that the dataset will 
already be in some spreadsheet or table format such that each of the 
100 people are in rows, and their data for each of those variables 
(weight, age, sex) are each column. As circumstance would have it, 
MATLAB has a build-in dataset just like this. Here I'll give some code 
to load and modify this dataset so it ultimatly contains just the info 
we are interested in, and adjust the values for demo purposes...


----------------------------------------------------------------
     Command Window
----------------------------------------------------------------

    clear
    load hospital
    DATA = dataset2table(hospital(:,[6 4 2 3]))
    DATA.Properties.RowNames={}
    DATA.BP = DATA.BloodPressure(:,1)
    DATA.BP(DATA.Weight<140) = DATA.BP(DATA.Weight<140).*.8
    DATA.BP(DATA.Weight>190) = DATA.BP(DATA.Weight>190).*1.2
    DATA.BP(DATA.Sex=='Male') = DATA.BP(DATA.Sex=='Male').*1.2;
    DATA.BP(DATA.Age>40) = DATA.BP(DATA.Age>40).*.9
    DATA.Sex = grp2idx(DATA.Sex).*2-3;
    DATA.BP = round(DATA.BP)
    DATA = DATA(:,[5 2 3 4])

    >>
----------------------------------------------------------------



Previewing the DATA table we just created...


  DATA =
    100×4 table
    BP     Weight     Sex      Age
    ___    ______    ______    ___
    149    176       Male      38 
    118    163       Male      43 
    100    131       Female    38 
     94    133       Female    40 
     88    119       Female    49 
    109    142       Female    46 
    130    142       Female    33 
    138    180       Male      40 
    138    183       Male      28 
    ...
    123    177       Male      48 


Everything should be self explainetory. (n.b. I would typically 
mean-center each IV before fitting a linear model, but it's
not a big deal to just keep them as raw values for this tutorial)

Anyway, back to the utility and purpose of matrix multiplication. So we
have our model:

    BP = Weight*w + Age*a + Sex*s

And we want it to make estimates for BP. Knowing that we could use matrix
multiplication to get quick esitmates, we might just try to guess a set
of theta coefficients for t=[w;a;s]. Let's try it...


----------------------------------------------------------------
     Command Window
----------------------------------------------------------------

    % Step-1:  separate the DV from the IVs.

    y = DATA.BP;
    X = [DATA.Weight, DATA.Age, DATA.Sex]


    % Step-2: make some random guesses for the theta coefficients:

    t = [1.1; 0.5; -4]


    % Step-3: matrix multiplication:

    h = X*t



    h =
        208.6
        196.8
        167.1
        170.3
        159.4
        183.2
        176.7

    >>
----------------------------------------------------------------


Ok so that worked. I mean, the matrix multiplication worked. But what 
about the quality of the predictions? We can assess the quality by 
examining how close the BP estimates (eBP) were to the actual values 
for BP. Specifically, we want the error, so we can calculate the 
sum of squared errors (sse):


----------------------------------------------------------------
     Command Window
----------------------------------------------------------------

    err = h - y

    sse = sum(err.^2)

>>
----------------------------------------------------------------


Why?.... Ta-daaa!

    m = length(y)

    J = sum((X * t - y).^2) ./ (2*m)


We just knocked out the main portion of the cost function. Put together
into a single step, we did...


    >> sse = sum((X * t - y).^2)


Now all that's left is to divide by 2*m (m being the number of people):


    >> J = sse ./ ( 2 * m)


Done. The cost of that 'round' of BP parameter estimates was:


    2471.1 = J


HOW CAN WE GET BETTER ESTIMATES OF THE THETA PARAMETER ESTIMATES?

That is what machine learning is all about! 

In practice, there are two ways:

1. normal equation
2. gradient descient


The normal equation allows one to arrive at the global minima in a single
step. Holyfuckingshit. The problem is that the normal equation is
computationally intense, and cannot be applied in all situations. It is
however, elegant. It will give you the coefficient for all thetas in a
single vectorized function. It's good to have a theta that accounts for
the intercept so before equation, I'll just add a column of ones to X:


x = [ones(length(X),1) X]

t = pinv(x' * x) * (x' * y)


Then our prediction for each y is simply given as:

h = sum( x .* t' ,2)



We can assess how good the predictions of this experimental model are
compared to a null model that simply predicts the mean of y for every
sample. Using statistical convention, we will compute the sum of squared
errors, or deviations between the prediction and the actual value of y:


SSE_hypo = sum(  (h - y).^2  )

SSE_null = sum(  (mean(y) - y).^2  )



From those sums of squared error we can actually calculate an F-statistic.
Before I give that equation though, we first need a few house-keeping 
values that will be input into the equation:


n = length(y)             % 100
df_hypo = numel(t)        % 4
df_null = numel(mean(y))  % 1


SSR = SSE_null - SSE_hypo



Where SSR represents the sums of squares reduced by the hypothesized model.
Now we have all we need to compute the F-statistic:



F = (SSR / (df_hypo - df_null)) / (SSE_hypo / (n-df_hypo))



Simplifying by plugging-in n and df values:

F = (SSR / 3) / (SSE_hypo / 96)     % 172.48  (that's a huge F!)



Something often lost in translation, that in my opinion is worth 
pointing out, is what the F-statistic actually represents. Here we have
in front of us an absolutely perfect situation to learn (or remind 
ourselves)...

F is a ratio. 

    (SSR/3) : (SSE_hypo/96)

But what does that ratio represent? At a high-level, we tried to find a
better way to estimate y based on some set of parameters we computed using
the dataset itself. Any time we use the sample data to compute a parameter
estimate (e.g. mean(y) ) those are considered 'burned' degrees of freedom,
or 'df'. We burned 4 degrees of freedom (of 100 available) to make these
estimates (i.e. we computed 4 thetas). The null model used just one, the 
mean of y. So the left-hand of the ratio represents the amount of SSE we 
got rid of, per estimated parameter above the null model. The right-hand 
side of the ratio represents the amount of SSE that still remains, per 
remaining degree of freedom. In essense, this ratio tells us how good our 
model parameters are compared to some random junk parameter estimated from 
the dataset. The ratio ends up being something like:

    21435 : 124


meaning that on average our parameters accounted for 21435 squared errors,
while any old junk parameter would accout for about 124 squared errors.
So, yea, I think we did pretty good choosing those independent variables
as predictors of the dependent variable.













SSN = sum((y0 - Y).^2);
SSA = sum((y1 - Y).^2);
SSR = (SSN - SSA);

r = SSR / SSN;

dfN = 1;
dfA = 3;


Fa = (SSR / (dfA - dfN)) / (SSA / (n-dfA));

Pa = fcdf(Fa,(dfA - dfN),(n-dfA),'upper');




lm = fitlm(DATA,'BP ~ Weight + Sex + Age')








Y     = Acceleration - mean(Acceleration);
Xa    = Model_Year   - mean(Model_Year);
Xb    = Horsepower   - mean(Horsepower);

%---
n = numel(Y);

X = [ones(n,1) Xa Xb];
XX = (X' * X)^-1;
XY = X' * Y;
b = XX * XY;
disp(b);


B0 = b(1);
B1 = b(2);
B2 = b(3);

y0 = B0;

y1 = B0 + B1.*Xa + B2.*Xb;





















To calculate the predictions for weight, for each
person, all at the same time is simply...

    weight = X*t

given...

    X is 100x3
    t is   3x1

    w = X*t

    verify input dims:   100x[3]  [3]x1

    verify output dims:  [100]x3  3x[1]


Everything looks good. 

        height  age  sex
        [72     22    1;     [1;
         45      6    1;  *   2;
         60     23   -1;      3]
         47      6   -1]


       sum( [0 1 2] .* [1 2 3] )  %= [ 8 ;
       sum( [1 1 1] .* [1 2 3] )  %=   6 ;
       sum( [1 5 1] .* [1 2 3] )  %=  14 ;
       sum( [2 2 4] .* [1 2 3] )  %=  18 ]

clear
load hospital
hospital


%}



%% LINEAR REGRESSION NOTES
%{
clc; close all; clear;
load carsmall

% Acceleration = B0 + B1.*Model_Year + B2.*Horsepower
Acceleration(isnan(Acceleration))  = nanmean(Acceleration);
Model_Year(isnan(Model_Year))      = nanmean(Model_Year);
Horsepower(isnan(Horsepower))      = nanmean(Horsepower);

Y     = Acceleration - mean(Acceleration);
Xa    = Model_Year   - mean(Model_Year);
Xb    = Horsepower   - mean(Horsepower);

%---
n = numel(Y);

X = [ones(n,1) Xa Xb];
XX = (X' * X)^-1;
XY = X' * Y;
b = XX * XY;
disp(b);


B0 = b(1);
B1 = b(2);
B2 = b(3);

y0 = B0;

y1 = B0 + B1.*Xa + B2.*Xb;


SSN = sum((y0 - Y).^2);
SSA = sum((y1 - Y).^2);
SSR = (SSN - SSA);

r = SSR / SSN;

dfN = 1;
dfA = 3;


Fa = (SSR / (dfA - dfN)) / (SSA / (n-dfA));

Pa = fcdf(Fa,(dfA - dfN),(n-dfA),'upper');



fprintf('\n     Beta       SSE      Fstat       pValue     \n')
fprintf('B0 %8.3f %8.0f %5s %11s  \n',B0,SSN,'-','-')
fprintf('B1 %8.3f %8.0f %10.3f %14.4g  \n',B1,SSA,Fa,Pa)
fprintf('B2 %8.3f %8.0f %10.3f %14.4g  \n',B2,SSA,Fa,Pa)
fprintf('r = %5.4f  \n\n',r)
disp(XY)


%}



% compute and display initial cost
J = computeCost(X, y, theta);

fprintf('With theta = [0 ; 0]\nCost computed = %f\n', J);
fprintf('Expected cost value (approx) 32.07\n');



% further testing of the cost function
J = computeCost(X, y, [-1 ; 2]);


fprintf('\nWith theta = [-1 ; 2]\nCost computed = %f\n', J);
fprintf('Expected cost value (approx) 54.24\n');




fprintf('Program paused. Press enter to continue.\n');
fprintf('\nRunning Gradient Descent ...\n')

[theta, Jh] = gradientDescent(X, y, theta, alpha, iterations);



% print theta to screen
fprintf('Theta found by gradient descent:\n');
fprintf('%f\n', theta);
fprintf('Expected theta values (approx)\n');
fprintf(' -3.6303\n  1.1664\n\n');

    % Plot the linear fit
    hold on; % keep previous plot visible
plot(X(:,2), X*theta, '-')
    legend('Training data', 'Linear regression')
    hold off % don't overlay any more plots on this figure

% Predict values for population sizes of 35,000 and 70,000

predict1 = [1, 3.5] *theta;

fprintf('For population = 35,000, we predict a profit of %f\n',predict1*10000);

predict2 = [1, 7] * theta;

fprintf('For population = 70,000, we predict a profit of %f\n',predict2*10000);


%% ============= Part 4: Visualizing J(theta_0, theta_1) =============
fprintf('Visualizing J(theta_0, theta_1) ...\n')

% Grid over which we will calculate J
theta0_vals = linspace(-10, 10, 100);
theta1_vals = linspace(-1, 4, 100);

% initialize J_vals to a matrix of 0's
J_vals = zeros(length(theta0_vals), length(theta1_vals));

% Fill out J_vals
for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];
	  J_vals(i,j) = computeCost(X, y, t);
    end
end


% Because of the way meshgrids work in the surf command, we need to
% transpose J_vals before calling surf, or else the axes will be flipped
J_vals = J_vals';
% Surface plot
figure;
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');

% Contour plot
figure;
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
contour(theta0_vals, theta1_vals, J_vals, logspace(-2, 3, 20))
xlabel('\theta_0'); ylabel('\theta_1');
hold on;
plot(theta(1), theta(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);



































%% Machine Learning Online Class
%  Exercise 1: Linear regression with multiple variables
%
%  Instructions
%  ------------
% 
%  This file contains code that helps you get started on the
%  linear regression exercise. 
%
%  You will need to complete the following functions in this 
%  exericse:
%
%     warmUpExercise.m
%     plotData.m
%     gradientDescent.m
%     computeCost.m
%     gradientDescentMulti.m
%     computeCostMulti.m
%     featureNormalize.m
%     normalEqn.m
%
%  For this part of the exercise, you will need to change some
%  parts of the code below for various experiments (e.g., changing
%  learning rates).
%
%% ================ Part 1: Feature Normalization ================

% Clear and Close Figures
clear ; close all; clc

fprintf('Loading data ...\n');

% Load Data
data = load('ex1data2.txt');
X = data(:, 1:2);
y = data(:, 3);
m = length(y);


% Scale features and set them to zero mean
fprintf('Normalizing Features ...\n');

[X, mu, sigma] = featureNormalize(X);

% Add intercept term to X
X = [ones(m, 1) X];


%% ================ Part 2: Gradient Descent ================

% ====================== YOUR CODE HERE ======================
% Instructions: We have provided you with the following starter
%               code that runs gradient descent with a particular
%               learning rate (alpha). 
%
%               Your task is to first make sure that your functions - 
%               computeCost and gradientDescent already work with 
%               this starter code and support multiple variables.
%
%               After that, try running gradient descent with 
%               different values of alpha and see which one gives
%               you the best result.
%
%               Finally, you should complete the code at the end
%               to predict the price of a 1650 sq-ft, 3 br house.
%
% Hint: By using the 'hold on' command, you can plot multiple
%       graphs on the same figure.
%
% Hint: At prediction, make sure you do the same feature normalization.


fprintf('Running gradient descent ...\n');


nloops = 400;
num_iters = 400;

ai = linspace(.001,.03,nloops);
Ji = zeros(num_iters, nloops);
Ti = zeros(3, nloops);

% Choose some alpha value
for nn = 1:nloops

    %alpha = .01;
    alpha = ai(nn);
    
    % Init Theta and Run Gradient Descent 
    theta = zeros(3, 1);
    [theta, J_history] = gradientDescentMulti(X, y, theta, alpha, num_iters);
    
    Ti(:,nn) = theta;
    Ji(:,nn) = J_history;

end


close all; clc;
fh1=figure('Units','normalized','OuterPosition',[.02 .1 .95 .8],'Color','w','MenuBar','none'); 
ax1 = axes('Position',[.05 .1 .40 .8],'Color','none');
    xlabel('Number of iterations');
    ylabel('Cost J');
    hold on;
ax2 = axes('Position',[.52 .1 .40 .8],'Color','none');
    xlabel('\theta_1'); 
    ylabel('\theta_2');
    zlabel('J');
    view(50,20); grid on;
    hold on;

axes(ax1) 
plot(  repmat([1:size(Ji,1)]',1,nloops)  , Ji, '-b', 'LineWidth', 2);
pause(.2)

axes(ax2)
surf(Ti(2,:), Ti(3,:), Ji')


% Display gradient descent's result
[MinJ , MinThetaInd] = min(min(Ji));

fprintf('\nLowest cost value found in test loop: %5.0f \n', MinThetaInd);
fprintf('Using alpha value: %13.3f \n', ai(MinThetaInd));

fprintf('\nTheta computed from gradient descent (on Z scaled data): \n')
fprintf('\n    %.3f ', Ti(:,MinThetaInd));
fprintf('\n\n')

Tnorm = pinv(X' * X) * X' * y;

fprintf('\nTheta computed from normal equation (on Z scaled data): \n')
fprintf('\n    %.3f ', Tnorm);
fprintf('\n\n')


theta = Ti(:,MinThetaInd);

% Estimate the price of a 1650 sq-ft, 3 br house
% ====================== YOUR CODE HERE ======================

x = ([1650 , 3] - mu) ./ sigma;

X = [1 x];

price = X * theta;

% ============================================================

fprintf(['\nPredicted price of a 1650 sq-ft, 3 br house ' ...
         '(using gradient descent):\n\n %4.s $%.2f\n\n\n'],' ', price);



%% ================ Part 3: Normal Equations ================

% ====================== YOUR CODE HERE ======================
% Instructions: The following code computes the closed form 
%               solution for linear regression using the normal
%               equations. You should complete the code in 
%               normalEqn.m
%
%               After doing so, you should complete this code 
%               to predict the price of a 1650 sq-ft, 3 br house.
%

%% Load Data
data = csvread('ex1data2.txt');
X = data(:, 1:2);
y = data(:, 3);
m = length(y);

% Add intercept term to X
X = [ones(m, 1) X];

% Calculate the parameters from the normal equation
theta = normalEqn(X, y);

% Display normal equation's result
fprintf('\nTheta computed from the normal equations (on raw data): \n');
fprintf(' %12.2f \n', theta);
fprintf('\n');


% Estimate the price of a 1650 sq-ft, 3 br house
% ====================== YOUR CODE HERE ======================

X = [1 , 1650 , 3];

price = X * theta;

% ============================================================

fprintf(['\nPredicted price of a 1650 sq-ft, 3 br house ' ...
         '(using normal equations):\n\n %4.s $%.2f\n\n\n'],' ', price);







%################################################################
%%  PERFORM MACHINE LEARNING USING HOMEBREW METHODS
%################################################################
% LOGISTIC REGRESSION NEURAL NET WITH GRADIENT DECENT

Niters = 20;
RES = zeros(Niters,5);

% szTx = size(TRAINX,1);
% nP = fliplr(round(linspace(100,szTx,Niters)))';



for nn = 1:Niters


clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
VLOCI VCASE VCTRL VTRCASE VTRCTRL VTECASE VTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN nP


% The *TRNN* matrix is formatted as follows
%
% ROW-1: VLOCI.VID of variant for TRNN(1,4:end)
% COL-1: PHEN.SRR participant ID
% COL-2: PHEN.AD participant AD CASE\CTRL status
% COL-3: bias input node 'ones'
%
% TRNN(2:end,3:end) ... the bias + the variant data
% TRNN(2:end,4:end) ... just the variant data



% DEAL OUT VARIANTS AND LABELS FOR TRAINING DATA
ADNN    =  TRNN;
TRAINX  =  ADNN(2:end,4:end);  
TRAINL  =  ADNN(2:end,2) + 1;


% DEAL OUT VARIANTS AND LABELS FOR TESTING DATA
TESTX   =  TENN(2:end,4:end);
TESTL   =  TENN(2:end,2) + 1;



% RANDOMIZE ROWS
randp        = randperm(size(TRAINX,1));
TRAINX       = TRAINX(randp,:);
TRAINL       = TRAINL(randp,:);


% % UNCOMMENT FOR SERIAL REDUCTION
% TRAINX = TRAINX(  1:nP(nn) ,:);
% TRAINL = TRAINL(  1:nP(nn) ,:);



% SET NUMBER OF L1 NEURONS EQUAL TO NUMBER OF VARIANTS
L1neurons  = size(TRAINX,2);
nLabels = 2;



% NEURAL NET PARAMETERS
%----------------------------------------------------------------------
lambda = 0.005;                 % .001 - .01 is a good operational window
epsInit = 0.22;                 % random initial theta weights
maxIters  = 50;                 % 20-50 iterations should be sufficient
L2neurons = 35;                 % 10 - 50 neurons should be sufficient
%L2neurons = 20 + randi(20);    % 10 - 50 neurons should be sufficient
%maxIters  = 50 + randi(20);    % 20-50 iterations should be sufficient
%----------------------------------------------------------------------



% ESTABLISH NEURON ACTIVATION FUNCTION
sig    = @(z) 1 ./ (1 + exp(-z) );     % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));    % sigmoid gradient function
sstf = @(n) 2 ./ (1 + exp(-2*n)) - 1;  % sigmoid symetric transfer func
% subplot(2,1,1), plot(sstf(-5:.1:5))
% subplot(2,1,2), plot(sig(-5:.1:5))



clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
VLOCI VCASE VCTRL VTRCASE VTRCTRL VTECASE VTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN ADNN TRAINX TRAINL TESTX TESTL...
nLabels maxIters lambda L1neurons L2neurons epsInit sig sigrad sstf nP





% PREP NEURAL NET SYNAPTIC WEIGHTS (THETAS) & GRADIENT DECENT COST FUNC
%----------------------------------------------------------------------

% INITIALIZE THETA WEIGHTS
initTheta1 = rand(L2neurons, L1neurons+1) * 2 * epsInit - epsInit;
initTheta2 = rand(nLabels, L2neurons+1) * 2 * epsInit - epsInit;


% UNROLL THETA WEIGHTS 
initial_Thetas = [initTheta1(:) ; initTheta2(:)];
nn_Thetas = initial_Thetas;


% ESTABLISH COST FUNCTION
J = NNCostF(nn_Thetas, L1neurons, L2neurons, nLabels, TRAINX, TRAINL, lambda);








% TRAIN NEURAL NETWORK USING TRAINING DATASET
%----------------------------------------------------------------------
options = optimset('MaxIter', maxIters);

GcostFun = @(T) NNCostF(T, L1neurons, L2neurons, nLabels, TRAINX, TRAINL, lambda);

[nn_Thetas, cost] = NNfmincg(GcostFun, initial_Thetas, options);
%----------------------------------------------------------------------




% REROLL THETA WEIGHTS
Theta1 = reshape(nn_Thetas(1:L2neurons * (L1neurons + 1)), L2neurons, (L1neurons + 1));

Theta2 = reshape(nn_Thetas((1 + (L2neurons * (L1neurons + 1))):end), nLabels, (L2neurons + 1));




clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
VLOCI VCASE VCTRL VTRCASE VTRCTRL VTECASE VTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN ADNN TRAINX TRAINL TESTX TESTL...
nLabels maxIters lambda L1neurons L2neurons epsInit sig sigrad sstf nP...
Theta1 Theta2




% EVALUATE PERFORMANCE ON **TRAINING** SET
%----------------------------------------------------------------------

[p , a , h] = NNpredict(Theta1, Theta2, TRAINX);

% p : prediction label {1,2}
% a : confidence level of p
% h : confidence level of p and min label



TRAINPCTCORRECT = mean(p == TRAINL);

disp('Percent accuracy on training data:')
disp( TRAINPCTCORRECT )





% EVALUATE NEURAL NETWORK PERFORMANCE ON **TEST** SET
%----------------------------------------------------------------------

[p , a , h] = NNpredict(Theta1, Theta2, TESTX);


TESTPCTCORRECT = mean(p == TESTL);
disp('Percent accuracy on testing data:')
disp( TESTPCTCORRECT )






% EXAMINE & REPORT ACCURACY FOR HIGH CONFIDENCE DATA
%----------------------------------------------------------------------
MidConfThresh = .70;
HiConfThresh  = .85;


y = TESTL;

numCorrect = p == y;

MidConf      = a >= MidConfThresh & a < HiConfThresh;
MidConfPCT   = (mean(p(MidConf) == y(MidConf)))*100;
MidConfNp    = numel(p(MidConf));
MidConfPCTp  = (numel(p(MidConf)) / numel(p))*100;

MidHiConf      = a >= MidConfThresh;
MidHiConfPCT   = (mean(p(MidHiConf) == y(MidHiConf)))*100;
MidHiConfNp    = numel(p(MidHiConf));
MidHiConfPCTp  = (numel(p(MidHiConf)) / numel(p))*100;

HiConf       = a >= HiConfThresh;                   % Logical index of high confidence predictions
HiConfPCT    = (mean(p(HiConf) == y(HiConf)))*100;  % Percent correct hi conf predictions
HiConfNp     = numel(p(HiConf));                    % Total number of hi conf predictions
HiConfPCTp   = (numel(p(HiConf)) / numel(p))*100;   % Percent of hi conf predictions
HiConfNCorr  = HiConfNp * (HiConfPCT / 100);        % Total number of correct hi conf predictions


fprintf('\nPercent accuracy on MID conf (%.2f < a < %.2f) testing data: >>> %.1f%% <<<\n',...
    MidConfThresh, HiConfThresh, MidConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] Subset:Total patients) \n\n',...
    MidConfNp , numel(p) , MidConfPCTp )


fprintf('\nPercent accuracy on MID+HI conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    MidConfThresh, MidHiConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] Subset:Total patients) \n\n',...
    MidHiConfNp , numel(p) , MidHiConfPCTp )


fprintf('\nPercent accuracy on HIGH conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    HiConfThresh, HiConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] Subset:Total patients) \n\n',...
    HiConfNp , numel(p) , HiConfPCTp )





RES(nn,:) = [TESTPCTCORRECT MidHiConfPCT MidHiConfPCTp HiConfPCT HiConfPCTp];

clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
VLOCI VCASE VCTRL VTRCASE VTRCTRL VTECASE VTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN ADNN TRAINX TRAINL TESTX TESTL...
nLabels maxIters lambda L1neurons L2neurons epsInit sig sigrad sstf nP...
Theta1 Theta2

end

clc;
RES(:,1) = RES(:,1) .*100;
disp(RES)
disp(RES(:,2:end))





%% HOW DOES NN SCORE FAKE PERSON WITH 1 SINGLE VARIANT... IN THE TOP GENE
% THE TOP GENE IS USUALLY APOE


sig    = @(z) 1 ./ (1 + exp(-z) );      % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));     % sigmoid gradient function
sstf   = @(n) 2 ./ (1 + exp(-2*n)) - 1; % sigmoid symetric transfer func

m = size(TRAINX, 1);
h1 = sig([ones(m, 1) TRAINX] * Theta1');
h2 = sig([ones(m, 1) h1] * Theta2');
[a, p] = max(h2, [], 2);

TRAINX(1,:)
h1(1,:)
h2(1,:)

FAKEPEOPLE = eye(size(TRAINX,2));

m = size(FAKEPEOPLE, 1);
h1 = sig([ones(m, 1) FAKEPEOPLE] * Theta1');
h2 = sig([ones(m, 1) h1] * Theta2');
[a, p] = max(h2, [], 2);

VLOCI.NNactiv  = h2(:,2);
VLOCI.NNguess  = p-1;
NNEYE = VLOCI(:,[7 26 27]);
disp(NNEYE(1:10,:))