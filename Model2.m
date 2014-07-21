%Matrices
    %Parameters
        %Outputs
%maximum expected reproductive success between t and T given state
f = zeros(101, 60);
%maximum expected reproductive success between t and T from visiting a
%specifed patch (1 through 3)
v = zeros(3, 101, 60);
%optimal patch to visit given state and time
d = zeros(101, 60);

        %Inputs
%Chance of mortality if patch i is visited
m = [0.05 0.04 0];
%Probability that food is found in patch
p = [0.2 0.5 0];
%Energetic cost of foraging in patch
a = [3 3 3];
%Energetic value of food in patch
y = [20 8 0];


%critical level of reserves
x_crit=0;
%maximum possible level of reserves
x_max=100;
%maximum length of breeding season
t_max=60;



        %Displays
disp ('Parameters are these:');
disp (['x_crit = ' num2str(x_crit) ', x_max = ' num2str(x_max) ', t_max = ' num2str(t_max)]);
for patch = 1:3
    disp (['For patch ' num2str(patch) ': m = ' num2str(m(patch)) ', p = ' num2str(p(patch)) ', a = ' num2str(a(patch)) ', y = ' num2str(y(patch))])
end
%Alternate way to do the same thing - sprintf
 %sprintf - at the beginning you specify the format of you string.
    %Previously we merely typed in what we wanted displayed. However
    %when you want to display a sequence of values you can enter
    %placeholders and refer to what they should contain later on in the
    %string function. E.g. %d is an integer, %f is a real number, and
    % \t is a tab. Numbers before f specify format in terms of decimal
    % places
for patch = 1:3
   disp (sprintf('For patch %d: m = %1.2f, p = %1.1f, a = %d, y = %d', patch, m(patch), p(patch), a(patch), y(patch)));
end

%End Conditions

%In the case of the end conditions,
%Phi is just the reproductive fitness on day T depending on your
%current levels of energy, x.
%A (acap) is the asymptotic value of phi(x,T)
%This is a limit, it is the upper bound of Phi given x and T.
%This end condition models that food is a satiable resource. Not linear
%(Squashed) diminishing returns.
acap = 200;
x_0 = 0.25 * x_max;

disp 'End conditions at T (60)'
disp (sprintf('x\tF(x,T)'))
%Because MATLAB starts from 1 (1 indexed arrays) and
%BASIC uses 0, when we access x we actually neeed to start
%at 1 so that the language understands where we want to write
%our data (correct row/column in the matrix).
for x = x_crit:x_max;
    x_surplus = x - x_crit;
    phi_x_T = acap * x_surplus / (x_surplus + x_0);
    f(x+1, t_max) = phi_x_T;
    disp (sprintf('%d\t%f', x, f(x+1, t_max)));
end

    %Constraints

    %MATLAB has you put zeros in when specifying space to put data.
    %As such this next step is pointless as we are just inserting zeroes on top of zeroes.
for t = 1:t_max-1;
    f(x_crit+1, t) = 0;
end

 
for t = t_max - 1:-1:1;
    disp (sprintf('t\tx\tF(x,t)\td(x,t)\tV(1,x,t)\tV(2,x,t)\tV(3,x,t)'))
    for x = x_crit+1:x_max;
        
        for i = 1:3;
            
            %Effect of visiting a foraging patch
            if (i == 1 || i == 2);
                foodyes = min(x - a(i) + y(i), x_max);
                foodno = max(x - a(i), x_crit);
                
                %Interpolation for non integers
                if rem(foodyes,1) == 0;
                    foodyes_f = f(foodyes+1, t+1);
                else
                    fl = floor(foodyes);
                    ce = ceil(foodyes);
                    w = foodyes - fl;
                    foodyes_f = ((1-w) * f(fl+1, t+1)) + (w * f(ce+1, t+1));
                end
                if rem(foodno,1) == 0;
                    foodno_f = f(foodno+1, t+1);
                else
                    fl = floor(foodno);
                    ce = ceil(foodno);
                    w = foodno - fl;
                    foodno_f = ((1-w) * f(fl+1, t+1)) + ((w * f(ce+1, t+1)));
                end
                
                %Patch Fitness
                v(i,x+1,t) = (1-m(i)) * ((p(i) * foodyes_f) + ((1-p(i)) * foodno_f));
                
                
            else
            %Effect of visiting a refuge
                refuge = max(x-a(3),x_crit);
                
                %Interpolation
                if rem(refuge,1) == 0;
                    refuge_f = f(refuge+1, t+1);
                else
                    fl = floor(refuge);
                    ce = ceil(refuge);
                    w = refuge - fl;
                    refuge_f = ((1-w) * f(fl+1, t+1)) + (w * f(ce+1, t+1));
                end
                    
                %Patch Fitness
                v(3,x+1,t) = (1-m(3)) * refuge_f; 
            end
        end
    
        vmax = max([v(1,x+1,t), v(2,x+1,t), v(3,x+1,t)]);
        %Record optimal patch choice
        if vmax == v(1,x+1,t)
            d(x+1, t) = 1;
        elseif vmax == v(2,x+1, t)
            d(x+1, t) = 2;
        elseif vmax == v(3,x+1, t)
            d(x+1, t) = 3;
        end
        %Record optimal fitness value
    	f(x+1, t) = vmax;
    
        disp(sprintf('%d\t%d\t%2.2f\t%d\t\t%2.2f\t\t%2.2f\t\t%2.2f', t, x, f(x+1,t), d(x+1,t), v(1,x+1,t), v(2,x+1,t), v(3,x+1,t)))
    end
end

%save patch matrix for graphing
cd('E:\Users\Jaggerous\Documents\MATLAB\Project\plot_csv\plot_csv')
csvwrite('mort_0.4_0.04.csv', d)


%Forward iteration probability distribution function

%Set up matrix and parameters
f_prob = zeros(101, 61);
init_state = 50;

%Input initial state
f_prob(init_state+1, 1) = 1;

for t = 1:t_max-1;
    %identify non zero states (possible states) within a particular time step
    possible_states = find(f_prob(:, t));
    %for each state within possible_states do x
    for state_i = transpose(possible_states);
        
        %extracts particular probabilities for state/time from matrix which we want to manipulate
        state_prob = f_prob(state_i,t);
        if state_i == 1;
            
            f_prob(state_i,t+1) = f_prob(state_i,t+1) + state_prob;
        else
            decision = d(state_i, t);
            m_prob = m(decision);
     %Chance of death from visiting patch adjusts state = 0 probability
            f_prob(1,t+1) = f_prob(1,t+1) + (state_prob * m_prob);
            
            p_success = p(decision);
            p_fail = 1 - p(decision);
            reward = y(decision);
            cost = a(decision);
     %Change prob for t+1 based on chance of finding food
            success_i = max(1, state_i + reward - cost);
            success_i = min(101,success_i);
            change = (state_prob * (p_success * (1 - m_prob)));
            
            %Interpolation
            if rem(success_i,1) == 0;
                f_prob(success_i,t+1) = f_prob(success_i,t+1) + change;   
            else
                fl = floor(success_i);
                ce = ceil(success_i);
                w = success_i - fl;
                    
                %Split change between two states at t+1
                f_prob(fl,t+1) = f_prob(fl,t+1) + ((1-w) * change);
                f_prob(ce,t+1) = f_prob(ce,t+1) + ((w) * change);
            end
            
     %Change prob for t+1 based on chance of not food
            fail_i = max(1, state_i - cost);
            fail_i = min(101,fail_i);
            change = (state_prob * (p_fail * (1 - m_prob)));
            
             if rem(fail_i,1) == 0;
                f_prob(fail_i,t+1) = f_prob(fail_i,t+1) + change;   
            else
                fl = floor(fail_i);
                ce = ceil(fail_i);
                w = fail_i - fl;
                    
                %Split change between two states at t+1
                f_prob(fl,t+1) = f_prob(fl,t+1) + ((1-w) * change);
                f_prob(ce,t+1) = f_prob(ce,t+1) + ((w) * change);
            end
            
        end
    end
end

%Extract mean state for forward iteration
mean_state = zeros(1,59);
    for t = 2:60;
        for i = 2:101;
            mean_state(1,t-1)= mean_state(1,t-1) + (f_prob(i,t)*(i-1));
        end
    mean_state(1,t-1)= mean_state(1,t-1)/(1-(f_prob(1,t)));
    end

%Calculate proportion of time spent in each patch
%name empty matrix
patch_time = zeros (3,58);
    %for times you arent dead
     for t = 2:59;
         %identify non-zero probabilities excluding dead state
        for state = transpose(find(f_prob(2:end,t)));
            %adjust state as excluding dead state minuses 1
            state = state + 1;
            patch = d(state, t);
            %increment patch_time accumulator
            patch_time(patch, t-1) = patch_time(patch, t-1) + f_prob(state, t);
            %normalise across non-dead states
        end
        for patch = 1:3;
            patch_time(patch, t-1)=patch_time(patch, t-1)/(1-(f_prob(1,t)));
        end
     
     end
        