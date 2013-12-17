function v=nfkbBasal(v,flag)
  
%% Call the ODE function to reset the persistent variables
nfkbOde([],[],[],v);

%% Run Phase 1 (equilibrium phase)
screen('Phase 1',v.D_FLAG);
v.PHASE     = 1;
static      = false;    % Set to 'true' after reaching equilibrium
count       = 1;        % Iteration Counter
threshold   = 1;        % Max % difference used by evaluate_phase1()
initial_values = v.INITVALUES{1};
v.BASAL_VALUES = [];

while ~static % Iterate through Phase 1 until equilibrium is reached
    options = odeset('RelTol', 1e-4);
    [t1, r1] = ode15s('nfkbOde', [v.START_TIME 0], initial_values,options,v);
%     plot(t1,r1(:,44));hold on ; 
    % Evaluate results and return true if at equilibrium
    
    if (evaluatePhase1(r1, threshold)) % values have converged
        v.START_TIME    = (count * v.START_TIME) - (count-1);
        static          = true;
        v.BASAL_VALUES  = r1(end,:);
        screen(['Equilibrium met at ' num2str(v.START_TIME) ...
            ' min'],v.D_FLAG);
    elseif count > 100 % values are not converging so stop running
        screen('===> max phase 1 steps reached (100 steps)',v.D_FLAG);
        static          = true;
        v.START_TIME    = (count * v.START_TIME) - (count-1);
    else % equilibrium not met, so run another round of Phase 1
        count           = count + 1;
        initial_values  = r1(end,:);
    end
end

%% subsubroutine static
    function static = evaluatePhase1(results,threshold)
        % threshold= percentage
        
        prev_values     = abs(results(1,:)'); %init value
        current_values  = abs(results(end,:)'); %end value
        
        max = current_values * (1 + (threshold/100));
        min = current_values * (1 - (threshold/100));
        
        
        for i = 1:length(prev_values) % examine all the state valuables (SVs).
            if ( (prev_values(i) >= max(i)) || (prev_values(i) <= min(i)) )
                if (current_values(i) > 1e-21)  % Resolves noise errors
                    static = false; % make sure all SVs statisfy the condition
                    return;
                end
            end
        end
        static = true;
    end
end