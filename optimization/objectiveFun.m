function obj = objectiveFun(params,target_population,alpha,beta,current_params)

   global errorValues;
    if isempty(errorValues)
        errorValues = [];
    end

    growth_rate = params(1);
    death_rate = params(2);
    initialTumorSize=params(3);
    vssl_bdary_value=params(4);
    
    %Run ABM from jar file
    system( 'java -jar /Users/4477116/Documents/projects/Ploidy_abm/Ploidy_abm.jar')

    % Load updated population from file
    %filePath = ['/Users/4477116/Documents/projects/Ploidy_abm/output/',num2str(current_params.tEnd), '_karyotype_Location.csv'];
     filePath = ['/Users/4477116/Documents/projects/Ploidy_abm/output/',num2str(1000),'_karyotype_Location.csv'];
    current_population = readCurrentPop(filePath);
    
    %Using standard deviation instead of pop size
    pop=readtable('/Users/4477116/Documents/projects/Ploidy_abm/output/PopSize.csv');
    pop=pop.Population;

   

     
    % Calculate population deviation
    abm_pop = size(current_population, 1);
    target_pop = size(target_population, 1);
    populationDeviation = abs(abm_pop - target_pop)/(40000-target_pop); %Scaling: Dividing by the maximum possible difference

    %Constraints
     rate_penalty = max(0, abs(growth_rate - death_rate) - 0.015); % Constrain the difference between the rates to within  reasonable range

      min_pop=(target_pop-100);
      max_pop=(target_pop+100);

     if abm_pop < min_pop
        pop_penalty =  (abs(min_pop - abm_pop))/min_pop;
     elseif abm_pop > max_pop
        pop_penalty  = (abs(abm_pop - max_pop))/abm_pop;
     else
         pop_penalty=0;
    end
    

     % Calculate location deviation using the K nearest neigbor search
    [Idx, D] = knnsearch(target_population, current_population);
    locationDeviation=sum(D)/9306;
  

    % Calculate combined objective
     n=current_params.initialTumorSize;

     %current=struct2cell(current_params);
     %valuesArray = cell2mat(current);
     disp(['populationDeviation ', num2str(populationDeviation)])
     disp(['pop_penalty ',num2str(pop_penalty)])
     disp([ 'rate_penalty ',num2str(rate_penalty )])
     disp(['locationDeviation ',num2str(locationDeviation)])

    obj =  std(pop); %(alpha * populationDeviation) %+ pop_penalty + rate_penalty  %+ (beta * locationDeviation); 


     %Saving the error at each iteration
    errorValues=[errorValues,obj];
    save('errorValues.mat', 'errorValues');
 

end