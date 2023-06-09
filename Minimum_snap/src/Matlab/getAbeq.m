 function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    n_poly_perseg = (n_order+1); % coef number of perseg
    %#####################################################
    % p,v,a,j constraint in start, 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    Aeq_start(:,1:n_poly_perseg) = [calc_tvec(0,n_order,0);
                                    calc_tvec(0,n_order,1);
                                    calc_tvec(0,n_order,2);
                                    calc_tvec(0,n_order,3)];
    beq_start = start_cond';
    %
    %
    %
    
    %#####################################################
    % p,v,a,j constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    Aeq_end(:,end-n_order:end) = [calc_tvec(ts(end),n_order,0);
                                  calc_tvec(ts(end),n_order,1);
                                  calc_tvec(ts(end),n_order,2);
                                  calc_tvec(ts(end),n_order,3)];
    beq_end = end_cond';
    %
    %
    %
    
    %#####################################################
    % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    for i = 1:n_seg-1
        Aeq_wp(i,1+(i-1)*(n_order+1):i*(n_order+1)) = calc_tvec(ts(i),n_order,0);
        beq_wp(i) = waypoints(i+1);
    end
    %
    %
    %
    %
    
    %#####################################################
    % position continuity constrain between each 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    for i = 1:n_seg-1
        tvec_p_start = calc_tvec(0,n_order,0);
        tvec_p_end = calc_tvec(ts(i),n_order,0);
        Aeq_con_p(i,1+(i-1)*(n_order+1):(i+1)*(n_order+1)) = [tvec_p_end, -tvec_p_start];
    end
    %
    %
    %
    
    %#####################################################
    % velocity continuity constrain between each 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    for i = 1:n_seg-1
        tvec_v_start = calc_tvec(0,n_order,1);
        tvec_v_end = calc_tvec(ts(i),n_order,1);
        Aeq_con_v(i,1+(i-1)*(n_order+1):(i+1)*(n_order+1)) = [tvec_v_end, -tvec_v_start];
    end
    %
    %
    %

    %#####################################################
    % acceleration continuity constrain between each 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    for i = 1:n_seg-1
        tvec_a_start = calc_tvec(0,n_order,2);
        tvec_a_end = calc_tvec(ts(i),n_order,2);
        Aeq_con_a(i,1+(i-1)*(n_order+1):(i+1)*(n_order+1)) = [tvec_a_end, -tvec_a_start];
    end
    %
    %
    %
    
    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    for i = 1:n_seg-1
        tvec_j_start = calc_tvec(0,n_order,3);
        tvec_j_end = calc_tvec(ts(i),n_order,3);
        Aeq_con_j(i,1+(i-1)*(n_order+1):(i+1)*(n_order+1)) = [tvec_j_end, -tvec_j_start];
    end
    %
    %
    %
    
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end