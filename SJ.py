import numpy as np
from scipy.spatial.distance import cdist

def sparse_jump(Y, n_states, max_features, jump_penalty=1e-5,
                max_iter=10, tol=1e-4, n_init=10, verbose=False):
    # Implementation of sparse jump model
    n_obs, n_features = Y.shape
    max_features = np.clip(max_features, a_min=1, a_max=np.sqrt(n_features)) #define upper bound for feature-weights
    feat_w = np.repeat(1 / np.sqrt(n_features), n_features) #init feature weights
    states = None

    for it in range(max_iter):
        states = jump(Y * np.sqrt(feat_w), # jump model is fitted on weighted features
                      n_states,
                      initial_states=states,
                      jump_penalty=jump_penalty,
                      n_init=n_init)
        if len(np.unique(states)) == 1:
            break
        else:
            new_w = get_weights(Y, states, max_features, n_states) # compute new feature weights
        if abs(new_w - feat_w).sum() / abs(feat_w).sum() < tol: # check if converge is reached
            break
        elif verbose:
            print('Iteration {}, w diff {:.6e}'.format(it, abs(new_w - feat_w).sum()))
        feat_w = new_w

    return states, feat_w


def jump(Y, n_states, jump_penalty=1e-5, initial_states=None,
         max_iter=10, n_init=10, tol=None, verbose=False):
    # Fit jump model using framework of Bemporad et al. (2018)

    # state initialization, can be provided as input or computed via init_states
    if initial_states is not None:
        initial_states = np.array(initial_states, dtype=np.int64)
        if len(np.unique(initial_states)) == n_states:
            s = initial_states.copy()
        else:
            s = init_states(Y, n_states)
    else:
        s = init_states(Y, n_states)

    n_obs, n_features = Y.shape
    Gamma = jump_penalty * (1 - np.eye(n_states)) 
    best_loss = None
    best_s = None
    # FC: counter for the number of iterations (both max_iter and n_init)
    #n_maxiter = 1
    #n_ninit = 1
    #

    for init in range(n_init):
        mu = np.zeros((n_states, n_features)) # mu is a matrix of zeros, one value for each feat and for each state
        loss_old = 1e10
        for it in range(max_iter):
            #
            #
            #FC: This part can be changed according to different loss functions##
            for i in np.unique(s):
                # Fit model by updating mean of observed states
                mu[i] = np.mean(Y[s==i], axis=0) # mu is a vector, one for each state
            # Fit state sequence
            s_old = s.copy()
            loss_by_state = cdist(mu, Y, 'euclidean').T**2
            #FC: #################################################################
            #
            #
            V = loss_by_state.copy()
            #FC: compute V for each t in T,...,1
            for t in range(n_obs-1, 0, -1):
                V[t-1] = loss_by_state[t-1] + (V[t] + Gamma).min(axis=1)
            #FC: compute state sequence as argmin(V+lambda*I)
            s[0] = V[0].argmin()
            for t in range(1, n_obs):
                s[t] = (Gamma[s[t-1]] + V[t]).argmin()
            #FC: Monitor convergence
            if len(np.unique(s)) == 1:
                break
            loss = min(V[0])
            if verbose:
                print('Iteration {}: {:.6e}'.format(it, loss))
            if tol:
                epsilon = loss_old - loss 
                if epsilon < tol:
                    break
            #FC: if s^i=s^i-1 break
            elif np.array_equal(s, s_old):
                break
            loss_old = loss
            #n_maxiter = n_maxiter+1
            
        if (best_s is None) or (loss_old < best_loss):
            best_loss = loss_old
            best_s = s.copy()
        s = init_states(Y, n_states)
        #n_ninit = n_ninit+1

    #return best_s, mu, n_maxiter, n_ninit
    return best_s


def init_states(Y, n_states):
    # Generate initial states using K-means++ (Arthur and Vassilvitskii, 2007) 
    n_obs, n_features = Y.shape
    centers = np.zeros((n_states, n_features))
    center_idx = np.random.randint(n_obs)
    centers[0] = Y[center_idx]
    n_local_trials = 2 + int(np.log(n_states))
    closest_dist_sq = cdist(centers[0, None], Y, 'euclidean')**2
    current_pot = closest_dist_sq.sum()
    
    for i in range(1, n_states):
        rand_vals = np.random.sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(np.cumsum(closest_dist_sq), 
                                        rand_vals)
        distance_to_candidates = cdist(Y[candidate_ids], Y, 'euclidean')**2
        # Decide which candidate is the best
        best_candidate = None
        best_pot = None
        best_dist_sq = None
        for trial in range(n_local_trials):
            # Compute potential when including center candidate
            new_dist_sq = np.minimum(closest_dist_sq,
                                     distance_to_candidates[trial])
            new_pot = new_dist_sq.sum()

            # Store result if it is the best local trial so far
            if (best_candidate is None) or (new_pot < best_pot):
                best_candidate = candidate_ids[trial]
                best_pot = new_pot
                best_dist_sq = new_dist_sq

        centers[i] = Y[best_candidate]
        current_pot = best_pot
        closest_dist_sq = best_dist_sq
    
    # Compute the state assignment
    states = cdist(centers, Y, 'euclidean').argmin(axis=0)
        
    return states


def get_weights(Y, states, max_features, n_states):
    # Find weights given a state sequence by maximizing the interstate distance
    BCSS = get_BCSS(Y, states)
    delta = binary_search(BCSS, max_features)
    w = calc_new_feature_weights(BCSS, delta)

    return w


def get_BCSS(Y, states):
    # Find BCSS given a state sequence
    WCSS = np.zeros(Y.shape[1])
    for i in np.unique(states):
        mask = (states == i)
        if mask.sum() > 1:
            WCSS += np.square(Y[mask] - np.mean(Y[mask], axis=0)).sum(axis=0)
    TSS = np.square(Y - np.mean(Y, axis=0)).sum(axis=0)

    return TSS - WCSS

def get_WCSS(Y, states):
    # Find WCSS given a state sequence
    WCSS = np.zeros(Y.shape[1])
    for i in np.unique(states):
        mask = (states == i)
        if mask.sum() > 1:
            WCSS += np.square(Y[mask] - np.mean(Y[mask], axis=0)).sum(axis=0)
    #TSS = np.square(Y - np.mean(Y, axis=0)).sum(axis=0)

    return WCSS

def get_TCSS(Y, states):
    # Find TCSS given a state sequence
    #WCSS = np.zeros(Y.shape[1])
    #for i in np.unique(states):
     #   mask = (states == i)
      #  if mask.sum() > 1:
       #     WCSS += np.square(Y[mask] - np.mean(Y[mask], axis=0)).sum(axis=0)
    TSS = np.square(Y - np.mean(Y, axis=0)).sum(axis=0)

    return TSS


def binary_search(objective, norm_constraint, max_iter=15):
    l2n_arg = np.linalg.norm(objective)
    if l2n_arg == 0 or abs(objective / l2n_arg).sum() <= norm_constraint:
        return 0
    lam1 = 0
    lam2 = abs(objective).max() - 1e-5
    for iter in range(max_iter):
        su = soft_threshold(objective, (lam1 + lam2) / 2)
        if abs(su / np.linalg.norm(su)).sum() < norm_constraint:
            lam2 = (lam1 + lam2) / 2
        else:
            lam1 = (lam1 + lam2) / 2
        if (lam2 - lam1) < 1e-4:
            break

    return (lam1 + lam2) / 2


def calc_new_feature_weights(objective, delta):
    # Calculate feature weights using soft thresholding
    soft = soft_threshold(objective, delta)
    w = soft / np.linalg.norm(soft)
    
    return w


def soft_threshold(x, delta):
    
    return np.sign(x) * np.maximum(0, np.abs(x) - delta)
