
def gelman_rubin(chain_ensemble):
    """
    returns: Gelman-Rubin scale reduction factor for set of chains
    Parameters:
    chain_ensemble: the ensemble of chains to use; this paramater expects that the first nsteps/2 of each chain have already been discarded;
    chain_ensemble should be an array of the format nchains x nsteps x ndim
    """
    nchains,nsteps,ndim=chain.shape
    
    #calculate the mean of each chain
    mean=np.mean(chain,axis=1)

    #calculate the variance of each chain
    wvar=np.var(chain,axis=1)

    #calculate the mean of the variances of each chain
    W=np.mean(wvar, axis=0)

    #calculate the variance of the chain means multiplied by n
    B=nstep*np.var(mean, axis=0)

    #calculate estimated variance:
    var_est=(1-1/nstep) * W + (1/nstep)*B

    #calculate the potential scale reduction factor
    R=np.sqrt(var_est/W)

    return R

