from scipy.stats import gaussian_kde
import numpy as np

def zFPKM(log_fpkm, resolution=100):
    """
    Should we consider a gene "expressed" or not?

    Implements algorithm from Hart et al 2013 (DOI:10.1186/1471-2164-14-778).
    In that paper, using human data they conclude that zFPKM < -3 is
    a reasonable cutoff since at that point, the frequency of repressive
    chromatin marks outweighs active marks.

    Parameters
    ----------
    log_fpkm : array-like
        FPKM (or TPM) on a log scale

    resolution : int
        Use this many evenly-spaced points between min and max log_fpkm

    Returns
    -------
    Dictionary of values as follows:

        xi: basis over which the kde was evaluated
        yi: the kde
        mu: maximum of the kde
        U: mean of log_fpkm values > mu
        sigma: std of the fitted half-gaussian
        zfpkm: z-transformed FPKM
    """
    kernel = gaussian_kde(log_fpkm)
    xi = np.linspace(log_fpkm.min(), log_fpkm.max(), resolution)
    yi = kernel.evaluate(xi)
    mu = xi[np.argmax(yi)]
    U = log_fpkm[log_fpkm > mu].mean()
    sigma = (U - mu) * np.sqrt(np.pi / 2)
    zlog_fpkm = (log_fpkm - mu) / sigma
    return {
        'xi': xi,
        'yi': yi,
        'mu': mu,
        'U': U,
        'z': zlog_fpkm}

