from kbutil.plotting import pylab_pretty_plot
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from numpy import array,median


def plot_R_delta(self,R,deltaij):
    '''
    Makes an R/delta (Reproducbility and its fluctuations) plot that lets you
    determine the dimension of the sparse subspace.  The median R and deltaij
    is a horizontal bar, and the actual samples are plotted as points.
    '''
    # median Reproducibility and deltaij
    medR = median(R,axis=0)
    meddij = median(deltaij,axis=0)
    signum = array(range(1,len(medR)+1))
    nSig = len(medR)

    fig,ax = plt.subplots()
    # sizing optional?
    fig.set_figwidth(10)
    fig.set_figheight(5)
    # plot
    pylab_pretty_plot()
    ax.plot(signum.T + 0.1, R.T,'r.', label = r'$R_{ij}$', alpha = 0.5, mec = 'none')
    ax.plot(signum.T-0.1, deltaij.T,'k.', label = r'$\delta_{ij}$', alpha = 0.5,mec = 'none')
    ax.scatter(signum.T + 0.1, medR.T, s = 400, c = 'r', marker = '_', edgecolor = 'r', alpha = 1., linewidth = 2, rasterized = True)
    ax.scatter(signum.T - 0.1, meddij.T, s = 400, c = 'k', marker = '_', edgecolor = 'k', alpha = 1., linewidth = 2, rasterized = True)

    # legend stuff
    red_patch = mpatches.Patch(ls = '-',color='red', label=r'$R_{ij}$', alpha = 0.7)
    black_patch = mpatches.Patch(ls = '-',color = 'black', label = r'$\delta_{ij}$', alpha = 0.7)
    plt.legend(handles=[red_patch,black_patch], bbox_to_anchor = (.31,0.63), prop = {'size':25})

    # additional formatting
    ax.set_xlabel(r'Independent components', size = 25)
    ax.set_ylabel(r'$\delta_{ij}$ and $R_{ij}$', size = 25)
    ax.set_xlim([0,nSig+1])
    ax.set_ylim([-0.2,1.2])
    plt.yticks([0,0.5,1])
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_tick_params(width=2)

    return fig
    #repro_mean = np.zeros((1,nSig))
    #repro_std = np.zeros((1,nSig))
    #repro_median = np.zeros((1,nSig))
    #for i in xrange(0,nSig):
    #    repro_mean[0,i] = np.mean(true_repro[:,i])
    #    repro_std[0,i] = np.std(true_repro[:,i])
    #    repro_median[0,i] = np.median(true_repro[:,i])

    #delta_ij = np.zeros((nSamp,nSig))
    #for i in xrange(0,nSamp):
    #    for j in xrange(0,nSig):
    #        delta_ij[i][j] = np.abs(true_repro[i][j]-true_repro[0][j])

    #delta_mean = np.zeros((1,nSig))
    #delta_std = np.zeros((1,nSig))
    #delta_median = np.zeros((1,nSig))

    #for i in xrange(0,nSig):
    #    delta_mean[0,i] = np.mean(delta_ij[:,i])
    #    delta_std[0,i] = np.std(delta_ij[:,i])
    #    delta_median[0,i] = np.median(delta_ij[:,i])
    #sig_no = np.linspace(1,nSig,nSig)

#     fig,ax = plt.subplots()
#     fig.set_figwidth(10)
#     fig.set_figheight(5)
#     kbplt.pylab_pretty_plot(lines = 4, width =4, size = 8, labelsize =25, markersize = 13, fontsize = 32, usetex = True)
#     ax.plot(sig_no.T + 0.1, true_repro.T,'r.', label = r'$R_{ij}$', alpha = 0.5, mec = 'none')
#     ax.plot(sig_no.T-0.1, delta_ij.T,'k.', label = r'$\delta_{ij}$', alpha = 0.5,mec = 'none')
#     ax.scatter(sig_no.T + 0.1, repro_median.T, s = 400, c = 'r', marker = '_', edgecolor = 'r', alpha = 1., linewidth = 2, rasterized = True)
#     ax.scatter(sig_no.T - 0.1, delta_median.T, s = 400, c = 'k', marker = '_', edgecolor = 'k', alpha = 1., linewidth = 2, rasterized = True)
#
#     red_patch = mpatches.Patch(ls = '-',color='red', label=r'$R_{ij}$', alpha = 0.7)
#     black_patch = mpatches.Patch(ls = '-',color = 'black', label = r'$\delta_{ij}$', alpha = 0.7)
#     plt.legend(handles=[red_patch,black_patch], bbox_to_anchor = (.31,0.63), prop = {'size':25})
#
#     ax.set_xlabel(r'Independent components', size = 25)
#     ax.set_ylabel(r'$\delta_{ij}$ and $R_{ij}$', size = 25)
#     ax.set_xlim([0,nSig+1])
#     ax.set_ylim([-0.2,1.2])
#     plt.yticks([0,0.5,1])
#     ax.xaxis.set_major_locator(plt.NullLocator())
#     ax.yaxis.set_tick_params(width=2)
#
#     fig.savefig('Fisher_delta_repro.pdf', bbox_inches = 'tight')
#
#     fig1, axes1 = plt.subplots(ncols=4, figsize = plt.figaspect(1.2/4))
#     kbplt.pylab_pretty_plot(lines = 4, width =4, size = 8, labelsize =25, markersize = 13, fontsize = 32, usetex = True)
#     for axes , i in zip(axes1,[0,1,2,3]):
#         kbplt.plot_hist(results[0]['raic_comp'][i,:].T, nbins = 15, ax = axes, kde = True)
# #    axes.set_xlim([-3.5,3.5])
# #    axes.set_ylim([-0.1,1.1])
#         axes.set_xlabel(r'$R_{%i}$'%(i+1), size = 25)
#         axes.xaxis.set_major_locator(plt.NullLocator())
#         axes.yaxis.set_major_locator(plt.NullLocator())
#     plt.tight_layout()
#     fig1.savefig('Fisher_hist.pdf', bbox_inches = 'tight')
#     plt.close('all')
#
