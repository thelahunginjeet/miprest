import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from numpy import array,median

def pylab_pretty_plot(lines=2,width=3,size=4,labelsize=16,markersize=10,fontsize=20,lfontsize=16,lframeon=False,usetex=True):
    """
    Changes pylab plot defaults to get nicer plots - frame size, marker size, etc.

    Parameters:
    ------------
    lines      : linewidth
    width      : width of framelines and tickmarks
    size       : tick mark length
    labelsize  : font size of ticklabels
    markersize : size of plotting markers
    fontsize   : size of font for axes labels
    lfontsize  : legend fontsize
    usetex     : use latex for labels/text?

    """
    pylab.rc("lines",linewidth=lines)
    pylab.rc("lines",markeredgewidth=size/3)
    pylab.rc("lines",markersize=markersize)
    pylab.rc("ytick",labelsize=labelsize)
    pylab.rc("ytick.major",pad=size)
    pylab.rc("ytick.minor",pad=size)
    pylab.rc("ytick.major",size=size*1.8)
    pylab.rc("ytick.minor",size=size)
    pylab.rc("xtick",labelsize=labelsize)
    pylab.rc("xtick.major",pad=size)
    pylab.rc("xtick.minor",pad=size)
    pylab.rc("xtick.major",size=size*1.8)
    pylab.rc("xtick.minor",size=size)
    pylab.rc("axes",linewidth=width)
    pylab.rc("text",usetex=usetex)
    pylab.rc("font",size=fontsize)
    pylab.rc("legend",fontsize=lfontsize)
    pylab.rc("legend",frameon=lframeon)


def plot_R_delta(R,deltaij):
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
