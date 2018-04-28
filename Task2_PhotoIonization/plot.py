import matplotlib.pyplot as plt

__all__ = ['plt_axis']


def plt_axis(fontsize, grid=False):
    
    ax=plt.gca()
    ax.xaxis.labelpad = 15
    ax.yaxis.labelpad = 15
    
    ax.tick_params('both', length=5, width=2, which='major')
    plt.minorticks_on()
    
    if grid:
        ax.grid(color='silver', linestyle='--')

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
