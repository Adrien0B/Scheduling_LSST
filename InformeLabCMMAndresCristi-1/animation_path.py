import matplotlib.pyplot as plt
from matplotlib import animation

import healpy as hpy
import numpy as np


def animation_path(sched,rot,title,*args, **kwargs):
    vid_title = "videos/"+title + ".mp4"
    hpy.mollview(title=title,rot=rot)
    hpy.graticule(verbose=0)
    points = hpy.projscatter(np.pi/2-sched[:,1],sched[:,0],lonlat=False,cmap="hsv",c=np.arange(np.size(sched,0)),*args, **kwargs)
    offto = points.get_offsets()
    arrayto = points.get_array()
    climto = points.get_clim()
    cmapto = points.get_cmap()

    plt.clf()
    hpy.mollview(title=title,rot=rot)
    hpy.graticule(verbose=0)
    points = hpy.projscatter([],[],lonlat=False,*args, **kwargs)
    fig = points.get_figure()
    points.set_clim(climto)
    points.set_cmap(cmapto)

    def animate(i,data,color):
        start_data = max(0,i-300)
        # start_data = 0
        points.set_offsets(data[start_data:i,:])
        points.set_array(color[start_data:i])
        return points,

    anim = animation.FuncAnimation(fig, animate,fargs=(offto,arrayto),frames=np.size(offto,0), interval=20, blit=True)
    anim.save(vid_title)
    # plt.show()
