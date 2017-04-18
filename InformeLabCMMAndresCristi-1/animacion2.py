import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class PolarScatAnim(object):
	def __init__(self,Data,colors,SS):
		self.Data=Data
		self.colors=colors
		self.fig=plt.figure()
		self.ax=self.fig.add_subplot(111,projection='polar')
		self.ax.set_rmax(np.pi/3)
		self.sc_anim=animation.FuncAnimation(self.fig,self.update_ScheduleAA,frames=SS-1,fargs=(self.Data,self.colors),init_func=self.setup_plot,interval=100,blit=True,repeat_delay=200)

	def update_ScheduleAA(self,num,dato,colors):
		self.points.set_offsets(dato[:num,:])
		self.points.set_array(colors[:num])
		self.points.set_clim(vmin=0,vmax=1)
		self.points.set_cmap('viridis')
		return self.points,

	def setup_plot(self):
		self.points=self.ax.scatter([],[],edgecolor='none',vmin=0,vmax=1,cmap='viridis',s=50,animated=True)
		self.ax.set_rmax(np.pi/3)
		return self.points,
		
	def show(self):
		plt.show()


if __name__ == '__main__':
	schedAA=np.load('schedAA.npy')
	SS=schedAA.size/2
	print SS
	schedAA[:,1]=np.pi/2-schedAA[:,1]
	colores=np.load('colores.npy')
	A=PolarScatAnim(schedAA,colores,SS)

	A.sc_anim.save("fact_sch_sesgo_TodosT1_1000it.mp4")#, codec='avi')#writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'])

	A.show()
