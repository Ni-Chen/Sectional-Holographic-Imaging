import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.io import loadmat



#z = np.linspace(0, 1, 100)
#x = z * np.sin(20 * z)
#y = z * np.cos(20 * z)
#
#c = x + y
#
#ax.scatter(x, y, z, c=c)
#
#ax.plot(x, y, z, '-b')


volume = loadmat('../data/SNUE_3d')['obj3d']

print(volume.shape)
print(volume)

x = 1-np.arange(volume.shape[1])[:, None, None]
y = 1-np.arange(volume.shape[0])[None, :, None]
z = 1-np.arange(volume.shape[2])[None, None, :]
x, y, z = np.broadcast_arrays(x, y, z)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.scatter(x, y, z, cmap='viridis', linewidth=0.5, alpha=0.05);
#color_bar.set_alpha(1)
#color_bar.draw_all()
#ax.plot(x, y, z, 'gray')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');

#plt.xlabel("x") #xlabel : 设置X轴的文字
#plt.ylabel("y") #ylabel : 设置Y轴的文字
#plt.title("PyPlot scatter Example") 

fig.savefig('test.pdf', transparent=True)