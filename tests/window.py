from tkinter import *
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk) 
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import matplotlib.pyplot as plt

hdulist = fits.open('data/rim_Ha_wcs.fits')
# hdulist.info()

# plot function is created for 
# plotting the graph in 
# tkinter window 
def plot(): 

  # the figure that will contain the plot 
  fig = Figure(figsize = (5, 5), 
        dpi = 100) 

  # # list of squares 
  # y = [i**2 for i in range(101)] 

  # # adding the subplot 
  # plot1 = fig.add_subplot(111) 

  wcs = WCS(hdulist[0].header)
  ax = WCSAxes(fig, [0, 0, 1, 1], wcs=wcs)
  # ax = plt.subplot(projection=wcs)
  # fig.add_subplot(ax)

  fig.add_axes(ax)



  ax.set_xlabel('RA')
  ax.set_ylabel('Dec')

  ax.grid(color='black', ls='solid')
  
  # plotting the graph 
  ax.imshow(hdulist[0].data) 

  # creating the Tkinter canvas 
  # containing the Matplotlib figure 
  canvas = FigureCanvasTkAgg(fig, 
              master = window) 
  canvas.draw() 

  # placing the canvas on the Tkinter window 
  canvas.get_tk_widget().pack() 

  # creating the Matplotlib toolbar 
  toolbar = NavigationToolbar2Tk(canvas, 
                window) 
  toolbar.update() 

  # placing the toolbar on the Tkinter window 
  canvas.get_tk_widget().pack() 


# the main Tkinter window 
window = Tk() 

# setting the title 
window.title('Plotting in Tkinter') 

# dimensions of the main window 
window.geometry("500x500") 

# button that displays the plot 
plot_button = Button(master = window, 
					command = plot, 
					height = 2, 
					width = 10, 
					text = "Plot") 

# place the button 
# in main window 
plot_button.pack() 

# run the gui 
window.mainloop() 
