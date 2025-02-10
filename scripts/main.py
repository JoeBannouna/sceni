import click
from astrosceni.image import Image
from astrosceni.mufinder import MuFinder
from astrosceni.stars_filter import StarsFilter

@click.command()
@click.option('--nb', type=click.Path(exists=True), required = True, help = 'Narrowband Image Path')
@click.option('--bb', type=click.Path(exists=True), required = True, help = 'Broadband Image Path')

def subtract(nb, bb):
    """
    Loads NB and BB images, performs subtraction via MuFinder, plots result
    """
    #Loads images
    nb_img = Image(nb)
    bb_img = Image(bb)
    #Creates MuFinder
    mufinder = MuFinder(nb_img, bb_img, mu_resolution = 0.05)

    #Get Result images
    images = mufinder.getResultImages()

    #Plot First result image
    images[0].plot(cmap = 'turbo')

if __name__ == '__main__':
    subtract()