import click
from astrosceni.image import Image
from astrosceni.mufinder import MuFinder
from astrosceni.stars_filter import StarsFilter
from matplotlib import pyplot as plt

@click.command()
@click.option('--nb', type=click.Path(exists=True), required = True, help = 'Narrowband Image Path')
@click.option('--bb', type=click.Path(exists=True), required = True, help = 'Broadband Image Path')
@click.option('--x-start', type=int, default=None, help='Starting x coordinate for cropping')
@click.option('--x-end', type=int, default=None, help='Ending x coordinate for cropping')
@click.option('--y-start', type=int, default=None, help='Starting y coordinate for cropping')
@click.option('--y-end', type=int, default=None, help='End y coordinate for cropping')
@click.option('--filter-stars', is_flag=True, help='Removes any visible stars within images')
@click.option('--export', type=click.Path(), default=None, help='Sets export path for plots')

#If running from sceni folder use following as test:
#python3 -m scripts.main --nb data/bs_h_ave_wcs.fits --bb data/bs_r_ave_wcs.fits --x-start=100 --x-end=-100 --y-start=100 --y-end=-100 --filter-stars --export scripts/subtracted_plot.png

def subtract(nb, bb, x_start, x_end, y_start, y_end, filter_stars, export):
    #Loads images
    nb_img = Image(nb)
    bb_img = Image(bb)

    #Crop Images
    nb_img.cropPixels(x_start, x_end, y_start, y_end)
    bb_img.cropPixels(x_start, x_end, y_start, y_end)

    #Saves images into MuFinder, gets optimal mu
    mufinder = MuFinder(nb_img, bb_img, mu_resolution = 0.05)
    result_img = None

    if not filter_stars:        
        #Get Result images
        images = mufinder.getResultImages()

        #Plot result image
        result_img = images[0]

    # If filtering, get optimal mu, remove stars, and then subtract with already found mu
    else:
        #Obtain optimal mus for images with stars
        optimal_mus = mufinder.getOptimalMus()

        #Remove visible stars from each image using hipparcus cataogue
        starsFilter = StarsFilter()
        filtered_nb = starsFilter.filterStars(nb_img)
        filtered_bb = starsFilter.filterStars(bb_img)

        #Subtract filtered images with optimal mu from before
        result_img = Image.subtract(filtered_nb, filtered_bb, optimal_mus[0])

    # If export path given, show plot and save to path
    if export:
        result_img.plot(export_path = export)
    else:
        result_img.plot()        
if __name__ == '__main__':
    subtract()