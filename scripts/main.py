import click
from astrosceni.image import Image
from astrosceni.mufinder import MuFinder
from astrosceni.stars_filter import StarsFilter

@click.command()
@click.option('--nb', type=click.Path(exists=True), required = True, help = 'Narrowband Image Path')
@click.option('--bb', type=click.Path(exists=True), required = True, help = 'Broadband Image Path')
@click.option('--x-start', type=int, default=None, help='Starting x pixel for cropping')
@click.option('--x-end', type=int, default=None, help='Ending x pixel for cropping')
@click.option('--y-start', type=int, default=None, help='Starting y pixel for cropping')
@click.option('--y-end', type=int, default=None, help='End y pixel for cropping')
@click.option('--filter-stars', is_flag=True, help='Removes any visible stars within images')
@click.option('--export', type=click.Path(), default=None, help='Sets export path for plots')

#If running from sceni folder use following as test: 
#python3 -m scripts.main --nb data/bs_h_ave_wcs.fits --bb data/bs_r_ave_wcs.fits --x-start=100 --x-end=-100 --y-start=100 --y-end=-100 --filter-stars --export scripts/subtracted_image.png

def subtract(nb, bb, x_start, x_end, y_start, y_end, filter_stars, export):
    """ Program that performs nb - mu*bb, where nb and bb are narrowband and broadband images respectively and mu is some optimal scaling factor"""
    #Loads images
    try:
        nb_img = Image(nb)
        bb_img = Image(bb)
    except Exception as e:
        click.echo(f"Error loading images: {e}")
        raise click.Abort()

    #Crop Images
    try:
        nb_img.cropPixels(x_start, x_end, y_start, y_end)
        bb_img.cropPixels(x_start, x_end, y_start, y_end)
    except Exception as e:
        click.echo(f"Error cropping images: {e}")
        raise click.Abort()

    #Initialize MuFinder with the cropped images
    try:
        mufinder = MuFinder(nb_img, bb_img, mu_resolution = 0.05)
    except Exception as e:
        click.echo(f"Error initializing MuFinder: {e}")
        raise click.Abort()

    result_img = None

    # Compute result image with or without filtering stars
    if not filter_stars: 
        try:       
            images = mufinder.getResultImages()
            result_img = images[0]
        except Exception as e:
            click.echo(f"Error generating result images: {e}")
            raise click.Abort()

    else:
        #Obtain optimal mu values
        try:
            optimal_mus = mufinder.getOptimalMus()
        except Exception as e:
            click.echo(f"Error computing optimal mu values: {e}")
            raise click.Abort()

        # Filter stars from the images using the star filter
        try:
            starsFilter = StarsFilter()
            # Uncomment line below to use bigger star catalogue
            # starsFilter.setCatalogue(catalogue_id="I/259/tyc2", ra_col_name="RA(ICRS)", dec_col_name="DE(ICRS)", app_mag_col_name="VTmag")
            filtered_nb = starsFilter.filterStars(nb_img)
            filtered_bb = starsFilter.filterStars(bb_img)
        except Exception as e:
            click.echo(f"Error filtering stars: {e}")
            raise click.Abort()

        #Subtract filtered images using optimal mu found earlier
        try:
            result_img = Image.subtract(filtered_nb, filtered_bb, optimal_mus[0])
        except Exception as e:
            click.echo(f"Error subtracting filtered images: {e}")
            raise click.Abort()

    # If export path given, show plot and save to path
    try:
        if export:
            result_img.plot(export_path = export, cmap = 'viridis')
        else:
            result_img.plot(cmap = 'viridis') 
    except Exception as e:
        click.echo(f"Error plotting image: {e}")
        raise click.Abort()
               
if __name__ == '__main__':
    subtract()